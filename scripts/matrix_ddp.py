__authors__ = "Ragousandirane RADJASANDIRANE"
__contact__ = "radja.ragou@gmail.com"
__date__ = "14/09/2021"
__version__= "1.0"


# Scripts import
from dict_sequence import *
from classes import *
import score_align
import main

# Modules import
import itertools
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, mahalanobis
import time
import multiprocessing as mp
from joblib import Parallel, delayed
import blosum



def matrix_distance(coords_pdb, sequence_pdb_ca):
	matrix = pd.DataFrame(cdist(coords_pdb, coords_pdb,metric="euclidean"))
	matrix.columns = sequence_pdb_ca.seq
	matrix.index = sequence_pdb_ca.seq
	return matrix

def Hmatrix(sequence_target,sequence_pdb_res,sequence_pdb_ca, matrix_dist, dope_score,options):
	n_col = sequence_pdb_res.length + 1
	n_row = sequence_target.length + 1
	print(sequence_target.seq)
	print(sequence_target.seq3)

	if options.blosum:
		if options.weight == -1:
			perc_id = score_align.first_align(sequence_target, sequence_pdb_res)
			weight = score_align.get_weigth(perc_id)
			print("Percent identity : {}%".format(perc_id))
			print("Weights used : {} for dope and {} for blosum\n".format(weight[0],weight[1]))
		else:
			weight = weight_value[options.weight]
	else:
		weight = (1,0)
	options.weight_used = weight

	if options.dssp:
		array_gap = score_align.ss_to_gap(score_align.get_ss(sequence_pdb_res.filename))
	else:
		array_gap = np.zeros((n_col-1))

	matrix_check = do_LM_parallel(sequence_target,sequence_pdb_res,
	sequence_pdb_ca, matrix_dist,array_gap, dope_score,weight,n_row, n_col, options)

	weight_dope, weight_blosum = weight
	matrix_high = np.zeros((n_row,n_col), dtype=np.float64)
	matrix_backtrace = np.empty((n_row,n_col),dtype=str)
	matrix_backtrace[0,1:] = "l"
	matrix_backtrace[1:,0] = "u"
	matrix_backtrace[0,0] = "x"

	for i in range(1, n_row):
		for j in range(1, n_col):
			left = matrix_high[i][j-1] + array_gap[j-1]
			up = matrix_high[i-1][j] + array_gap[j-1]
			best_score = matrix_check[i-1][j-1]
			diag = matrix_high[i-1][j-1]  + (weight_dope*best_score) + (weight_blosum*blosum.blosum62(sequence_target.seq[i-1].res,sequence_pdb_res.seq[j-1].res))
			matrix_high[i][j] =  max(diag, left, up)  
			if matrix_high[i][j] == diag:
				matrix_backtrace[i][j] = "d"
			elif matrix_high[i][j] == left:
				matrix_backtrace[i][j] = "l"
			else:
				matrix_backtrace[i][j] = "u"
	return matrix_high, matrix_backtrace

def do_LM_parallel(sequence_target,sequence_pdb_res,sequence_pdb_ca, matrix_dist, array_gap,dope_score,weight,n_row,n_col,options):
	coords_ij = itertools.product(range(1,n_row), range(1,n_col))
	print("Filling all low matrix in paralell...")
	matrix_check = Parallel(n_jobs= options.cpu, verbose = 0, prefer="processes")(delayed(Lmatrix)
	(sequence_target,sequence_pdb_res,sequence_pdb_ca,matrix_dist, array_gap,dope_score,weight,(fixed_i,fixed_j),
	n_row,n_col) for fixed_i,fixed_j in coords_ij)
	matrix_check = np.array(matrix_check, dtype=np.float64)
	matrix_check = np.reshape(matrix_check, (n_row-1,n_col-1))
	return matrix_check


def Lmatrix(sequence_target,sequence_pdb_res,sequence_pdb_ca, matrix_dist,array_gap, dope_score,weight, residue_fixed, n_row, n_col):

	weight_dope, weight_blosum = weight
	matrix_low = np.full([n_row,n_col], -100, dtype=np.float64)
	matrix_low[0][:] = np.zeros((1,n_col))
	matrix_low[:,0] = np.zeros((1,n_row))
	residue_fixed_i = residue_fixed[0]
	residue_fixed_j = residue_fixed[1]
	for i in range(1,n_row):		
		for j in range(1,n_col):
			if (i,j) == (residue_fixed_i,residue_fixed_j):
				matrix_low[residue_fixed_i][residue_fixed_j] = matrix_low[i-1][j-1]
			if (j < residue_fixed_j and i < residue_fixed_i) or (j > residue_fixed_j and i > residue_fixed_i):
				pairs_residues = [(residue_fixed_i,residue_fixed_j),(i,j)]
				score = score_align.get_score(sequence_target, sequence_pdb_ca, dope_score, matrix_dist, pairs_residues)
				left = matrix_low[i][j-1] + array_gap[j-1]
				up = matrix_low[i-1][j] + array_gap[j-1]
				diag = matrix_low[i-1][j-1] 
				matrix_low[i][j] = (weight_dope*score) + max(diag,left,up)  + (weight_blosum*blosum.blosum62(sequence_target.seq[i-1].res,sequence_pdb_res.seq[j-1].res))

	if residue_fixed_j == n_col-1:
		time_x_line = time.time() - main.GLOBAL_START_TIME
		x = residue_fixed_i
		y = n_row - residue_fixed_i
		time_remaining = y * time_x_line / x
		print("Time remaining (estimated): {:.2f}s".format(time_remaining))

	if matrix_low[-1][-1] == -100:
		return 0
	return matrix_low[-1][-1] 

def get_alignment(sequence_target,sequence_pdb_res,matrix_backtrace):
	j = matrix_backtrace.shape[1] - 1 
	i = matrix_backtrace.shape[0] - 1
	align_res = []
	align_pdb = []
	while i > 0 and j > 0:
		if matrix_backtrace[i][j] == "d":
			align_pdb.append(sequence_pdb_res.seq[j-1].res)
			align_res.append(sequence_target.seq[i-1].res)
			i -= 1
			j -= 1			
		elif matrix_backtrace[i][j] == "u":
			align_pdb.append("-")
			align_res.append(sequence_target.seq[i-1].res)
			i -= 1
		elif matrix_backtrace[i][j] == "l":			
			align_pdb.append(sequence_pdb_res.seq[j-1].res)
			align_res.append("-")
			j -= 1
	return align_pdb, align_res


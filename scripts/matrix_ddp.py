"""All functions for dealing with matrix"""

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
from joblib import Parallel, delayed
import blosum



def matrix_distance(coords_pdb, sequence_pdb_ca):
	"""Compute distance matrix from 3D coordinates of alpha carbon

	Parameters
	----------
    coords_pdb : list
		3D coordinates of alpha carbon

    sequence_pdb_ca : Sequence Object
        Sequence of alpha carbon from PDB file (template)

	Returns
	-------
	matrix_dist : Numpy matrix
		Matrix of pair distance between all alpha carbon from PDB
	"""
	
	matrix_dist = pd.DataFrame(cdist(coords_pdb, coords_pdb,metric="euclidean"))
	matrix_dist.columns = sequence_pdb_ca.seq
	matrix_dist.index = sequence_pdb_ca.seq
	return matrix_dist

def Hmatrix(sequence_target,sequence_pdb_res,sequence_pdb_ca, matrix_dist, dope_score,options):
	"""Fill the High matrix

	Each position of High matrix is a fixed target residue on the template
	and to know if this residue is in a good position, Low matrix is
	filled using dynamic programming considering the fixed residue and computing
	a global score
	This global score from Low matrix is retrieve in High matrix for that particular
	position

	To reduce time execution, all Low matrix are filled in parallel before
	dealing with High matrix

	Parameters
	----------
    sequence_target : Sequence_Residue Object
        Target sequence to align with template

    sequence_pdb_res : Sequence_Residue Object
        Sequence from PDB file (template)

    sequence_pdb_ca : Sequence Object
        Sequence of alpha carbon from PDB file (template)

	matrix_dist : Numpy matrix
		All distances computed between alpha carbon of template

	dope_score : DataFrame
		All dope score in a dataframe

	options : namespace
		all argument gave by user with argparse

	Returns
	-------
	matrix_high : Numpy matrix
		High matrix filled

	matrix_backtrack : Numpy matrix
		Backtracking matrix to find the final alignment
	"""
	
	n_col = sequence_pdb_res.length + 1
	n_row = sequence_target.length + 1
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
	matrix_backtrack = np.empty((n_row,n_col),dtype=str)
	matrix_backtrack[0,1:] = "l"
	matrix_backtrack[1:,0] = "u"
	matrix_backtrack[0,0] = "x"

	for i in range(1, n_row):
		for j in range(1, n_col):
			left = matrix_high[i][j-1] + array_gap[j-1]
			up = matrix_high[i-1][j] + array_gap[j-1]
			best_score = matrix_check[i-1][j-1]
			diag = matrix_high[i-1][j-1]  + (weight_dope*best_score) + (weight_blosum*blosum.blosum62(sequence_target.seq[i-1].res,sequence_pdb_res.seq[j-1].res))
			matrix_high[i][j] =  max(diag, left, up)  
			if matrix_high[i][j] == diag:
				matrix_backtrack[i][j] = "d"
			elif matrix_high[i][j] == left:
				matrix_backtrack[i][j] = "l"
			else:
				matrix_backtrack[i][j] = "u"
	return matrix_high, matrix_backtrack

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
	"""Fill the High matrix

	For each position in High matrix, a full Low matrix is filled considering 
	a fixed residue

	Parameters
	----------
    sequence_target : Sequence_Residue Object
        Target sequence to align with template

    sequence_pdb_res : Sequence_Residue Object
        Sequence from PDB file (template)

    sequence_pdb_ca : Sequence Object
        Sequence of alpha carbon from PDB file (template)

	matrix_dist : Numpy matrix
		All distances computed between alpha carbon of template

	array_gap : Numpy array
		Array of gap score

	dope_score : DataFrame
		All dope score in a dataframe
	
	weight : Tuple
		Weights used 

	residue_fixed : Tuple
		Coordinates of the fixed residue

	n_row : int
		Number of row for Low matrix (same as High matrix)
	
	n_col : int
		Number of column for Low matrix (same as High matrix)

	Returns
	-------
	float
		Best score of Low matrix for a particular fixed residue
	"""

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

def get_alignment(sequence_target,sequence_pdb_res,matrix_backtrack):
	"""Redo the path drawn in the backtrack matrix to get the final alignment
	
	Parameters
	----------
	sequence_target : Sequence_Residue Object
        Target sequence to align with template

    sequence_pdb_res : Sequence_Residue Object
        Sequence from PDB file (template)
	
	matrix_backtrack : Numpy matrix
		Matrix that save all the direction taken to get the best scores
		of high matrix

	Returns
	-------
	align_pdb : list
		Sequence of residue from PDB aligned

	align_res : list
		Sequence of residue from Fasta aligned
	"""
	j = matrix_backtrack.shape[1] - 1 
	i = matrix_backtrack.shape[0] - 1
	align_res = []
	align_pdb = []
	while i > 0 and j > 0:
		if matrix_backtrack[i][j] == "d":
			align_pdb.append(sequence_pdb_res.seq[j-1].res)
			align_res.append(sequence_target.seq[i-1].res)
			i -= 1
			j -= 1			
		elif matrix_backtrack[i][j] == "u":
			align_pdb.append("-")
			align_res.append(sequence_target.seq[i-1].res)
			i -= 1
		elif matrix_backtrack[i][j] == "l":			
			align_pdb.append(sequence_pdb_res.seq[j-1].res)
			align_res.append("-")
			j -= 1
	return align_pdb, align_res


__authors__ = "Ragousandirane RADJASANDIRANE"
__contact__ = "radja.ragou@gmail.com"
__date__ = "14/09/2021"
__version__= "1.0"

# Scripts import
import matrix_ddp


# Modules import
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, mahalanobis
import multiprocessing as mp
from dict_sequence import *
from joblib import Parallel, delayed
import blosum
import copy
import os


def first_align(sequence_target, sequence_pdb_res):
	n_col = sequence_pdb_res.length + 1
	n_row = sequence_target.length + 1
	matrix_align = np.zeros((n_row,n_col), dtype=np.float64)
	matrix_backtrace = np.empty((n_row,n_col),dtype=str)
	matrix_backtrace[0,1:] = "l"
	matrix_backtrace[1:,0] = "u"
	matrix_backtrace[0,0] = "x"
	gap = 0
	for i in range(1, n_row):
		for j in range(1, n_col):
			left = matrix_align[i][j-1] + gap
			up = matrix_align[i-1][j] + gap
			diag = matrix_align[i-1][j-1] + blosum.blosum62(sequence_target.seq[i-1].res,sequence_pdb_res.seq[j-1].res)
			matrix_align[i][j] =  max(diag, left, up)  
			if matrix_align[i][j] == diag:
				matrix_backtrace[i][j] = "d"
			elif matrix_align[i][j] == left:
				matrix_backtrace[i][j] = "l"
			else:
				matrix_backtrace[i][j] = "u"
	first_align_pdb, first_align_res = matrix_ddp.get_alignment(sequence_target, sequence_pdb_res, matrix_backtrace)
	return get_perc_id(first_align_pdb, first_align_res)

def get_perc_id(first_align_pdb, first_align_res):
	match = 0
	total = 0
	for i in range(len(first_align_pdb)):
		if first_align_pdb[i] != "-" or first_align_res[i] != '-':
			if first_align_pdb[i] == first_align_res[i]:
				match += 1
			total += 1
	return int(round(match / total,2) * 100)

def get_weigth(perc_id):
	perc_id = float(perc_id/100)
	if perc_id <= 0.3:
		return (0.9,0.1)
	elif perc_id <= 0.6:
		return (0.7,0.3)
	elif perc_id <= 0.9:
		return (0.5,0.5)
	return (0.3,0.7)


def get_ss(pdb_file):
	dssp_filename = "dssp.txt"
	os.system("dssp {} -o {}".format(pdb_file,dssp_filename))
	structure = []

	with open(dssp_filename, 'r') as stride_file:
		read = False
		for line in stride_file:
			if read:          
				SS = line[16]
				if SS == "H" or SS == "E":
					structure.append(SS)
				else:
					structure.append("C")
			else:
				if line.strip().startswith("#"):
					read = True
	os.remove(dssp_filename)

	return structure

def ss_to_gap(secondary_structure):
	gap_score = {"C" : 0, "H" : -2, "E" : -2}
	return [gap_score[ss] for ss in secondary_structure]


def binary_search(col, beg, end, dist):
    if beg > end:
        return False
    m_index = (beg + end) // 2
    mid = col[m_index]
    values = mid.split("_")
    first = float(values[0])
    second = float(values[1])
    if first <= dist <= second:
        return mid
    if dist < first:
        return binary_search(col, beg,m_index-1,dist)
    return binary_search(col, m_index+1,end,dist)


def get_score(sequence_target,sequence_pdb_ca,dope_score, matrix_dist ,pairs_residues):

	target_res1 = sequence_target.seq3[pairs_residues[0][0] - 1].res
	target_res2 = sequence_target.seq3[pairs_residues[1][0] - 1].res
	calpha1 = sequence_pdb_ca.seq[pairs_residues[0][1] - 1]
	calpha2 = sequence_pdb_ca.seq[pairs_residues[1][1] - 1]
	distance_observed = matrix_dist[str(calpha1)][str(calpha2)]

	all_col =  dope_score.columns[4:]
	if distance_observed > 15.25:
		return 0
	if distance_observed < 0.25:
		return -10
	colname = binary_search(all_col,0,all_col.size - 1,distance_observed)
	index_res = "_".join([target_res1,target_res2])
	return  float(dope_score[colname][index_res])
		

def distribution(sequence_target, sequence_pdb_res,sequence_pdb_ca, matrix_dist, dope_score,options):
	shuffled_sequence = copy.deepcopy(sequence_target)
	value = []
	print("\n\nStarting Zcore computation")
	for step in range(options.zscore):
		print("\nRun {} / {}".format(step+1,options.zscore))	
		print("Sequence shuffling ...")
		shuffled_sequence.shuffle_seq()

		HM, _  = matrix_ddp.Hmatrix(shuffled_sequence, sequence_pdb_res,sequence_pdb_ca, matrix_dist, dope_score,options)
		value.append(HM[-1][-1])
	return value

def compute_zscore(value, obs):
	print(value,obs,np.mean(value),np.std(value))
	return ((obs - np.mean(value)) / np.std(value))
	

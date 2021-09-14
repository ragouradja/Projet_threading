"All functions needed to deal with scoring"

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
	"""Compute first alignment to get initial percent identity
	
	It is compute using simple dynamic programming and a Blosum
	only score (without Dope) if the option "--weight -1" is activated

    Parameters
    ----------
    sequence_target : Sequence_Residue Object
        Target sequence to align with template

	sequence_pdb_res : Sequence_Residue Object
		Sequence of residue from PDB file (template)

    Returns
    -------
    Call
        Call get_perf_id() to determinate the percent identity from this alignment
    """
	# Dimensions of Matrix with extra row and column for gap (Needleman & Wunsch)
	n_col = sequence_pdb_res.length + 1
	n_row = sequence_target.length + 1
	matrix_align = np.zeros((n_row,n_col), dtype=np.float64)
	matrix_backtrack = np.empty((n_row,n_col), dtype=str)
	matrix_backtrack[0,1:] = "l"
	matrix_backtrack[1:,0] = "u"
	matrix_backtrack[0,0] = "x"
	gap = 0
	for i in range(1, n_row):
		for j in range(1, n_col):
			left = matrix_align[i][j-1] + gap
			up = matrix_align[i-1][j] + gap
			diag = matrix_align[i-1][j-1] + blosum.blosum62(sequence_target.seq[i-1].res, sequence_pdb_res.seq[j-1].res)
			matrix_align[i][j] =  max(diag, left, up)  
			if matrix_align[i][j] == diag:
				matrix_backtrack[i][j] = "d"
			elif matrix_align[i][j] == left:
				matrix_backtrack[i][j] = "l"
			else:
				matrix_backtrack[i][j] = "u"
	first_align_pdb, first_align_res = matrix_ddp.get_alignment(sequence_target, sequence_pdb_res, matrix_backtrack)
	return get_perc_id(first_align_pdb, first_align_res)

def get_perc_id(first_align_pdb, first_align_res):
	"""Determinate the initial percent identity from first alignment
	
	The initial percent identity will give the weights to use for Dope
	and Blosum score.

    Parameters
    ----------
    first_align_pdb : List
        PDB sequence resulting from first alignment

	first_align_res : List
		Target sequence resulting from first alignment

    Returns
    -------
    int
        Inital percent identity
    """
	match = 0
	total = 0
	for i in range(len(first_align_pdb)):
		if first_align_pdb[i] != "-" or first_align_res[i] != '-':
			if first_align_pdb[i] == first_align_res[i]:
				match += 1
			total += 1
	return int(round(match / total,2) * 100)

def get_weigth(perc_id):
	"""Determinate weights to use to balance Dope and Blosum scores

	Blosum score is favored if the inital percent identity between
	sequences is high, otherwise, Dope score is favored for very
	low inital percent identity.

    Parameters
    ----------
    perc_id : int
        Inital percent identity

    Returns
    -------
    Tuple
        Weight to use for (Dope,Blosum)
    """
	perc_id = float(perc_id/100)
	if perc_id <= 0.3:
		return (0.9,0.1)
	elif perc_id <= 0.6:
		return (0.7,0.3)
	elif perc_id <= 0.9:
		return (0.5,0.5)
	return (0.3,0.7)


def get_ss(pdb_file):
	"""Using DSSP to assign secondary structure of template

	The secondary structure can be used to adjust gap score. This score
	is more penalizing if gap occurs in alpha helix or strand.

    Parameters
    ----------
    pdb_file : String
        Name of the PDB file (template)

    Returns
    -------
    structure : List
        Secondary structure for each residue in a list
    """
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
	"""Determinate score gap from secondary structure

    Parameters
    ----------
    secondary_structure : list
        Secondary structure for each residue in a list

    Returns
    -------
    list
        List of gap score for each residue of template
    """
	gap_score = {"C" : 0, "H" : -2, "E" : -2}
	return [gap_score[ss] for ss in secondary_structure]


def binary_search(col, beg, end, dist):
	"""Binary search to find the lower and upper bound where observed distance fit in

	Lower and upper bound of distance are the column name of the Dope score matrix
	It is needed to get a Dope score from a pair of residue and a distance
	The colname of Dope matrix are formated as : lowerbound_upperbound

	Parameters
	----------
	col : List
		Colname of Dope Matrix
	beg : int
		first index of binary search
	end : int
		last index of binary search
	dist : float
		Observed distance between two residue

	Returns
	-------
	Recursive call
		Return the name of the column to get the Dope score according to
		a pair of residue and a observed distance
	"""
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
	"""Extract Dope score between two target residue using 3D information from template 

	Observed distance between residue of a pair of residue is extracted from distance matrix 
	Binary search determine in which column we need to get the dope score according to
	the observed distance

	Parameters
	----------
    sequence_target : Sequence_Residue Object
        Target sequence to align with template

    sequence_pdb_ca : Sequence Object
        Sequence of alpha carbon from PDB file (template)

	dope_score : DataFrame
		All dope score in a dataframe
	matrix_dist : Numpy matrix
		All distances computed between alpha carbon of template

	pairs_residues : Tuple
		Pair of target residue with a alpha carbon from template
			((i_res1,j_carbon1), (i_res2,j_carbon2))

	Returns
	-------
	float
		Dope score of interaction between pair of target residue
	"""

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
		

def distribution(sequence_target, sequence_pdb_res,sequence_pdb_ca,
 matrix_dist, dope_score,options):
	"""Compute shuffling runs for zscore

	Each shuffling run give a score from the best score of high matrix

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
	value : list
		List of all best score from High matrix with random target sequence
	"""
	
	shuffled_sequence = copy.deepcopy(sequence_target)
	value = []
	print("\n\nStarting Zcore computation")
	for step in range(options.zscore):
		print("\nRun {} / {}".format(step+1,options.zscore))	
		print("Sequence shuffling ...")
		shuffled_sequence.shuffle_seq()

		HM, _  = matrix_ddp.Hmatrix(shuffled_sequence, sequence_pdb_res,
		sequence_pdb_ca, matrix_dist, dope_score,options)
		value.append(HM[-1][-1])
	return value

def compute_zscore(value, obs):
	"""Compute zscore using random score and observed score from high matrix

	Parameters
	----------
	value : list
		 List of score of random target sequence

	obs : float
		Observed best score from high matrix used for zscore

	Returns
	-------
	float
		Zscore
	"""
	return ((obs - np.mean(value)) / np.std(value))
	

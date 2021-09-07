import itertools
from Bio import SeqIO
import numpy as np
from numpy.core.records import array
from numpy.lib.shape_base import _make_along_axis_idx
import pandas as pd
from scipy.spatial.distance import cdist
import time
import random
import copy
from multiprocessing import Process,Pool
from joblib import Parallel, delayed
import sys

"""
Installation CONDA et autre dependances
Tester iterative vs parallel : pb arrondi --> aligment different !
Marche pas avec WINDOWS !
ne pas oublier dope_filter.py
"""

def pdb_to_seq(file_pdb):
	with open(file_pdb, "r") as pdb:
		for record in SeqIO.parse(pdb, 'pdb-atom'):
			return record.seq


def pdb_to_df(file_pdb):
	content = []
	header_pdb = ["ATOM","NumATOM","NameATOM","Residue",
				"Chain","NumResidue", "x","y","z","Occ","Temp","Atom_Element"]
	with open(file_pdb,"r") as pdb:
		for line in pdb:
			if line.startswith("ATOM"):
				items = line.split()
				if items[2] == "CA":
					content.append(items)
	df = pd.DataFrame(content, columns= header_pdb)
	return df

def sequence_one_to_three(sequence):
	amino_acid = {'C' : 'CYS', 'D' : 'ASP', 'S' : 'SER', 'Q' : 'GLN', 'K' : 'LYS',
     'I' : 'ILE', 'P' : 'PRO', 'T' : 'THR', 'F' : 'PHE', 'N' : 'ASN', 
     'G' : 'GLY', 'H' : 'HIS', 'L' : 'LEU', 'R' : 'ARG', 'W' : 'TRP', 
     'A' : 'ALA', 'V' : 'VAL', 'E' : 'GLU', 'Y' : 'TYR', 'M' : 'MET'}

	sequence_three_letter = [amino_acid[sequence[i]]
	for i in range(len(sequence))]

	return sequence_three_letter

def sequence_three_to_one(sequence):
	amino_acid = {'CYS' :'C', 'ASP' :'D', 'SER' :'S', 'GLN' :'Q', 'LYS' :'K',
     'ILE' :'I', 'PRO' :'P', 'THR' :'T', 'PHE' :'F', 'ASN' :'N', 
     'GLY' :'G', 'HIS' :'H', 'LEU' :'L', 'ARG' :'R', 'TRP' :'W', 
     'ALA' :'A', 'VAL' :'V', 'GLU' :'E', 'TYR' :'Y', 'MET' :'M'}

	sequence_one_letter = [amino_acid[sequence[i]]
	for i in range(len(sequence))]

	return sequence_one_letter



def fasta_to_seq(file_fasta):
	amino_acid = {'C' : 'CYS', 'D' : 'ASP', 'S' : 'SER', 'Q' : 'GLN', 'K' : 'LYS',
     'I' : 'ILE', 'P' : 'PRO', 'T' : 'THR', 'F' : 'PHE', 'N' : 'ASN', 
     'G' : 'GLY', 'H' : 'HIS', 'L' : 'LEU', 'R' : 'ARG', 'W' : 'TRP', 
     'A' : 'ALA', 'V' : 'VAL', 'E' : 'GLU', 'Y' : 'TYR', 'M' : 'MET'}
	sequence_one_letter = []
	with open(file_fasta, "r") as fasta:
		for line in fasta:
			if not line.startswith(">"):
				sequence_one_letter += line.strip()
	sequence_one_letter = list(sequence_one_letter)
	sequence_three_letter = sequence_one_to_three(sequence_one_letter)
	return sequence_one_letter,sequence_three_letter


def get_coords(df):
	coords = df[["x","y","z"]].to_numpy()
	sequence_pdb_CA = df["NumResidue"].to_numpy()
	sequence_pdb_res = sequence_three_to_one(df["Residue"].to_numpy())
	return coords,sequence_pdb_CA,sequence_pdb_res


def matrix_distance(coords_pdb, sequence_pdb_ca):
	matrix = pd.DataFrame(cdist(coords_pdb, coords_pdb,metric="euclidean"))
	matrix.columns = sequence_pdb_ca
	matrix.index = sequence_pdb_ca
	return matrix

def get_score(dope_score, matrix_dist, pairs_residues):
	target_res1 = pairs_residues[0][:3]
	target_res2 = pairs_residues[1][:3]
	calpha1 = pairs_residues[0][3:]	
	calpha2 = pairs_residues[1][3:]
	distance_observed = matrix_dist[str(calpha1)][str(calpha2)]
	for dist_col in dope_score.columns[4:]:
		distance_values = dist_col.split("_")
		first_distance = float(distance_values[0])
		second_distance = float(distance_values[1])
		if distance_observed > first_distance and distance_observed < second_distance:
			return  float(dope_score[(dope_score["res1"] == target_res1) 
			& (dope_score["res2"] == target_res2)][dist_col])
		if distance_observed > 15.25:
			return 0
		if distance_observed < 0.25:
			return -10

def Hmatrix(sequence_target, matrix_dist, dope_score):
	n_col = matrix_dist.shape[1] + 1
	n_row = len(sequence_target) + 1
	matrix_high = np.zeros([n_row,n_col], dtype=np.float16)
	matrix_backtrace = np.empty((n_row,n_col),dtype=np.str)
	matrix_backtrace[0,1:] = "l"
	matrix_backtrace[1:,0] = "u"
	matrix_backtrace[0,0] = "x"
	beg_start = time.time()

	#matrix_check = np.zeros((n_row-1,n_col-1))
	matrix_check = do_LM_parallel(sequence_target, matrix_dist, dope_score, n_row, n_col)
	print("PARALLEL DONE IN : ",time.time() - beg_start)
	# FOR col : 0.8 - 0.9s --> iterative
	# FOR col : 0.7 - 0.8s --> rec
	for i in range(1, n_row):
		for j in range(1, n_col):
			left = matrix_high[i][j-1]
			up = matrix_high[i][j-1]
			residue_fixed = (i,j)
			best_score = matrix_check[i-1][j-1]
			#best_score = Lmatrix(sequence_target, matrix_dist, dope_score, residue_fixed, n_row, n_col)
			#matrix_check[i-1][j-1] = best_score
			diag = matrix_high[i-1][j-1] + best_score 
			matrix_high[i][j] =  max(diag, left, up)
			if matrix_high[i][j] == diag:
				matrix_backtrace[i][j] = "d"
			elif matrix_high[i][j] == left:
				matrix_backtrace[i][j] = "l"
			else:
				matrix_backtrace[i][j] = "u"
			
	return matrix_high, matrix_check,matrix_backtrace

def do_LM_parallel(sequence_target,matrix_dist, dope_score,n_row,n_col):
	coords = list(itertools.product(range(1,n_row), range(1,n_col)))
	matrix_check = Parallel(n_jobs=6, verbose = 0)(delayed(Lmatrix)(sequence_target,
	 matrix_dist, dope_score,[i,j],n_row,n_col) for i,j in coords)
	matrix_check = np.array(matrix_check)
	matrix_check = np.reshape(matrix_check, (n_row-1,n_col-1))
	return matrix_check
   
def Lmatrix(sequence_target, matrix_dist, dope_score, residue_fixed, n_row, n_col):
	print("DOING : ", residue_fixed, n_row,n_col)
	list_calpha = matrix_dist.columns 
	matrix_low = np.full([n_row,n_col], -100,dtype=np.float16)

	matrix_low[0][:] = np.zeros((1,n_col))
	matrix_low[:,0] = np.zeros((1,n_row))

	residue_fixed_i = residue_fixed[0]
	residue_fixed_j = residue_fixed[1]
	name_res_fixed = sequence_target[residue_fixed_i-1] + list_calpha[residue_fixed_j-1]  # SER2
	for i in range(1,n_row):
		for j in range(1,n_col):
			if i == residue_fixed_i and j == residue_fixed_j:
				matrix_low[residue_fixed_i][residue_fixed_j] = matrix_low[i-1][j-1]
			if (j < residue_fixed_j and i < residue_fixed_i) or (j > residue_fixed_j and i > residue_fixed_i):
				actual_residue = sequence_target[i-1] + list_calpha[j-1]
				pairs_residues = [name_res_fixed,actual_residue]
				score = get_score(dope_score,matrix_dist, pairs_residues)
				left = matrix_low[i][j-1]
				up = matrix_low[i-1][j]
				diag = matrix_low[i-1][j-1] 
				matrix_low[i][j] = score + max(diag,left,up)
	return matrix_low[-1][-1]


def get_alignement(sequence_target,list_pdb_res,matrix_backtrace):
	j = matrix_backtrace.shape[1] - 1 
	i = matrix_backtrace.shape[0] - 1
	align_res = []
	align_pdb = []
	while i > 0 and j > 0:
		if matrix_backtrace[i][j] == "d":
			align_pdb.append(list_pdb_res[j-1])
			align_res.append(sequence_target[i-1])
			i -= 1
			j -= 1			
		elif matrix_backtrace[i][j] == "u":
			align_pdb.append("-")
			align_res.append(sequence_target[i-1])
			i -= 1
		elif matrix_backtrace[i][j] == "l":			
			align_pdb.append(list_pdb_res[j-1])
			align_res.append("-")
			j -= 1
	return align_pdb, align_res


def print_alignement(align_pdb, align_res,informations):
	seq = informations[0]
	pdb = informations[1]
	fasta = informations[2]
	last_value_HM = informations[3]
	shape = informations[4]
	time = informations[5]
	align = open("align.txt","a")
	align.write('# Sequence (Length) : {} ({})\n'.format("".join(seq),len(seq)))
	align.write("# Fasta file : {}\n".format(fasta))
	align.write("# PDB file : {}\n".format(pdb))
	align.write("# Last value of High Matrix : {:.3f}\n".format(last_value_HM))
	align.write("# Shape of High Matrix : {}\n".format(shape))
	align.write("# Total time : {:.3f}s\n\n".format(time))
	# A optimiser !
	for i in range(len(align_pdb)-1,-1,-1):
		align.write("{}".format(align_pdb[i]))
		print("{}".format(align_pdb[i]), end = "")
	align.write("\n")
	print()
	for i in range(len(align_res)-1,-1,-1):
		align.write("{}".format(align_res[i]))
		print("{}".format(align_res[i]), end = "")
	align.write("\n\n\n")
	align.close()
	print()


def main(file_pdb, file_fasta):

	start = time.time()


	# Fasta
	sequence_one_letter,sequence_three_letter = fasta_to_seq(file_fasta)
	print(sequence_one_letter,sequence_three_letter)
	print("FASTA : ",time.time() - start)
	# PDB
	df = pdb_to_df(file_pdb)
	coords_pdb,sequence_pdb_ca, sequence_pdb_res = get_coords(df)
	matrix_dist = matrix_distance(coords_pdb, sequence_pdb_ca)
	matrix_check = np.zeros((20,20))
	# Dope
	file_dope = "../data/dope_clean.txt"
	dope_score = pd.read_table(file_dope)
	# Matrix
	matrix_high, matrix_check,matrix_backtrace = Hmatrix(sequence_three_letter, matrix_dist, dope_score)

	informatiions_to_write = [sequence_one_letter,file_pdb,file_fasta,
	matrix_high[-1][-1],matrix_high.shape, time.time()-start]
	print(matrix_high)
	align_pdb, align_res = get_alignement(sequence_one_letter, sequence_pdb_res, matrix_backtrace)
	print_alignement(align_pdb, align_res,informatiions_to_write)
	print(time.time()-start, matrix_high.shape,matrix_high[-1][-1])



if __name__ == "__main__":
	file_pdb = sys.argv[1]
	file_fasta = sys.argv[2]

	main(file_pdb,file_fasta)


	# (78,78) --> 91minutes : 5482.439617872238s
	# Total time : 3.514s --> WSL
	# Total time : 3.831s --> W10 shell
	# Total time : 4.195s --> VSC
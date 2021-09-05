from Bio import SeqIO
from Bio.SeqUtils import seq3
import numpy as np
from numpy.lib.shape_base import _make_along_axis_idx
import pandas as pd
from scipy.spatial.distance import cdist
import time
import random
import copy
from multiprocessing import Process


"""
Remplacer FOR par reccursive !
Debug chaque fonction pour choper le bon align
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


def fasta_to_seq(file_fasta):
	amino_acid = {'C' : 'CYS', 'D' : 'ASP', 'S' : 'SER', 'Q' : 'GLN', 'K' : 'LYS',
     'I' : 'ILE', 'P' : 'PRO', 'T' : 'THR', 'F' : 'PHE', 'N' : 'ASN', 
     'G' : 'GLY', 'H' : 'HIS', 'L' : 'LEU', 'R' : 'ARG', 'W' : 'TRP', 
     'A' : 'ALA', 'V' : 'VAL', 'E' : 'GLU', 'Y' : 'TYR', 'M' : 'MET'}
	
	with open(file_fasta, "r") as fasta:
		for record in SeqIO.parse(fasta,"fasta"):
			sequence_one_letter = list(record.seq)
			if sequence_one_letter[0] == "M":
				sequence_one_letter = sequence_one_letter[1:]

	sequence_three_letter = [amino_acid[sequence_one_letter[i]]
	for i in range(len(sequence_one_letter))]
	return sequence_one_letter,sequence_three_letter


def get_coords(df):
	return df[["x","y","z"]].to_numpy(),df["NumResidue"].to_numpy()


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
	beg_start = time.time()
	# FOR col : 0.8 - 0.9s --> iterative
	# FOR col : 0.7 - 0.8s --> rec
	for i in range(1, n_row):
		print(i, matrix_high.shape)
		start = time.time()
		for j in range(1, n_col):
			print(f"{j}")
			left = matrix_high[i][j-1]
			up = matrix_high[i][j-1]
			residue_fixed = (i,j)
			best_score = Lmatrix(sequence_target, matrix_dist, dope_score, residue_fixed, n_row, n_col)
			diag = matrix_high[i-1][j-1] 
			matrix_high[i][j] =  best_score + max(diag, left, up)
		end = time.time() - start
		print("TIME FOR : ", end)
		global_end = time.time() - beg_start
		print("TIME REMAINING : ", end*(n_row-1) - global_end)
	return matrix_high

   
def Lmatrix(sequence_target, matrix_dist, dope_score, residue_fixed, n_row, n_col):
	list_calpha = matrix_dist.columns 
	matrix_low = np.zeros([n_row,n_col], dtype=np.float16)
	residue_fixed_i = residue_fixed[0]
	residue_fixed_j = residue_fixed[1]
	name_res_fixed = sequence_target[residue_fixed_i-1] + list_calpha[residue_fixed_j-1]  # SER2
	for i in range(1,n_row):
		for j in range(1,n_col):
			if i == residue_fixed_i and j == residue_fixed_j:
				matrix_low[residue_fixed_i][residue_fixed_j] = matrix_low[i-1][j-1]
			elif (j < residue_fixed_j and i < residue_fixed_i) or (j > residue_fixed_j and i > residue_fixed_i):
				actual_residue = sequence_target[i-1] + list_calpha[j-1]
				pairs_residues = [name_res_fixed,actual_residue]
				score = get_score(dope_score,matrix_dist, pairs_residues)
				left = matrix_low[i][j-1]
				up = matrix_low[i-1][j]
				diag = matrix_low[i-1][j-1] 
				matrix_low[i][j] = score + max(diag,left,up)
	return matrix_low[-1][-1]

def sum_low(matrix_low):
	i = matrix_low.shape[0] - 1
	j = matrix_low.shape[1] - 1 
	som = 0
	while i > 0 and j > 0:
		diag = matrix_low[i-1][j-1]
		left = matrix_low[i][j-1]
		up = matrix_low[i-1][j]
		if diag >= left and diag >= up:
			som += diag
			i -= 1
			j -= 1
		elif left > diag and left > diag:
			som += left

			j -= 1
		elif up > diag and up > left:
			som += up
			i -= 1
	return som


def distribution(sequence_target, matrix_dist, dope_score):
	step = 10
	shufflled_sequence = copy.deepcopy(sequence_target)
	value = []
	for i in range(step):
		random.shuffle(shufflled_sequence)
		print(shufflled_sequence)
		HM = Hmatrix(shufflled_sequence, matrix_dist, dope_score)
		value.append(HM)
		print(HM)
	return value


def zscore(value, obs):
	print((obs - np.mean(value)) / np.std(value))

def get_alignement(sequence_target,list_calpha,matrix_high):
	i = matrix_high.shape[0] - 1
	j = matrix_high.shape[1] - 1 
	align_res = []
	align_ca = []
	while i > 0 and j > 0:
		diag = matrix_high[i-1][j-1]
		left = matrix_high[i][j-1]
		up = matrix_high[i-1][j]
		if (diag >= left and diag >= up) or (up > diag and up > left):
			align_ca.append(list_calpha[j-1])
			align_res.append(sequence_target[i-1])
			i -= 1
		else:
			align_ca.append(list_calpha[j-1])
			align_res.append("-")
		j -= 1
	return align_ca, align_res


def print_alignement(align_ca, align_res,informations):
	seq = informations[0]
	pdb = informations[1]
	last_value_HM = informations[2]
	shape = informations[3]
	time = informations[4]
	align = open("align.txt","a")
	align.write('# Sequence (Length) : {} ({})\n'.format("".join(seq),len(seq)))
	align.write("# PDB file : {}\n".format(pdb))
	align.write("# Last value of High Matrix : {:.3f}\n".format(last_value_HM))
	align.write("# Shape of High Matrix : {}\n".format(shape))
	align.write("# Total time : {:.3f}s\n\n".format(time))
	for i in range(len(align_ca)-1,-1,-1):
		align.write("{:^5s}".format(align_ca[i]))
		print("{:^5s}".format(align_ca[i]), end = "")
	align.write("\n")
	print()
	for i in range(len(align_res)-1,-1,-1):
		align.write("{:^5s}".format(align_res[i]))
		print("{:^5s}".format(align_res[i]), end = "")
	align.write("\n\n\n")
	align.close()
	print()


def main(file_pdb, file_fasta):

	start = time.time()
	# Fasta
	sequence_one_letter,sequence_three_letter = fasta_to_seq(file_fasta)
	print(time.time() - start)
	# PDB
	df = pdb_to_df(file_pdb)
	coords_pdb,sequence_pdb_ca = get_coords(df)
	matrix_dist = matrix_distance(coords_pdb, sequence_pdb_ca)

	# Dope
	file_dope = "../data/dope_clean.txt"
	dope_score = pd.read_table(file_dope)

	# Matrix
	matrix_high = Hmatrix(sequence_three_letter, matrix_dist, dope_score)

	#value = distribution(sequence_three_letter, matrix_dist, dope_score)
	#print(value)
	#zscore(value, matrix_high)

	informatiions_to_write = [sequence_one_letter,file_pdb,
	matrix_high[-1][-1],matrix_high.shape, time.time()-start]

	align_ca, align_res = get_alignement(sequence_one_letter, sequence_pdb_ca, matrix_high)
	print_alignement(align_ca, align_res,informatiions_to_write)

	print(matrix_high)
	print(time.time()-start, matrix_high.shape,matrix_high[-1][-1])

 
if __name__ == "__main__":
	file_fasta = "../data/1dfn.fasta"
	file_pdb = "../data/1bnb.pdb"

	main(file_pdb, file_fasta)


	# (78,78) --> 91minutes : 5482.439617872238s
	# Total time : 3.514s --> WSL
	# Total time : 3.831s --> W10 shell
	# Total time : 4.195s --> VSC
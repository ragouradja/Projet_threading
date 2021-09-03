from Bio import SeqIO
from Bio.SeqUtils import seq3
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import time

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
	with open(file_fasta, "r") as fasta:
		for record in SeqIO.parse(fasta,"fasta"):
			fasta_sequence_one_letter = record.seq
			if fasta_sequence_one_letter[0] == "M":
				fasta_sequence_one_letter = fasta_sequence_one_letter[1:]
			fasta_sequence_three_letter = seq3(fasta_sequence_one_letter)

	final_sequence_three_letter = [fasta_sequence_three_letter.upper()[i:i+3]
	for i in range(0,len(fasta_sequence_three_letter),3)]
	return fasta_sequence_one_letter,final_sequence_three_letter


def get_coords(df):
	return df[["x","y","z"]].to_numpy(),df["NumResidue"].to_numpy()


def matrix_distance(coords_pdb, sequence_pdb):
	matrix = pd.DataFrame(cdist(coords_pdb, coords_pdb,metric="euclidean"))
	matrix.columns = sequence_pdb
	matrix.index = sequence_pdb
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
		if distance_observed > 15.25:
			return 0
		if distance_observed > first_distance and distance_observed < second_distance:
			return  float(dope_score[(dope_score["res1"] == target_res1) 
			& (dope_score["res2"] == target_res2)][dist_col])

def Hmatrix(sequence_target, matrix_dist, dope_score):
	n_col = len(matrix_dist.columns)
	n_row = len(sequence_target)
	matrix_high = np.zeros([n_row,n_col])
	fixed = np.full((1,n_col),np.nan)
	for i in range(1, n_row):
		for j in range(1, n_col):
			left = matrix_high[i][j-1]
			up = matrix_high[i][j-1]
			residue_fixed = [i,j]
			best_score = Lmatrix(sequence_target, matrix_dist, dope_score, residue_fixed)
			diag = matrix_high[i-1][j-1] + best_score
			matrix_high[i][j] = max(diag, left, up)

	return matrix_high



def Lmatrix(sequence_target, matrix_dist, dope_score, residue_fixed):
	list_calpha = matrix_dist.columns 
	n_col = len(list_calpha) + 1
	n_row = len(sequence_target) + 1
	matrix_low = np.zeros([n_row,n_col])
	fixed = np.full((1,n_col-1),np.nan)
	residue_fixed_i = residue_fixed[0]
	residue_fixed_j = residue_fixed[1]
	#matrix_low[residue_fixed_i][1:] = fixed
	#matrix_low[1:,residue_fixed_j] = fixed
	target_res_fixed = sequence_target[residue_fixed[0] -1] + list_calpha[residue_fixed_j-1]  # SER2

	#print(target_res_fixed)
	#print(matrix_low)
	#print(matrix_dist)
	#print(dope_score)
	
	for i in range(1,n_row):
		for j in range(1,n_col):
			if i == residue_fixed_i and j == residue_fixed_j:
				matrix_low[residue_fixed_i][residue_fixed_j] = matrix_low[i-1][j-1]
			elif matrix_low[i][j] != np.nan:
				if (j < residue_fixed_j and i < residue_fixed_i) or (j > residue_fixed_j and i > residue_fixed_i):
					actual_residue = sequence_target[i-1] + list_calpha[j-1]
					pairs_residues = [target_res_fixed,actual_residue]
					score = get_score(dope_score,matrix_dist, pairs_residues)
					left = matrix_low[i][j-1]
					up = matrix_low[i-1][j]
					diag = matrix_low[i-1][j-1] + score
					# max entre valeur de GAUCHE + gap, HAUT + gap ou DIAG + score DOPE
					matrix_low[i][j] = max(diag,left,up)
	return matrix_low[n_row-1][n_col-1]

def get_alignement(sequence_target,list_calpha,matrix_high):
	j = matrix_high.shape[1] - 1 
	i = matrix_high.shape[0] - 1
	align_res = []
	align_ca = []
	while i > 0 and j > 0:
		diag = matrix_high[i-1][j-1]
		left = matrix_high[i][j-1]
		up = matrix_high[i-1][j]
		if diag >= left and diag >= up:
			align_ca.append(list_calpha[j-1])
			align_res.append(sequence_target[i-1])
			i -= 1
			j -= 1
		elif left > diag and left > diag:
			align_ca.append(list_calpha[j-1])
			align_res.append("-")
			j -= 1
		elif up > diag and up > left:
			align_ca.append("-")
			align_res.append(sequence_target[i-1])
			i -= 1
	return align_ca, align_res


def print_alignement(align_ca, align_res):
	for i in range(len(align_ca)-1,-1,-1):
		print("{:^3s}".format(align_ca[i]), end = "")

	print()
	for i in range(len(align_res)-1,-1,-1):
		print("{:^3s}".format(align_res[i]), end = "")
	print()

if __name__ == "__main__":
	start = time.time()
	file_pdb = "../data/prot.pdb"
	file_fasta = "../data/2d0a.fasta"
	file_dope = "../data/dope_clean.txt"
	n_res = 20
	sequence_one_letter,sequence_three_letter = fasta_to_seq(file_fasta)
	sequence_target = sequence_three_letter[30:30+n_res]
	print(sequence_three_letter[30:30+n_res])
	dope_score = pd.read_table(file_dope)
	df = pdb_to_df(file_pdb)
	coords_pdb,sequence_pdb = get_coords(df.head(50))
	matrix_dist = matrix_distance(coords_pdb,sequence_pdb)
	#Hmatrix(sequence_target, matrix_dist, dope_score)
	Hmatrice = Hmatrix(sequence_target, matrix_dist, dope_score)
	print(Hmatrice)
	align_ca, align_res = get_alignement(sequence_one_letter[30:30+n_res],matrix_dist.columns,Hmatrice)
	print_alignement(align_ca, align_res)
	print(time.time() - start) #1.35s #160s pour 20x50
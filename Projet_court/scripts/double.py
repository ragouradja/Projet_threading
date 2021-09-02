from Bio import SeqIO
from Bio.SeqUtils import seq3
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import time

"""
Remplacer FOR par reccursive !

A FAIRE :

Remplir la matrice H !
matrice L : OK

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
			fasta_sequence_three_letter = seq3(fasta_sequence_one_letter)

	final_sequence = [fasta_sequence_three_letter.upper()[i:i+3]
	for i in range(0,len(fasta_sequence_three_letter),3)]

	return final_sequence


def get_coords(df):
	return df[["x","y","z"]].to_numpy(),df["NumResidue"].to_numpy()


def matrix_distance(coords_pdb, sequence_pdb):
	matrix = pd.DataFrame(cdist(coords_pdb, coords_pdb,metric="euclidean"))
	matrix.columns = sequence_pdb
	matrix.index = sequence_pdb
	return matrix

def get_score(dope_score, matrix_dist, pairs_residues = ["SER2","VAL5"]):
	target_res1 = pairs_residues[0][:3]
	target_res2 = pairs_residues[1][:3]
	calpha1 = pairs_residues[0][-1]	
	calpha2 = pairs_residues[1][-1]	
	distance_observed = matrix_dist[calpha1][calpha2]
	for dist_col in dope_score.columns[4:]:
		distance_values = dist_col.split("_")
		first_distance = float(distance_values[0])
		second_distance = float(distance_values[1])
		if distance_observed > first_distance and distance_observed < second_distance:
			return float(dope_score[(dope_score["res1"] == target_res1) 
			& (dope_score["res2"] == target_res2)][dist_col])

def Hmatrix(sequence_target, matrix_dist, dope_score):
	n_col = len(matrix_dist.columns)
	n_row = len(sequence_target)
	matrix_high = np.zeros([n_row,n_col])
	fixed = np.full((1,n_col),np.nan)
	print(matrix_high)



def Lmatrix(sequence_target, matrix_dist, dope_score, residue_fixed = [4,6]):
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
					score = get_score(dope_score,matrix_dist)
					left = matrix_low[i][j-1]
					up = matrix_low[i-1][j]
					diag = matrix_low[i-1][j-1] + score
					# max entre valeur de GAUCHE + gap, HAUT + gap ou DIAG + score DOPE
					matrix_low[i][j] = max(diag,left,up)
	print(matrix_low)

if __name__ == "__main__":

	file_pdb = "prot.pdb"
	file_fasta = "2d0a.fasta"
	sequence_target = fasta_to_seq("2d0a.fasta")[:9]
	print(sequence_target)
	dope_score = pd.read_table("dope_clean.txt")
	df = pdb_to_df(file_pdb)
	coords_pdb,sequence_pdb = get_coords(df.head(9))
	matrix_dist = matrix_distance(coords_pdb,sequence_pdb)
	#Hmatrix(sequence_target, matrix_dist, dope_score)
	get_score(dope_score,matrix_dist)
	Lmatrix(sequence_target, matrix_dist, dope_score)


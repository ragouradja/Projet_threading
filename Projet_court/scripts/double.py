import itertools
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import time
import multiprocessing as mp
from joblib import Parallel, delayed
import sys
import blosum


"""
Installation CONDA et autre dependances
Tester iterative vs parallel : pb arrondi --> aligment different !
Marche pas avec WINDOWS !
ne pas oublier dope_filter.py
"""
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
	matrix = pd.DataFrame(np.around(cdist(coords_pdb, coords_pdb,metric="euclidean"),decimals=2))
	matrix.columns = sequence_pdb_ca
	matrix.index = sequence_pdb_ca
	return matrix

def create_dict():
	array = ['0.25_0.75', '0.75_1.25', '1.25_1.75', '1.75_2.25', '2.25_2.75', '2.75_3.25', '3.25_3.75', '3.75_4.25',
	'4.25_4.75', '4.75_5.25', '5.25_5.75', '5.75_6.25', '6.25_6.75', '6.75_7.25', '7.25_7.75', '7.75_8.25', '8.25_8.75',
	'8.75_9.25', '9.25_9.75', '9.75_10.25', '10.25_10.75', '10.75_11.25', '11.25_11.75', '11.75_12.25', '12.25_12.75',
	'12.75_13.25', '13.25_13.75', '13.75_14.25', '14.25_14.75', '14.75_15.25']

	value = []
	for i in range(len(array)):
		items = array[i].split("_")
		new_list = [(n,array[i]) for n in np.round(np.arange(float(items[0]),float(items[1]),0.01),2)]
		value += new_list

	all = dict(value)
	return all


def binary_search(col, beg, end, dist):
    if beg > end:
        return False
    m_index = (beg + end) // 2
    mid = col[m_index]
    values = mid.split("_")
    first = float(values[0])
    second = float(values[1])
    if first <= dist < second:
        return mid
    if dist < first:
        return binary_search(col, beg,m_index-1,dist)
    return binary_search(col, m_index+1,end,dist)


def get_score(dope_score, matrix_dist, pairs_residues,dope_dict):
	target_res1 = pairs_residues[0][:3]
	target_res2 = pairs_residues[1][:3]
	calpha1 = pairs_residues[0][3:]	
	calpha2 = pairs_residues[1][3:]
	distance_observed = matrix_dist[str(calpha1)][str(calpha2)]

	all_col =  dope_score.columns[4:]

	if distance_observed > 15.25:
		return 0
	col_distance = binary_search(all_col,0,len(all_col) - 1,distance_observed)
	# col_distance = dope_dict[distance_observed]
	if col_distance:
		return  float(dope_score[(dope_score["res1"] == target_res1) 
			& (dope_score["res2"] == target_res2)][col_distance])
	return -10
		


def get_score0(dope_score, matrix_dist, pairs_residues):
	target_res1 = pairs_residues[0][:3]
	target_res2 = pairs_residues[1][:3]
	calpha1 = pairs_residues[0][3:]	
	calpha2 = pairs_residues[1][3:]
	distance_observed = matrix_dist[str(calpha1)][str(calpha2)]
	all_col = set(dope_score.columns[4:].tolist())
	for dist_col in all_col:
		distance_values = dist_col.split("_")
		first_distance = float(distance_values[0])
		second_distance = float(distance_values[1])
		if distance_observed >= first_distance and distance_observed <= second_distance:

			return  float(dope_score[(dope_score["res1"] == target_res1) 
			& (dope_score["res2"] == target_res2)][dist_col])
		if distance_observed > 15.25:

			return 0
		if distance_observed < 0.25:
			return -10


def Hmatrix(sequence_target,list_pdb_res, matrix_dist, dope_score,nproc):

	seq_one = sequence_three_to_one(sequence_target)
	dico_dope = create_dict()
	n_col = matrix_dist.shape[1] + 1
	n_row = len(sequence_target) + 1
	matrix_high = np.zeros([n_row,n_col], dtype=np.float16)
	matrix_backtrace = np.empty((n_row,n_col),dtype=np.str)
	matrix_backtrace[0,1:] = "l"
	matrix_backtrace[1:,0] = "u"
	matrix_backtrace[0,0] = "x"
	beg_start = time.time()
	gap = 0
	#matrix_check = np.zeros((n_row-1,n_col-1))
	matrix_check = do_LM_parallel(sequence_target,seq_one,list_pdb_res, matrix_dist, dope_score, n_row, n_col,nproc,dico_dope)
	# FOR col : 0.8 - 0.9s --> iterative
	# FOR col : 0.7 - 0.8s --> rec
	for i in range(1, n_row):
		for j in range(1, n_col):
			left = matrix_high[i][j-1] + gap
			up = matrix_high[i-1][j] + gap
			best_score = matrix_check[i-1][j-1]
			#best_score = Lmatrix(sequence_target, matrix_dist, dope_score, residue_fixed, n_row, n_col)
			#matrix_check[i-1][j-1] = best_score
			diag = matrix_high[i-1][j-1]  + best_score + blosum.blosum62(seq_one[i-1],list_pdb_res[j-1])
			matrix_high[i][j] =  max(diag, left, up)
			if matrix_high[i][j] == diag:
				matrix_backtrace[i][j] = "d"
			elif matrix_high[i][j] == left:
				matrix_backtrace[i][j] = "l"
			else:
				matrix_backtrace[i][j] = "u"
	return matrix_high, matrix_check,matrix_backtrace

def do_LM_parallel(sequence_target,seq_one ,list_pdb_res,matrix_dist, dope_score,n_row,n_col,nproc,dico_dope):
	coords = list(itertools.product(range(1,n_row), range(1,n_col)))
	matrix_check = Parallel(n_jobs=nproc - 1, verbose = 0, prefer="processes")(delayed(Lmatrix)
	(sequence_target,seq_one ,list_pdb_res,matrix_dist, dope_score,[fixed_i,fixed_j],
	n_row,n_col,dico_dope) for fixed_i,fixed_j in coords)
	matrix_check = np.array(matrix_check, dtype=np.float16)
	matrix_check = np.reshape(matrix_check, (n_row-1,n_col-1))
	return matrix_check
   
def Lmatrix(sequence_target, seq_one, list_pdb_res, matrix_dist, dope_score, residue_fixed, n_row, n_col,dico_dope):
	list_calpha = matrix_dist.columns 
	matrix_low = np.full([n_row,n_col], -100,dtype=np.float16)
	start = time.time()
	matrix_low[0][:] = np.zeros((1,n_col))
	matrix_low[:,0] = np.zeros((1,n_row))
	gap = 0
	residue_fixed_i = residue_fixed[0]
	residue_fixed_j = residue_fixed[1]
	name_res_fixed = sequence_target[residue_fixed_i-1] + list_calpha[residue_fixed_j-1]  # SER2
	for i in range(1,n_row):
		for j in range(1,n_col):
			if (i,j) == (residue_fixed_i,residue_fixed_j):
				matrix_low[residue_fixed_i][residue_fixed_j] = matrix_low[i-1][j-1]
			if (j < residue_fixed_j and i < residue_fixed_i) or (j > residue_fixed_j and i > residue_fixed_i):
				actual_residue = sequence_target[i-1] + list_calpha[j-1]
				pairs_residues = [name_res_fixed,actual_residue]
				score = get_score(dope_score,matrix_dist, pairs_residues,dico_dope)
				left = matrix_low[i][j-1] + gap
				up = matrix_low[i-1][j] + gap
				diag = matrix_low[i-1][j-1] 
				matrix_low[i][j] = score + max(diag,left,up) + blosum.blosum62(seq_one[i-1],list_pdb_res[j-1])

	if residue_fixed_j == n_col-1:
		time_x_line = time.time() - GLOBAL_START_TIME
		x = residue_fixed_i
		y = n_row - residue_fixed_i
		time_remaining = y * time_x_line / x
		print("Time remaining (estimated): {:.2f}s".format(time_remaining))
				
	if matrix_low[-1][-1] == -100:
		return 0
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
	nproc = informations[6]
	align = open("align.txt","a")
	align.write('# Sequence (Length) : {} ({})\n'.format("".join(seq),len(seq)))
	align.write("# Fasta file : {}\n".format(fasta))
	align.write("# PDB file : {}\n".format(pdb))
	align.write("# Last value of High Matrix : {:.3f}\n".format(last_value_HM))
	align.write("# Shape of High Matrix : {}\n".format(shape))
	align.write("# Total time : {:.3f}s\n".format(time))
	align.write("# Nb of processes used : {}\n\n".format(nproc))

	align_pdb.reverse()
	align_res.reverse()

	align.writelines(align_pdb)
	align.write("\n")	
	align.writelines(align_res)
	align.write("\n\n")	

	print(*align_pdb)
	print(*align_res)


def main(file_pdb, file_fasta,nproc):

	start = time.time()


	# Fasta
	fasta_one_letter,fasta_three_letter = fasta_to_seq(file_fasta)
	# PDB
	pdb_content = pdb_to_df(file_pdb)

	coords_res_pdb,list_pdb_ca, list_pdb_res = get_coords(pdb_content)

	matrix_dist_pdb = matrix_distance(coords_res_pdb, list_pdb_ca)
	# Dope
	file_dope = "../data/dope_clean.txt"
	matrix_dope = pd.read_table(file_dope)

	# Matrix
	matrix_high, matrix_check_parallel,matrix_backtrace = Hmatrix(fasta_three_letter, list_pdb_res,matrix_dist_pdb, matrix_dope,nproc)
	informatiions_to_write = [fasta_one_letter,file_pdb,file_fasta,
	matrix_high[-1][-1],matrix_high.shape, time.time()-start,nproc]
	align_pdb, align_res = get_alignement(fasta_one_letter, list_pdb_res, matrix_backtrace)
	print_alignement(align_pdb, align_res,informatiions_to_write)
	print(time.time()-start, matrix_high.shape,matrix_high[-1][-1])


GLOBAL_START_TIME = time.time()
if __name__ == "__main__":
	file_pdb = sys.argv[1]
	file_fasta = sys.argv[2]
	nproc = mp.cpu_count() - 1

	main(file_pdb,file_fasta,nproc)

	# (78,78) --> 91minutes : 5482.439617872238s
	# Total time : 3.514s --> WSL
	# Total time : 3.831s --> W10 shell
	# Total time : 4.195s --> VSC
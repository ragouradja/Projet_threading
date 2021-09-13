"""Main script for compute Threading using Double Dynamics Programming

Usage:
======
	$ python main.py -p file.pdb -f file.fasta 
	
	file.pdb: path to the PDB file (template)
	file.fasta: path to the Fasta file (target for threading)
"""

__authors__ = "Ragousandirane RADJASANDIRANE"
__contact__ = "radja.ragou@gmail.com"
__date__ = "14/09/2021"
__version__= "1.0"

# Scripts import
from matrix_ddp import *
from classes import *
from visual_stuff import *
from dict_sequence import *

# Modules import
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, mahalanobis
import time
import multiprocessing as mp
from joblib import Parallel, delayed



"""
INSTALLER DSSP !!
Enlever les prints !

1.6.2
Installation CONDA et autre dependances
Tester iterative vs parallel : pb arrondi --> aligment different !
Marche pas avec WINDOWS !
ne pas oublier dope_filter.py
"""
GLOBAL_START_TIME = time.time()

def main():
	"""Main function that run all tasks
	
	- Save the user's arguments
	- Extracting informations from the PDB file : 
		- Sequence of residue
		- Sequence of alpha carbon's positions
		- 3D coordinates of alpha carbon
	- Calculating matrix distance from 3D Coordinates of alpha carbon
	- Extracting target sequence from fasta file
	- Load Dope score matrix
	- Fill the High matrix
	- Write all results and options in a log file
	"""

	start = time.time()
	options = user_param()
	# Sequences of template from PDB file
	pdb = PDB(options.pdb, options.chain)
	seq_pdb_res = pdb.get_seq()
	seq_pdb_CA = pdb.get_CA()
	coords_pdb = pdb.get_coords()
	# Pair distance between residue of PDB
	matrix_dist_pdb = matrix_distance(coords_pdb, seq_pdb_CA)
	# Target sequence from Fasta file
	seq_fasta = Fasta(options.fasta).get_seq()
	start = time.time()
	score_align.first_align(seq_fasta,seq_pdb_res)
	# Loading Dope matrix
	file_dope = "../data/dope_clean.txt"
	matrix_dope = pd.read_table(file_dope, index_col=0)
	# Fill the High Matrix
	matrix_high, matrix_backtrace = Hmatrix(seq_fasta,seq_pdb_res,seq_pdb_CA,
	matrix_dist_pdb, matrix_dope, options)
	# Writing alignement results in log file
	informatiions_to_write = []
	if options.zscore > 0:
		list_value = score_align.distribution(seq_fasta, seq_pdb_res,seq_pdb_CA,
		matrix_dist_pdb, matrix_dope, options)
		zscore = round(score_align.compute_zscore(list_value,
		matrix_high[-1][-1]),2)
		informatiions_to_write += [zscore]
	informatiions_to_write = [seq_fasta,options.pdb,options.fasta,
	matrix_high[-1][-1],matrix_high.shape, time.time()-start,
	options.cpu] + informatiions_to_write
	align_pdb, align_res = get_alignment(seq_fasta, seq_pdb_res, matrix_backtrace)
	print_alignment(align_pdb, align_res,informatiions_to_write,options)
	print(time.time()-start, matrix_high.shape,matrix_high[-1][-1])

if __name__ == "__main__":
	main()
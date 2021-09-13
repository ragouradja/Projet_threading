"""Functions for save user's arguements and write results in log file"""

__authors__ = "Ragousandirane RADJASANDIRANE"
__contact__ = "radja.ragou@gmail.com"
__date__ = "14/09/2021"
__version__= "1.0"


# Scripts import
from classes import *
from dict_sequence import *
import score_align

# Modules import
import multiprocessing as mp
import sys
import os
import argparse
import textwrap






def user_param():
	"""Get the user's argument using argparse"""
		
	parser = argparse.ArgumentParser(description="Compute alignment between a target sequence and a template using DOPE score and Double dynamic programming",
	formatter_class=argparse.RawTextHelpFormatter)
	
	parser.add_argument("-p","--pdb", help="Path to the PDB file (template).", type=str, required=True, metavar="")
	parser.add_argument("-f","--fasta", help="Path to the Fasta file (target).", type=str, required=True, metavar="")

	parser.add_argument("--cpu", help=textwrap.dedent('''\
		Number of CPUs to use (Default -1 : All CPUs but one are used).
		If the number entered exceed your PC capacity, max CPU will be use.'''),type=int, metavar="", default=-1)

	parser.add_argument("--blosum", help="Add Blosum62 score to DOPE score for alignment.", action = "store_true")
	parser.add_argument("--weight", help=textwrap.dedent('''\
		Define weights to use for DOPE and Blosum score. 
		Default 0 : No weights is used (DOPE and Blosum are used equally).
		-1 : Weights are computed according to the initial percent identity
		between target and template sequences using blosum only
		1 : Weights = (0.9,0.1) --> if option is -1, these weights are chosen if 
					     the sequence identity is <= 30
		2 : Weights = (0.7,0.3) --> if sequence identity is <= 60
		3 : Weights = (0.5,0.5) --> if sequence identity is <= 90
		4 : Weights = (0.3,0.7) --> if sequence identity is > 90
		(Option used with --blosum)'''),  type = int, default = 0, metavar="")

	parser.add_argument("--zscore", help=textwrap.dedent('''\
		Number of iteration (run) to compute zscore.
		Default 0 : Simple run without computing zscore.
		Execution time will be multiplicated by this number of run.'''),default=0, type=int, metavar="")
	parser.add_argument("--dssp", help="Add gap penalty based on the secondary structure of the template.", action = "store_true")
	parser.add_argument("--log", help="Name of the log file. Default `results.log`.",default="results.log", type = str,metavar="")
	parser.add_argument("--chain", help="Chain to select in PDB file.",default=None, type = str,metavar="")
		

	options = parser.parse_args()
	pdb = options.pdb
	fasta = options.fasta
	cpu = options.cpu
	blosum = options.blosum
	weight = options.weight
	zscore = options.zscore
	dssp = options.dssp
	log = options.log
	chain = options.chain

	if not os.path.isfile(pdb):
		sys.exit(f"The file {pdb} does not exist. Please enter a valid PDB file.")

	if not os.path.isfile(fasta):
		sys.exit(f"The file {fasta} does not exist. Please enter a valid Fasta file.")

	if cpu != -1 and cpu <= 0:
		sys.exit(f"Can't use {cpu} CPU. Please enter a valid number (>=1 or -1) of CPU to use.")
	else:
		max_cpu = mp.cpu_count()
		if cpu == -1:
			n_cpu = mp.cpu_count() - 1
			options.cpu = n_cpu
		elif cpu >= 1:			
			if cpu >= max_cpu:
				options.cpu = max_cpu
	print(f"{options.cpu} CPUs will be use.")
	
	if blosum:
		print("Blosum62 will be use for computing score addtionally to DOPE score.")
		if not weight:
			print("No weight is used for Dope and Blosum score since --weight = 0.")

	if weight != 0:
		if not blosum:
			print("--wieght is to be use with --blosum. This will have no effect on results.")
		elif weight == -1:
			print("Weights for Dope and Blosum are going to be calculated according to percent identity between sequences.")

	if zscore:
		print(f"Zscore will be computed on {zscore} runs.")
	options.zrun = False

	if dssp:
		print("Using secondary structure of the template to compute gap penalty.")
	if chain:
		print(f"Chain selected : {chain}")

	options.log = f"../log/{options.log}"
	print(f"All results are going to be written in {options.log}")
	print()
	return options

def print_alignment(align_pdb, align_res,informations,options):
	"""Print and write results and options used in a log file"""
	seq = informations[0]
	pdb = informations[1]
	fasta = informations[2]
	last_value_HM = informations[3]
	shape = informations[4]
	time = informations[5]
	nproc = informations[6]
	perc_id = score_align.get_perc_id(align_pdb, align_res)
	logfile_path = options.log
	align = open(logfile_path,"a")
	align.write('########### Options used ###########\n')
	align.write("# PDB file : {}\n".format(pdb))
	align.write("# Fasta file : {}\n".format(fasta))
	align.write("# Nb of processes used : {}\n".format(nproc))

	if options.blosum:
		align.write('# Blosum matrix\n')
		if options.weight != 0:
			align.write('# Weights used : {} Dope {} Blosum\n'.format(options.weight_used[0],options.weight_used[1]))
	
	if options.zscore != 0:
		align.write("# Zscore ({} runs): {}\n".format(options.zscore,informations[-1]))

	if options.dssp:
		align.write('# Secondary structure used for gap penalty\n')

	if options.chain:
		align.write('# Chain selected for template : {}\n'.format(options.chain))
		
	align.write('######## Other informations ########\n')
	align.write('# Sequence (Length) : {} ({})\n'.format(seq(),seq.length))
	align.write("# Final percent identity : {}%\n".format(perc_id))
	align.write("# Last value of High Matrix : {:.3f}\n".format(last_value_HM))
	align.write("# Shape of High Matrix : {}\n".format(shape))
	align.write("# Total time : {:.3f}s\n".format(time))
	align.write("####################################\n\n")

	align_pdb.reverse()
	align_res.reverse()
	align.writelines(align_pdb)
	align.write("\n")	
	align.writelines(align_res)
	align.write("\n\n\n\n")	
	print(*align_pdb)
	print(*align_res)

if __name__ == "__main__":
    user_param()
# Threading by double dynamics programming

This program is an implementation of the method of Threading using double dynamics programming.


Requirements
------------
All dependencies can be installed with the env/threading.yml file.

Optionally, DSSP can be used. 

Installation
------------
Clone this repository : 

    $ git clone git@github.com:ragouradja/Projet_threading.git
Or : 

    $ git clone https://github.com/ragouradja/Projet_threading.git
Then :

    $ cd Projet_threading

Load the conda environment : 
    

    $ conda env create -f env/threading.yml

Install DSPP : 

    $ sudo apt install dssp

Usage
-----

For a simple usage, give a PDB and Fasta file : 

    $ python3 main.py -p ../data/1ard_1.pdb -f ../data/1znf_1.fasta
    7 CPUs will be use.
    All results are going to be written in ../log/results.log

    Filling all low matrix in paralell...
    Time remaining (estimated): 6.05s
    Time remaining (estimated): 4.37s
    Time remaining (estimated): 1.23s
    Time remaining (estimated): 2.65s
    Time remaining (estimated): 1.02s
    Time remaining (estimated): 1.31s
    Time remaining (estimated): 1.11s
    Time remaining (estimated): 1.38s
    Time remaining (estimated): 1.29s
    Time remaining (estimated): 0.93s
    Time remaining (estimated): 1.08s
    Time remaining (estimated): 1.11s
    Time remaining (estimated): 0.99s
    Time remaining (estimated): 0.76s
    Time remaining (estimated): 1.23s
    Time remaining (estimated): 0.68s
    Time remaining (estimated): 0.51s
    Time remaining (estimated): 0.58s
    Time remaining (estimated): 0.83s
    Time remaining (estimated): 0.33s
    Time remaining (estimated): 0.43s
    Time remaining (estimated): 0.13s
    Time remaining (estimated): 0.22s
    Time remaining (estimated): 0.35s
    Time remaining (estimated): 0.07s
    R S F V C E V C T R A F A R Q E H L K R H Y R S H T N E K
    Y K - C - G L C E R S F V E K S A - L S R - H Q R V H K N

    Check the ../log/results.log file
    
Help
----
    $ python3 main.py -h
    usage: main.py [-h] -p pdb -f fasta [--cpu int] [--blosum] [--weight int] [--zscore int] [--dssp] [--log str]
               [--chain str]

    Compute alignment between a target sequence and a template using DOPE score and Double dynamic programming

    optional arguments:
      -h, --help    show this help message and exit
      -p pdb        Path to the PDB file (template).
      -f fasta      Path to the Fasta file (target).
      --cpu int     Number of CPUs to use (Default -1 : All CPUs but one are used).
                    If the number entered exceed your PC capacity, max CPU will be use.
      --blosum      Add Blosum62 score to DOPE score for alignment.
      --weight int  Define weights to use for DOPE and Blosum score.
                    Default 0 : No weights is used (DOPE and Blosum are used equally).
                    -1 : Weights are computed according to the initial percent identity
                    between target and template sequences using blosum only
                    1 : Weights = (0.9,0.1) --> if option is -1, these weights are chosen if
                                                 the sequence identity is <= 30
                    2 : Weights = (0.7,0.3) --> if sequence identity is <= 60
                    3 : Weights = (0.5,0.5) --> if sequence identity is <= 90
                    4 : Weights = (0.3,0.7) --> if sequence identity is > 90
                    (Option used with --blosum)
      --zscore int  Number of iteration (run) to compute zscore.
                    Default 0 : Simple run without computing zscore.
                    Execution time will be multiplicated by this number of run.
      --dssp        Add gap penalty based on the secondary structure of the template.
      --log str     Name of the log file. Default `results.log`.
      --chain str   Chain to select in PDB file.

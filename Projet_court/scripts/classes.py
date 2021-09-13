__authors__ = "Ragousandirane RADJASANDIRANE"
__contact__ = "radja.ragou@gmail.com"
__date__ = "14/09/2021"
__version__= "1.0"

# Scripts import
from dict_sequence import *

# Modules import
import random
import pandas as pd
import time
import sys


class Residue:
    """Class Residue to represent each residue in a sequence"""
    def __init__(self, residue):
        """Transforming a residue """
        self.res = residue
        if isinstance(residue, Residue):
            self.state = len(residue.res)
        else:
            self.state = len(residue)
        
    def to_three(self):
        if self.state == 1:
            return Residue(amino_acid_one_three[self.res])
        return self

    def to_one(self):
        if self.state == 3:
            return Residue(amino_acid_three_one[self.res])
        return self

    def __call__(self):
        print(self.res)

    def __str__(self):
        return self.res

    def __repr__(self):
        return self.res

class Sequence:
    """Class Sequence"""
    def __init__(self,sequence,filename):
        self.seq = []
        if not isinstance(self,Sequence_Residue):
            self.seq = sequence
        self.length = len(sequence)
        self.filename = filename

    def reverse_sequence(self):
        self.seq.reverse()
        
    def shuffle_seq(self):
        random.shuffle(self.seq)

    def __call__(self):
        return "".join([str(x) for x in self.seq])

class Sequence_Residue(Sequence):

    def __init__(self, sequence, filename): # sequence : list
        Sequence.__init__(self,sequence,filename)
        self.seq3 = []

        for amino_acid in sequence:
            res = Residue(amino_acid)
            self.seq.append(res.to_one())
            self.seq3.append(res.to_three())

    def shuffle_seq(self):        
        super().shuffle_seq()
        self.seq3 = []
        for res in self.seq:
            self.seq3.append(res.to_three())

    def __call__(self, state = 1):
        if state == 1:
            return super().__call__()
        else:
            return "".join([str(x) for x in self.seq3])
            

class PDB:
    def __init__(self,file_pdb,chain = None):
        self.filename = file_pdb
        self.chain = chain
        self.pdb_to_df()
    
    def pdb_to_df(self):
        content = []
        header_pdb = ["ATOM","NumATOM","NameATOM","Residue",
                    "Chain","NumResidue", "x","y","z","Occ","Temp","Atom_Element"]
        with open(self.filename,"r") as pdb:
            for line in pdb:
                if line.startswith("ATOM"):
                    items = line.split()
                    if items[2] == "CA":
                        if self.chain:
                            if items[4] == self.chain:
                                content.append(items)
                        else:
                            content.append(items)
        if len(content) == 0:
            sys.exit("Please, enter a valid file in PDB format or a valid chain")
        self.df = pd.DataFrame(content, columns = header_pdb)
    
    def get_df(self):
        return self.df

    def get_coords(self):
        return self.df[["x","y","z"]].to_numpy()

    def get_seq(self):
        seq_res = self.df["Residue"].to_numpy()
        return Sequence_Residue(seq_res, self.filename)

    def get_CA(self):
        seq_CA = self.df["NumResidue"].to_numpy()
        return Sequence(seq_CA, self.filename)

class Fasta:
    def __init__(self, file_fasta):
        self.filename = file_fasta
        self.fasta_to_seq()

    def fasta_to_seq(self):
        self.sequence = []
        with open(self.filename, "r") as fasta:
            for line in fasta:
                if not line.startswith(">"):
                    for aa in line.strip():
                        if aa not in amino_acid_one_three:
                            print("Deleting from fasta sequence : ", aa)
                            line = line.replace(aa, "")
                        else:
                            self.sequence += aa
        self.sequence = list(self.sequence)

    def get_seq(self):
        return Sequence_Residue(self.sequence, self.filename)        

if __name__ =="__main__":
    file = "../data/immunoglobulin/6fab.pdb"
    pdb = PDB(file,None).get_seq().seq
    print(pdb)
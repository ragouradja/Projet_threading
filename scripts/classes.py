"""Definition of object"""

__authors__ = "Ragousandirane RADJASANDIRANE"
__contact__ = "radja.ragou@gmail.com"
__date__ = "14/09/2021"
__version__= "1.0"

# Scripts import
from dict_sequence import *

# Modules import
import random
import pandas as pd
import sys


class Residue:
    """Class Residue to represent each residue in a sequence"""
    def __init__(self, residue):
        """Constructor of Residue
        
        Parameters
        ----------
        residue : string
            String of one residue
        """

        # Save the residue as a string
        self.res = residue

        # If residue is already an Residue Object
        # Determine length of the attribute that contains the residue
        if isinstance(residue, Residue):
            self.state = len(residue.res)
        else:
            self.state = len(residue)
        
    def to_three(self):
        """Converting one letter amino acid to three letter amino acid

        The conversion is enable only if the actual residue is in one letter state
        """
        if self.state == 1:
            return Residue(amino_acid_one_three[self.res])
        return self

    def to_one(self):
        """Converting three letter amino acid to one letter amino acid

        The conversion is enable only if the actual residue is in three letter state
        """
        if self.state == 3:
            return Residue(amino_acid_three_one[self.res])
        return self

    def __call__(self):
        """Redefining call method to print the residue"""
        print(self.res)

    def __str__(self):
        """Redefining str method to return the residue using the name of the object"""
        return self.res

    def __repr__(self):
        """Redefining repr method to modify the behavior of object when it is in a list"""
        return self.res

class Sequence:
    """Class Sequence to represent all sequences"""

    def __init__(self,sequence,filename):
        """Constructor of Sequence

        The attribute self.seq is either for one letter sequence of amino acid or sequence of alpha carbon
        from PDB        

        Parameters
        ----------
        sequence : list
            List of residue
        filename : string
            Name of the file where comes the sequence
        """
        self.seq = []
        if not isinstance(self,Sequence_Residue):
            self.seq = sequence
        self.length = len(sequence)
        self.filename = filename

    def reverse_sequence(self):
        """Reverse one letter sequence"""
        self.seq.reverse()
        
    def shuffle_seq(self):
        """Shuffling one letter sequence"""
        random.shuffle(self.seq)

    def __call__(self):
        """Redefining call method to return a string of sequence"""
        return "".join([str(x) for x in self.seq])

class Sequence_Residue(Sequence):
    """Class Sequence_Residue inherited from Sequence class to represent all sequence of amino acids"""

    def __init__(self, sequence, filename):
        """Constructor of Sequence_Residue

        The attribute self.seq3 is for three letter sequence of amino acid

        Parameters
        ----------
        sequence : list
            List of residue
        filename : string
            Name of the file where comes the sequence
        """
        Sequence.__init__(self,sequence,filename)
        self.seq3 = []

        # Filling attributes for one and three letter amino acids
        for amino_acid in sequence:
            res = Residue(amino_acid)
            self.seq.append(res.to_one())
            self.seq3.append(res.to_three())

    def shuffle_seq(self):        
        """Update three letter sequence after shuffling one letter sequence by calling superclass method"""
        super().shuffle_seq()
        self.seq3 = []
        for res in self.seq:
            self.seq3.append(res.to_three())

    def __call__(self, state = 1):
        """Redefining call method to return a string of sequence
        
        Parameters
        ----------
        state : int
            state of the sequence we want to return to string
            If the state wanted is 1, supercall method will be called to deal with 
            one letter sequence (self.seq)
            Otherwise, three letter sequence (self.seq3) is transformed in string
        """
        if state == 1:
            return super().__call__()
        else:
            return "".join([str(x) for x in self.seq3])
            

class PDB:
    """Class PDB for PDB file"""

    def __init__(self,file_pdb,chain = None):
        """Constructor of PDB
        
        Parameters
        ----------
        file_pdb : string
            file name of PDB
        chain : string
            Letter of chain to extract if there is one
            Otherwise, extract all residue
        """
        self.filename = file_pdb
        self.chain = chain
        self.pdb_to_df()
    
    def pdb_to_df(self):
        """Converting the PDB file to a Dataframe for extracting informations"""

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
        """Return the dataframe"""
        return self.df

    def get_coords(self):
        """Return the 3D coordinates of alpha carbon"""
        return self.df[["x","y","z"]].to_numpy()

    def get_seq(self):
        """Return the sequence of residue from PDB in a Sequence_Residue object"""
        seq_res = self.df["Residue"].to_numpy()
        return Sequence_Residue(seq_res, self.filename)

    def get_CA(self):
        """Return the sequence of alpha carbon from PDB in a Sequence object"""
        seq_CA = self.df["NumResidue"].to_numpy()
        return Sequence(seq_CA, self.filename)

class Fasta:
    """Class Fasta for Fasta file"""
    def __init__(self, file_fasta):
        """Constructor of Fasta"""
        self.filename = file_fasta
        self.fasta_to_seq()

    def fasta_to_seq(self):
        """Extracting the sequence of residue from fasta file"""
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
        """Return the sequence of residue in a Sequence_Residue object"""
        return Sequence_Residue(self.sequence, self.filename)        

if __name__ =="__main__":
    file = "../data/immunoglobulin/6fab.pdb"
    pdb = PDB(file,None).get_seq().seq
    print(pdb)
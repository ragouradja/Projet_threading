from typing import ContextManager
from Bio import SeqIO
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import time

def pdbtoseq(file):
	with open(file, "r") as pdb:
		for record in SeqIO.parse(pdb, 'pdb-atom'):
			return record.seq

def pdb2df(file):
	content = []
	with open(file,"r") as pdb:
		for line in pdb:
			if line.startswith("ATOM"):
				items = line.split()
				if items[2] == "CA":
					content.append(items)
	df = pd.DataFrame(content)
	return df

def get_coords(df):
	return df[[6,7,8]].to_numpy(),(df[3]+df[5]).to_numpy()

def matrix_distance(coords, colnames):
	matrix = pd.DataFrame(cdist(coords, coords,metric="euclidean"))
	matrix.columns = colnames
	matrix.index = colnames
	print(matrix)




start = time.time()

file = "prot.pdb"
df = pdb2df(file)
coords,colnames = get_coords(df)

matrix_distance(coords,colnames)

print(time.time() - start)
#print(sequence)

"""
a = [[1,1,1],[2,1,2],[3,3,2]]
b = [[1,1,1],[2,1,2],[3,3,2]]

start = time.time()
results = distance.cdist(a,b,metric="euclidean")
print(results)
print(time.time() - start)
"""
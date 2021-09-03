import pandas as pd
import numpy as np

# Data initiale
file = "dope.txt"

# NOMS DE COLONNES
# Creation des noms de colonnes --> distance à partir de 0.25 puis pas de 0.5
# 30 distances
distance = np.arange(0.25,0.25+0.5*31,0.5)
dist_name = ["{}_{}".format(distance[i],distance[i+1]) for i in range(len(distance)-1)]
col = ["res1","atom1","res2","atom2"]
colname = col + dist_name
print(colname)

# CHARGEMENT DES DONNÉES
dtf = pd.read_table(file, delimiter=" ",header = None)

# FILTRE CARBONE ALPHA
CA = dtf[(dtf[1] == "CA") & (dtf[3] == "CA")]
# MAJ des noms de colonnes
CA.columns = colname

# Selection des colonnes numériques
# On multiplie par -1 toutes les valeurs de scores
negative_dtf = -CA.drop(col,axis = 1)

# On Concat avec les colonnes non numériques
final_dtf = pd.concat([CA[col], negative_dtf], axis = 1)
print(final_dtf)

# Save le dope clean --> Que les CA avec les noms de colonnes
final_dtf.to_csv("dope_clean.txt", header = True, sep="\t", index=False)



"""
DOPE clean :
* Uniquement les CA
* Noms de colonnes avec les fourchettes de distances
* Valeurs négatives devenues positives
"""
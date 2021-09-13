import pandas as pd
import numpy as np

# Data initiale
file = "../data/dope.txt"

# NOMS DE COLONNES
# Creation des noms de colonnes --> distance à partir de 0.25 puis pas de 0.5
# 30 distances
distance = np.arange(0.25,0.25+0.5*31,0.5)
dist_name = ["{}_{}".format(distance[i],distance[i+1]) for i in range(len(distance)-1)]
col = ["res1","atom1","res2","atom2"]
colname = col + dist_name
# print(colname)

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
dtf_wo_ca = pd.concat([CA[col], negative_dtf], axis = 1)
residues = dtf_wo_ca[["res1","res2"]]
res1_list = residues["res1"].tolist()
res2_list = residues["res2"].tolist()

new_index_res = ["_".join([res1_list[i],res2_list[i]]) for i in range(len(res2_list))]
dtf_wo_ca.index = new_index_res
dtf_final = dtf_wo_ca.drop(["res1","atom1","res2","atom2"], axis = 1)
print(dtf_final)
# Save le dope clean --> Que les CA avec les noms de colonnes
# dtf_final.to_csv("dope_clean.txt", header = True, sep="\t")



"""
DOPE clean :
* Uniquement les CA
* Noms de colonnes avec les fourchettes de distances
* Valeurs négatives devenues positives
* Noms des lignes : directement les noms des résidus
"""
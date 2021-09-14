"""Reformating the dope file"""

import pandas as pd
import numpy as np

# Inital data
file_dope = "../data/dope.par"

# COLNAME
# Reformating colname to : lowerbound_upperbound (0.25_0.75 for example)
# Step of 0.5 for each column
distance = np.arange(0.25,0.25+0.5*31,0.5)
dist_name = ["{}_{}".format(distance[i],distance[i+1]) for i in range(len(distance)-1)]
col = ["res1","atom1","res2","atom2"]
colname = col + dist_name

# LOAD DATA
dtf = pd.read_table(file_dope, delimiter=" ",header = None)

# Keep only alpha carbon
CA = dtf[(dtf[1] == "CA") & (dtf[3] == "CA")]
# Update colnames
CA.columns = colname

# Multiply by -1 all numerical column
negative_dtf = -CA.drop(col,axis = 1)

# Construct of the new dataframe
dtf_wo_ca = pd.concat([CA[col], negative_dtf], axis = 1)
residues = dtf_wo_ca[["res1","res2"]]
res1_list = residues["res1"].tolist()
res2_list = residues["res2"].tolist()

# Deleting first fourth column to put the residues names in index as Res1_Res2 (ALA_ALA for example)
new_index_res = ["_".join([res1_list[i],res2_list[i]]) for i in range(len(res2_list))]
dtf_wo_ca.index = new_index_res
dtf_final = dtf_wo_ca.drop(["res1","atom1","res2","atom2"], axis = 1)
print(dtf_final)
#dtf_final.to_csv("dope_clean.txt", header = True, sep="\t")

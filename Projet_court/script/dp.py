import numpy as np
from numpy.core.numerictypes import maximum_sctype

# Check le pb d'une sequnce plus grande que l'autre

seq2 = "ATGCT"
seq1 = "AGCT"

n1 = len(seq1)
n2 = len(seq2)
matrix = np.zeros([n2+1,n1+1])

for i in range(1,n2+1):
    for j in range(1,n1+1):
        left = matrix[i][j-1]
        up = matrix[i-1][j]
        if seq1[j-1] == seq2[i-1]:
            match = 1
        else:
            match = -1
        matrix[i][j] = max(left,up,matrix[i-1][j-1] + match)
print(matrix)

i = n2

align_seq1 = ""
align_seq2 = ""

for j in range(n1,0,-1):
    matrix[i][j] = -1
    print(matrix[i][j])
    if seq1[j-1] == seq2[i-1]:
        align_seq1 = seq1[j-1] + align_seq1
        align_seq2 = seq2[i-1] + align_seq2
        i -= 1
    else:
        if matrix[i][j-1] > matrix[i-1][j]:
            align_seq1 = seq1[j-1] + align_seq1
            align_seq2 = "-" + align_seq2
        else:
            align_seq1 = "-" + align_seq1       
            align_seq2 = seq2[i-1] + align_seq2
            i -= 1
            j += 1
print(matrix)

print(align_seq1)
print(align_seq2)

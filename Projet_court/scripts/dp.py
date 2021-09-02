import numpy as np
import time

# Check le pb d'une sequnce plus grande que l'autre
# TEST AVEC RECCURENCE PUIS time.time() !

seq2 = "ATGCT"
seq1 = "ACTG"

n1 = len(seq1)
n2 = len(seq2)

# Une colonne de plus pour les gaps
matrix = np.zeros([n2+1,n1+1])
start = time.time()

for i in range(1,n2+1):
    for j in range(1,n1+1):
        left = matrix[i][j-1]
        up = matrix[i-1][j]
        if seq1[j-1] == seq2[i-1]:
            match = 1
        else:
            match = -1
        # max entre valeur de GAUCHE + gap, HAUT + gap ou DIAG + match/missmatch
        matrix[i][j] = max(left,up,matrix[i-1][j-1] + match)
print(matrix)

print(time.time() - start)

align_seq1 = ""
align_seq2 = ""

# Coordonnées de la dernière case = best score
i = n2
j = n1

# On remonte le long de la matrice tant qu'on est pas sur les valeurs des gaps
# --> i ou j == 0
while(j > 0 and i > 0):
    # Match
    if seq1[j-1] == seq2[i-1]:
        align_seq1 = seq1[j-1] + align_seq1
        align_seq2 = seq2[i-1] + align_seq2
        # On recule en diagonale
        i -= 1

    # Missmatch
    else:
        # Si GAUCHE == MAX ; donc > HAUT
        # On va à GAUCHE
        if matrix[i][j-1] > matrix[i-1][j]:
            # gap pour la séquence 2
            align_seq1 = seq1[j-1] + align_seq1
            align_seq2 = "-" + align_seq2
        else:
            # Si HAUT == max
            # On va en HAUT
            # gap pour la séquence 1
            align_seq1 = "-" + align_seq1       
            align_seq2 = seq2[i-1] + align_seq2
            # Pour aller en haut --> i -= 1
            i -= 1
            # Pour contrer le j -= 1 car le j ne doit pas bouger dans ce cas (meilleur moyen ?)
            j += 1
    j -= 1
    
print(align_seq1)
print(align_seq2)
print(time.time() - start)

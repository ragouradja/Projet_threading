import itertools
import numpy as np
from multiprocessing import Pool
import time
from joblib import Parallel, delayed

def fill_array(start_val):
    return list(range(start_val, start_val+10))

def fill(a,x):
    print("DOING : ",x)
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            if i == x[0] and j == x[1]:
                a[i][j] = -100
            else :
                a[i][j] = 0

    print(a)
    return a
if __name__=='__main__':
    start = time.time()

    a = np.zeros((5,5))
    coords = list(itertools.product(range(1,5), repeat = 2))
    
    Parallel(n_jobs=1, verbose = 0)(delayed(fill)(a,[x,y]) for x,y in coords)
    print(time.time()-  start)
    
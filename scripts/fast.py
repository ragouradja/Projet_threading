from os import stat
from sys import getsizeof
import time
import numpy as np
from joblib import Parallel, delayed, parallel
from numba import jit
from dict_sequence import *


# def binary_search(array, beg, end,x):
#     if beg > end:
#         return False
#     m = (beg + end) // 2
#     mid = array[m]
#     values = mid.split("_")
#     first = float(values[0])
#     second = float(values[1])

#     if first <= x < second:
#         return mid
#     if x < first:
#         return binary_search(array, beg,m-1,x)
#     else:
#         return binary_search(array, m+1,end,x)


# def get_score(obs,dico):
#     return dico[obs]

# @jit(nopython=True)
# def score_DOPE(obs):
#     for name in array:
#         items = name.split("_")
#         val1 = float(items[0])
#         val2 = float(items[1])
#         if val1 <= obs < val2:
#             return (name)

# @jit(nopython=True)
# def search(obs,distances_array):
#     for line in range(len(distances_array)):
#         first = distances_array[line][0]
#         second = distances_array[line][1]
#         if first <= obs < second:
#             return first,second

# array = ['0.25_0.75', '0.75_1.25', '1.25_1.75', '1.75_2.25', '2.25_2.75', '2.75_3.25', '3.25_3.75', '3.75_4.25',
#  '4.25_4.75', '4.75_5.25', '5.25_5.75', '5.75_6.25', '6.25_6.75', '6.75_7.25', '7.25_7.75', '7.75_8.25', '8.25_8.75',
#   '8.75_9.25', '9.25_9.75', '9.75_10.25', '10.25_10.75', '10.75_11.25', '11.25_11.75', '11.75_12.25', '12.25_12.75',
#   '12.75_13.25', '13.25_13.75', '13.75_14.25', '14.25_14.75', '14.75_15.25']

# stealth_check = dict(
#                     [(n, 'You are about as stealthy as thunderstorm.')
#                         for n in range(1, 6)] +
#                     [(n, 'You tip-toe through the crowd of walkers, while loudly calling them names.')
#                         for n in range(6, 11)] +
#                     [(n, 'You are quiet, and deliberate, but still you smell.')
#                         for n in range(11, 16)] +
#                     [(n, 'You move like a ninja, but attracting a handful of walkers was inevitable.')
#                         for n in range(16, 20)]
#                     )


# # value = []
# # for i in range(len(array)):
# #     items = array[i].split("_")
# #     new_list = [(n,array[i]) for n in np.round(np.arange(float(items[0]),float(items[1]),0.01),2)]
# #     value += new_list

# # all_dico = dict(value)
# import pandas as pd
# import numpy as np
# start = time.time()
# for i in range(100000):
#     score_DOPE(i)
# print(time.time() - start)

# start = time.time()

# a = np.arange(100000)
# list(map(score_DOPE,a))
# print(time.time() - start)

# # for i in range(78*78*78*78):
# #     first,second = search(6,a)
# #     col = "_".join([str(first),str(second)])
# # print(time.time() - start)

# # dtf = pd.DataFrame(value, dtype = np.float16)

# # start = time.time()
# # for i in range(14):
# #     obs = 12
# #     inf = dtf[dtf <= obs]
# #     inf_index = inf.stack().index[-1][0]
# #     distances = dtf.values[inf_index]
# #     colname = "_".join([str(distances[0]),str(distances[1])])
# # print(time.time() - start)

seq = "ARKLKXXJK"
for aa in seq:
    if aa not in amino_acid_one_three:
        print(aa)
        seq = seq.replace(aa,"")
print(seq)
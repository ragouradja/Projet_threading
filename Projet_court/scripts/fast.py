import time
import numpy as np
def add(n):
    return n + 5

start = time.time()
a = np.zeros((1000,1000))

b = np.array(list(map(add,a[1]))).reshape((1000,1000))

print(b)
print(time.time() - start)
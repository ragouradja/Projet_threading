import time
start = time.time()
for i in range(1000):
    print(i*i)

print(time.time() - start)
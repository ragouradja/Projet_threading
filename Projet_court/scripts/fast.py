import time

def binary_search(array, beg, end,x):
    if beg > end:
        return False
    m = (beg + end) // 2
    mid = array[m]
    values = mid.split("_")
    first = float(values[0])
    second = float(values[1])

    if first <= x < second:
        return mid
    if x < first:
        return binary_search(array, beg,m-1,x)
    else:
        return binary_search(array, m+1,end,x)

array = ['0.25_0.75', '0.75_1.25', '1.25_1.75', '1.75_2.25', '2.25_2.75', '2.75_3.25', '3.25_3.75', '3.75_4.25',
 '4.25_4.75', '4.75_5.25', '5.25_5.75', '5.75_6.25', '6.25_6.75', '6.75_7.25', '7.25_7.75', '7.75_8.25', '8.25_8.75',
  '8.75_9.25', '9.25_9.75', '9.75_10.25', '10.25_10.75', '10.75_11.25', '11.25_11.75', '11.75_12.25', '12.25_12.75',
  '12.75_13.25', '13.25_13.75', '13.75_14.25', '14.25_14.75', '14.75_15.25']

start = time.time()




print(time.time()  - start)
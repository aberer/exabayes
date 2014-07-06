#! /usr/bin/python

import sys 
import re 

if len(sys.argv) != 2  : 
    print(sys.argv[0] + " file")
    sys.exit()



fh = open(sys.argv[1], "r")
fh.readline()

lines = fh.readlines()


lens = []
for line in lines: 
    num = line.strip().split()[0]
    lens.append(int(num))

thesum = sum(lens)
# print("total %d" % thesum )

# 99 percentile 
n99 = 0 
sofar = 0
goal = float(thesum) * 0.99
for elem in lens: 
    n99 += 1 
    sofar += elem
    if sofar >= goal : 
        break
# print("99-th perc: %d " % n99)

n50 = 0 
sofar = 0
goal = float(thesum) * 0.50
for elem in lens: 
    n50 += 1 
    sofar += elem
    if sofar >= goal : 
        break
# print("50-th perc: %d " % n50)

print("%d\t%d\t%d"  % (thesum , n50, n99))

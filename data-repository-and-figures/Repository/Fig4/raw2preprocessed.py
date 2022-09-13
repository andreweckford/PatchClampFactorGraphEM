# this module implements the decimation filter
# takes the raw data and averages over non-overlapping blocks of 50 samples

import numpy as np
import matplotlib.pyplot as plt
import csv
import sys

def readCSV(f,encoding=None):

    if encoding is None:
        with open(f) as csvfile:
            r = csv.reader(csvfile)
            m = []
            for i in r:
                m.append(list(i))
    else:
        with open(f,encoding=encoding) as csvfile:
            r = csv.reader(csvfile)
            m = []
            for i in r:
                m.append(list(i))
                
    return m

filename = sys.argv[1]
a = np.array(readCSV(filename,encoding='utf-8-sig'))

data_string = a[:,1]
data = np.zeros(len(data_string))
for i in range(0,len(data_string)):
  data[i] = float(data_string[i])

dec = []
current = 0
f = 50
while (current+1)*f < len(data):
  dec.append(np.sum(data[(current*f):((current+1)*f)])/f)
  current += 1
  
for i in range(0,len(dec)):
  if i < (len(dec)-1):
    print(str(dec[i])+',',end='')
  else:
    print(str(dec[i]))

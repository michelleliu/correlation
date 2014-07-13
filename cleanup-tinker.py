import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import math
import time

# usage: python cleanup.py inp_file

if len(sys.argv)<2:
    print "ERROR: not enough arguments"
inputfile = sys.argv[len(sys.argv)-1]
outputfile = inputfile + "-cleanup"

print "Reading {0}, saving to {1}".format(inputfile,outputfile)

data = open(inputfile)
out = open(outputfile, "w+")
for timestep in np.arange(200000):
    line = data.readline()
    split=line.split()
    for i in np.arange(len(split)-1):
        out.write("{0} ".format(float(split[i])))
    out.write("{0}\n".format(float(split[len(split)-1])))
    #line = data.readline()
    #out.write("{0}".format(line))
    #for atom in np.arange(5184):
    #    line = data.readline()
    #    split=line.split()
    #    for i in np.arange(2):
    #        out.write("{0} ".format(int(split[i])))
    #    for i in np.arange(5):
    #        out.write("{0} ".format(float(split[i+2])))
    #    out.write("{0}\n".format(float(split[7])))

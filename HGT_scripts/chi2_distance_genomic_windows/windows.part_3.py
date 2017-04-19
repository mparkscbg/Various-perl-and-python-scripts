
#This script takes output from the two companion scripts, windows.part_1.py and windows.part_2.py, and outputs chi2 values
#for each useable window in the genome output from windows.part_2.py. Output is in two, tab-delimited columns - the name
#of each genomic window, and its chi2 value against the whole genome.
#
#Syntax is:
#
#python windows.part_3.py fiG_values.txt fig_values.txt > chi2_values.txt
#
#NOTE 1: fiG, fig and chi2 values/calculations are as defined in Becq J, Churlaud C, Deschavanne P. 2010. A benchmark of parametric
#methods for horizontal transfers detection. PLoSONE 5(4):e9989.
#NOTE 2: this script is pretty fast, in part because the input (i.e., output from windows.part_2.py is only typically in the thousands
#to hundreds of thousands of lines, depending on genome size and window/step size used.



from __future__ import division

import re
import sys
import math
from Bio import SeqIO
import ast


fiG=sys.argv[1]
fig=sys.argv[2]


with open (fiG,'r') as f:
    fiG_var=f.readline().rstrip()       #read fiG values into fiG_var variable
    fiG_list=fiG_var.split(',')         #split fiG_var variable into list with comma delimiter
    fiG_list.pop(0)                     #pop off first element of fiG_list (first element is just 'fiG_list')

#print fiG_list


with open (fig,'r') as f:               #open fig file and process
    for line in f:                      #one line at a time
        fig_var=line.rstrip()
        fig_list=fig_var.split(',')
        if len(fig_list)>2:             #this line filters out windows with warning for too many masked positions, output from windows_2.py
            tetra_values=[]             #set empty tetra_values list, resets to empty for each line
            fig_name=fig_list[0]        #window name is first element
            fig_list.pop(0)             #pop off first element for calculations below
            for g,G in zip(fig_list,fiG_list):      #below zip fiG and fig lists together, calculate chi2 values
                if float(G)>0:
                    value=((float(g)-float(G))*(float(g)-float(G)))/float(G)
                    tetra_values.append(value)
                elif float(g)==0 and float(G)==0:
                    value=0
                    tetra_values.append(value)
                elif float(G)==0:
                    value=((float(g)-0.00000001)*(float(g)-0.00000001))/0.00000001
                    tetra_values.append(value)

            chi2_dist=sum(tetra_values)                     #sum individual comparison values to get final chi2 value for window
            print fig_name+"\t"+("%.5f" % chi2_dist)     #print window name<tab>chi2 value to output file specified in command line
            chi2_dist=0




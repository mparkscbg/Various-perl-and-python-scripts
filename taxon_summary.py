#!/usr/bin/python
from __future__ import division
import re
import os
import sys
from os import listdir
import sys
from collections import defaultdict
import operator
import numpy as np

#Author: Matthew Parks 2016

#this script takes a collection of gene trees in separate .tre files, and a list of the names of species found in those trees (see below),
#and returns the min, max and average total taxon count for the trees, as well as the min, max and average unique taxon count for the collection of trees
#
#syntax:
#python taxon_summary.py <taxon_names.txt> <tree_directory/>
#
#sys.argv[1] is file with all taxon names from trees, with each name followed by a colon
#sys.argv[2] is the directory with all of the trees of interest


with open(sys.argv[1]) as file:
    taxon_list=file.read().splitlines()                         #initialize list of taxon names from filenames.txt input file


counts=defaultdict(list)                #initialize empty dictionary called counts, will be populated first with taxon names, then with values of taxon presence per tree



file.close()



tree_taxon_count=0
tree_unique_taxon_count=0

taxon_count=[]
unique_taxon_count=[]


suffix=".tre"                       #will use so that script will only read files ending in ".tre"



f=open("tree_count_summaries.txt", "w")
for filename in listdir(sys.argv[2]):
    if filename.endswith(suffix):                           #only read files ending in ".tre"
        with open(sys.argv[2] + filename) as currentFile:   #open tree files one at a time
            text = currentFile.read()                       #read currently open tree file

            for name in taxon_list:                         #set up iteration of file over list of names
                ind_count=text.count(name)                  #count all instances of name, store count in variable ind_count
                counts[name].append(int(ind_count))         #append count of taxon in tree to its key's values in counts dictionary, for later use below

                tree_taxon_count+=text.count(name)          #add all instances of name to variable tree_taxon_count

                if (name in text):                          #test if name present in tree
                    tree_unique_taxon_count+=1              #if so, add to count of unique taxa for tree
            taxon_count.append(tree_taxon_count)
            unique_taxon_count.append(tree_unique_taxon_count)          #add number of unique taxa from tree to unique_taxon_count list

    taxon_count_sorted=sorted(taxon_count)
    unique_taxon_count_sorted=sorted(unique_taxon_count)

    tree_taxon_count=0
    tree_unique_taxon_count=0





f.write("Numerical summary for trees in: "+sys.argv[2]+"\n\n")

f.write("Min taxon count for all trees is: "+str(taxon_count_sorted[0])+"\n")
f.write("Max taxon count for all trees is: "+str(taxon_count_sorted[-1])+"\n")
f.write("Average taxon count per tree is: "+str(round(reduce(lambda x, y: x + y, taxon_count_sorted)/len(taxon_count_sorted),2))+"\n\n")
    
f.write("Min unique taxon count for all trees is: "+str(unique_taxon_count_sorted[0])+"\n")
f.write("Max unique taxon count for all trees is: "+str(unique_taxon_count_sorted[-1])+"\n")
f.write("Average unique taxon count per tree is: "+str(round(reduce(lambda x, y: x + y, unique_taxon_count_sorted)/len(unique_taxon_count_sorted),2))+"\n\n")

print "Numerical summary for trees in: "+sys.argv[2]+"\n\n"

print "Min taxon count for all trees is: "+str(taxon_count_sorted[0])
print "Max taxon count for all trees is: "+str(taxon_count_sorted[-1])
print "Average taxon count per tree is: "+str(round(reduce(lambda x, y: x + y, taxon_count_sorted)/len(taxon_count_sorted),2))+"\n"

print "Min unique taxon count for all trees is: "+str(unique_taxon_count_sorted[0])
print "Max unique taxon count for all trees is: "+str(unique_taxon_count_sorted[-1])
print "Average unique taxon count per tree is: "+str(round(reduce(lambda x, y: x + y, unique_taxon_count_sorted)/len(unique_taxon_count_sorted),2))+"\n"





averages={}

for name, numbers in counts.items():
    average=float(sum(numbers))/len(numbers)
    averages[name]=average

#print averages

minimums={}
for name, numbers in counts.items():
    minimum=(min(numbers))
    minimums[name]=minimum

#print minimums

maximums={}
for name, numbers in counts.items():
    maximum=(max(numbers))
    maximums[name]=maximum





print "transcripts per tree, for all "+str(len(taxon_count_sorted))+" trees in directory:\n(sample: average<tab>minimum<tab>maximum<tab>total trees occupied by sample (proportion of all trees))"
f.write("transcripts per tree, for all "+str(len(taxon_count_sorted))+" trees in directory:\n(sample: average<tab>minimum<tab>maximum<tab>total trees occupied by sample (proportion of all trees))\n")

sorted_averages=sorted(averages.items())
for k,v in sorted_averages:
    tree_count=np.count_nonzero(counts[k])
    tree_proportion=tree_count/len(counts[k])
    print k,"\t",round(v,2),"\t",minimums[k],"\t",maximums[k],"\t",tree_count,"(",round(tree_proportion,2),")"
    f.write(k+"\t"+str(round(v,2))+"\t"+str(minimums[k])+"\t"+str(maximums[k])+"\t"+str(tree_count)+" ("+str(round(tree_proportion,2))+")"+"\n")



print "\n\nThank you for using Python and have a great day!\n"
f.write("\n\nThank you for using Python and have a great day!\n")



f.close()


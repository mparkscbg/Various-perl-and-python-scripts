#!/usr/bin/python
from __future__ import division
import re
import os
import sys

#Author: Matthew Parks 2016

#usage: python filter_ortho_output.py orthofinder_output.txt <lower orthogroup size cutoff> <upper orthogroup size cutoff> <per-species transcript count cutoff>
#this script is used to filter orthofinder output for a specific size of orthologous cluster, as well
#as a certain level of uniqueness in the list of species contributing to the cluster

ortho_file=open(sys.argv[1],"r")

for line in ortho_file:
    lines=line.rstrip().split()
    transcript_count=len(lines)-1
    OG=lines[0]                                 #isolate orthogroup identifier
    transcript_array=lines[1:]                  #isolate transcripts from orthogroup identifier
    species_names=[]                            #initialize array for species names only
    for transcript in transcript_array:
        species_id=transcript.split(".")        #cut transcript names at "."
        species_names.append(species_id[0])     #append first field (species names) to species_names list
    species_names.sort()                        #sort species_names
    #print species_names                        #check for correctness of sorted species names
        
    species_counts=[]                           #initialize list for species counts


    uniq_list=list(set(species_names))          #create new list with only unique species names from original list
    #print str(len(uniq_list))                  #check for correctness of unique list

    for species in uniq_list:
        species_counts.append(species_names.count(species))     #count numbers of each species in original list
    species_counts.sort()                                       #sort list numbers so highest number is last
    if len(uniq_list)>=int(sys.argv[2]) and len(uniq_list)<=int(sys.argv[3]) and int(species_counts[-1])<=int(sys.argv[4]):               #compare highest species count to command line input argv[4]
        print OG







ortho_file.close()

















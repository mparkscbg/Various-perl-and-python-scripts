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
import pprint
from collections import Counter

#Author: Matthew Parks, 2016

#this script takes a collection of multiple sequence DNA alignments in fasta format in separate .fa, .fasta or .aln-cln files, and a list of the names of species found
#in those trees (see below), and returns the following:
#For each alignment:
#number of taxa
#alignment length
#total number of characters, including gaps
#average number of unique (non-gap) characters per alignment position
#total number of characters, excluding gaps
#total number of gap characters
#percentage of gap characters in alignment
#average gaps per alignment position
#number of variable sites in alignment
#number of parsimony informative sites in alignment
#
#For all files:
#min, max and average alignment length across all files
#min, max and average percent missing data/gaps across all files
#min, max and average proportion of variable positions in alignments across all files
#min, max and average proportion of parsimony informative positions in alignments across all files
#min, max and average percent missing data/gaps across all files for each taxon in alignments
#
#syntax:
#python alignments_summary.py <taxon_names.txt> <alignment_directory/>
#
#sys.argv[1] is file with all taxon names from alignments
#sys.argv[2] is the directory with all of the alignments of interest


with open(sys.argv[1]) as file:
    taxon_list=file.read().splitlines()                         #initialize list of taxon names from filenames.txt input file

alignment_names=defaultdict(list)                               #initialize empty dictionary called alignment_names, will be populated with alignment names and summary numbers from alignments
alignment_taxa=dict.fromkeys(taxon_list,[])                       #initialize alignment_taxa dictionary, populate with tax names from input taxon_names.txt (sys.argv[1])
alignment_collection={}

file.close()

suffix1=".fa"                   #will use these so that script will only read files ending in ".fa", ".fasta", ".fsa", or ".aln-cln"
suffix2=".fasta"
suffix3=".fsa"
suffix4=".aln-cln"

for_lengths=[]                  #initialize for_lengths list
alignment_dict={}               #initialize 'temporary' alignment dictionary
alignment_dict_separate={}      #initialize alignment dictionary for to hold values as character-separated sequence (i.e., each value is a character in the protein sequence corresponding to a key)
for_vars=[]                     #initialize list for calculating variable and parsimony informative counts
var_site_count=0                #initialize count for variable sites in an alignment at 0
pars_inf_site_count=0           #initilaize count for parsimony informative sites in an alignment at 0
uniq_count_list=[]              #initialize list for unique character counts at alignment positions, values to be averaged
al_col=[]                       #initialize list for alignment columns
tax_count=0                     #initialize taxon count for individual taxa in alignments
seq_length=0                    #initialize sequence lengths for individual taxa in alignments
tax_gap_count=0                 #initialize gap counts for individual taxa in alignments
taxon_summaries=[]              #initialze list to hold taxon summary numbers tuples

f=open("alignment_summaries.txt", "w")          #open text file for results to be written to


print"\nALIGNMENT_NAME\t#_TAXA_IN_ALIGNMENT\tALIGNMENT_LENGTH\tTOTAL_#_CHARS_INCL_GAPS\tAVG_#_UNIQ_(NON-GAP)_CHARS_PER_ALN_POS\tTOTAL_#_CHARS_EXCL_GAPS\tTOTAL_#_GAP_CHARS\tPERC_GAP_POS_IN_ALN\tAVG_GAPS_PER_ALN_POS\t#_VAR_POS_IN_ALN\tPROP_VAR_POS_IN_ALN\t#_PARS_INF_POS_IN_ALN\tPROP_PARS_INF_POS_IN_ALN"

for filename in listdir(sys.argv[2]):                                                                                                   #for files in directory
    if filename.endswith(suffix1) or filename.endswith(suffix2) or filename.endswith(suffix3) or filename.endswith(suffix4):            #only read files ending in ".fa" or ".fasta"
        with open(sys.argv[2] + filename) as currentFile:                                                                               #and open alignment files one at a time
            for line in currentFile:
                if line.startswith(">"):                                                                #separate fasta headers from sequences
                    sequence_name=line.rstrip().lstrip(">")
                else:
                    alignment_dict.setdefault(sequence_name, []).append(line.rstrip())
                    line_list=list(line.rstrip("\n"))                                                   #create character-separated list of sequence
                    for i, val in enumerate(line_list):                                                 #build new dictionary with key=taxa and values=character-separated list of the taxon's sequence
                        alignment_dict_separate.setdefault(sequence_name, []).append(val)
        
            alignment_collection[filename]=alignment_dict

            num_taxa=len(alignment_dict.keys())                                                         #get number of taxa in alignment
        
            val4len=str(alignment_dict.values()[0]).replace("'","").replace("[","").replace("]","")     #turn first key's value into a string, with no ', [ or ]
            aln_length=len(val4len)                                                                     #get length of first key's value, i.e., length of alignment
        
            for key in alignment_dict:
                for_lengths.append(alignment_dict[key])                                                 #append all sequences to for_lengths list
            

            all_seqs=''.join(map(str,for_lengths)).replace("'","").replace("[","").replace("]","").replace(",","").replace(" ","")      #concatenate all sequences, get rid of special characters associated with dictionary
    
            dash_count=all_seqs.count("-")                  #count dashes in all concatenated sequences in alignment
            total_length=len(all_seqs)                      #determine total length of all concatenated sequences in alignment
            gap_perc=round(dash_count/total_length*100,1)   #determine percentage of dashes in total length of all concatenated sequences in alignment

            alignment_dict_list=alignment_dict_separate.values()    #convert alignment_dict_separate values to list
            columns=zip(*alignment_dict_list)                       #convert list of alignment_dict_separate values into new list by position of lists

            for column in columns:                                  #iterate through columns
                al_col=list(column)                                 #convert tuple 'column' into list 'al_col'
                while "-" in al_col:                                #remove dashes from alignment position before calculating variable and parsimony sites
                    try:
                        al_col.remove("-")
                    except ValueError:
                        break
                #print al_col                                       #print check for al_col
                c=Counter(al_col)                                   #set counter for counting membership numbers in al_col
                k=c.most_common()                                   #order membership numbers from most to least common
                char_count=len(al_col)                                              #count number of characters in alignment position after dashes have been removed
                uniq_count=len(set(al_col))                                         #count number of unique characters in alignment position after dashes have been removed
                if uniq_count>1:                                                    #condition for a variable site
                    var_site_count+=1
                if uniq_count<len(al_col) and uniq_count>=2 and k[1][1]>=2:         #conditions for a parsimony informative site: 1) not every site is different, 2) at least two
                    pars_inf_site_count+=1                                          #different characters, 3) at least two characters present two or more times
                uniq_count_list.append(uniq_count)
                uniq_count=0
            avg_uniq_non_gap_chars=round(reduce(lambda x, y: x+y,uniq_count_list)/len(uniq_count_list),2)


            print filename+"\t"+str(num_taxa)+"\t"+str(aln_length)+"\t"+str(total_length)+"\t"+str(round(reduce(lambda x, y: x+y,uniq_count_list)/len(uniq_count_list),2))+"\t"+str(total_length-dash_count)+"\t"+str(dash_count)+"\t"+str(gap_perc)+"\t"+str(round((dash_count/aln_length),2))+"\t"+str(var_site_count)+"\t"+str(round(var_site_count/aln_length,2))+"\t"+str(pars_inf_site_count)+"\t"+str(round(pars_inf_site_count/aln_length,2))
            
            

            alignment_names[filename].append(num_taxa)            #append values to new dictionary alignment_names, using input filenames for keys
            alignment_names[filename].append(aln_length)
            alignment_names[filename].append(total_length)
            alignment_names[filename].append(avg_uniq_non_gap_chars)
            alignment_names[filename].append(total_length-dash_count)
            alignment_names[filename].append(dash_count)
            alignment_names[filename].append(gap_perc)
            alignment_names[filename].append(round(dash_count/total_length,2))
            alignment_names[filename].append(round(var_site_count/aln_length,2))
            alignment_names[filename].append(round(pars_inf_site_count/aln_length,2))



            alignment_dict={}                   #reset alignment dictionary to empty
            alignment_dict_separate={}          #reset alignment_dict_separate to empty


            for_lengths=[]                      #reset list of all sequences to empty
            line_list=[]                        #reset character-separated list of sequence
            for_vars=[]                         #reset list for variable and parsimony informative calculations
            var_site_count=0                    #reset variable site count to 0
            pars_inf_site_count=0               #reset parsimony informative site count to 0
            diverse_prop_list=[]                #reset div_prop_list to empty
            uniq_count_list=[]                  #reset uniq_count_list



alignment_final=alignment_names.values()                    #convert alignment_dict_separate values to list
final_columns=zip(*alignment_final)                         #convert list of alignment_dict_separate values into new list by position of lists

aln_length_final=list(final_columns[1])
perc_gaps_final=list(final_columns[6])
prop_var_pos_final=list(final_columns[8])
prop_pars_inf_pos_final=list(final_columns[9])



aln_length_final.sort()
perc_gaps_final.sort()
prop_var_pos_final.sort()
prop_pars_inf_pos_final.sort()


avg_aln_length_final=round(reduce(lambda x, y: x+y,aln_length_final)/len(aln_length_final),2)
avg_perc_gaps_final=round(reduce(lambda x, y: x+y,perc_gaps_final)/len(perc_gaps_final),2)
avg_prop_var_pos_final=round(reduce(lambda x, y: x+y,prop_var_pos_final)/len(prop_var_pos_final),2)
avg_prop_pars_inf_pos_final=round(reduce(lambda x, y: x+y,prop_pars_inf_pos_final)/len(prop_pars_inf_pos_final),2)


print "\n"
print "CATEGORY\tAVERAGE\tMINIMUM\tMAXIMUM"
print "aln_length\t"+str(avg_aln_length_final)+"\t"+str(aln_length_final[0])+"\t"+str(aln_length_final[-1])
print "perc_gaps\t"+str(avg_perc_gaps_final)+"\t"+str(perc_gaps_final[0])+"\t"+str(perc_gaps_final[-1])
print "prop_var_pos\t"+str(avg_prop_var_pos_final)+"\t"+str(prop_var_pos_final[0])+"\t"+str(prop_var_pos_final[-1])
print "prop_pars_inf_pos\t"+str(avg_prop_pars_inf_pos_final)+"\t"+str(prop_pars_inf_pos_final[0])+"\t"+str(prop_pars_inf_pos_final[-1])+"\n"



print "TAXON\tALN_COUNT\tTOTAL_SEQ_LENGTH\tGAP_COUNT\tGAP_PROP"

for taxon in taxon_list:
    for alignment_name in alignment_collection:
        try:
            seq=alignment_collection[alignment_name][taxon]
            #print seq[0]
            tax_count+=1
            seq_length+=len(seq[0])
            tax_gap_count+=seq[0].count('-')
        except KeyError:
            pass

    try:
        tax_gap_prop=round(tax_gap_count/seq_length,2)
    except ZeroDivisionError:
        tax_gap_prop=float('NaN')

    print taxon+"\t"+str(tax_count)+"\t"+str(seq_length)+"\t"+str(tax_gap_count)+"\t"+str(tax_gap_prop)

    tax_count=0
    seq_length=0
    tax_gap_count=0

print "\n\n\nHoller!\n"


f.close()







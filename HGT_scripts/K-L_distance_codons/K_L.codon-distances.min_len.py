

#This script takes a multi-sequence fasta file as input, with each fasta sequence representing a
#coding sequence (must be in frame) from a transcriptome. The *.cds output from TransDecoder is a
#good example.
#The script calculates the Kullback-Leibler (K-L) distance for codon usage for each coding sequence
#against the entire collection of coding sequences in the supplied fasta file, given a minimum sequence length cutoff.
#Syntax is:

#python K_L.codon_distances.py input_cds.fasta min_length_cutoff > output.txt

#Output is two tab-separated columns - the name of each coding sequence (from its fasta header) and
#its K-L distance. Ouput K-L values can be plotted separately for outlier analysis.
#
#NOTE: this script works pretty fast, it typically gets through ca. 20,000 transcripts in around 10-20 minutes on a single processor.




from __future__ import division

import re
import sys
import math
from Bio import SeqIO


min_len=sys.argv[2]
seq_dict={}

CodonsDict = {'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0, 'ACA': 0,
                'ACC': 0, 'ACG': 0, 'ACT': 0, 'AGA': 0, 'AGC': 0,
                'AGG': 0, 'AGT': 0, 'ATA': 0, 'ATC': 0, 'ATG': 0,
                'ATT': 0, 'CAA': 0, 'CAC': 0, 'CAG': 0, 'CAT': 0,
                'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0, 'CGA': 0,
                'CGC': 0, 'CGG': 0, 'CGT': 0, 'CTA': 0, 'CTC': 0,
                'CTG': 0, 'CTT': 0, 'GAA': 0, 'GAC': 0, 'GAG': 0,
                'GAT': 0, 'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0,
                'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0, 'GTA': 0,
                'GTC': 0, 'GTG': 0, 'GTT': 0, 'TAA': 0, 'TAC': 0,
                'TAG': 0, 'TAT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0,
                'TCT': 0, 'TGA': 0, 'TGC': 0, 'TGG': 0, 'TGT': 0,
                'TTA': 0, 'TTC': 0, 'TTG': 0, 'TTT': 0}


for record in SeqIO.parse(sys.argv[1],"fasta"):
    name=record.id
    sequence=record.seq
    #print sequence
    if len(record.seq)>int(min_len)-1:
        #print record.seq
        temp_dict= {'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0, 'ACA': 0,
                    'ACC': 0, 'ACG': 0, 'ACT': 0, 'AGA': 0, 'AGC': 0,
                    'AGG': 0, 'AGT': 0, 'ATA': 0, 'ATC': 0, 'ATG': 0,
                    'ATT': 0, 'CAA': 0, 'CAC': 0, 'CAG': 0, 'CAT': 0,
                    'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0, 'CGA': 0,
                    'CGC': 0, 'CGG': 0, 'CGT': 0, 'CTA': 0, 'CTC': 0,
                    'CTG': 0, 'CTT': 0, 'GAA': 0, 'GAC': 0, 'GAG': 0,
                    'GAT': 0, 'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0,
                    'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0, 'GTA': 0,
                    'GTC': 0, 'GTG': 0, 'GTT': 0, 'TAA': 0, 'TAC': 0,
                    'TAG': 0, 'TAT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0,
                    'TCT': 0, 'TGA': 0, 'TGC': 0, 'TGG': 0, 'TGT': 0,
                    'TTA': 0, 'TTC': 0, 'TTG': 0, 'TTT': 0}




        for i in xrange(0,len(sequence),3):
            codon=sequence[i:i+3]
            #print codon

            if codon.upper() == 'AAA':
                CodonsDict['AAA']+=1
                temp_dict['AAA']+=1
            elif codon.upper() == 'AAC':
                CodonsDict['AAC']+=1
                temp_dict['AAC']+=1
            elif codon.upper() == 'AAG':
                CodonsDict['AAG']+=1
                temp_dict['AAG']+=1
            elif codon.upper() == 'AAT':
                CodonsDict['AAT']+=1
                temp_dict['AAT']+=1
            elif codon.upper() == 'ACA':
                CodonsDict['ACA']+=1
                temp_dict['ACA']+=1
            elif codon.upper() == 'ACC':
                CodonsDict['ACC']+=1
                temp_dict['ACC']+=1
            elif codon.upper() == 'ACG':
                CodonsDict['ACG']+=1
                temp_dict['ACG']+=1
            elif codon.upper() == 'ACT':
                CodonsDict['ACT']+=1
                temp_dict['ACT']+=1
            elif codon.upper() == 'AGA':
                CodonsDict['AGA']+=1
                temp_dict['AGA']+=1
            elif codon.upper() == 'AGC':
                CodonsDict['AGC']+=1
                temp_dict['AGC']+=1
            elif codon.upper() == 'AGG':
                CodonsDict['AGG']+=1
                temp_dict['AGG']+=1
            elif codon.upper() == 'AGT':
                CodonsDict['AGT']+=1
                temp_dict['AGT']+=1
            elif codon.upper() == 'ATA':
                CodonsDict['ATA']+=1
                temp_dict['ATA']+=1
            elif codon.upper() == 'ATC':
                CodonsDict['ATC']+=1
                temp_dict['ATC']+=1
            elif codon.upper() == 'ATG':
                CodonsDict['ATG']+=1
                temp_dict['ATG']+=1
            elif codon.upper() == 'ATT':
                CodonsDict['ATT']+=1
                temp_dict['ATT']+=1
            elif codon.upper() == 'CAA':
                CodonsDict['CAA']+=1
                temp_dict['CAA']+=1
            elif codon.upper() == 'CAC':
                CodonsDict['CAC']+=1
                temp_dict['CAC']+=1
            elif codon.upper() == 'CAG':
                CodonsDict['CAG']+=1
                temp_dict['CAG']+=1
            elif codon.upper() == 'CAT':
                CodonsDict['CAT']+=1
                temp_dict['CAT']+=1
            elif codon.upper() == 'CCA':
                CodonsDict['CCA']+=1
                temp_dict['CCA']+=1
            elif codon.upper() == 'CCC':
                CodonsDict['CCC']+=1
                temp_dict['CCC']+=1
            elif codon.upper() == 'CCG':
                CodonsDict['CCG']+=1
                temp_dict['CCG']+=1
            elif codon.upper() == 'CCT':
                CodonsDict['CCT']+=1
                temp_dict['CCT']+=1
            elif codon.upper() == 'CGA':
                CodonsDict['CGA']+=1
                temp_dict['CGA']+=1
            elif codon.upper() == 'CGC':
                CodonsDict['CGC']+=1
                temp_dict['CGC']+=1
            elif codon.upper() == 'CGG':
                CodonsDict['CGG']+=1
                temp_dict['CGG']+=1
            elif codon.upper() == 'CGT':
                CodonsDict['CGT']+=1
                temp_dict['CGT']+=1
            elif codon.upper() == 'CTA':
                CodonsDict['CTA']+=1
                temp_dict['CTA']+=1
            elif codon.upper() == 'CTC':
                CodonsDict['CTC']+=1
                temp_dict['CTC']+=1
            elif codon.upper() == 'CTG':
                CodonsDict['CTG']+=1
                temp_dict['CTG']+=1
            elif codon.upper() == 'CTT':
                CodonsDict['CTT']+=1
                temp_dict['CTT']+=1
            elif codon.upper() == 'GAA':
                CodonsDict['GAA']+=1
                temp_dict['GAA']+=1
            elif codon.upper() == 'GAC':
                CodonsDict['GAC']+=1
                temp_dict['GAC']+=1
            elif codon.upper() == 'GAG':
                CodonsDict['GAG']+=1
                temp_dict['GAG']+=1
            elif codon.upper() == 'GAT':
                CodonsDict['GAT']+=1
                temp_dict['GAT']+=1
            elif codon.upper() == 'GCA':
                CodonsDict['GCA']+=1
                temp_dict['GCA']+=1
            elif codon.upper() == 'GCC':
                CodonsDict['GCC']+=1
                temp_dict['GCC']+=1
            elif codon.upper() == 'GCG':
                CodonsDict['GCG']+=1
                temp_dict['GCG']+=1
            elif codon.upper() == 'GCT':
                CodonsDict['GCT']+=1
                temp_dict['GCT']+=1
            elif codon.upper() == 'GGA':
                CodonsDict['GGA']+=1
                temp_dict['GGA']+=1
            elif codon.upper() == 'GGC':
                CodonsDict['GGC']+=1
                temp_dict['GGC']+=1
            elif codon.upper() == 'GGG':
                CodonsDict['GGG']+=1
                temp_dict['GGG']+=1
            elif codon.upper() == 'GGT':
                CodonsDict['GGT']+=1
                temp_dict['GGT']+=1
            elif codon.upper() == 'GTA':
                CodonsDict['GTA']+=1
                temp_dict['GTA']+=1
            elif codon.upper() == 'GTC':
                CodonsDict['GTC']+=1
                temp_dict['GTC']+=1
            elif codon.upper() == 'GTG':
                CodonsDict['GTG']+=1
                temp_dict['GTG']+=1
            elif codon.upper() == 'GTT':
                CodonsDict['GTT']+=1
                temp_dict['GTT']+=1
            elif codon.upper() == 'TAA':
                CodonsDict['TAA']+=1
                temp_dict['TAA']+=1
            elif codon.upper() == 'TAC':
                CodonsDict['TAC']+=1
                temp_dict['TAC']+=1
            elif codon.upper() == 'TAG':
                CodonsDict['TAG']+=1
                temp_dict['TAG']+=1
            elif codon.upper() == 'TAT':
                CodonsDict['TAT']+=1
                temp_dict['TAT']+=1
            elif codon.upper() == 'TCA':
                CodonsDict['TCA']+=1
                temp_dict['TCA']+=1
            elif codon.upper() == 'TCC':
                CodonsDict['TCC']+=1
                temp_dict['TCC']+=1
            elif codon.upper() == 'TCG':
                CodonsDict['TCG']+=1
                temp_dict['TCG']+=1
            elif codon.upper() == 'TCT':
                CodonsDict['TCT']+=1
                temp_dict['TCT']+=1
            elif codon.upper() == 'TGA':
                CodonsDict['TGA']+=1
                temp_dict['TGA']+=1
            elif codon.upper() == 'TGC':
                CodonsDict['TGC']+=1
                temp_dict['TGC']+=1
            elif codon.upper() == 'TGG':
                CodonsDict['TGG']+=1
                temp_dict['TGG']+=1
            elif codon.upper() == 'TGT':
                CodonsDict['TGT']+=1
                temp_dict['TGT']+=1
            elif codon.upper() == 'TTA':
                CodonsDict['TTA']+=1
                temp_dict['TTA']+=1
            elif codon.upper() == 'TTC':
                CodonsDict['TTC']+=1
                temp_dict['TTC']+=1
            elif codon.upper() == 'TTG':
                CodonsDict['TTG']+=1
                temp_dict['TTG']+=1
            elif codon.upper() == 'TTT':
                CodonsDict['TTT']+=1
                temp_dict['TTT']+=1


        z={name:temp_dict}  #use variable key name to avoid variable variables
        seq_dict.update(z)  #update sequence-specific dictionary by adding temp_dict dictionary to it with sequence name as key


total_codons=sum(CodonsDict.values())

fiG_AAA=round(CodonsDict['AAA']/total_codons,5)
fiG_AAC=round(CodonsDict['AAC']/total_codons,5)
fiG_AAG=round(CodonsDict['AAG']/total_codons,5)
fiG_AAT=round(CodonsDict['AAT']/total_codons,5)
fiG_ACA=round(CodonsDict['ACA']/total_codons,5)
fiG_ACC=round(CodonsDict['ACC']/total_codons,5)
fiG_ACG=round(CodonsDict['ACG']/total_codons,5)
fiG_ACT=round(CodonsDict['ACT']/total_codons,5)
fiG_AGA=round(CodonsDict['AGA']/total_codons,5)
fiG_AGC=round(CodonsDict['AGC']/total_codons,5)
fiG_AGG=round(CodonsDict['AGG']/total_codons,5)
fiG_AGT=round(CodonsDict['AGT']/total_codons,5)
fiG_ATA=round(CodonsDict['ATA']/total_codons,5)
fiG_ATC=round(CodonsDict['ATC']/total_codons,5)
fiG_ATG=round(CodonsDict['ATG']/total_codons,5)
fiG_ATT=round(CodonsDict['ATT']/total_codons,5)
fiG_CAA=round(CodonsDict['CAA']/total_codons,5)
fiG_CAC=round(CodonsDict['CAC']/total_codons,5)
fiG_CAG=round(CodonsDict['CAG']/total_codons,5)
fiG_CAT=round(CodonsDict['CAT']/total_codons,5)
fiG_CCA=round(CodonsDict['CCA']/total_codons,5)
fiG_CCC=round(CodonsDict['CCC']/total_codons,5)
fiG_CCG=round(CodonsDict['CCG']/total_codons,5)
fiG_CCT=round(CodonsDict['CCT']/total_codons,5)
fiG_CGA=round(CodonsDict['CGA']/total_codons,5)
fiG_CGC=round(CodonsDict['CGC']/total_codons,5)
fiG_CGG=round(CodonsDict['CGG']/total_codons,5)
fiG_CGT=round(CodonsDict['CGT']/total_codons,5)
fiG_CTA=round(CodonsDict['CTA']/total_codons,5)
fiG_CTC=round(CodonsDict['CTC']/total_codons,5)
fiG_CTG=round(CodonsDict['CTG']/total_codons,5)
fiG_CTT=round(CodonsDict['CTT']/total_codons,5)
fiG_GAA=round(CodonsDict['GAA']/total_codons,5)
fiG_GAC=round(CodonsDict['GAC']/total_codons,5)
fiG_GAG=round(CodonsDict['GAG']/total_codons,5)
fiG_GAT=round(CodonsDict['GAT']/total_codons,5)
fiG_GCA=round(CodonsDict['GCA']/total_codons,5)
fiG_GCC=round(CodonsDict['GCC']/total_codons,5)
fiG_GCG=round(CodonsDict['GCG']/total_codons,5)
fiG_GCT=round(CodonsDict['GCT']/total_codons,5)
fiG_GGA=round(CodonsDict['GGA']/total_codons,5)
fiG_GGC=round(CodonsDict['GGC']/total_codons,5)
fiG_GGG=round(CodonsDict['GGG']/total_codons,5)
fiG_GGT=round(CodonsDict['GGT']/total_codons,5)
fiG_GTA=round(CodonsDict['GTA']/total_codons,5)
fiG_GTC=round(CodonsDict['GTC']/total_codons,5)
fiG_GTG=round(CodonsDict['GTG']/total_codons,5)
fiG_GTT=round(CodonsDict['GTT']/total_codons,5)
fiG_TAA=round(CodonsDict['TAA']/total_codons,5)
fiG_TAC=round(CodonsDict['TAC']/total_codons,5)
fiG_TAG=round(CodonsDict['TAG']/total_codons,5)
fiG_TAT=round(CodonsDict['TAT']/total_codons,5)
fiG_TCA=round(CodonsDict['TCA']/total_codons,5)
fiG_TCC=round(CodonsDict['TCC']/total_codons,5)
fiG_TCG=round(CodonsDict['TCG']/total_codons,5)
fiG_TCT=round(CodonsDict['TCT']/total_codons,5)
fiG_TGA=round(CodonsDict['TGA']/total_codons,5)
fiG_TGC=round(CodonsDict['TGC']/total_codons,5)
fiG_TGG=round(CodonsDict['TGG']/total_codons,5)
fiG_TGT=round(CodonsDict['TGT']/total_codons,5)
fiG_TTA=round(CodonsDict['TTA']/total_codons,5)
fiG_TTC=round(CodonsDict['TTC']/total_codons,5)
fiG_TTG=round(CodonsDict['TTG']/total_codons,5)
fiG_TTT=round(CodonsDict['TTT']/total_codons,5)

fiG_list=(fiG_AAA,fiG_AAC,fiG_AAG,fiG_AAT,fiG_ACA,fiG_ACC,fiG_ACG,fiG_ACT,fiG_AGA,fiG_AGC,fiG_AGG,fiG_AGT,fiG_ATA,fiG_ATC,fiG_ATG,fiG_ATT,fiG_CAA,fiG_CAC,fiG_CAG,fiG_CAT,fiG_CCA,fiG_CCC,fiG_CCG,fiG_CCT,fiG_CGA,fiG_CGC,fiG_CGG,fiG_CGT,fiG_CTA,fiG_CTC,fiG_CTG,fiG_CTT,fiG_GAA,fiG_GAC,fiG_GAG,fiG_GAT,fiG_GCA,fiG_GCC,fiG_GCG,fiG_GCT,fiG_GGA,fiG_GGC,fiG_GGG,fiG_GGT,fiG_GTA,fiG_GTC,fiG_GTG,fiG_GTT,fiG_TAA,fiG_TAC,fiG_TAG,fiG_TAT,fiG_TCA,fiG_TCC,fiG_TCG,fiG_TCT,fiG_TGA,fiG_TGC,fiG_TGG,fiG_TGT,fiG_TTA,fiG_TTC,fiG_TTG,fiG_TTT)
#print fiG_list

for key in sorted(seq_dict.iterkeys()):
    seq_name=key        #set name for sequence
    #print seq_name
    new_dict=seq_dict[key]      #create new dictionary from nested dictionary
    #print new_dict
    key_list=[k for k in new_dict]      #create new list of keys from new dictionary from nested dictionary
    value_list=[v for v in new_dict.values()]   #create new list of values from new dictionary from nested dictionary
    total_codons_seq=sum(new_dict.values())     #calculate total codon count for sequence
    #print total_codons_seq
    #print key_list
    #print value_list
    value_list_sorted=[y for (x,y) in sorted(zip(key_list,value_list))]     #sort value_list by alphabetically sorted list of codons from key_list
    #print value_list_sorted
    fig_AAA=round(value_list_sorted[0]/total_codons_seq,5)
    fig_AAC=round(value_list_sorted[1]/total_codons_seq,5)
    fig_AAG=round(value_list_sorted[2]/total_codons_seq,5)
    fig_AAT=round(value_list_sorted[3]/total_codons_seq,5)
    fig_ACA=round(value_list_sorted[4]/total_codons_seq,5)
    fig_ACC=round(value_list_sorted[5]/total_codons_seq,5)
    fig_ACG=round(value_list_sorted[6]/total_codons_seq,5)
    fig_ACT=round(value_list_sorted[7]/total_codons_seq,5)
    fig_AGA=round(value_list_sorted[8]/total_codons_seq,5)
    fig_AGC=round(value_list_sorted[9]/total_codons_seq,5)
    fig_AGG=round(value_list_sorted[10]/total_codons_seq,5)
    fig_AGT=round(value_list_sorted[11]/total_codons_seq,5)
    fig_ATA=round(value_list_sorted[12]/total_codons_seq,5)
    fig_ATC=round(value_list_sorted[13]/total_codons_seq,5)
    fig_ATG=round(value_list_sorted[14]/total_codons_seq,5)
    fig_ATT=round(value_list_sorted[15]/total_codons_seq,5)
    fig_CAA=round(value_list_sorted[16]/total_codons_seq,5)
    fig_CAC=round(value_list_sorted[17]/total_codons_seq,5)
    fig_CAG=round(value_list_sorted[18]/total_codons_seq,5)
    fig_CAT=round(value_list_sorted[19]/total_codons_seq,5)
    fig_CCA=round(value_list_sorted[20]/total_codons_seq,5)
    fig_CCC=round(value_list_sorted[21]/total_codons_seq,5)
    fig_CCG=round(value_list_sorted[22]/total_codons_seq,5)
    fig_CCT=round(value_list_sorted[23]/total_codons_seq,5)
    fig_CGA=round(value_list_sorted[24]/total_codons_seq,5)
    fig_CGC=round(value_list_sorted[25]/total_codons_seq,5)
    fig_CGG=round(value_list_sorted[26]/total_codons_seq,5)
    fig_CGT=round(value_list_sorted[27]/total_codons_seq,5)
    fig_CTA=round(value_list_sorted[28]/total_codons_seq,5)
    fig_CTC=round(value_list_sorted[29]/total_codons_seq,5)
    fig_CTG=round(value_list_sorted[30]/total_codons_seq,5)
    fig_CTT=round(value_list_sorted[31]/total_codons_seq,5)
    fig_GAA=round(value_list_sorted[32]/total_codons_seq,5)
    fig_GAC=round(value_list_sorted[33]/total_codons_seq,5)
    fig_GAG=round(value_list_sorted[34]/total_codons_seq,5)
    fig_GAT=round(value_list_sorted[35]/total_codons_seq,5)
    fig_GCA=round(value_list_sorted[36]/total_codons_seq,5)
    fig_GCC=round(value_list_sorted[37]/total_codons_seq,5)
    fig_GCG=round(value_list_sorted[38]/total_codons_seq,5)
    fig_GCT=round(value_list_sorted[39]/total_codons_seq,5)
    fig_GGA=round(value_list_sorted[40]/total_codons_seq,5)
    fig_GGC=round(value_list_sorted[41]/total_codons_seq,5)
    fig_GGG=round(value_list_sorted[42]/total_codons_seq,5)
    fig_GGT=round(value_list_sorted[43]/total_codons_seq,5)
    fig_GTA=round(value_list_sorted[44]/total_codons_seq,5)
    fig_GTC=round(value_list_sorted[45]/total_codons_seq,5)
    fig_GTG=round(value_list_sorted[46]/total_codons_seq,5)
    fig_GTT=round(value_list_sorted[47]/total_codons_seq,5)
    fig_TAA=round(value_list_sorted[48]/total_codons_seq,5)
    fig_TAC=round(value_list_sorted[49]/total_codons_seq,5)
    fig_TAG=round(value_list_sorted[50]/total_codons_seq,5)
    fig_TAT=round(value_list_sorted[51]/total_codons_seq,5)
    fig_TCA=round(value_list_sorted[52]/total_codons_seq,5)
    fig_TCC=round(value_list_sorted[53]/total_codons_seq,5)
    fig_TCG=round(value_list_sorted[54]/total_codons_seq,5)
    fig_TCT=round(value_list_sorted[55]/total_codons_seq,5)
    fig_TGA=round(value_list_sorted[56]/total_codons_seq,5)
    fig_TGC=round(value_list_sorted[57]/total_codons_seq,5)
    fig_TGG=round(value_list_sorted[58]/total_codons_seq,5)
    fig_TGT=round(value_list_sorted[59]/total_codons_seq,5)
    fig_TTA=round(value_list_sorted[60]/total_codons_seq,5)
    fig_TTC=round(value_list_sorted[61]/total_codons_seq,5)
    fig_TTG=round(value_list_sorted[62]/total_codons_seq,5)
    fig_TTT=round(value_list_sorted[63]/total_codons_seq,5)

    fig_list=(fig_AAA,fig_AAC,fig_AAG,fig_AAT,fig_ACA,fig_ACC,fig_ACG,fig_ACT,fig_AGA,fig_AGC,fig_AGG,fig_AGT,fig_ATA,fig_ATC,fig_ATG,fig_ATT,fig_CAA,fig_CAC,fig_CAG,fig_CAT,fig_CCA,fig_CCC,fig_CCG,fig_CCT,fig_CGA,fig_CGC,fig_CGG,fig_CGT,fig_CTA,fig_CTC,fig_CTG,fig_CTT,fig_GAA,fig_GAC,fig_GAG,fig_GAT,fig_GCA,fig_GCC,fig_GCG,fig_GCT,fig_GGA,fig_GGC,fig_GGG,fig_GGT,fig_GTA,fig_GTC,fig_GTG,fig_GTT,fig_TAA,fig_TAC,fig_TAG,fig_TAT,fig_TCA,fig_TCC,fig_TCG,fig_TCT,fig_TGA,fig_TGC,fig_TGG,fig_TGT,fig_TTA,fig_TTC,fig_TTG,fig_TTT)

    #print fig_list

    seq_values=[]

    for g,G in zip(fig_list,fiG_list):
        if g>0 and G>0:
            value=g*(math.log(g/G))
            seq_values.append(value)
        elif g==0 and G>0:
            value=g*(math.log(0.00000004/G))
            seq_values.append(value)
        elif g>0 and G==0:
            value=g*(math.log(g/0.00000004))
            seq_values.append(value)
        elif g==0 and G==0:
            value=0
            seq_values.append(value)

    dist=sum(seq_values)
    print seq_name+"\t"+("%.5f" % dist)



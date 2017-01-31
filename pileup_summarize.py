#!/usr/bin/python
from __future__ import division
import sys
import re


#Author: Matthew Parks, 2016

inFile = open(sys.argv[1],'r')


print "chrom\tpos\tref_base\tA_count\tG_count\tC_count\tT_count\tN_count\tA_qual_avg\tG_qual_avg\tC_qual_avg\tT_qual_avg\tN_qual_avg\t\tread_count_insertion\tread_count_deletion\tprop_reads_insertion\tprop_reads_deletion"

for line in inFile:
    data = line.rstrip("\n").split('\t')
    if "-" not in data[4] and "+" not in data[4]:
        chrom = data[0]
        bp = data[1]
        ref_base = data[2].upper()
        bases = data[4].upper()
        base_calls = re.sub('[$]|[\^].','',bases)
        base_calls = base_calls.replace('.',ref_base).replace(',',ref_base).replace('*',ref_base)
        #print base_calls
        qual_scores=data[5]
        #print qual_scores
        qual_dict={"A":[],"T":[],"C":[],"G":[],"N":[]}

        for call in xrange(len(base_calls)):
            qual_dict[base_calls[call]].append(ord(qual_scores[call])-33)

            if (len(qual_dict["A"]) != 0):
                A=round(sum(qual_dict["A"])/len(qual_dict["A"]),3)
            else:
                A="0"

            if (len(qual_dict["G"]) != 0):
                G=round(sum(qual_dict["G"])/len(qual_dict["G"]),3)
            else:
                G="0"

            if (len(qual_dict["C"]) != 0):
                C=round(sum(qual_dict["C"])/len(qual_dict["C"]),3)
            else:
                C="0"

            if (len(qual_dict["T"]) != 0):
                T=round(sum(qual_dict["T"])/len(qual_dict["T"]),3)
            else:
                T="0"

            if (len(qual_dict["N"]) != 0):
                N=round(sum(qual_dict["N"])/len(qual_dict["N"]),3)
            else:
                N="0"



        print chrom+"\t"+bp+"\t"+ref_base+"\t"+str(len(qual_dict["A"]))+"\t"+str(len(qual_dict["G"]))+"\t"+str(len(qual_dict["C"]))+"\t"+str(len(qual_dict["T"]))+"\t"+str(len(qual_dict["N"]))+"\t"+str(A)+"\t"+str(G)+"\t"+str(C)+"\t"+str(T)+"\t"+str(N)+"\t0\t0\t0.0\t0.0"



    else:
        chrom = data[0]
        bp = data[1]
        bases = data[4].upper()
        base_calls = re.sub('[$]|[\^].','',bases)
        base_calls = base_calls.replace('.',ref_base).replace(',',ref_base).replace('*',ref_base)
        base_calls_list=list(base_calls)
        #print base_calls_list
        qual_scores=data[5]
        
        insertion_prop=str(round(bases.count('+')/len(qual_scores),3))
        deletion_prop=str(round(bases.count('-')/len(qual_scores),3))
        
        #print chrom+"\t"+bp+"\t"+str(bases.count('+'))+"\t"+"and"+"\t"+str(bases.count("-"))+"\t"+"read(s) out of "+str(len(qual_scores))+" reads have/has an insertion and deletion, respectively"

        numbers=re.findall(r'[+-]\d+',base_calls)
        numbers=[abs(int(i)) for i in numbers]
        #print str(numbers)+"\tindel sizes"
        #print numbers[::-1]
        numbers_rev=numbers[::-1]
        #print str(numbers_rev)+"\tindel sizes reversed"
        #numbers_rev=numbers.reverse()
        #print numbers.reverse()


        indel_starts=[]
        for starts in re.finditer('[+-]',base_calls):
            indel_starts.append(starts.start())
        #print str(indel_starts)+"\tindel start positions"


        for element in sorted(indel_starts, reverse=True):
            position=indel_starts.index(element)
            number_length=len(str(numbers[position]))
            #print numbers[position]
            del base_calls_list[element:element+int(numbers[position])+number_length+1]
            #del base_calls_list[element:element+3]
            
        base_calls_no_indels=''.join(base_calls_list)
        #print base_calls
        #print base_calls_no_indels
        
        qual_dict={"A":[],"T":[],"C":[],"G":[],"N":[]}

        for call in xrange(len(base_calls_no_indels)):
            qual_dict[base_calls_no_indels[call]].append(ord(qual_scores[call])-33)
        
            if (len(qual_dict["A"]) != 0):
                A=round(sum(qual_dict["A"])/len(qual_dict["A"]),3)
            else:
                A="0"
    
            if (len(qual_dict["G"]) != 0):
                G=round(sum(qual_dict["G"])/len(qual_dict["G"]),3)
            else:
                G="0"
        
            if (len(qual_dict["C"]) != 0):
                C=round(sum(qual_dict["C"])/len(qual_dict["C"]),3)
            else:
                C="0"

            if (len(qual_dict["T"]) != 0):
                T=round(sum(qual_dict["T"])/len(qual_dict["T"]),3)
            else:
                T="0"
        
            if (len(qual_dict["N"]) != 0):
                N=round(sum(qual_dict["N"])/len(qual_dict["N"]),3)
            else:
                N="0"



        print chrom+"\t"+bp+"\t"+ref_base+"\t"+str(len(qual_dict["A"]))+"\t"+str(len(qual_dict["G"]))+"\t"+str(len(qual_dict["C"]))+"\t"+str(len(qual_dict["T"]))+"\t"+str(len(qual_dict["N"]))+"\t"+str(A)+"\t"+str(G)+"\t"+str(C)+"\t"+str(T)+"\t"+str(N)+"\t"+str(bases.count('+'))+"\t"+str(bases.count("-"))+"\t"+insertion_prop+"\t"+deletion_prop









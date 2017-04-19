

#This script is modified from a set of posts on StackOverflow (http://stackoverflow.com/questions/31757876/python-find-longest-orf-in-dna-sequence)
#The script reads sequences from a multi-fasta file (in both directions? - not sure), and does the following:
#1) identifies the longest open reading frame that starts with a start codon (ATG) and ends
#   either at a stop codon (TAG, TGA, TAA) or the end of the sequence
#2) shortens the fasta header at the first whitespace
#3) prints out the new headers and corresponding sequences
#
#NOTE 1: if no stop codon is found, the script will print out the ORF until either the end of the sequence of the last multiple of three
#nucleotides (and so ignoring the final one or two nucleotides)
#
#NOTE 2: this script also filters for a specified minimum ORF length, to 'ignore' just set this at 1
#
#syntax
# python get_longest_orf.py input.fasta min_ORF_length > output.fasta



import sys
import math
from Bio import SeqIO
from Bio import Seq
import regex as re

min_len=sys.argv[2]


for record in SeqIO.parse(sys.argv[1],"fasta"):
    startP = re.compile('ATG')
    name=record.description
    sequence = record.seq
    input_seq=str(sequence)
    nuc = input_seq.replace('\n','')
    longest = (0,)
    for m in startP.finditer(nuc, overlapped=True):
        if len(Seq.Seq(nuc)[m.start():].translate(to_stop=True)) > longest[0]:
            pro = Seq.Seq(nuc)[m.start():].translate(to_stop=True)
            longest = (nuc[m.start():m.start()+len(pro)*3+3])
            if len(longest)>int(min_len)-1 and len(longest) % 3 == 0:
                print ">"+name
                print longest
            elif len(longest)>int(min_len)-1 and len(longest) % 3 != 0:
                longest_cut = (nuc[m.start():m.start()+len(pro)*3])
                print ">"+name
                print longest_cut






#longest = (len(pro),
#           m.start(),
#           str(pro),
#           nuc[m.start():m.start()+len(pro)*3+3])







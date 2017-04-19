

#This script adjusts transcript (nucleotide sequence) output from the maker annotation process.
#Maker outputs a transcript file which contains nucleotide sequences for all of the identified
#transcripts. However, in many cases the transcripts have an 'offset' of DNA sequence prior to
#the identified start codon; this information is identified in sequence headers as 'transcript
#offset: #'. If the # is zero, there is no offset and the sequence can be used as-is, and is in-frame
#If the number is something other than zero, the leading offset must be deleted to obtain the actual CDS.
#
#This script reads the transcript fasta headers and identifies whether there is an offset.
#If there is no offset ('transcript offset:0'), it prints the shortened transcript header (trimmed
#after the first space) and corresponding full sequence. If there is an offset ('transcript offset:5',
#for example), the leading offset sequence is trimmed off, and the shortened header and remaining
#sequence are printed. There is also a minimum transcript length that can be set, corresponding to
#the final sequences.
#
#Syntax is:
# python K_L.codon_distances.py maker_transcripts.fasta min_length_cutoff > cleaned_maker_transcripts.fasta






from __future__ import division

import re
import sys
import math
from Bio import SeqIO


min_len=sys.argv[2]


for record in SeqIO.parse(sys.argv[1],"fasta"):
    short_name=record.id
    name=record.description
    sequence=record.seq
    if len(record.seq)>int(min_len)-1:
        if "transcript offset:0" in name:
            print ">"+short_name
            print sequence
        else:
            offset_list=re.split(r'[ :]', name)
            offset=offset_list[3]
            clean_sequence=sequence[int(offset):]
            if len(clean_sequence)>int(min_len)-1:
                print ">"+short_name
                print clean_sequence

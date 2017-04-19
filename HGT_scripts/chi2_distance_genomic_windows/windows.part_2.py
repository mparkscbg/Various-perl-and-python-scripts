
#this script calculates the tetranucleotide proportions in sliding windows for all sequences in a fasta file, and outputs them as comma-delimited lists starting with fasta header_windowstart_windowstop,
#and followed by proportions of each tetranucleotide in the window
#
#syntax:
#
#python windows.part_2.py input.fasta window_size step_size min_#_tetranucleotides_in_window > output.txt
#
#NOTE 1: the min_num_tetranucleotides_in_window is to minimize the impact of hard-masked sequences (N'd out) taking up significant portions of a window. In a 5kbp window, there are 4997 tetranucleotides.
#If min_#_tetranucleotides_in_window is set to 3998, for example, this would mean that no more than ~20% of any given window is masked (i.e., 80% of tetranucleotides do NOT have any 'N' characters in
#their sequence.
#
#NOTE 2: fig values/calculations are as defined in Becq J, Churlaud C, Deschavanne P. 2010. A benchmark of parametric methods for horizontal transfers detection. PLoSONE 5(4):e9989.
#
#NOTE 3: this script is definitely slower than windows.part_1.py. I have run this on a ca. 90Mbp genome on a single processor with a run
#time of a day or so.

from __future__ import division

import re
import sys
import math
from Bio import SeqIO

win=sys.argv[2]
step=sys.argv[3]
occupancy=sys.argv[4]



for record in SeqIO.parse(sys.argv[1],"fasta"):
    name=record.id
    sequence=record.seq
    seqlen=len(sequence)
    #print name
    #print sequence
    #print seqlen
	
    for i in range (0,int(seqlen)-int(win)+1,int(step)):        #set up sliding window parameters
        
        subsequence=sequence[i:i+int(win)]                      #isolate sequence in window
        subseq_name=name+"_"+str(i+1)+"_"+str(i+int(win))       #name subsequence from original fasta header<tab>starting position<tab>end position
        #print subseq_name
        #initialize tetranucleotide dictionary for window
        tetra_dict =   {'AAAA': 0, 'AAAC': 0, 'AAAG': 0, 'AAAT': 0, 'AACA': 0, 'AACC': 0, 'AACG': 0, 'AACT': 0, 'AAGA': 0, 'AAGC': 0, 'AAGG': 0, 'AAGT': 0, 'AATA': 0, 'AATC': 0, 'AATG': 0, 'AATT': 0,
                        'ACAA': 0, 'ACAC': 0, 'ACAG': 0, 'ACAT': 0, 'ACCA': 0, 'ACCC': 0, 'ACCG': 0, 'ACCT': 0, 'ACGA': 0, 'ACGC': 0, 'ACGG': 0, 'ACGT': 0, 'ACTA': 0, 'ACTC': 0, 'ACTG': 0, 'ACTT': 0,
                        'AGAA': 0, 'AGAC': 0, 'AGAG': 0, 'AGAT': 0, 'AGCA': 0, 'AGCC': 0, 'AGCG': 0, 'AGCT': 0, 'AGGA': 0, 'AGGC': 0, 'AGGG': 0, 'AGGT': 0, 'AGTA': 0, 'AGTC': 0, 'AGTG': 0, 'AGTT': 0,
                        'ATAA': 0, 'ATAC': 0, 'ATAG': 0, 'ATAT': 0, 'ATCA': 0, 'ATCC': 0, 'ATCG': 0, 'ATCT': 0, 'ATGA': 0, 'ATGC': 0, 'ATGG': 0, 'ATGT': 0, 'ATTA': 0, 'ATTC': 0, 'ATTG': 0, 'ATTT': 0,
                        'CAAA': 0, 'CAAC': 0, 'CAAG': 0, 'CAAT': 0, 'CACA': 0, 'CACC': 0, 'CACG': 0, 'CACT': 0, 'CAGA': 0, 'CAGC': 0, 'CAGG': 0, 'CAGT': 0, 'CATA': 0, 'CATC': 0, 'CATG': 0, 'CATT': 0,
                        'CCAA': 0, 'CCAC': 0, 'CCAG': 0, 'CCAT': 0, 'CCCA': 0, 'CCCC': 0, 'CCCG': 0, 'CCCT': 0, 'CCGA': 0, 'CCGC': 0, 'CCGG': 0, 'CCGT': 0, 'CCTA': 0, 'CCTC': 0, 'CCTG': 0, 'CCTT': 0,
                        'CGAA': 0, 'CGAC': 0, 'CGAG': 0, 'CGAT': 0, 'CGCA': 0, 'CGCC': 0, 'CGCG': 0, 'CGCT': 0, 'CGGA': 0, 'CGGC': 0, 'CGGG': 0, 'CGGT': 0, 'CGTA': 0, 'CGTC': 0, 'CGTG': 0, 'CGTT': 0,
                        'CTAA': 0, 'CTAC': 0, 'CTAG': 0, 'CTAT': 0, 'CTCA': 0, 'CTCC': 0, 'CTCG': 0, 'CTCT': 0, 'CTGA': 0, 'CTGC': 0, 'CTGG': 0, 'CTGT': 0, 'CTTA': 0, 'CTTC': 0, 'CTTG': 0, 'CTTT': 0,
                        'GAAA': 0, 'GAAC': 0, 'GAAG': 0, 'GAAT': 0, 'GACA': 0, 'GACC': 0, 'GACG': 0, 'GACT': 0, 'GAGA': 0, 'GAGC': 0, 'GAGG': 0, 'GAGT': 0, 'GATA': 0, 'GATC': 0, 'GATG': 0, 'GATT': 0,
                        'GCAA': 0, 'GCAC': 0, 'GCAG': 0, 'GCAT': 0, 'GCCA': 0, 'GCCC': 0, 'GCCG': 0, 'GCCT': 0, 'GCGA': 0, 'GCGC': 0, 'GCGG': 0, 'GCGT': 0, 'GCTA': 0, 'GCTC': 0, 'GCTG': 0, 'GCTT': 0,
                        'GGAA': 0, 'GGAC': 0, 'GGAG': 0, 'GGAT': 0, 'GGCA': 0, 'GGCC': 0, 'GGCG': 0, 'GGCT': 0, 'GGGA': 0, 'GGGC': 0, 'GGGG': 0, 'GGGT': 0, 'GGTA': 0, 'GGTC': 0, 'GGTG': 0, 'GGTT': 0,
                        'GTAA': 0, 'GTAC': 0, 'GTAG': 0, 'GTAT': 0, 'GTCA': 0, 'GTCC': 0, 'GTCG': 0, 'GTCT': 0, 'GTGA': 0, 'GTGC': 0, 'GTGG': 0, 'GTGT': 0, 'GTTA': 0, 'GTTC': 0, 'GTTG': 0, 'GTTT': 0,
                        'TAAA': 0, 'TAAC': 0, 'TAAG': 0, 'TAAT': 0, 'TACA': 0, 'TACC': 0, 'TACG': 0, 'TACT': 0, 'TAGA': 0, 'TAGC': 0, 'TAGG': 0, 'TAGT': 0, 'TATA': 0, 'TATC': 0, 'TATG': 0, 'TATT': 0,
                        'TCAA': 0, 'TCAC': 0, 'TCAG': 0, 'TCAT': 0, 'TCCA': 0, 'TCCC': 0, 'TCCG': 0, 'TCCT': 0, 'TCGA': 0, 'TCGC': 0, 'TCGG': 0, 'TCGT': 0, 'TCTA': 0, 'TCTC': 0, 'TCTG': 0, 'TCTT': 0,
                        'TGAA': 0, 'TGAC': 0, 'TGAG': 0, 'TGAT': 0, 'TGCA': 0, 'TGCC': 0, 'TGCG': 0, 'TGCT': 0, 'TGGA': 0, 'TGGC': 0, 'TGGG': 0, 'TGGT': 0, 'TGTA': 0, 'TGTC': 0, 'TGTG': 0, 'TGTT': 0,
                        'TTAA': 0, 'TTAC': 0, 'TTAG': 0, 'TTAT': 0, 'TTCA': 0, 'TTCC': 0, 'TTCG': 0, 'TTCT': 0, 'TTGA': 0, 'TTGC': 0, 'TTGG': 0, 'TTGT': 0, 'TTTA': 0, 'TTTC': 0, 'TTTG': 0, 'TTTT': 0}


        for i in range (0,int(win)-4+1,1):                      #set up tetra-nucleotide sliding windows within larger window
            tetra=subsequence[i:i+4]
            #print tetra

            if tetra.upper() == 'AAAA':
                tetra_dict['AAAA']+=1

            elif tetra.upper() == 'AAAC':
                tetra_dict['AAAC']+=1

            elif tetra.upper() == 'AAAG':
                tetra_dict['AAAG']+=1

            elif tetra.upper() == 'AAAT':
                tetra_dict['AAAT']+=1

            elif tetra.upper() == 'AACA':
                tetra_dict['AACA']+=1

            elif tetra.upper() == 'AACC':
                tetra_dict['AACC']+=1

            elif tetra.upper() == 'AACG':
                tetra_dict['AACG']+=1

            elif tetra.upper() == 'AACT':
                tetra_dict['AACT']+=1

            elif tetra.upper() == 'AAGA':
                tetra_dict['AAGA']+=1

            elif tetra.upper() == 'AAGC':
                tetra_dict['AAGC']+=1

            elif tetra.upper() == 'AAGG':
                tetra_dict['AAGG']+=1

            elif tetra.upper() == 'AAGT':
                tetra_dict['AAGT']+=1

            elif tetra.upper() == 'AATA':
                tetra_dict['AATA']+=1

            elif tetra.upper() == 'AATC':
                tetra_dict['AATC']+=1

            elif tetra.upper() == 'AATG':
                tetra_dict['AATG']+=1

            elif tetra.upper() == 'AATT':
                tetra_dict['AATT']+=1


            elif tetra.upper() == 'ACAA':
                tetra_dict['ACAA']+=1

            elif tetra.upper() == 'ACAC':
                tetra_dict['ACAC']+=1

            elif tetra.upper() == 'ACAG':
                tetra_dict['ACAG']+=1

            elif tetra.upper() == 'ACAT':
                tetra_dict['ACAT']+=1

            elif tetra.upper() == 'ACCA':
                tetra_dict['ACCA']+=1

            elif tetra.upper() == 'ACCC':
                tetra_dict['ACCC']+=1

            elif tetra.upper() == 'ACCG':
                tetra_dict['ACCG']+=1

            elif tetra.upper() == 'ACCT':
                tetra_dict['ACCT']+=1

            elif tetra.upper() == 'ACGA':
                tetra_dict['ACGA']+=1

            elif tetra.upper() == 'ACGC':
                tetra_dict['ACGC']+=1

            elif tetra.upper() == 'ACGG':
                tetra_dict['ACGG']+=1

            elif tetra.upper() == 'ACGT':
                tetra_dict['ACGT']+=1

            elif tetra.upper() == 'ACTA':
                tetra_dict['ACTA']+=1

            elif tetra.upper() == 'ACTC':
                tetra_dict['ACTC']+=1

            elif tetra.upper() == 'ACTG':
                tetra_dict['ACTG']+=1

            elif tetra.upper() == 'ACTT':
                tetra_dict['ACTT']+=1


            elif tetra.upper() == 'AGAA':
                tetra_dict['AGAA']+=1

            elif tetra.upper() == 'AGAC':
                tetra_dict['AGAC']+=1

            elif tetra.upper() == 'AGAG':
                tetra_dict['AGAG']+=1

            elif tetra.upper() == 'AGAT':
                tetra_dict['AGAT']+=1

            elif tetra.upper() == 'AGCA':
                tetra_dict['AGCA']+=1

            elif tetra.upper() == 'AGCC':
                tetra_dict['AGCC']+=1

            elif tetra.upper() == 'AGCG':
                tetra_dict['AGCG']+=1

            elif tetra.upper() == 'AGCT':
                tetra_dict['AGCT']+=1

            elif tetra.upper() == 'AGGA':
                tetra_dict['AGGA']+=1

            elif tetra.upper() == 'AGGC':
                tetra_dict['AGGC']+=1

            elif tetra.upper() == 'AGGG':
                tetra_dict['AGGG']+=1

            elif tetra.upper() == 'AGGT':
                tetra_dict['AGGT']+=1

            elif tetra.upper() == 'AGTA':
                tetra_dict['AGTA']+=1

            elif tetra.upper() == 'AGTC':
                tetra_dict['AGTC']+=1

            elif tetra.upper() == 'AGTG':
                tetra_dict['AGTG']+=1

            elif tetra.upper() == 'AGTT':
                tetra_dict['AGTT']+=1


            elif tetra.upper() == 'ATAA':
                tetra_dict['ATAA']+=1

            elif tetra.upper() == 'ATAC':
                tetra_dict['ATAC']+=1

            elif tetra.upper() == 'ATAG':
                tetra_dict['ATAG']+=1

            elif tetra.upper() == 'ATAT':
                tetra_dict['ATAT']+=1

            elif tetra.upper() == 'ATCA':
                tetra_dict['ATCA']+=1

            elif tetra.upper() == 'ATCC':
                tetra_dict['ATCC']+=1

            elif tetra.upper() == 'ATCG':
                tetra_dict['ATCG']+=1

            elif tetra.upper() == 'ATCT':
                tetra_dict['ATCT']+=1

            elif tetra.upper() == 'ATGA':
                tetra_dict['ATGA']+=1

            elif tetra.upper() == 'ATGC':
                tetra_dict['ATGC']+=1

            elif tetra.upper() == 'ATGG':
                tetra_dict['ATGG']+=1

            elif tetra.upper() == 'ATGT':
                tetra_dict['ATGT']+=1

            elif tetra.upper() == 'ATTA':
                tetra_dict['ATTA']+=1

            elif tetra.upper() == 'ATTC':
                tetra_dict['ATTC']+=1

            elif tetra.upper() == 'ATTG':
                tetra_dict['ATTG']+=1

            elif tetra.upper() == 'ATTT':
                tetra_dict['ATTT']+=1



            elif tetra.upper() == 'CAAA':
                tetra_dict['CAAA']+=1

            elif tetra.upper() == 'CAAC':
                tetra_dict['CAAC']+=1

            elif tetra.upper() == 'CAAG':
                tetra_dict['CAAG']+=1

            elif tetra.upper() == 'CAAT':
                tetra_dict['CAAT']+=1

            elif tetra.upper() == 'CACA':
                tetra_dict['CACA']+=1

            elif tetra.upper() == 'CACC':
                tetra_dict['CACC']+=1

            elif tetra.upper() == 'CACG':
                tetra_dict['CACG']+=1

            elif tetra.upper() == 'CACT':
                tetra_dict['CACT']+=1

            elif tetra.upper() == 'CAGA':
                tetra_dict['CAGA']+=1

            elif tetra.upper() == 'CAGC':
                tetra_dict['CAGC']+=1

            elif tetra.upper() == 'CAGG':
                tetra_dict['CAGG']+=1

            elif tetra.upper() == 'CAGT':
                tetra_dict['CAGT']+=1

            elif tetra.upper() == 'CATA':
                tetra_dict['CATA']+=1

            elif tetra.upper() == 'CATC':
                tetra_dict['CATC']+=1

            elif tetra.upper() == 'CATG':
                tetra_dict['CATG']+=1

            elif tetra.upper() == 'CATT':
                tetra_dict['CATT']+=1


            elif tetra.upper() == 'CCAA':
                tetra_dict['CCAA']+=1

            elif tetra.upper() == 'CCAC':
                tetra_dict['CCAC']+=1

            elif tetra.upper() == 'CCAG':
                tetra_dict['CCAG']+=1

            elif tetra.upper() == 'CCAT':
                tetra_dict['CCAT']+=1

            elif tetra.upper() == 'CCCA':
                tetra_dict['CCCA']+=1

            elif tetra.upper() == 'CCCC':
                tetra_dict['CCCC']+=1

            elif tetra.upper() == 'CCCG':
                tetra_dict['CCCG']+=1

            elif tetra.upper() == 'CCCT':
                tetra_dict['CCCT']+=1

            elif tetra.upper() == 'CCGA':
                tetra_dict['CCGA']+=1

            elif tetra.upper() == 'CCGC':
                tetra_dict['CCGC']+=1

            elif tetra.upper() == 'CCGG':
                tetra_dict['CCGG']+=1

            elif tetra.upper() == 'CCGT':
                tetra_dict['CCGT']+=1

            elif tetra.upper() == 'CCTA':
                tetra_dict['CCTA']+=1

            elif tetra.upper() == 'CCTC':
                tetra_dict['CCTC']+=1

            elif tetra.upper() == 'CCTG':
                tetra_dict['CCTG']+=1

            elif tetra.upper() == 'CCTT':
                tetra_dict['CCTT']+=1


            elif tetra.upper() == 'CGAA':
                tetra_dict['CGAA']+=1

            elif tetra.upper() == 'CGAC':
                tetra_dict['CGAC']+=1

            elif tetra.upper() == 'CGAG':
                tetra_dict['CGAG']+=1

            elif tetra.upper() == 'CGAT':
                tetra_dict['CGAT']+=1

            elif tetra.upper() == 'CGCA':
                tetra_dict['CGCA']+=1

            elif tetra.upper() == 'CGCC':
                tetra_dict['CGCC']+=1

            elif tetra.upper() == 'CGCG':
                tetra_dict['CGCG']+=1

            elif tetra.upper() == 'CGCT':
                tetra_dict['CGCT']+=1

            elif tetra.upper() == 'CGGA':
                tetra_dict['CGGA']+=1

            elif tetra.upper() == 'CGGC':
                tetra_dict['CGGC']+=1

            elif tetra.upper() == 'CGGG':
                tetra_dict['CGGG']+=1

            elif tetra.upper() == 'CGGT':
                tetra_dict['CGGT']+=1

            elif tetra.upper() == 'CGTA':
                tetra_dict['CGTA']+=1

            elif tetra.upper() == 'CGTC':
                tetra_dict['CGTC']+=1

            elif tetra.upper() == 'CGTG':
                tetra_dict['CGTG']+=1

            elif tetra.upper() == 'CGTT':
                tetra_dict['CGTT']+=1


            elif tetra.upper() == 'CTAA':
                tetra_dict['CTAA']+=1

            elif tetra.upper() == 'CTAC':
                tetra_dict['CTAC']+=1

            elif tetra.upper() == 'CTAG':
                tetra_dict['CTAG']+=1

            elif tetra.upper() == 'CTAT':
                tetra_dict['CTAT']+=1

            elif tetra.upper() == 'CTCA':
                tetra_dict['CTCA']+=1

            elif tetra.upper() == 'CTCC':
                tetra_dict['CTCC']+=1

            elif tetra.upper() == 'CTCG':
                tetra_dict['CTCG']+=1

            elif tetra.upper() == 'CTCT':
                tetra_dict['CTCT']+=1

            elif tetra.upper() == 'CTGA':
                tetra_dict['CTGA']+=1

            elif tetra.upper() == 'CTGC':
                tetra_dict['CTGC']+=1

            elif tetra.upper() == 'CTGG':
                tetra_dict['CTGG']+=1

            elif tetra.upper() == 'CTGT':
                tetra_dict['CTGT']+=1

            elif tetra.upper() == 'CTTA':
                tetra_dict['CTTA']+=1

            elif tetra.upper() == 'CTTC':
                tetra_dict['CTTC']+=1

            elif tetra.upper() == 'CTTG':
                tetra_dict['CTTG']+=1

            elif tetra.upper() == 'CTTT':
                tetra_dict['CTTT']+=1




            elif tetra.upper() == 'GAAA':
                tetra_dict['GAAA']+=1

            elif tetra.upper() == 'GAAC':
                tetra_dict['GAAC']+=1

            elif tetra.upper() == 'GAAG':
                tetra_dict['GAAG']+=1

            elif tetra.upper() == 'GAAT':
                tetra_dict['GAAT']+=1

            elif tetra.upper() == 'GACA':
                tetra_dict['GACA']+=1

            elif tetra.upper() == 'GACC':
                tetra_dict['GACC']+=1

            elif tetra.upper() == 'GACG':
                tetra_dict['GACG']+=1

            elif tetra.upper() == 'GACT':
                tetra_dict['GACT']+=1

            elif tetra.upper() == 'GAGA':
                tetra_dict['GAGA']+=1

            elif tetra.upper() == 'GAGC':
                tetra_dict['GAGC']+=1

            elif tetra.upper() == 'GAGG':
                tetra_dict['GAGG']+=1

            elif tetra.upper() == 'GAGT':
                tetra_dict['GAGT']+=1

            elif tetra.upper() == 'GATA':
                tetra_dict['GATA']+=1

            elif tetra.upper() == 'GATC':
                tetra_dict['GATC']+=1

            elif tetra.upper() == 'GATG':
                tetra_dict['GATG']+=1

            elif tetra.upper() == 'GATT':
                tetra_dict['GATT']+=1


            elif tetra.upper() == 'GCAA':
                tetra_dict['GCAA']+=1

            elif tetra.upper() == 'GCAC':
                tetra_dict['GCAC']+=1

            elif tetra.upper() == 'GCAG':
                tetra_dict['GCAG']+=1

            elif tetra.upper() == 'GCAT':
                tetra_dict['GCAT']+=1

            elif tetra.upper() == 'GCCA':
                tetra_dict['GCCA']+=1

            elif tetra.upper() == 'GCCC':
                tetra_dict['GCCC']+=1

            elif tetra.upper() == 'GCCG':
                tetra_dict['GCCG']+=1

            elif tetra.upper() == 'GCCT':
                tetra_dict['GCCT']+=1

            elif tetra.upper() == 'GCGA':
                tetra_dict['GCGA']+=1

            elif tetra.upper() == 'GCGC':
                tetra_dict['GCGC']+=1

            elif tetra.upper() == 'GCGG':
                tetra_dict['GCGG']+=1

            elif tetra.upper() == 'GCGT':
                tetra_dict['GCGT']+=1

            elif tetra.upper() == 'GCTA':
                tetra_dict['GCTA']+=1

            elif tetra.upper() == 'GCTC':
                tetra_dict['GCTC']+=1

            elif tetra.upper() == 'GCTG':
                tetra_dict['GCTG']+=1

            elif tetra.upper() == 'GCTT':
                tetra_dict['GCTT']+=1


            elif tetra.upper() == 'GGAA':
                tetra_dict['GGAA']+=1

            elif tetra.upper() == 'GGAC':
                tetra_dict['GGAC']+=1

            elif tetra.upper() == 'GGAG':
                tetra_dict['GGAG']+=1

            elif tetra.upper() == 'GGAT':
                tetra_dict['GGAT']+=1

            elif tetra.upper() == 'GGCA':
                tetra_dict['GGCA']+=1

            elif tetra.upper() == 'GGCC':
                tetra_dict['GGCC']+=1

            elif tetra.upper() == 'GGCG':
                tetra_dict['GGCG']+=1

            elif tetra.upper() == 'GGCT':
                tetra_dict['GGCT']+=1

            elif tetra.upper() == 'GGGA':
                tetra_dict['GGGA']+=1

            elif tetra.upper() == 'GGGC':
                tetra_dict['GGGC']+=1

            elif tetra.upper() == 'GGGG':
                tetra_dict['GGGG']+=1

            elif tetra.upper() == 'GGGT':
                tetra_dict['GGGT']+=1

            elif tetra.upper() == 'GGTA':
                tetra_dict['GGTA']+=1

            elif tetra.upper() == 'GGTC':
                tetra_dict['GGTC']+=1

            elif tetra.upper() == 'GGTG':
                tetra_dict['GGTG']+=1

            elif tetra.upper() == 'GGTT':
                tetra_dict['GGTT']+=1


            elif tetra.upper() == 'GTAA':
                tetra_dict['GTAA']+=1

            elif tetra.upper() == 'GTAC':
                tetra_dict['GTAC']+=1

            elif tetra.upper() == 'GTAG':
                tetra_dict['GTAG']+=1

            elif tetra.upper() == 'GTAT':
                tetra_dict['GTAT']+=1

            elif tetra.upper() == 'GTCA':
                tetra_dict['GTCA']+=1

            elif tetra.upper() == 'GTCC':
                tetra_dict['GTCC']+=1

            elif tetra.upper() == 'GTCG':
                tetra_dict['GTCG']+=1

            elif tetra.upper() == 'GTCT':
                tetra_dict['GTCT']+=1

            elif tetra.upper() == 'GTGA':
                tetra_dict['GTGA']+=1

            elif tetra.upper() == 'GTGC':
                tetra_dict['GTGC']+=1

            elif tetra.upper() == 'GTGG':
                tetra_dict['GTGG']+=1

            elif tetra.upper() == 'GTGT':
                tetra_dict['GTGT']+=1

            elif tetra.upper() == 'GTTA':
                tetra_dict['GTTA']+=1

            elif tetra.upper() == 'GTTC':
                tetra_dict['GTTC']+=1

            elif tetra.upper() == 'GTTG':
                tetra_dict['GTTG']+=1

            elif tetra.upper() == 'GTTT':
                tetra_dict['GTTT']+=1





            elif tetra.upper() == 'TAAA':
                tetra_dict['TAAA']+=1

            elif tetra.upper() == 'TAAC':
                tetra_dict['TAAC']+=1

            elif tetra.upper() == 'TAAG':
                tetra_dict['TAAG']+=1

            elif tetra.upper() == 'TAAT':
                tetra_dict['TAAT']+=1

            elif tetra.upper() == 'TACA':
                tetra_dict['TACA']+=1

            elif tetra.upper() == 'TACC':
                tetra_dict['TACC']+=1

            elif tetra.upper() == 'TACG':
                tetra_dict['TACG']+=1

            elif tetra.upper() == 'TACT':
                tetra_dict['TACT']+=1

            elif tetra.upper() == 'TAGA':
                tetra_dict['TAGA']+=1

            elif tetra.upper() == 'TAGC':
                tetra_dict['TAGC']+=1

            elif tetra.upper() == 'TAGG':
                tetra_dict['TAGG']+=1

            elif tetra.upper() == 'TAGT':
                tetra_dict['TAGT']+=1

            elif tetra.upper() == 'TATA':
                tetra_dict['TATA']+=1

            elif tetra.upper() == 'TATC':
                tetra_dict['TATC']+=1

            elif tetra.upper() == 'TATG':
                tetra_dict['TATG']+=1

            elif tetra.upper() == 'TATT':
                tetra_dict['TATT']+=1


            elif tetra.upper() == 'TCAA':
                tetra_dict['TCAA']+=1

            elif tetra.upper() == 'TCAC':
                tetra_dict['TCAC']+=1

            elif tetra.upper() == 'TCAG':
                tetra_dict['TCAG']+=1

            elif tetra.upper() == 'TCAT':
                tetra_dict['TCAT']+=1

            elif tetra.upper() == 'TCCA':
                tetra_dict['TCCA']+=1

            elif tetra.upper() == 'TCCC':
                tetra_dict['TCCC']+=1

            elif tetra.upper() == 'TCCG':
                tetra_dict['TCCG']+=1

            elif tetra.upper() == 'TCCT':
                tetra_dict['TCCT']+=1

            elif tetra.upper() == 'TCGA':
                tetra_dict['TCGA']+=1

            elif tetra.upper() == 'TCGC':
                tetra_dict['TCGC']+=1

            elif tetra.upper() == 'TCGG':
                tetra_dict['TCGG']+=1

            elif tetra.upper() == 'TCGT':
                tetra_dict['TCGT']+=1

            elif tetra.upper() == 'TCTA':
                tetra_dict['TCTA']+=1

            elif tetra.upper() == 'TCTC':
                tetra_dict['TCTC']+=1

            elif tetra.upper() == 'TCTG':
                tetra_dict['TCTG']+=1

            elif tetra.upper() == 'TCTT':
                tetra_dict['TCTT']+=1


            elif tetra.upper() == 'TGAA':
                tetra_dict['TGAA']+=1

            elif tetra.upper() == 'TGAC':
                tetra_dict['TGAC']+=1

            elif tetra.upper() == 'TGAG':
                tetra_dict['TGAG']+=1

            elif tetra.upper() == 'TGAT':
                tetra_dict['TGAT']+=1

            elif tetra.upper() == 'TGCA':
                tetra_dict['TGCA']+=1

            elif tetra.upper() == 'TGCC':
                tetra_dict['TGCC']+=1

            elif tetra.upper() == 'TGCG':
                tetra_dict['TGCG']+=1

            elif tetra.upper() == 'TGCT':
                tetra_dict['TGCT']+=1

            elif tetra.upper() == 'TGGA':
                tetra_dict['TGGA']+=1

            elif tetra.upper() == 'TGGC':
                tetra_dict['TGGC']+=1

            elif tetra.upper() == 'TGGG':
                tetra_dict['TGGG']+=1

            elif tetra.upper() == 'TGGT':
                tetra_dict['TGGT']+=1

            elif tetra.upper() == 'TGTA':
                tetra_dict['TGTA']+=1

            elif tetra.upper() == 'TGTC':
                tetra_dict['TGTC']+=1

            elif tetra.upper() == 'TGTG':
                tetra_dict['TGTG']+=1

            elif tetra.upper() == 'TGTT':
                tetra_dict['TGTT']+=1


            elif tetra.upper() == 'TTAA':
                tetra_dict['TTAA']+=1

            elif tetra.upper() == 'TTAC':
                tetra_dict['TTAC']+=1

            elif tetra.upper() == 'TTAG':
                tetra_dict['TTAG']+=1

            elif tetra.upper() == 'TTAT':
                tetra_dict['TTAT']+=1

            elif tetra.upper() == 'TTCA':
                tetra_dict['TTCA']+=1

            elif tetra.upper() == 'TTCC':
                tetra_dict['TTCC']+=1

            elif tetra.upper() == 'TTCG':
                tetra_dict['TTCG']+=1

            elif tetra.upper() == 'TTCT':
                tetra_dict['TTCT']+=1

            elif tetra.upper() == 'TTGA':
                tetra_dict['TTGA']+=1

            elif tetra.upper() == 'TTGC':
                tetra_dict['TTGC']+=1

            elif tetra.upper() == 'TTGG':
                tetra_dict['TTGG']+=1

            elif tetra.upper() == 'TTGT':
                tetra_dict['TTGT']+=1

            elif tetra.upper() == 'TTTA':
                tetra_dict['TTTA']+=1

            elif tetra.upper() == 'TTTC':
                tetra_dict['TTTC']+=1

            elif tetra.upper() == 'TTTG':
                tetra_dict['TTTG']+=1

            elif tetra.upper() == 'TTTT':
                tetra_dict['TTTT']+=1



        total_tetras=sum(tetra_dict.values())               #specify value of total_tetras as sum of all tetranucleotides measured in all windows

        if int(total_tetras)>int(occupancy)-1:
            fig_AAAA=round(tetra_dict['AAAA']/total_tetras,5)
            fig_AAAC=round(tetra_dict['AAAC']/total_tetras,5)
            fig_AAAG=round(tetra_dict['AAAG']/total_tetras,5)
            fig_AAAT=round(tetra_dict['AAAT']/total_tetras,5)
            fig_AACA=round(tetra_dict['AACA']/total_tetras,5)
            fig_AACC=round(tetra_dict['AACC']/total_tetras,5)
            fig_AACG=round(tetra_dict['AACG']/total_tetras,5)
            fig_AACT=round(tetra_dict['AACT']/total_tetras,5)
            fig_AAGA=round(tetra_dict['AAGA']/total_tetras,5)
            fig_AAGC=round(tetra_dict['AAGC']/total_tetras,5)
            fig_AAGG=round(tetra_dict['AAGG']/total_tetras,5)
            fig_AAGT=round(tetra_dict['AAGT']/total_tetras,5)
            fig_AATA=round(tetra_dict['AATA']/total_tetras,5)
            fig_AATC=round(tetra_dict['AATC']/total_tetras,5)
            fig_AATG=round(tetra_dict['AATG']/total_tetras,5)
            fig_AATT=round(tetra_dict['AATT']/total_tetras,5)

            fig_ACAA=round(tetra_dict['ACAA']/total_tetras,5)
            fig_ACAC=round(tetra_dict['ACAC']/total_tetras,5)
            fig_ACAG=round(tetra_dict['ACAG']/total_tetras,5)
            fig_ACAT=round(tetra_dict['ACAT']/total_tetras,5)
            fig_ACCA=round(tetra_dict['ACCA']/total_tetras,5)
            fig_ACCC=round(tetra_dict['ACCC']/total_tetras,5)
            fig_ACCG=round(tetra_dict['ACCG']/total_tetras,5)
            fig_ACCT=round(tetra_dict['ACCT']/total_tetras,5)
            fig_ACGA=round(tetra_dict['ACGA']/total_tetras,5)
            fig_ACGC=round(tetra_dict['ACGC']/total_tetras,5)
            fig_ACGG=round(tetra_dict['ACGG']/total_tetras,5)
            fig_ACGT=round(tetra_dict['ACGT']/total_tetras,5)
            fig_ACTA=round(tetra_dict['ACTA']/total_tetras,5)
            fig_ACTC=round(tetra_dict['ACTC']/total_tetras,5)
            fig_ACTG=round(tetra_dict['ACTG']/total_tetras,5)
            fig_ACTT=round(tetra_dict['ACTT']/total_tetras,5)

            fig_AGAA=round(tetra_dict['AGAA']/total_tetras,5)
            fig_AGAC=round(tetra_dict['AGAC']/total_tetras,5)
            fig_AGAG=round(tetra_dict['AGAG']/total_tetras,5)
            fig_AGAT=round(tetra_dict['AGAT']/total_tetras,5)
            fig_AGCA=round(tetra_dict['AGCA']/total_tetras,5)
            fig_AGCC=round(tetra_dict['AGCC']/total_tetras,5)
            fig_AGCG=round(tetra_dict['AGCG']/total_tetras,5)
            fig_AGCT=round(tetra_dict['AGCT']/total_tetras,5)
            fig_AGGA=round(tetra_dict['AGGA']/total_tetras,5)
            fig_AGGC=round(tetra_dict['AGGC']/total_tetras,5)
            fig_AGGG=round(tetra_dict['AGGG']/total_tetras,5)
            fig_AGGT=round(tetra_dict['AGGT']/total_tetras,5)
            fig_AGTA=round(tetra_dict['AGTA']/total_tetras,5)
            fig_AGTC=round(tetra_dict['AGTC']/total_tetras,5)
            fig_AGTG=round(tetra_dict['AGTG']/total_tetras,5)
            fig_AGTT=round(tetra_dict['AGTT']/total_tetras,5)

            fig_ATAA=round(tetra_dict['ATAA']/total_tetras,5)
            fig_ATAC=round(tetra_dict['ATAC']/total_tetras,5)
            fig_ATAG=round(tetra_dict['ATAG']/total_tetras,5)
            fig_ATAT=round(tetra_dict['ATAT']/total_tetras,5)
            fig_ATCA=round(tetra_dict['ATCA']/total_tetras,5)
            fig_ATCC=round(tetra_dict['ATCC']/total_tetras,5)
            fig_ATCG=round(tetra_dict['ATCG']/total_tetras,5)
            fig_ATCT=round(tetra_dict['ATCT']/total_tetras,5)
            fig_ATGA=round(tetra_dict['ATGA']/total_tetras,5)
            fig_ATGC=round(tetra_dict['ATGC']/total_tetras,5)
            fig_ATGG=round(tetra_dict['ATGG']/total_tetras,5)
            fig_ATGT=round(tetra_dict['ATGT']/total_tetras,5)
            fig_ATTA=round(tetra_dict['ATTA']/total_tetras,5)
            fig_ATTC=round(tetra_dict['ATTC']/total_tetras,5)
            fig_ATTG=round(tetra_dict['ATTG']/total_tetras,5)
            fig_ATTT=round(tetra_dict['ATTT']/total_tetras,5)

            fig_CAAA=round(tetra_dict['CAAA']/total_tetras,5)
            fig_CAAC=round(tetra_dict['CAAC']/total_tetras,5)
            fig_CAAG=round(tetra_dict['CAAG']/total_tetras,5)
            fig_CAAT=round(tetra_dict['CAAT']/total_tetras,5)
            fig_CACA=round(tetra_dict['CACA']/total_tetras,5)
            fig_CACC=round(tetra_dict['CACC']/total_tetras,5)
            fig_CACG=round(tetra_dict['CACG']/total_tetras,5)
            fig_CACT=round(tetra_dict['CACT']/total_tetras,5)
            fig_CAGA=round(tetra_dict['CAGA']/total_tetras,5)
            fig_CAGC=round(tetra_dict['CAGC']/total_tetras,5)
            fig_CAGG=round(tetra_dict['CAGG']/total_tetras,5)
            fig_CAGT=round(tetra_dict['CAGT']/total_tetras,5)
            fig_CATA=round(tetra_dict['CATA']/total_tetras,5)
            fig_CATC=round(tetra_dict['CATC']/total_tetras,5)
            fig_CATG=round(tetra_dict['CATG']/total_tetras,5)
            fig_CATT=round(tetra_dict['CATT']/total_tetras,5)

            fig_CCAA=round(tetra_dict['CCAA']/total_tetras,5)
            fig_CCAC=round(tetra_dict['CCAC']/total_tetras,5)
            fig_CCAG=round(tetra_dict['CCAG']/total_tetras,5)
            fig_CCAT=round(tetra_dict['CCAT']/total_tetras,5)
            fig_CCCA=round(tetra_dict['CCCA']/total_tetras,5)
            fig_CCCC=round(tetra_dict['CCCC']/total_tetras,5)
            fig_CCCG=round(tetra_dict['CCCG']/total_tetras,5)
            fig_CCCT=round(tetra_dict['CCCT']/total_tetras,5)
            fig_CCGA=round(tetra_dict['CCGA']/total_tetras,5)
            fig_CCGC=round(tetra_dict['CCGC']/total_tetras,5)
            fig_CCGG=round(tetra_dict['CCGG']/total_tetras,5)
            fig_CCGT=round(tetra_dict['CCGT']/total_tetras,5)
            fig_CCTA=round(tetra_dict['CCTA']/total_tetras,5)
            fig_CCTC=round(tetra_dict['CCTC']/total_tetras,5)
            fig_CCTG=round(tetra_dict['CCTG']/total_tetras,5)
            fig_CCTT=round(tetra_dict['CCTT']/total_tetras,5)

            fig_CGAA=round(tetra_dict['CGAA']/total_tetras,5)
            fig_CGAC=round(tetra_dict['CGAC']/total_tetras,5)
            fig_CGAG=round(tetra_dict['CGAG']/total_tetras,5)
            fig_CGAT=round(tetra_dict['CGAT']/total_tetras,5)
            fig_CGCA=round(tetra_dict['CGCA']/total_tetras,5)
            fig_CGCC=round(tetra_dict['CGCC']/total_tetras,5)
            fig_CGCG=round(tetra_dict['CGCG']/total_tetras,5)
            fig_CGCT=round(tetra_dict['CGCT']/total_tetras,5)
            fig_CGGA=round(tetra_dict['CGGA']/total_tetras,5)
            fig_CGGC=round(tetra_dict['CGGC']/total_tetras,5)
            fig_CGGG=round(tetra_dict['CGGG']/total_tetras,5)
            fig_CGGT=round(tetra_dict['CGGT']/total_tetras,5)
            fig_CGTA=round(tetra_dict['CGTA']/total_tetras,5)
            fig_CGTC=round(tetra_dict['CGTC']/total_tetras,5)
            fig_CGTG=round(tetra_dict['CGTG']/total_tetras,5)
            fig_CGTT=round(tetra_dict['CGTT']/total_tetras,5)

            fig_CTAA=round(tetra_dict['CTAA']/total_tetras,5)
            fig_CTAC=round(tetra_dict['CTAC']/total_tetras,5)
            fig_CTAG=round(tetra_dict['CTAG']/total_tetras,5)
            fig_CTAT=round(tetra_dict['CTAT']/total_tetras,5)
            fig_CTCA=round(tetra_dict['CTCA']/total_tetras,5)
            fig_CTCC=round(tetra_dict['CTCC']/total_tetras,5)
            fig_CTCG=round(tetra_dict['CTCG']/total_tetras,5)
            fig_CTCT=round(tetra_dict['CTCT']/total_tetras,5)
            fig_CTGA=round(tetra_dict['CTGA']/total_tetras,5)
            fig_CTGC=round(tetra_dict['CTGC']/total_tetras,5)
            fig_CTGG=round(tetra_dict['CTGG']/total_tetras,5)
            fig_CTGT=round(tetra_dict['CTGT']/total_tetras,5)
            fig_CTTA=round(tetra_dict['CTTA']/total_tetras,5)
            fig_CTTC=round(tetra_dict['CTTC']/total_tetras,5)
            fig_CTTG=round(tetra_dict['CTTG']/total_tetras,5)
            fig_CTTT=round(tetra_dict['CTTT']/total_tetras,5)

            fig_GAAA=round(tetra_dict['GAAA']/total_tetras,5)
            fig_GAAC=round(tetra_dict['GAAC']/total_tetras,5)
            fig_GAAG=round(tetra_dict['GAAG']/total_tetras,5)
            fig_GAAT=round(tetra_dict['GAAT']/total_tetras,5)
            fig_GACA=round(tetra_dict['GACA']/total_tetras,5)
            fig_GACC=round(tetra_dict['GACC']/total_tetras,5)
            fig_GACG=round(tetra_dict['GACG']/total_tetras,5)
            fig_GACT=round(tetra_dict['GACT']/total_tetras,5)
            fig_GAGA=round(tetra_dict['GAGA']/total_tetras,5)
            fig_GAGC=round(tetra_dict['GAGC']/total_tetras,5)
            fig_GAGG=round(tetra_dict['GAGG']/total_tetras,5)
            fig_GAGT=round(tetra_dict['GAGT']/total_tetras,5)
            fig_GATA=round(tetra_dict['GATA']/total_tetras,5)
            fig_GATC=round(tetra_dict['GATC']/total_tetras,5)
            fig_GATG=round(tetra_dict['GATG']/total_tetras,5)
            fig_GATT=round(tetra_dict['GATT']/total_tetras,5)

            fig_GCAA=round(tetra_dict['GCAA']/total_tetras,5)
            fig_GCAC=round(tetra_dict['GCAC']/total_tetras,5)
            fig_GCAG=round(tetra_dict['GCAG']/total_tetras,5)
            fig_GCAT=round(tetra_dict['GCAT']/total_tetras,5)
            fig_GCCA=round(tetra_dict['GCCA']/total_tetras,5)
            fig_GCCC=round(tetra_dict['GCCC']/total_tetras,5)
            fig_GCCG=round(tetra_dict['GCCG']/total_tetras,5)
            fig_GCCT=round(tetra_dict['GCCT']/total_tetras,5)
            fig_GCGA=round(tetra_dict['GCGA']/total_tetras,5)
            fig_GCGC=round(tetra_dict['GCGC']/total_tetras,5)
            fig_GCGG=round(tetra_dict['GCGG']/total_tetras,5)
            fig_GCGT=round(tetra_dict['GCGT']/total_tetras,5)
            fig_GCTA=round(tetra_dict['GCTA']/total_tetras,5)
            fig_GCTC=round(tetra_dict['GCTC']/total_tetras,5)
            fig_GCTG=round(tetra_dict['GCTG']/total_tetras,5)
            fig_GCTT=round(tetra_dict['GCTT']/total_tetras,5)

            fig_GGAA=round(tetra_dict['GGAA']/total_tetras,5)
            fig_GGAC=round(tetra_dict['GGAC']/total_tetras,5)
            fig_GGAG=round(tetra_dict['GGAG']/total_tetras,5)
            fig_GGAT=round(tetra_dict['GGAT']/total_tetras,5)
            fig_GGCA=round(tetra_dict['GGCA']/total_tetras,5)
            fig_GGCC=round(tetra_dict['GGCC']/total_tetras,5)
            fig_GGCG=round(tetra_dict['GGCG']/total_tetras,5)
            fig_GGCT=round(tetra_dict['GGCT']/total_tetras,5)
            fig_GGGA=round(tetra_dict['GGGA']/total_tetras,5)
            fig_GGGC=round(tetra_dict['GGGC']/total_tetras,5)
            fig_GGGG=round(tetra_dict['GGGG']/total_tetras,5)
            fig_GGGT=round(tetra_dict['GGGT']/total_tetras,5)
            fig_GGTA=round(tetra_dict['GGTA']/total_tetras,5)
            fig_GGTC=round(tetra_dict['GGTC']/total_tetras,5)
            fig_GGTG=round(tetra_dict['GGTG']/total_tetras,5)
            fig_GGTT=round(tetra_dict['GGTT']/total_tetras,5)

            fig_GTAA=round(tetra_dict['GTAA']/total_tetras,5)
            fig_GTAC=round(tetra_dict['GTAC']/total_tetras,5)
            fig_GTAG=round(tetra_dict['GTAG']/total_tetras,5)
            fig_GTAT=round(tetra_dict['GTAT']/total_tetras,5)
            fig_GTCA=round(tetra_dict['GTCA']/total_tetras,5)
            fig_GTCC=round(tetra_dict['GTCC']/total_tetras,5)
            fig_GTCG=round(tetra_dict['GTCG']/total_tetras,5)
            fig_GTCT=round(tetra_dict['GTCT']/total_tetras,5)
            fig_GTGA=round(tetra_dict['GTGA']/total_tetras,5)
            fig_GTGC=round(tetra_dict['GTGC']/total_tetras,5)
            fig_GTGG=round(tetra_dict['GTGG']/total_tetras,5)
            fig_GTGT=round(tetra_dict['GTGT']/total_tetras,5)
            fig_GTTA=round(tetra_dict['GTTA']/total_tetras,5)
            fig_GTTC=round(tetra_dict['GTTC']/total_tetras,5)
            fig_GTTG=round(tetra_dict['GTTG']/total_tetras,5)
            fig_GTTT=round(tetra_dict['GTTT']/total_tetras,5)

            fig_TAAA=round(tetra_dict['TAAA']/total_tetras,5)
            fig_TAAC=round(tetra_dict['TAAC']/total_tetras,5)
            fig_TAAG=round(tetra_dict['TAAG']/total_tetras,5)
            fig_TAAT=round(tetra_dict['TAAT']/total_tetras,5)
            fig_TACA=round(tetra_dict['TACA']/total_tetras,5)
            fig_TACC=round(tetra_dict['TACC']/total_tetras,5)
            fig_TACG=round(tetra_dict['TACG']/total_tetras,5)
            fig_TACT=round(tetra_dict['TACT']/total_tetras,5)
            fig_TAGA=round(tetra_dict['TAGA']/total_tetras,5)
            fig_TAGC=round(tetra_dict['TAGC']/total_tetras,5)
            fig_TAGG=round(tetra_dict['TAGG']/total_tetras,5)
            fig_TAGT=round(tetra_dict['TAGT']/total_tetras,5)
            fig_TATA=round(tetra_dict['TATA']/total_tetras,5)
            fig_TATC=round(tetra_dict['TATC']/total_tetras,5)
            fig_TATG=round(tetra_dict['TATG']/total_tetras,5)
            fig_TATT=round(tetra_dict['TATT']/total_tetras,5)

            fig_TCAA=round(tetra_dict['TCAA']/total_tetras,5)
            fig_TCAC=round(tetra_dict['TCAC']/total_tetras,5)
            fig_TCAG=round(tetra_dict['TCAG']/total_tetras,5)
            fig_TCAT=round(tetra_dict['TCAT']/total_tetras,5)
            fig_TCCA=round(tetra_dict['TCCA']/total_tetras,5)
            fig_TCCC=round(tetra_dict['TCCC']/total_tetras,5)
            fig_TCCG=round(tetra_dict['TCCG']/total_tetras,5)
            fig_TCCT=round(tetra_dict['TCCT']/total_tetras,5)
            fig_TCGA=round(tetra_dict['TCGA']/total_tetras,5)
            fig_TCGC=round(tetra_dict['TCGC']/total_tetras,5)
            fig_TCGG=round(tetra_dict['TCGG']/total_tetras,5)
            fig_TCGT=round(tetra_dict['TCGT']/total_tetras,5)
            fig_TCTA=round(tetra_dict['TCTA']/total_tetras,5)
            fig_TCTC=round(tetra_dict['TCTC']/total_tetras,5)
            fig_TCTG=round(tetra_dict['TCTG']/total_tetras,5)
            fig_TCTT=round(tetra_dict['TCTT']/total_tetras,5)

            fig_TGAA=round(tetra_dict['TGAA']/total_tetras,5)
            fig_TGAC=round(tetra_dict['TGAC']/total_tetras,5)
            fig_TGAG=round(tetra_dict['TGAG']/total_tetras,5)
            fig_TGAT=round(tetra_dict['TGAT']/total_tetras,5)
            fig_TGCA=round(tetra_dict['TGCA']/total_tetras,5)
            fig_TGCC=round(tetra_dict['TGCC']/total_tetras,5)
            fig_TGCG=round(tetra_dict['TGCG']/total_tetras,5)
            fig_TGCT=round(tetra_dict['TGCT']/total_tetras,5)
            fig_TGGA=round(tetra_dict['TGGA']/total_tetras,5)
            fig_TGGC=round(tetra_dict['TGGC']/total_tetras,5)
            fig_TGGG=round(tetra_dict['TGGG']/total_tetras,5)
            fig_TGGT=round(tetra_dict['TGGT']/total_tetras,5)
            fig_TGTA=round(tetra_dict['TGTA']/total_tetras,5)
            fig_TGTC=round(tetra_dict['TGTC']/total_tetras,5)
            fig_TGTG=round(tetra_dict['TGTG']/total_tetras,5)
            fig_TGTT=round(tetra_dict['TGTT']/total_tetras,5)

            fig_TTAA=round(tetra_dict['TTAA']/total_tetras,5)
            fig_TTAC=round(tetra_dict['TTAC']/total_tetras,5)
            fig_TTAG=round(tetra_dict['TTAG']/total_tetras,5)
            fig_TTAT=round(tetra_dict['TTAT']/total_tetras,5)
            fig_TTCA=round(tetra_dict['TTCA']/total_tetras,5)
            fig_TTCC=round(tetra_dict['TTCC']/total_tetras,5)
            fig_TTCG=round(tetra_dict['TTCG']/total_tetras,5)
            fig_TTCT=round(tetra_dict['TTCT']/total_tetras,5)
            fig_TTGA=round(tetra_dict['TTGA']/total_tetras,5)
            fig_TTGC=round(tetra_dict['TTGC']/total_tetras,5)
            fig_TTGG=round(tetra_dict['TTGG']/total_tetras,5)
            fig_TTGT=round(tetra_dict['TTGT']/total_tetras,5)
            fig_TTTA=round(tetra_dict['TTTA']/total_tetras,5)
            fig_TTTC=round(tetra_dict['TTTC']/total_tetras,5)
            fig_TTTG=round(tetra_dict['TTTG']/total_tetras,5)
            fig_TTTT=round(tetra_dict['TTTT']/total_tetras,5)


#build fig_list from fig calculations
            fig_list=[fig_AAAA,fig_AAAC,fig_AAAG,fig_AAAT,fig_AACA,fig_AACC,fig_AACG,fig_AACT,fig_AAGA,fig_AAGC,fig_AAGG,fig_AAGT,fig_AATA,fig_AATC,fig_AATG,fig_AATT,fig_ACAA,fig_ACAC,fig_ACAG,fig_ACAT,fig_ACCA,fig_ACCC,fig_ACCG,fig_ACCT,fig_ACGA,fig_ACGC,fig_ACGG,fig_ACGT,fig_ACTA,fig_ACTC,fig_ACTG,fig_ACTT,fig_AGAA,fig_AGAC,fig_AGAG,fig_AGAT,fig_AGCA,fig_AGCC,fig_AGCG,fig_AGCT,fig_AGGA,fig_AGGC,fig_AGGG,fig_AGGT,fig_AGTA,fig_AGTC,fig_AGTG,fig_AGTT,fig_ATAA,fig_ATAC,fig_ATAG,fig_ATAT,fig_ATCA,fig_ATCC,fig_ATCG,fig_ATCT,fig_ATGA,fig_ATGC,fig_ATGG,fig_ATGT,fig_ATTA,fig_ATTC,fig_ATTG,fig_ATTT,fig_CAAA,fig_CAAC,fig_CAAG,fig_CAAT,fig_CACA,fig_CACC,fig_CACG,fig_CACT,fig_CAGA,fig_CAGC,fig_CAGG,fig_CAGT,fig_CATA,fig_CATC,fig_CATG,fig_CATT,fig_CCAA,fig_CCAC,fig_CCAG,fig_CCAT,fig_CCCA,fig_CCCC,fig_CCCG,fig_CCCT,fig_CCGA,fig_CCGC,fig_CCGG,fig_CCGT,fig_CCTA,fig_CCTC,fig_CCTG,fig_CCTT,fig_CGAA,fig_CGAC,fig_CGAG,fig_CGAT,fig_CGCA,fig_CGCC,fig_CGCG,fig_CGCT,fig_CGGA,fig_CGGC,fig_CGGG,fig_CGGT,fig_CGTA,fig_CGTC,fig_CGTG,fig_CGTT,fig_CTAA,fig_CTAC,fig_CTAG,fig_CTAT,fig_CTCA,fig_CTCC,fig_CTCG,fig_CTCT,fig_CTGA,fig_CTGC,fig_CTGG,fig_CTGT,fig_CTTA,fig_CTTC,fig_CTTG,fig_CTTT,fig_GAAA,fig_GAAC,fig_GAAG,fig_GAAT,fig_GACA,fig_GACC,fig_GACG,fig_GACT,fig_GAGA,fig_GAGC,fig_GAGG,fig_GAGT,fig_GATA,fig_GATC,fig_GATG,fig_GATT,fig_GCAA,fig_GCAC,fig_GCAG,fig_GCAT,fig_GCCA,fig_GCCC,fig_GCCG,fig_GCCT,fig_GCGA,fig_GCGC,fig_GCGG,fig_GCGT,fig_GCTA,fig_GCTC,fig_GCTG,fig_GCTT,fig_GGAA,fig_GGAC,fig_GGAG,fig_GGAT,fig_GGCA,fig_GGCC,fig_GGCG,fig_GGCT,fig_GGGA,fig_GGGC,fig_GGGG,fig_GGGT,fig_GGTA,fig_GGTC,fig_GGTG,fig_GGTT,fig_GTAA,fig_GTAC,fig_GTAG,fig_GTAT,fig_GTCA,fig_GTCC,fig_GTCG,fig_GTCT,fig_GTGA,fig_GTGC,fig_GTGG,fig_GTGT,fig_GTTA,fig_GTTC,fig_GTTG,fig_GTTT,fig_TAAA,fig_TAAC,fig_TAAG,fig_TAAT,fig_TACA,fig_TACC,fig_TACG,fig_TACT,fig_TAGA,fig_TAGC,fig_TAGG,fig_TAGT,fig_TATA,fig_TATC,fig_TATG,fig_TATT,fig_TCAA,fig_TCAC,fig_TCAG,fig_TCAT,fig_TCCA,fig_TCCC,fig_TCCG,fig_TCCT,fig_TCGA,fig_TCGC,fig_TCGG,fig_TCGT,fig_TCTA,fig_TCTC,fig_TCTG,fig_TCTT,fig_TGAA,fig_TGAC,fig_TGAG,fig_TGAT,fig_TGCA,fig_TGCC,fig_TGCG,fig_TGCT,fig_TGGA,fig_TGGC,fig_TGGG,fig_TGGT,fig_TGTA,fig_TGTC,fig_TGTG,fig_TGTT,fig_TTAA,fig_TTAC,fig_TTAG,fig_TTAT,fig_TTCA,fig_TTCC,fig_TTCG,fig_TTCT,fig_TTGA,fig_TTGC,fig_TTGG,fig_TTGT,fig_TTTA,fig_TTTC,fig_TTTG,fig_TTTT]
        #print total_tetras

        #print tetra_dict

            fig_list_var=','.join(map(str,fig_list))

        else:
            fig_list_var="masked positions over specified limit"

        print str(subseq_name) + "," + str(fig_list_var)     #print window_name and list of tetranucleotide frequencies for window, all comma-delimited





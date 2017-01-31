#!/usr/bin/perl -w
# parsing tab-delimited fastq files based on biases at "start"-1 and "end"+1 positions
# recommend the following command line:
# shuf $file | perl start_end_bias.pl - | head -n # < outfile
# or
# shuf $file | perl start_end_bias.pl - | head -n # | awk 'BEGIN {OFS="\t"}  {print $1,substr($2,2,30),$3,substr($4,2,30)}' > outfile


use strict;

# get filename and read count from command line argument

my $fqtd =$ARGV[0];

#my $readcount = $ARGV[1];


open (FILE, $fqtd) or die ("Could not open file.\n");




while (my $line = <FILE>) {
    chomp $line;                                    # chomp line
    my @line = split("\t",$line);                   # explode line into array
    #print @line,"\n";
    
    my $seq = $line[1];
    my @seqarray = split ("",$seq);
    #print @seqarray,"\n";
    
    my @begseq;
    my @endseq;
    
    
    
    # add misincorporations to 5' end of sequence

    
    if ($seqarray[0] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <4001) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
        
        }
    }
    else {
        push (@begseq,$seqarray[0]);
    }
    
    
    if ($seqarray[1] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <3001) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[1]);
    }


    if ($seqarray[2] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <2251) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[2]);
    }

    
    if ($seqarray[3] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <1689) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[3]);
    }

    
    if ($seqarray[4] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <1267) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[4]);
    }

    
    if ($seqarray[5] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <950) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[5]);
    }

    
    if ($seqarray[6] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <713) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[6]);
    }

    
    if ($seqarray[7] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <535) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[7]);
    }

    
    if ($seqarray[8] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <401) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[8]);
    }

    
    if ($seqarray[9] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <301) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[9]);
    }

    
    if ($seqarray[10] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <226) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[10]);
    }

    
    if ($seqarray[11] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <170) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[11]);
    }

   
    if ($seqarray[12] =~ m/C/) {
        my $rand = int(rand(10000))+1;
        if ($rand <128) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[12]);
    }

    
    
    
    
    
    # add misincorporations to 3' end of sequence
    
    
    if ($seqarray[-1] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <4001) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-1]);
    }
    
    
    if ($seqarray[-2] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <3001) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-2]);
    }
    
    
    if ($seqarray[-3] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <2251) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-3]);
    }
    
    
    if ($seqarray[-4] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <1689) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-4]);
    }
    
    
    if ($seqarray[-5] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <1267) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-5]);
    }
    
    
    if ($seqarray[-6] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <950) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-6]);
    }
    
    
    if ($seqarray[-7] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <713) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-7]);
    }
    
    
    if ($seqarray[-8] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <535) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-8]);
    }
    
    
    if ($seqarray[-9] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <401) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-9]);
    }
    
    
    if ($seqarray[-10] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <301) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-10]);
    }
    
    
    if ($seqarray[-11] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <226) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-11]);
    }
    
    
    if ($seqarray[-12] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <170) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-12]);
    }

    
    if ($seqarray[-13] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <128) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-13]);
    }

    
    splice (@seqarray,0,13,@begseq);
    splice (@seqarray,-13,13,@endseq);
    
    
    #print @begseq,"\n";
    #print @endseq,"\n";
    #print @seqarray,"\n";
    #print $line[1],"\n";
    
    my $newseq = join("",@seqarray);
    #print $newseq,"\n";
    splice (@line,1,1,$newseq);
    print $line[0],"\t",$line[1],"\t",$line[2],"\t",$line[3],"\n";
    
    
    
}






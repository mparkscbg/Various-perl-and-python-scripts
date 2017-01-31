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
        if ($rand <1201) {
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
        if ($rand <901) {
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
        if ($rand <676) {
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
        if ($rand <507) {
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
        if ($rand <381) {
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
        if ($rand <286) {
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
        if ($rand <215) {
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
        if ($rand <161) {
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
        if ($rand <121) {
            push (@begseq,"T");
            
        } else {
            push (@begseq,"C");
            
        }
    }
    else {
        push (@begseq,$seqarray[8]);
    }


    
    
    
    
    
    # add misincorporations to 3' end of sequence
    
    
    if ($seqarray[-1] =~ m/G/) {
        my $rand = int(rand(10000))+1;
        if ($rand <1201) {
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
        if ($rand <901) {
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
        if ($rand <676) {
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
        if ($rand <507) {
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
        if ($rand <381) {
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
        if ($rand <286) {
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
        if ($rand <215) {
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
        if ($rand <161) {
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
        if ($rand <121) {
            unshift (@endseq,"A");
            
        } else {
            unshift (@endseq,"G");
            
        }
    }
    else {
        unshift (@endseq,$seqarray[-9]);
    }
    

    
    
    
    
    splice (@seqarray,0,9,@begseq);
    splice (@seqarray,-9,9,@endseq);
    
    
    #print @begseq,"\n";
    #print @endseq,"\n";
    #print @seqarray,"\n";
    #print $line[1],"\n";
    
    my $newseq = join("",@seqarray);
    #print $newseq,"\n";
    splice (@line,1,1,$newseq);
    print $line[0],"\t",$line[1],"\t",$line[2],"\t",$line[3],"\n";
    
    
    
}






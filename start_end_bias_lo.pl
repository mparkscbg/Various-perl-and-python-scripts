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
    my $rand = int(rand(100))+1;                    # generate random number inclusive of range 1-100

        

    if ($line[1] =~ m/^[AG]/ && $line[1] =~ m/[CT]$/ && $rand>0) {
        print $line[0],"\t",$line[1],"\t",$line[2],"\t",$line[3],"\n";
    }
    elsif ($line[1] =~ m/^[AG]/ && $line[1] =~ m/[AG]$/ && $rand>25) {
        print $line[0],"\t",$line[1],"\t",$line[2],"\t",$line[3],"\n";
    }
    elsif ($line[1] =~ m/^[CT]/ && $line[1] =~ m/[CT]$/ && $rand>25) {
        print $line[0],"\t",$line[1],"\t",$line[2],"\t",$line[3],"\n";
    }
    elsif ($line[1] =~ m/^[CT]/ && $line[1] =~ m/[AG]$/ && $rand>50) {
        print $line[0],"\t",$line[1],"\t",$line[2],"\t",$line[3],"\n";
    }
}



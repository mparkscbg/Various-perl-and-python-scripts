#!/usr/bin/perl -w
# partitioning tab-delimited fastq files based into 54-74bp length partitions (average 64bp)
# recommend the following command line:
# perl read_partitioning.pl FILE < OUTFILE

#Author: Matthew Parks, 2014



use strict;

# get filename and read count from command line argument

my $fqtd =$ARGV[0];



open (FILE, $fqtd) or die ("Could not open file.\n");




while (my $line = <FILE>) {
    chomp $line;                                    # chomp line
    my @line = split("\t",$line);                   # explode line into array
    
    
    
    # determine which length pool read goes into and trim read to that length

    
        my $rand = int(rand(10000))+1;

    if ($rand <=1500) {
        print $line[0],"\n",substr($line[1],0,64),"\n",$line[2],"\n",substr($line[3],0,64),"\n";
    }
    
    elsif ($rand >= 1501 && $rand <= 2625){
        print $line[0],"\n",substr($line[1],0,65),"\n",$line[2],"\n",substr($line[3],0,65),"\n";
        
    }

    elsif ($rand >= 2626 && $rand <= 3750){
        print $line[0],"\n",substr($line[1],0,63),"\n",$line[2],"\n",substr($line[3],0,63),"\n";
        
    }

    elsif ($rand >= 3751 && $rand <= 4594){
        print $line[0],"\n",substr($line[1],0,66),"\n",$line[2],"\n",substr($line[3],0,66),"\n";
        
    }

    elsif ($rand >= 4595 && $rand <= 5438){
        print $line[0],"\n",substr($line[1],0,62),"\n",$line[2],"\n",substr($line[3],0,62),"\n";
        
    }

    elsif ($rand >= 5439 && $rand <= 6070){
        print $line[0],"\n",substr($line[1],0,67),"\n",$line[2],"\n",substr($line[3],0,67),"\n";
        
    }

    elsif ($rand >= 6071 && $rand <= 6703){
        print $line[0],"\n",substr($line[1],0,61),"\n",$line[2],"\n",substr($line[3],0,61),"\n";
        
    }

    elsif ($rand >= 6704 && $rand <= 7178){
        print $line[0],"\n",substr($line[1],0,68),"\n",$line[2],"\n",substr($line[3],0,68),"\n";
        
    }

    elsif ($rand >= 7179 && $rand <= 7652){
        print $line[0],"\n",substr($line[1],0,60),"\n",$line[2],"\n",substr($line[3],0,60),"\n";
        
    }

    elsif ($rand >= 7653 && $rand <= 8008){
        print $line[0],"\n",substr($line[1],0,69),"\n",$line[2],"\n",substr($line[3],0,69),"\n";
        
    }

    elsif ($rand >= 8009 && $rand <= 8364){
        print $line[0],"\n",substr($line[1],0,59),"\n",$line[2],"\n",substr($line[3],0,59),"\n";
        
    }

    elsif ($rand >= 8365 && $rand <= 8631){
        print $line[0],"\n",substr($line[1],0,70),"\n",$line[2],"\n",substr($line[3],0,70),"\n";
        
    }

    elsif ($rand >= 8632 && $rand <= 8898){
        print $line[0],"\n",substr($line[1],0,58),"\n",$line[2],"\n",substr($line[3],0,58),"\n";
        
    }

    elsif ($rand >= 8899 && $rand <= 9098){
        print $line[0],"\n",substr($line[1],0,71),"\n",$line[2],"\n",substr($line[3],0,71),"\n";
        
    }

    elsif ($rand >= 9099 && $rand <= 9299){
        print $line[0],"\n",substr($line[1],0,57),"\n",$line[2],"\n",substr($line[3],0,57),"\n";
        
    }

    elsif ($rand >= 9300 && $rand <= 9449){
        print $line[0],"\n",substr($line[1],0,72),"\n",$line[2],"\n",substr($line[3],0,72),"\n";
        
    }

    elsif ($rand >= 9450 && $rand <= 9599){
        print $line[0],"\n",substr($line[1],0,56),"\n",$line[2],"\n",substr($line[3],0,56),"\n";
        
    }

    elsif ($rand >= 9600 && $rand <= 9712){
        print $line[0],"\n",substr($line[1],0,73),"\n",$line[2],"\n",substr($line[3],0,73),"\n";
        
    }

    elsif ($rand >= 9713 && $rand <= 9824){
        print $line[0],"\n",substr($line[1],0,55),"\n",$line[2],"\n",substr($line[3],0,55),"\n";
        
    }

    elsif ($rand >= 9825 && $rand <= 9912){
        print $line[0],"\n",substr($line[1],0,74),"\n",$line[2],"\n",substr($line[3],0,74),"\n";
        
    }

    elsif ($rand >= 9913 && $rand <= 10000){
        print $line[0],"\n",substr($line[1],0,54),"\n",$line[2],"\n",substr($line[3],0,54),"\n";
        
    }

}




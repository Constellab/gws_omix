#!/usr/bin/perl

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

use strict;
use warnings;

# Reading genes <-> EC number(s) table file

my$EC_nb_file=$ARGV[0]; # First parameter 
my$blast_output_file=$ARGV[1];
my%hash_table_EC; # Hash table initialisation

open(ECFILE,$EC_nb_file);

while(<ECFILE>){
    chomp;
    if($_=~/^#/){
        next;
    }
    else{
        my@t=split/\t/;
        if(defined($t[7])){
            $hash_table_EC{$t[0]}=$t[7]; # Filling the hash table : key = gene id ; value = EC number(s)
        }
        else{
            $hash_table_EC{$t[0]}="NA";
        }
    }
}

close(ECFILE);

# Going through blast output results (read through standard input stream aka STDIN)

open(BLFILE,$blast_output_file);

while(<BLFILE>){
    chomp;
    if($_=~/^#/){
        next;
    }
    else{
        my@t=split/\t/;
        my@tGeneID=split("#", $t[1]);
        print $_,"\t",$hash_table_EC{$tGeneID[1]},"\n"; # Printing in the STDOUT blast file lines + corresponding EC Nb
    }
} 
 
close(BLFILE);


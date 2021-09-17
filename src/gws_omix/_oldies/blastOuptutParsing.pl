#!/usr/bin/perl

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

use strict;
use warnings;

my$identity=$ARGV[0];
my$blast_parsed=$ARGV[1];

open(BLFILE,$blast_parsed);
while(<BLFILE>){
    chomp;
    my@t=split/\t/;
    if($t[2] >= $identity){
        print $_,"\n";
    }
    else{
        next;
    }
}

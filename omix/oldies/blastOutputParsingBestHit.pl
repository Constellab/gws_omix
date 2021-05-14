#!/usr/bin/perl

# This software is the exclusive property of Gencovery SAS. 
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

use strict;
use warnings;

my%hash_table;

while(<STDIN>){
    chomp;
    my@t=split/\t/;
    print ++$hash_table{$t[0]},"\t",$_,"\n";
}

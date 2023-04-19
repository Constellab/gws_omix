#!/bin/bash

gmt_file=$1

# Extract gene universe from GMT file

echo "First two lines of the gmt file : " ;
head -2 $gmt_file ;
cat $gmt_file  | perl -ne 'chomp; @t=split/\t/; $cpt=0;  foreach(@t){ $cpt++; if($cpt>2){ $h{$_}{$t[0]}++;}} END{ foreach $k (sort keys %h){ print $k."\t"; foreach $l (sort keys %{$h{$k}}){ print ",".$l; } print "\n";  } } ' | sed 's/\t\,/\t/' > TopGO2_gene_universe.csv ;

echo "GMT convertion DONE" ;
echo  "First two lines of the gmt file : " ;
head -2 TopGO2_gene_universe.csv ;


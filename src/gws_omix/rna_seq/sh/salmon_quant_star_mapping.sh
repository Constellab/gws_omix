
gffread -w transcripto.tmp.fa -g $genomeFasta $annotation ;

cat transcripto.tmp.fa  | cut -d " " -f 1 > transcripto.tmp.2.fa ;

rm transcripto.tmp.fa ;

for i in $bamFolder
do
    salmon quant  -p $threads  -t transcripto.tmp.2.fa -l $experiment_opt -a $i -o $i.count.output
done

rm transcripto.tmp.2.fa

salmon quantmerge --quants *.count.output --column tpm -o salmon_quantmerge.tpm.txt ;
salmon quantmerge --quants *.count.output --column numreads -o salmon_quantmerge.raw_count.txt ;




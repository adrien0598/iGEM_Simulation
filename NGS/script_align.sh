#!/bin/bash
cd X204SC21093372-Z01-F001/raw_data
for i in A*; do
    name=$i.txt
    cd $i
    for j in *1.fq.gz; do
        doc1=$j
    done
    for j in *2.fq.gz; do
        doc2=$j
    done
    cd ../../..
    doc1=X204SC21093372-Z01-F001/raw_data/$i/$doc1
    doc2=X204SC21093372-Z01-F001/raw_data/$i/$doc2

    minimap2 -t 8 -a -x sr ref.fa $doc1 $doc2  | \
    samtools fixmate -O bam,level=0 -m - - | \
    samtools sort -l 1 -@8 -o sort.bam -T /tmp/example_prefix
    samtools markdup -O bam,level=1 sort.bam final.bam
    samtools mpileup -f ref.fa final.bam >myPileup.pileup
    java -jar VarScan.v2.3.9.jar pileup2snp myPileup.pileup --min-coverage 1 --min-reads2 1 --min-avg-qual 40 --min-var-freq 0.0000000001 --p-value 0.99 >$name
    cd X204SC21093372-Z01-F001/raw_data
done


#!/bin/bash

for i in $(cat /mnt/raid/projects/99Lives_analysis/scripts/chr.list); do

sleep 1s
(
INDIR=/mnt/raid/projects/99Lives_analysis/result/annotation/snp_impact
INFILE=$i.snp.impact.ens.tab
INVCFDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp_filter
INVCFFILE=$i.filter.vcf.gz
OUTDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp_impact
OUTFILE=$i.snp.impact.vcf
TMPDIR=/mnt/raid/projects/99Lives_analysis/TMP

cat $INDIR/$INFILE | cut -f2 | sed "s/:/\t/g" > $TMPDIR/$i.sites

vcftools --gzvcf $INVCFDIR/$INVCFFILE --positions $TMPDIR/$i.sites --stdout --recode --recode-INFO-all | bgzip -c > $OUTDIR/$OUTFILE.gz
)&

done

wait


for i in $(cat /mnt/raid/projects/99Lives_analysis/scripts/chr.list | grep -v UNMAPPED); do

sleep 1s
(
INDIR=/mnt/raid/projects/99Lives_analysis/result/annotation/indel_impact
INFILE=$i.indel.impact.ens.tab
INVCFDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel_filter
INVCFFILE=$i.filter.vcf.gz
OUTDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel_impact
OUTFILE=$i.indel.impact.vcf
TMPDIR=/mnt/raid/projects/99Lives_analysis/TMP

cat $INDIR/$INFILE | cut -f1 | sed "s/_/\t/g" | awk '{print $1"\t"$2-1}' | grep -v "#" > $TMPDIR/$i.sites

vcftools --gzvcf $INVCFDIR/$INVCFFILE --positions $TMPDIR/$i.sites --stdout --recode --recode-INFO-all | bgzip -c > $OUTDIR/$OUTFILE.gz
)&
done

i=UNMAPPED

INDIR=/mnt/raid/projects/99Lives_analysis/result/annotation/indel_impact
INFILE=$i.indel.impact.ens.tab
INVCFDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel_filter
INVCFFILE=$i.filter.vcf.gz
OUTDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel_impact
OUTFILE=$i.indel.impact.vcf
TMPDIR=/mnt/raid/projects/99Lives_analysis/TMP

cat $INDIR/UNMAPPED.indel.impact.ens.tab | cut -f1 | sed "s/_/\n/g" | grep -v [[:alpha:]] | grep -v "\." | awk '{print $1-1}' > $TMPDIR/$i.pos
cat $INDIR/UNMAPPED.indel.impact.ens.tab | cut -f2 | cut -f1 -d":" | grep -v "#" | grep -v "Location" > $TMPDIR/$i.chr

paste $TMPDIR/$i.chr $TMPDIR/$i.pos > $TMPDIR/$i.sites
vcftools --gzvcf $INVCFDIR/$INVCFFILE --positions $TMPDIR/$i.sites --stdout --recode --recode-INFO-all | bgzip -c > $OUTDIR/$OUTFILE.gz




#!/bin/bash


INDIR="/mnt/raid/projects/99Lives_analysis/vcf_200117/thin_snps"
OUTDIR="/mnt/raid/projects/99Lives_analysis/result/thin_snp_stats"



for i in $(cat /mnt/raid/projects/99Lives_analysis/99Lives_scripts/chr.list | grep -v UNMAPPED); do
sleep 1s
(  
vcftools --gzvcf $INDIR/combined.thin_10kb.vcf.gz --out $OUTDIR/$i.combined.10kb --LROH --chr $i
)&
done


sleep 1s
vcftools --gzvcf $INDIR/combined.thin_10kb.vcf.gz --out $OUTDIR/combined.10kb --geno-r2 &
sleep 1s
vcftools --gzvcf $INDIR/combined.thin_1kb.vcf.gz --out $OUTDIR/combined.1kb --relatedness2 &
sleep 1s
vcftools --gzvcf $INDIR/combined.thin_1kb.vcf.gz --out $OUTDIR/combined.1kb --geno-r2 &

wait


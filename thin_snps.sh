#!/bin/bash

INDIR="/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp_filter"
OUTDIR="/mnt/raid/projects/99Lives_analysis/vcf_200117/thin_snps"


for i in $(cat /mnt/raid/projects/99Lives_analysis/99Lives_scripts/chr.list); do
sleep 1s
(  
vcftools --gzvcf $INDIR/$i.filter.vcf.gz --stdout --non-ref-af 0.01 --max-non-ref-af 0.99 --recode --recode-INFO-all |
vcftools --vcf - --thin 1000 --recode --recode-INFO-all --stdout | bgzip -c > $OUTDIR/$i.thin.1maf.vcf.gz
tabix $OUTDIR/$i.thin.1maf.vcf.gz
)&
done

wait

vcfcombine $(echo $(ls $OUTDIR/*thin.1maf.vcf.gz)) | bgzip -c > $OUTDIR/combined.thin_1kb.1maf.vcf.gz
tabix $OUTDIR/combined.thin_1kb.1maf.vcf.gz

vcftools --gzvcf $OUTDIR/combined.thin_1kb.1maf.vcf.gz --thin 10000 --recode --recode-INFO-all --stdout | bgzip -c > $OUTDIR/combined.thin_10kb.1maf.vcf.gz
tabix $OUTDIR/combined.thin_10kb.1maf.vcf.gz


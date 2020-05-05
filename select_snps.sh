#/bin/bash

INDIR="/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp_filter"
OUTDIR="/mnt/raid/projects/99Lives_analysis/vcf_200117/thin_snps"

SELECT="/mnt/raid/projects/99Lives_analysis/accessory_data/DNA_chip_sites.tsv"


for i in $(cat /mnt/raid/projects/99Lives_analysis/99Lives_scripts/chr.list); do
sleep 1s
(  
vcftools --gzvcf $INDIR/$i.filter.vcf.gz --stdout --positions $SELECT --recode --recode-INFO-all | bgzip -c > $OUTDIR/$i.select.chip.vcf.gz
tabix $OUTDIR/$i.select.chip.vcf.gz
)&
done

wait

vcfcombine $(echo $(ls $OUTDIR/*select.chip.vcf.gz | grep -v combined)) | bgzip -c > $OUTDIR/combined.select.chip.vcf.gz
tabix $OUTDIR/combined.select.chip.vcf.gz



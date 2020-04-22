#!/bin/bash


INFILT=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp_filter
INFULL=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp/200117_domestics_WGS
TMP=/mnt/raid/projects/99Lives_analysis/TMP
OUTDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp_rm


for i in $(cat /mnt/raid/projects/99Lives_analysis/99Lives_scripts/chr.list); do
sleep 1s
(

vcftools --gzvcf $INFULL.$i.snps.recalibrated.filtered.vcf.gz --gzdiff $INFILT/$i.filter.vcf.gz --diff-site --stdout | grep "\.$" | cut -f1,2 > $TMP/$i.snp.diff-site

vcftools --gzvcf $INFULL.$i.snps.recalibrated.filtered.vcf.gz --positions $TMP/$i.snp.diff-site --recode --recode-INFO-all --stdout |
vcfremovesamples - felCat.CVB13769.20679 felCat.Fcat20796.BIRA003CT189 felCat.SAMEA104694019.K510 felCat.SAMEA1061687.K499 felCat.SAMEA1061689.K515 felCat.SAMEA1061799.K511 felCat.SAMEA2785958.K450 felCat.Zivziv.Zivziv felCat.cat_07.Tesla felCat.cat_09.Moses felCat.cat_10.Louie felCat.cat_17.Sundance felCat.cat_24.Gumbo felCat.cat_33.D05-0399 felCat.Fcat12682.Cinnamon |
vcffixup - |
vcffilter -f "AC > 0" | bgzip -c > $OUTDIR/$i.snp.rm.vcf.gz

)&
done

wait

INFILT=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel_filter
INFULL=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel/200117_domestics_WGS
TMP=/mnt/raid/projects/99Lives_analysis/TMP
OUTDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel_rm


for i in $(cat /mnt/raid/projects/99Lives_analysis/99Lives_scripts/chr.list); do
sleep 1s
(

vcftools --gzvcf $INFULL.$i.indels.recalibrated.filtered.vcf.gz --gzdiff $INFILT/$i.filter.vcf.gz --diff-site --stdout | grep "\.$" | cut -f1,2 > $TMP/$i.indel.diff-site

vcftools --gzvcf $INFULL.$i.indels.recalibrated.filtered.vcf.gz --positions $TMP/$i.indel.diff-site --recode --recode-INFO-all --stdout |
vcfremovesamples - felCat.CVB13769.20679 felCat.Fcat20796.BIRA003CT189 felCat.SAMEA104694019.K510 felCat.SAMEA1061687.K499 felCat.SAMEA1061689.K515 felCat.SAMEA1061799.K511 felCat.SAMEA2785958.K450 felCat.Zivziv.Zivziv felCat.cat_07.Tesla felCat.cat_09.Moses felCat.cat_10.Louie felCat.cat_17.Sundance felCat.cat_24.Gumbo felCat.cat_33.D05-0399 felCat.Fcat12682.Cinnamon |
vcffixup - |
vcffilter -f "AC > 0" | bgzip -c > $OUTDIR/$i.indel.rm.vcf.gz

)&
done

wait
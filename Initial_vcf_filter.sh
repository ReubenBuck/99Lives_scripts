#!/bin/bash


INDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp
INFILE=200117_domestics_WGS.chrE3.snps.recalibrated.filtered.vcf.gz
OUTDIR=/mnt/raid/projects/99Lives_analysis/result/first_pass
OUTFILE=chrE3_vcf_200117


#genotypeSummary --type PL --target $(seq -s , 0 283) --file $INDIR/$INFILE > $OUTDIR/$OUTFILE.gt_sum.tsv

vcftools --gzvcf $INDIR/$INFILE --depth --out $OUTDIR/$OUTFILE




for i in $(cat /mnt/raid/projects/99Lives_analysis/scripts/chr.list); do
sleep 1s
(
INDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp
INFILE=200117_domestics_WGS.$i.snps.recalibrated.filtered.vcf.gz
OUTDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp_filter
OUTFILE=$i.filter.vcf

vcfremovesamples $INDIR/$INFILE felCat.CVB13769.20679 felCat.Fcat20796.BIRA003CT189 felCat.SAMEA104694019.K510 felCat.SAMEA1061687.K499 felCat.SAMEA1061689.K515 felCat.SAMEA1061799.K511 felCat.SAMEA2785958.K450 felCat.Zivziv.Zivziv felCat.cat_07.Tesla felCat.cat_09.Moses felCat.cat_10.Louie felCat.cat_17.Sundance felCat.cat_24.Gumbo felCat.cat_33.D05-0399 felCat.Fcat12682.Cinnamon |
/home/buckley/software/vcflib/scripts/vcfbiallelic | vcffilter -f "FILTER = PASS" | vcffixup - |
vcffilter -f "AC > 0" | bgzip -c > $OUTDIR/$OUTFILE.gz
)&
done


for i in $(cat /mnt/raid/projects/99Lives_analysis/scripts/chr.list); do
sleep 1s
(
INDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel
INFILE=200117_domestics_WGS.$i.indels.recalibrated.filtered.vcf.gz
OUTDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel_filter
OUTFILE=$i.filter.vcf

vcfremovesamples $INDIR/$INFILE felCat.CVB13769.20679 felCat.Fcat20796.BIRA003CT189 felCat.SAMEA104694019.K510 felCat.SAMEA1061687.K499 felCat.SAMEA1061689.K515 felCat.SAMEA1061799.K511 felCat.SAMEA2785958.K450 felCat.Zivziv.Zivziv felCat.cat_07.Tesla felCat.cat_09.Moses felCat.cat_10.Louie felCat.cat_17.Sundance felCat.cat_24.Gumbo felCat.cat_33.D05-0399 felCat.Fcat12682.Cinnamon |
/home/buckley/software/vcflib/scripts/vcfbiallelic | vcffilter -f "FILTER = PASS" | vcffixup - |
vcffilter -f "AC > 0" | bgzip -c > $OUTDIR/$OUTFILE.gz
)&
done




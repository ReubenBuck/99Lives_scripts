#!/bin/bash

declare -a DIRarr=(
	"indel_rm"
	"indel_filter"
	"snp_rm"
	"snp_filter"
	)

declare -a INpre=(
	""
	""
	""
	""
	)

declare -a INsuf=(
	"indel.rm.vcf.gz"
	"filter.vcf.gz"
	"snp.rm.vcf.gz"
	"filter.vcf.gz"
	)

declare -a OUTsuf=(
	"indel.rm.stats"
	"indel.filter.stats"
	"snp.rm.stats"
	"snp.filter.stats"
	)

declare -a SAMPnum=(
	268
	268
	268
	268
	)


INDIR="/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file"

OUTDIR=/mnt/raid/projects/99Lives_analysis/result/vcfStats_advanced


for i in $(seq 0 3); do
for j in $(cat /mnt/raid/projects/99Lives_analysis/99Lives_scripts/chr.list); do
sleep 1s
(
INFILE=$INDIR/${DIRarr[$i]}/${INpre[$i]}$j.${INsuf[$i]}
OUTFILE=$OUTDIR/${DIRarr[$i]}/$j.${OUTsuf[$i]}

genotypeSummary --type PL --target $(echo $(seq -s"," 0 ${SAMPnum[$i]})) --file $INFILE > $OUTFILE.gtSum

vcftools --gzvcf $INFILE --out $OUTFILE --depth

vcftools --gzvcf $INFILE --out $OUTFILE --TajimaD 10000

vcftools --gzvcf $INFILE --out $OUTFILE --SNPdensity 10000

vcftools --gzvcf $INFILE --out $OUTFILE --singletons

vcftools --gzvcf $INFILE --out $OUTFILE --mendel /mnt/raid/projects/99Lives_analysis/accessory_data/99Lives.ped.csv

vcftools --gzvcf $INFILE --out $OUTFILE --indv-freq-burden

vcftools --gzvcf $INFILE --out $OUTFILE --freq

vcftools --gzvcf $INFILE --out $OUTFILE --counts

vcftools --gzvcf $INFILE --out $OUTFILE --hardy

vcftools --gzvcf $INFILE --out $OUTFILE --missing-site

vcftools --gzvcf $INFILE --out $OUTFILE --missing-indv
)&
done
wait
done
wait




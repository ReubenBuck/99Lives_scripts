#!/bin/bash

declare -a DIRarr=(
	"indel"
	"indel_filter"
	"indel_impact"
	"snp"
	"snp_filter"
	"snp_impact"
	)

declare -a INpre=(
	"200117_domestics_WGS."
	""
	""
	"200117_domestics_WGS."
	""
	""
	)

declare -a INsuf=(
	"indels.recalibrated.filtered.vcf.gz"
	"filter.vcf.gz"
	"indel.impact.vcf.gz"
	"snps.recalibrated.filtered.vcf.gz"
	"filter.vcf.gz"
	"snp.impact.vcf.gz"
	)

declare -a OUTsuf=(
	"indels.all.stats"
	"indel.filter.stats"
	"indel.impact.stats"
	"snps.all.stats"
	"snp.filter.stats"
	"snp.impact.stats"
	)


INDIR="/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file"

OUTDIR=/mnt/raid/projects/99Lives_analysis/result/vcfStats



echo ${INpre[$i]}$j.${INsuf[$i]}

for i in $(seq 0 5); do
for j in $(cat /mnt/raid/projects/99Lives_analysis/99Lives_scripts/chr.list); do
sleep 1s
(
vcfstats $INDIR/${DIRarr[$i]}/${INpre[$i]}$j.${INsuf[$i]} > $OUTDIR/${DIRarr[$i]}/$j.${OUTsuf[$i]}
)&
done
wait
done
wait
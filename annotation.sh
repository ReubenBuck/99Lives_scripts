#!/bin/bash

# lets do the annotation
#


for i in $(cat /mnt/raid/projects/99Lives_analysis/scripts/chr.list); do
sleep 1s
(
INDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp_filter
INFILE=$i.filter.vcf.gz
OUTDIR=/mnt/raid/projects/99Lives_analysis/result/annotation/snp
OUTFILE=$i.snp.ens.tab

~/src/ensembl-vep/vep \
-i $INDIR/$INFILE \
-o $OUTDIR/$OUTFILE.gz \
--species felis_catus \
--cache \
--offline \
--compress_output gzip \
--symbol \
--canonical \
--tab
)&
done

wait

for i in $(cat /mnt/raid/projects/99Lives_analysis/scripts/chr.list); do
sleep 1s
(
INDIR=/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel_filter
INFILE=$i.filter.vcf.gz
OUTDIR=/mnt/raid/projects/99Lives_analysis/result/annotation/indel
OUTFILE=$i.indel.ens.tab

~/src/ensembl-vep/vep \
-i $INDIR/$INFILE \
-o $OUTDIR/$OUTFILE.gz \
--species felis_catus \
--cache \
--offline \
--compress_output gzip \
--symbol \
--canonical \
--tab
)&
done

wait



for i in $(cat /mnt/raid/projects/99Lives_analysis/scripts/chr.list); do
sleep 1s
(
INDIR=/mnt/raid/projects/99Lives_analysis/result/annotation/snp
INFILE=$i.snp.ens.tab.gz
OUTDIR=/mnt/raid/projects/99Lives_analysis/result/annotation/snp_impact
OUTFILE=$i.snp.impact.ens.tab


gzip -cd $INDIR/$INFILE | grep -v "MODIFIER" > $OUTDIR/$OUTFILE

)&
done

wait





for i in $(cat /mnt/raid/projects/99Lives_analysis/scripts/chr.list); do
sleep 1s
(
INDIR=/mnt/raid/projects/99Lives_analysis/result/annotation/indel
INFILE=$i.indel.ens.tab.gz
OUTDIR=/mnt/raid/projects/99Lives_analysis/result/annotation/indel_impact
OUTFILE=$i.indel.impact.ens.tab


gzip -cd $INDIR/$INFILE | grep -v "MODIFIER" > $OUTDIR/$OUTFILE

)&
done

wait



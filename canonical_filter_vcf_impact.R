rm(list = ls())
options(stringsAsFactors = FALSE)
library(vcfR)

files.vcf <- list.files("/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp_impact/")
files.vcf <- files.vcf[grep("gz$", files.vcf)]
files.impact <- list.files("/mnt/raid/projects/99Lives_analysis/result/annotation/snp_impact/")


for(i in 1:length(files.vcf)){
  impact0 <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/annotation/snp_impact/", files.impact[i], sep = ""),
                        sep = "\t", header = FALSE,
                        col.names = c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
                                      "Consequence", "cDNA_position", "CDS_position", "Protein_position", 
                                      "Amino_acids", "Codons", "Existing_variation", "IMPACT", "DISTANCE",
                                      "STRAND", "FLAGS", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL"))
  impact0 <- impact0[impact0$CANONICAL == "YES",]
  
  vcf0 <- read.vcfR(paste("/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp_impact/", files.vcf[i], sep = ""))
  var.name <- paste(getCHROM(vcf0), "_",getPOS(vcf0), "_",getREF(vcf0), "/", getALT(vcf0), sep = "")
  vcf0 <- vcf0[var.name %in% impact0$Uploaded_variation]
  if(i == 1){
    vcf <- vcf0
    impact <- impact0
  }else{
    vcf <- rbind2(vcf, vcf0)
    impact <- rbind(impact, impact0)
  }
}

write.vcf(vcf, "/mnt/raid/projects/99Lives_analysis/vcf_200117/combined_file/canonical_gene_impact_allchr.snp.vcf.gz")
write.table(impact, "/mnt/raid/projects/99Lives_analysis/result/annotation/snp_impact/canonical_gene_impact_allchr.snp.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# now to stitch together indel files
rm(list = ls())


files.vcf <- list.files("/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel_impact/")
files.vcf <- files.vcf[grep("gz$", files.vcf)]
files.impact <- list.files("/mnt/raid/projects/99Lives_analysis/result/annotation/indel_impact/")
files.impact <- files.impact[grep("indel.impact.ens.tab", files.impact)]


for(i in 1:length(files.vcf)){
  impact0 <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/annotation/indel_impact/", files.impact[i], sep = ""),
                        sep = "\t", header = FALSE,
                        col.names = c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
                                      "Consequence", "cDNA_position", "CDS_position", "Protein_position", 
                                      "Amino_acids", "Codons", "Existing_variation", "IMPACT", "DISTANCE",
                                      "STRAND", "FLAGS", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL"))
  if(nrow(impact0) == 0){
    next()
  }
  impact0 <- impact0[impact0$CANONICAL == "YES",]
  
  vcf0 <- read.vcfR(paste("/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/indel_impact/", files.vcf[i], sep = ""))
  df.ref.alt <- data.frame(ref = getREF(vcf0), alt = getALT(vcf0))
  df.indel <- data.frame(ins = nchar(df.ref.alt$ref) == 1, del = nchar(df.ref.alt$alt) == 1)
  df.new.ref.alt <- df.ref.alt
  df.new.ref.alt$ref[df.indel$ins] <- "-"
  df.new.ref.alt$alt[df.indel$ins] <- substring(text = df.new.ref.alt$alt[df.indel$ins], first = 2)
  df.new.ref.alt$alt[df.indel$del] <- "-"
  df.new.ref.alt$ref[df.indel$del] <- substring(text = df.new.ref.alt$ref[df.indel$del], first = 2)
  
  var.name <- paste(getCHROM(vcf0), "_",getPOS(vcf0) + 1, "_",df.new.ref.alt$ref, "/", df.new.ref.alt$alt, sep = "")
  var.name.orig <- paste(getCHROM(vcf0), "_",getPOS(vcf0), "_",getREF(vcf0), "/", getALT(vcf0), sep = "")
  names(var.name.orig) <- var.name
  
  impact0$Uploaded_variation <- var.name.orig[impact0$Uploaded_variation]
  
  vcf0 <- vcf0[var.name.orig %in% impact0$Uploaded_variation]
  
  
  
  if(i == 1){
    vcf <- vcf0
    impact <- impact0
  }else{
    vcf <- rbind2(vcf, vcf0)
    impact <- rbind(impact, impact0)
  }
}

write.vcf(vcf, "/mnt/raid/projects/99Lives_analysis/vcf_200117/combined_file/canonical_gene_impact_allchr.indel.vcf.gz")
write.table(impact, "/mnt/raid/projects/99Lives_analysis/result/annotation/indel_impact/canonical_gene_impact_allchr.indel.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


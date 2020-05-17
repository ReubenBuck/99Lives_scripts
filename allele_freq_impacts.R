rm(list = ls())
options(stringsAsFactors = FALSE)
library(GenomicRanges)
library(vcfR)

## allele frequency time

# tables

# counts of snps! and ts/tv ratio


vcf.snp <- read.vcfR("/mnt/raid/projects/99Lives_analysis/vcf_200117/combined_file/filtered_allchr.snp.vcf.gz")
gt.snp <- extract.gt(vcf.snp)
rownames(gt.snp) <- paste(getCHROM(vcf.snp), "_",getPOS(vcf.snp), "_",getREF(vcf.snp), "/", getALT(vcf.snp), sep = "")



alleles.snp <- data.frame(ref = rowSums(gt.snp == "0/0", na.rm = TRUE),
                          het = rowSums(gt.snp == "0/1", na.rm = TRUE),
                          hom = rowSums(gt.snp == "1/1", na.rm = TRUE),
                          row.names = rownames(gt.snp))
alleles.snp$ac <- alleles.snp$het + (alleles.snp$hom*2)
alleles.snp$an <- (alleles.snp$het*2) + (alleles.snp$hom*2) + (alleles.snp$ref*2)
alleles.snp$af <- alleles.snp$ac/alleles.snp$an
alleles.snp$maf <- alleles.snp$af
alleles.snp$maf[alleles.snp$maf > 0.5] <- 1 - alleles.snp$maf[alleles.snp$maf > 0.5]
alleles.snp$singleton <- NA
alleles.snp$singleton[alleles.snp$ac == 1] <- colnames(gt.snp)[t(col(gt.snp[alleles.snp$ac == 1,]))[t(gt.snp[alleles.snp$ac == 1,] == "0/1" & !is.na(gt.snp[alleles.snp$ac == 1,]))]]
alleles.snp$doubleton <- NA
alleles.snp$doubleton[alleles.snp$ac == 2 & alleles.snp$hom == 1] <- colnames(gt.snp)[t(col(gt.snp[alleles.snp$ac == 2 & alleles.snp$hom == 1,]))[t(gt.snp[alleles.snp$ac == 2 & alleles.snp$hom == 1,] == "1/1" & !is.na(gt.snp[alleles.snp$ac == 2 & alleles.snp$hom == 1,]))]]

impact.snp <- read.table("/mnt/raid/projects/99Lives_analysis/result/annotation/snp_impact/filtered_allchr.snp.tsv",
                         sep = "\t", header = TRUE)
impact.snp.gr <- GRanges(seqnames = Rle(sapply(strsplit(impact.snp$Location, ":"), "[[", 1)),
                         ranges = IRanges(start = as.integer(sapply(strsplit(impact.snp$Location, ":"), "[[", 2)), 
                                          width = 1),
                         impact = impact.snp$IMPACT, symbol = impact.snp$SYMBOL,
                         gene = impact.snp$Gene, var.name = impact.snp$Uploaded_variation)
elementMetadata(impact.snp.gr) <- cbind(elementMetadata(impact.snp.gr), alleles.snp[impact.snp$Uploaded_variation,])



vcf.indel <- read.vcfR("/mnt/raid/projects/99Lives_analysis/vcf_200117/combined_file/canonical_gene_impact_allchr.indel.vcf.gz")
gt.indel <- extract.gt(vcf.indel)
rownames(gt.indel) <- paste(getCHROM(vcf.indel), "_",getPOS(vcf.indel), "_",getREF(vcf.indel), "/", getALT(vcf.indel), sep = "")

alleles.indel <- data.frame(ref = rowSums(gt.indel == "0/0", na.rm = TRUE),
                            het = rowSums(gt.indel == "0/1", na.rm = TRUE),
                            hom = rowSums(gt.indel == "1/1", na.rm = TRUE),
                            row.names = rownames(gt.indel))
alleles.indel$ac <- alleles.indel$het + (alleles.indel$hom*2)
alleles.indel$an <- (alleles.indel$het*2) + (alleles.indel$hom*2) + (alleles.indel$ref*2)
alleles.indel$af <- alleles.indel$ac/alleles.indel$an
alleles.indel$maf <- alleles.indel$af
alleles.indel$maf[alleles.indel$maf > 0.5] <- 1 - alleles.indel$maf[alleles.indel$maf > 0.5]
alleles.indel$singleton <- NA
alleles.indel$singleton[alleles.indel$ac == 1] <- colnames(gt.indel)[t(col(gt.indel[alleles.indel$ac == 1,]))[t(gt.indel[alleles.indel$ac == 1,] == "0/1" & !is.na(gt.indel[alleles.indel$ac == 1,]))]]
alleles.indel$doubleton <- NA
alleles.indel$doubleton[alleles.indel$ac == 2 & alleles.indel$hom == 1] <- colnames(gt.indel)[t(col(gt.indel[alleles.indel$ac == 2 & alleles.indel$hom == 1,]))[t(gt.indel[alleles.indel$ac == 2 & alleles.indel$hom == 1,] == "1/1" & !is.na(gt.indel[alleles.indel$ac == 2 & alleles.indel$hom == 1,]))]]


impact.indel <- read.table("/mnt/raid/projects/99Lives_analysis/result/annotation/indel_impact/canonical_gene_impact_allchr.indel.tsv",
                           sep = "\t", header = TRUE)
indel.1bp <- !grepl("-",impact.indel$Location)
impact.indel[indel.1bp,"Location"] <- paste(impact.indel[indel.1bp,"Location"], "-",
                                            sapply(strsplit(impact.indel[indel.1bp,"Location"], ":"), "[[", 2),
                                            sep = "")
impact.indel.gr <- GRanges(seqnames = Rle(sapply(strsplit(impact.indel$Location, ":|-"), "[[", 1)),
                           ranges = IRanges(start = as.integer(sapply(strsplit(impact.indel$Location, ":|-"), "[[", 2)), 
                                            end = as.integer(sapply(strsplit(impact.indel$Location, ":|-"), "[[", 3))))


var.type <- c("snp", "indel")
df.all <- NULL
for(t in var.type){
  
  impact <- get(paste("impact.",t,sep = ""))
  alleles <- get(paste("alleles.",t,sep = ""))
  
  
  s.spl <- strsplit(impact$Consequence, split = ",")
  s.len <- sapply(s.spl,length)
  s.impact <- rep(impact$IMPACT, s.len)
  s.var <- rep(impact$Uploaded_variation, s.len)
  
  df <- data.frame(var = s.var,
                   impact = s.impact,
                   consequence = unlist(s.spl))
  df <- cbind(df,alleles[df$var,])
  df <- unique(df)
  df$type <- t
  
  df.all <- rbind(df.all, df)
  
}


# now we can do the summerising table
# just need to circle through the maf ranges 

maf_low <- c(0,0.01,0.10,0)
maf_high <- c(0.01,0.1,0.5,0.5)
imp <- c("LOW", "MODERATE", "HIGH")
cons <- unique(df.all$consequence)
df.all$consequence <-factor(df.all$consequence, levels = cons)
df.all$impact <- factor(df.all$impact, levels = imp)
df.table <- NULL
for(maf in 1:length(maf_low)){
  df0 <- df.all[df.all$maf > maf_low[maf] & df.all$maf <= maf_high[maf],]
  df1 <- unique(df0[,c("var","type","impact","consequence")])
  df1 <- df1[,c("type","impact","consequence")]
  df2 <- data.frame(table(df1))
  df1 <- unique(df0[,c("var","type","impact")])
  df1 <- df1[,c("impact","type")]
  df1$consequence <- "All"
  df2 <- rbind(df2,data.frame(table(df1)))
  df1 <- unique(df0[,c("var","type")])
  df1 <- df1[,c("type")]
  df1 <- data.frame(type = df1, impact = "All", consequence = "All")
  df2 <- rbind(df2,data.frame(table(df1)))
  df2$low_maf <- maf_low[maf]
  df2$high_maf <- maf_high[maf]
  df.table <- rbind(df.table, df2)
}

df.table$consequence <-factor(df.table$consequence, 
                              levels = c(as.character(cons), "All"))
df.table$impact <- factor(df.table$impact, 
                          levels = c(as.character(imp), "All"))
df.table$type <- factor(df.table$type, level = var.type)

df.table <- df.table[order(df.table$type, df.table$impact, df.table$high_maf, df.table$low_maf, df.table$consequence),]

df.filt <- df.table[,c("type", "impact", "consequence", "low_maf", "high_maf", "Freq")]
df.filt <- df.filt[df.filt$Freq > 0,]

write.table(df.filt,
            "/mnt/raid/projects/99Lives_analysis/result/allele_frequencies/impact_af.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


af <- read.table("/mnt/raid/projects/99Lives_analysis/result/vcfStats_advanced/snp_filter/all.frq",
                 sep = "\t", comment.char = "F")
af <- af$V1
maf <- af
maf[maf > 0.5] <- 1-maf[maf > 0.5]

layout(1)
pdf(file = "/mnt/raid/projects/99Lives_analysis/result/allele_frequencies/MAF_dist.pdf")
h <- hist(maf,breaks = seq(0,.5,0.01), xlab = "MAF", main = "")
dev.off()

write.table(data.frame(mids = h$mids, counts = h$counts, pct = h$counts/sum(h$counts) * 100),
            "/mnt/raid/projects/99Lives_analysis/result/allele_frequencies/MAF_dist.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)









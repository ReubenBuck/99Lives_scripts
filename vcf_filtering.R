rm(list = ls())
options(stringsAsFactors = FALSE)


setwd("/mnt/raid/projects/99Lives_analysis/result/vcfStats")

snp_all_files <- paste("snp/",list.files("snp/"), sep = "")
snp_filt_files <- paste("snp_filter/",list.files("snp_filter/"), sep = "")
snp_impact_files <- paste("snp_impact/",list.files("snp_impact/"), sep = "")

snp.tstv <- NULL
snp.sites <- NULL
snp_files <- c("snp_all_files", "snp_filt_files", "snp_impact_files")
df.snp.tstv <- NULL
df.snp.site <- NULL
for(i in snp_files){
  snp_file <- get(i)
  
  snp.tstv <- NULL
  snp.sites <- NULL
  
  for(j in snp_file){
    snp <- read.table(j, sep = "\t", fill = TRUE)
    snp.sites <- c(snp.sites, snp[1,2])
    snp.tstv <- c(snp.tstv, snp[10,2])
  }
  df.snp.site <- cbind(df.snp.site, as.numeric(snp.sites))
  df.snp.tstv <- cbind(df.snp.tstv, as.numeric(snp.tstv))
}

colnames(df.snp.site) <- snp_files
colnames(df.snp.tstv) <- snp_files
rownames(df.snp.tstv) <- rownames(df.snp.site) <- gsub(".snps.all.stats", "", gsub("^.+/","",snp_all_files))

wt.mean <- NULL
for(i in 1:3){
wt.mean <- c(wt.mean, weighted.mean(df.snp.tstv[,i], 
                                    w = df.snp.site[,i]/sum(df.snp.site[,i])))
}
df.snp.tstv <- rbind(df.snp.tstv, wt.mean)
df.snp.site <- rbind(df.snp.site, sum = colSums(df.snp.site))


write.csv(df.snp.site, "/mnt/raid/projects/99Lives_analysis/result/vcfStats/collected_snp_result/snp_site_filtering_stats.csv")
write.csv(df.snp.tstv, "/mnt/raid/projects/99Lives_analysis/result/vcfStats/collected_snp_result/snp_tstv_filtering_stats.csv")








indel_all_files <- paste("indel/",list.files("indel/"), sep = "")
indel_filt_files <- paste("indel_filter/",list.files("indel_filter/"), sep = "")
indel_impact_files <- paste("indel_impact/",list.files("indel_impact/"), sep = "")

indel.tstv <- NULL
indel.sites <- NULL
indel_files <- c("indel_all_files", "indel_filt_files", "indel_impact_files")
df.indel.tstv <- NULL
df.indel.site <- NULL
for(i in indel_files){
  indel_file <- get(i)
  
  indel.tstv <- NULL
  indel.sites <- NULL
  
  for(j in indel_file){
    indel <- read.table(j, sep = "\t", fill = TRUE)
    indel.sites <- c(indel.sites, indel[1,2])
    indel.tstv <- c(indel.tstv, indel[10,2])
  }
  df.indel.site <- cbind(df.indel.site, as.numeric(indel.sites))
  df.indel.tstv <- cbind(df.indel.tstv, as.numeric(indel.tstv))
}

colnames(df.indel.site) <- indel_files
colnames(df.indel.tstv) <- indel_files
rownames(df.indel.tstv) <- rownames(df.indel.site) <- gsub(".indels.all.stats", "", gsub("^.+/","",indel_all_files))

wt.mean <- NULL
for(i in 1:3){
  wt.mean <- c(wt.mean, weighted.mean(df.indel.tstv[,i], 
                                      w = df.indel.site[,i]/sum(df.indel.site[,i])))
}
df.indel.tstv <- rbind(df.indel.tstv, wt.mean)
df.indel.site <- rbind(df.indel.site, sum = colSums(df.indel.site))


write.csv(df.indel.site, "/mnt/raid/projects/99Lives_analysis/result/vcfStats/collected_snp_result/indel_site_filtering_stats.csv")
write.csv(df.indel.tstv, "/mnt/raid/projects/99Lives_analysis/result/vcfStats/collected_snp_result/indel_tstv_filtering_stats.csv")





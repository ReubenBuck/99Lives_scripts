rm(list = ls())
options(stringsAsFactors = FALSE)
library(vcfR)

args = commandArgs(trailingOnly=TRUE)


vcf_file = paste("/mnt/raid/projects/99Lives_analysis/vcf_200117/split_file/snp_filter/",args,".filter.vcf.gz", sep = "")

samp.vcf <- read.vcfR(vcf_file)

gt <- extract.gt(samp.vcf)
gt.var <- gt != "0/0" & !is.na(gt)

new_var_count <- NULL
total_var_count <- NULL
indiv_var_count <- NULL
indiv_var_med_count <- NULL

set.seed(100)

obj.list <- NULL
for( i in seq(10,ncol(gt.var),by = 10)){
  obj <- NULL
  for(j in 1:10){
    obj <- rbind(obj,
                 sample(x=1:ncol(gt.var),size = i, replace = FALSE))
  }
  obj.list <- c(obj.list, list(obj))
}

for(i in 1:length(obj.list)){
  new_vars1 <- NULL
  total_vars1 <- NULL
  indiv_vars1 <- NULL
  indiv_vars_med1 <- NULL
  for(j in 1:10){
    samp.gt.var <- gt.var[,obj.list[[i]][j,]]
    
    new_vars0 <- sum(rowSums(samp.gt.var, na.rm = TRUE) == 1, na.rm = TRUE)
    total_vars0 <- sum(rowSums(samp.gt.var, na.rm = TRUE) > 0, na.rm = TRUE)
    indiv_vars0 <- mean(colSums(samp.gt.var,na.rm = TRUE), na.rm = TRUE)
    indiv_vars_med0 <- mean(colSums(samp.gt.var,na.rm = TRUE), na.rm = TRUE)
    
    new_vars1 <- c(new_vars1, new_vars0)
    total_vars1 <- c(total_vars1, total_vars0)
    indiv_vars1 <- c(indiv_vars1, indiv_vars0)
    indiv_vars_med1 <- c(indiv_vars_med1, indiv_vars_med0)
  }
  new_var_count <- rbind(new_var_count, new_vars1)
  total_var_count <- rbind(total_var_count, total_vars1)
  indiv_var_count <- rbind(indiv_var_count, indiv_vars1)
  indiv_var_med_count <- rbind(indiv_var_med_count, indiv_vars_med1)
}


outfile_new <- paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/",args,".new_var.tsv", sep = "")

write.table(x = new_var_count,
            file = outfile_new, 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE)


outfile_total <- paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/",args,".total_var.tsv", sep = "")

write.table(x = total_var_count,
            file = outfile_total, 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE)



outfile_indiv <- paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/",args,".indiv_var.tsv", sep = "")

write.table(x = indiv_var_count,
            file = outfile_indiv, 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE)

outfile_indiv_med <- paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/",args,".indiv_var_med.tsv", sep = "")

write.table(x = indiv_var_med_count,
            file = outfile_indiv_med, 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE)

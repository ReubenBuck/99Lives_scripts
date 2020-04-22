rm(list = ls())

files <- list.files("/mnt/raid/projects/99Lives_analysis/result/snp_counts/")

files_new <- files[grep("new", files)]
files_indv <- files[grep("indiv", files)]
files_tot <- files[grep("total", files)]

new_var <- matrix(0, nrow = 26, ncol = 10)
indv_var <- matrix(0, nrow = 26, ncol = 10)
tot_var <- matrix(0, nrow = 26, ncol = 10)

for(i in 1:length(files_new)){
  new_var <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/",files_new[i] ,sep = ""),
                        sep = "\t", header = FALSE) + new_var
  indv_var <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/",files_indv[i] ,sep = ""),
                         sep = "\t", header = FALSE) + indv_var
  tot_var <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/",files_tot[i] ,sep = ""),
                        sep = "\t", header = FALSE) + tot_var
}



pdf(file = "/mnt/raid/projects/99Lives_analysis/plots/snp_count/snp_count_genome.pdf",
    width = 6, height = 4)

plot(seq(10,by = 10, length.out = nrow(new_var)),
     rowMeans(new_var)/seq(10,by = 10, length.out = nrow(new_var)), 
     type = "n", ylim = c(0,140000),
     xlab = "Samples", ylab = "New SNVs per sample")
for(i in 1:10){
  lines(seq(10,by = 10, length.out = nrow(new_var)),
        new_var[,i]/seq(10,by = 10, length.out = nrow(new_var)),
        col=8)
}
lines(seq(10,by = 10, length.out = nrow(new_var)),
      rowMeans(new_var)/seq(10,by = 10, length.out = nrow(new_var)),
      col=1, lwd = 3)
grid()

dev.off()


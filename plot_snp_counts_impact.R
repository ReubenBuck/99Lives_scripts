rm(list = ls())

files <- list.files("/mnt/raid/projects/99Lives_analysis/result/snp_counts/snp/")



impact <- c("HIGH", "MODERATE", "LOW")
cols <- c("red", "darkorange", "blue")
names(cols) <- impact


layout(matrix(1:3, ncol = 1))
for(m in impact){

files_new <- files[grepl("new", files) & grepl(m, files) & !grepl("med", files)]
files_indv <- files[grepl("indiv", files) & grepl(m, files)]
files_tot <- files[grepl("total", files) & grepl(m, files)]
files_med <- files[grepl("med", files) & grepl(m, files)]

new_var <- matrix(0, nrow = 26, ncol = 10)
indv_var <- matrix(0, nrow = 26, ncol = 10)
tot_var <- matrix(0, nrow = 26, ncol = 10)
med_var <- matrix(0, nrow = 26, ncol = 10)



new_var <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/snp/",files_new ,sep = ""),
                      sep = "\t", header = FALSE)
indv_var <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/snp/",files_indv ,sep = ""),
                       sep = "\t", header = FALSE)
tot_var <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/snp/",files_tot ,sep = ""),
                      sep = "\t", header = FALSE)
med_var <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/snp/",files_med ,sep = ""),
                      sep = "\t", header = FALSE)




#pdf(file = "/mnt/raid/projects/99Lives_analysis/plots/snp_count/snp_count_genome.pdf",
#    width = 6, height = 4)

plot(seq(10,by = 10, length.out = nrow(new_var)),
     rowMeans(med_var)/seq(10,by = 10, length.out = nrow(new_var)), 
     type = "n", ylim = c(0,rev(rowMeans(med_var))[1] * 6),
     xlab = "Samples", ylab = "New SNVs per sample")
for(i in 1:10){
   lines(seq(10,by = 10, length.out = nrow(new_var)),
         new_var[,i]/seq(10,by = 10, length.out = nrow(new_var)),
         col=8)
  #lines(seq(10,by = 10, length.out = nrow(new_var)),
  #      med_var[,i],
  #      col=8)
}
lines(seq(10,by = 10, length.out = nrow(new_var)),
      rowMeans(med_var),
      col=cols[m], lwd = 3)
grid()


}

dev.off()


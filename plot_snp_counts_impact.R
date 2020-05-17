rm(list = ls())

files <- list.files("/mnt/raid/projects/99Lives_analysis/result/snp_counts/snp/")



impact <- c("HIGH", "MODERATE", "LOW")
cols <- c("red", "darkorange", "blue")
names(cols) <- impact
type <- c("snp", "indel", "both")


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
  
  
  
  new_var.snp <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/snp/",files_new ,sep = ""),
                            sep = "\t", header = FALSE)
  indv_var.snp <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/snp/",files_indv ,sep = ""),
                             sep = "\t", header = FALSE)
  tot_var.snp <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/snp/",files_tot ,sep = ""),
                            sep = "\t", header = FALSE)
  med_var.snp <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/snp/",files_med ,sep = ""),
                            sep = "\t", header = FALSE)
  
  
  new_var.indel <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/indel/",files_new ,sep = ""),
                              sep = "\t", header = FALSE)
  indv_var.indel <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/indel/",files_indv ,sep = ""),
                               sep = "\t", header = FALSE)
  tot_var.indel <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/indel/",files_tot ,sep = ""),
                              sep = "\t", header = FALSE)
  med_var.indel <- read.table(paste("/mnt/raid/projects/99Lives_analysis/result/snp_counts/indel/",files_med ,sep = ""),
                              sep = "\t", header = FALSE)
  
  new_var.both <- new_var.snp + new_var.indel
  indv_var.both <- indv_var.snp + indv_var.indel
  tot_var.both <- tot_var.snp + tot_var.indel
  med_var.both <- med_var.snp + med_var.indel
  
  
  
  
  
  
  
  
  for(t in type){
    
    indv_var <- get(paste("indv_var.",t , sep = ""))
    med_var <- get(paste("med_var.",t , sep = ""))
    new_var <- get(paste("new_var.",t , sep = ""))
    tot_var <- get(paste("tot_var.",t , sep = ""))
    
    
    pdf(file = paste("/mnt/raid/projects/99Lives_analysis/plots/snp_count/mean_",
                     t,"_",m,"_count_genome.pdf", sep = ""),
        width = 6, height = 4)
    
    par(mar= c(5,5,5,5))
    
    
    mean_indiv <- mean(as.matrix(indv_var))
    
    
    plot(seq(10,by = 10, length.out = nrow(new_var)),
         rowMeans(new_var)/seq(10,by = 10, length.out = nrow(new_var))/mean_indiv * 100, 
         type = "n", ylim = c(0,1),
         xlab = "Cohort size",
         ylab = "Mean unique sites per sample (%)", las = 1)
    
    
    for(i in 1:10){
      lines(seq(10,by = 10, length.out = nrow(new_var)),
            new_var[,i]/seq(10,by = 10, length.out = nrow(new_var)) / indv_var[,i] * 100,
            col=8)
    }
    lines(seq(10,by = 10, length.out = nrow(new_var)),
          rowMeans(new_var)/seq(10,by = 10, length.out = nrow(new_var)) / mean_indiv * 100,
          col=cols[m], lwd = 3)
    
    
    grid()
    
    par(new = TRUE)
    plot(seq(10,by = 10, length.out = nrow(new_var)),
         rowMeans(new_var)/seq(10,by = 10, length.out = nrow(new_var))/mean_indiv * 100, 
         type = "n", ylim = c(0,mean_indiv/100), axes = FALSE, xlab = "", ylab = "")
    
    
    ticks <- axTicks(2)
    den <- 10^median(nchar(round(ticks))) / 10
    if(den >= 100000){
      labels <- sapply(ticks/den, function(i,ex) as.expression(bquote(.(i) %*% 10^.(ex))), ex = log10(den))
    }else{
      labels = ticks
    }
    axis(side = 4, at=ticks, labels=labels, las = 2)
    mtext(side = 4, text = "Mean unique sites per sample (n)", line = 3.75)
    
    dev.off()
    
    
    
    
    
    
    pdf(file = paste("/mnt/raid/projects/99Lives_analysis/plots/snp_count/median_",
                     t,"_",m,"_count_genome.pdf", sep = ""),
        width = 6, height = 4)
    
    par(mar = c(5,5,5,5))
    
    mean_indiv <- rowMeans(as.matrix(indv_var))
    
    plot(seq(10,by = 10, length.out = nrow(new_var)),
         rowMeans(med_var)/mean_indiv * 100, 
         type = "n", ylim = c(0,1),
         xlab = "Cohort size",
         ylab = "Median unique sites per sample (%)", las = 1)
    
    
    for(i in 1:10){
      lines(seq(10,by = 10, length.out = nrow(new_var)),
            (med_var[,i])/indv_var[,i] *100,
            col=8)
    }
    lines(seq(10,by = 10, length.out = nrow(new_var)),
          rowMeans(med_var)/mean_indiv * 100,
          col=cols[m], lwd = 3)
    grid()
    
    par(new = TRUE)
    plot(seq(10,by = 10, length.out = nrow(new_var)),
         rowMeans(new_var)/seq(10,by = 10, length.out = nrow(new_var))/mean_indiv * 100, 
         type = "n", ylim = c(0,mean(mean_indiv * 0.01)), axes = FALSE, xlab = "", ylab = "")
    
    
    ticks <- axTicks(2)
    den <- 10^median(nchar(round(ticks))) / 10
    if(den >= 100000){
      labels <- sapply(ticks/den, function(i,ex) as.expression(bquote(.(i) %*% 10^.(ex))), ex = log10(den))
    }else{
      labels = ticks
    }
    axis(side = 4, at=ticks, labels=labels, las = 2)
    mtext(side = 4, text = "Median unique sites per sample (n)", line = 3.75)
    
    
    dev.off()
    
    
    
    pdf(file = paste("/mnt/raid/projects/99Lives_analysis/plots/snp_count/total_",
                     t,"_",m,"_count_genome.pdf", sep = ""),
        width = 6, height = 4)
    
    par(mar = c(5,7,5,3))
    
    
    plot(seq(10,by = 10, length.out = nrow(tot_var)),
         rowMeans(tot_var), 
         type = "n", yaxt = "n",
         xlab = "Cohort size",
         ylab = "", las = 1)
    
    
    for(i in 1:10){
      lines(seq(10,by = 10, length.out = nrow(new_var)),
            (tot_var[,i]),
            col=8)
    }
    lines(seq(10,by = 10, length.out = nrow(new_var)),
          rowMeans(tot_var),
          col=cols[m], lwd = 3)
    grid()
    
    ticks <- axTicks(2)
    den <- 10^median(nchar(as.integer(round(ticks)))) / 10
    if(den >= 100000){
      labels <- sapply(ticks/den, function(i,ex) as.expression(bquote(.(i) %*% 10^.(ex))), ex = log10(den))
    }else{
      labels = ticks
    }
    axis(side = 2, at=ticks, labels=labels, las = 2)
    mtext(side = 2, text = "Total sites", line = 4.25)
    
    
    dev.off()
    
    
    
    
    pdf(file = paste("/mnt/raid/projects/99Lives_analysis/plots/snp_count/new_",
                     t,"_",m,"_count_genome.pdf", sep = ""),
        width = 6, height = 4)
    
    par(mar = c(5,7,5,3))
    
    
    plot(c(10, nrow(new_var)*10),
         range(as.matrix(new_var)), 
         type = "n", yaxt = "n",
         xlab = "Cohort size",
         ylab = "", las = 1)
    
    
    for(i in 1:10){
      lines(seq(10,by = 10, length.out = nrow(new_var)),
            (new_var[,i]),
            col=8)
    }
    lines(seq(10,by = 10, length.out = nrow(new_var)),
          rowMeans(new_var),
          col=cols[m], lwd = 3)
    grid()
    
    ticks <- axTicks(2)
    den <- 10^median(nchar(as.integer(round(ticks)))) / 10
    if(den >= 100000){
      labels <- sapply(ticks/den, function(i,ex) as.expression(bquote(.(i) %*% 10^.(ex))), ex = log10(den))
    }else{
      labels = ticks
    }
    axis(side = 2, at=ticks, labels=labels, las = 2)
    mtext(side = 2, text = "Total unique sites", line = 4.25)
    
    
    dev.off()
    
    
    df <- data.frame(samples = (1:nrow(med_var)) * 10,
                     total_site_mean = rowMeans(tot_var),
                     total_site_sd = apply(tot_var, MARGIN = 1, FUN = sd),
                     total_unique_site_mean = rowMeans(new_var),
                     total_unique_site_sd = apply(new_var, MARGIN = 1, FUN = sd),
                     sites_per_indv_mean = rowMeans(indv_var),
                     sites_per_indv_sd = apply(indv_var, MARGIN = 1, FUN = sd),
                     mean_unique_sites_mean = rowMeans(new_var/(1:nrow(med_var)) * 10),
                     mean_unique_sites_sd = apply(new_var/(1:nrow(med_var)) * 10, MARGIN = 1, FUN = sd),
                     median_unique_sites_mean = rowMeans(med_var),
                     median_unique_sites_sd = apply(med_var, MARGIN = 1, FUN = sd)
    )
    
    
    
    
    write.table(df, file = paste("/mnt/raid/projects/99Lives_analysis/plots/snp_count/coding",
                     t,"_",m,"_counts.tsv", sep = ""),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    
    
  }
  
  
  
  
}











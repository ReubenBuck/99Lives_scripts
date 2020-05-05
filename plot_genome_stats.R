rm(list = ls())
options(stringsAsFactors = FALSE)
library(gdsfmt)
library(SNPRelate)



genofile <-snpgdsOpen("~/Desktop/TMP/data.gds")

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.5, autosome.only = FALSE)
str(snpset)
snpset <- snpset[!grepl("_", names(snpset))]
snpset <- snpset[!grepl("chrX", names(snpset))]
str(snpset)
snpset.id <- unlist(unname(snpset))
length(snpset.id)

inbreed <- snpgdsIndInb(genofile, autosome.only = FALSE, snp.id = snpset.id, method = "mle")



chr <- read.table("/mnt/raid/projects/99Lives_analysis/99Lives_scripts/chr.list")[,1]
samples <- read.table("/mnt/raid/projects/99Lives_analysis/accessory_data/samples_filtered.txt")[,1]


file.dir <- "/mnt/raid/projects/99Lives_analysis/result/vcfStats_advanced/snp_filter"

singleton <- doubleton <- rep(0,length(samples))
for(i in chr){
  file <- paste(file.dir, "/", i, ".snp.filter.stats.singletons", sep = "")
  dat <- read.table(file, sep = "\t", header = TRUE)
  singleton0 <- table(factor(dat[dat$SINGLETON.DOUBLETON == "S","INDV"], levels = samples))
  singleton <- singleton + singleton0
  doubleton0 <- table(factor(dat[dat$SINGLETON.DOUBLETON == "D","INDV"], levels = samples))
  doubleton <- doubleton + doubleton0
  print(i)
  }


gt_sum <- data.frame(sample_ID = samples,
                     n.nocall = 0, n.hom.ref = 0, 
                     n.het = 0, n.hom.alt = 0)
for(i in chr){
  file <- paste(file.dir, "/", i, ".snp.filter.stats.gtSum", sep = "")
  dat <- read.table(file, sep = "\t", header = TRUE, comment.char = "")
  gt_sum[,2:ncol(gt_sum)] <- gt_sum[,2:ncol(gt_sum)] + dat[,2:ncol(dat)]
}



fam_info <- read.table("/mnt/raid/projects/99Lives_analysis/result/cat_families/cat_family_id.tsv",
                       sep= "\t", header = TRUE, comment.char = "")

meta_data <- read.table("/mnt/raid/projects/99Lives_analysis/result/PCA_result/meta_data.tsv", 
                        sep = "\t", header = TRUE, comment.char = "")
rownames(meta_data) <- meta_data$smaple_ID

meta_data$singleton <- singleton
meta_data$doubleton <- doubleton
meta_data$unique_site <- singleton + doubleton


meta_data <- data.frame(meta_data,gt_sum[,2:ncol(gt_sum)])
meta_data$sites <- meta_data$n.het + meta_data$n.hom.alt

meta_data$inbreeding <- inbreed$inbreeding
meta_data$inbreed_iter <- inbreed$out.num.iter

meta_data <- data.frame(meta_data, fam_info[,2:ncol(fam_info)])




pdf(file = "~/Desktop/stats.pdf",width = 10, height = 4,onefile = TRUE)
par(mar = c(5,5,5,1))

plots <-  c("singleton", "sites", "inbreeding")
trans <- c(TRUE, FALSE, FALSE)
names(trans) <- plots
ylab <- c("Singletons per sample", "Sites per sample", "Inbreeding coefficient (F)")
names(ylab) <- plots
expon <- c(FALSE,TRUE,FALSE)
names(expon) <- plots

for(i in plots){

meta_data2 <- meta_data[order(meta_data$cluster,meta_data[,i]),]

cluster_names <- c(1:(length(unique(meta_data2$cluster_col))-1), "Outlier")
names(cluster_names) <- c(unique(meta_data2$cluster)[grep("G", unique(meta_data2$cluster))],
                          unique(meta_data2$cluster)[median(grep("O", unique(meta_data2$cluster)))])
cluster_position <- aggregate(1:nrow(meta_data2), by = list(meta_data2$cluster), median)[,2]
names(cluster_position) <- unique(meta_data2$cluster)
cluster_position <- cluster_position[names(cluster_names)]

meta_data2$cluster_col <- factor(meta_data2$cluster_col,levels = unique(meta_data2$cluster_col))


x <- meta_data2[,i]
if(trans[i]){
  x <- log10(x)
}
barplot(x,
        col = scales::alpha(as.character(meta_data2$cluster_col), .5), 
        border = FALSE, yaxt= "n",
        space = 0, xlab = "Samples",
        ylab = ylab[i])

points(x = seq(.5,nrow(meta_data2),1)[meta_data2$col != "white"],
       y = rep(0,sum(meta_data2$col != "white")),
       pch = 15, cex = .5,
       col = meta_data2$col[meta_data2$col != "white"])


axis(side = 1, at = cluster_position, labels = NA, 
     tick = TRUE,lwd = 1, lwd.ticks = 1)  
axis(side = 1, at = cluster_position[seq(1,length(cluster_position),2)], labels = NA, 
     tick = TRUE,lwd = 0, lwd.ticks = 1, line = .6)  
mtext(text = cluster_names, at = cluster_position, side = 1, 
      line = c(1,.4), col = levels(meta_data2$cluster_col), cex = .9, font = 2)
bndry <- c(0,cumsum(table(meta_data2$cluster_col)))
#abline(v = bndry,lty = 2)

if(trans[i]){
  ticks <- axTicks(2)
  labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  axis(side = 2, at=ticks, labels=labels, las = 2)
}else if(expon[i]){
  ticks <- axTicks(2)
  labels <- sapply(ticks/1e6, function(i) as.expression(bquote(.(i)^ .(6))))
  axis(side = 2, at=ticks, labels=labels, las = 2)
}else{
  axis(side = 2,las = 2)
}

abline(h = axTicks(2), col = "grey90", lty = 3)


box <- boxplot(meta_data2[,i] ~ meta_data2$cluster_col, plot = FALSE)

if(trans[i]){
  box$stats <- log10(box$stats)
}

segments(x0 = bndry[1:(length(bndry)-1)] + .5,
         x1 = bndry[2:(length(bndry))] - .5,
         y0 = box$stats[3,],
         y1 = box$stats[3,],
         col = levels(meta_data2$cluster_col),
         lwd = 3)
segments(x0 = bndry[1:(length(bndry)-1)] + .5,
         x1 = bndry[2:(length(bndry))] - .5,
         y0 = box$stats[2,],
         y1 = box$stats[2,],
         col = levels(meta_data2$cluster_col),
         lwd = 1, lty = 2)
segments(x0 = bndry[1:(length(bndry)-1)] + .5,
         x1 = bndry[2:(length(bndry))] - .5,
         y0 = box$stats[4,],
         y1 = box$stats[4,],
         col = levels(meta_data2$cluster_col),
         lwd = 1, lty = 2)

}


dev.off()
# we can look at the shift between publications!


# put together these facts in a way to make sense


### what other things can we plot?



aov_res <- aov(meta_data$doubleton ~ meta_data$cluster)
summary(aov_res)
tuk_res <- TukeyHSD(aov_res)

barplot(-log10(sort(tuk_res$`meta_data$cluster`[,4])))
abline(h = -log10(.05), col = 2)


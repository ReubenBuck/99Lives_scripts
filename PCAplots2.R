rm(list = ls())
options(stringsAsFactors = FALSE)
library(gdsfmt)
library(SNPRelate)
library(RColorBrewer)
library(dendextend)

#system("rm ~/Desktop/TMP/data.gds")
vcf.fn <- "/mnt/raid/projects/99Lives_analysis/vcf_200117/thin_snps/combined.select.chip.vcf.gz"
snpgdsVCF2GDS(vcf.fn,"~/Desktop/TMP/data.select.gds",method ="biallelic.only")


genofile <-snpgdsOpen("~/Desktop/TMP/data.select.gds")
set.seed(100)


fam_info <- read.table("/mnt/raid/projects/99Lives_analysis/result/cat_families/cat_family_id.tsv",
                       sep= "\t", header = TRUE, comment.char = "")


snpset <- snpgdsLDpruning(genofile, ld.threshold=0.5, autosome.only = FALSE, num.thread = 10)
str(snpset)
snpset <- snpset[!grepl("_", names(snpset))]
snpset <- snpset[!grepl("chrX", names(snpset))]
str(snpset)
snpset.id <- unlist(unname(snpset))
length(snpset.id)


ibs <- snpgdsIBS(genofile,snp.id = snpset.id, autosome.only=FALSE)
ibs.hc<-snpgdsHCluster(ibs)
rv <- snpgdsCutTree(ibs.hc, outlier.n = 2, label.H = FALSE, label.Z = FALSE, col.list = NA)




cols2 <- rep("grey", sum(grepl("Out", levels(rv$samp.group))))
nb.cols <-  sum(!grepl("Out", levels(rv$samp.group)))
mycolors <- rev(colorRampPalette(brewer.pal(8, "Set1"))(nb.cols))
cols <- c(mycolors, cols2)
names(cols) <- levels(rv$samp.group)



meta_data <- read.table("/mnt/raid/projects/99Lives_analysis/accessory_data/signalment.tsv.csv", sep= "\t", header = TRUE)
rownames(meta_data) <- meta_data$LabID
meta_data <- meta_data[rv$sample.id,]


df.meta <- data.frame(sample_ID = meta_data$LabID,
                 breed = as.character(meta_data$breed_simplified),
                 breed_specifc = meta_data$breed,
                 Institution = meta_data$Institution,
                 Institution_symbol = meta_data$Institution_symbol,
                 Country = meta_data$Country,
                 Continent = meta_data$Continent,
                 cluster = rv$samp.group)



pdf(file = "/mnt/raid/projects/99Lives_analysis/result/PCA_result/breed_member.pdf",
    height = 8,width = 8)

main.names <- 1:sum(!duplicated(cols))
main.names[length(main.names)] <- "Out"


# breeds 


df <- data.frame(breed = as.character(meta_data$breed_simplified),
                 cluster = rv$samp.group)
df$cluster[grep("Out", df$cluster)] <- "Outlier001"
df.tab <- table(df)
df.tab[df.tab == 0] <- NA

cat.no <- table(df$cluster)
cat.no <- c(cat.no, Out = sum(cat.no[grepl("Out",names(cat.no))]))
cat.no <- cat.no[!grepl("Outlier",names(cat.no))]


layout(matrix(1:(sum(!duplicated(cols)) + 1), nrow = 1), 
       widths = c(3,rep(1,(sum(!duplicated(cols)) + 1))))

par(mar = c(5,0,5,0), oma = c(0,5,0,2))
mat.na <- t(df.tab[,1])
mat.na[1,] <- 0
image(mat.na , y = 1:nrow(df.tab),
      yaxt = "n", xaxt = "n", axes = FALSE,
      col = "white", 
      main = "")
mtext("Cluster:",side = 3, font = 2, line = 3, adj = 1)
mtext("Samples:",side = 3, line = 1.5, adj = 1, cex = .7)
mtext("Breeds:",side = 3, line = 0.2, adj = 1, cex = .7)

mtext(text = paste(rev(rownames(df.tab)), ":", sep = ""), 
      side = 4, at = 1:nrow(df.tab), las = 2,adj = 1,
      cex = .7, line = -1.4)

mtext(text = rev(rowSums(df.tab, na.rm = TRUE)), 
      side = 4, at = 1:nrow(df.tab), las = 2,adj = 1,
      cex = .7, line = -.3)


for(i in 1:(sum(!duplicated(cols)))){
  image(t(rev(df.tab[,i])), y = 1:nrow(df.tab),
        yaxt = "n", xaxt = "n",
        main = "", 
        col = scales::alpha(cols[i], seq(.3,1,.1))
  )
  grid(ny = nrow(df.tab), nx = 0)
  mtext(main.names[i],side = 3, col = cols[i], font = 2, line = 3)
  mtext(cat.no[i],side = 3, line = 1.5, cex = .7)
  mtext(sum(!is.na(df.tab[,i])),side = 3, line = 0.2, cex = .7)
}

dev.off()

# institution
pdf(file = "/mnt/raid/projects/99Lives_analysis/result/PCA_result/institute_member.pdf",
    height = 8,width = 8)


df <- data.frame(inst = as.character(meta_data$Institution_symbol),
                 cluster = rv$samp.group)
df$cluster[grep("Out", df$cluster)] <- "Outlier001"
df.tab <- table(df)
df.tab[df.tab == 0] <- NA



layout(matrix(1:(sum(!duplicated(cols)) + 1), nrow = 1), 
       widths = c(3,rep(1,(sum(!duplicated(cols)) + 1))))

par(mar = c(5,0,5,0), oma = c(0,5,0,2))
mat.na <- t(df.tab[,1])
mat.na[1,] <- 0
image(mat.na , y = 1:nrow(df.tab),
      yaxt = "n", xaxt = "n", axes = FALSE,
      col = "white", 
      main = "")
mtext("Cluster:",side = 3, font = 2, line = 3, adj = 1)
mtext("Samples:",side = 3, line = 1.5, adj = 1, cex = .7)
mtext("Institutions:",side = 3, line = 0.2, adj = 1, cex = .7)

mtext(text = paste(rev(rownames(df.tab)), ":", sep = ""), 
      side = 4, at = 1:nrow(df.tab), las = 2,adj = 1,
      cex = .7, line = -1.4)

mtext(text = rev(rowSums(df.tab, na.rm = TRUE)), 
      side = 4, at = 1:nrow(df.tab), las = 2,adj = 1,
      cex = .7, line = -.3)

for(i in 1:(sum(!duplicated(cols)))){
  image(t(rev(df.tab[,i])), y = 1:nrow(df.tab),
        yaxt = "n", xaxt = "n",
        main = "", 
        col = scales::alpha(cols[i], seq(.3,1,.1))
  )
  grid(ny = nrow(df.tab), nx = 0)
  mtext(main.names[i],side = 3, col = cols[i], font = 2, line = 3)
  mtext(cat.no[i],side = 3, line = 1.5, cex = .7)
  mtext(sum(!is.na(df.tab[,i])),side = 3, line = 0.2, cex = .7)
}

dev.off()



df.tab.breed <- table(df.meta[,c("breed","cluster")])
df.tab.breed[df.tab.breed == 0] <- NA

most_common_breed = NULL
for(i in 1:ncol(df.tab)){
  most_common_breed <- c(most_common_breed, names(which.max(df.tab.breed[,i])))
}


df.tab.country <- table(df.meta[,c("Country","cluster")])
df.tab.country[df.tab.country == 0] <- NA

most_common_country = NULL
for(i in 1:ncol(df.tab)){
  most_common_country <- c(most_common_country, names(which.max(df.tab.country[,i])))
}

df.tab.continent <- table(df.meta[df.meta$breed == "Random bred",c("Continent","cluster")])
df.tab.continent[df.tab.continent == 0] <- NA

most_common_continent = rep(NA,ncol(df.tab.continent))
names(most_common_continent) <- colnames(df.tab.continent)
for(i in 1:ncol(df.tab)){
  if(all(is.na(df.tab.continent[,i]))){
    next()
  }
  most_common_continent[i] <- names(which.max(df.tab.continent[,i]))
}

most_common <- most_common_breed

most_common[most_common_breed == "Random bred"] <- paste("Random Bred", 
                                                   most_common_continent[most_common_breed == "Random bred"],
                                                   sep = "/")

names(most_common) <- colnames(df.tab.breed)


df.stats <- data.frame(df.meta,
                       most_common_breed = most_common[df$cluster],
                       cluster_col = cols[df.meta$cluster])


pdf(file = "/mnt/raid/projects/99Lives_analysis/result/PCA_result/dendrogram.pdf",
    height = 9, width = 8)
layout(1)
par(mar = c(5,5,5,10))
plot(rv$dendrogram, horiz = TRUE, edge.root = TRUE, leaflab = "none", 
     type = "triangle", xlim = c(.35,0),yaxt= "n", xlab = "Height", ylim = c(269,1))
axis(side = 1, at = seq(0,.3, by = .05))
df.lab <- data.frame(cluster = rv$samp.group[rv$samp.order], order = 1:269)
agg.lab <- aggregate(df.lab$order, by = list(df.lab$cluster), FUN = median)
agg.lab <- agg.lab[grep("G", x = agg.lab$Group.1),]
mtext(text = substr(agg.lab$Group.1, 3,4), side = 4, at = agg.lab$x, las = 2, 
      col = cols[agg.lab$Group.1], cex = 0.7)
mtext(text = most_common[agg.lab$Group.1], side = 4, at = agg.lab$x, las = 2, 
      col = cols[agg.lab$Group.1], line = 2, cex = 0.7)

mtext(text = "Cluster", side = 4, at = -4, las = 2, cex = 0.7,padj = 0, line = -.5)
mtext(text = "Most common breed", side = 4, at = -4, las = 2, cex = 0.7,padj = 0, line = 2)

agg.box <- aggregate(df.lab$order, by = list(df.lab$cluster), FUN = range)
agg.box <- agg.box[grep("G", x = agg.box$Group.1),]
rect(xleft = .3, xright = 0,
     ybottom = agg.box$x[,1], ytop = agg.box$x[,2], border = FALSE,
     col = scales::alpha(cols[agg.box$Group.1],.2))



ins_ID <- unique(df.meta$Institution_symbol)




fam_ids <- unique(fam_info$fam.id)
fam_ids <- fam_ids[fam_ids > 0]
for(i in fam_ids){
  fam_info0 <- fam_info[fam_info$fam.id == i,]
  df.ord <- data.frame(samp = rv$sample.id[rv$samp.order], order = 1:269)
  
  b_heights <- get_leaves_attr(rv$dendrogram, attribute = "height")
  
  y = df.ord$order[df.ord$samp %in% fam_info0$sample.ID]
  x = b_heights[df.ord$samp %in% fam_info0$sample.ID]
  points(x,y, cex = 1, pch = 16, col = fam_info0$col[1])

}

fam_col <- unique(fam_info[,c("fam.id","col")])
fam_col <- fam_col[!(fam_col$col == "black" | fam_col$col == "white"),]
legend("topleft",pch = 16 ,col = fam_col$col,
       legend = fam_col$fam.id, cex = .7, 
       title = "Family ID", bty = "n",pt.cex = 1)

dev.off()



# colour cats according to breed


pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = FALSE)

pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  stringsAsFactors = FALSE)
head(tab)

pdf(file = "/mnt/raid/projects/99Lives_analysis/result/PCA_result/biplot.pdf", width = 10.5, height = 7)
layout(matrix(1:2, nrow = 1), widths = c(10,7))
par(mar = c(5,4,3,0))
plot(tab$EV1, tab$EV2, 
     xlab=paste("Eigenvector 1", " (",signif(pc.percent[1],3)," %)", sep = ""),  
     ylab=paste("Eigenvector 2", " (",signif(pc.percent[2],3)," %)", sep = ""), 
     col = scales::alpha(cols[rv$samp.group], .7), pch = 16)

cluster_names <- c(1:sum(grepl("G", unique(rv$samp.group))), "Out")
cluster_breed <- c(most_common[grepl("G", names(table(rv$samp.group)))], "Outliers")
cluster_col <- cols[!duplicated(cols)]

par(mar = c(5,0,3,2))
plot.new()
legend("topleft", legend = paste(cluster_names,cluster_breed,sep= ", "), 
       pch = 16, col = cluster_col, cex = 1, bty = "n",
       title = "Cluster, Most common breed")

dev.off()


write.table(df.stats, "/mnt/raid/projects/99Lives_analysis/result/PCA_result/meta_data.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)





fst.clust <- snpgdsFst(genofile,population = factor(rv$samp.group),snp.id = snpset.id ,autosome.only = FALSE,
                 method = "W&H02")

fst.breed <- snpgdsFst(genofile,population = factor(df.stats$breed),snp.id = snpset.id ,autosome.only = FALSE,
                 method = "W&H02")

fst.breed <- snpgdsFst(genofile,population = factor(df.stats$breed_specifc[df.stats$breed != "Random bred" & df.stats$breed != "Cross-breed colony"]),
                       snp.id = snpset.id ,autosome.only = FALSE,
                       method = "W&H02", sample.id = df.stats$sample_ID[df.stats$breed != "Random bred" & df.stats$breed != "Cross-breed colony"])


fst.institution <- snpgdsFst(genofile,population = factor(df.stats$Institution_symbol),snp.id = snpset.id ,autosome.only = FALSE,
                       method = "W&H02")

fst.devon <- snpgdsFst(genofile,population = factor(df.stats$cluster == "G013"),snp.id = snpset.id ,autosome.only = FALSE,
                       method = "W&H02")

plot(fst.devon$FstSNP, pch = 16, cex = .5, 
     col = cols[rep(1:length(snpset), sapply(snpset,length))])


fst.persian <- snpgdsFst(genofile,population = factor(df.stats$cluster == "G010"),
                         snp.id = snpset.id,autosome.only = FALSE, method = "W&H02")

plot(fst.persian$FstSNP, pch = 16, cex = .5, 
     col = cols[rep(1:length(snpset), sapply(snpset,length))])



fst.clust$MeanFst
fst.breed$MeanFst
fst.institution$MeanFst


hist(fst.clust$FstSNP)
hist(fst.breed$FstSNP, add = TRUE)







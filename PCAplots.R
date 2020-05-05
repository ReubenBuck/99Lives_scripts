rm(list = ls())
options(stringsAsFactors = FALSE)
library(vcfR)
library(RColorBrewer)

# get genotypes
vcf <- read.vcfR("/mnt/raid/projects/99Lives_analysis/vcf_200117/thin_snps/combined.thin_10kb.vcf.gz")
gt <- extract.gt(vcf)

# get breed info
meta_data <- read.table("Documents/projects/project_files/WGS_database/sample_table.tsv", sep= "\t", header = TRUE)
rownames(meta_data) <- meta_data$LabID
meta_data <- meta_data[colnames(gt),]


# convert genotypes to numbers
gt.num <- gt
gt.num[gt.num == "0/0"] <- 0
gt.num[gt.num == "0/1"] <- 1
gt.num[gt.num == "1/1"] <- 2
gt.num <- data.frame(gt.num)
for(i in 1:ncol(gt.num)){
  gt.num[,i] <- as.integer(gt.num[,i])
}
gt.num <- gt.num - 1
gt.num.cc <- gt.num[complete.cases(gt.num),]

# perform pca
pca.cc <- prcomp(t(gt.num.cc),center = TRUE)



# work out how many dimensions will be used
pdf(file = "/mnt/raid/projects/99Lives_analysis/result/PCA_result/param_select.pdf", 
    onefile = TRUE, height = 4, width = 6)

#layout(matrix(1:2,ncol =1))
pca.cc.sum <- summary(pca.cc)
dim.choice <- 10

par(mar = c(5,5,2,5))
plot(1:269,pca.cc.sum$importance[2,], pch = 16, type="b", 
     xlim = c(0,90), 
     xlab = "Principal component", 
     ylab = "Proportion of variance")
grid()
text(dim.choice, pca.cc.sum$importance[2,][dim.choice], labels =dim.choice, pos = 1, col = 2)
points(dim.choice, pca.cc.sum$importance[2,][dim.choice], pch = 16, col = 2)

plot(1:269,pca.cc.sum$importance[3,], pch = 16, type="b", 
     xlim = c(0,90), ylim = c(0,1),
     xlab = "Principal component", 
     ylab = "Cumulative variance")
grid()
text(dim.choice, pca.cc.sum$importance[3,][dim.choice], dim.choice, pos = 1, col = 2)
points(dim.choice, pca.cc.sum$importance[3,][dim.choice], pch = 16, col = 2)


par(mar = c(5,5,2,5))
# decide on the number of clusters
k.explained <- NULL
for(i in 1:40){
  k.clust <- kmeans(pca.cc$x[,1:dim.choice], centers = i)
  k.explained <- c(k.explained,k.clust$betweenss/k.clust$totss)
}
plot(1:40,k.explained,xlab = "k", 
     ylab = "Variance explained\n(Between SS / Total SS)",
     type = "b", pch = 16
     )
#k.choice <- identify(1:100,k.explained)
k.choice = 11
text(k.choice, k.explained[k.choice], k.choice, pos = 3, col = 2)
points(k.choice, k.explained[k.choice], pch = 16, col = 2)
dev.off()

k.clust <- kmeans(pca.cc$x[,1:dim.choice], centers = k.choice)
k.clust <- meanShift(pca.cc$x[,1:2])


h.clust <- hclust(dist(pca.cc$x[,1:dim.choice]), method = "average")
h.clust <- hclust(dist(gt.num.cc), method = "average")
plot(as.dendrogram(h.clust))
h.cut <- cutree(h.clust, k = 11)
plot(as.dendrogram(h.clust), col = h.cut)


cols <- c(brewer.pal(8, "Set2"),brewer.pal(8, "Dark2"))
plot(pca.cc$x[,c(1,2)],
     col = scales::alpha(cols[h.cut], .7), 
     pch = 16)

text(agg[,c("PC1","PC2")],
     labels = 1:k.choice, col= 1, cex = 2.3, font = 2 
     #pos = sample(1:4, size = k.choice, replace = TRUE)
)
# we can use some sort of iterativer k means to choose the number of clusters, but ultimately use dendrogram to cut.

# at each cut height we can look at within ss and between ss


agg <- aggregate(pca.cc$x[,1:dim.choice], median, by = list(h.cut))






# plot results
layout(1)

pdf(file = "/mnt/raid/projects/99Lives_analysis/result/PCA_result/biplot.pdf")
par(mar = c(5,5,5,2))
cols <- c(brewer.pal(8, "Set2"),brewer.pal(8, "Dark2"))
plot(pca.cc$x[,c(1,2)],
     col = scales::alpha(cols[h.cut], .7), 
     pch = 16)

text(k.clust$centers[,1:2],
     labels = 1:k.choice, col= 1, cex = 2.3, font = 2 
     #pos = sample(1:4, size = k.choice, replace = TRUE)
)


df <- data.frame(as.character(meta_data$breed), k.clust$cluster)
df.tab <- table(df)
df.tab[df.tab == 0] <- NA
dev.off()


pdf(file = "/mnt/raid/projects/99Lives_analysis/result/PCA_result/breed_cluster.pdf")
layout(matrix(1:(k.choice + 1), nrow = 1), widths = c(3,rep(1,(k.choice + 1))))
par(mar = c(5,0,5,0), oma = c(0,0,0,2))
mat.na <- t(df.tab[,1])
mat.na[1,] <- 0
image(mat.na , y = 1:nrow(df.tab),
      yaxt = "n", xaxt = "n", axes = FALSE,
      col = "white", 
      main = "Cluster:")
mtext(text = rev(rownames(df.tab)), 
      side = 4, at = 1:nrow(df.tab), las = 2,adj = 1,
      cex = .6, line = -.3)
for(i in 1:k.choice){
  image(t(rev(df.tab[,i])), y = 1:nrow(df.tab),
        yaxt = "n", xaxt = "n",
        main = i, col = scales::alpha(cols[i], seq(.3,1,.1))
        )
  grid(ny = nrow(df.tab), nx = 0)
}
dev.off()


most_common = NULL
for(i in 1:ncol(df.tab)){
  most_common <- c(most_common, names(which.max(df.tab[,i])))
}


df.stats <- data.frame(smaple_ID = rownames(df), 
                 breed = df$as.character.meta_data.breed., 
                 cluster = df$k.clust.cluster,
                 most_common_breed = most_common[df$k.clust.cluster],
                 meta_data[c("biomaterial_provider","collected_by","disease","disease_stage")])

write.table(df.stats, "/mnt/raid/projects/99Lives_analysis/result/PCA_result/PCA_metadata.tsv",
            sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

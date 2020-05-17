
rm(list = ls())
library(gdsfmt)
library(SNPRelate)
library(RColorBrewer)


#system("rm ~/Desktop/TMP/data.gds")
#vcf.fn <- "/mnt/raid/projects/99Lives_analysis/vcf_200117/thin_snps/combined.thin_1kb.1maf.vcf.gz"
#snpgdsVCF2GDS(vcf.fn,"~/Desktop/TMP/data.1maf.gds",method ="biallelic.only")



genofile <-snpgdsOpen("~/Desktop/TMP/data.1maf.gds")
set.seed(100)




snpset <- snpgdsLDpruning(genofile, ld.threshold=0.5, autosome.only = FALSE, num.thread = 5, method = "corr")
str(snpset)
snpset <- snpset[!grepl("_", names(snpset))]
snpset <- snpset[!grepl("chrX", names(snpset))]
str(snpset)
snpset.id <- unlist(unname(snpset))
length(snpset.id)
# 319815 markers after LD pruning at 0.5 with method corr


ibd.robust <- snpgdsIBDKING(genofile, snp.id = snpset.id, autosome.only = FALSE)
ibd.mat <- ibd.robust$kinship
hist(ibd.mat, breaks = 1000, xlim = c(-.1,.5))
ibd.mat[ibd.mat < 0] <- NA
layout(1)
par(mar = c(5,5,5,5))
hist(ibd.mat, breaks = 100)


df.fam <- data.frame(row = row(ibd.mat)[!is.na(ibd.mat) & ibd.mat < .5 & ibd.mat > .125 - (.065/2)],
                     col = col(ibd.mat)[!is.na(ibd.mat) & ibd.mat < .5 & ibd.mat > .125 - (.065/2)])
length(unique(df.fam$row))

# maybe we can pull out family clusters using single linkage clustering.

df.fam0 <- df.fam[2:nrow(df.fam),]
df.fam1 <- df.fam[1,]
df.fam.list <- NULL

while(nrow(df.fam0) > 0 ){
  
  # put matches in match pile
  df.fam1 <- rbind(df.fam1, df.fam0[df.fam0$row %in% unlist(df.fam1),])
  # then remove them
  df.fam0 <- df.fam0[!(df.fam0$row %in% unlist(df.fam1)),]
  
  # put matches in match pile
  df.fam1 <- rbind(df.fam1, df.fam0[df.fam0$col %in% unlist(df.fam1),])
  # then remove them
  df.fam0 <- df.fam0[!(df.fam0$col %in% unlist(df.fam1)),]
  
  
  if(sum(unlist(df.fam0) %in% unlist(df.fam1)) == 0){
    df.fam.list <- c(df.fam.list, list(df.fam1))
    df.fam1 <- df.fam0[1,]
    df.fam0 <- df.fam0[2:nrow(df.fam0),]
    df.fam0 <- df.fam0[complete.cases(df.fam0),]
    print(nrow(df.fam0))
    
  }
}

df.fam.list
fam.list <- lapply(df.fam.list, function(x){unique(unlist(x))})

# now we have all the cat families!!!!

pdf(file = "/mnt/raid/projects/99Lives_analysis/result/cat_families/cat_families.pdf",
    height = 7.5, width = 7.5)

hist(ibd.robust$kinship, breaks = 100, xlim = c(-.1,.5),
     xlab = "Kinship estimate (Phi)", main  = "Relatedness in 99 Lives",
     yaxs = "i")
abline(v = .125 - (.065/2), col = 2, lwd = 2)

colRamp <- colorRampPalette(c("white", "red"))
par(mar = c(0,0,0,0), oma = c(5,5,5,5))

for(i in 1:length(fam.list)){
  
  
  if( length(fam.list[[i]]) > 2 ){
    
    
    layout(matrix(1:4,byrow = TRUE, nrow = 2),
           widths = c(20 - length(fam.list[[i]]), 
                      length(fam.list[[i]])),
           heights  = c(length(fam.list[[i]]),
                        20 - length(fam.list[[i]]))
    )
    
    
    plot.new()
    
    image(ibd.mat[fam.list[[i]], fam.list[[i]]],
          zlim = c(0,.5), col = colRamp(20),
          xaxt = "n", yaxt = "n")
    
    axis(side = 1, 
         at = seq(0,1, length.out = length(fam.list[[i]])),
         ibd.robust$sample.id[fam.list[[i]]], las = 2)
    
    axis(side = 2, 
         at = seq(0,1, length.out = length(fam.list[[i]])),
         ibd.robust$sample.id[fam.list[[i]]], las = 2)
    
    title(main = paste("Family", i), outer = TRUE)
  }
}

par(mar = c(0,0,0,0))
layout(matrix(1:9,byrow = TRUE, nrow = 3),
       widths = c(1,.2,1),
       heights  = c(1,1,1))
plot.new()
plot.new()
plot.new()
plot.new()
par(mar = c(0,0,5,0))
image(t(matrix(seq(0,.5,.1))), col = colRamp(20),
      xaxt = "n", yaxt = "n", main = "Kinship\nestimate\n(Phi)")
axis(side = 2,at = seq(0,1,length.out = 6), 
     labels = seq(0,.5,.1), las = 2)
dev.off()
# now we have heat maps for all the cat families.


# write a table of all cat relationships that can be used else where
# may be we could add a second mark to the histogram for the cat family ids.

cat_fam_id <- data.frame(sample.ID = ibd.robust$sample.id,
                         fam.id = 0,
                         fam.pair = FALSE)
for(i in 1:length(fam.list)){
  cat_fam_id$fam.id[fam.list[[i]]] <- i
  if(length(fam.list[[i]]) == 2){
    cat_fam_id$fam.pair[fam.list[[i]]] <- TRUE
  }
}


non_pair_families <- unique(cat_fam_id[!cat_fam_id$fam.pair & cat_fam_id$fam.id > 0,"fam.id"])
pair_families <- unique(cat_fam_id[cat_fam_id$fam.pair & cat_fam_id$fam.id > 0,"fam.id"])

cols <- RColorBrewer::brewer.pal(5,"Set1")
cols <- colorRampPalette(colors = cols)
cols1 <- cols(length(non_pair_families))
cols2 <- rep("black",length(pair_families))
colsall <- c(cols1,cols2)
colsall <- colsall[order(c(non_pair_families,pair_families))]

cat_fam_id$col <- "white"
cat_fam_id$col[cat_fam_id$fam.id > 0] <- colsall[cat_fam_id$fam.id]

cat_fam_id_col <- cat_fam_id$fam.id
names(cat_fam_id_col) <- cat_fam_id$sample.ID

kinship.col <- data.frame(row = cat_fam_id_col[ibd.robust$sample.id[row(ibd.robust$kinship)]],
                          col = cat_fam_id_col[ibd.robust$sample.id[col(ibd.robust$kinship)]])
kinship.col[kinship.col$row != kinship.col$col,] <- NA
kinship.col[kinship.col==0] <- NA

cat_fam_id_col <- cat_fam_id$col
names(cat_fam_id_col) <- cat_fam_id$sample.ID
kinship.col$row[!is.na(kinship.col$row)] <- colsall[kinship.col$row[!is.na(kinship.col$row)]]
kinship.col <- kinship.col$row


pdf(file = "/mnt/raid/projects/99Lives_analysis/result/cat_families/phi_ibs.pdf")
plot(ibd.robust$IBS0, ibd.robust$kinship, 
     pch = 16, col = scales::alpha(as.character(kinship.col), .5),
     xlim = c(0,0.05), ylim = c(-.15,.5),
     xlab = "Proportion of Zero IBS",
     ylab = "Estimated Kinship Coefficient",
     main = "Familial relationships identified in 99 Lives")
abline(h = .25 + (.1), lty = 2)
abline(h = .25 - (.125/2), lty = 2)
abline(h = .125 - (.065/2), lty = 2)
abline(h = .065 - (.0325/3), lty = 2)
abline(h = 0, lty = 2, lwd = 3)
legend("topright",pch = 16 ,col = colsall[colsall!="black"],
       legend = (1:length(colsall))[colsall!="black"], cex = 1, 
       title = "Family ID", bg = "white")
dev.off()

write.table(x = cat_fam_id,
            "/mnt/raid/projects/99Lives_analysis/result/cat_families/cat_family_id.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


write.table(x = ibd.robust$kinship,
            "/mnt/raid/projects/99Lives_analysis/result/cat_families/cat_kinship.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)





snpgdsClose(genofile)




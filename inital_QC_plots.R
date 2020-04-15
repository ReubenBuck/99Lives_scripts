rm(list = ls())
options(stringsAsFactors = FALSE)

dat <- read.table("/mnt/raid/projects/99Lives_analysis/result/first_pass/chrE3_vcf_200117.gt_sum.tsv",
                  sep = "\t", comment.char = "", header = TRUE)
dep <- read.table("/mnt/raid/projects/99Lives_analysis/result/first_pass/chrE3_vcf_200117.idepth",
                  sep = "\t", comment.char = "", header = TRUE)

# if we do a coverage and missing call 

nocall <- dat$n.nocall
names(nocall) <- dat$X.sample.id

cov <- dep$MEAN_DEPTH
names(cov) <- dep$INDV

nocall_perc <- (nocall/ (dat$n.nocall + dat$n.hom.ref + dat$n.het + dat$n.hom.alt))


pdf(file = "/mnt/raid/projects/99Lives_analysis/result/first_pass/qulaity_control_plots.pdf", onefile = TRUE)
plot(cov, 
     nocall_perc * 100,
     ylim =c(0,10),
     xlim = c(0,50),
     ylab = "sites with missing genotype (%)",
     xlab = "mean site depth",
     main = "missing genotypes")
grid()



plot(cov,dat$n.het/dat$n.hom.alt,
     ylim = c(0,5),
     xlim = c(0,50),
     ylab = "het/hom ratio",
     xlab = "mean site depth",
     main = "het/hom ratio")
grid()


plot(cov,dat$n.het+dat$n.hom.alt,
     xlim = c(0,50),
     ylab = "sites (n)",
     xlab = "mean site depth",
     main = "Variant sites")
grid()


dev.off()


names(cov)[cov > 30 & dat$n.het/dat$n.hom.alt < .5]








a <- identify(cov,dat$n.het/dat$n.hom.alt)
names(cov)[a]



cov[dat$n.het/dat$n.hom.alt > 3]



abline(h = 2)
abline(v = 10)

sum(nocall_perc < .025)
sum(cov > 15)
sum(cov > 10)
sum(cov > 20)

sum(nocall_perc < .025 & cov > 15)



plot(cov, (dat$n.hom.ref + dat$n.het + dat$n.hom.alt))

nocall_perc[cov >= 10]
cov["felCat.Zivziv.Zivziv"]


names(cov)[cov < 10]



options(stringsAsFactors = FALSE)
samp <- read.table("/mnt/raid/projects/99Lives_analysis/vcf_200117/combined_file/samples.txt")[,1]

samp_rm <- c("felCat.CVB13769.20679", 
             "felCat.Fcat20796.BIRA003CT189",
             "felCat.SAMEA104694019.K510",
             "felCat.SAMEA1061687.K499",
             "felCat.SAMEA1061689.K515",
             "felCat.SAMEA1061799.K511",
             "felCat.SAMEA2785958.K450",
             "felCat.Zivziv.Zivziv",
             "felCat.cat_07.Tesla",
             "felCat.cat_09.Moses",
             "felCat.cat_10.Louie",
             "felCat.cat_17.Sundance",
             "felCat.cat_24.Gumbo",
             "felCat.cat_33.D05-0399",
             "felCat.Fcat12682.Cinnamon")

samp <- samp[!(samp %in% samp_rm)]

dat <- read.table("Documents/projects/project_files/WGS_database/sample_table.tsv", sep= "\t", header = TRUE)

rownames(dat) <- dat$LabID


samp_dat <- dat[samp,]
no_samp_dat <- dat[!(rownames(dat) %in% samp),]

write.table(samp_dat,"/mnt/raid/projects/99Lives_analysis/result/signalment/kept_samples.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(no_samp_dat,"/mnt/raid/projects/99Lives_analysis/result/signalment/removed_samples.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)




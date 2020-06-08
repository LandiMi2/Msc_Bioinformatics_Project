setwd("~/Documents/IMGT/Simulated_bovine_IMGT")
#read repertpire vdj recombination file 
rep_vdj <- read.table("rep_vdj_recombination.txt", header = F, sep = ";")
colnames(rep_vdj) <- c("Antibidy_ID", "V_call", "D_call", "J_call")

View(rep_vdj)
#remove _ on the antibody name 
rep_vdj$Antibidy_ID <- gsub("\\_", "", rep_vdj$Antibidy_ID)

#read imgt annotation
imgt <- read.delim("./4_IMGT-gapped-AA-sequences.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

imgt_specific <- data.frame(imgt$Sequence.ID, imgt$V.GENE.and.allele, imgt$D.GENE.and.allele, imgt$J.GENE.and.allele) 
View(imgt_specific)
View(imgt)


######################################## V call analysis ##########################################################

imgt_Vcall <- data.frame(imgt_specific$imgt.Sequence.ID, imgt_specific$imgt.V.GENE.and.allele)
head(imgt_Vcall)

#process imgt columns before analysis 
#1. remove Bostau, F & P at the end of the gene names

# remove Bostau
imgt_Vcall$imgt_specific.imgt.V.GENE.and.allele <-  gsub("Bostau ", "", imgt_Vcall$imgt_specific.imgt.V.GENE.and.allele)
#consider the first gene in the column 
imgt_Vcall$imgt_specific.imgt.V.GENE.and.allele <-  gsub("^(.*?),.*", "\\1", imgt_Vcall$imgt_specific.imgt.V.GENE.and.allele)
#remove the last character P and F
imgt_Vcall$imgt_specific.imgt.V.GENE.and.allele <-  gsub(' .{1}$', "", imgt_Vcall$imgt_specific.imgt.V.GENE.and.allele)
#rename column names for Vcall dataframe 
colnames(imgt_Vcall) <- c("Antibody_ID", "IMGT_Vcall")



View(imgt_Vcall)

rep_vdj_Vcall <- data.frame(rep_vdj$Antibidy_ID, rep_vdj$V_call)
head(rep_vdj_Vcall)
head(imgt_Vcall)
colnames(rep_vdj_Vcall) <- c("Antibody_ID", "True_Vcall")
#merger to see matches and mismatches

imgt_compare <- merge(data.frame(imgt_Vcall), data.frame(rep_vdj_Vcall), by = "Antibody_ID", all = TRUE)

head(imgt_compare)


#create a hit column
imgt_compare$HIT <- NA
head(imgt_compare)
#rename the column names to be more informative
colnames(imgt_compare) <- c("Antibodty_ID", "IMGT_Vcall", "True_Vcall", "Hit_or_Mishit")

head(imgt_compare)
#where there is a match give TRUE , whereas a mishit give false on the HIT column

imgt_compare$Hit_or_Mishit <- as.character(imgt_compare$IMGT_Vcall) == as.character(imgt_compare$True_Vcall)
View(imgt_compare)

#drop NA rows 
install.packages("tidyr")
library(tidyr)
#removing NA
imgt_compare_without_Na <- imgt_compare %>% drop_na()
#blank rows example 
imgt_compare_without_Na[212,]
#remove blank rows
imgt_compare_wo_NA_Space <- imgt_compare_without_Na[-which(imgt_compare_without_Na$IMGT_Vcall == ""), ]

View(imgt_compare_without_Na)
nrow(imgt_compare_without_Na)

imgt_table_hit <- table(imgt_compare$Hit_or_Mishit)
imgt_table_hit

imgt_table_modified <- table(imgt_compare_wo_NA_Space$Hit_or_Mishit)
imgt_table_modified

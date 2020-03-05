setwd("/home/cofia/Documents/Msc_Bioinformatics_Project/Data/Annotatio_Output")

#mixcr <- read.delim("./MIXCR/Simulated_data/simulated_clones_0.5M.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
igblast <- read.delim("./Igblast/simulated_igblast_merged.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

View(igblast)

#######Read VDJ recombination; we are using this file for comaprision because it contains information about V(D)J recombination for each reads in merged file.
read_vdj <- read.table("read_vdj.txt", header = F, sep = ";")
#rename column names 
colnames(read_vdj) <- c("Antibody_ID", "V_call", "D_call", "J_call")
View(read_vdj)

nrow(read_vdj)#420668
nrow(igblast)#420668

########Take 4 columns from annotation - (seq ID , Vcall, Dcall and J call)
igblast_specific <- data.frame(igblast$sequence_id, igblast$v_call, igblast$d_call, igblast$j_call) 
head(igblast_specific)


############################################### V call analysis for Igblast ############################################################################## 
igblast_Vcall <- data.frame(igblast_specific$igblast.sequence_id, igblast_specific$igblast.v_call)
head(igblast_Vcall)
#rename column name for Igblast Vcall data frame 
colnames(igblast_Vcall) <- c("Antibody_ID", "IgBlast_Vcall")
read_vdj_Vcall <- data.frame(read_vdj$Antibody_ID, read_vdj$V_call)
head(read_vdj_Vcall)
#rename column names for read_vdj_Vcall dataframe
colnames(read_vdj_Vcall) <- c("Antibody_ID", "True_Vcall")


#IgBLAST will give you miltiple genes for the same sequence, when it is undecided.
#e.g. this one
igblast_Vcall[4,]

#it is convention in these cases, to consider only the first one before matching with the true genes

igblast_Vcall$IgBlast_Vcall <- gsub("^(.*?),.*", "\\1", igblast_Vcall$IgBlast_Vcall)
#confirm 
igblast_Vcall[4,]


#merger to see matches and mismatches
comparison_Vcall_Igblast <- merge(data.frame(igblast_Vcall), data.frame(read_vdj_Vcall), 
                                  by = "Antibody_ID" , all = TRUE)

View(comparison_Vcall_Igblast)
#create a hit column
comparison_Vcall_Igblast$HIT <- NA

#where there is a match give TRUE , whereas a mishit give false on the HIT column
comparison_Vcall_Igblast$HIT <- as.character(comparison_Vcall_Igblast$IgBlast_Vcall) == as.character(comparison_Vcall_Igblast$True_Vcall)
head(comparison_Vcall_Igblast)
nrow(comparison_Vcall_Igblast) #420668
#create frequency table of hits and mishits for Vgene
igblast_Vcall_hit_table <- table(comparison_Vcall_Igblast$HIT)
igblast_Vcall_hit_table





#############################################   D  call analysis for Igblast ###################################################################
igblast_Dcall <- data.frame(igblast_specific$igblast.sequence_id, igblast_specific$igblast.d_call)
head(igblast_Dcall)
#rename column name for Igblast Dcall data frame 
colnames(igblast_Dcall) <- c("Antibody_ID", "IgBlast_Dcall")

read_vdj_Dcall <- data.frame(read_vdj$Antibody_ID, read_vdj$D_call)
head(read_vdj_Dcall)
#rename column names for read_vdj_Dcall dataframe
colnames(read_vdj_Dcall) <- c("Antibody_ID", "True_Dcall")


#IgBLAST will give you miltiple genes for the same sequence, when it is undecided.
#e.g. this one
igblast_Dcall[3,]

#it is convention in these cases, to consider only the first one before matching with the true genes

igblast_Dcall$IgBlast_Dcall <- gsub("^(.*?),.*", "\\1", igblast_Dcall$IgBlast_Dcall)
#confirm 
igblast_Dcall[3,]


#merger to see matches and mismatches
comparison_Dcall_Igblast <- merge(data.frame(igblast_Dcall), data.frame(read_vdj_Dcall), 
                                  by = "Antibody_ID" , all = TRUE)

head(comparison_Dcall_Igblast)
View(comparison_Dcall_Igblast)
#create a hit column
comparison_Dcall_Igblast$HIT <- NA

#where there is a match give TRUE , whereas a mishit give false on the HIT column
comparison_Dcall_Igblast$HIT <- as.character(comparison_Dcall_Igblast$IgBlast_Dcall) == as.character(comparison_Dcall_Igblast$True_Dcall)
head(comparison_Dcall_Igblast)

#create frequency table of hits and mishits for Dgene
igblast_Dcall_hit_table <- table(comparison_Dcall_Igblast$HIT)
igblast_Dcall_hit_table

#Notice there are rows where IgBlast doesn't call any allele  e.g.
comparison_Dcall_Igblast[9, ]

#lets remove the blank rows because they don't tell us anything concerning hit or mishit 
comparison_Dcall_Igblast_wo_space <- comparison_Dcall_Igblast[-which(comparison_Dcall_Igblast$IgBlast_Dcall == ""), ]
#confirm

comparison_Dcall_Igblast_wo_space[9, ]

#create a frequency table of hits witout space 
igblast_Dcall_hit_table_wo_space <- table(comparison_Dcall_Igblast_wo_space$HIT)
igblast_Dcall_hit_table_wo_space

#############################################   J  call analysis for Igblast ###################################################################
igblast_Jcall <- data.frame(igblast_specific$igblast.sequence_id, igblast_specific$igblast.j_call)
head(igblast_Jcall)
#rename column name for Igblast Jcall data frame 
colnames(igblast_Jcall) <- c("Antibody_ID", "IgBlast_Jcall")
read_vdj_Jcall <- data.frame(read_vdj$Antibody_ID, read_vdj$J_call)
head(read_vdj_Jcall)
#rename column names for read_vdj_Vcall dataframe
colnames(read_vdj_Jcall) <- c("Antibody_ID", "True_Jcall")


#IgBLAST will give you miltiple genes for the same sequence, when it is undecided.
#e.g. this one
igblast_Jcall[3,]

#it is convention in these cases, to consider only the first one before matching with the true genes

igblast_Jcall$IgBlast_Jcall <- gsub("^(.*?),.*", "\\1", igblast_Jcall$IgBlast_Jcall)
#confirm 
igblast_Jcall[3,]


#merger to see matches and mismatches
comparison_Jcall_Igblast <- merge(data.frame(igblast_Jcall), data.frame(read_vdj_Jcall), 
                                  by = "Antibody_ID" , all = TRUE)

View(comparison_Jcall_Igblast)
#create a hit column
comparison_Jcall_Igblast$HIT <- NA

#where there is a match give TRUE , whereas a mishit give false on the HIT column
comparison_Jcall_Igblast$HIT <- as.character(comparison_Jcall_Igblast$IgBlast_Jcall) == as.character(comparison_Jcall_Igblast$True_Jcall)
head(comparison_Jcall_Igblast)

#create frequency table of hits and mishits for Jgene
igblast_Jcall_hit_table <- table(comparison_Jcall_Igblast$HIT)
igblast_Jcall_hit_table

View(comparison_Jcall_Igblast)

#Notice there are rows where IgBlast doesn't call any allele  e.g.
comparison_Jcall_Igblast[6, ]

#lets remove the blank rows because they don't tell us anything concerning hit or mishit 
comparison_Jcall_Igblast_wo_space <- comparison_Jcall_Igblast[-which(comparison_Jcall_Igblast$IgBlast_Jcall == ""), ]
#confirm

comparison_Jcall_Igblast_wo_space[6, ]

#create a frequency table of hits witout space 
igblast_Jcall_hit_table_wo_space <- table(comparison_Jcall_Igblast_wo_space$HIT)
igblast_Jcall_hit_table_wo_space



########################################### combine V, D, J calls mishit data ##################################################################
igblast_Vcall_hit_table #False 38003  ,  True 382665 
igblast_Dcall_hit_table #False 178229 ,  True 242439 
igblast_Jcall_hit_table #False 203314 ,  True 217354 

#create a martrix for the frequencies 
igblast_hit_mishit <- matrix(c(38003, 382665, 178229, 242439, 203314, 217354), ncol = 2, byrow = TRUE)
colnames(igblast_hit_mishit) <- c("Mishits", "Hits")
rownames(igblast_hit_mishit) <- c("V call", "D call", "J call")
igblast_hit_mishit

#convert the martrix to a table so that we can start to analysis these counts 
igblast_hits_table <- as.table(igblast_hit_mishit)

igblast_hits_table

barplot(igblast_hits_table,legend.text = T,beside = TRUE, col = c("red", "green", "yellow"), 
        main = "Igblast annotation analysis of hits and mishits")

#plot only hits 
library(ggplot2)

ggplot(igblast_hits_df, aes(x= rownames(igblast_hits_df) ,y=igblast_hits_df$Hits)) + geom_bar(stat = "identity")

#for unassigned calls in percentage for D call and J call ..

percent_unassigned_Dcall <- ((nrow(comparison_Dcall_Igblast) - nrow(comparison_Dcall_Igblast_wo_space))) / 
  nrow(comparison_Dcall_Igblast) * 100
# 2.177251
percent_unassigned_Jcall <- ((nrow(comparison_Jcall_Igblast) - nrow(comparison_Jcall_Igblast_wo_space))) / 
  nrow(comparison_Jcall_Igblast) * 100
#26.93288


########################################## Do IMGT analysis ######################################################
#read imgt annotation 
imgt <- read.delim("./IMGT/Bovine_Simulated_Annotation/4_IMGT-gapped-AA-sequences.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

View(imgt)

#read read_recombination file and process the sequence ID to match imgt specifi sequence ID column 

read_vdj_imgt <- read.table("read_vdj.txt", header = F, sep = ";")
head(read_vdj_imgt)
colnames(read_vdj_imgt) <- c("Sequence_ID", "True_Vcall", "True_Dcall", "True_Jcall")

read_vdj_imgt$Sequence_ID <- sub("_multiplicity.*", "", read_vdj_imgt$Sequence_ID)

########Take 4 columns from annotation of imgt - (seq ID , Vcall, Dcall and J call)
imgt_specific <- data.frame(imgt$Sequence.ID, imgt$V.GENE.and.allele, imgt$D.GENE.and.allele, imgt$J.GENE.and.allele) 
View(imgt_specific)
#rename the column names of imgt specific 
colnames(imgt_specific) <- c("Sequence_ID", "IMGT_Vcall", "IMGT_Dcall", "IMGT_Jcall")
head(imgt_specific)
#process Sequence ID column of imgt_specific file;
imgt_specific$Sequence_ID <-  sub("_multiplicity.*", "", imgt_specific$Sequence_ID)


######################################## V call analysis ##########################################################
imgt_Vcall <- data.frame(imgt_specific$Sequence_ID, imgt_specific$IMGT_Vcall)
colnames(imgt_Vcall) <- c("Sequence_ID", "IMGT_Vcall")
head(imgt_Vcall)

#process imgt V call columns before analysis 
#1. remove Bostau, F & P at the end of the gene names

# remove Bostau
imgt_Vcall$IMGT_Vcall <-  gsub("Bostau ", "", imgt_Vcall$IMGT_Vcall)
#consider the first gene in the column 
imgt_Vcall$IMGT_Vcall <-  gsub("^(.*?),.*", "\\1", imgt_Vcall$IMGT_Vcall)
#remove the last character P and F
imgt_Vcall$IMGT_Vcall <-  gsub(' .{1}$', "", imgt_Vcall$IMGT_Vcall)
head(imgt_Vcall)
read_vdj_imgt_Vcall <- data.frame(read_vdj_imgt$Sequence_ID, read_vdj_imgt$True_Vcall)
head(read_vdj_imgt_Vcall)
#merger to see matches and mismatches
comparison_Vcall_imgt <- merge(data.frame(imgt_Vcall), data.frame(read_vdj_imgt_Vcall), by.x = "Sequence_ID",
                          by.y = "read_vdj_imgt.Sequence_ID", all = TRUE)
head(comparison_Vcall_imgt)
#create a hit column
comparison_Vcall_imgt$HIT <- NA
#rename the column names to be more informative
colnames(comparison_Vcall_imgt) <- c("Antibodty_ID", "IMGT_Vcall", "True_Vcall", "HIT")
head(comparison_Vcall_imgt)
#where there is a match give TRUE , whereas a mishit give false on the HIT column
comparison_Vcall_imgt$HIT <- as.character(comparison_Vcall_imgt$IMGT_Vcall) == as.character(comparison_Vcall_imgt$True_Vcall)
head(comparison_Vcall_imgt)
#create frequency table of hits and mishits for Vgene
imgt_Vcall_hit_table <- table(comparison_Vcall_imgt$HIT)
imgt_Vcall_hit_table

########################################## D call analysis #####################################
imgt_Dcall <- data.frame(imgt_specific$Sequence_ID, imgt_specific$IMGT_Dcall)
colnames(imgt_Dcall) <- c("Sequence_ID", "IMGT_Dcall")
head(imgt_Dcall)

#process imgt V call columns before analysis 
#1. remove Bostau, F & P at the end of the gene names

# remove Bostau
imgt_Dcall$IMGT_Dcall <-  gsub("Bostau ", "", imgt_Dcall$IMGT_Dcall)
#consider the first gene in the column 
imgt_Dcall$IMGT_Dcall <-  gsub("^(.*?),.*", "\\1", imgt_Dcall$IMGT_Dcall)
#remove the last character P and F
imgt_Dcall$IMGT_Dcall <-  gsub(' .{1}$', "", imgt_Dcall$IMGT_Dcall)
View(imgt_Dcall)
#remove ORF at the end of some gene calls 
imgt_Dcall$IMGT_Dcall <- gsub(" ORF", "", imgt_Dcall$IMGT_Dcall)

read_vdj_imgt_Dcall <- data.frame(read_vdj_imgt$Sequence_ID, read_vdj_imgt$True_Dcall)
head(read_vdj_imgt_Dcall)
head(imgt_Dcall)
#merger to see matches and mismatches
comparison_Dcall_imgt <- merge(data.frame(imgt_Dcall), data.frame(read_vdj_imgt_Dcall), by.x = "Sequence_ID",
                               by.y = "read_vdj_imgt.Sequence_ID", all = TRUE)
head(comparison_Dcall_imgt)
#create a hit column
comparison_Dcall_imgt$HIT <- NA
#rename the column names to be more informative
colnames(comparison_Dcall_imgt) <- c("Antibodty_ID", "IMGT_Dcall", "True_Dcall", "HIT")
head(comparison_Dcall_imgt)
#where there is a match give TRUE , whereas a mishit give false on the HIT column
comparison_Dcall_imgt$HIT <- as.character(comparison_Dcall_imgt$IMGT_Dcall) == as.character(comparison_Dcall_imgt$True_Dcall)

head(comparison_Dcall_imgt)
#create frequency table of hits and mishits for Vgene
imgt_Dcall_hit_table <- table(comparison_Dcall_imgt$HIT)
imgt_Dcall_hit_table

#Notice there are rows where IMGT doesn't call any allele  e.g.
comparison_Dcall_imgt[6, ]

#lets remove the blank rows because they don't tell us anything concerning hit or mishit 
comparison_Dcall_imgt_wo_space <- comparison_Dcall_imgt[-which(comparison_Dcall_imgt$IMGT_Dcall == ""), ]
#confirm
comparison_Dcall_imgt_wo_space[6, ]

#create a frequency table of hits witout space 
imgt_Dcall_hit_table_wo_space <- table(comparison_Dcall_imgt_wo_space$HIT)
imgt_Dcall_hit_table_wo_space


########################################## J call analysis #####################################
imgt_Jcall <- data.frame(imgt_specific$Sequence_ID, imgt_specific$IMGT_Jcall)
colnames(imgt_Jcall) <- c("Sequence_ID", "IMGT_Jcall")
head(imgt_Jcall)

#process imgt V call columns before analysis 
#1. remove Bostau, F & P at the end of the gene names

# remove Bostau
imgt_Jcall$IMGT_Jcall <-  gsub("Bostau ", "", imgt_Jcall$IMGT_Jcall)
#consider the first gene in the column 
imgt_Jcall$IMGT_Jcall <-  gsub("^(.*?),.*", "\\1", imgt_Jcall$IMGT_Jcall)
#remove the last character P and F
imgt_Jcall$IMGT_Jcall <-  gsub(' .{1}$', "", imgt_Jcall$IMGT_Jcall)
#remove ORF at the end of some gene calls 
imgt_Jcall$IMGT_Jcall <- gsub(" ORF", "", imgt_Jcall$IMGT_Jcall)
head(imgt_Jcall)

read_vdj_imgt_Jcall <- data.frame(read_vdj_imgt$Sequence_ID, read_vdj_imgt$True_Jcall)
head(read_vdj_imgt_Jcall)
head(imgt_Dcall)
#merger to see matches and mismatches
comparison_Jcall_imgt <- merge(data.frame(imgt_Jcall), data.frame(read_vdj_imgt_Jcall), by.x = "Sequence_ID",
                               by.y = "read_vdj_imgt.Sequence_ID", all = TRUE)
head(comparison_Jcall_imgt)
#create a hit column
comparison_Jcall_imgt$HIT <- NA
#rename the column names to be more informative
colnames(comparison_Jcall_imgt) <- c("Antibodty_ID", "IMGT_Jcall", "True_Jcall", "HIT")
head(comparison_Jcall_imgt)
#where there is a match give TRUE , whereas a mishit give false on the HIT column
comparison_Jcall_imgt$HIT <- as.character(comparison_Jcall_imgt$IMGT_Jcall) == as.character(comparison_Jcall_imgt$True_Jcall)

head(comparison_Jcall_imgt)
#create frequency table of hits and mishits for Vgene
imgt_Jcall_hit_table <- table(comparison_Jcall_imgt$HIT)
imgt_Jcall_hit_table

#Notice there are rows where IMGT doesn't call any allele  e.g.
comparison_Jcall_imgt[6, ]

#lets remove the blank rows because they don't tell us anything concerning hit or mishit 
comparison_Jcall_imgt_wo_space <- comparison_Jcall_imgt[-which(comparison_Jcall_imgt$IMGT_Jcall == ""), ]
#confirm
comparison_Jcall_imgt_wo_space[6, ]

#create a frequency table of hits witout space 
imgt_Jcall_hit_table_wo_space <- table(comparison_Jcall_imgt_wo_space$HIT)
imgt_Jcall_hit_table_wo_space





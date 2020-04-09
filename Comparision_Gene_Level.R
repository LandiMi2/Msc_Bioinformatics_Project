setwd("/home/cofia/Documents/Msc_Bioinformatics_Project/Data/Annotatio_Output")

#mixcr <- read.delim("./MIXCR/Simulated_data/simulated_clones_0.5M.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
igblast <- read.delim("./Igblast/simulated_igblast_merged.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

#######Read VDJ recombination; we are using this file for comaprision because it contains information about V(D)J recombination for each reads in merged file.
read_vdj <- read.table("read_vdj.txt", header = F, sep = ";")
#rename column names 
colnames(read_vdj) <- c("Antibody_ID", "V_call", "D_call", "J_call")
View(read_vdj)

nrow(read_vdj)#420668
nrow(igblast)#420668

########Take 4 columns from annotation - (seq ID , Vcall, Dcall and J call)
igblast_specific <- data.frame(igblast$sequence_id, igblast$v_call, igblast$d_call, igblast$j_call) 


############################################### V call analysis for Igblast ############################################################################## 
igblast_Vcall <- data.frame(igblast_specific$igblast.sequence_id, igblast_specific$igblast.v_call)
head(igblast_Vcall)
#rename column name for Igblast Vcall data frame 
colnames(igblast_Vcall) <- c("Antibody_ID", "IgBlast_Vcall")
read_vdj_Vcall <- data.frame(read_vdj$Antibody_ID, read_vdj$V_call)
View(read_vdj_Vcall)
#rename column names for read_vdj_Vcall dataframe
colnames(read_vdj_Vcall) <- c("Antibody_ID", "True_Vcall")

#before comparing remove the allele naming to enable you compare by gene level. eg IGHV1-30*03 - IGHV1-30
#remove *03 ...the last three character that define the allele 
read_vdj_Vcall$True_Vcall <- gsub(".{3}$","",read_vdj_Vcall$True_Vcall)


#IgBLAST will give you miltiple genes for the same sequence, when it is undecided.
#e.g. this one
igblast_Vcall[4,]

#it is convention in these cases, to consider only the first one before matching with the true genes
igblast_Vcall$IgBlast_Vcall <- gsub("^(.*?),.*", "\\1", igblast_Vcall$IgBlast_Vcall)
#confirm 
igblast_Vcall[4,]
#remove *03 ...the last three character that define the allele 
igblast_Vcall$IgBlast_Vcall <- gsub(".{3}$","", igblast_Vcall$IgBlast_Vcall)
View(igblast_Vcall)
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
#remove *03 ...the last three character that define the allele 
read_vdj_Dcall$True_Dcall <- gsub(".{3}$","",read_vdj_Dcall$True_Dcall)

#IgBLAST will give you miltiple genes for the same sequence, when it is undecided.
#e.g. this one
igblast_Dcall[3,]

#it is convention in these cases, to consider only the first one before matching with the true genes

igblast_Dcall$IgBlast_Dcall <- gsub("^(.*?),.*", "\\1", igblast_Dcall$IgBlast_Dcall)
#confirm 
igblast_Dcall[3,]

#remove *03 ...the last three character that define the allele 
igblast_Dcall$IgBlast_Dcall <- gsub(".{3}$","", igblast_Dcall$IgBlast_Dcall)

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
View(comparison_Dcall_Igblast)
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

#remove *03 ...the last three character that define the allele 
read_vdj_Jcall$True_Jcall <- gsub(".{3}$","",read_vdj_Jcall$True_Jcall)
#IgBLAST will give you miltiple genes for the same sequence, when it is undecided.
#e.g. this one
igblast_Jcall[3,]

#it is convention in these cases, to consider only the first one before matching with the true genes

igblast_Jcall$IgBlast_Jcall <- gsub("^(.*?),.*", "\\1", igblast_Jcall$IgBlast_Jcall)
#confirm 
igblast_Jcall[3,]

#remove *03 ...the last three character that define the allele 
igblast_Jcall$IgBlast_Jcall <- gsub(".{3}$","", igblast_Jcall$IgBlast_Jcall)
View(igblast_Jcall)
#merger to see matches and mismatches
comparison_Jcall_Igblast <- merge(data.frame(igblast_Jcall), data.frame(read_vdj_Jcall), 
                                  by = "Antibody_ID" , all = TRUE)

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
#remove eg *03 ...the last three character that define the allele 
imgt_Vcall$IMGT_Vcall <- gsub(".{3}$","", imgt_Vcall$IMGT_Vcall)

read_vdj_imgt_Vcall <- data.frame(read_vdj_imgt$Sequence_ID, read_vdj_imgt$True_Vcall)
head(read_vdj_imgt_Vcall)
#remove *03 ...the last three character that define the allele 
read_vdj_imgt_Vcall$read_vdj_imgt.True_Vcall <- gsub(".{3}$","",read_vdj_imgt_Vcall$read_vdj_imgt.True_Vcall)
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
#remove ORF at the end of some gene calls 
imgt_Dcall$IMGT_Dcall <- gsub(" ORF", "", imgt_Dcall$IMGT_Dcall)
#remove *03 ...the last three character that define the allele 
imgt_Dcall$IMGT_Dcall <- gsub(".{3}$","", imgt_Dcall$IMGT_Dcall)

read_vdj_imgt_Dcall <- data.frame(read_vdj_imgt$Sequence_ID, read_vdj_imgt$True_Dcall)
#remove *03 ...the last three character that define the allele 
read_vdj_imgt_Dcall$read_vdj_imgt.True_Dcall <- gsub(".{3}$","",read_vdj_imgt_Dcall$read_vdj_imgt.True_Dcall)

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
#remove *03 ...the last three character that define the allele 
imgt_Jcall$IMGT_Jcall <- gsub(".{3}$","", imgt_Jcall$IMGT_Jcall)
head(imgt_Jcall)

read_vdj_imgt_Jcall <- data.frame(read_vdj_imgt$Sequence_ID, read_vdj_imgt$True_Jcall)
#remove *03 ...the last three character that define the allele 
read_vdj_imgt_Jcall$read_vdj_imgt.True_Jcall <- gsub(".{3}$","",read_vdj_imgt_Jcall$read_vdj_imgt.True_Jcall)
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


############################################## mixcr analysis ##############################
#read mixcr annotation output 
mixcr <-read.table("./MIXCR/Simulated_data/simulated_alignments.txt", header = TRUE, sep = "\t")
#rename column names 
colnames(mixcr) <- c("Sequence_ID", "Vcall", "Dcall", "Jcall")

#read read_vdj_recombination file 
read_vdj_mixcr <- read.table("./read_vdj.txt", header = F, sep = ";")
colnames(read_vdj_mixcr) <- c("Antibidy_ID", "V_call", "D_call", "J_call")

head(read_vdj_mixcr)
#modify the first column to match mixcr sequence ID 
read_vdj_mixcr$Antibidy_ID <- gsub("_merged_read_", "", read_vdj_mixcr$Antibidy_ID)
read_vdj_mixcr$Antibidy_ID <- gsub("^[0-9]+", "", read_vdj_mixcr$Antibidy_ID)

head(read_vdj_mixcr)
head(mixcr)


###################################### V call analysis ########################################
mixcr_Vcall <- data.frame(mixcr$Sequence_ID, mixcr$Vcall)
head(mixcr_Vcall)

read_vdj_mixcr_Vcall <- data.frame(read_vdj_mixcr$Antibidy_ID, read_vdj_mixcr$V_call)
head(read_vdj_mixcr_Vcall)
#remove *03 ...the last three character that define the allele 
read_vdj_mixcr_Vcall$read_vdj_mixcr.V_call <- gsub(".{3}$","",read_vdj_mixcr_Vcall$read_vdj_mixcr.V_call)

#modify V call column 
mixcr_Vcall$mixcr.Vcall <- gsub("Bostau_", "", mixcr_Vcall$mixcr.Vcall)
head(mixcr_Vcall)

#remove *03 ...the last three character that define the allele 
mixcr_Vcall$mixcr.Vcall <- gsub(".{3}$","", mixcr_Vcall$mixcr.Vcall)
#merge the two files
compare_mixcr_Vcall <- merge(data.frame(mixcr_Vcall), data.frame(read_vdj_mixcr_Vcall), by.x = "mixcr.Sequence_ID", 
                             by.y = "read_vdj_mixcr.Antibidy_ID", all = TRUE)


#create a HIT column 
compare_mixcr_Vcall$HIT <- NA

#where there is match print TRUE and vice versa
compare_mixcr_Vcall$HIT <- as.character(compare_mixcr_Vcall$mixcr.Vcall) == as.character(compare_mixcr_Vcall$read_vdj_mixcr.V_call)

#create frequency table of hits and mishits for Vgene
mixcr_Vcall_hit_table <- table(compare_mixcr_Vcall$HIT)
mixcr_Vcall_hit_table

######################################## D call analysis ########################################
mixcr_Dcall <- data.frame(mixcr$Sequence_ID, mixcr$Dcall)
head(mixcr_Dcall)

read_vdj_mixcr_Dcall <- data.frame(read_vdj_mixcr$Antibidy_ID, read_vdj_mixcr$D_call)
head(read_vdj_mixcr_Dcall)
#remove *03 ...the last three character that define the allele 
read_vdj_mixcr_Dcall$read_vdj_mixcr.D_call <- gsub(".{3}$","",read_vdj_mixcr_Dcall$read_vdj_mixcr.D_call)

#modify V call column 
mixcr_Dcall$mixcr.Dcall <- gsub("Bostau_", "", mixcr_Dcall$mixcr.Dcall)

#remove *03 ...the last three character that define the allele 
mixcr_Dcall$mixcr.Dcall <- gsub(".{3}$","", mixcr_Dcall$mixcr.Dcall)
head(mixcr_Dcall)
#merge the two files
compare_mixcr_Dcall <- merge(data.frame(mixcr_Dcall), data.frame(read_vdj_mixcr_Dcall), by.x = "mixcr.Sequence_ID", 
                             by.y = "read_vdj_mixcr.Antibidy_ID", all = TRUE)

head(compare_mixcr_Dcall)
#create a HIT column 
compare_mixcr_Dcall$HIT <- NA

#where there is match print TRUE and vice versa
compare_mixcr_Dcall$HIT <- as.character(compare_mixcr_Dcall$mixcr.Dcall) == as.character(compare_mixcr_Dcall$read_vdj_mixcr.D_call)

#create frequency table of hits and mishits for Vgene
mixcr_Dcall_hit_table <- table(compare_mixcr_Dcall$HIT)
mixcr_Dcall_hit_table

######################################## J call analysis ########################################
mixcr_Jcall <- data.frame(mixcr$Sequence_ID, mixcr$Jcall)


read_vdj_mixcr_Jcall <- data.frame(read_vdj_mixcr$Antibidy_ID, read_vdj_mixcr$J_call)
head(read_vdj_mixcr_Jcall)
#remove *03 ...the last three character that define the allele 
read_vdj_mixcr_Jcall$read_vdj_mixcr.J_call <- gsub(".{3}$","",read_vdj_mixcr_Jcall$read_vdj_mixcr.J_call)

#modify V call column 
mixcr_Jcall$mixcr.Jcall <- gsub("Bostau_", "", mixcr_Jcall$mixcr.Jcall)
#remove *03 ...the last three character that define the allele 
mixcr_Jcall$mixcr.Jcall <- gsub(".{3}$","", mixcr_Jcall$mixcr.Jcall)
head(mixcr_Jcall)
#merge the two files
compare_mixcr_Jcall <- merge(data.frame(mixcr_Jcall), data.frame(read_vdj_mixcr_Jcall), by.x = "mixcr.Sequence_ID", 
                             by.y = "read_vdj_mixcr.Antibidy_ID", all = TRUE)

head(compare_mixcr_Jcall)
#create a HIT column 
compare_mixcr_Jcall$HIT <- NA

#where there is match print TRUE and vice versa
compare_mixcr_Jcall$HIT <- as.character(compare_mixcr_Jcall$mixcr.Jcall) == as.character(compare_mixcr_Jcall$read_vdj_mixcr.J_call)
head(compare_mixcr_Jcall)
#create frequency table of hits and mishits for Vgene
mixcr_Jcall_hit_table <- table(compare_mixcr_Jcall$HIT)
mixcr_Jcall_hit_table

########################################### combine V, D, J calls mishit data for all tools ###################################################
#Igblast VDJ call 
igblast_Vcall_hit_table #False 15826  ,  True 404842 
igblast_Dcall_hit_table #False 178229 ,  True 242439 
igblast_Jcall_hit_table #False 170090 ,  True 250578  
#IMGT VDJ call 
imgt_Vcall_hit_table  #False  16346 , True 404322 
imgt_Dcall_hit_table #False 291328 , True 129340
imgt_Jcall_hit_table #False 181796 , True 238872 
#MiXCR VDJ call 
mixcr_Vcall_hit_table #False 33754 , True 235918 
mixcr_Dcall_hit_table #False 115539 , True 154133 
mixcr_Jcall_hit_table #False 40760 , True 228912

#calculate the percentage frequencies for each tool True and False counts 
#Igblast percentages
#percentage igblast V call
igblast_Vcall_hit_df <- as.data.frame(igblast_Vcall_hit_table)#convert table to dataframe 
igblast_True_V_count <- igblast_Vcall_hit_df$Freq[2]
igblast_False_V_count <- igblast_Vcall_hit_df$Freq[1]

percent_igblast_V_call_hit <- round((igblast_True_V_count / (igblast_True_V_count + igblast_False_V_count)) * 100)
percent_igblast_V_call_mishit <- round((igblast_False_V_count / (igblast_True_V_count + igblast_False_V_count)) * 100)

#percentage igblast D call
igblast_Dcall_hit_df <- as.data.frame(igblast_Dcall_hit_table)#convert table to dataframe 
igblast_True_D_count <- igblast_Dcall_hit_df$Freq[2]
igblast_False_D_count <- igblast_Dcall_hit_df$Freq[1]

percent_igblast_D_call_hit <- round((igblast_True_D_count / (igblast_True_D_count + igblast_False_D_count)) * 100)
percent_igblast_D_call_mishit <- round((igblast_False_D_count / (igblast_True_D_count + igblast_False_D_count)) * 100)

#percantage igblast J call 
igblast_Jcall_hit_df <- as.data.frame(igblast_Jcall_hit_table)#convert table to dataframe 
igblast_True_J_count <- igblast_Jcall_hit_df$Freq[2]
igblast_False_J_count <- igblast_Jcall_hit_df$Freq[1]

percent_igblast_J_call_hit <- round((igblast_True_J_count / (igblast_True_J_count + igblast_False_J_count)) * 100)
percent_igblast_J_call_mishit <- round((igblast_False_J_count / (igblast_True_J_count + igblast_False_J_count)) * 100)

#imgt percentages
#percentages imgt V call 
imgt_Vcall_hit_df <- as.data.frame(imgt_Vcall_hit_table)#convert table to dataframe 
imgt_True_V_count <- imgt_Vcall_hit_df$Freq[2]
imgt_False_V_count <- imgt_Vcall_hit_df$Freq[1]

percent_imgt_V_call_hit <- round((imgt_True_V_count / (imgt_True_V_count + imgt_False_V_count)) * 100)
percent_imgt_V_call_mishit <- round((imgt_False_V_count / (imgt_True_V_count + imgt_False_V_count)) * 100)

#percentages imgt D call 
imgt_Dcall_hit_df <- as.data.frame(imgt_Dcall_hit_table)#convert table to dataframe 
imgt_True_D_count <- imgt_Dcall_hit_df$Freq[2]
imgt_False_D_count <- imgt_Dcall_hit_df$Freq[1]

percent_imgt_D_call_hit <- round((imgt_True_D_count / (imgt_True_D_count + imgt_False_D_count)) * 100)
percent_imgt_D_call_mishit <- round((imgt_False_D_count / (imgt_True_D_count + imgt_False_D_count)) * 100)

#percentages imgt J call 
imgt_Jcall_hit_df <- as.data.frame(imgt_Jcall_hit_table)#convert table to dataframe 
imgt_True_J_count <- imgt_Jcall_hit_df$Freq[2]
imgt_False_J_count <- imgt_Jcall_hit_df$Freq[1]

percent_imgt_J_call_hit <- round((imgt_True_J_count / (imgt_True_J_count + imgt_False_J_count)) * 100)
percent_imgt_J_call_mishit <- round((imgt_False_J_count / (imgt_True_J_count + imgt_False_J_count)) * 100)

#mixcr percentages
#percentages mixcr V call 
mixcr_Vcall_hit_df <- as.data.frame(mixcr_Vcall_hit_table)#convert table to dataframe 
mixcr_True_V_count <- mixcr_Vcall_hit_df$Freq[2]
mixcr_False_V_count <- mixcr_Vcall_hit_df$Freq[1]

percent_mixcr_V_call_hit <- round((mixcr_True_V_count / (mixcr_True_V_count + mixcr_False_V_count)) * 100)
percent_mixcr_V_call_mishit <- round((mixcr_False_V_count / (mixcr_True_V_count + mixcr_False_V_count)) * 100)

#percentages mixcr D call 
mixcr_Dcall_hit_df <- as.data.frame(mixcr_Dcall_hit_table)#convert table to dataframe 
mixcr_True_D_count <- mixcr_Dcall_hit_df$Freq[2]
mixcr_False_D_count <- mixcr_Dcall_hit_df$Freq[1]

percent_mixcr_D_call_hit <- round((mixcr_True_D_count / (mixcr_True_D_count + mixcr_False_D_count)) * 100)
percent_mixcr_D_call_mishit <- round((mixcr_False_D_count / (mixcr_True_D_count + mixcr_False_D_count)) * 100)

#percentages mixcr J call 
mixcr_Jcall_hit_df <- as.data.frame(mixcr_Jcall_hit_table)#convert table to dataframe 
mixcr_True_J_count <- mixcr_Jcall_hit_df$Freq[2]
mixcr_False_J_count <- mixcr_Jcall_hit_df$Freq[1]

percent_mixcr_J_call_hit <- round((mixcr_True_J_count / (mixcr_True_J_count + mixcr_False_J_count)) * 100)
percent_mixcr_J_call_mishit <- round((mixcr_False_J_count / (mixcr_True_J_count + mixcr_False_J_count)) * 100)


#Make a plot for hits and mishits for the three tools

#create a martrix for the percentage hits 
percent_hit_matrix <- matrix(c(percent_igblast_V_call_hit, percent_mixcr_V_call_hit, percent_imgt_V_call_hit,
                               percent_igblast_D_call_hit, percent_mixcr_D_call_hit, percent_imgt_D_call_hit,
                               percent_igblast_J_call_hit, percent_mixcr_J_call_hit, percent_imgt_J_call_hit), ncol = 3, byrow = TRUE)

percent_hit_matrix
colnames(percent_hit_matrix) <- c("IgBlast", "MiXCR", "IMGT")
rownames(percent_hit_matrix) <- c("Vcall", "Dcall", "Jcall")

#convert the martrix to a table so that we can start to analysis these counts 
percent_hit_table <- as.table(percent_hit_matrix)
percent_hit_table

hit_plot <- barplot(percent_hit_table, legend.text = T, beside = TRUE, ylim=c(0,100), col = c("red", "green", "blue"), 
                    main = "Predicted percentage hit of annotation tools")

#add text ontop of the bar plots 
text(x = hit_plot, y = percent_hit_table, labels = percent_hit_table,  pos = 3,)

#create a martrix for the percentage mishits 
percent_mishit_matrix <- matrix(c(percent_igblast_V_call_mishit, percent_mixcr_V_call_mishit, percent_imgt_V_call_mishit,
                                  percent_igblast_D_call_mishit, percent_mixcr_D_call_mishit, percent_imgt_D_call_mishit,
                                  percent_igblast_J_call_mishit, percent_mixcr_J_call_mishit, percent_imgt_J_call_mishit), ncol = 3, byrow = TRUE)

percent_mishit_matrix
colnames(percent_mishit_matrix) <- c("IgBlast", "MiXCR", "IMGT")
rownames(percent_mishit_matrix) <- c("Vcall", "Dcall", "Jcall")

#convert the martrix to a table so that we can start to analysis these counts 
percent_mishit_table <- as.table(percent_mishit_matrix)
percent_mishit_table

mishit_plot <- barplot(percent_mishit_table, legend.text =T,args.legend = list(x = 'top'), beside = TRUE, ylim=c(0,100), col = c("red", "green", "blue"), 
                       main = "Predicted percentage error rate of annotation tools")

#add text ontop of the bar plots 
text(x = mishit_plot, y = percent_mishit_table, labels = percent_mishit_table,  pos = 3,)


#unassigned calls in percentage for D call and J call for IgBlast

percent_unassigned_Dcall_igblast <- round(((nrow(comparison_Dcall_Igblast) - nrow(comparison_Dcall_Igblast_wo_space))) / 
                                            nrow(comparison_Dcall_Igblast) * 100)
# 2
percent_unassigned_Jcall_igblast <- round(((nrow(comparison_Jcall_Igblast) - nrow(comparison_Jcall_Igblast_wo_space))) / 
                                            nrow(comparison_Jcall_Igblast) * 100)
#27

#unassigned calls in percentage for D and J calls for IMGT
percent_unassigned_Dcall_imgt <- round(((nrow(comparison_Dcall_imgt) - nrow(comparison_Dcall_imgt_wo_space))) / 
                                         nrow(comparison_Dcall_imgt) * 100)
#45
percent_unassigned_Jcall_imgt <- round(((nrow(comparison_Jcall_imgt) - nrow(comparison_Jcall_imgt_wo_space))) / 
                                         nrow(comparison_Jcall_imgt) * 100)
#26

percent_unassigned_matrix <- matrix(c(percent_unassigned_Dcall_igblast, percent_unassigned_Dcall_imgt,
                                      percent_unassigned_Jcall_igblast, percent_unassigned_Jcall_imgt), ncol = 2, byrow = TRUE)

colnames(percent_unassigned_matrix) <- c("IgBlast", "IMGT")
rownames(percent_unassigned_matrix) <- c("Dcall", "Jcall")
#convert matrix to table 
percent_unassigned_table <- as.table(percent_unassigned_matrix)

unassigned_plot <-  barplot(percent_unassigned_table, legend.text =T, beside = TRUE, ylim=c(0,100), col = c("green", "blue"), 
                            main = "Percentage unassigned D and J call") 
#add text ontop of the bar plots 
text(x = unassigned_plot, y = percent_unassigned_table, labels = percent_unassigned_table,  pos = 3,)

############################# Creating confusion matrixs ############################################
data_Igblast_Vgene <- data.frame(comparison_Vcall_Igblast$IgBlast_Vcall,comparison_Vcall_Igblast$True_Vcall,
                                 comparison_Vcall_Igblast$HIT)
#columnames
colnames(data_Igblast_Vgene) <- c("IgBlast_Vcall", "True_Vcall", "HIT")
View(data_Igblast_Vgene)
#add IGHV1-33 as a level in IgBlast Vcall
#below is not my function....
addLevel <- function(x, newlevel=NULL) {
  if(is.factor(x)) {
    if (is.na(match(newlevel, levels(x))))
      return(factor(x, levels=c(levels(x), newlevel)))
  }
  return(x)
}
data_Igblast_Vgene$IgBlast_Vcall <- addLevel(data_Igblast_Vgene$IgBlast_Vcall, "IGHV1-33")
str(data_Igblast_Vgene$IgBlast_Vcall)


data_mixcr_Vgene <- data.frame(compare_mixcr_Vcall$mixcr.Vcall, compare_mixcr_Vcall$read_vdj_mixcr.V_call,
                               compare_mixcr_Vcall$HIT)
colnames(data_mixcr_Vgene) <- c("MiXCR_Vcall","True_Vcall", "HIT")

data_mixcr_Vgene$MiXCR_Vcall <- addLevel(data_mixcr_Vgene$MiXCR_Vcall, "IGHV1S1")
data_mixcr_Vgene$MiXCR_Vcall <- addLevel(data_mixcr_Vgene$MiXCR_Vcall, "IGHV1S1")


data_imgt_Vgene <- data.frame(comparison_Vcall_imgt$IMGT_Vcall, comparison_Vcall_imgt$True_Vcall,
                              comparison_Vcall_imgt$HIT)
colnames(data_imgt_Vgene) <- c("IMGT_Vcall","True_Vcall", "HIT")
#to remove 466 empty rows from this datafram IMGT
data_imgt_Vgene <- data_imgt_Vgene[-which(data_imgt_Vgene$IMGT_Vcall == ""), ]
sum(data_imgt_Vgene$IMGT_Vcall == "") #confirm 
#### draw a hit map for V gene calls #############################
#construct a table for comaparison for Igblast V gene 
Igblast_Vgene_CrossTable <- data.frame(table(data_Igblast_Vgene$IgBlast_Vcall, data_Igblast_Vgene$True_Vcall))
#change column names 
colnames(Igblast_Vgene_CrossTable) <- c("Igblast_Vcall", "True_Vcall", "Counts")
View(Igblast_Vgene_CrossTable)

#log transform
Igblast_Vgene_CrossTable$Counts <- log(Igblast_Vgene_CrossTable$Counts)
library(ggplot2)
library(gridExtra)

heatmap_Vgene_Igblast <- ggplot(Igblast_Vgene_CrossTable, aes(Igblast_Vcall, True_Vcall)) +
  geom_tile(aes(fill = Counts), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red", name = "log(Counts)") + 
  ggtitle("Frequency of hits and mishits of IgBlast V gene annotation") +
  xlab("IgBlast V calls") + ylab("True V calls")

#construct a table for comaparison for IMGT V gene 
imgt_Vgene_CrossTable <- data.frame(table(data_imgt_Vgene$IMGT_Vcall, data_imgt_Vgene$True_Vcall))
#change column names 
colnames(imgt_Vgene_CrossTable) <- c("IMGT_Vcall", "True_Vcall", "Counts")
#log transformation
imgt_Vgene_CrossTable$Counts <- log(imgt_Vgene_CrossTable$Counts)
#remove empty rows 
imgt_Vgene_CrossTable <- imgt_Vgene_CrossTable[-which(imgt_Vgene_CrossTable$IMGT_Vcall == ""), ]

heatmap_Vgene_imgt <- ggplot(imgt_Vgene_CrossTable, aes(IMGT_Vcall, True_Vcall)) +
  geom_tile(aes(fill = Counts), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red", name = "log(Counts)") +
  ggtitle("Frequency of hits and mishits of IMGT V gene annotation") +
  xlab("IMGT V calls") + ylab("True V calls")

#construct a table for comaparison for MiXCR V gene 
mixcr_Vgene_CrossTable <- data.frame(table(data_mixcr_Vgene$MiXCR_Vcall, data_mixcr_Vgene$True_Vcall))

#change column names 
colnames(mixcr_Vgene_CrossTable) <- c("MiXCR_Vcall", "True_Vcall", "Counts")

#log transformation
mixcr_Vgene_CrossTable$Counts <- log(mixcr_Vgene_CrossTable$Counts)
#View(mixcr_Vgene_CrossTable)

heatmap_Vgene_mixcr <- ggplot(mixcr_Vgene_CrossTable, aes(MiXCR_Vcall, True_Vcall)) +
  geom_tile(aes(fill = Counts), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red", name = "log(Counts)") +
  ggtitle("Frequency of hits and mishits of MiXCR V gene annotation") +
  ylab("True V calls") + xlab("MiXCR V calls")

##### combine the three heat maps ################
heat_list <- list(heatmap_Vgene_Igblast, heatmap_Vgene_imgt, heatmap_Vgene_mixcr)
heat_list[["ncol"]] <- 1
do.call(grid.arrange, heat_list)


###################################### J gene analysis ########################################

data_Igblast_Jgene <- data.frame(comparison_Jcall_Igblast$IgBlast_Jcall,comparison_Jcall_Igblast$True_Jcall)
#columnames
colnames(data_Igblast_Jgene) <- c("IgBlast_Jcall", "True_Jcall")
#remove blank spaces 
data_Igblast_Jgene <- data_Igblast_Jgene [-which(data_Igblast_Jgene$IgBlast_Jcall == ""), ]

str(data_Igblast_Jgene)
# True column has only 2 levels ....add more 9 levels to fit a square matrix
#create a list of missing levels 
addlist_Jgene <- list("IGHJ1-1", "IGHJ1-2", "IGHJ1-3", "IGHJ1-4", "IGHJ1-5", "IGHJ2-3",
                      "IGHJ2-5", "IGHJ2-6")
#do a for loop to add the genes not present
for (a in addlist_Jgene){
  data_Igblast_Jgene$True_Jcall <- addLevel(data_Igblast_Jgene$True_Jcall, a)
}

#construct a table for comaparison for Igblast J gene 
Igblast_Jgene_CrossTable <- data.frame(table(data_Igblast_Jgene$IgBlast_Jcall, data_Igblast_Jgene$True_Jcall))
#change column names 
colnames(Igblast_Jgene_CrossTable) <- c("Igblast_Jcall", "True_Jcall", "Counts")
#remove the two emtpy spaces
Igblast_Jgene_CrossTable <- Igblast_Jgene_CrossTable[-which(Igblast_Jgene_CrossTable$Igblast_Jcall== ""),]

#log transform
Igblast_Jgene_CrossTable$Counts <- log(Igblast_Jgene_CrossTable$Counts)
head(Igblast_Jgene_CrossTable)

heatmap_Jgene_Igblast <- ggplot(Igblast_Jgene_CrossTable, aes(Igblast_Jcall, True_Jcall)) +
  geom_tile(aes(fill = Counts), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red", name = "log(Counts)") + 
  ggtitle("Frequency of hits and mishits of IgBlast J gene annotation") +
  ylab("True J calls") + xlab("IgBlast J calls")


#################### mixcr J gene analysis ############################
data_mixcr_Jgene <- data.frame(compare_mixcr_Jcall$mixcr.Jcall, compare_mixcr_Jcall$read_vdj_mixcr.J_call)
colnames(data_mixcr_Jgene) <- c("MiXCR_Jcall","True_Jcall")
str(data_mixcr_Jgene)
addlistmixc_Jgene <- list("IGHJ2-6", "IGHJ1-5")
#add the two genes 
for (a in addlistmixc_Jgene){
  data_mixcr_Jgene$True_Jcall <- addLevel(data_mixcr_Jgene$True_Jcall, a)
}

#construct a table for comaparison for MiXCCR J gene 
mixcr_Jgene_CrossTable <- data.frame(table(data_mixcr_Jgene$MiXCR_Jcall, data_mixcr_Jgene$True_Jcall))
#change column names 
colnames(mixcr_Jgene_CrossTable) <- c("MiXCR_Jcall", "True_Jcall", "Counts")
#log transformation
mixcr_Jgene_CrossTable$Counts <- log(mixcr_Jgene_CrossTable$Counts)

heatmap_Jgene_mixcr <- ggplot(mixcr_Jgene_CrossTable, aes(True_Jcall,MiXCR_Jcall)) +
  geom_tile(aes(fill = Counts), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red", name = "log(Counts)") + 
  ggtitle("Frequency of hits and mishits of MiXCR J gene annotation") +
  xlab("True J calls") + ylab("MiXCR J calls")


########################## imgt J gene coss tables ##################
data_imgt_Jgene <- data.frame(comparison_Jcall_imgt$IMGT_Jcall, comparison_Jcall_imgt$True_Jcall)
colnames(data_imgt_Jgene) <- c("IMGT_Jcall","True_Jcall")
#remove empty spaces 
data_imgt_Jgene <- data_imgt_Jgene[-which(data_imgt_Jgene$IMGT_Jcall == ""), ]
data_imgt_Jgene <- data_imgt_Jgene[-which(data_imgt_Jgene$IMGT_Jcall == "less than 6 nucleotides are alig"), ]

unique(data_imgt_Jgene$IMGT_Jcall)

#create a list of missing levels 
addlist_Jgene_imgt <- list("IGHJ1-1", "IGHJ1-2", "IGHJ1-3", "IGHJ1-4", "IGHJ1-5", "IGHJ2-3",
                      "IGHJ2-5", "IGHJ2-6")
#do a for loop to add the genes not present
for (a in addlist_Jgene_imgt){
  data_imgt_Jgene$True_Jcall <- addLevel(data_imgt_Jgene$True_Jcall, a)
}

#construct a table for comaparison for imgt J gene 
imgt_Jgene_CrossTable <- data.frame(table(data_imgt_Jgene$IMGT_Jcall, data_imgt_Jgene$True_Jcall))
#change column names 
colnames(imgt_Jgene_CrossTable) <- c("IMGT_Jcall", "True_Jcall", "Counts")
#log transformation
imgt_Jgene_CrossTable$Counts <- log(imgt_Jgene_CrossTable$Counts)
#remove emtpty space 
imgt_Jgene_CrossTable <- imgt_Jgene_CrossTable[-which(imgt_Jgene_CrossTable$IMGT_Jcall == ""), ]
imgt_Jgene_CrossTable <- imgt_Jgene_CrossTable[-which(imgt_Jgene_CrossTable$IMGT_Jcall == "less than 6 nucleotides are alig"), ]

heatmap_Jgene_imgt <- ggplot(imgt_Jgene_CrossTable, aes(True_Jcall,IMGT_Jcall)) +
  geom_tile(aes(fill = Counts), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red", name = "log(Counts)") + 
  ggtitle("Frequency of hits and mishits of IMGT J gene annotation") +
  xlab("True J calls") + ylab("MiXCR J calls")


##### combine the three heat maps ################
heat_list_Jgene <- list(heatmap_Jgene_Igblast, heatmap_Jgene_imgt, heatmap_Jgene_mixcr)
heat_list[["ncol"]] <- 1
do.call(grid.arrange, heat_list_Jgene)



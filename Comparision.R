setwd("/home/cofia/Documents/Msc_Bioinformatics_Project/Data/Annotatio_Output")

mixcr <- read.delim("./MIXCR/Simulated_data/simulated_clones_0.5M.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
igblast <- read.delim("./Igblast/simulated_igblast.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
imgt <- read.delim("./IMGT/Simulated_data/Simulated_Bovine_Annotation/4_IMGT-gapped-AA-sequences.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")


#######Read VDJ recombination; we are using this file for comaprision because it contains information about V(D)J recombination for each constructed antibody.
#replace tab separator with ; using sed 's/\t/;/g' repertoire_vdj_recombination.txt > rep_vdj_recombination.txt
#read the file to have separate columns 
repertoire_vdj_recombination <- read.table("rep_vdj_recombination.txt", header = FALSE, sep = ";")
#add column names to the dataframe 
colnames(repertoire_vdj_recombination) <- c("Sequence ID", "V call", "D call", "J calls")
#modify your antibody name by removing _ character 
repertoire_vdj_recombination$`Sequence ID`=gsub("\\_", "", repertoire_vdj_recombination$`Sequence ID`)
View(repertoire_vdj_recombination)
nrow(repertoire_vdj_recombination) # 163688 

########Take 4 columns from annotation - (seq ID , Vcall, Dcall and J call)
igblast_specific <- data.frame(igblast$sequence_id, igblast$v_call, igblast$d_call, igblast$j_call) 
head(igblast_specific)
nrow(igblast_specific)  # 422390

############################################### V call analysis for Igblast ############################################################################## 
igblast_Vcall <- data.frame(igblast_specific$igblast.sequence_id, igblast_specific$igblast.v_call)
head(igblast_Vcall)
rep_vdj_Vcall <- data.frame(repertoire_vdj_recombination$`Sequence ID`, repertoire_vdj_recombination$`V call`)
head(rep_vdj_Vcall)
# observe IgBLAST will give you miltiple genes for the same sequence, when it is undecided.
# e.g. this one
View(igblast_Vcall)
#IgBLAST will give you miltiple genes for the same sequence, when it is undecided.
#e.g. this one
igblast_Vcall[3,]

#it is convention in these cases, to consider only the first one before matching with the true genes

igblast_Vcall$igblast_specific.igblast.v_call <- gsub("^(.*?),.*", "\\1", igblast_Vcall$igblast_specific.igblast.v_call)
#confirm 
igblast_Vcall[3,]

#merger to see matches and mismatches
comparison_Vcall <- merge(data.frame(igblast_Vcall), data.frame(rep_vdj_Vcall), by.x = "igblast_specific.igblast.sequence_id",
                          by.y = "repertoire_vdj_recombination..Sequence.ID.", all = TRUE)

View(comparison_Vcall)
#create a hit column
comparison_Vcall$HIT <- NA
#rename the column names to be more informative
colnames(comparison_Vcall) <- c("Antibodty_ID", "Igblast_Vcall", "True_V_call", "Hit_or_Mishit")

#where there is a match give TRUE , whereas a mishit give false on the HIT column
comparison_Vcall$Hit_or_Mishit <- as.character(comparison_Vcall$Igblast_Vcall) == as.character(comparison_Vcall$True_V_call)
head(comparison_Vcall)

#create frequency table of hits and mishits for Vgene
igblast_Vcall_hit_mishit_table <- table(comparison_Vcall$Hit_or_Mishit)
igblast_Vcall_hit_mishit_table


#############################################   D  call analysis for Igblast ###################################################################

igblast_Dcall <- data.frame(igblast_specific$igblast.sequence_id, igblast_specific$igblast.d_call)
rep_vdj_Dcall <- data.frame(repertoire_vdj_recombination$`Sequence ID`, repertoire_vdj_recombination$`D call`)
head(igblast_Dcall)
head(rep_vdj_Dcall)

igblast_Dcall$igblast_specific.igblast.d_call <- gsub("^(.*?),.*", "\\1", igblast_Dcall$igblast_specific.igblast.d_call)


comparison_Dcall <- merge(data.frame(igblast_Dcall), data.frame(rep_vdj_Dcall), by.x = "igblast_specific.igblast.sequence_id",
                          by.y = "repertoire_vdj_recombination..Sequence.ID.", all = TRUE)
head(comparison_Dcall)
#create a hit column
comparison_Dcall$HIT <- NA

#rename the column names to be more informative
colnames(comparison_Dcall) <- c("Antibodty_ID", "Igblast_Dcall", "True_D_call", "Hit_or_Mishit")

#where there is a match give TRUE , whereas a mishit give false on the HIT column
comparison_Dcall$Hit_or_Mishit <- as.character(comparison_Dcall$Igblast_Dcall) == as.character(comparison_Dcall$True_D_call)

#create frequency table of hits and mishits for Vgene
igblast_Dcall_hit_mishit_table <- table(comparison_Dcall$Hit_or_Mishit)
igblast_Dcall_hit_mishit_table

#############################################   J  call analysis for Igblast ###################################################################

igblast_Jcall <- data.frame(igblast_specific$igblast.sequence_id, igblast_specific$igblast.j_call)
rep_vdj_Jcall <- data.frame(repertoire_vdj_recombination$`Sequence ID`, repertoire_vdj_recombination$`J calls`)

head(rep_vdj_Jcall)

igblast_Jcall$igblast_specific.igblast.j_call <- gsub("^(.*?),.*", "\\1", igblast_Jcall$igblast_specific.igblast.j_call)

#merger to see matches and mismatches
comparison_Jcall <- merge(data.frame(igblast_Jcall), data.frame(rep_vdj_Jcall), by.x = "igblast_specific.igblast.sequence_id",
                          by.y = "repertoire_vdj_recombination..Sequence.ID.", all = TRUE)
#create a hit column
comparison_Jcall$HIT <- NA

#rename the column names to be more informative
colnames(comparison_Jcall) <- c("Antibodty_ID", "Igblast_Jcall", "True_J_call", "Hit_or_Mishit")


#where there is a match give TRUE , whereas a mishit give false on the HIT column
comparison_Jcall$Hit_or_Mishit <- as.character(comparison_Jcall$Igblast_Jcall) == as.character(comparison_Jcall$True_J_call)

View(comparison_Jcall)
#Note that there are alot of blank row examples
comparison_Jcall[8,]
#remove blank to see what frequency I get for hit and mishit 

comparison_Jcall_wo_space <- comparison_Jcall[-which(comparison_Jcall$Igblast_Jcall == ""), ]

#create a hit and mishit table for Jgene
igblast_Jcall_hit_mishit_table <- table(comparison_Jcall$Hit_or_Mishit)
igblast_Jcall_hit_mishit_table

#create a modified table for J genes
igblast_Jcall_hit_mishit_table_modified <- table(comparison_Jcall_wo_space$Hit_or_Mishit)
igblast_Jcall_hit_mishit_table_modified
########################################### combine V, D, J calls mishit data ##################################################################
igblast_Vcall_hit_mishit_table #False 157864 , True 5723 
igblast_Dcall_hit_mishit_table #False 156863 , True 6724
igblast_Jcall_hit_mishit_table #False 156981 , True 6606 

nrow(comparison_Vcall)

#create a martrix for the frequencies 
igblast_hit_mishit <- matrix(c(157864, 5723, 156863, 6724, 156981, 6606), ncol = 2, byrow = TRUE)
colnames(igblast_hit_mishit) <- c("Mishits", "Hits")
rownames(igblast_hit_mishit) <- c("V call", "D call", "J call")
igblast_hit_mishit
str(igblast_hit_mishit)

#convert the martrix to a table so that we can start to analysis these counts 
igblast_hits_table <- as.table(igblast_hit_mishit)

igblast_hits_table

igblast_hits_df <- as.data.frame(igblast_hit_mishit)
igblast_hits_df

barplot(igblast_hits_table,legend.text = T,beside = TRUE, col = c("red", "green", "yellow"), 
        main = "Igblast annotation analysis of hits and mishits")

#plot only hits 
library(ggplot2)

ggplot(igblast_hits_df, aes(x= rownames(igblast_hits_df) ,y=igblast_hits_df$Hits)) + geom_bar(stat = "identity")


########################################## Do IMGT analysis ######################################################
########Take 4 columns from annotation of imgt - (seq ID , Vcall, Dcall and J call)
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
imgt_Vcall$imgt_specific.imgt.V.GENE.and.allele <-  gsub('.{1}$', "", imgt_Vcall$imgt_specific.imgt.V.GENE.and.allele)
head(imgt_Vcall)
head(rep_vdj_Vcall)

#merger to see matches and mismatches

comparison_Vcall_imgt <- merge(data.frame(imgt_Vcall), data.frame(rep_vdj_Vcall), by.x = "imgt_specific.imgt.Sequence.ID",
                          by.y = "repertoire_vdj_recombination..Sequence.ID.", all = TRUE)


#create a hit column
comparison_Vcall_imgt$HIT <- NA

#rename the column names to be more informative
colnames(comparison_Vcall_imgt) <- c("Antibodty_ID", "IMGT_Vcall", "True_Vcall", "Hit_or_Mishit")

View(comparison_Vcall_imgt)
#where there is a match give TRUE , whereas a mishit give false on the HIT column


comparison_Vcall_imgt$Hit_or_Mishit <- as.character(comparison_Vcall_imgt$IMGT_Vcall) == as.character(comparison_Vcall_imgt$True_Vcall)

View(comparison_Vcall_imgt)







#This gives you your expected mishit rate 
#using the mishit rate you calculate of the expected number of mistakes for all the tools , then compare with tool is better for you 



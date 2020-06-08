setwd("/home/cofia/Documents/Msc_Bioinformatics_Project/TIgGER/Results/")
library(tidyr)
############################ Boran ###################################
Boran <- read.table("./Bovine_Novel_alleles/NovelBoran.txt", sep = "\t")
View(Boran)
#filter out novel allele that have 0 counts of unique CDR3 and J's
Boran_filtered <- dplyr::filter(Boran, grepl("^[1-9]" , Boran$novel_imgt_unique_cdr3))
#create a dataframe of novel allele name and novel sequences discovered
Boran_seq_df <- data.frame(Boran_filtered$germline_call, Boran_filtered$novel_imgt)
View(Boran_seq_df)
#change tab file to fasta file by first writing into tab separated file then use awk
write.table(Boran_seq_df, "Boran_sequences.tab", sep = "\t")
#on the terminal use awk '{print ">"$2"\n"$3}'

######################## Friesian ###############################
Friesian <- read.table("./Bovine_Novel_alleles/NovelFriesian.txt", sep = "\t")
#create a datafame of sequences and id
Friesian_seq_df <- data.frame(Friesian$germline_call, Friesian$novel_imgt)
#change tab file to fasta file by first writing into tab separated file then use awk
write.table(Friesian_seq_df, "Friesian_sequences.tab", sep = "\t")

####################### Ndama ######################################
Ndama <- read.table("./Bovine_Novel_alleles/NovelNdama.txt", sep = "\t")
#View(Ndama)
#create a datafame of sequences and id
Ndama_seq_df <- data.frame(Ndama$germline_call, Ndama$novel_imgt)
#change tab file to fasta file by first writing into tab separated file then use awk
write.table(Ndama_seq_df, "Ndama_sequences.tab", sep = "\t")

####################### Ankole ######################################
Ankole <- read.table("./Bovine_Novel_alleles/NovelAnkole.txt", sep = "\t")
#View(Ankole)
#create a datafame of sequences and id
Ankole_seq_df <- data.frame(Ankole$germline_call, Ankole$novel_imgt)
#change tab file to fasta file by first writing into tab separated file then use awk
write.table(Ankole_seq_df, "Ankole_sequences.tab", sep = "\t")

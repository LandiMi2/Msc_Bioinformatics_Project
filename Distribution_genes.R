setwd("/home/cofia/Documents/Msc_Bioinformatics_Project/Data/Annotatio_Output")
#open igblast annotation of diverse antibodty
igblast <- read.delim("./Igblast/simulated_igblast_diverse.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#######################################################################################################
#extract the V call column and remove the last three characters that define the alleles 
igblast_Vgenes <- data.frame(igblast$v_call)
#consider the first allele called 
igblast_Vgenes$igblast.v_call <- gsub("^(.*?),.*", "\\1", igblast_Vgenes$igblast.v_call)
#make a gene column tables 
igblast_Vgenes$igblast.v_call <- gsub(".{3}$","", igblast_Vgenes$igblast.v_call)
#rename the column name 
colnames(igblast_Vgenes) <- c("Igblast_V_genes")
head(igblast_Vgenes)
#create a dataframe having unique values and their frequencies 
igblast_Vgenes_counts <- as.data.frame(table(igblast_Vgenes))
View(igblast_Vgenes_counts)
#convert the frequencies in percentages to make a frequencies barplot and a piechart 
igblast_Vgenes_counts$Percent <- round((igblast_Vgenes_counts$Freq / sum(igblast_Vgenes_counts$Freq)) * 100)
#create a pie chart
require(RColorBrewer) #color options for the pie chart
#modify labels
lbls_percent_V <- paste(igblast_Vgenes_counts$Percent,"%", sep="") 
lbls_V <- paste(igblast_Vgenes_counts$igblast_Vgenes,lbls_percent_V,sep="  ") 
pie(igblast_Vgenes_counts$Percent, labels = lbls_V,
    col = brewer.pal(length(igblast_Vgenes_counts$igblast_Vgenes), "Paired"), main = "Igblast Vgene distribution")

##############################################################################################
#extract the D call column and remove the last three characters that define the alleles 
igblast_Dgenes <- data.frame(igblast$d_call)
#consider the first allele called 
igblast_Dgenes$igblast.d_call <- gsub("^(.*?),.*", "\\1", igblast_Dgenes$igblast.d_call)
#make a gene column tables 
igblast_Dgenes$igblast.d_call <- gsub(".{3}$","", igblast_Dgenes$igblast.d_call)
#rename the column name 
colnames(igblast_Dgenes) <- c("Igblast_D_genes")
head(igblast_Dgenes)
#create a dataframe having unique values and their frequencies 
igblast_Dgenes_counts <- as.data.frame(table(igblast_Dgenes))
View(igblast_Dgenes_counts)
#convert the frequencies in percentages to make a frequencies barplot and a piechart 
igblast_Dgenes_counts$Percent <- round((igblast_Dgenes_counts$Freq / sum(igblast_Dgenes_counts$Freq)) * 100)
#create a pie chart
require(RColorBrewer) #color options for the pie chart
#modify labels
lbls_percent_D <- paste(igblast_Dgenes_counts$Percent,"%", sep="") 
lbls_D <- paste(igblast_Dgenes_counts$igblast_Dgenes,lbls_percent_D,sep="  ") 

pie(igblast_Dgenes_counts$Percent, labels = lbls_D,
    col = brewer.pal(length(igblast_Dgenes_counts$igblast_Dgenes), "Paired"), main = "Igblast Dgene distribution")


##############################################################################################
#extract the J call column and remove the last three characters that define the alleles 
igblast_Jgenes <- data.frame(igblast$j_call)
#consider the first allele called 
igblast_Jgenes$igblast.j_call <- gsub("^(.*?),.*", "\\1", igblast_Jgenes$igblast.j_call)
#make a gene column tables 
igblast_Jgenes$igblast.j_call <- gsub(".{3}$","", igblast_Jgenes$igblast.j_call)
#rename the column name 
colnames(igblast_Jgenes) <- c("Igblast_J_genes")
head(igblast_Jgenes)
#create a dataframe having unique values and their frequencies 
igblast_Jgenes_counts <- as.data.frame(table(igblast_Jgenes))
View(igblast_Jgenes_counts)
#convert the frequencies in percentages to make a frequencies barplot and a piechart 
igblast_Jgenes_counts$Percent <- round((igblast_Jgenes_counts$Freq / sum(igblast_Jgenes_counts$Freq)) * 100)
#create a pie chart
require(RColorBrewer) #color options for the pie chart
#modify labels
lbls_percent_J <- paste(igblast_Jgenes_counts$Percent,"%", sep="") 
lbls_J <- paste(igblast_Jgenes_counts$igblast_Jgenes,lbls_percent_J,sep="  ") 

pie(igblast_Jgenes_counts$Percent, labels = lbls_J,
    col = brewer.pal(length(igblast_Jgenes_counts$igblast_Jgenes), "Paired"), main = "Igblast Jgene distribution")

############################## mixcr distribution ###################################################
#open mixcr annotation of diverse antibodty
mixcr <- read.delim("./MIXCR/Simulated_data/diverse_simulated_alignments.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#######################################################################################################
#extract the V call column and remove the last three characters that define the alleles 
mixcr_Vgenes <- data.frame(mixcr$bestVHit)
head(mixcr_Vgenes)
mixcr_Vgenes$mixcr.bestVHit <- gsub("Bostau_", "", mixcr_Vgenes$mixcr.bestVHit)
#make a gene column tables 
mixcr_Vgenes$mixcr.bestVHit <- gsub(".{3}$","", mixcr_Vgenes$mixcr.bestVHit)
colnames(mixcr_Vgenes) <- c("mixcr_V_genes")

#create a dataframe having unique values and their frequencies 
mixcr_Vgenes_counts <- as.data.frame(table(mixcr_Vgenes))

#convert the frequencies in percentages to make a frequencies barplot and a piechart 
mixcr_Vgenes_counts$Percent <- round((mixcr_Vgenes_counts$Freq / sum(mixcr_Vgenes_counts$Freq)) * 100)
#create a pie chart
require(RColorBrewer) #color options for the pie chart
#modify labels
lbls_percent_mixcr_V <- paste(mixcr_Vgenes_counts$Percent,"%", sep="") 
lbls_mixcr_V <- paste(mixcr_Vgenes_counts$mixcr_Vgenes,lbls_percent_mixcr_V,sep="  ") 
pie(mixcr_Vgenes_counts$Percent, labels = lbls_mixcr_V,
    col = brewer.pal(length(mixcr_Vgenes_counts$mixcr_Vgenes), "Paired"), main = "MiXCR Vgene distribution")

##############################################################################################
#extract the D call column and remove the last three characters that define the alleles 
mixcr_Dgenes <- data.frame(mixcr$bestDHit)
head(mixcr_Dgenes)
mixcr_Dgenes$mixcr.bestDHit <- gsub("Bostau_", "", mixcr_Dgenes$mixcr.bestDHit)
#make a gene column tables 
mixcr_Dgenes$mixcr.bestDHit <- gsub(".{3}$","", mixcr_Dgenes$mixcr.bestDHit)
colnames(mixcr_Dgenes) <- c("mixcr_D_genes")

#create a dataframe having unique values and their frequencies 
mixcr_Dgenes_counts <- as.data.frame(table(mixcr_Dgenes))

#convert the frequencies in percentages to make a frequencies barplot and a piechart 
mixcr_Dgenes_counts$Percent <- round((mixcr_Dgenes_counts$Freq / sum(mixcr_Dgenes_counts$Freq)) * 100)
#create a pie chart
require(RColorBrewer) #color options for the pie chart
#modify labels
lbls_percent_mixcr_D<- paste(mixcr_Dgenes_counts$Percent,"%", sep="") 
lbls_mixcr_D <- paste(mixcr_Dgenes_counts$mixcr_Dgenes,lbls_percent_mixcr_D,sep="  ") 
pie(mixcr_Dgenes_counts$Percent, labels = lbls_mixcr_D,
    col = brewer.pal(length(mixcr_Dgenes_counts$mixcr_Dgenes), "Paired"), main = "MiXCR Dgene distribution")

#################################################################################################
#extract the J call column and remove the last three characters that define the alleles 
mixcr_Jgenes <- data.frame(mixcr$bestJHit)
head(mixcr_Jgenes)
mixcr_Jgenes$mixcr.bestJHit <- gsub("Bostau_", "", mixcr_Jgenes$mixcr.bestJHit)
#make a gene column tables 
mixcr_Jgenes$mixcr.bestJHit <- gsub(".{3}$","", mixcr_Jgenes$mixcr.bestJHit)
colnames(mixcr_Jgenes) <- c("mixcr_J_genes")

#create a dataframe having unique values and their frequencies 
mixcr_Jgenes_counts <- as.data.frame(table(mixcr_Jgenes))

#convert the frequencies in percentages to make a frequencies barplot and a piechart 
mixcr_Jgenes_counts$Percent <- round((mixcr_Jgenes_counts$Freq / sum(mixcr_Jgenes_counts$Freq)) * 100)
#create a pie chart
require(RColorBrewer) #color options for the pie chart
#modify labels
lbls_percent_mixcr_J <- paste(mixcr_Jgenes_counts$Percent,"%", sep="") 
lbls_mixcr_J <- paste(mixcr_Jgenes_counts$mixcr_Jgenes,lbls_percent_mixcr_J,sep="  ") 
pie(mixcr_Jgenes_counts$Percent, labels = lbls_mixcr_J,
    col = brewer.pal(length(mixcr_Jgenes_counts$mixcr_Jgenes), "Paired"), main = "MiXCR J gene distribution")


############################## IMGT distribution ###################################################
#open mixcr annotation of diverse antibodty
imgt <- read.delim("./IMGT/Bovine_Simulated_Annotation_diverse/4_IMGT-gapped-AA-sequences.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
#######################################################################################################
#extract the V call column and remove the last three characters that define the alleles 
imgt_Vgenes <- data.frame(imgt$V.GENE.and.allele)
head(imgt_Vgenes)
imgt_Vgenes$imgt.V.GENE.and.allele <-  gsub("Bostau ", "", imgt_Vgenes$imgt.V.GENE.and.allele)
#consider the first gene in the column 
imgt_Vgenes$imgt.V.GENE.and.allele <-  gsub("^(.*?),.*", "\\1", imgt_Vgenes$imgt.V.GENE.and.allele)
#remove the last character P and F
imgt_Vgenes$imgt.V.GENE.and.allele <-  gsub(' .{1}$', "", imgt_Vgenes$imgt.V.GENE.and.allele)
#remove eg *03 ...the last three character that define the allele 
imgt_Vgenes$imgt.V.GENE.and.allele <- gsub(".{3}$","", imgt_Vgenes$imgt.V.GENE.and.allele)

#create a dataframe having unique values and their frequencies 
imgt_Vgenes_counts <- as.data.frame(table(imgt_Vgenes))

#convert the frequencies in percentages to make a frequencies barplot and a piechart 
imgt_Vgenes_counts$Percent <- round((imgt_Vgenes_counts$Freq / sum(imgt_Vgenes_counts$Freq)) * 100)
#create a pie chart
require(RColorBrewer) #color options for the pie chart
#modify labels
lbls_percent_imgt_V <- paste(imgt_Vgenes_counts$Percent,"%", sep="") 
lbls_imgt_V <- paste(imgt_Vgenes_counts$imgt_Vgenes,lbls_percent_imgt_V,sep="  ") 
pie(imgt_Vgenes_counts$Percent, labels = lbls_imgt_V,
    col = brewer.pal(length(imgt_Vgenes_counts$imgt_Vgenes), "Paired"), main = "IMGT Vgene distribution")

###################################################################################################

#extract the D call column and remove the last three characters that define the alleles 
imgt_Dgenes <- data.frame(imgt$D.GENE.and.allele)
View(imgt_Dgenes)
imgt_Dgenes$imgt.D.GENE.and.allele <-  gsub("Bostau ", "", imgt_Dgenes$imgt.D.GENE.and.allele)
#consider the first gene in the column 
imgt_Dgenes$imgt.D.GENE.and.allele <-  gsub("^(.*?),.*", "\\1", imgt_Dgenes$imgt.D.GENE.and.allele)
#remove the last character P and F
imgt_Dgenes$imgt.D.GENE.and.allele <-  gsub(' .{1}$', "", imgt_Dgenes$imgt.D.GENE.and.allele)
#remove ORF at the end of some gene calls 
imgt_Dgenes$imgt.D.GENE.and.allele <- gsub(" ORF", "", imgt_Dgenes$imgt.D.GENE.and.allele )
#remove eg *03 ...the last three character that define the allele 
imgt_Dgenes$imgt.D.GENE.and.allele <- gsub(".{3}$","", imgt_Dgenes$imgt.D.GENE.and.allele)
View(imgt_Dgenes)
#create a dataframe having unique values and their frequencies 
imgt_Dgenes_counts <- as.data.frame(table(imgt_Dgenes))

#convert the frequencies in percentages to make a frequencies barplot and a piechart 
imgt_Dgenes_counts$Percent <- round((imgt_Dgenes_counts$Freq / sum(imgt_Dgenes_counts$Freq)) * 100)
#create a pie chart
require(RColorBrewer) #color options for the pie chart
#modify labels
lbls_percent_imgt_D <- paste(imgt_Dgenes_counts$Percent,"%", sep="") 
lbls_imgt_D <- paste(imgt_Dgenes_counts$imgt_Dgenes,lbls_percent_imgt_D,sep="  ") 
pie(imgt_Dgenes_counts$Percent, labels = lbls_imgt_D,
    col = brewer.pal(length(imgt_Dgenes_counts$imgt_Dgenes), "Paired"), main = "IMGT Dgene distribution")


####################################################################################################
#extract the V call column and remove the last three characters that define the alleles 
imgt_Jgenes <- data.frame(imgt$J.GENE.and.allele)
View(imgt_Jgenes)
imgt_Jgenes$imgt.J.GENE.and.allele <-  gsub("Bostau ", "", imgt_Jgenes$imgt.J.GENE.and.allele)
#consider the first gene in the column 
imgt_Jgenes$imgt.J.GENE.and.allele <-  gsub("^(.*?),.*", "\\1", imgt_Jgenes$imgt.J.GENE.and.allele)
#remove the last character P and F
imgt_Jgenes$imgt.J.GENE.and.allele <-  gsub(' .{1}$', "", imgt_Jgenes$imgt.J.GENE.and.allele)
#remove ORF at the end of some gene calls 
imgt_Jgenes$imgt.J.GENE.and.allele <- gsub(" ORF", "", imgt_Jgenes$imgt.J.GENE.and.allele )
#remove eg *03 ...the last three character that define the allele 
imgt_Jgenes$imgt.J.GENE.and.allele <- gsub(".{3}$","", imgt_Jgenes$imgt.J.GENE.and.allele)

#create a dataframe having unique values and their frequencies 
imgt_Jgenes_counts <- as.data.frame(table(imgt_Jgenes))

#convert the frequencies in percentages to make a frequencies barplot and a piechart 
imgt_Jgenes_counts$Percent <- round((imgt_Jgenes_counts$Freq / sum(imgt_Jgenes_counts$Freq)) * 100)
#create a pie chart
require(RColorBrewer) #color options for the pie chart
#modify labels
lbls_percent_imgt_J <- paste(imgt_Jgenes_counts$Percent,"%", sep="") 
lbls_imgt_J <- paste(imgt_Jgenes_counts$imgt_Jgenes,lbls_percent_imgt_J,sep="  ") 
pie(imgt_Jgenes_counts$Percent, labels = lbls_imgt_J,
    col = brewer.pal(length(imgt_Jgenes_counts$imgt_Jgenes), "Paired"), main = "IMGT Jgene distribution")


mixcr_Dgenes_counts
imgt_Dgenes_counts

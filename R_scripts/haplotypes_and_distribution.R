# haplotype network analysis with R and the ape and pegas packages
#install.packages("devtools", dependencies = TRUE)
#library(devtools)
#install.packages("pegas", dependencies = T)
#install_github("thibautjombart/adegenet")
#install_github("emmanuelparadis/pegas/pegas")
library(pegas) #ape v5.3 and pegas v0.13
library(DECIPHER)
library(magrittr) 
library(dplyr)
library(ape)
library(cultevo)
library(ggplot2)
library(readxl)
#####
setwd("~/Documents/Msc_Bioinformatics_Project/Combined_Igdiscover_Tigger/")
setwd("~/Documents/Msc_Bioinformatics_Project/IgDiscover/")
setwd("~/Documents/Msc_Bioinformatics_Project/TIgGER/Results/Bovine_Novel_fasta/")
#setwd("~/Documents/Msc_Bioinformatics_Project/IgDiscover/Bovine_expt_10_iteration_Ankole/database")
#setwd("~/Documents/Msc_Bioinformatics_Project/TIgGER/Results/Bovine_Novel_fasta/")

#### ALIGN SEQUENCES FIRST ####
seqs <- readDNAStringSet("./", format = "fasta")
# nucleotide sequences need to be in the same orientation. if they are not, then they can be reoriented (optional)
seqs <- OrientNucleotides(seqs)
# perform the alignment
aligned <- AlignSeqs(seqs)
# write the alignment to a new FASTA file
writeXStringSet(aligned,file="./Align_allele/Boran_tigger_align.fasta")

#### READING THE INPUT DATA (ALIGNMENT) ####
# We first read the alignment into a DNAbin object:
alignment <- read.dna("./Align_allele/friesian_align.fasta", format = "fasta") 
class(alignment) # it's an object of class "DNAbin"
str(alignment) # 45 sequences

#### CREATE A HAPLOTYPE NETWORK ######
h <- haplotype(alignment)


# create a table with sequence name and haplotype names 
## this table will have 27 rows:
haplotype_t <- tibble(sequence_name = character(), haplotype_name = character())
#concatenate sequence name and the type of haplotype 
name <- paste(rownames(alignment), rownames(h), sep = " - ")


haplotype_t %<>%bind_rows(list(sequence_name = rownames(alignment), haplotype_name = rownames(h)))
haplotype_table <- table(haplotype_t)
nrow(haplotype_table) # ok, 27

#check out just a chuck of the table 1-10 row and column 
haplotype_table[1:10, 1:10]  

#create a network
net <- haploNet(h)
str(name)

#plot the network
pdf("./Networks/Friesian_net1.pdf", width = 14, height = 10)
plot(net, size= 4, scale.ratio = 1,cex = 1.2, threshold = 0,
     show.mutation = 3, asp = 0.5,legend= F, main = "Friesian Breed VH alleles discovered by IgDiscover")
legend(x =-100, y = 0, ncol = 2, legend = name, bty = "n", pch = 20, title = "VH alleles")
dev.off()

############## igraph ####################
library(igraph)

plot(as.igraph(net, directed = F, use.labels = T, altlinks = T), 
     main = "Friesian VH alleles discovered by IgDiscover")



####### hamming distances of African and western breeds novel alelles discovered by IgDiscover 
alignment_african <- read.dna("./Align_allele/combine_african_align.fasta", format = "fasta")

hamming_distances_african <- hammingdists(alignment_african)


hamming_distances_african_matrix <- as.matrix(hamming_distances_african)
#extracting distances for pairs defined by the upper triangle of the distance matrix
xy <- t(combn(colnames(hamming_distances_african_matrix), 2))
hamming_distances_african_df <- data.frame(xy, dist=hamming_distances_african_matrix[xy])
View(hamming_distances_african_df)
colnames(hamming_distances_african_df) <- c("allele_X1", "allele_X2", "African_Dist")
#add a column of breed type
hamming_distances_african_df$Breed <- "African_IgDiscover"
#calculate p-value
summary(hamming_distances_african_df$African_Dist)
t.test(hamming_distances_african_df$African_Dist)

#export
write.csv(hamming_distances_african_df, "./Hamming_Dists/african1.csv")

######### friesian novel alleles
alignment_friesian <- read.dna("./Align_allele/friesian_align.fasta", format = "fasta")
hamm_dist_friesian <- hammingdists(alignment_friesian)
str(hamm_dist_friesian)


hamm_dist_friesian_matrix <- as.matrix(hamm_dist_friesian)
#extracting distances for pairs defined by the upper triangle of the distance matrix
ab <- t(combn(colnames(hamm_dist_friesian_matrix), 2))
hamm_dist_friesian_df <- data.frame(ab, dist=hamm_dist_friesian_matrix[ab])
View(hamm_dist_friesian_df)
colnames(hamm_dist_friesian_df) <- c("allele_X1", "allele_X2", "Friesian_Dist")
hamm_dist_friesian_df$Breed <- "Friesian_IgDiscover"

t.test(hamming_distances_african_df$African_Dist, hamm_dist_friesian_df$Friesian_Dist)
#export 
write.csv(hamm_dist_friesian_df, "./Hamming_Dists/friesian1.csv")

##plotting a box plot

hamm_dists <- read_excel("./Hamming_Dists/african_friesian_hamm_dists.xlsx")
View(hamm_dists)
head(hamm_dists)


ggplot(data = hamm_dists, aes(x = Breed, y = Hamm_Dist)) +  geom_jitter(aes(size = Hamm_Dist, colour =Breed)) + 
  geom_boxplot(size = 0.7, alpha=0.7, outlier.colour = "black", varwidth = T) +
  xlab("African breeds vs Friesian breed")+ ylab("Hamming distances") +
  ggtitle("Pairwise hamming distances of novel alelles discovered by TIgGER")
  
#extra calculations 
sd(hamming_distances_african_df$African_Dist) #10.03485-IgDiscover; 13.45994-TIgGER
sd(hamm_dist_friesian_df$Friesian_Dist) #8.437102-IgDiscover; 10.54795 - TIgGER
median(hamming_distances_african_df$African_Dist) #42-IgDiscover; 59-TIgGER
mean(hamming_distances_african_df$African_Dist)#41.93162-IgDiscover; 56.47692-TIgGER
mean(hamm_dist_friesian_df$Friesian_Dist)#31.05445-IgDiscover; 29-TIgGER
median(hamm_dist_friesian_df$Friesian_Dist)#31-IgDiscover; 27.5-TIgGER

###ploting histogram of the distribution of pairwise hamming distances

ggplot(hamm_dists, aes(Hamm_Dist, fill = Breed)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..),position = 'identity', bins = 20) +
  geom_density(linetype = "dashed", alpha = 0.1, adjust = 2) +
  xlab("Pairwise hamming distances") + ylab("Density") + 
  ggtitle("Distribution of pairwise distances of novel alleles discovered by TIgGER")


###### plot a boxplot combining all novel alleles discovered by IgDiscover and TIgGER 

hamm_dists_combined <- read_excel("../Combined_Igdiscover_Tigger/hamming_dist_plots/combined_african_friesian_tigger_igdiscover.xlsx")

head(hamm_dists_combined)

ggplot(data = hamm_dists_combined, aes(x = Breed, y = Hamm_Dist)) +  geom_jitter(aes(colour =Breed)) + 
  geom_boxplot(size = 0.7, alpha=0.7, outlier.colour = "black", varwidth = T) +
  xlab("")+ ylab("Hamming distances") +
  ggtitle("Pairwise hamming distances of novel alelles discovered by IgDiscover & TIgGER")+
  theme(text = element_text(size = 18))

####create a table showing the sd, mean and medians of pairwise distances within breeds
############ Boran IgDiscover 
alignment_Boran_Igdiscover <- read.dna("./Align_allele/boran_align.fasta", format = "fasta")

hamm_dist_Boran_IgDiscover <- hammingdists(alignment_Boran_Igdiscover)

hamm_dist_Boran_IgDiscover_matrix <- as.matrix(hamm_dist_Boran_IgDiscover)
#extracting distances for pairs defined by the upper triangle of the distance matrix
xy <- t(combn(colnames(hamm_dist_Boran_IgDiscover_matrix), 2))
hamm_dist_Boran_IgDiscover_df <- data.frame(xy, dist=hamm_dist_Boran_IgDiscover_matrix[xy])
View(hamm_dist_Boran_IgDiscover_df)
colnames(hamm_dist_Boran_IgDiscover_df) <- c("allele_X1", "allele_X2", "Boran_Dist")
hamm_dist_Boran_IgDiscover_df$Breed <- "Boran_IgDiscover"

write.csv(hamm_dist_Boran_IgDiscover_df, "./Hamming_Dists/Boran_new.csv")

sd(hamm_dist_Boran_IgDiscover_df$Boran_Dist) # 9.308315
mean(hamm_dist_Boran_IgDiscover_df$Boran_Dist)# 39.66667
median(hamm_dist_Boran_IgDiscover_df$Boran_Dist)# 39
############ Ndama IgDiscover 
alignment_ndama_Igdiscover <- read.dna("./Align_allele/ndama_align.fasta", format = "fasta")

hamm_dist_ndama_IgDiscover <- hammingdists(alignment_ndama_Igdiscover)

hamm_dist_ndama_IgDiscover_matrix <- as.matrix(hamm_dist_ndama_IgDiscover)
#extracting distances for pairs defined by the upper triangle of the distance matrix
xy <- t(combn(colnames(hamm_dist_ndama_IgDiscover_matrix), 2))
hamm_dist_ndama_IgDiscover_df <- data.frame(xy, dist=hamm_dist_ndama_IgDiscover_matrix[xy])
View(hamm_dist_Boran_IgDiscover_df)
colnames(hamm_dist_ndama_IgDiscover_df) <- c("allele_X1", "allele_X2", "Ndama_Dist")
hamm_dist_ndama_IgDiscover_df$Breed <- "Ndama_IgDiscover"

write.csv(hamm_dist_ndama_IgDiscover_df,"./Hamming_Dists/Ndama_new.csv")

sd(hamm_dist_ndama_IgDiscover_df$Ndama_Dist) # 10.71581
mean(hamm_dist_ndama_IgDiscover_df$Ndama_Dist)# 50.6
median(hamm_dist_ndama_IgDiscover_df$Ndama_Dist)# 52

############### Ankole IgDiscover 
alignment_ankole_Igdiscover <- read.dna("./Align_allele/ankole_align.fasta", format = "fasta")

hamm_dist_ankole_IgDiscover <- hammingdists(alignment_ankole_Igdiscover)

hamm_dist_ankole_IgDiscover_matrix <- as.matrix(hamm_dist_ankole_IgDiscover)
#extracting distances for pairs defined by the upper triangle of the distance matrix
xy <- t(combn(colnames(hamm_dist_ankole_IgDiscover_matrix), 2))
hamm_dist_ankole_IgDiscover_df <- data.frame(xy, dist=hamm_dist_ankole_IgDiscover_matrix[xy])
View(hamm_dist_ankole_IgDiscover_df)
colnames(hamm_dist_ankole_IgDiscover_df) <- c("allele_X1", "allele_X2", "Ankole_Dist")
hamm_dist_ankole_IgDiscover_df$Breed <- "Ankole_IgDiscover"
write.csv(hamm_dist_ankole_IgDiscover_df, "./Hamming_Dists/ankole_new.csv")

sd(hamm_dist_ankole_IgDiscover_df$Ankole_Dist) # 12.49
mean(hamm_dist_ankole_IgDiscover_df$Ankole_Dist)# 34
median(hamm_dist_ankole_IgDiscover_df$Ankole_Dist)# 38


################ Friesian IgDiscover
alignment_friesian_Igdiscover <- read.dna("./Align_allele/friesian_align.fasta", format = "fasta")
######
hamm_dist_friesian_IgDiscover <- hammingdists(alignment_friesian_Igdiscover)

hamm_dist_friesian_IgDiscover_matrix <- as.matrix(hamm_dist_friesian_IgDiscover)
#extracting distances for pairs defined by the upper triangle of the distance matrix
xy <- t(combn(colnames(hamm_dist_friesian_IgDiscover_matrix), 2))
hamm_dist_friesian_IgDiscover_df <- data.frame(xy, dist=hamm_dist_friesian_IgDiscover_matrix[xy])
colnames(hamm_dist_friesian_IgDiscover_df) <- c("allele_X1", "allele_X2", "Friesian_Dist")


sd(hamm_dist_friesian_IgDiscover_df$Friesian_Dist) # 8.437102
mean(hamm_dist_friesian_IgDiscover_df$Friesian_Dist)# 31.05445
median(hamm_dist_friesian_IgDiscover_df$Friesian_Dist)# 31


###################### plot different boxplot per breed ################
hamm_dists_per_breed <- read.csv("./Hamming_Dists/combined_Ankole_Boran_Friesian_Ndama.csv")
View(hamm_dists_per_breed)
head(hamm_dists_per_breed)


ggplot(data = hamm_dists_per_breed, aes(x = Breed, y = Dist)) +  geom_jitter(aes(colour =Breed)) + 
  geom_boxplot(size = 0.7, alpha=0.7, outlier.colour = "black", varwidth = T) +
  xlab("")+ ylab("Hamming distances") +
  ggtitle("Pairwise hamming distances of novel alelles discovered by IgDiscover per Breed")

############ Boran TIgGER
alignment_Boran_tigger <- read.dna("./Align_allele/Boran_tigger_align.fasta", format = "fasta")

hamm_dist_Boran_tigger <- hammingdists(alignment_Boran_tigger)

hamm_dist_Boran_tigger_matrix <- as.matrix(hamm_dist_Boran_tigger)
#extracting distances for pairs defined by the upper triangle of the distance matrix
xy <- t(combn(colnames(hamm_dist_Boran_tigger_matrix), 2))
hamm_dist_Boran_tigger_df <- data.frame(xy, dist=hamm_dist_Boran_tigger_matrix[xy])
colnames(hamm_dist_Boran_tigger_df) <- c("allele_X1", "allele_X2", "Dist")
hamm_dist_Boran_tigger_df$Breed <- "Boran_TIgGER"

write.csv(hamm_dist_Boran_tigger_df, "./hamming_dist/Boran_new.csv")

sd(hamm_dist_Boran_tigger_df$Boran_Dist) # 16.13573
mean(hamm_dist_Boran_tigger_df$Boran_Dist)# 50.19048
median(hamm_dist_Boran_tigger_df$Boran_Dist)# 50

############ Ndama TIgGER
alignment_ndama_tigger <- read.dna("./Align_allele/Ndama_tigger_align.fasta", format = "fasta")

hamm_dist_ndama_tigger <- hammingdists(alignment_ndama_tigger)

hamm_dist_ndama_tigger_matrix <- as.matrix(hamm_dist_ndama_tigger)
#extracting distances for pairs defined by the upper triangle of the distance matrix
xy <- t(combn(colnames(hamm_dist_ndama_tigger_matrix), 2))
hamm_dist_ndama_tigger_df <- data.frame(xy, dist=hamm_dist_ndama_tigger_matrix[xy])
colnames(hamm_dist_ndama_tigger_df) <- c("allele_X1", "allele_X2", "Dist")
hamm_dist_ndama_tigger_df$Breed <- "Ndama_TIgGER"
write.csv(hamm_dist_ndama_tigger_df, "./hamming_dist/Ndama_new.csv")

sd(hamm_dist_ndama_tigger_df$Ndama_Dist) # 15.28431
mean(hamm_dist_ndama_tigger_df$Ndama_Dist)# 57.79739
median(hamm_dist_ndama_tigger_df$Ndama_Dist)# 60

################ Friesian TIgGER
alignment_friesian_tigger <- read.dna("./Align_allele/Friesian_tigger_align.fasta", format = "fasta")
######
hamm_dist_friesian_tigger <- hammingdists(alignment_friesian_tigger)

hamm_dist_friesian_tigger_matrix <- as.matrix(hamm_dist_friesian_tigger)
#extracting distances for pairs defined by the upper triangle of the distance matrix
xy <- t(combn(colnames(hamm_dist_friesian_tigger_matrix), 2))
hamm_dist_friesian_tigger_df <- data.frame(xy, dist=hamm_dist_friesian_tigger_matrix[xy])
colnames(hamm_dist_friesian_tigger_df) <- c("allele_X1", "allele_X2", "Dist")
hamm_dist_friesian_tigger_df$Breed <- "Friesian_TIgGER"

write.csv(hamm_dist_friesian_tigger_df, "./hamming_dist/fri_new.csv")

sd(hamm_dist_friesian_tigger_df$Friesian_Dist) # 10.54795
mean(hamm_dist_friesian_tigger_df$Friesian_Dist)# 29
median(hamm_dist_friesian_tigger_df$Friesian_Dist)# 27.5

###################### plot different boxplot per breed ################
hamm_dists_per_breed <- read.csv("./hamming_dist/combined_boran_ndama_friesian.csv")
View(hamm_dists_per_breed)
head(hamm_dists_per_breed)


ggplot(data = hamm_dists_per_breed, aes(x = Breed, y = Dist)) +  geom_jitter(aes(colour =Breed)) + 
  geom_boxplot(size = 0.7, alpha=0.7, outlier.colour = "black", varwidth = T) +
  xlab("")+ ylab("Hamming distances") +
  ggtitle("Pairwise hamming distances of novel alelles discovered by TIgGER per Breed")


######## plot different boxplot per breed two tools combined #####

hamm_dists_per_breed_tigger_igdiscover <- read.csv("../Combined_Igdiscover_Tigger/hamming_dist_plots/combined_ankole_boran_ndama_friesian_tigger_igdiscover.csv")
head(hamm_dists_per_breed_tigger_igdiscover)

ggplot(data = hamm_dists_per_breed_tigger_igdiscover, aes(x = Breed, y = Dist)) +  geom_jitter(aes(colour =Breed)) + 
  geom_boxplot(size = 0.7, alpha=0.7, outlier.colour = "black", varwidth = T) +
  xlab("")+ ylab("Hamming distances") +
  ggtitle("Pairwise hamming distances of novel alelles discovered per Breed")+
  theme(text = element_text(size = 14))



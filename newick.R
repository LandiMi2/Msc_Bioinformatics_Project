setwd("/home/cofia/Documents/Msc_Bioinformatics_Project/IgDiscover/")
setwd("~/Documents/Msc_Bioinformatics_Project/TIgGER/Results/Bovine_Novel_fasta/")
setwd("/home/cofia/Documents/Msc_Bioinformatics_Project/Combined_Igdiscover_Tigger/")

library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)

seqs <- readDNAStringSet("./Ankole_combined.fasta", format = "fasta")
# look at some of the sequences (optional)
# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs)

# write the alignment to a new FASTA file
writeXStringSet(aligned,
                file="combined_ankole.fasta")

# read in the aligned data
align <- read.alignment("combined_ankole.fasta", format = "fasta")

# create a distance matrix for the alignment 
D <- dist.alignment(align, matrix = "similarity")

h_cluster <- hclust(D, method = "complete", members = NULL)



class(h_cluster) # must be hclust class

my_tree <- as.phylo(h_cluster) 

#create newick file 
write.tree(phy=my_tree, file="combined_ankole.newick") # look for the file in your working directory
                  
                  
                  

library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)
setwd("~/Documents/Msc_Bioinformatics_Project/IgDiscover/")
############################ J gene analysis #######################################
# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet("./combined_V_new_germline_boran_ndama_ankole.fasta", format = "fasta")
# look at some of the sequences (optional)
# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs <- OrientNucleotides(seqs)

# perform the alignment
aligned <- AlignSeqs(seqs)

# view the alignment in a browser (optional)
BrowseSeqs(aligned, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned,
                file="combined_aligned.fasta")

# read in the aligned data
combined_align <- read.alignment("combined_aligned.fasta", format = "fasta")

# create a distance matrix for the alignment 
D <- dist.alignment(combined_align, matrix = "similarity")

#heatmap(as.matrix(D))
temp <- as.data.frame(as.matrix(D))
View(temp)
# visualizing the distance matrix option 1
D_matrix <- data.matrix(temp)
nrow(D_matrix)

image(1:nrow(D_matrix), 1:nrow(D_matrix), D_matrix, axes = F, xlab="", ylab="", main = "Germline J gene distance matrix")  
axis(1, 1:nrow(D_matrix), colnames(temp), cex.axis = 0.5, las=3) 
axis(2, 1:nrow(D_matrix),colnames(temp), cex.axis = 0.5, las=2) 
text(expand.grid(1:nrow(D_matrix), 1:nrow(D_matrix)), sprintf("%0.1f", D_matrix), cex=0.6)

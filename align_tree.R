# install necessary packges 
install.packages("remotes")
library(remotes)
remotes::install_github("GuangchuangYu/treeio", force = TRUE)
BiocManager::install("DECIPHER")
BiocManager::install(version = "3.10")

library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)
setwd("~/Documents/Msc_Bioinformatics_Project/Data/Annotatio_Output/Alignments_Tree/")
############################ J gene analysis #######################################
# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readDNAStringSet("./J.fasta", format = "fasta")
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
                file="J_aligned.fasta")

# read in the aligned data
J_align <- read.alignment("J_aligned.fasta", format = "fasta")

# create a distance matrix for the alignment 
D <- dist.alignment(J_align, matrix = "similarity")

#heatmap(as.matrix(D))
temp <- as.data.frame(as.matrix(D))
# visualizing the distance matrix option 1
D_matrix <- data.matrix(temp)
nrow(D_matrix)

image(1:18, 1:18, D_matrix, axes = F, xlab="", ylab="", main = "Germline J gene distance matrix")  
axis(1, 1:18, colnames(temp), cex.axis = 0.5, las=3) 
axis(2, 1:18,colnames(temp), cex.axis = 0.5, las=2) 
text(expand.grid(1:18, 1:18), sprintf("%0.1f", D_matrix), cex=0.6)

#visualizing the distance matrix option 2 - using networks 
dist_mi <- 1/D_matrix # one over, as qgraph takes similarity matrices as input
library(qgraph)
j_gene_network <- qgraph(dist_mi, layout='spring', vsize=7,
                         title = "Germline J gene distance matrix network")

############################# V gene analysis ##############################################
seqs_V <- readDNAStringSet("./V.fasta", format = "fasta")
# look at some of the sequences (optional)
# nucleotide sequences need to be in the same orientation
# if they are not, then they can be reoriented (optional)
seqs_V <- OrientNucleotides(seqs_V)

# perform the alignment
aligned_V <- AlignSeqs(seqs_V)

# view the alignment in a browser (optional)
BrowseSeqs(aligned_V, highlight=0)

# write the alignment to a new FASTA file
writeXStringSet(aligned_V,
                file="V_aligned.fasta")

# read in the aligned data
V_align <- read.alignment("V_aligned.fasta", format = "fasta")

# create a distance matrix for the alignment 
D_V <- dist.alignment(V_align, matrix = "similarity")

#heatmap(as.matrix(D))
temp_V <- as.data.frame(as.matrix(D_V))
View(temp_V)
# visualizing the distance matrix option 1
D_matrix_V <- data.matrix(temp_V)

nrow(D_matrix_V)

image(1:nrow(D_matrix_V), 1:nrow(D_matrix_V), D_matrix_V, axes = F, xlab="", ylab="", main = "Germline V gene distance matrix")  
axis(1, 1:nrow(D_matrix_V), colnames(temp_V), cex.axis = 0.5, las=3) 
axis(2, 1:nrow(D_matrix_V),colnames(temp_V), cex.axis = 0.5, las=2) 
text(expand.grid(1:nrow(D_matrix_V), 1:nrow(D_matrix_V)), sprintf("%0.1f", D_matrix_V), cex=0.6)

#visualizing the distance matrix option 2 - using networks 
dist_mi_v <- 1/D_matrix_V # one over, as qgraph takes similarity matrices as input
library(qgraph)
V_gene_network <- qgraph(dist_mi_v, layout='spring', vsize=7,
                         title = "Germline V gene distance matrix network")



##### Draw a phylogenetic tree ########


##### dendrogramm for germline J gene
h_cluster <- hclust(D, method = "average", members = NULL) # method = average is used for UPGMA, members can be equal to NULL or a vector with a length of size D
plot(h_cluster, cex = 0.8, main = "Germline J gene dendrogram")

######## dendrogram for gemline V gene 

h_cluster_V <- hclust(D_V, method = "average", members = NULL) # method = average is used for UPGMA, members can be equal to NULL or a vector with a length of size D
plot(h_cluster_V, cex = 0.8, main = "Germline V gene dendrogram")

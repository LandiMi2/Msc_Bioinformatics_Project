

############################# Creating confusion matrixs ############################################
####################### IgBlast V matrix ###########################
data_Igblast_Vgene <- data.frame(comparison_Vcall_Igblast$IgBlast_Vcall,comparison_Vcall_Igblast$True_Vcall,
                                 comparison_Vcall_Igblast$HIT)
#columnames
colnames(data_Igblast_Vgene) <- c("IgBlast_Vcall", "True_Vcall", "HIT")
#View(data_Igblast_Vgene)
str(data_Igblast_Vgene)
#add IGHV1-33 as a level in IgBlast Vcall
#below is not my function....

addLevel <- function(x, newlevel=x[10]) {
  if(is.factor(x)) {
    if (is.na(match(newlevel, levels(x))))
      return(factor(x, levels=c(levels(x), newlevel)))
  }
  return(x)
}

data_Igblast_Vgene$IgBlast_Vcall <- addLevel(data_Igblast_Vgene$IgBlast_Vcall, "IGHV1-33")
data_Igblast_Vgene$True_Vcall <- addLevel(data_Igblast_Vgene$True_Vcall, "IGHV1-33")
#### draw a hit map for V gene calls #############################
#construct a table for comaparison for Igblast V gene 
Igblast_Vgene_CrossTable <- data.frame(table(data_Igblast_Vgene$IgBlast_Vcall, data_Igblast_Vgene$True_Vcall))
#change column names 
colnames(Igblast_Vgene_CrossTable) <- c("Igblast_Vcall", "True_Vcall", "Counts")


#log transform
Igblast_Vgene_CrossTable$Counts <- log(Igblast_Vgene_CrossTable$Counts)
str(Igblast_Vgene_CrossTable)
levels(Igblast_Vgene_CrossTable$Igblast_Vcall)
#force the levels to be ordered the same as other matrix
Igblast_Vgene_CrossTable$True_Vcall <- factor(Igblast_Vgene_CrossTable$True_Vcall, 
                                                 levels = c("IGHV1-10", "IGHV1-14", "IGHV1-17", "IGHV1-20", "IGHV1-21", "IGHV1-25", "IGHV1-27", "IGHV1-30",
                                                            "IGHV1-32","IGHV1-33", "IGHV1-37", "IGHV1-39", "IGHV1-7" , "IGHV1S1" , "IGHV2S1"))
 
library(ggplot2)
library(gridExtra)

heatmap_Vgene_Igblast <- ggplot(Igblast_Vgene_CrossTable, aes(Igblast_Vcall, True_Vcall)) +
  geom_tile(aes(fill = Counts), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red", name = "log(Counts)") + 
  ggtitle("Distribution of hits and mishits of IgBlast") +
  xlab("IgBlast V calls") + ylab("True V calls") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(aspect.ratio=1)


################################################ MiXCR V matrix ############################

data_mixcr_Vgene <- data.frame(compare_mixcr_Vcall$mixcr.Vcall, compare_mixcr_Vcall$read_vdj_mixcr.V_call,
                               compare_mixcr_Vcall$HIT)
colnames(data_mixcr_Vgene) <- c("MiXCR_Vcall","True_Vcall", "HIT")
str(data_mixcr_Vgene)
#add levels to make it a square matrix
levels(data_mixcr_Vgene$True_Vcall)
levels(data_mixcr_Vgene$MiXCR_Vcall)

data_mixcr_Vgene$MiXCR_Vcall <- addLevel(data_mixcr_Vgene$MiXCR_Vcall, "IGHV1S1")
data_mixcr_Vgene$MiXCR_Vcall <- addLevel(data_mixcr_Vgene$MiXCR_Vcall, "IGHV2S1")
data_mixcr_Vgene$True_Vcall <- addLevel(data_mixcr_Vgene$True_Vcall, "IGHV1-33")
str(data_mixcr_Vgene)
#construct a table for comaparison for MiXCR V gene 
mixcr_Vgene_CrossTable <- data.frame(table(data_mixcr_Vgene$MiXCR_Vcall, data_mixcr_Vgene$True_Vcall))

#change column names 
colnames(mixcr_Vgene_CrossTable) <- c("MiXCR_Vcall", "True_Vcall", "Counts")

#log transformation
mixcr_Vgene_CrossTable$Counts <- log(mixcr_Vgene_CrossTable$Counts)
#View(mixcr_Vgene_CrossTable)
levels(mixcr_Vgene_CrossTable$MiXCR_Vcall)
#modify levels
mixcr_Vgene_CrossTable$True_Vcall <- factor(mixcr_Vgene_CrossTable$True_Vcall, 
                                                 levels = c("IGHV1-10", "IGHV1-14", "IGHV1-17", "IGHV1-20", "IGHV1-21", "IGHV1-25", "IGHV1-27", "IGHV1-30",
                                                            "IGHV1-32","IGHV1-33", "IGHV1-37", "IGHV1-39", "IGHV1-7" , "IGHV1S1" , "IGHV2S1"))

heatmap_Vgene_mixcr <- ggplot(mixcr_Vgene_CrossTable, aes(MiXCR_Vcall, True_Vcall)) +
  geom_tile(aes(fill = Counts), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red", name = "log(Counts)") +
  ggtitle("Frequency of hits and mishits of MiXCR") +
  ylab("True V calls") + xlab("MiXCR V calls") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(aspect.ratio=1)

############################ IMGT V matrix ..............................
data_imgt_Vgene <- data.frame(comparison_Vcall_imgt$IMGT_Vcall, comparison_Vcall_imgt$True_Vcall,
                              comparison_Vcall_imgt$HIT)
colnames(data_imgt_Vgene) <- c("IMGT_Vcall","True_Vcall", "HIT")

str(data_imgt_Vgene)
levels(data_imgt_Vgene$True_Vcall)
#data_imgt_Vgene$IMGT_Vcall <- addLevel(data_imgt_Vgene$IMGT_Vcall, "IGHV1-33")
#to remove 466 empty rows from this datafram IMGT
data_imgt_Vgene <- data_imgt_Vgene[-which(data_imgt_Vgene$IMGT_Vcall == ""), ]
sum(data_imgt_Vgene$IMGT_Vcall == "") #confirm 
#add level
data_imgt_Vgene$IMGT_Vcall <- addLevel(data_imgt_Vgene$IMGT_Vcall, "IGHV1-33")
data_imgt_Vgene$True_Vcall <- addLevel(data_imgt_Vgene$True_Vcall, "IGHV1-33")

#construct a table for comaparison for IMGT V gene 
imgt_Vgene_CrossTable <- data.frame(table(data_imgt_Vgene$IMGT_Vcall, data_imgt_Vgene$True_Vcall))
#change column names 
colnames(imgt_Vgene_CrossTable) <- c("IMGT_Vcall", "True_Vcall", "Counts")
#log transformation
imgt_Vgene_CrossTable$Counts <- log(imgt_Vgene_CrossTable$Counts)
#remove empty rows 
imgt_Vgene_CrossTable <- imgt_Vgene_CrossTable[-which(imgt_Vgene_CrossTable$IMGT_Vcall == ""), ]

levels(imgt_Vgene_CrossTable$True_Vcall) 
#rearrange levels 
imgt_Vgene_CrossTable$True_Vcall <-  factor(imgt_Vgene_CrossTable$True_Vcall , 
                                            levels = c("IGHV1-10", "IGHV1-14", "IGHV1-17", "IGHV1-20", "IGHV1-21", "IGHV1-25", "IGHV1-27", "IGHV1-30",
                                                       "IGHV1-32","IGHV1-33", "IGHV1-37", "IGHV1-39", "IGHV1-7" , "IGHV1S1" , "IGHV2S1"))

heatmap_Vgene_imgt <- ggplot(imgt_Vgene_CrossTable, aes(IMGT_Vcall, True_Vcall)) +
  geom_tile(aes(fill = Counts), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red", name = "log(Counts)") +
  ggtitle("Frequency of hits and mishits of IMGT") +
  xlab("IMGT V calls") + ylab("True V calls") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(aspect.ratio=1)

##### combine the three heat maps ################
grid.arrange(heatmap_Vgene_Igblast,heatmap_Vgene_mixcr, heatmap_Vgene_imgt, ncol = 3)

grid.arrange(heatmap_Vgene_Igblast,heatmap_Vgene_mixcr, heatmap_Vgene_imgt, ncol = 2)



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

heatmap_Jgene_mixcr <- ggplot(mixcr_Jgene_CrossTable, aes(MiXCR_Jcall,True_Jcall,)) +
  geom_tile(aes(fill = Counts), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red", name = "log(Counts)") + 
  ggtitle("Frequency of hits and mishits of MiXCR J gene annotation") +
  ylab("True J calls") + xlab("MiXCR J calls")


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

heatmap_Jgene_imgt <- ggplot(imgt_Jgene_CrossTable, aes(IMGT_Jcall, True_Jcall)) +
  geom_tile(aes(fill = Counts), colour = "black") +
  scale_fill_gradient(low = "blue", high = "red", name = "log(Counts)") + 
  ggtitle("Frequency of hits and mishits of IMGT J gene annotation") +
  ylab("True J calls") + xlab("IMGT J calls")


##### combine the three heat maps ################
heat_list_Jgene <- list(heatmap_Jgene_Igblast, heatmap_Jgene_imgt, heatmap_Jgene_mixcr)
heat_list[["ncol"]] <- 1
do.call(grid.arrange, heat_list_Jgene)

setwd("./Documents/Msc_Bioinformatics_Project/IgDiscover/Novel_allele_Igdiscover/")

library(readxl)
library(gridExtra)
###### Boran
Boran_alleles <- read_excel("./Boran_alleles.xlsx")
#order using number of Jgenes counts 
Boran_alleles_ordered <- Boran_alleles[order(-Boran_alleles[,9]), ]
View(Boran_alleles_ordered)
write.csv(Boran_alleles_ordered, "Boran_alleles_ordered.csv")
#save as png table 
png(filename = "Boran_alleles_ordered_table.png", width=1200,height=500,bg = "white")
grid.table(Boran_alleles_ordered)
dev.off()

##### Ndama
ndama_alleles <- read_excel("./Ndama_alleles.xlsx")
#order
ndama_alleles_ordered <- ndama_alleles[order(-ndama_alleles[,5]), ]
write.csv(ndama_alleles_ordered, "ndama_alleles_ordered.csv")
#save as png table 
png(filename = "Ndama_alleles_ordered_table.png", width=1200,height=200,bg = "white")
grid.table(ndama_alleles)
dev.off()

#### Ankole 
ankole_alleles <- read_excel("./Ankole_alleles.xlsx")
#order
ankole_alleles_ordered <- ankole_alleles[order(-ankole_alleles[,9]), ]
write.csv(ankole_alleles_ordered, "Ankole_alleles_ordered.csv")
#save as png table 
png(filename = "Ankole_alleles_ordered_table.png", width=1200,height=200,bg = "white")
grid.table(ankole_alleles)
dev.off()



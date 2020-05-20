setwd("~/Documents/Msc_Bioinformatics_Project/IgDiscover")
# Boran annotated_V_germline file 
annotated_V_germline_boran_10 <- read.table("./Bovine_expt_10_iterations_Boran/iteration-10/annotated_V_germline.tab",
                                     header = TRUE, sep = "\t")

#filter the novel alleles identified by IgDiscover
library(tidyr)
novel_boran_10 <- dplyr::filter(annotated_V_germline_boran_10, grepl("^0" , annotated_V_germline_boran_10$is_filtered))
View(novel_boran_10)
#write.table(novel_boran_10, "Novel_Boran.tab", sep = "\t")

# Ndama annotated_V_germline file 
annotated_V_germline_ndama_10 <- read.table("./Bovine_expt_10_iteration_Ndama/iteration-10/annotated_V_germline.tab",
                                            header = TRUE, sep = "\t")
#filter the novel alleles identified by IgDiscover
novel_ndama_10 <- dplyr::filter(annotated_V_germline_ndama_10, grepl("^0" , annotated_V_germline_ndama_10$is_filtered))

View(novel_ndama_10)
#write.table(novel_ndama_10, "Novel_Ndama.tab", sep = "\t")
# Friesian annotated_V_germline file 
annotated_V_germline_friesian_10 <- read.table("./Bovine_expt_10_iteration_Freisian/iteration-10/annotated_V_germline.tab",
                                            header = TRUE, sep = "\t")
#filter the novel alleles identified by IgDiscover
novel_friesian_10 <- dplyr::filter(annotated_V_germline_friesian_10, grepl("^0" , annotated_V_germline_friesian_10$is_filtered))

View(novel_friesian_10)
#write.table(novel_friesian_10, "Novel_friesian.tab", sep = "\t", col.names = TRUE)

# Ankole annoatated_V_germline file 
annotated_V_germline_ankole_10 <- read.table("./Bovine_expt_10_iteration_Ankole/iteration-10/annotated_V_germline.tab",
                                               header = TRUE, sep = "\t")
#filter the novel alleles identified by IgDiscover
novel_ankole_10 <- dplyr::filter(annotated_V_germline_ankole_10, grepl("^0" , annotated_V_germline_ankole_10$is_filtered))

View(novel_ankole_10)
#write.table(novel_ankole_10, "Novel_Ankole.tab", sep = "\t")

setwd("~/Documents/Msc_Bioinformatics_Project/IgDiscover")
#open annotated V germline 4th iteration 

annotated_V_germline_10 <- read.table("./Bovine_expt_10_iterations/iteration-10/annotated_V_germline.tab",
                                     header = TRUE, sep = "\t")

#annotated_final <- read.table("./Bovine_expt_4_iterations/final/assigned.tab.gz",
 #                             header = TRUE, sep = "\t")

View(annotated_V_germline_10)

#filter the novel alleles identified by IgDiscover
library(tidyr)
novel_10 <- dplyr::filter(annotated_V_germline_10, grepl("^0" , annotated_V_germline_10$is_filtered))
View(novel_10)

#filter based on database difference
novel_10_database_change <- dplyr::filter(novel_10, grepl("^[1-9]" , novel_10$database_diff))
View(novel_10_database_change)


##### check out marcels file 
annotated_V_germline_marcel <- read.table("./annotated_V_germline_marcel.tab",
                                          header = TRUE, sep = "\t")
novel_marcel <- dplyr::filter(annotated_V_germline_marcel, grepl("^0", annotated_V_germline_marcel$is_filtered))



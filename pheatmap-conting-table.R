### Reed Woyda
# reedwoyda@gmail.com
# 11/17/2020
# Script to make pheatmaps with metadata

library(pheatmap)
library(viridis)

# data must be in the format: rownames = isolate, colnames = genes
# metadata must be in the format: rownames = isolates (*have to be in same order as data rownames), colnames = metadata categories
data <- as.matrix(read.csv("db_element_ID_vfdb-conting-table.csv", row.names = 1))
metadata <- as.data.frame(read.csv("metadata_794.csv", row.names = 1))

# for(i in 1:length(colnames(data))){
#   for(j in 1:length(rownames(data))){
#     if(data[j,i] == ""){
#       data[j,i] <- 0
#     }  
#     else{
#       data[j,i] <- 1
#     }
#       
#     }
# }
#data <- as.data.frame(data)
# str(data)
# sapply(data[1:794, 1:252], as.integer)
# str(data)
# class(data) <- "numeric"

# data.t.copy <- as.data.frame(data.t)
# allzero = c()
# j=1
# for(i in 5260:1){
#   if(dim(table(data.t[,i])) != 1)
#   {
#     print(i)
#     allzero[j] <- i
#     j = j + 1
#     data.t.copy = subset(data.t.copy, select = i)
#   }
#     
# }

#rownames(data.t) <- rownames(metadata)
colnames(metadata)[3] <- "Isolation Source"
#colnames(metadata)[2] <- "Combined Plasmids"
heatmap(data)
pheatmap(data)

pheatmap(
  mat               = data,
  color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  drop_levels       = FALSE,
  fontsize          = 14,
  annotation_row    = metadata,
  main              = "Virulence Finder Database Presence/Absence Heatmap",
  legend_breaks     = c(0,1),
  legend_labels     = c("0", "1")

)





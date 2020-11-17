### Reed Woyda
# 11/17/2020
# Script for doing hierarchical clustering on contingency tables
# Mainly uses: Distance() for distance matrix calculation, fviz_nbclust() for silhouette optimal cluster number estimation, 
#    and hclust() for hierarchical clustering

library("fastcluster")
library("dendextend")
library("ggplot2")
library("ggdendro")
library("IntClust")
library("ggtree")
library("factoextra")
library(cluster)

### Read in data in format: rownames = isolate, colnames = genes
# metadata must be in the format: rownames = isolates (*have to be in same order as data rownames), colnames = metadata categories
data <- as.matrix(read.csv("example_contingency_table.csv", row.names = 1))
metadata <- as.data.frame(read.csv("metadata_example.csv", row.names = 1))


### Calculate distance matrix using the coefficiant of Rogers and Tanimoto
# http://www.imsbio.co.jp/RGM/R_rdfile?f=IntClust/man/Distance.Rd&d=R_CC
# "euclidean" for abundance
# "tanimoto" for presence/absence
distance.mat <- Distance(data,distmeasure = "tanimoto", normalize = FALSE, method = NULL )

### hclust through fastcluster
# https://cran.r-project.org/web/packages/fastcluster/fastcluster.pdf
# "complete" for abundance
# "ward.D" for pres/absence
hclust.results <- hclust(as.dist(distance.mat), method = "ward.D", members = NULL)
hclust.results.den <- as.dendrogram(hclust.results)


### Estimate optimal number of clusters using the Silhouette method
fviz_nbclust(data, kmeans, method='silhouette', k.max = 10)

### Creating dendrogram
# https://www.researchgate.net/post/How_do_i_perform_a_cluster_analysis_on_a_very_large_data_set_in_R_Which_will_be_the_best_complete_or_single_linkage_method
# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
# https://rdrr.io/cran/dendextend/man/rect.dendrogram.html
# replace k = x with number of optimal clusters obtained from silhouette method above
sub_grp1 <- cutree(hclust.results, k = 2)
table(sub_grp1)
plot(hclust.results, cex = 0.6, main = "RGI and VFDB cluster dendrogram", xlab = "Isolates colored by cluster") 
rect.hclust(hclust.results, k = 2, border = 2:5)
hclust.results <- data.frame(sub_grp1)
# Optional: output table of isolates with their respective cluster #
write.csv(sub_grp1, "RGI_VFDB_cluster-2-11172020.csv")

### Better dendrogram uisn

hclust.results.den %>% color_branches(k=2) %>% plot
# Change color and size for labels
hclust.results.den %>%  color_branches(k=2) %>%# Change size
  plot(main = "RGI and VFDB cluster dendrogram") # plot



### The below code is to add metadata labels onto the dendrogram
# To do this you need a metadata CSV with the rownames as the isolate names, in the same order as the data, and the colnames can contain metadata categories
# Here for example it is using the "isolation_source" metadata
# can edit "labels_col" to color the metadata categories
# To be honest, the labeling is not optimal and needs to be worked on to get a better representation of the metadata distribution along the x-axis of the dendrogram
hclust.results.den %>% set("labels", c(metadata$Host)) %>% set("labels_col", c("blue","darkorange4")) %>% 
  set("labels_cex", 0.5) %>% 
  set("labels", c(metadata$Host)) %>% color_branches(k=2) %>%
  plot(main = "RGI and VFDB cluster dendrogram")


# Some example colors
# c("green", "blue","red","black","aquamarine","chartreuse","chocolate1","darkgoldenrod1","darkmagenta","darkorange4","firebrick4","darksalmon","darkslateblue")



























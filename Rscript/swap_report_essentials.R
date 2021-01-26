library(reshape2)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(viridis)


#functions to load jaccard matrix

getSampleRunData <- function(jaccard_matrix){
  #"input: jaccard matrix dataframe"
  #"output: dataframe for sample, donor info"
  head(jaccard_matrix)
  jm <- jaccard_matrix[,1:ncol(jaccard_matrix)-1]
  
  # create a sample id map
  sample_miso_donor_map <- data.frame(cbind(names(jm)[2:length(names(jm))], stringr::str_split_fixed(stringr::str_replace(stringr::str_replace(names(jm)[2:length(names(jm))], ".annotated",""),"[.]", "-"), "_",10)[,3:10]))
  sample_miso_donor_map$miso_donorID <- paste(sample_miso_donor_map$X2, sample_miso_donor_map$X3, sep = "_")
  sample_miso_donor_map$effective_donorID <- stringr::str_split_fixed(sample_miso_donor_map$X9, "_",3)[,1]
  
  donor_maps <- unique(sample_miso_donor_map[,c("miso_donorID", "effective_donorID")])
  donor_maps <- donor_maps[grepl("MATS", donor_maps$effective_donorID),]
  row.names(donor_maps) <- donor_maps$miso_donorID
  
  sample_miso_donor_map$effective_donorID <- donor_maps[sample_miso_donor_map$miso_donorID,]$effective_donorID
  
  sample_miso_donor_map$workingSampleID <- stringr::str_split_fixed(sample_miso_donor_map$X9, "-", 1)[,1]
  sample_miso_donor_map$Expt <- ifelse(grepl("WG", sample_miso_donor_map$X1), "WG", "EX/TS")
  
  return (sample_miso_donor_map)
}



# loadJaccardMatrix <- function(jaccard_matrix_file){
#   # function loads jaccard matrix as a data frame
#   jaccard_matrix <- jaccard_matrix_file
#   laneData <- read.delim(file=jaccard_matrix, header = T)
#   sample_miso_donor_map <- getSampleRunData(laneData)
#   # check for the effective_DonorId
#   sample_miso_donor_map$effective_donorID <- sample_miso_donor_map$miso_donorID
#   sample_miso_donor_map$workingSampleID <- paste0(sample_miso_donor_map$miso_donorID, "_", sample_miso_donor_map$workingSampleID)
# 
#   # number of donors
#   tot.donors <- length(unique(sample_miso_donor_map$effective_donorID))
#   tot.libs <- dim(sample_miso_donor_map)[1]
#   
#   laneData <- laneData[,!(colnames(laneData) %in% c("X", "SNPs"))]
#   row.names(laneData) <- sample_miso_donor_map$workingSampleID
#   colnames(laneData) <- sample_miso_donor_map$workingSampleID
#   
#   return(laneData)
# }


drawHeatMap <- function(data_frame, 
                        annotation_df, 
                        tree_cuts, 
                        filename, 
                        breaksList, 
                        height = 30, 
                        width = 40, 
                        custom_colors, 
                        showAnnotations = F, 
                        treeheight = 0, 
                        fontsize = 12){
  # draw and save heatmap 
  pheatmap(as.matrix(data_frame), 
           color = colorRampPalette(rev(brewer.pal(n = 10, 
                                                   name ="RdYlBu")))(length(breaksList)),
           breaks = breaksList,
           annotation_row = annotation_df,
           annotation_col = annotation_df,
           cluster_rows = T,
           cluster_cols = T,
           show_rownames = showAnnotations,
           show_colnames = showAnnotations,
           cutree_cols = tree_cuts,
           cutree_rows = tree_cuts,
           treeheight_row = treeheight, 
           treeheight_col = treeheight,
           filename = filename,
           height = height, width = width,
           annotation_colors = custom_colors,
           fontsize = fontsize)
}




loadRelevantJMData <- function(jaccard_matrix_file){
  # relevantJMData <- c()
  jaccard_matrix <- jaccard_matrix_file
  laneData <- read.delim(file=jaccard_matrix, header = T)
  sample_miso_donor_map <- getSampleRunData(laneData)
  # check for the effective_DonorId
  sample_miso_donor_map$effective_donorID <- sample_miso_donor_map$miso_donorID
  sample_miso_donor_map$workingSampleID <- paste0(sample_miso_donor_map$miso_donorID, "_", sample_miso_donor_map$workingSampleID)
  # head(sample_miso_donor_map)
  # knitr::kable(sample_miso_donor_map)
  
  # number of donors
  tot.donors <- length(unique(sample_miso_donor_map$effective_donorID))
  tot.libs <- dim(sample_miso_donor_map)[1]
  
  laneData <- laneData[,!(colnames(laneData) %in% c("X", "SNPs"))]
  row.names(laneData) <- sample_miso_donor_map$workingSampleID
  colnames(laneData) <- sample_miso_donor_map$workingSampleID
  relevantJMData <- c("laneData" = laneData,
                      "total.donors" = tot.donors,
                      "total.libraries" = tot.libs,
                      "sample.miso.map" = sample_miso_donor_map)
  return(relevantJMData)
  
}

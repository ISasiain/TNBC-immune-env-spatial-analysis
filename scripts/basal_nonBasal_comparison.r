rm(list = ls())
setwd("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots")

# Installing packages if requitred
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("gtools", quietly = TRUE))
  install.packages("gtools")

if (!require("circlize", quietly = TRUE))
  install.packages("circlize")


#Loading required packages
library(ggplot2)
library(ComplexHeatmap)
library(gtools)
library(circlize)


#
# Reading and transforming the annotation file
#

annotation <- read.csv2("/home/Illumina/Iñaki_Sasiain/immune_spatial/annotation/supplData_withimages.csv")
rownames(annotation) <- annotation$uid

#Grouping non-basal PAM subtypes into a single one
annotation$"PAM50_basal_nonBasal" <- annotation$"PAM50_NCN"
annotation$"PAM50_basal_nonBasal"[annotation$"PAM50_basal_nonBasal" == "unclassified" ] <- NA
annotation$"PAM50_basal_nonBasal"[annotation$"PAM50_basal_nonBasal" != "Basal" ] <- "NonBasal"

#Grouping non LAR molecular subtypes into a single one
annotation$"LAR_nonLAR" <- annotation$"TNBCtype4_n235_notPreCentered"
annotation$"LAR_nonLAR"[annotation$"LAR_nonLAR" != "LAR" ] <- "nonLAR"


#
# Plottting subgroups versus TILs
#
annotation$"PAM50_basal_nonBasal" <- as.factor(annotation$"PAM50_basal_nonBasal")

basal_nonBasal <- annotation$"PAM50_basal_nonBasal"[!is.na(annotation$"PAM50_basal_nonBasal")]
TILs <- annotation$"TILs"[!is.na(annotation$"PAM50_basal_nonBasal")]


# Difference in TILs between basal and non basal

to_plot <- as.data.frame(
  cbind(
    "basal_nonBasal" <- basal_nonBasal,
    "TILs" <- TILs
  )
)

ggplot(data=to_plot,
       aes(x=basal_nonBasal,
           y=TILs,
           fill=basal_nonBasal)) +
       geom_violin() +
       geom_boxplot(width=0.1, fill=c("grey")) +
       xlab("PAM50 subtype") +
       ylab("TILs (%)") +
       theme_classic() +
       scale_fill_discrete(name="PAM50 subtype")


# Difference in TILs within basal TNBC samples

IM_status <- as.factor(annotation[annotation$"PAM50_basal_nonBasal" == "Basal" ,"TNBCtype_IMpositive"])
TILs <- annotation[annotation$"PAM50_basal_nonBasal" == "Basal" ,"TILs"]

to_plot <- data.frame(
  
    IM_status,
    TILs
)

to_plot <- to_plot[!is.na(to_plot$"IM_status"), ]

ggplot(data=to_plot,
       aes(x=IM_status,
           y=TILs,
           fill=IM_status)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill=c("grey")) +
  xlab("PAM50 subtype") +
  ylab("TILs (%)") +
  theme_classic() +
  scale_x_discrete(labels=c("IM-", "IM+"))
  scale_fill_discrete(name="PAM50 subtype")



  # Difference in TILs within non basal TNBC samples: IM- and IM+
  
  IM_status <- as.factor(annotation[annotation$"PAM50_basal_nonBasal" == "NonBasal" ,"TNBCtype_IMpositive"])
  TILs <- annotation[annotation$"PAM50_basal_nonBasal" == "NonBasal" ,"TILs"]
  
  to_plot <- data.frame(
    
    IM_status,
    TILs
  )
  
  to_plot <- to_plot[!is.na(to_plot$"IM_status"), ]
  
  ggplot(data=to_plot,
         aes(x=IM_status,
             y=TILs,
             fill=IM_status)) +
    geom_violin() +
    geom_boxplot(width=0.1, fill=c("grey")) +
    xlab("PAM50 subtype") +
    ylab("TILs (%)") +
    theme_classic()

  # Difference in TILs within non basal TNBC samples: LAR and nonLAR
  
  LAR_nonLAR <- annotation[annotation$"PAM50_basal_nonBasal" == "NonBasal" ,"TNBCtype4_n235_notPreCentered"]
  LAR_nonLAR[LAR_nonLAR != "LAR"] <- "notLAR"
  
  TILs <- annotation[annotation$"PAM50_basal_nonBasal" == "NonBasal" ,"TILs"]
  
  to_plot <- data.frame(
    
    LAR_nonLAR,
    TILs
  )
  
  to_plot <- to_plot[!is.na(to_plot$"LAR_nonLAR"), ]
  
  ggplot(data=to_plot,
         aes(x=LAR_nonLAR,
             y=TILs,
             fill=LAR_nonLAR)) +
    geom_violin() +
    geom_boxplot(width=0.1, fill=c("grey")) +
    xlab("Molecular subtypes") +
    ylab("TILs (%)") +
    theme_classic()
  
  
  
  # Difference in TILs within basal TNBC samples: IM- and IM+
  
  IM_status <- as.factor(annotation[annotation$"PAM50_basal_nonBasal" == "NonBasal" ,"TNBCtype_IMpositive"])
  TILs <- annotation[annotation$"PAM50_basal_nonBasal" == "NonBasal" ,"TILs"]
  
  to_plot <- data.frame(
    
    IM_status,
    TILs
  )
  
  to_plot <- to_plot[!is.na(to_plot$"IM_status"), ]
  
  ggplot(data=to_plot,
         aes(x=IM_status,
             y=TILs,
             fill=IM_status)) +
    geom_violin() +
    geom_boxplot(width=0.1, fill=c("grey")) +
    xlab("PAM50 subtype") +
    ylab("TILs (%)") +
    theme_classic()
  



## Reading metrics
  
  reference_phen <- c("CD3", "CD8", "CD20", "p53", "FOXP3", "CD68", "CD4")
  target_phen <- c("CD3", "CD8", "CD20", "p53", "FOXP3", "CD68", "CD4")
  
  path_to_metrics <- "/home/Illumina/Iñaki_Sasiain/immune_spatial/spatial_metrics/"
  
  list_of_metrics <- list()



  #Store the metrics of the markers of interest in a list
  for (marker1 in reference_phen) {
    for (marker2 in target_phen[!target_phen == marker1]) {
      path_to_file <- paste0(path_to_metrics,  marker1, "_", marker2, "_spatial_stats.RData")

      list_of_metrics[[paste(marker1, marker2, sep="_")]] <- readRDS(path_to_file)
      
    }
  }
  


#
# COMPARING CELL COUNTS
#
  
  ## Basal TNBC samples
  
  #Defining rows and columns of the heat map, not of the heatmap_data matrix!
  markers <- reference_phen
  groups <- c("Total", "Basal", "Basal-IM+", "Basal-IM-", "Basal-HRD: high", "Basal-HRD: low/inter")
  
  heatmap_data <- matrix(NA, 
                         ncol = length(markers),
                         nrow = length(groups))
  
  #Defining colnames of the heatmap_data matrix
  colnames(heatmap_data) <- markers
  rownames(heatmap_data) <- groups

  for (marker in markers) {
    
    #Generating a named vector with the cell counts of interest
    counts <- unlist(sapply(list_of_metrics[[paste(marker, markers[!markers == marker][1], sep="_")]], 
                             function(my_metrics) unname(my_metrics$"Cell_Counts"[marker])))
    

    for(group in groups) {

      #Generate median cell count for each subtype 
      heatmap_data[group, marker] <- 
             
        if (group=="Total") { 
          median(na.omit(counts))
          
        } else if (group=="Basal") { 
            filter <- (annotation[names(counts), "PAM50_basal_nonBasal"] == "Basal")
            median(na.omit(counts[filter]))
            
        } else if (group=="Basal-IM+") {
        
          filter <- (annotation[names(counts), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(counts), "TNBCtype_IMpositive"] == "1")
          median(na.omit(counts[filter]))
          
        } else if (group == "Basal-IM-") {
          filter <- (annotation[names(counts), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(counts), "TNBCtype_IMpositive"] == "0")
          median(na.omit(counts[filter]))
        
        } else if (group == "Basal-HRD: high") {
          filter <- (annotation[names(counts), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(counts), "HRD.2.status"] == "high")
          median(na.omit(counts[filter]))
        
        } else if (group == "Basal-HRD: low/inter") {
          filter <- (annotation[names(counts), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(counts), "HRD.2.status"] == "low/inter")
          median(na.omit(counts[filter]))}
    }  
  }

  heatmap_data_log <- log(heatmap_data)
 
pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/basal_cell_count_with_p53.pdf")

    Heatmap(heatmap_data,
            col = colorRamp2(c(min(heatmap_data), max(heatmap_data)), c("red", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=markers,
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Cell counts")
            )

  dev.off()

change_heatmap_data <- t(apply(heatmap_data, 
                            1,
                            FUN = function(row) {
                              row / heatmap_data["Total", ]
                            }))


pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/basal_cell_count_change_with_p53.pdf")

    Heatmap(change_heatmap_data,
            col = colorRamp2(c(min(change_heatmap_data), 1, max(change_heatmap_data)), c("red", "white", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=markers,
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Change in cell counts")
            )

  dev.off()

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/basal_cell_count_change_without_p53.pdf")

    Heatmap(change_heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")],
            col = colorRamp2(c(min(change_heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")]), 1, max(change_heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")])), c("red", "white", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4"),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Change in cell counts")
            )

  dev.off()

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/basal_cell_count_without_p53.pdf")

    Heatmap(heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")],
            col = colorRamp2(c(min(heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")]), max(heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")])), c("red", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4"),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Cell counts")
            )
dev.off()

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/basal_log_cell_count_with_p53.pdf")

    Heatmap(heatmap_data_log,
            col = colorRamp2(c(min(heatmap_data_log), max(heatmap_data_log)), c("red", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=markers,
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="log(Cell counts)")
            )

  dev.off()


pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/basal_log_cell_count_without_p53.pdf")

    Heatmap(heatmap_data_log[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")],
            col = colorRamp2(c(min(heatmap_data_log[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")]), max(heatmap_data_log[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")])), c("red", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4"),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="log(Cell counts)")
            )
dev.off()


## Non-Basal TNBC samples
    
  
  # Generate a heatmap with the log(cell_counts) of the reference marker
  
  #Defining rows and columns of the heat map, not of the heatmap_data matrix!
  markers <- reference_phen
  groups <- c("Total", "Non Basal", "Non Basal-IM+", "Non Basal-IM-", "Non Basal-LAR", "Non Basal-nonLAR")
  
  heatmap_data <- matrix(NA, 
                         ncol = length(markers),
                         nrow = length(groups))
  
  #Defining colnames of the heatmap_data matrix
  colnames(heatmap_data) <- markers
  rownames(heatmap_data) <- groups

  for (marker in markers) {
    
    #Generating a named vector with the cell counts of interest
    counts <- unlist(sapply(list_of_metrics[[paste(marker, markers[!markers == marker][1], sep="_")]], 
                             function(my_metrics) unname(my_metrics$"Cell_Counts"[marker])))
    

    for(group in groups) {

      #Generate median cell count for each subtype analysed using switch()
      heatmap_data[group, marker] <- 
             
        if (group=="Total") { 
          median(na.omit(counts))
          
        } else if (group=="Non Basal") { 
            filter <- (annotation[names(counts), "PAM50_basal_nonBasal"] == "NonBasal")
            median(na.omit(counts[filter]))
            
        } else if (group=="Non Basal-IM+") {
        
          filter <- (annotation[names(counts), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(counts), "TNBCtype_IMpositive"] == "1")
          median(na.omit(counts[filter]))
          
        } else if (group == "Non Basal-IM-") {
          filter <- (annotation[names(counts), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(counts), "TNBCtype_IMpositive"] == "0")
          median(na.omit(counts[filter]))
        
        } else if (group == "Non Basal-LAR") {
          filter <- (annotation[names(counts), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(counts), "LAR_nonLAR"] == "LAR")
          median(na.omit(counts[filter]))
        
        } else if (group == "Non Basal-nonLAR") {
          filter <- (annotation[names(counts), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(counts), "LAR_nonLAR"] == "nonLAR")
          median(na.omit(counts[filter]))}
    }  
  }

  heatmap_data_log <- log(heatmap_data)

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/nonBasal_cell_count_with_p53.pdf")

    Heatmap(heatmap_data,
            col = colorRamp2(c(min(heatmap_data), max(heatmap_data)), c("red", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=markers,
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Cell counts")
            )

  dev.off()

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/nonBasal_cell_count_without_p53.pdf")

    Heatmap(heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")],
            col = colorRamp2(c(min(heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")]), max(heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")])), c("red", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4"),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Cell counts")
            )
dev.off()

change_heatmap_data <- t(cbind(apply(heatmap_data, 
                            1,
                            FUN = function(row) {
                              row / heatmap_data["Total", ]
                            })))

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/nonBasal_cell_count_change_with_p53.pdf")

    Heatmap(change_heatmap_data,
            col = colorRamp2(c(min(change_heatmap_data),1, max(change_heatmap_data)), c("red", "white", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=markers,
            row_split=c(1,2,3,3,4,4)
            )

  dev.off()

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/nonBasal_cell_count_change_without_p53.pdf")

    Heatmap(change_heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")],
            col = colorRamp2(c(min(change_heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")]), 1,  max(change_heatmap_data[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")])), c("red", "white", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4"),
            row_split=c(1,2,3,3,4,4)
            )

  dev.off()

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/nonBasal_log_cell_count_with_p53.pdf")

    Heatmap(heatmap_data_log,
            col = colorRamp2(c(min(heatmap_data_log), max(heatmap_data_log)), c("red", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=markers,
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="log(Cell counts)")
            )

  dev.off()

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/nonBasal_log_cell_count_without_p53.pdf")

    Heatmap(heatmap_data_log[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")],
            col = colorRamp2(c(min(heatmap_data_log[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")]), max(heatmap_data_log[,c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4")])), c("red", "green")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=c("CD3", "CD8", "CD20", "FOXP3", "CD68", "CD4"),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="log(Cell counts)")
            )
dev.off()


#
# COMPARING CELL DISTANCES
#

  #Determine possible combinations of markers and getting AMD matrix

  marker_permutations <- permutations(reference_phen, r=2, n=length(reference_phen))
  marker_permutations <- sapply(1:nrow(marker_permutations), function(row) paste(marker_permutations[row,], collapse="_"))

  groups <- c("Total", "Basal", "Basal-IM+", "Basal-IM-", "Basal-HRD: high", "Basal-HRD: low/inter")
  
  heatmap_data <- matrix(NA, 
                         ncol = length(marker_permutations),
                         nrow = length(groups))
  
  #Defining colnames of the heatmap_data matrix
  colnames(heatmap_data) <- marker_permutations
  rownames(heatmap_data) <- groups



for (marker in marker_permutations) {
    
    #Generating a named vector with the cell counts of interest
    AMD <- unlist(sapply(list_of_metrics[[marker]], 
                             function(my_metrics) {
                              tryCatch(
                                {my_metrics[["AMD"]][1,"Median"]},
                                error = function(e) {NA}
                              )
                            }
            ))


    for(group in groups) {

      #Generate median cell count for each subtype analysed
      heatmap_data[group, marker] <- 


        if (group=="Total") { 
          median(na.omit(AMD))
          
        } else if (group=="Basal") { 
            filter <- (annotation[names(AMD), "PAM50_basal_nonBasal"] == "Basal")
            median(na.omit(AMD[filter]))
            
        } else if (group=="Basal-IM+") {
        
          filter <- (annotation[names(AMD), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(AMD), "TNBCtype_IMpositive"] == "1")
          median(na.omit(AMD[filter]))
          
        } else if (group == "Basal-IM-") {
          filter <- (annotation[names(AMD), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(AMD), "TNBCtype_IMpositive"] == "0")
          median(na.omit(AMD[filter]))
        
        } else if (group == "Basal-HRD: high") {
          filter <- (annotation[names(AMD), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(AMD), "HRD.2.status"] == "high")
          median(na.omit(AMD[filter]))
        
        } else if (group == "Basal-HRD: low/inter") {
          filter <- (annotation[names(AMD), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(AMD), "HRD.2.status"] == "low/inter")
          median(na.omit(AMD[filter]))}
    }  
  }


pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/basal_distances_with_p53.pdf")

    Heatmap(heatmap_data,
            col = colorRamp2(c(min(heatmap_data), max(heatmap_data)), c("red", "blue")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="AMD")
            )

  dev.off()

change_heatmap_data <- t(apply(heatmap_data, 
                            1,
                            FUN = function(row) {
                              row / heatmap_data["Total", ]
                            }))

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/basal_distances.change_with_p53.pdf")

    Heatmap(change_heatmap_data,
            col = colorRamp2(c(min(change_heatmap_data), 1,  max(change_heatmap_data)), c("red", "white", "blue")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Change in AMD")
            )

  dev.off()


  #PAM 5O non Basal

   # Generate a data frame to store distances between cells

  #Determine possible combinations of markers

  marker_permutations <- permutations(reference_phen, r=2, n=length(reference_phen))
  marker_permutations <- sapply(1:nrow(marker_permutations), function(row) paste(marker_permutations[row,], collapse="_"))

  groups <- c("Total", "Non Basal", "Non Basal-IM+", "Non Basal-IM-", "Non Basal-LAR", "Non Basal-non LAR")
  
  heatmap_data <- matrix(NA, 
                         ncol = length(marker_permutations),
                         nrow = length(groups))
  
  #Defining colnames of the heatmap_data matrix
  colnames(heatmap_data) <- marker_permutations
  rownames(heatmap_data) <- groups


for (marker in marker_permutations) {
    
    #Generating a named vector with the cell counts of interest
    AMD <- unlist(sapply(list_of_metrics[[marker]], 
                             function(my_metrics) {
                              tryCatch(
                                {my_metrics[["AMD"]][1,"Median"]},
                                error = function(e) {NA}
                              )
                            }
            ))


    for(group in groups) {

      #Generate median cell count for each subtype analysed using switch()
      heatmap_data[group, marker] <- 


        if (group=="Total") { 
          median(na.omit(AMD))
          
        } else if (group=="Non Basal") { 
            filter <- (annotation[names(AMD), "PAM50_basal_nonBasal"] == "NonBasal")
            median(na.omit(AMD[filter]))
            
        } else if (group=="Non Basal-IM+") {
        
          filter <- (annotation[names(AMD), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(AMD), "TNBCtype_IMpositive"] == "1")
          median(na.omit(AMD[filter]))
          
        } else if (group == "Non Basal-IM-") {
          filter <- (annotation[names(AMD), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(AMD), "TNBCtype_IMpositive"] == "0")
          median(na.omit(AMD[filter]))
        
        } else if (group == "Non Basal-LAR") {
          filter <- (annotation[names(AMD), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(AMD), "LAR_nonLAR"] == "LAR")
          median(na.omit(AMD[filter]))
        
        } else if (group == "Non Basal-non LAR") {
          filter <- (annotation[names(AMD), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(AMD), "LAR_nonLAR"] == "nonLAR")
          median(na.omit(AMD[filter]))}
    }  
  }


pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/NonBasal_distances_with_p53.pdf")

    Heatmap(heatmap_data,
            col = colorRamp2(c(min(heatmap_data), max(heatmap_data)), c("red", "blue")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="AMD")
            )

  dev.off()

change_heatmap_data <- t(apply(heatmap_data, 
                            1,
                            FUN = function(row) {
                              row / heatmap_data["Total", ]
                            }))

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/NonBasal_distances.change_with_p53.pdf")

    Heatmap(change_heatmap_data,
            col = colorRamp2(c(min(change_heatmap_data), 1,  max(change_heatmap_data)), c("red", "white", "blue")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Change in AMD")
            )

  dev.off()


#
# COMPARING DENSITY IN THE NEIGHBORHOOD
#

## PAM50 Basal samples

# Generate matrix with median of median Density in the Neighborhood values



  #Determine possible combinations of markers and getting DIN matrix

  marker_permutations <- permutations(reference_phen, r=2, n=length(reference_phen))
  marker_permutations <- sapply(1:nrow(marker_permutations), function(row) paste(marker_permutations[row,], collapse="_"))

  groups <- c("Total", "Basal", "Basal-IM+", "Basal-IM-", "Basal-HRD: high", "Basal-HRD: low/inter")
  
  heatmap_data <- matrix(NA, 
                         ncol = length(marker_permutations),
                         nrow = length(groups))

  #Defining colnames of the heatmap_data matrix
  colnames(heatmap_data) <- marker_permutations
  rownames(heatmap_data) <- groups


  for (marker in marker_permutations) {
    
   DIN <- unlist(sapply(list_of_metrics[[marker]], 
                             function(my_metrics) {
                              tryCatch(
                                {mean(my_metrics[["DIN"]])},
                                error = function(e) {NA}
                              )
                            }
            ))


    for(group in groups) {

      #Generate median cell count for each subtype analysed
      heatmap_data[group, marker] <- 


        if (group=="Total") { 
          median(na.omit(DIN))
          
        } else if (group=="Basal") { 
            filter <- (annotation[names(DIN), "PAM50_basal_nonBasal"] == "Basal")
            median(na.omit(DIN[filter]))
            
        } else if (group=="Basal-IM+") {
        
          filter <- (annotation[names(DIN), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(DIN), "TNBCtype_IMpositive"] == "1")
          median(na.omit(DIN[filter]))
          
        } else if (group == "Basal-IM-") {
          filter <- (annotation[names(DIN), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(DIN), "TNBCtype_IMpositive"] == "0")
          median(na.omit(DIN[filter]))
        
        } else if (group == "Basal-HRD: high") {
          filter <- (annotation[names(DIN), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(DIN), "HRD.2.status"] == "high")
          median(na.omit(DIN[filter]))
        
        } else if (group == "Basal-HRD: low/inter") {
          filter <- (annotation[names(DIN), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(DIN), "HRD.2.status"] == "low/inter")
          median(na.omit(DIN[filter]))}
    }  
  }

  # Generating and saving plots

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/Basal_DIN_with_p53.pdf")

    Heatmap(heatmap_data,
            col = colorRamp2(c(min(heatmap_data), max(heatmap_data)), c("blue", "red")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="DIN (cells/pixel)")
            )
print(heatmap_data)
  dev.off()

change_heatmap_data <- t(apply(heatmap_data, 
                            1,
                            FUN = function(row) {
                              row / heatmap_data["Total", ]
                            }))



pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/Basal_ChangeDIN_with_p53.pdf")

    Heatmap(change_heatmap_data,
            col = colorRamp2(c(min(change_heatmap_data), 1,  max(change_heatmap_data)), c("red", "white", "blue")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Change in DIN)")
            )

  dev.off()


## PAM50 Non Basal samples


  marker_permutations <- permutations(reference_phen, r=2, n=length(reference_phen))
  marker_permutations <- sapply(1:nrow(marker_permutations), function(row) paste(marker_permutations[row,], collapse="_"))

  groups <- c("Total", "Non Basal", "Non Basal-IM+", "Non Basal-IM-", "Non Basal-LAR", "Non Basal-non LAR")
  
  heatmap_data <- matrix(NA, 
                         ncol = length(marker_permutations),
                         nrow = length(groups))
  
  #Defining colnames of the heatmap_data matrix
  colnames(heatmap_data) <- marker_permutations
  rownames(heatmap_data) <- groups


for (marker in marker_permutations) {
    
    #Generating a named vector with the cell counts of interest
    DIN <- unlist(sapply(list_of_metrics[[marker]], 
                             function(my_metrics) {
                              tryCatch(
                                {mean(my_metrics[["DIN"]])},
                                error = function(e) {NA}
                              )
                            }
            ))


    for(group in groups) {

      #Generate median cell count for each subtype analysed
      heatmap_data[group, marker] <- 


        if (group=="Total") { 
          median(na.omit(DIN))
          
        } else if (group=="Non Basal") { 
            filter <- (annotation[names(DIN), "PAM50_basal_nonBasal"] == "NonBasal")
            median(na.omit(DIN[filter]))
            
        } else if (group=="Non Basal-IM+") {
        
          filter <- (annotation[names(DIN), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(DIN), "TNBCtype_IMpositive"] == "1")
          median(na.omit(DIN[filter]))
          
        } else if (group == "Non Basal-IM-") {
          filter <- (annotation[names(DIN), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(DIN), "TNBCtype_IMpositive"] == "0")
          median(na.omit(DIN[filter]))
        
        } else if (group == "Non Basal-LAR") {
          filter <- (annotation[names(DIN), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(DIN), "LAR_nonLAR"] == "LAR")
          median(na.omit(DIN[filter]))
        
        } else if (group == "Non Basal-non LAR") {
          filter <- (annotation[names(DIN), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(DIN), "LAR_nonLAR"] == "nonLAR")
          median(na.omit(DIN[filter]))}
    }  
  }


pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/NonBasal_DIN_with_p53.pdf")

    Heatmap(heatmap_data,
            col = colorRamp2(c(min(heatmap_data), max(heatmap_data)), c("blue", "red")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="DIN (cells/pixel)")
            )

  dev.off()


change_heatmap_data <- t(apply(heatmap_data, 
                            1,
                            FUN = function(row) {
                              row / heatmap_data["Total", ]
                            }))



pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/NonBasal_ChangeDIN_with_p53.pdf")

    Heatmap(change_heatmap_data,
            col = colorRamp2(c(min(change_heatmap_data), 1,  max(change_heatmap_data)), c("red", "white", "blue")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Change in DIN)")
            )

  dev.off()


#
# COMPARING ATTRACTION
#


## PAM50 Basal samples

# Generate matrix with median of attraction values

  #Determine possible combinations of markers and getting DIN matrix

  marker_permutations <- permutations(reference_phen, r=2, n=length(reference_phen))
  marker_permutations <- sapply(1:nrow(marker_permutations), function(row) paste(marker_permutations[row,], collapse="_"))

  groups <- c("Total", "Basal", "Basal-IM+", "Basal-IM-", "Basal-HRD: high", "Basal-HRD: low/inter")
  
  heatmap_data <- matrix(NA, 
                         ncol = length(marker_permutations),
                         nrow = length(groups))

  #Defining colnames of the heatmap_data matrix
  colnames(heatmap_data) <- marker_permutations
  rownames(heatmap_data) <- groups

  for (marker in marker_permutations) {


   Attr <- unlist(sapply(list_of_metrics[[marker]], 
                             function(my_metrics) {
                              tryCatch(
                                {(mean(my_metrics[["DIN"]])) / (unname(my_metrics[["Cell_Counts"]][strsplit(marker, "_")[[1]][2]]) / (pi * 1500^2) )},
                                error = function(e) {NA}
                              )
                            }
            ))

    for(group in groups) {

      #Generate median cell count for each subtype analysed
      heatmap_data[group, marker] <- 


        if (group=="Total") { 
          median(na.omit(Attr))
          
        } else if (group=="Basal") { 
            filter <- (annotation[names(Attr), "PAM50_basal_nonBasal"] == "Basal")
            median(na.omit(Attr[filter]))
            
        } else if (group=="Basal-IM+") {
        
          filter <- (annotation[names(Attr), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(Attr), "TNBCtype_IMpositive"] == "1")
          median(na.omit(Attr[filter]))
          
        } else if (group == "Basal-IM-") {
          filter <- (annotation[names(Attr), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(Attr), "TNBCtype_IMpositive"] == "0")
          median(na.omit(Attr[filter]))
        
        } else if (group == "Basal-HRD: high") {
          filter <- (annotation[names(Attr), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(Attr), "HRD.2.status"] == "high")
          median(na.omit(Attr[filter]))
        
        } else if (group == "Basal-HRD: low/inter") {
          filter <- (annotation[names(Attr), "PAM50_basal_nonBasal"] == "Basal" &
                     annotation[names(Attr), "HRD.2.status"] == "low/inter")
          median(na.omit(Attr[filter]))}
    }  
  }

  # Generating and saving plots

pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/Basal_Attr_with_p53.pdf")

    Heatmap(heatmap_data,
            col = colorRamp2(c(min(heatmap_data), max(heatmap_data)), c("blue", "red")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Attr (DIN/density)")
            )

  dev.off()

change_heatmap_data <- t(apply(heatmap_data, 
                            1,
                            FUN = function(row) {
                              row / heatmap_data["Total", ]
                            }))



pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/Basal_ChangeAttr_with_p53.pdf")

    Heatmap(change_heatmap_data,
            col = colorRamp2(c(min(change_heatmap_data), 1,  max(change_heatmap_data)), c("red", "white", "blue")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Change in attraction")
            )

dev.off()


## PAM50 NonBasal samples

# Generate matrix with median of attraction values

  #Determine possible combinations of markers and getting DIN matrix

  marker_permutations <- permutations(reference_phen, r=2, n=length(reference_phen))
  marker_permutations <- sapply(1:nrow(marker_permutations), function(row) paste(marker_permutations[row,], collapse="_"))

  groups <- c("Total", "Non Basal", "Non Basal-IM+", "Non Basal-IM-", "Non Basal-LAR", "Non Basal-non LAR")
  
  heatmap_data <- matrix(NA, 
                         ncol = length(marker_permutations),
                         nrow = length(groups))

  
  #Defining colnames of the heatmap_data matrix
  colnames(heatmap_data) <- marker_permutations
  rownames(heatmap_data) <- groups

  for (marker in marker_permutations) {
    
   Attr <- unlist(sapply(list_of_metrics[[marker]], 
                             function(my_metrics) {
                              tryCatch(
                                {(mean(my_metrics[["DIN"]])) / (unname(my_metrics[["Cell_Counts"]][strsplit(marker, "_")[[1]][2]]) / (pi * 1500^2) )},
                                error = function(e) {NA}
                              )
                            }
            ))

    for(group in groups) {

      #Generate median cell count for each subtype analysed
      heatmap_data[group, marker] <- 

        if (group=="Total") {
          median(na.omit(Attr))
          
        } else if (group=="Non Basal") { 
            filter <- (annotation[names(Attr), "PAM50_basal_nonBasal"] == "NonBasal")
            median(na.omit(Attr[filter]))
            
        } else if (group=="Non Basal-IM+") {
        
          filter <- (annotation[names(Attr), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(Attr), "TNBCtype_IMpositive"] == "1")
          median(na.omit(Attr[filter]))
          
        } else if (group == "Non Basal-IM-") {
          filter <- (annotation[names(Attr), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(Attr), "TNBCtype_IMpositive"] == "0")
          median(na.omit(Attr[filter]))
        
        } else if (group == "Non Basal-LAR") {
          filter <- (annotation[names(Attr), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(Attr), "LAR_nonLAR"] == "LAR")
          median(na.omit(Attr[filter]))
        
        } else if (group == "Non Basal-non LAR") {
          filter <- (annotation[names(Attr), "PAM50_basal_nonBasal"] == "NonBasal" &
                     annotation[names(Attr), "LAR_nonLAR"] == "nonLAR")
          median(na.omit(Attr[filter]))}
    }  
  }

  # Generating and saving plots



pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/NonBasal_Attr_with_p53.pdf")

    Heatmap(heatmap_data,
            col = colorRamp2(c(min(heatmap_data), max(heatmap_data)), c("blue", "red")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Attr (DIN/density)")
            )

dev.off()

change_heatmap_data <- t(apply(heatmap_data, 
                            1,
                            FUN = function(row) {
                              row / heatmap_data["Total", ]
                            }))



pdf("/home/Illumina/Iñaki_Sasiain/immune_spatial/plots/NonBasal_ChangeAttr_with_p53.pdf")

    Heatmap(change_heatmap_data,
            col = colorRamp2(c(min(change_heatmap_data), 1,  max(change_heatmap_data)), c("red", "white", "blue")),
            cluster_rows=FALSE,
            cluster_columns=FALSE,
            row_order=groups,
            column_order=marker_permutations,
            column_names_gp = gpar(fontsize = 6),
            row_split=c(1,2,3,3,4,4),
            heatmap_legend_param=list(title="Change in attraction")
            )

dev.off()

rm(list = ls())
setwd("/home/isc/Spatial_immune_env")


# Installing packages if requitred
if (!require("survival", quietly = TRUE))
  install.packages("survivalr")

if (!require("survminer", quietly = TRUE))
  install.packages("survminer")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("vioplot", quietly = TRUE))
  install.packages("vioplot")


library(ggplot2)
library(survival)
library(survminer)
library(vioplot)

#Reading annotation file and generating survival plots
annotation <- read.csv2("data_from_suze/data/supplData_withimages.csv")
rownames(annotation) <- annotation$uid



#
# Read single marker statistics file
#

single_marker_stats <- readRDS("CD3_props_and_CIN.RData")
cell_type <- c("CD3")

#
# Plotting TILs VS molecular subgroups (6)
#

to_plot2 <- na.omit(annotation[,c("TILs", "TNBCtype_org")])
vioplot(as.numeric(to_plot2$TILs) ~ to_plot2$"TNBCtype_org",
        xlab="Molecular subtypes",
        ylab="TILs (%)")

#
# Plotting TILs VS molecular subgroups (4)
#

to_plot2 <- na.omit(annotation[,c("TILs", "TNBCtype4_n235_notPreCentered")])
vioplot(as.numeric(to_plot2$TILs) ~ to_plot2$"TNBCtype4_n235_notPreCentered",
        xlab="Molecular subtypes",
        ylab="TILs (%)")


#
# Plotting TILs VS molecular subgroups (IM)
#

to_plot2 <- na.omit(annotation[,c("TILs", "TNBCtype_IMpositive")])
vioplot(as.numeric(to_plot2$TILs) ~ to_plot2$"TNBCtype_IMpositive", 
        xlab="Immunomodulatory signature",
        ylab="TILs (%)",
        names=c("IM-", "IM+")
        )

# Cell counts in groups

p53 <- unlist(sapply(single_marker_stats, function(item) unname(item$"Cell_Counts"["p53"])))
CD3 <- unlist(sapply(single_marker_stats, function(item) unname(item$"Cell_Counts"["CD3"])))
CD20 <- unlist(sapply(single_marker_stats, function(item) unname(item$"Cell_Counts"["CD20"])))
CD8 <- unlist(sapply(single_marker_stats, function(item) unname(item$"Cell_Counts"["CD8"])))

names <- names(p53)

# Create dataframe to plot values

to_plot <- as.data.frame(cbind(
 
  "ID"=names,
  "p53"=p53,
  "CD3"=CD3,
  "CD8"=CD8,
  "CD20"=CD20,
  "Subtypes_6"=annotation[names,"TNBCtype_org"],
  "Subtypes_4"=annotation[names, "TNBCtype4_n235_notPreCentered"],
  "Subtypes_IM"=annotation[names, "TNBCtype_IMpositive"],
   
))


#
# Analyising the impact of molecular subgroups in APD
#

APD <- lapply(single_marker_stats, function(item) item$APD$"Median")
APD[sapply(APD, is.null)] <- NA
APD <- unlist(APD)

APD_names <- names(lapply(single_marker_stats, function(item) item$APD$"Median"))


subgroups <- annotation[APD_names, "TNBCtype_org"]

to_plot <- as.data.frame(cbind(
  "APD"=APD,
  "Subgroups"=subgroups,
  "Sample"=APD_names
))

to_plot$Subgroups <- as.factor(to_plot$Subgroups)
to_plot$APD <- as.numeric(to_plot$APD)


boxplot(APD~Subgroups, data=na.omit(to_plot))

summary(as.factor(annotation$TNBCtype_org))


#
# Analyising the impact of molecular subgroups in attraction
#

obs_vs_expected <- sapply(single_marker_stats, function(item) unname(item$"CIN" / item$"Cell_Proportions"[cell_type]))
obs_vs_expected[sapply(obs_vs_expected, is.null)] <- NA

attraction <- unlist(obs_vs_expected)
attraction_names <- names(attraction)

subgroups <- annotation[attraction_names, "TNBCtype_org"]

subgroups

to_plot <- as.data.frame(cbind(
  "attraction"=attraction,
  "Subgroups"=subgroups,
  "Sample"=attraction_names
))

to_plot$Subgroups <- as.factor(to_plot$Subgroups)
to_plot$attraction <- as.numeric(to_plot$attraction)


boxplot(attraction~Subgroups, data=na.omit(to_plot), ylim=c(0,10))


#
# Analysing the impact of molecular subgroups in clustering (ANNi)
#




clustering <- unlist(sapply(single_marker_stats, function(item) unname(item$"ANNI"$"ANN_index")))
names_clustering <- names(clustering)

subgroups <- annotation[names_clustering, "TNBCtype_org"]


to_plot <- as.data.frame(cbind(
  "ANN_index"=clustering,
  "Subgroups"=subgroups,
  "Sample"=names_clustering
))


to_plot$Subgroups <- as.factor(to_plot$Subgroups)
to_plot$ANN_index <- as.numeric(to_plot$ANN_index)


boxplot(ANN_index~Subgroups, data=na.omit(to_plot))

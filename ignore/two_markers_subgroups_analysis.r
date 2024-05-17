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
# Read two markers statistics file
#

two_marker_stats <- readRDS("p53_CD20_spatial_stats.RData")
target_cell_type <- "CD20"

#groups <- "TNBCtype_org"
groups <- "TNBCtype4_n235_notPreCentered"
#groups <- "TNBCtype_IMpositive"


#
# Analysis the impact of AMD
#

AMD <- unlist(lapply(two_marker_stats, function(item) item$AMD[2,"Median"]))

AMD_names <- names(AMD)


subgroups <- annotation[AMD_names, groups]

to_plot <- as.data.frame(cbind(
  "AMD"=AMD,
  "Subgroups"=subgroups,
  "Sample"=AMD_names
))

to_plot$Subgroups <- as.factor(to_plot$Subgroups)
to_plot$AMD <- as.numeric(to_plot$AMD)


vioplot(AMD~Subgroups, data=na.omit(to_plot)
 #       , names=c("IM-", "IM+")
        )



#
# Analysis the impact of CIN
#


CIN <- unlist(lapply(two_marker_stats, function(item) item$CIN / item$"Cell_Proportions"[target_cell_type]))

CIN_names <- names(CIN)


subgroups <- annotation[CIN_names, groups]

to_plot <- as.data.frame(cbind(
  "CIN"=CIN,
  "Subgroups"=subgroups,
  "Sample"=CIN_names
))

to_plot$Subgroups <- as.factor(to_plot$Subgroups)
to_plot$CIN <- as.numeric(to_plot$CIN)


vioplot(CIN~Subgroups, data=na.omit(to_plot)
#      , names=c("IM-", "IM+")
        )


#
# Analysis the impact of attraction
#


CIN <- unlist(lapply(two_marker_stats, function(item) item$CIN / unname(item$"Cell_Proportions"[target_cell_type])))

CIN_names <- names(CIN)



subgroups <- annotation[CIN_names, groups]

to_plot <- as.data.frame(cbind(
  "CIN"=CIN,
  "Subgroups"=subgroups,
  "Sample"=CIN_names
))

to_plot$Subgroups <- as.factor(to_plot$Subgroups)
to_plot$CIN <- as.numeric(to_plot$CIN)


vioplot(CIN~Subgroups, data=na.omit(to_plot)
             , names=c("IM-", "IM+")
        #,ylim=c(0,4)
        ,ylab = "Obs./Exp. CIN"
)




#
# Analysis the impact of AUC of cross Kcross function
#


Kcross <- lapply(two_marker_stats, function(item) item$AUC)
AUC_Kcross <- unlist(lapply(Kcross, function(sample) AUC_of_cross_function(sample)))
names_AUC <- names(AUC_Kcross)


subgroups <- annotation[names_AUC, groups]

to_plot <- as.data.frame(cbind(
  "AUC"=AUC_Kcross,
  "Subgroups"=subgroups,
  "Sample"=names_AUC
))

to_plot$Subgroups <- as.factor(to_plot$Subgroups)
to_plot$AUC <- as.numeric(to_plot$AUC)


vioplot(AUC~Subgroups, data=na.omit(to_plot))


rm(list = ls())
setwd("/home/isc/Spatial_immune_env")

annotation <- read.csv2("data_from_suze/data/supplData_withimages.csv")
rownames(annotation) <- annotation$uid

# Installing packages if requitred
if (!require("survival", quietly = TRUE))
  install.packages("survivalr")

if (!require("survminer", quietly = TRUE))
  install.packages("survminer")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")


library(ggplot2)
library(survival)
library(survminer)

#
# TILs survival plot
#

# Transform TILs into a numeric variable
annotation$TILs <- as.numeric(annotation$TILs)

# Group TILs values in group
annotation$grouped_TILs <- cut(annotation$TILs, 
                               breaks = c(0, 5, 15, 30, 60, 100),
                               labels = c("0-5", "5-15", "15-30", "30-60", "60-100"),
                               include.lowest = TRUE)


#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(as.numeric(time), event) ~ grouped_TILs, data = annotation)
fit

png(filename = "TILs_survival_plot.png")

ggsurvplot(fit, 
           data=annotation,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="TILs (%)",
           legend.labs=c("0-5", "5-15", "15-30", "30-60", "60-100"))

dev.off()

#
# Moleculat subgroupssurvival plot
#

## IM+ and IM-

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(as.numeric(time), event) ~ TNBCtype_IMpositive, data = annotation)

ggsurvplot(fit, 
           data=annotation,
           size=2,
           conf.int = FALSE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="IM signature", 
           legend.labs=c("IM-", "IM+"))


## Gene expression profiles

annotation$TNBCtype4_n235_notPreCentered

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(as.numeric(time), event) ~ TNBCtype4_n235_notPreCentered, data = annotation)

ggsurvplot(fit, 
           data=annotation,
           size=2,
           conf.int = FALSE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="Molecular subgroup",
           legend.labs=c("BL1", "BL2", "LAR", "M")
           )

#
# Generate survival plot for the attraction between p53 cells. 2 CAT
#

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(as.numeric(time), event) ~ TNBCtype_org, data = annotation)

ggsurvplot(fit, 
           data=annotation,
           size=2,
           conf.int = FALSE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="Molecular subgroup",
           legend.labs=c("BL1", "BL2", "IM", "LAR", "M", "USL", "UNS"))

p53_props_and_CIN <- readRDS("p53_props_and_CIN.RData")

obs_vs_expected <- sapply(p53_props_and_CIN, function(item) item$"CIN" / item$"Cell_Proportions"["p53"])

#Remove samples over 1.45. They are normally IHC artefacts
obs_vs_expected <- obs_vs_expected[obs_vs_expected <= 1.45]



names_to_plot <- unlist(strsplit(names(obs_vs_expected), ".p53"))

analyse_p53_overrepresentation <- as.data.frame(cbind(
  
  "name"=names_to_plot,
  "time"=as.numeric(annotation[names_to_plot, "time"]),
  "event"= annotation[names_to_plot, "event"],
  "p53_overrepresentation"= cut(obs_vs_expected,
                              breaks = c(1, 1.1, 1.5),
                              labels = c("1-1.1", "1.1-1.5"),
                              include_lowest = TRUE)
  
  
))

analyse_p53_overrepresentation$time <- as.numeric(analyse_p53_overrepresentation$time)
analyse_p53_overrepresentation$event <- as.numeric(analyse_p53_overrepresentation$event)
summary(analyse_p53_overrepresentation)

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ p53_overrepresentation, data = analyse_p53_overrepresentation)

png(filename = "p53_attraction_survival_plot_2cat.png")

ggsurvplot(fit, 
           data=analyse_p53_overrepresentation,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="Obs./Ex. probability",
           legend.labs=c("1-1.1", "1.1-1.5"))
dev.off()


#
# Generate survival plot for the attraction between p53 cells. 5 CAT
#

p53_props_and_CIN <- readRDS("p53_props_and_CIN.RData")

obs_vs_expected <- sapply(p53_props_and_CIN, function(item) item$"CIN" / item$"Cell_Proportions"["p53"])

#Remove samples over 1.45. They are normally IHC artefacts
obs_vs_expected <- obs_vs_expected[obs_vs_expected <= 1.45]



names_to_plot <- unlist(strsplit(names(obs_vs_expected), ".p53"))

analyse_p53_overrepresentation <- as.data.frame(cbind(
  
  "name"=names_to_plot,
  "time"=as.numeric(annotation[names_to_plot, "time"]),
  "event"= annotation[names_to_plot, "event"],
  "p53_overrepresentation"= cut(obs_vs_expected,
                                breaks = c(1, 1.1, 1.2, 1.3, 1.4, 1.5),
                                labels = c("1-1.1", "1.1-1.2", "1.2-1.3", "1.3-1.4", "1.4-1.5"),
                                include_lowest = TRUE)
  
  
))

analyse_p53_overrepresentation$time <- as.numeric(analyse_p53_overrepresentation$time)
analyse_p53_overrepresentation$event <- as.numeric(analyse_p53_overrepresentation$event)
summary(analyse_p53_overrepresentation)

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ p53_overrepresentation, data = analyse_p53_overrepresentation)

png(filename = "p53_attraction_survival_plot_5cat.png")

ggsurvplot(fit, 
           data=analyse_p53_overrepresentation,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="Obs./Ex. probability",
           legend.labs=c("1-1.1", "1.1-1.2", "1.2-1.3", "1.3-1.4", "1.4-1.5"))
dev.off()



#
# Generate survival plot for the clustering between p53 cells
#


p53_clusters_info <- readRDS("clustering_metrics_p53.RData")

analyse_p53_clustering <- as.data.frame(cbind(
  
  "name"=rownames(p53_clusters_info),
  "time"=as.numeric(annotation[p53_clusters_info$uid, "time"]),
  "event"= annotation[p53_clusters_info$uid, "event"],
  "p53_clustered"= p53_clusters_info$Clustered
  
))

analyse_p53_clustering$time <- as.numeric(analyse_p53_clustering$time)
analyse_p53_clustering$event <- as.numeric(analyse_p53_clustering$event)

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ p53_clustered, data = analyse_p53_clustering)

png(filename = "p53_clustering_survival_plot.png")

ggsurvplot(fit, 
           data=analyse_p53_clustering,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",           
           legend.title="p53 clusters",
           legend.labs=c("FALSE", "TRUE"))

dev.off()

#
# Generate survival plot for the attraction between CD3 cells. 2 CAT
#

CD3_props_and_CIN <- readRDS("CD3_props_and_CIN.RData")

obs_vs_expected <- sapply(CD3_props_and_CIN, function(item) item$"CIN" / item$"Cell_Proportions"["CD3"])

#Removc("1-1.1", "1.1-1.3", "1.3-1.5")e samples over 1.45. They are normally IHC artefacts
#obs_vs_expected <- obs_vs_expected[obs_vs_expected <= 1.45]



names_to_plot <- unlist(strsplit(names(obs_vs_expected), ".CD3"))

analyse_CD3_overrepresentation <- as.data.frame(cbind(
  
  "name"=names_to_plot,
  "time"=as.numeric(annotation[names_to_plot, "time"]),
  "event"= annotation[names_to_plot, "event"],
  "CD3_overrepresentation"= cut(obs_vs_expected,
                                breaks = c(1, 1.5, 5),
                                labels = c("1-1.5", "1.5-5"),
                                include_lowest = TRUE)
  
  
))

analyse_CD3_overrepresentation$time <- as.numeric(analyse_CD3_overrepresentation$time)
analyse_CD3_overrepresentation$event <- as.numeric(analyse_CD3_overrepresentation$event)
summary(analyse_CD3_overrepresentation)

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ CD3_overrepresentation, data = analyse_CD3_overrepresentation)

png(filename = "CD3_attraction_survival_plot_2cat.png")

ggsurvplot(fit, 
           data=analyse_CD3_overrepresentation,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="Obs./Ex. probability",
           legend.labs=c("1-5", "1.5-5"))
dev.off()


#
# Generate survival plot for the attraction between CD3 cells. 5 CAT
#

CD3_props_and_CIN <- readRDS("CD3_props_and_CIN.RData")

obs_vs_expected <- sapply(CD3_props_and_CIN, function(item) item$"CIN" / item$"Cell_Proportions"["CD3"])

#Removc("1-1.1", "1.1-1.3", "1.3-1.5")e samples over 1.45. They are normally IHC artefacts
#obs_vs_expected <- obs_vs_expected[obs_vs_expected <= 1.45]



names_to_plot <- unlist(strsplit(names(obs_vs_expected), ".CD3"))

analyse_CD3_overrepresentation <- as.data.frame(cbind(
  
  "name"=names_to_plot,
  "time"=as.numeric(annotation[names_to_plot, "time"]),
  "event"= annotation[names_to_plot, "event"],
  "CD3_overrepresentation"= cut(obs_vs_expected,
                                breaks = c(1, 1.35, 1.5, 2, 3, 5),
                                labels = c("1-1.35", "1.35-1.5", "1.5-2", "2-3", "3-5"),
                                include_lowest = TRUE)
  
  
))



analyse_CD3_overrepresentation$time <- as.numeric(analyse_CD3_overrepresentation$time)
analyse_CD3_overrepresentation$event <- as.numeric(analyse_CD3_overrepresentation$event)
summary(analyse_CD3_overrepresentation)

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ CD3_overrepresentation, data = analyse_CD3_overrepresentation)

png(filename = "CD3_attraction_survival_plot_5cat.png")

ggsurvplot(fit, 
           data=analyse_CD3_overrepresentation,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="Obs./Ex. probability",
           legend.labs=c("1-1.35", "1.35-1.5", "1.5-2", "2-3", "3-5"))
dev.off()



#
# Generate survival plot for the clustering between CD3 cells
#


CD3_clusters_info <- readRDS("clustering_metrics_CD3.RData")

analyse_CD3_clustering <- as.data.frame(cbind(
  
  "name"=rownames(CD3_clusters_info),
  "time"=as.numeric(annotation[CD3_clusters_info$uid, "time"]),
  "event"= annotation[CD3_clusters_info$uid, "event"],
  "CD3_clustered"= CD3_clusters_info$Clustered
  
))

analyse_CD3_clustering$time <- as.numeric(analyse_CD3_clustering$time)
analyse_CD3_clustering$event <- as.numeric(analyse_CD3_clustering$event)
summary(as.factor(CD3_clusters_info))

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ CD3_clustered, data = analyse_CD3_clustering)

png(filename = "CD3_clustering_survival_plot.png")

ggsurvplot(fit, 
           data=analyse_CD3_clustering,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",           
           legend.title="CD3 clusters",
           legend.labs=c("FALSE", "TRUE")
           )

dev.off()


#
# Generate survival plot for the attraction between CD8 cells. 2 CAT
#

CD8_props_and_CIN <- readRDS("CD8_props_and_CIN.RData")

obs_vs_expected <- sapply(CD8_props_and_CIN, function(item) item$"CIN" / item$"Cell_Proportions"["CD8"])

#Removc("1-1.1", "1.1-1.3", "1.3-1.5")e samples over 1.45. They are normally IHC artefacts
#obs_vs_expected <- obs_vs_expected[obs_vs_expected <= 1.45]



names_to_plot <- unlist(strsplit(names(obs_vs_expected), ".CD8"))

analyse_CD8_overrepresentation <- as.data.frame(cbind(
  
  "name"=names_to_plot,
  "time"=as.numeric(annotation[names_to_plot, "time"]),
  "event"= annotation[names_to_plot, "event"],
  "CD8_overrepresentation"= cut(obs_vs_expected,
                                breaks = c(1, 1.5, 5),
                                labels = c("1-1.5", "1.5-5"),
                                include_lowest = TRUE)
  
  
))

analyse_CD8_overrepresentation$time <- as.numeric(analyse_CD8_overrepresentation$time)
analyse_CD8_overrepresentation$event <- as.numeric(analyse_CD8_overrepresentation$event)
summary(analyse_CD8_overrepresentation)

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ CD8_overrepresentation, data = analyse_CD8_overrepresentation)

png(filename = "CD8_attraction_survival_plot_2cat.png")

ggsurvplot(fit, 
           data=analyse_CD8_overrepresentation,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="Obs./Ex. probability",
           legend.labs=c("1-5", "1.5-5"))
dev.off()


#
# Generate survival plot for the attraction between CD8 cells. 5 CAT
#

CD8_props_and_CIN <- readRDS("CD8_props_and_CIN.RData")

obs_vs_expected <- sapply(CD8_props_and_CIN, function(item) item$"CIN" / item$"Cell_Proportions"["CD8"])

#Removc("1-1.1", "1.1-1.3", "1.3-1.5")e samples over 1.45. They are normally IHC artefacts
#obs_vs_expected <- obs_vs_expected[obs_vs_expected <= 1.45]



names_to_plot <- unlist(strsplit(names(obs_vs_expected), ".CD8"))

analyse_CD8_overrepresentation <- as.data.frame(cbind(
  
  "name"=names_to_plot,
  "time"=as.numeric(annotation[names_to_plot, "time"]),
  "event"= annotation[names_to_plot, "event"],
  "CD8_overrepresentation"= cut(obs_vs_expected,
                                breaks = c(1, 1.35, 1.5, 2, 3, 5),
                                labels = c("1-1.35", "1.35-1.5", "1.5-2", "2-3", "3-5"),
                                include_lowest = TRUE)
  
  
))



analyse_CD8_overrepresentation$time <- as.numeric(analyse_CD8_overrepresentation$time)
analyse_CD8_overrepresentation$event <- as.numeric(analyse_CD8_overrepresentation$event)
summary(analyse_CD8_overrepresentation)

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ CD8_overrepresentation, data = analyse_CD8_overrepresentation)

png(filename = "CD8_attraction_survival_plot_5cat.png")

ggsurvplot(fit, 
           data=analyse_CD8_overrepresentation,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="Obs./Ex. probability",
           legend.labs=c("1-1.35", "1.35-1.5", "1.5-2", "2-3", "3-5"))
dev.off()



#
# Generate survival plot for the clustering between CD8 cells
#


CD8_clusters_info <- readRDS("clustering_metrics_CD8.RData")

analyse_CD8_clustering <- as.data.frame(cbind(
  
  "name"=rownames(CD8_clusters_info),
  "time"=as.numeric(annotation[CD8_clusters_info$uid, "time"]),
  "event"= annotation[CD8_clusters_info$uid, "event"],
  "CD8_clustered"= CD8_clusters_info$Clustered
  
))

analyse_CD8_clustering$time <- as.numeric(analyse_CD8_clustering$time)
analyse_CD8_clustering$event <- as.numeric(analyse_CD8_clustering$event)
summary(as.factor(CD8_clusters_info))

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ CD8_clustered, data = analyse_CD8_clustering)

png(filename = "CD8_clustering_survival_plot.png")

ggsurvplot(fit, 
           data=analyse_CD8_clustering,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",           
           legend.title="CD8 clusters",
           legend.labs=c("FALSE", "TRUE")
)

dev.off()


#
# Generate survival plot for the attraction between CD20 cells. 2 CAT
#

CD20_props_and_CIN <- readRDS("CD20_props_and_CIN.RData")

obs_vs_expected <- sapply(CD20_props_and_CIN, function(item) item$"CIN" / item$"Cell_Proportions"["CD20"])

#Removc("1-1.1", "1.1-1.3", "1.3-1.5")e samples over 1.45. They are normally IHC artefacts
#obs_vs_expected <- obs_vs_expected[obs_vs_expected <= 1.45]



names_to_plot <- unlist(strsplit(names(obs_vs_expected), ".CD20"))

analyse_CD20_overrepresentation <- as.data.frame(cbind(
  
  "name"=names_to_plot,
  "time"=as.numeric(annotation[names_to_plot, "time"]),
  "event"= annotation[names_to_plot, "event"],
  "CD20_overrepresentation"= cut(obs_vs_expected,
                                breaks = c(1, 1.5, 5),
                                labels = c("1-1.5", "1.5-5"),
                                include_lowest = TRUE)
  
  
))

analyse_CD20_overrepresentation$time <- as.numeric(analyse_CD20_overrepresentation$time)
analyse_CD20_overrepresentation$event <- as.numeric(analyse_CD20_overrepresentation$event)
summary(analyse_CD20_overrepresentation)

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ CD20_overrepresentation, data = analyse_CD20_overrepresentation)

png(filename = "CD20_attraction_survival_plot_2cat.png")

ggsurvplot(fit, 
           data=analyse_CD20_overrepresentation,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="Obs./Ex. probability",
           legend.labs=c("1-5", "1.5-5"))
dev.off()


#
# Generate survival plot for the attraction between CD20 cells. 5 CAT
#

CD20_props_and_CIN <- readRDS("CD20_props_and_CIN.RData")

obs_vs_expected <- sapply(CD20_props_and_CIN, function(item) item$"CIN" / item$"Cell_Proportions"["CD20"])
obs_vs_expected <- sapply(CD20_props_and_CIN, function(item) item$"CIN" / item$"Cell_Proportions"["CD20"])


#Removc("1-1.1", "1.1-1.3", "1.3-1.5")e samples over 1.45. They are normally IHC artefacts
#obs_vs_expected <- obs_vs_expected[obs_vs_expected <= 1.45]



names_to_plot <- unlist(strsplit(names(obs_vs_expected), ".CD20"))

analyse_CD20_overrepresentation <- as.data.frame(cbind(
  
  "name"=names_to_plot,
  "time"=as.numeric(annotation[names_to_plot, "time"]),
  "event"= annotation[names_to_plot, "event"],
  "CD20_overrepresentation"= cut(obs_vs_expected,
                                breaks = c(1, 1.35, 1.5, 2, 3, 5),
                                labels = c("1-1.35", "1.35-1.5", "1.5-2", "2-3", "3-5"),
                                include_lowest = TRUE)
  
  
))



analyse_CD20_overrepresentation$time <- as.numeric(analyse_CD20_overrepresentation$time)
analyse_CD20_overrepresentation$event <- as.numeric(analyse_CD20_overrepresentation$event)
summary(analyse_CD20_overrepresentation)

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ CD20_overrepresentation, data = analyse_CD20_overrepresentation)

png(filename = "CD20_attraction_survival_plot_5cat.png")

ggsurvplot(fit, 
           data=analyse_CD20_overrepresentation,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="Obs./Ex. probability",
           legend.labs=c("1-1.35", "1.35-1.5", "1.5-2", "2-3", "3-5"))
dev.off()



#
# Generate survival plot for the clustering between CD20 cells
#


CD20_clusters_info <- readRDS("clustering_metrics_CD20.RData")

analyse_CD20_clustering <- as.data.frame(cbind(
  
  "name"=rownames(CD20_clusters_info),
  "time"=as.numeric(annotation[CD20_clusters_info$uid, "time"]),
  "event"= annotation[CD20_clusters_info$uid, "event"],
  "CD20_clustered"= CD20_clusters_info$Clustered
  
))

analyse_CD20_clustering$time <- as.numeric(analyse_CD20_clustering$time)
analyse_CD20_clustering$event <- as.numeric(analyse_CD20_clustering$event)
summary(as.factor(CD20_clusters_info))

#Generating a Kapland-Meier survival plot
fit <- survfit(Surv(time, event) ~ CD20_clustered, data = analyse_CD20_clustering)

png(filename = "CD20_clustering_survival_plot.png")

ggsurvplot(fit, 
           data=analyse_CD20_clustering,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",           
           legend.title="CD20 clusters",
           legend.labs=c("FALSE", "TRUE")
           )

dev.off()

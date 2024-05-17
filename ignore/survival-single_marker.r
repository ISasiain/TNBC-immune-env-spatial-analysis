rm(list = ls())
setwd("/home/isc/Spatial_immune_env")


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

#Reading annotation file and generating survival plots
annotation <- read.csv2("data_from_suze/data/supplData_withimages.csv")
rownames(annotation) <- annotation$uid



#
# Read single marker statistics file
#

single_marker_stats <- readRDS("p53_props_and_CIN.RData")
cell_type <- c("p53")

#
# Plot the effect of APD in survival
#

APD_Mean <- unlist(lapply(single_marker_stats, function(item) item$APD$Mean))
#names_APD <- names(unlist(lapply(single_marker_stats, function(item) item$APD$Mean)))

#plot(sort(APD_Mean))

#analyse_APD <- as.data.frame(cbind(
  
#  "name"=names_APD,
#  "time"=as.numeric(annotation[names_APD, "time"]),
#  "event"= annotation[names_APD, "event"],
#  "Mean_APD"= cut(APD_Mean,
#                  breaks = c(0, 1100, 1250, 1300, 2500),
#                  labels = c("0-1100", "1100-1200", "1200-1300", "1300-2500"),
#                  include_lowest = TRUE)
  
#))

#analyse_APD$time <- as.numeric(analyse_APD$time)
#analyse_APD$event <- as.numeric(analyse_APD$event)
#analyse_APD$Mean_APD <- as.factor(analyse_APD$Mean_APD)


#fit <- survfit(Surv(time, event) ~ Mean_APD, data = analyse_APD)

#ggsurvplot(fit, 
#           data=analyse_APD,
#           size=2,
#           conf.int = TRUE,
#           censor.shape="|",
#           censor.size=4,
#           risk.table = TRUE,
#           xlab="Time (years)",
#           ylab="Survival probability (%)")


#
# Plot the effect of CIN_200 in survival
#

CIN <- unlist(lapply(single_marker_stats, function(item) item$CIN))
names_CIN <- names(unlist(lapply(single_marker_stats, function(item) item$CIN)))

hist(sort(CIN))
sort(CIN)

  analyse_CIN <- as.data.frame(cbind(
  
  "name"=names_CIN,
  "time"=as.numeric(annotation[names_CIN, "time"]),
  "event"= annotation[names_CIN, "event"],
  "CIN_200"= cut(CIN,
                  breaks = c(0, 50, 65, 80, 90, 100),
                  labels = c("0-55", "50-65", "65-80", "80-90", "90-100"),
                  include_lowest = TRUE)
  
))

analyse_CIN$time <- as.numeric(analyse_CIN$time)
analyse_CIN$event <- as.numeric(analyse_CIN$event)
analyse_CIN$CIN_200 <- as.factor(analyse_CIN$CIN_200)


fit <- survfit(Surv(time, event) ~ CIN, data = analyse_CIN)

ggsurvplot(fit, 
           data=analyse_CIN,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)", 
           label.title="CIN (%)",
           legend.labs=c("0-55", "50-65", "65-80", "80-90", "90-100")
           )

#
# Plot the effect of attraction in survival
#

obs_vs_expected <- sapply(single_marker_stats, function(item) unname(item$"CIN" / item$"Cell_Proportions"[cell_type]))

#Removc("1-1.1", "1.1-1.3", "1.3-1.5")e samples over 1.45. They are normally IHC artefacts
#obs_vs_expected <- obs_vs_expected[obs_vs_expected <= 1.45]


attraction <- unlist(unname(obs_vs_expected))
names_attraction <- names(obs_vs_expected)

plot(sort(attraction))

analyse_attraction<- as.data.frame(cbind(
  
  "name"=names_attraction,
  "time"=as.numeric(annotation[names_attraction, "time"]),
  "event"= annotation[names_attraction, "event"],
  "attraction"= cut(attraction,
                                breaks = c(1, 1.35, 1.5, 2, 3, 5),
                                labels = c("1-1.35", "1.35-1.5", "1.5-2", "2-3", "3-5"),
                                include_lowest = TRUE)
  
  
))



analyse_attraction$time <- as.numeric(analyse_attraction$time)
analyse_attraction$event <- as.numeric(analyse_attraction$event)


fit <- survfit(Surv(time, event) ~ attraction, data = analyse_attraction)

ggsurvplot(fit, 
           data=analyse_attraction,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)")


#
# Plot the effect of ANN index (clustering) in survival
#


clustering <- unname(unlist(sapply(single_marker_stats, function(item) unname(item$"ANNI"$"pattern"))))
names_clustering <- names(unlist(sapply(single_marker_stats, function(item) unname(item$"ANNI"$"pattern"))))


analyse_clustering<- as.data.frame(cbind(
  
  "name"=names_clustering,
  "time"=as.numeric(annotation[names_clustering, "time"]),
  "event"= annotation[names_clustering, "event"],
  "clusters"= clustering)
  
  
)


analyse_clustering$time <- as.numeric(analyse_clustering$time)
analyse_clustering$event <- as.numeric(analyse_clustering$event)


fit <- survfit(Surv(time, event) ~ clusters, data = analyse_clustering)


ggsurvplot(fit, 
           data=analyse_clustering,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)")

ANNi <- unlist(sapply(single_marker_stats, function(item) unname(item$"ANNI"$"ANN_index")))
APD <- unlist(sapply(single_marker_stats, function(item) unname(item$"APD"$"Mean")))
CIN <- unlist(sapply(single_marker_stats, function(item) unname(item$"CIN")))
attraction <- unlist(sapply(single_marker_stats, function(item) unname(item$"CIN"))) / unlist(sapply(single_marker_stats, function(item) unname(item$"Cell_Proportions"[cell_type])))



#
# Cell proportions and survival
#


proportions <- unname(unlist(sapply(single_marker_stats, function(item) unname(item$"Cell_Proportions"[cell_type]))))
names_proportions <- names(unlist(sapply(single_marker_stats, function(item) unname(item$"Cell_Proportions"[cell_type]))))


analyse_proportions<- as.data.frame(cbind(
  
  "name"=names_proportions,
  "time"=as.numeric(annotation[names_proportions, "time"]),
  "event"= annotation[names_proportions, "event"],
  "proportions"= cut(proportions,
                     breaks = c(0, 40, 60, 80, 100),
                     labels = c("0-40", "40-60", "60-80", "80-100"),
                     include_lowest = TRUE))
  
  
)


analyse_proportions$time <- as.numeric(analyse_proportions$time)
analyse_proportions$event <- as.numeric(analyse_proportions$event)


fit <- survfit(Surv(time, event) ~ proportions, data = analyse_proportions)


ggsurvplot(fit, 
           data=analyse_proportions,
           size=2,
           conf.int = TRUE,
           censor.shape="|",
           censor.size=4,
           risk.table = TRUE,
           xlab="Time (years)",
           ylab="Survival probability (%)",
           legend.title="Cell proportions (%)",
           legend.labs=c("0-40", "40-60", "60-80", "80-100"))

ANNi <- unlist(sapply(single_marker_stats, function(item) unname(item$"ANNI"$"ANN_index")))
APD <- unlist(sapply(single_marker_stats, function(item) unname(item$"APD"$"Mean")))
CIN <- unlist(sapply(single_marker_stats, function(item) unname(item$"CIN")))
attraction <- unlist(sapply(single_marker_stats, function(item) unname(item$"CIN"))) / unlist(sapply(single_marker_stats, function(item) unname(item$"Cell_Proportions"[cell_type])))

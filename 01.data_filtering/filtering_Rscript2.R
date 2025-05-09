library(tidyverse)
library(cowplot)
library(stringi)
library(ggpubr)

QD <- read.delim("QD_values.txt") 
colnames(QD)[1]<-"QD"
QUAL <- read.delim("QUAL_values.txt")  
colnames(QUAL)[1]<-"QUAL" 
SOR  <- read.delim("SOR_values.txt")  
colnames(SOR)[1]<-"SOR" 
FS <- read.delim("FS_values.txt")  
colnames(FS)[1]<-"FS" 
MQ <- read.delim("MQ_values.txt")  
colnames(MQ)[1]<-"MQ" 
MQRankSum <- read.delim("MQRankSum_values.txt")
colnames(MQRankSum)[1]<-"MQRankSum" 
RPRS <- read.delim("RPRS_values.txt")  
colnames(RPRS)[1]<-"RPRS"  

QD_filtered <- read.delim("QD_values_filter.txt") 
colnames(QD_filtered)[1]<-"QD" 
QUAL_filtered  <- read.delim("QUAL_values_filter.txt") 
colnames(QUAL_filtered)[1]<-"QUAL" 
SOR_filtered   <- read.delim("SOR_values_filter.txt") 
colnames(SOR_filtered)[1]<-"SOR"
FS_filtered  <- read.delim("FS_values_filter.txt") 
colnames(FS_filtered)[1]<-"FS"
MQ_filtered  <- read.delim("MQ_values_filter.txt")
colnames(MQ_filtered)[1]<-"MQ"
MQRankSum_filtered  <- read.delim("MQRankSum_values_filter.txt") 
colnames(MQRankSum_filtered)[1]<-"MQRankSum"
RPRS_filtered  <- read.delim("RPRS_values_filter.txt") 
colnames(RPRS_filtered)[1]<-"RPRS"

q1 <- ggplot(QD) + geom_histogram(aes(as.numeric(x=QD), y=..density..), binwidth = 1, color="black", fill = "#455A64") + 
  geom_histogram(data=QD_filtered, aes(x=as.numeric(QD), y=..density..), binwidth = 1, color="black", fill = "red") +
  xlab("Quailty Depth") + ylab("Density") + 
  theme_classic() + 
  geom_vline(xintercept =2, linetype ="dashed", color = "red", size = 1.5)

q2 <- ggplot(QUAL) + geom_histogram(aes(x=as.numeric(QUAL), y=..density..), binwidth = 10, color="black", fill = "#455A64") +
  geom_histogram(data=QUAL_filtered, aes(x=as.numeric(QUAL), y=..density..), binwidth = 10, color="red", fill = "red") +
  xlab("Quailty") + ylab("Density") + 
  theme_classic() + xlim(0,5000) +
  geom_vline(xintercept =30, linetype ="dashed", color = "red", size = 1.5)

q3 <- ggplot(SOR) + geom_histogram(aes(x=as.numeric(SOR), y=..density..), binwidth = 1, color="black", fill = "#455A64") + 
  geom_histogram(data=SOR_filtered, aes(x=as.numeric(SOR), y=..density..), binwidth = 1, color="black", fill = "red") +
  xlab("Strand Odds Ratio") + ylab("Density") + 
  theme_classic() + 
  geom_vline(xintercept =3, linetype ="dashed", color = "red", size = 1.5)

q4 <- ggplot(FS) + geom_histogram(aes(x=as.numeric(FS), y=..density..), binwidth = 2, color="black", fill = "#455A64") +
  geom_histogram(data=FS_filtered, aes(x=as.numeric(FS), y=..density..), binwidth = 2, color="black", fill = "red") +
  xlab("Disher Strand Bias") + ylab("Density") + xlim(0,100) +
  theme_classic() + 
  geom_vline(xintercept =60, linetype ="dashed", color = "red", size = 1.5)

q5 <- ggplot(MQRankSum) + geom_histogram(aes(x=as.numeric(MQRankSum), y=..density..), binwidth = 2, color="black", fill = "#455A64") + 
  geom_histogram(data=MQRankSum_filtered, aes(x=as.numeric(MQRankSum), y=..density..), binwidth = 2, color="black", fill = "red") + 
  xlab("Mapping Quailty Rank Sum") + ylab("Density") + 
  theme_classic() + 
  geom_vline(xintercept =-12.5, linetype ="dashed", color = "red", size = 1.5)

q6 <- ggplot(RPRS) + geom_histogram(aes(x=as.numeric(RPRS), y=..density..), binwidth = 2, color="black", fill = "#455A64") + 
  geom_histogram(data=RPRS_filtered, aes(x=as.numeric(RPRS), y=..density..), binwidth = 2, color="black", fill = "red") +
  xlab("Rank Sum Test of Alt vs. Ref read position") + ylab("Density") + 
  theme_classic() + 
  geom_vline(xintercept =-8, linetype ="dashed", color = "red", size = 1.5) 

plot_grid(q1,q2,q3,q4,q5,q6, nrow=2, labels =c("QD","QUAL","SOR","FS","MQRankSum","ReadPosRankSum"))
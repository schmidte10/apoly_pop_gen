library("pcadapt")
library("ggplot2")

apoly_bed <- "data/apoly_popualtions.bed"
apoly_pcadapt <-  read.pcadapt(apoly_bed, type = "bed")

#apoly_populations_kplot <- pcadapt(input = apoly_pcadapt, K=20)
#pdf("analysis/pcadapt/apoly_populations_kplot.pdf")
#plot(apoly_populations_kplot, option = "screeplot")
#dev.off() 

poplist.names <- read.delim("data/poppop.txt", header=FALSE)
apoly_populations_pca <- pcadapt(input =apoly_pcadapt, K=3)  
write.table(apoly_populations_pca$scores, "analysis/pcadapt/pcapcadapt_scores.txt")
summary(apoly_populations_pca)

pdf("analysis/pcadapt/pcadapt_apoly_projection1v2.pdf")
plot(apoly_populations_pca, option = "scores", i=1, j=2, pop = poplist.names$V3)
dev.off() 

pdf("analysis/pcadapt/pcadapt_apoly_projection1v3.pdf")
plot(apoly_populations_pca, option = "scores", i=1, j=3, pop = poplist.names$V3)
dev.off() 

pdf("analysis/pcadapt/pcadapt_apoly_projection2v3.pdf")
plot(apoly_populations_pca, option = "scores", i=2, j=3, pop = poplist.names$V3)
dev.off() 

pdf("analysis/pcadapt/pcadapt_apoly_manhattan.pdf")
plot(apoly_populations_pca, option= "manhattan") 
dev.off() 

pdf("analysis/pcadapt/pcadapt_apoly_qqplot.pdf") 
plot(apoly_populations_pca, option ="qqplot") 
dev.off()

apoly_pcadapt_pvalues <- as.data.frame(apoly_populations_pca$pvalues) 

pdf("analysis/pcadapt/apoly_populations_pvalues.pdf")
hist(apoly_populations_pca$pvalues, xlab ="p-values", main =NULL, breaks =50, col ="skyblue") 
dev.off() 

apoly_populations_padj <- p.adjust(apoly_populations_pca$pvalues, method ="bonferroni") 
alpha <- 0.1 
outliers <- which(apoly_populations_padj < alpha) 
length(outliers) 

write.table(outliers, file="analysis/pcadapt/apoly_populations_pcaoutliers.txt") 

q()



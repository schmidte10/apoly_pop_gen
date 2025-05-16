library("pcadapt")
library("ggplot2")

apoly_bed <- "../data/gatk.filtered.relaxed_studywide.bicolor.biallelic.test.bed"
apoly_pcadapt <-  read.pcadapt(apoly_bed, type = "bed")

apoly_populations_kplot <- pcadapt(input = apoly_pcadapt, K=10)
pdf("outputs/apoly_populations_kplot.pdf")
plot(apoly_populations_kplot, option = "screeplot")
dev.off() 

apoly_populations_pca <- pcadapt(input =apoly_pcadapt, K =3)
write.table(apoly_populations_pca$scores, "outputs/pcadapt_scores.txt")
summary(apoly_populations_pca)

poplist.names <- read.delim("../data/population_list.txt", header =FALSE)

pdf("outputs/apoly_populations_pca1v2.pdf")
plot(apoly_populations_pca, option = "scores", i=1, j=2, pop = poplist.names$V3)
dev.off()

pdf("outputs/apoly_populations_pca1v3.pdf")
plot(apoly_populations_pca, option = "scores", i=1, j=3, pop = poplist.names$V3)
dev.off()

pdf("outputs/apoly_populations_pca2v3.pdf")
plot(apoly_populations_pca, option = "scores", i=2, j=3, pop = poplist.names$V3)
dev.off()

png("outputs/pcadapt_apoly_manhattan.png")
plot(apoly_populations_pca, option= "manhattan") 
dev.off() 

png("outputs/pcadapt_apoly_qqplot.png") 
plot(apoly_populations_pca, option ="qqplot") 
dev.off()

apoly_pcadapt_pvalues <- as.data.frame(apoly_populations_pca$pvalues) 

pdf("outputs/apoly_populations_pvalues.pdf")
hist(apoly_populations_pca$pvalues, xlab ="p-values", main =NULL, breaks =50, col ="skyblue") 
dev.off() 

apoly_populations_padj <- p.adjust(apoly_populations_pca$pvalues, method ="bonferroni") 
alpha <- 0.1 #[adjust as necessary]
outliers <- which(apoly_populations_padj < alpha) 
length(outliers) 

write.table(outliers, file="outputs/apoly_populations_pcaoutliers.txt") 

q()



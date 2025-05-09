
library(tidyverse)
library(cowplot)
library(stringi)
library(ggpubr)

idepth <- read_tsv("./out.idepth")

imiss <- read_tsv("./out.imiss")

region.list <- factor(rep(c("CORE", "CHAUVEL","LEADING","HERON","LEADING","TORRES STRAIT","REEF HQ"), c(25, 8, 8, 8, 16, 15, 3)))

df <- left_join(idepth,imiss) %>%  
  tidyr::extract(INDV,into="site",regex="([^_]*)",remove = FALSE) %>% 
  mutate(region = region.list) %>%
  arrange(F_MISS) %>% 
  dplyr::mutate(rn = row_number())

mycol=c("#455A64","#D32F2F","#00BCD4","dodgerblue", "#212121","#FF5722")

p1 = ggplot(df) + geom_col(aes(x=reorder(INDV, rn), y=F_MISS, fill=region)) + 
  scale_fill_manual(values =mycol) + 
  xlab("individuals") +
  ylab("Proportion of genotypes missing") +
  theme_pubr() + 
  theme(legend.position = "none", axis.text.x = element_text(angle=90, size = 8))

p2 <- ggplot(df) + geom_col(aes(x=reorder(INDV, rn), y=MEAN_DEPTH, fill=region)) + 
  scale_fill_manual(values =mycol) + 
  xlab("") + 
  ylab("Mean Sequencing Depth") + 
  theme_pubr() + 
  theme(axis.text.x =element_blank(), legend.position ="none")

plot_grid(p2,p1,nrow=2, labels = c("A","B")) 

snpden <- read_tsv("./out.snpden") %>% 
  filter(SNP_COUNT>10) 

ggplot(snpden) + geom_histogram(aes(x=SNP_COUNT, y=..density..), binwidth = 5, color ="black", fill="#455A64") + 
  xlab("SNP count per 10000kb") + ylab("Density") + 
  theme_classic()

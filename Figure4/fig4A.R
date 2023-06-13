# Kmer difference enriched in statistic significant group

library(ggplot2)
library(ggpubr)

# load data
k_fisher <- read.csv("fig4A_diffKmerStat.csv")

# subset mono nucleotide result
k_one <- 
  transform(subset(k_fisher, kmer%in%c("A", "C", "G", "T")),
            kmer = factor(kmer, levels = c("A", "C", "G", "T")),
            utr = factor(utr, levels = c("5'UTR", "3'UTR")))

# plot
ggplot(k_one)+
  geom_errorbar(aes(x = kmer, y = odds, ymax = conf_hi, ymin = conf_lo), width = 0.2, linetype = 1)+
  geom_point(aes(x = kmer, y = odds, color = kmer), size = 3)+
  geom_hline(yintercept = 1, linetype = 3, linewidth = 0.9, color = "red4")+
  facet_wrap(~utr)+
  scale_color_manual(values = c("#1f77b4","#ff7f0e","#2ca02c","#d62728"))+
  theme_pubclean()+
  labs(y = "odds ratio w. 95%CI", subtitle = "Kmer difference enrich in significance group")+
  theme(legend.position = "bottom")

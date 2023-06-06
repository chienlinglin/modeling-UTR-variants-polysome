# library
library(ggplot2)
library(ggpubr)
library(ggrepel)

# load raw data
plotTbl <- read.csv("/work/temp/polysome/fig4aTbl.csv")

# focus on motif, remove none significance & duplicated record
tbl <- unique(transform(plotTbl, utr = factor(utr, levels = c("5'UTR", "3'UTR")), RBP = NULL, matrixName = NULL, motif = NULL, pnt_col=as.character(pnt_col)))

# plot
  ggplot(tbl)+
  geom_jitter(aes(x = log2(risk.r), y = -log10(pval), color = pnt_col), show.legend = FALSE, height = 0.02, width = 0.02, size = 0.5)+
  scale_color_manual(values = c("#d3d3d3", "#ff3300"))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#001e43", size = 0.2)+
  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "#422400", size = 0.2)+
  #geom_text_repel(aes(x = log2(risk.r), y = -log10(pval), label = label), size = 2.75, max.overlaps = 30, segment.size = 0.1)+  ##ver1
  geom_text_repel(aes(x = log2(risk.r), y = -log10(pval), label = label), size = 2.6, max.overlaps = 30, segment.size = 0.1)+
  facet_grid(~utr)+ 
  xlim(-2.1,4.8)+
  labs(x = "log2 risk ratio", y = "-log10(p.value)", title = "Motif enriched in UTRs with significantly changed translatability")+
  theme_bw()

# Compare WT/mt UTR predicted-free energy

library(coin)
library(ggplot2)
library(ggpubr)

# load table
dG_tbl <- read.csv("fig3C_deltaG.csv")

# subset 5'UTR result
# permutation test
dg5U <- subset(dG_tbl, utr=="5utr")
dg5U_diffTest <- oneway_test(dg5U$dG~as.factor(dg5U$exp_results))
dg5U_absTest <- oneway_test(dg5U$abs_dG~as.factor(dg5U$exp_results))

# plot
ggplot(dg5U)+
  geom_boxplot(aes(x = exp_results, y = dG, color = exp_results), fill = NA, size = 0.8)+
  scale_color_manual(values = c("#4e79a7","#f28e2c"))+
  labs(subtitle = sprintf("Fisher-Pitman permutation test p.value = %.3e", coin::pvalue(dg5U_diffTest)), x = NULL, y = NULL)+
  theme_pubclean()+
  coord_flip()+
  theme(legend.position = "bottom", axis.text.x = element_text(size = 12), plot.subtitle = element_text(size = 12))

ggplot(dg5U)+
  geom_boxplot(aes(x = exp_results, y = abs_dG, color = exp_results), fill = NA, size = 0.8)+
  scale_color_manual(values = c("#4e79a7","#f28e2c"))+
  labs(subtitle = sprintf("Fisher-Pitman permutation test p.value = %.3e", coin::pvalue(dg5U_absTest)), x = NULL, y = NULL)+
  theme_pubclean()+
  coord_flip()+
  theme(legend.position = "bottom", axis.text.x = element_text(size = 12), plot.subtitle = element_text(size = 12))


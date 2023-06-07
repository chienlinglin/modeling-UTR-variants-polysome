library(ggplot2)
library(ggpubr)

# load raw table
tbl <- read.csv("/work/temp/polysome/fig3A_cslevel.csv")

# Fisher's Exact Test
cslevel <- table(tbl$pathogenic, tbl$exp_results)
fTest <- fisher.test(cslevel)

# build plot table
csTbl <- 
  data.frame(sig = rep(c("HC_nonsig", "HC_sig"), each=2), diseaseCausing = rep(c("cause", "else")), Data = c(cslevel))
csTbl <- transform(csTbl, 
                   pos = unlist(lapply(split(c(cslevel), f = rep(c("HC_nonsig", "HC_sig"), each=2)), FUN = function(x){tmp <- rev(cumsum(rev(x))); tmp <- (tmp-0.5*x)/sum(x); return(tmp)}))) # calculate label 

# plot
ggplot(csTbl)+
  geom_bar(aes(x = sig, y = Data, fill = diseaseCausing), stat = "identity", position = "fill")+
  geom_text(aes(x = sig, y = pos, label = Data), size = 9, color = "white")+
  scale_fill_manual(values = c("#f28e2c","#4e79a7"))+
  labs(subtitle = sprintf("Fisher's exact test p.value = %.3f", fTest$p.value), x = NULL, y = NULL)+
  theme_pubclean()+
  theme(legend.position = "bottom", axis.text.x = element_text(size = 12), plot.subtitle = element_text(size = 12))


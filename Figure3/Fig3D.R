# Difference maximum stem-loop coverage

library(ggplot2)
library(ggpubr)

# load data
slCov <- read.csv("/work/temp/polysome/fig3D_stemLoopCov.csv")

# calculate mt-WT stem-loop coverage
diffTbl <- subset(slCov, utr=="utr5")
diffTbl$type <- factor(diffTbl$type, levels = c("WT", "mt"))

diffTbl <- 
  do.call(rbind, lapply(split(diffTbl, diffTbl$grpTbl), FUN = function(x){
    x <- x[order(x$type),]
    rst <- data.frame(grpNo = unique(x$grpTbl), diffCov = diff(x$covRate), utr = unique(x$utr), exp_results = unique(x$exp_results))
    return(rst)
  }))


# calculate summary value
plotTbl <- 
  do.call(rbind, tapply(diffTbl$diffCov, diffTbl$exp_results, simplify = FALSE, FUN = function(x){
    return(data.frame(oriMean = mean(x), oriSd = sd(x), oriSE=sd(x)/(sqrt(length(x))), absMean = mean(abs(x)), absSd = sd(abs(x)), absSE=sd(abs(x))/(sqrt(length(x)))))
  }))

plotTbl <- transform(plotTbl, group = c("HC_nonsig", "HC_sig"))

# plot
ggplot(plotTbl)+
  geom_hline(yintercept = 0, linewidth = 0.9)+
  geom_bar(aes(x = group, y = oriMean, fill = group, color = group), stat = "identity", width = 0.7)+
  geom_errorbar(aes(x = group, ymin = oriMean-oriSE, ymax = oriMean+oriSE), linewidth = 0.5, width = 0.2)+
  scale_fill_manual(values = c("#4e79a7","#f28e2c"))+
  scale_color_manual(values = c("#4e79a7","#f28e2c"))+
  theme_pubclean()+
  labs(x = "significance", y = "mean w. se", subtitle = "5'UTR mt-WT stem_loop coverage")+
  theme(legend.position = "bottom")

ggplot(plotTbl)+
  geom_hline(yintercept = 0, linewidth = 0.9)+
  geom_bar(aes(x = group, y = absMean, fill = group, color = group), stat = "identity", width = 0.7)+
  geom_errorbar(aes(x = group, ymin = absMean-absSE, ymax = absMean+absSE), linewidth = 0.5, width = 0.2)+
  scale_fill_manual(values = c("#4e79a7","#f28e2c"))+
  scale_color_manual(values = c("#4e79a7","#f28e2c"))+
  theme_pubclean()+
  labs(x = "significance", y = "abs mean w. se", subtitle = "5'UTR mt-WT stem_loop coverage\nabsolute value")+
  theme(legend.position = "bottom")

# library
library(coin)
library(ggplot2)
library(ggrepel)

# load raw data
# significanse group
poly_seq.list <- read.csv("fig4ab_seqList.csv")
poly_seq.list <- split(poly_seq.list$ID~poly_seq.list$group)

# diff motif
rbp_ptrn <- read.csv("fig4a_diffMotif.csv")

# diff gene symbol
rbp_gs <- read.csv("fig4b_diffGeneSymbol.csv")

# motif table
mMotif <- read.delim("ATtRACT_motif.tsv")

# analysis
enrich_stat <- 
  sapply(c("rbp_gs", "rbp_ptrn"), simplify = FALSE, FUN = function(var){
    tbl <- get(var); tbl <- split(tbl, tbl$UTR_Region)                                              # load data and seperate by UTR
    stat_rst <- do.call(rbind, lapply(tbl, FUN = function(utr){                                     # enrichment analysis, sig vs all stat grp
      utr_p <- unique(utr$UTR_Region)
      all.tbl <- utr[utr$Group%in%poly_seq.list$all&utr$RBP_Type!="constant", 5]                    # extract each feature
      sig.tbl <- utr[utr$Group%in%poly_seq.list$sig&utr$RBP_Type!="constant", 5]
      all.tbl <- table(strsplit(paste(all.tbl, collapse = ","), split = ","))
      all.tbl <- all.tbl[nchar(names(all.tbl))>0]
      sig.tbl <- table(strsplit(paste(sig.tbl, collapse = ","), split = ","))
      sig.tbl <- sig.tbl[nchar(names(sig.tbl))>0]
      sig_name <- names(sig.tbl)
      dt.array <- array(data = c(rbind(sig.tbl, rep(sum(sig.tbl), length(sig.tbl)), all.tbl[sig_name], rep(sum(all.tbl), length(sig.tbl)))), 
                        dim = c(2,2, length(sig_name)))
      rst <- do.call(rbind, lapply(seq_along(sig_name), FUN = function(idx){
        stat_rst <- fisher.test(dt.array[,,idx])                          
        return(c(var = var, utr = utr_p, feature = sig_name[idx], pval = as.numeric(stat_rst$p.value), risk.r = as.numeric(stat_rst$estimate)))
      }))
      return(rst)
    }))
    stat_rst <- transform(as.data.frame(stat_rst), pval = as.numeric(pval), risk.r = as.numeric(risk.r))
    return(stat_rst)
  })


# fig4B
tbl4A <- merge(enrich_stat$rbp_ptrn, subset(mMotif, select = c("Matrix_id", "Motif")), by.x = "feature", by.y = "Matrix_id", all.x = TRUE)
tbl4A <- 
  unique(transform(tbl4A, 
                   utr = factor(utr, levels = c("5'UTR", "3'UTR")), 
                   label = ifelse(pval<0.05&risk.r>1.5, Motif, NA),
                   pnt_col= ifelse(pval<0.05&risk.r>1, "1", "0"),
                   feature = NULL, Motif = NULL))
  
ggplot(tbl4A)+
  geom_jitter(aes(x = log2(risk.r), y = -log10(pval), color = pnt_col), show.legend = FALSE, height = 0.02, width = 0.02, size = 0.5)+
  scale_color_manual(values = c("#d3d3d3", "#ff3300"))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#001e43", size = 0.2)+
  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "#422400", size = 0.2)+
  geom_text_repel(aes(x = log2(risk.r), y = -log10(pval), label = label), size = 2.6, max.overlaps = 30, segment.size = 0.1)+
  facet_grid(~utr)+ 
  xlim(-2.1,4.8)+
  labs(x = "log2 risk ratio", y = "-log10(p.value)", title = "Motif enriched in UTRs with significantly changed translatability")+
  theme_bw()  

# fig4C
tbl4B <- 
  transform(enrich_stat$rbp_gs, 
                   utr = factor(utr, levels = c("5'UTR", "3'UTR")), 
                   label = ifelse(pval<0.05&risk.r>1.5, feature, NA),
                   pnt_col= ifelse(pval<0.05&risk.r>1, "1", "0"))

ggplot(tbl4B)+
  geom_jitter(aes(x = log2(risk.r), y = -log10(pval), color = pnt_col), show.legend = FALSE, height = 0.02, width = 0.02, size = 0.5)+
  scale_color_manual(values = c("#d3d3d3", "#ff3300"))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#001e43", size = 0.2)+
  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "#422400", size = 0.2)+
  geom_text_repel(aes(x = log2(risk.r), y = -log10(pval), label = label), size = 2.75, segment.size = 0.1, max.overlaps = 30)+
  facet_grid(~utr)+ 
  labs(x = "log2 risk ratio", y = "-log10(p.value)", title = "RBP binding enriched in UTRs with significantly changed translatability")+
  theme_bw()

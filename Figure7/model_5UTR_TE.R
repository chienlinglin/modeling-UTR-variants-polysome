# =============================================================================
# model_5UTR_TE.R
#
# The 5' UTR translation-efficiency model, and the four figures that show it:
# the coefficients, the fit on held-out sequences, and the two luciferase
# comparisons. Everything is written to OUT_DIR, which is separate from the
# working analysis, so nothing here overwrites it.
#
#   TE = log2( Poly / (Free + 40S) ),  Poly = F2..F8,  40S = F040
#
# WHAT IT DOES
#   1. TE of every wild-type 5' UTR sequence, with a weight for how precisely it
#      was measured.
#   2. Features reduced to weakly correlated representatives: single-linkage
#      clustering on 1 - |Spearman|, the member most strongly associated with TE
#      as the representative of each cluster, and the tree cut until every
#      representative has a variance inflation factor below 5.
#   3. A single split of the sequences into a training two thirds and a testing
#      third, fixed by SEED.
#   4. Stability selection on the training part: NUM_SS half-samples, a weighted
#      LASSO in each with the penalty chosen by ten-fold cross-validation, and
#      the features selected in more than STABLE_THRESHOLD of them.
#   5. Those features refitted by weighted linear regression on the same
#      training part. Its coefficients are the model; the LASSO only selects.
#   6. The testing third predicted, and the luciferase comparisons made.
#
#   Every random draw is governed by SEED, so a rerun reproduces the figures
#   exactly.
#
# WHAT IT WRITES
#   figures/coefficients(_with_probability)  the features of the model with their
#                         95% intervals and selection probabilities (Figure 7C),
#                         and the sequence logos of its RBP motifs
#   model.rds             the fitted model, for export_figure_data.R to read
#   sessionInfo.txt
#
#   Figures 7B, 7D and 7E are drawn by Fig7B.R, Fig7D.R and Fig7E.R from the
#   small tables that export_figure_data.R writes, so that those panels can be
#   redrawn without the raw data.
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse); library(glmnet); library(ggrepel); library(patchwork)
})

## ---- where the data is, and where the results go ---------------------------------
ROOT     <- "/mnt/cloudBackup/GoogleDrive_thejovialjoy/Projects/rna_utranstab"
DATA     <- file.path(ROOT, "polysome_redo", "prepared_data_v2.RData")
FEAT     <- file.path(ROOT, "HEK_fullDataNoCut_221118.rda")
PWM      <- file.path(ROOT, "human_RBP_pwm.rdata")
OUT_DIR  <- file.path(ROOT, "polysome_redo", "wt_lasso", "public_5UTR_TE_model")
FIG_DIR  <- file.path(OUT_DIR, "figures")
for (dd in c(OUT_DIR, FIG_DIR)) if (!dir.exists(dd)) dir.create(dd, recursive = TRUE)

# The training and testing sets used for the published figures are supplied as a
# file rather than drawn here, so that the figures reproduce exactly. Point this
# at nothing, or delete the file, and the script falls back to its own split,
# fixed by SEED.
SPLIT <- local({
  cands <- c("split.rds", file.path(OUT_DIR, "split.rds"),
             file.path(ROOT, "polysome_redo", "wt_lasso", "FINAL_5UTR_TE_stability",
                       "data", "data_reported_split.rds"))
  hit <- cands[file.exists(cands)]
  if (length(hit)) normalizePath(hit[1]) else ""
})


# the wild-type luciferase measurements; the file is looked for in the places it
# has lived, so a different checkout does not need the path edited
LUC_WT_FILE <- local({
  cands <- c("poly_lucy_validate.csv",
             file.path(ROOT, "poly_lucy_validate.csv"),
             file.path(ROOT, "polysome_redo", "poly_lucy_validate.csv"),
             file.path(ROOT, "polysome_redo", "wt_lasso", "poly_lucy_validate.csv"),
             file.path(OUT_DIR, "poly_lucy_validate.csv"))
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop("poly_lucy_validate.csv not found; set LUC_WT_FILE by hand")
  normalizePath(hit[1])
})


## ---- settings ---------------------------------------------------------------------
# the outcome
MIN_COUNT        <- 1       # reads required on each side of the ratio, per replicate
MIN_BATCH        <- 2       # replicates a sequence must pass the filter in
WEIGHT_POW       <- 0.5     # weight of a sequence is its mean read count to this power

# the candidate features
MAX_NA           <- 0.10    # dropped when missing in more than this share of groups
MIN_NONZERO      <- 0.05    # dropped when non-zero in this share of groups or fewer
CUTS_HCT         <- 50:200  # clusters tried, from the largest number down
VIF_THRESH       <- 5       # every representative must fall below this

# the split and the selection
TEST_FRACTION    <- 1/3     # a 2:1 split into training and testing
NUM_SS           <- 100     # half-samples of the training part
NFOLD            <- 10      # folds of cv.glmnet inside each half-sample
LAMBDA_RULE      <- "lambda.min"
STABLE_THRESHOLD <- 0.7     # a feature is kept above this selection probability
SEED             <- 41      # prime; every draw below is derived from it

# the figures
COEF_ROW_MM  <- 3.3      # height of one feature row, mm
COEF_BAR_MM  <- 62       # width of the coefficient panel, mm
COEF_FREQ_MM <- 20       # width of the selection-probability panel, mm
COEF_GAP_MM  <- 9        # space between the two, mm
COEF_BAND    <- 0.52     # thickness of an interval, as a share of the row
COEF_DOT     <- 0.40     # size of an estimate, as a share of the row
COEF_LEGEND  <- FALSE    # the colours are explained by the axis
CLASS_PT     <- 10       # size of the class names beside the ribbon
TOP_N_FEAT   <- 30       # features in the shorter coefficient figure
NEIGHBOUR_R  <- 0.02     # density radius in the observed-against-predicted figure
DENS_MIN     <- 6
DENS_MAX     <- 60

THR     <- STABLE_THRESHOLD
THR_TXT <- sub("0\\.", ".", format(THR))
TAG     <- "log2(Poly/(Free+F040)), 5' UTR"

GENE_OF <- c(G3219 = "IRF6", G3320 = "DMD",  G3387 = "NF1",  G3935 = "MUTYH",
             G4119 = "TP53", G4928 = "FTL",  G5598 = "LGI1", G6355 = "RFXAP")

TAG <- "[5'] log2(Poly/(Free+F040))"

cat("[Load]", basename(DATA), "\n"); load(DATA)

set.seed(SEED)
cat("[Load]", basename(DATA), "\n"); load(DATA)
cat("[Load]", basename(FEAT), "\n"); load(FEAT)
for (o in c("combine_all2", "groupTbl", "ftTbl", "mean_df"))
  if (!exists(o)) stop("object not found in the data files: ", o)


if (!"Free_F040" %in% names(combine_all2))
  combine_all2 <- combine_all2 %>% dplyr::mutate(Free_F040 = Free + F040)

utr_of <- function(s) dplyr::case_when(substr(s,1,1)=="5" ~ "5'",
                                       substr(s,1,1)=="3" ~ "3'", TRUE ~ NA_character_)

## =============================================================================
## 1. the data: outcome, features, and the representatives of each cluster
## =============================================================================

## ---- 1. bridge and RBP names (identical to the pipeline) -----------------------
build_bridge <- function(gt) {
  nm <- names(gt)
  if (all(c("grpTbl","seqName","type") %in% nm)) {
    gt %>% dplyr::filter(type == "WT") %>% dplyr::select(Group = grpTbl, seqName = seqName)
  } else if (all(c("grpNo","WT") %in% nm)) {
    gt %>% dplyr::select(Group = grpNo, seqName = WT)
  } else stop("Unrecognised groupTbl columns: ", paste(nm, collapse=", "))
}
gt_wt <- build_bridge(groupTbl) %>%
  dplyr::filter(!is.na(seqName), !is.na(Group)) %>% dplyr::distinct()

lab_of <- (function() {
  if (!exists("ATtRACT_dbInfo")) return(function(ids) ids)
  db <- as.data.frame(ATtRACT_dbInfo)
  db$rbp_ID  <- paste("RBP", db$Matrix_id, sep = "_")
  db$motif   <- gsub("T", "U", as.character(db$Motif))
  db$gene_nm <- as.character(db$Gene_name)
  ann <- db %>% dplyr::group_by(rbp_ID) %>%
    dplyr::summarise(mm = motif[which.max(nchar(motif))],
                     gn = paste(unique(na.omit(gene_nm[gene_nm != ""])), collapse = ";"),
                     .groups = "drop") %>%
    dplyr::mutate(dn = ifelse(gn != "", paste0(mm, " (", gn, ")"), mm))
  ann <- dplyr::bind_rows(ann, ann %>% dplyr::mutate(rbp_ID = paste0(rbp_ID, "_WT"))) %>%
    dplyr::distinct(rbp_ID, .keep_all = TRUE)
  function(ids) { m <- ann$dn[match(ids, ann$rbp_ID)]; ifelse(is.na(m), ids, m) }
})()

# miRNA features carry the family in their name; the seed sequence comes from the
# TargetScan family file, so that a site is shown the way an RBP motif is:
#   miRNA_m7.1a_let-7-5p/98-5p  ->  GAGGUAG 7mer-1A (let-7-5p/98-5p)
SITE_NAME <- c(m6 = "6mer", m7.1a = "7mer-1A", m7.m8 = "7mer-m8", m8 = "8mer")
mir_lab <- (function() {
  ts <- list.files(ROOT, pattern = "targetscan|miRFamily", recursive = TRUE,
                   full.names = TRUE, ignore.case = TRUE)
  seeds <- NULL
  if (length(ts)) {
    x <- tryCatch(readr::read_delim(ts[1], delim = "\t", show_col_types = FALSE),
                  error = function(e) NULL)
    if (!is.null(x)) {
      fc <- names(x)[grepl("family", names(x), ignore.case = TRUE)][1]
      sc <- names(x)[grepl("seed", names(x), ignore.case = TRUE)][1]
      if (!is.na(fc) && !is.na(sc))
        seeds <- tibble::tibble(key = gsub("[/.-]", "_", as.character(x[[fc]])),
                                seed = as.character(x[[sc]])) %>%
        dplyr::distinct(key, .keep_all = TRUE)
    }
  }
  function(ids) {
    p <- sub("^miRNA_", "", ids)
    site <- sub("_.*$", "", p); fam <- sub("^[^_]*_", "", p)
    site <- ifelse(is.na(SITE_NAME[site]), site, SITE_NAME[site])
    sd <- if (is.null(seeds)) rep(NA_character_, length(ids))
    else seeds$seed[match(gsub("[/.-]", "_", fam), seeds$key)]
    ifelse(is.na(sd), sprintf("%s %s", fam, site), sprintf("%s %s (%s)", sd, site, fam))
  }
})()

## ---- 2. outcomes and features (identical to the pipeline) ---------------------
frac_b <- combine_all2 %>%
  dplyr::mutate(Mono = F1, Poly = F2+F3+F4+F5+F6+F7+F8,
                F040c = F040, Freec = Free)
outcome_def <- list(
  D1 = list(num = quote(Poly), den = quote(F040c),         lab = "log2(Poly/F040)"),
  D2 = list(num = quote(Poly), den = quote(Freec + F040c), lab = "log2(Poly/(Free+F040))"))

seq_outcome <- function(oc, which_type) {
  d <- outcome_def[[oc]]; ne <- d$num; de <- d$den
  frac_b %>% dplyr::filter(type == which_type) %>%
    dplyr::mutate(num = !!ne, den = !!de) %>%
    dplyr::filter(num >= MIN_COUNT, den >= MIN_COUNT) %>%
    dplyr::mutate(r = log2(num/den), M = num + den) %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(y = mean(r, na.rm = TRUE), Mbar = mean(M, na.rm = TRUE),
                     B_kept = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(is.finite(y), B_kept >= MIN_BATCH)
}
mt_lab <- setdiff(unique(combine_all2$type), "WT")[1]
grp_seq_utr <- combine_all2 %>% dplyr::filter(type == "WT") %>%
  dplyr::distinct(Group, seqName) %>% dplyr::mutate(utr = utr_of(seqName))

# mutant oligos of each variant group (for predicting mutant TE)
gt_mt <- if (all(c("grpTbl", "seqName", "type") %in% names(groupTbl))) {
  groupTbl %>% dplyr::filter(type != "WT") %>%
    dplyr::select(Group = grpTbl, seqName) %>%
    dplyr::filter(!is.na(seqName), !is.na(Group)) %>% dplyr::distinct()
} else NULL

# Features are transformed, averaged per variant group, filtered, imputed and
# standardized on the WILD-TYPE sequences. The mutant sequences are put on the
# same scale with the wild-type medians, means and standard deviations, so that
# one set of coefficients can be applied to both.
# Upstream AUG counts, parsed from the uAUG(KOZAK) table: "=" is none, anything
# else is semicolon-separated positions, so the count is how many positions are
# listed. Only counts are used; a sequence without a uAUG has no position, and any
# filler value would be read as a real distance. The "optimal" class is present
# exactly once in every sequence, which is the annotated start codon rather than
# an upstream one, so the total is taken over the other three classes; the column
# itself is kept so that the filter log records it being dropped for no variance.
uaug_tbl <- local({
  ka <- if (exists("uAUG(KOZAK)")) get("uAUG(KOZAK)") else NULL   # exists() searches the calling environments
  if (is.null(ka)) { cat("[Features] uAUG(KOZAK) not found; upstream AUGs not added\n"); return(NULL) }
  cls <- c("optimal", "strong", "moderate", "weak")
  cls <- cls[cls %in% names(ka)]
  cnt <- vapply(ka[cls], function(z) {
    z <- as.character(z)
    vapply(z, function(v) if (is.na(v) || v %in% c("=", "")) 0L else
      sum(nzchar(trimws(strsplit(v, ";", fixed = TRUE)[[1]]))), integer(1))
  }, integer(nrow(ka)))
  out <- as.data.frame(cnt); names(out) <- paste0("uAUG_n_", cls)
  out$seqName <- ka$seqName
  # Every sequence carries exactly one "optimal" AUG, which is what the start
  # codon looks like rather than an upstream one, so its count has no variance
  # and is dropped. Its position is kept in its place.
  pos <- vapply(as.character(ka$optimal), function(v) {
    if (is.na(v) || v %in% c("=", "")) return(NA_real_)
    p <- suppressWarnings(as.numeric(strsplit(v, ";", fixed = TRUE)[[1]]))
    if (any(is.finite(p))) max(p[is.finite(p)]) else NA_real_
  }, numeric(1))
  out$uAUG_n_optimal <- NULL
  out$uAUG_optimal_pos <- unname(pos)
  out
})
# the aggregated miRNA site counts, in place of the 749 single-site columns
mir_tbl <- local({
  mi <- if (exists("miRNA")) miRNA else NULL
  if (is.null(mi)) { cat("[Features] miRNA table not found; miRNA counts not added\n"); return(NULL) }
  cols <- grep("_sum$", names(mi), value = TRUE)
  if (!length(cols)) return(NULL)
  out <- as.data.frame(lapply(mi[cols], function(z) as.numeric(as.character(z))))
  names(out) <- paste0("miR_", cols); out$seqName <- mi$seqName
  out
})

prep_ft <- function(u) {
  raw <- as.data.frame(ftTbl[[u]]); raw$seqName <- rownames(ftTbl[[u]])
  # upstream AUGs and miRNA sites join here, so they go through every filter,
  # the square-root transform and the standardisation with everything else
  if (!is.null(uaug_tbl)) raw <- dplyr::left_join(raw, uaug_tbl, by = "seqName")
  if (!is.null(mir_tbl))  raw <- dplyr::left_join(raw, mir_tbl,  by = "seqName")
  all_feat <- setdiff(names(raw), "seqName")
  # Nothing is excluded by kind: every column of the feature table is offered,
  # including the single miRNA sites and the G-quadruplex features, and the
  # filters below decide. Counts take a square root, conservation and the
  # position take a log, free energies are left as they are.
  cnt  <- grep("^RBP_|^ARE_|^kmer|stemLoopAmt|^uAUG_n_|^miR_|^miRNA|^RG4$", names(raw), value = TRUE)
  wayc <- grep("way$|Rate$", names(raw), value = TRUE)
  logc <- grep("^uAUG_optimal_pos$", names(raw), value = TRUE)
  df <- raw %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(cnt), sqrt)) %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(wayc), ~log(.x + 0.02))) %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(logc), ~log(.x + 1)))
  by_group <- function(bridge) df %>% dplyr::inner_join(bridge, by = "seqName") %>%
    dplyr::select(-seqName) %>% dplyr::group_by(Group) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x, na.rm = TRUE)),
                     .groups = "drop")
  wt <- by_group(gt_wt)
  g <- wt$Group; m <- dplyr::select(wt, -Group)
  
  # three filters, applied in this order and recorded one by one
  na_rate <- vapply(m, function(x) mean(is.na(x)), numeric(1))
  nz_rate <- vapply(m, function(x) sum(x != 0, na.rm = TRUE)/length(x), numeric(1))
  sd_val  <- vapply(m, function(x) sd(x, na.rm = TRUE), numeric(1))
  drop_na <- na_rate > MAX_NA
  drop_nz <- !drop_na & nz_rate <= MIN_NONZERO
  drop_sd <- !drop_na & !drop_nz & (!is.finite(sd_val) | sd_val == 0)
  keep <- !(drop_na | drop_nz | drop_sd)
  
  log_tbl <- tibble::tibble(
    feature = all_feat,
    missing_rate = unname(na_rate[all_feat]),
    non_zero_rate = unname(nz_rate[all_feat]),
    sd = unname(sd_val[all_feat]),
    dropped_at = dplyr::case_when(
      all_feat %in% names(m)[drop_na]    ~ sprintf("1. missing in more than %.0f%% of groups", 100*MAX_NA),
      all_feat %in% names(m)[drop_nz]    ~ sprintf("2. non-zero in %.0f%% of groups or fewer", 100*MIN_NONZERO),
      all_feat %in% names(m)[drop_sd]    ~ "3. no variance",
      TRUE                               ~ "kept"))
  
  m <- m[, keep, drop = FALSE]
  med <- vapply(m, median, numeric(1), na.rm = TRUE)
  m <- as.data.frame(Map(function(x, md) { x[is.na(x)] <- md; x }, m, med))
  mu <- vapply(m, mean, numeric(1)); sdv <- vapply(m, sd, numeric(1))
  out_wt <- as.data.frame(Map(function(x, a, b) (x - a)/b, m, mu, sdv)); out_wt$Group <- g
  out_mt <- NULL
  if (!is.null(gt_mt)) {
    mt <- by_group(gt_mt)
    mm <- as.data.frame(mt[, names(keep)[keep], drop = FALSE])
    mm <- as.data.frame(Map(function(x, md, a, b) { x[is.na(x)] <- md; (x - a)/b },
                            mm, med, mu, sdv))
    names(mm) <- names(out_wt)[names(out_wt) != "Group"]; mm$Group <- mt$Group
    out_mt <- mm
  }
  # the standardized name, so the log can be matched to the model later
  log_tbl$feature_id <- make.names(log_tbl$feature)
  log_tbl$source <- dplyr::case_when(grepl("^uAUG_n_", log_tbl$feature) ~ "uAUG(KOZAK)",
                                     grepl("^miR_", log_tbl$feature)    ~ "miRNA (aggregated)",
                                     TRUE                               ~ "ftTbl")
  list(wt = out_wt, mt = out_mt, log = log_tbl)
}
ft5 <- prep_ft("U5")
feat    <- list("5'" = ft5$wt)
feat_mt <- ft5$mt
feat_log <- ft5$log
cat("[Features] from", nrow(feat_log), "columns (feature table plus upstream AUGs and miRNA counts):\n")
print(as.data.frame(feat_log %>% dplyr::count(source, name = "columns")), row.names = FALSE)
print(as.data.frame(feat_log %>% dplyr::count(dropped_at, name = "features") %>%
                      dplyr::arrange(dropped_at)), row.names = FALSE)
# where a feature of interest ends up
look_for <- "uATG|uAUG|[Kk]ozak"
hits <- feat_log %>% dplyr::filter(grepl(look_for, feature))
cat(sprintf("  features matching '%s': %s\n", look_for,
            if (nrow(hits)) paste(sprintf("%s (%s)", hits$feature, hits$dropped_at), collapse = "; ")
            else "none in the feature table"))
cat(sprintf("[Features] %d wild-type and %s mutant variant groups, %d features\n",
            nrow(feat[["5'"]]), if (is.null(feat_mt)) "no" else nrow(feat_mt),
            ncol(feat[["5'"]]) - 1))

fast_vif_select <- function(dm, w = NULL) {
  y <- dm$outcome; Xf <- dm %>% dplyr::select(-outcome)
  if (is.null(w)) w <- rep(1, length(y))
  est <- vapply(Xf, function(x) tryCatch(abs(coef(lm(y ~ x, weights = w))[2]),
                                         error = function(e) NA_real_), numeric(1))
  cc <- suppressWarnings(cor(Xf, method = "spearman", use = "pairwise.complete.obs"))
  cc[!is.finite(cc)] <- 0
  hc <- hclust(as.dist(1 - abs(cc)), method = "single")
  best <- NULL; best_cl <- NULL
  for (k in rev(CUTS_HCT)) {
    if (k >= ncol(Xf)) next
    cl <- cutree(hc, k = k)
    reps <- unname(unlist(tapply(seq_along(cl), cl,
                                 function(i) names(est)[i][which.max(est[i])])))
    reps <- reps[!is.na(reps)]; if (length(reps) < 2) next
    v <- tryCatch(car::vif(lm(reformulate(reps, "y"),
                              data = cbind(y = y, Xf[reps]), weights = w)),
                  error = function(e) NA)
    if (all(is.finite(v)) && max(v) < VIF_THRESH) { best <- reps; best_cl <- cl; break }
  }
  if (is.null(best)) {
    best_cl <- cutree(hc, k = min(min(CUTS_HCT), ncol(Xf) - 1))
    best <- unname(unlist(tapply(seq_along(best_cl), best_cl,
                                 function(i) names(est)[i][which.max(est[i])])))
    best <- best[!is.na(best)]
  }
  list(sel = best, cormat = cc,
       membership = tibble::tibble(feature = names(best_cl),
                                   cluster = as.integer(best_cl)) %>%
         dplyr::mutate(is_rep = feature %in% best))
}

build_d <- function(oc, ut) {
  seq_outcome(oc, "WT") %>% dplyr::inner_join(grp_seq_utr, by = "Group") %>%
    dplyr::filter(utr == ut) %>%
    dplyr::arrange(dplyr::desc(Mbar)) %>% dplyr::distinct(seqName, .keep_all = TRUE) %>%
    dplyr::select(Group, seqName, y, Mbar) %>%
    dplyr::inner_join(feat[[ut]], by = "Group") %>%
    { mu <- mean(.$y); sg <- sd(.$y)
    dplyr::filter(., y > mu - 3*sg, y < mu + 3*sg) }
}

# the labels the coefficient table and the figures show
pretty_feat <- function(lbl, id) {
  out <- as.character(lbl); id <- as.character(id)
  mi <- grepl("^miRNA_", id);      out[mi] <- mir_lab(id[mi])
  ms <- grepl("^miR_", id)
  out[ms] <- sprintf("%s sites, all families",
                     ifelse(is.na(SITE_NAME[sub("_sum$", "", sub("^miR_", "", id[ms]))]),
                            sub("^miR_", "", id[ms]),
                            SITE_NAME[sub("_sum$", "", sub("^miR_", "", id[ms]))]))
  ua <- grepl("^uAUG_n_", id)
  out[ua] <- sprintf("uAUG, %s Kozak", sub("^uAUG_n_", "", id[ua]))
  out[id == "uAUG_optimal_pos"] <- "Position of the optimal AUG"
  out[id == "RG4"]   <- "G-quadruplex count"
  out[id == "RG4FE"] <- "G-quadruplex free energy"
  km <- grepl("^kmer_", id) & lbl == id
  kseq <- gsub("T", "U", sub("^kmer_", "", id[km]))
  out[km] <- sprintf("%d-mer %s", nchar(kseq), kseq)    # 3-mer CGG, 2-mer CG, 1-mer C
  out[id == "GCcontent"] <- "GC content"
  out
}

## =============================================================================
## 3. the model: one split, stability selection, weighted linear regression
## =============================================================================
cat(sprintf("\n[Model] %s\n", TAG))
d <- build_d("D2", "5'")
w <- d$Mbar^WEIGHT_POW; w <- w/mean(w)
dmod <- d %>% dplyr::select(-Group, -seqName, -Mbar) %>% dplyr::rename(outcome = y)
vs  <- fast_vif_select(dmod, w = w); sel <- vs$sel
X <- as.matrix(d[, sel, drop = FALSE]); y <- d$y
N <- nrow(d)
cat(sprintf("  n = %d sequences | %d representatives of %d features\n",
            N, length(sel), nrow(vs$membership)))

# the split: read from SPLIT when it is there, otherwise drawn here and fixed by
# SEED. Rows are matched on the sequence name, which identifies them whatever
# order the data arrive in.
if (file.exists(SPLIT)) {
  sp <- readRDS(SPLIT)
  tr_rows <- which(d$seqName %in% sp$train_seqName)
  te_rows <- which(d$seqName %in% sp$test_seqName)
  unmatched <- length(sp$train_seqName) + length(sp$test_seqName) -
    (length(tr_rows) + length(te_rows))
  if (unmatched) cat(sprintf("  note: %d sequences of the supplied split are not in this data\n",
                             unmatched))
  if (!length(tr_rows) || !length(te_rows)) stop("the supplied split matches nothing here")
  cat(sprintf("  split read from %s\n", basename(SPLIT)))
} else {
  set.seed(SEED)
  te_rows <- sort(sample(N, round(N*TEST_FRACTION)))
  tr_rows <- setdiff(seq_len(N), te_rows)
  cat(sprintf("  split drawn here, fixed by seed %d\n", SEED))
}
cat(sprintf("  %d sequences to fit on, %d held out\n", length(tr_rows), length(te_rows)))

# stability selection on the training part alone. The seed is the one the split
# file carries when there is one, so that the same half-samples are drawn and the
# same features come out; otherwise it follows SEED.
STAB_SEED <- if (exists("sp") && !is.null(sp$stability_seed)) sp$stability_seed else SEED + 1
cat(sprintf("  half-samples drawn with seed %d\n", STAB_SEED))
set.seed(STAB_SEED)
h  <- floor(length(tr_rows)/2)
pb <- utils::txtProgressBar(min = 0, max = NUM_SS, style = 3, width = 36)
sel_mat <- vapply(seq_len(NUM_SS), function(b) {
  sub    <- sample(tr_rows, h)
  foldid <- sample(rep(seq_len(NFOLD), length.out = h))
  fit <- cv.glmnet(X[sub, , drop = FALSE], y[sub], weights = w[sub], alpha = 1,
                   family = "gaussian", foldid = foldid)
  utils::setTxtProgressBar(pb, b)
  as.matrix(coef(fit, s = LAMBDA_RULE))[colnames(X), 1] != 0
}, logical(ncol(X)))
close(pb)
rownames(sel_mat) <- colnames(X); colnames(sel_mat) <- paste0("SS_", seq_len(NUM_SS))
sel_freq <- rowMeans(sel_mat)
stable   <- names(sel_freq)[sel_freq > THR]
if (!length(stable)) stop("no feature reached the threshold")

# the reported model: those features, refitted by weighted linear regression
fit_stable <- function(rows, feats) {
  df <- data.frame(y = y[rows], X[rows, feats, drop = FALSE], check.names = FALSE)
  lm(y ~ ., data = df, weights = w[rows])
}
fit_lm <- fit_stable(tr_rows, stable)
cm <- coef(fit_lm); names(cm) <- gsub("`", "", names(cm)); cm[is.na(cm)] <- 0
intercept <- unname(cm[1])
beta <- setNames(numeric(ncol(X)), colnames(X)); beta[names(cm)[-1]] <- cm[-1]
ci <- suppressWarnings(confint(fit_lm)); rownames(ci) <- gsub("`", "", rownames(ci))
sm <- summary(fit_lm)$coefficients; rownames(sm) <- gsub("`", "", rownames(sm))
vif <- if (length(stable) > 1 && requireNamespace("car", quietly = TRUE))
  car::vif(fit_lm) else setNames(rep(NA_real_, length(stable)), stable)
names(vif) <- gsub("`", "", names(vif))

cf <- tibble::tibble(
  feature_id = stable,
  label   = pretty_feat(lab_of(stable), stable),   # RBP ids become motif and protein
  coef    = unname(beta[stable]),
  se      = unname(sm[stable, "Std. Error"]),
  ci_low  = unname(ci[stable, 1]),
  ci_high = unname(ci[stable, 2]),
  p_value = unname(sm[stable, "Pr(>|t|)"]),
  vif     = unname(vif[stable]),
  sel_freq = unname(sel_freq[stable])) %>%
  dplyr::mutate(direction = ifelse(coef >= 0, "increases", "decreases"))

# predictions, from these coefficients alone
predict_TE <- function(newdata)
  as.numeric(intercept + as.matrix(newdata[, names(beta), drop = FALSE]) %*% beta)
pred_te <- as.numeric(intercept + X[te_rows, , drop = FALSE] %*% beta)
best <- tibble::tibble(Group = d$Group[te_rows], seqName = d$seqName[te_rows],
                       observed = y[te_rows], predicted = pred_te)
rs <- suppressWarnings(cor.test(best$predicted, best$observed, method = "spearman"))
rp <- suppressWarnings(cor.test(best$predicted, best$observed))
cat(sprintf("  %d features above %s of the %d half-samples; max VIF %.2f\n",
            length(stable), THR_TXT, NUM_SS, max(cf$vif, na.rm = TRUE)))
cat(sprintf("  held out: Spearman %.3f (P = %.3g), Pearson %.3f\n",
            rs$estimate, rs$p.value, rp$estimate))
fig_cf <- cf; fig_best <- best

## ---- the luciferase comparisons --------------------------------------------------
te_wt <- seq_outcome("D2", "WT")   %>% dplyr::select(Group, TE_WT = y)
te_mt <- seq_outcome("D2", mt_lab) %>% dplyr::select(Group, TE_mutant = y)
te_5 <- dplyr::inner_join(te_wt, te_mt, by = "Group") %>%
  dplyr::mutate(delta_TE = TE_mutant - TE_WT) %>%
  dplyr::inner_join(grp_seq_utr, by = "Group") %>%
  dplyr::filter(utr == "5'", is.finite(delta_TE))
luc_mt <- purrr::map_dfr(c(Library = "A", `Full-length` = "B"), function(k) {
  mean_df[[k]] %>% dplyr::filter(Type == "5'") %>%
    dplyr::select(Group, Gene, luciferase_log2_mt_over_wt = Value)
}, .id = "construct") %>%
  dplyr::filter(is.finite(luciferase_log2_mt_over_wt))
mtwt <- luc_mt %>% dplyr::inner_join(te_5 %>% dplyr::select(Group, delta_TE), by = "Group")


## =============================================================================
## 3b. the luciferase measurements, predicted with these coefficients
## =============================================================================
cat("\n[Luciferase WT]\n")
luc_raw <- readr::read_csv(LUC_WT_FILE, show_col_types = FALSE)
names(luc_raw)[1] <- "label"
rep_cols <- grep("^repeat", names(luc_raw), value = TRUE)
lucW <- luc_raw %>%
  dplyr::mutate(Group = stringr::str_extract(label, "G[0-9]+"),
                utr = ifelse(grepl("5'", label), "5'", ifelse(grepl("3'", label), "3'", NA)),
                luc_mean = rowMeans(dplyr::across(dplyr::all_of(rep_cols)), na.rm = TRUE),
                luc_log2 = log2(luc_mean),
                Gene = unname(GENE_OF[Group]),
                Gene = ifelse(is.na(Gene), Group, Gene)) %>%
  dplyr::filter(!is.na(Group), is.finite(luc_mean), utr == "5'")
pred_set <- lucW %>% dplyr::select(Group) %>%
  dplyr::inner_join(feat[["5'"]], by = "Group")
lucWT <- lucW %>%
  dplyr::select(label, Group, Gene, dplyr::all_of(rep_cols), luc_mean, luc_log2) %>%
  dplyr::left_join(tibble::tibble(Group = pred_set$Group, predicted_TE = predict_TE(pred_set)),
                   by = "Group") %>%
  dplyr::left_join(d %>% dplyr::select(Group, observed_TE = y), by = "Group") %>%
  dplyr::arrange(Group)
lw <- lucWT %>% dplyr::filter(is.finite(predicted_TE))
r3  <- cor(lw$predicted_TE, lw$luc_log2); rho3 <- cor(lw$predicted_TE, lw$luc_log2, method = "spearman")
p3_pearson  <- cor.test(lw$predicted_TE, lw$luc_log2)$p.value
p3_spearman <- suppressWarnings(cor.test(lw$predicted_TE, lw$luc_log2, method = "spearman"))$p.value
cat(sprintf("  %d constructs | Pearson %.3f, Spearman %.3f (P = %.3g)\n",
            nrow(lw), r3, rho3, p3_spearman))

cat("\n[Luciferase mutant/WT]\n")

stats_of <- function(df, x) df %>% dplyr::summarise(
  n = dplyr::n(),
  pearson = cor(.data[[x]], luciferase_log2_mt_over_wt),
  pearson_p = cor.test(.data[[x]], luciferase_log2_mt_over_wt)$p.value,
  spearman = cor(.data[[x]], luciferase_log2_mt_over_wt, method = "spearman"),
  spearman_p = suppressWarnings(cor.test(.data[[x]], luciferase_log2_mt_over_wt,
                                         method = "spearman"))$p.value,
  .groups = "drop")
mt_stats <- mtwt %>% dplyr::group_by(construct) %>% stats_of("delta_TE")
print(as.data.frame(mt_stats %>% dplyr::mutate(dplyr::across(c(pearson, spearman), ~round(.x, 3)))),
      row.names = FALSE)
mt_pool <- stats_of(mtwt, "delta_TE")
luc_stats <- tibble::tibble(n = nrow(lw), spearman = rho3, spearman_p = p3_spearman, pearson = r3)


## =============================================================================
## 4. the figures
## =============================================================================


## =============================================================================
## 6. figures  --  journal style
## =============================================================================
# Sized to real journal column widths (single 89 mm, double 183 mm) with 7-8 pt
# type, hairline axes and no grid, so the files can go into a manuscript as they
# are. They are vector (pdf, svg) and 600 dpi (png), so they also scale up
# cleanly onto a slide.

MM <- 1/25.4                                   # millimetres to inches

# the figure font, chosen once and used by every figure below
FONT <- "sans"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  fams <- unique(systemfonts::system_fonts()$family)
  for (ff in c("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"))
    if (ff %in% fams) { FONT <- ff; break }
}

BASE <- 8                                      # base font size, pt
AXIS_TITLE_PT <- BASE + 1.5                    # axis titles, which carry long labels
# and read small at the base size
cat("\n[Figures] ->", FIG_DIR, "| font:", FONT, "\n")

# palette: low saturation, harmonised, legible in greyscale and for the common
# forms of colour blindness (blue against terracotta rather than red/green)
C_UP   <- "#2F5D8A"   # steel blue    - increases TE
C_DN   <- "#C2593F"   # terracotta    - decreases TE
C_WT   <- "#1F7A6D"   # deep teal     - wild-type luciferase
C_LIB  <- "#2F5D8A"   # library constructs
C_FL   <- "#8A5A9E"   # full-length constructs
INK    <- "#1F1F1F"; AXIS <- "#333333"; SOFT <- "#5A5A5A"
LINE   <- "#1F1F1F"; BAND <- "#C9D1DC"; IDL <- "#A0A0A0"; QUAD <- "#EEF2F6"

theme_pub <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line  = element_line(linewidth = .3, colour = AXIS),
          axis.ticks = element_line(linewidth = .3, colour = AXIS),
          axis.ticks.length = unit(1.1, "mm"),
          axis.text  = element_text(size = base - .5, colour = AXIS),
          axis.title = element_text(size = AXIS_TITLE_PT, colour = INK),
          axis.title.x = element_text(margin = margin(t = 4)),
          axis.title.y = element_text(margin = margin(r = 4)),
          plot.title = element_text(size = base - .5, colour = SOFT, hjust = 0,
                                    margin = margin(b = 5)),
          plot.title.position = "plot",
          legend.background = element_blank(), legend.key = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = base - .5, colour = INK),
          legend.key.size = unit(3, "mm"),
          plot.margin = margin(5, 6, 4, 4))
}
TXT <- (BASE - .5)/.pt                          # geom text size matching axis text

# P values written as a journal would: three decimals down to 0.001, then
# mantissa x 10^exponent
p_expr <- function(p) {
  if (!is.finite(p)) return(quote(italic(P) == "NA"))
  if (p <= .Machine$double.xmin) return(quote(italic(P) < 10^-300))
  if (p >= 1e-3) return(bquote(italic(P) == .(sprintf("%.3f", p))))
  e <- floor(log10(p)); m <- p/10^e
  if (round(m, 1) >= 10) { m <- m/10; e <- e + 1 }
  bquote(italic(P) == .(sprintf("%.1f", m)) %*% 10^.(e))
}
title_spearman <- function(rho, p)
  bquote("Spearman" ~ rho == .(sprintf("%.3f", rho)) * "," ~~ .(p_expr(p)))
TE_lab <- function(prefix) bquote(.(prefix) ~ log[2] * "(Poly / (Free + 40S))")

# feature labels: k-mers named by their length in the RNA alphabet (3-mer CGG),
# like the RBP motifs beside them


save_fig <- function(p, name, w_mm, h_mm) {
  base <- file.path(FIG_DIR, name); w <- w_mm*MM; h <- h_mm*MM
  ggsave(paste0(base, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 600, bg = "white")
  svg_ok <- FALSE
  if (requireNamespace("svglite", quietly = TRUE)) {
    svg_ok <- tryCatch({
      ggsave(paste0(base, ".svg"), p, width = w, height = h, device = svglite::svglite)
      TRUE
    }, error = function(e) {
      if (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE)
      message("  svglite failed (", conditionMessage(e), "); using grDevices::svg")
      FALSE
    })
  }
  if (!svg_ok)
    ggsave(paste0(base, ".svg"), p, width = w, height = h, device = grDevices::svg)
  cat(sprintf("  %-36s %3.0f x %3.0f mm  (.pdf .png .svg)\n", name, w_mm, h_mm))
}

## ---- coefficients, grouped by feature class ------------------------------------
# Features are grouped by what they describe, so the figure says what KIND of
# sequence information the model uses. Each class is marked by a colour ribbon
# with its name on the left and a faint tint behind its rows; class colours are
# earth tones so they never compete with the blue/terracotta that encodes the
# direction of the effect.
# One coefficient (the CG k-mer) is several times larger than the rest and would
# compress every other bar into a sliver, so the axis is cut just past the
# second-largest effect and the cut bar carries a break mark.
feat_class <- function(id) dplyr::case_when(
  grepl("^kmer_", id)                    ~ "k-mer",
  grepl("^RBP_",  id)                    ~ "RBP motif",
  grepl("^ARE_",  id)                    ~ "ARE",
  grepl("^miRNA_|^miR_", id)             ~ "miRNA site",
  grepl("^uAUG_", id)                    ~ "Upstream AUG",
  grepl("^RG4|stemLoop|struc|psudoK|pseudoK|KISS|freeEnergy", id) ~ "Structure",
  grepl("way$", id)                      ~ "Conservation",
  TRUE                                   ~ "Composition")
CLASS_LEVELS <- c("Composition", "k-mer", "RBP motif", "ARE", "miRNA site",
                  "Upstream AUG", "Structure", "Conservation")
CLASS_COL  <- c(Composition = "#9A7F55", `k-mer` = "#5E8571",
                `RBP motif` = "#6E6A96", ARE = "#9A6653",
                `miRNA site` = "#4E7B9C", `Upstream AUG` = "#B08347",
                Structure = "#7D8B6A", Conservation = "#8A6E8F")
CLASS_TINT <- c(Composition = "#F6F1E8", `k-mer` = "#EEF4F0",
                `RBP motif` = "#F1F0F6", ARE = "#F7EFEC",
                `miRNA site` = "#EDF2F6", `Upstream AUG` = "#F8F2E9",
                Structure = "#F1F3EE", Conservation = "#F4F0F5")

prep_coef <- function(tbl, n_top) {
  tbl %>% dplyr::arrange(dplyr::desc(abs(coef))) %>%
    { if (is.finite(n_top)) dplyr::slice_head(., n = n_top) else . } %>%
    dplyr::mutate(label = pretty_feat(label, feature_id),
                  label = ifelse(duplicated(label) | duplicated(label, fromLast = TRUE),
                                 paste0(label, " [", feature_id, "]"), label),
                  class = factor(feat_class(feature_id), levels = CLASS_LEVELS),
                  direction = factor(direction, levels = c("increases", "decreases"))) %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::arrange(class, coef) %>%                   # within a class, by effect
    dplyr::mutate(label = factor(label, levels = unique(label)))
}

# Widths are set in millimetres from the content itself -- the longest feature
# name and the longest class name -- so no label can be clipped however many rows
# the figure has. Text width is estimated from the font size: about 0.55 em per
# character for Arial, with a margin.
char_mm <- function(pt, bold = FALSE) pt*0.3528*(if (bold) .62 else .56)

# The shape of the coefficient figure. A tall row height with narrow panels
# gives a figure close to square, which sits beside the other panels of a
# composite figure; the wide, short default reads better on its own.
# Two shapes are written for each coefficient figure: the near-square one, which
# sits beside another panel in a composite figure, and the long one, which reads
# better on its own. Each is a row height and three widths, in mm.
COEF_SHAPES <- list(
  square = list(row = 4.6, bar = 46, freq = 13, gap = 4, suffix = ""),
  wide   = list(row = 3.3, bar = 62, freq = 20, gap = 9, suffix = "_wide"))

COEF_ROW_MM  <- 4.6      # height of one feature row, mm
COEF_BAR_MM  <- 46       # width of the coefficient panel, mm
COEF_FREQ_MM <- 13       # width of the selection-probability panel, mm
COEF_GAP_MM  <- 4        # space between the two, mm
COEF_BAND    <- 0.52     # thickness of an interval, as a share of the row
COEF_DOT     <- 0.40     # size of an estimate, as a share of the row
COEF_LEGEND  <- FALSE    # the colours are explained by the axis, so the legend
# under the figure is off by default

make_coef_fig <- function(tbl, with_freq, row_mm = COEF_ROW_MM, txt = BASE - .5,
                          bar_mm = COEF_BAR_MM, freq_mm = COEF_FREQ_MM,
                          gap_mm = COEF_GAP_MM) {
  lo <- min(0, min(tbl$ci_low)); hi <- max(0, max(tbl$ci_high))
  xlim <- c(lo - .04*(hi - lo), hi + .06*(hi - lo))
  
  present <- levels(droplevels(tbl$class))
  facet <- facet_grid(class ~ ., scales = "free_y", space = "free_y")
  no_strip <- theme(strip.background = element_blank(), strip.text = element_blank(),
                    panel.spacing.y = unit(1.8, "mm"))
  tints <- lapply(present, function(k)
    geom_rect(data = data.frame(class = factor(k, levels = CLASS_LEVELS)),
              aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = CLASS_TINT[[k]]))
  
  ## ---- the class ribbon, sized to its longest name ----
  # The class names are written the same way for every class: along the ribbon
  # when every one of them fits, beside it otherwise. Deciding this once rather
  # than class by class keeps a short class from being turned while its
  # neighbour lies flat.
  rib_lab <- tbl %>% dplyr::group_by(class) %>%
    dplyr::summarise(n = dplyr::n(),
                     mid = levels(droplevels(label))[ceiling(dplyr::n()/2)],
                     .groups = "drop") %>%
    dplyr::mutate(name = as.character(class),
                  len_mm = nchar(name)*char_mm(CLASS_PT, bold = TRUE),
                  nudge = ifelse(n %% 2 == 0, .5, 0))
  rib_lab$vertical <- all(rib_lab$n*row_mm >= rib_lab$len_mm + 2)
  TILE <- 1.6                                             # ribbon width, mm
  horiz_mm <- if (!rib_lab$vertical[1]) max(rib_lab$len_mm) else 0
  rib_mm <- max(8, horiz_mm + TILE + 3.5)
  ribbon <- ggplot(tbl, aes(y = label)) + facet +
    lapply(present, function(k)
      geom_tile(data = tbl %>% dplyr::filter(class == k),
                aes(x = rib_mm - TILE/2, y = label),
                width = TILE, height = 1, fill = CLASS_COL[[k]])) +
    lapply(seq_len(nrow(rib_lab)), function(i) {
      r <- rib_lab[i, ]
      geom_text(data = r, aes(y = mid), inherit.aes = FALSE,
                x = if (r$vertical) rib_mm - TILE - 2.4 else rib_mm - TILE - 1.4,
                label = r$name, angle = if (r$vertical) 90 else 0,
                hjust = if (r$vertical) .5 else 1, vjust = .5, nudge_y = r$nudge,
                size = CLASS_PT/.pt, family = FONT, fontface = "bold",
                colour = CLASS_COL[[r$name]])
    }) +
    scale_x_continuous(limits = c(0, rib_mm), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = FONT) + no_strip +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  ## ---- the coefficients: estimate and 95% interval ----
  # The interval is drawn as a band that fades outwards rather than as a whisker:
  # each interval is cut into NSLICE pieces whose opacity falls with the distance
  # from the estimate, so the eye is drawn to the middle of the interval.
  NSLICE <- 48
  slices <- purrr::map_dfr(seq_len(nrow(tbl)), function(i) {
    r <- tbl[i, ]
    e <- seq(r$ci_low, r$ci_high, length.out = NSLICE + 1)
    mid <- (utils::head(e, -1) + e[-1])/2
    half <- max((r$ci_high - r$ci_low)/2, 1e-12)
    tibble::tibble(label = r$label, class = r$class, direction = r$direction,
                   x0 = utils::head(e, -1), x1 = e[-1],
                   a = (1 - abs(mid - r$coef)/half)^1.6)
  })
  bar <- ggplot(tbl, aes(coef, label)) + facet + tints +
    geom_vline(xintercept = 0, colour = AXIS, linewidth = .3) +
    geom_segment(data = slices, aes(x = x0, xend = x1, y = label, yend = label,
                                    colour = direction, alpha = a),
                 inherit.aes = FALSE, linewidth = row_mm*COEF_BAND, lineend = "butt") +
    geom_point(aes(fill = direction), shape = 21, size = row_mm*COEF_DOT,
               colour = "white", stroke = .45) +
    scale_colour_manual(values = c(increases = C_UP, decreases = C_DN), guide = "none") +
    scale_alpha_continuous(range = c(.06, .55), guide = "none") +
    scale_fill_manual(values = c(increases = C_UP, decreases = C_DN),
                      labels = c(increases = "Increases TE", decreases = "Decreases TE")) +
    scale_x_continuous(breaks = scales::breaks_pretty(4), expand = c(0, 0)) +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(x = "Coefficient (standardised)  ·  95% interval", y = NULL) +
    theme_pub() + no_strip +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(colour = INK, size = txt, margin = margin(r = 2)),
          axis.title.x = element_text(margin = margin(t = 1.5)),
          legend.position = if (COEF_LEGEND) "bottom" else "none",
          legend.justification = "center", legend.margin = margin(2, 0, 0, 0))
  
  ## ---- absolute widths, so the names always have room ----
  lab_mm <- max(nchar(as.character(tbl$label)))*char_mm(txt) + 3
  panels <- list(ribbon, bar); w <- c(rib_mm, bar_mm)
  if (with_freq) {
    # from zero, so the length of each bar is the selection probability itself
    # the value is written on the bar: inside it when the bar is long enough to
    # hold the text, just outside it when it is not
    freq_lab <- tbl %>%
      dplyr::mutate(txt = sprintf("%.0f%%", 100*sel_freq),
                    inside = sel_freq > .28,
                    x = ifelse(inside, sel_freq - .02, sel_freq + .02),
                    hj = ifelse(inside, 1, 0),
                    col = ifelse(inside, "white", SOFT))
    freq <- ggplot(tbl, aes(sel_freq, label)) + facet + tints +
      geom_col(width = .68, fill = "#8E9AA8") +
      geom_text(data = freq_lab, aes(x = x, y = label, label = txt, hjust = hj,
                                     colour = col),
                inherit.aes = FALSE, size = (txt - 1.2)/.pt, family = FONT,
                fontface = "bold", show.legend = FALSE) +
      scale_colour_identity() +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1),
                         labels = scales::percent_format(accuracy = 1),
                         expand = expansion(mult = c(0, .02))) +
      labs(x = "Selection\nprobability", y = NULL) +
      theme_pub() + no_strip +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_blank())
    # an explicit gap, so the tick labels of the two x-axes cannot meet
    panels <- c(panels, list(plot_spacer(), freq)); w <- c(w, gap_mm, freq_mm)
  }
  p <- wrap_plots(panels, nrow = 1, widths = unit(w, "mm")) &
    theme(plot.margin = margin(4, 3, 3, 2))
  list(plot = p,
       width_mm  = sum(w) + lab_mm + (length(w) - 1)*3 + 14,
       height_mm = nrow(tbl)*row_mm + (length(present) - 1)*1.8 +
         if (COEF_LEGEND) 30 else 22)
}

# Every feature of the reported model is drawn. A second, shorter figure showing
# only the largest coefficients is made only when there are more features than
# fit comfortably in one panel; below that the two would be the same figure.
cf_all <- prep_coef(fig_cf, Inf)
for (with_freq in c(TRUE, FALSE)) {
  g <- make_coef_fig(cf_all, with_freq)
  save_fig(g$plot, paste0("coefficients", if (with_freq) "_with_probability" else ""),
           g$width_mm, g$height_mm)
}
if (nrow(cf_all) > TOP_N_FEAT) {
  cf_top <- prep_coef(fig_cf, TOP_N_FEAT)
  save_coef(cf_top, "coefficients_top30")
  cat(sprintf("  (%d features in all; coefficients_top30 shows the %d largest)\n",
              nrow(cf_all), nrow(cf_top)))
} else {
  cat(sprintf("  (%d features, so one figure only)\n", nrow(cf_all)))
}

## ---- sequence logos of the RBP motifs in the model --------------------------
# Each RBP feature is a count of matches to an ATtRACT position weight matrix, so
# a logo shows what sequence the feature actually scores. PWMs come from
# human_RBP_pwm.rdata (RBP.matrix), as in the original stability analysis.
# Feature ids and PWM names are matched after turning every non-alphanumeric
# character into "_" on both sides and dropping a "_WT" suffix: the ids went
# through make.names() and the PWM names contain "/", "-" and ".".
suppressPackageStartupMessages(library(ggseqlogo))
RBP_PWM_FILE <- file.path(ROOT, "human_RBP_pwm.rdata")
if (file.exists(RBP_PWM_FILE)) {
  pe <- new.env(); load(RBP_PWM_FILE, envir = pe)
  pwm <- get("RBP.matrix", envir = pe); rm(pe)
  key <- function(x) gsub("[^A-Za-z0-9]+", "_", sub("_WT$", "", x))
  names(pwm) <- key(paste0("RBP_", names(pwm)))
  # RNA alphabet, and each position as probabilities so bits are computed properly
  pwm <- lapply(pwm, function(m) {
    m <- as.matrix(m)
    if (nrow(m) != 4 && ncol(m) == 4) m <- t(m)
    if (nrow(m) != 4) return(NULL)
    rownames(m) <- c("A", "C", "G", "U")
    sweep(m, 2, pmax(colSums(m), 1e-12), "/")
  })
  pwm <- Filter(Negate(is.null), pwm)
  
  NT_COL <- make_col_scheme(chars = c("A", "C", "G", "U"),
                            cols = c("#3A8E5C", "#2F5D8A", "#D08C2E", "#C2593F"))
  
  make_logo_fig <- function(tbl, name, ncol = 4) {
    rb <- tbl %>% dplyr::filter(grepl("^RBP_", feature_id)) %>%
      dplyr::arrange(dplyr::desc(abs(coef))) %>%
      dplyr::mutate(k = key(feature_id))
    miss <- rb$feature_id[!rb$k %in% names(pwm)]
    if (length(miss))
      cat("  no PWM found for:", paste(miss, collapse = ", "), "\n")
    rb <- rb %>% dplyr::filter(k %in% names(pwm))
    if (!nrow(rb)) { cat("  no RBP features with a PWM; logos skipped\n"); return(invisible()) }
    # panel titles: motif (protein), then the coefficient with its direction
    # labels were already made readable in prep_coef(); use them as they are
    ttl <- sprintf("%s\n\u03b2 = %+.3f", as.character(rb$label), rb$coef)
    ttl <- make.unique(ttl, sep = " ")
    mats <- setNames(pwm[rb$k], ttl)
    nr <- ceiling(length(mats)/ncol)
    # ggseqlogo still calls guides(... = FALSE) and aes_string(), which newer
    # ggplot2 flags as deprecated; the notices come from the package, not this code
    p <- suppressWarnings(ggseqlogo(mats, method = "bits", col_scheme = NT_COL, ncol = ncol)) +
      labs(y = "Bits") +
      theme_classic(base_size = BASE, base_family = FONT) +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = BASE - .5, colour = INK, lineheight = 1.05),
            axis.line = element_line(linewidth = .3, colour = AXIS),
            axis.ticks = element_line(linewidth = .3, colour = AXIS),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = BASE - 1.5, colour = AXIS),
            axis.title.y = element_text(size = BASE - .5, colour = INK),
            legend.position = "none",
            panel.spacing = unit(3, "mm"),
            plot.margin = margin(4, 6, 4, 4))
    withCallingHandlers(
      save_fig(p, name, 42*min(ncol, length(mats)) + 8, 26*nr + 8),
      warning = function(w) if (grepl("deprecated|aes_string|guides", conditionMessage(w)))
        invokeRestart("muffleWarning"))
  }
  make_logo_fig(cf_all, "RBP_motif_logos")
} else {
  cat("  human_RBP_pwm.rdata not found at", RBP_PWM_FILE, "- motif logos skipped\n")
}


## ---- observed against predicted, coloured by point density ---------------------
r2 <- cor(fig_best$observed, fig_best$predicted)
rho2 <- cor(fig_best$observed, fig_best$predicted, method = "spearman")
p2_pearson  <- cor.test(fig_best$observed, fig_best$predicted)$p.value
p2_spearman <- suppressWarnings(cor.test(fig_best$observed, fig_best$predicted,
                                         method = "spearman"))$p.value
# Point density as in ggpointdensity: every sequence is drawn, coloured by how
# many OTHER sequences lie within a small circle around it. Both axes are scaled
# by their own range first, so the circle is round on the plot; its radius is
# NEIGHBOUR_R of that range. The count is computed here rather than taken from a
# package, so the definition is explicit and does not change between versions.
NEIGHBOUR_R <- 0.05
n_fig2 <- nrow(fig_best)
leg_theme <- theme(legend.position = "right",
                   legend.title = element_text(size = BASE - .5, colour = INK,
                                               margin = margin(b = 6), lineheight = 1.1),
                   legend.text = element_text(size = BASE - 1, colour = SOFT),
                   legend.margin = margin(0, 0, 0, 3))
bar <- guide_colourbar(barwidth = unit(2.4, "mm"), barheight = unit(20, "mm"),
                       ticks.colour = "white", title.position = "top")
# sample size above the word Density
leg_title <- bquote(atop(italic(n) == .(format(n_fig2, big.mark = ",")), "Density"))

sx <- (fig_best$observed  - min(fig_best$observed))/diff(range(fig_best$observed))
sy <- (fig_best$predicted - min(fig_best$predicted))/diff(range(fig_best$predicted))
d2 <- outer(sx, sx, "-")^2 + outer(sy, sy, "-")^2
fb <- fig_best %>%
  dplyr::mutate(neighbours = rowSums(d2 <= NEIGHBOUR_R^2) - 1L) %>%   # not itself
  dplyr::arrange(neighbours)                 # densest drawn last, on top
cat(sprintf("  point density: neighbours within %.0f%% of the axis range, %d to %d\n",
            100*NEIGHBOUR_R, min(fb$neighbours), max(fb$neighbours)))
p2 <- ggplot(fb, aes(predicted, observed)) +
  geom_point(aes(colour = neighbours), size = 1.1, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black",
              fill = "grey70", alpha = .45, linewidth = .5) +
  scale_colour_viridis_c(option = "magma", end = .97, name = leg_title,
                         breaks = scales::breaks_pretty(3), guide = bar) +
  labs(title = title_spearman(rho2, p2_spearman),
       x = TE_lab("Predicted"), y = TE_lab("Observed")) +
  theme_pub() + leg_theme
save_fig(p2, "observed_vs_predicted_heldout", 98, 80)   # wider by the colour bar

TE_X      <- bquote(Delta * log[2] * "(Poly / (Free + 40S))")

mk_mtwt <- function(df, xvar, xlab, stats, subtitle = NULL,
                    colour_by = NULL, one_colour = C_LIB, labs_ = NULL) {
  p <- ggplot(df, aes(.data[[xvar]], luciferase_log2_mt_over_wt)) +
    # the two quadrants where the MPRA and the luciferase agree on the direction
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = QUAD) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = QUAD) +
    geom_hline(yintercept = 0, colour = "#B8BEC6", linewidth = .3) +
    geom_vline(xintercept = 0, colour = "#B8BEC6", linewidth = .3) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = LINE,
                fill = BAND, alpha = .45, linewidth = .5)
  p <- p + if (is.null(colour_by))
    geom_point(shape = 21, size = 2.6, fill = one_colour, colour = "white", stroke = .45)
  else
    geom_point(aes(fill = .data[[colour_by]]), shape = 21, size = 2.6,
               colour = "white", stroke = .45, alpha = .72)
  p <- p +
    ggrepel::geom_text_repel(aes(label = Gene), size = TXT, family = FONT,
                             fontface = "italic", colour = INK, box.padding = .45,
                             point.padding = .3, min.segment.length = Inf, seed = 41,
                             max.overlaps = Inf) +
    labs(title = title_spearman(stats$spearman, stats$spearman_p), subtitle = subtitle,
         x = xlab, y = expression(log[2] ~ "translation efficiency (mutant / WT)")) +
    theme_pub() + theme(plot.subtitle = element_text(size = BASE - 1, colour = SOFT))
  if (is.null(colour_by)) return(p)
  p <- p + scale_fill_manual(values = c(Library = C_LIB, `Full-length` = C_FL),
                             labels = if (is.null(labs_)) ggplot2::waiver() else labs_)
  if (is.null(labs_)) p + guides(fill = "none")
  else p + theme(legend.position = "top", legend.justification = "right",
                 legend.direction = "vertical", legend.title = element_blank(),
                 legend.background = element_blank(), legend.key = element_blank(),
                 legend.margin = margin(0, 0, 1, 0), legend.spacing.y = unit(0, "mm"),
                 legend.key.height = unit(3.0, "mm"), legend.key.width = unit(3.0, "mm"),
                 legend.text = element_text(size = BASE - 1, colour = INK),
                 plot.title = element_text(margin = margin(b = 0))) +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 2.4, alpha = 1)))
}

## ---- wild-type luciferase, with the two replicates as a range ----------------------
lw <- lw %>% dplyr::mutate(
  rep_lo = log2(pmin(.data[[rep_cols[1]]], .data[[rep_cols[length(rep_cols)]]])),
  rep_hi = log2(pmax(.data[[rep_cols[1]]], .data[[rep_cols[length(rep_cols)]]])))
p3 <- ggplot(lw, aes(predicted_TE, luc_log2)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = LINE,
              fill = BAND, alpha = .45, linewidth = .5) +
  geom_errorbar(aes(ymin = rep_lo, ymax = rep_hi), width = 0, linewidth = .35,
                colour = C_WT) +
  geom_point(shape = 21, size = 2.6, fill = C_WT, colour = "white", stroke = .45) +
  ggrepel::geom_text_repel(aes(label = Gene), size = TXT, family = FONT,
                           fontface = "italic", colour = INK, box.padding = .45,
                           point.padding = .3, min.segment.length = Inf, seed = 41) +
  labs(title = title_spearman(rho3, p3_spearman), x = TE_lab("Predicted"),
       y = expression(log[2] ~ "luciferase activity")) +
  theme_pub()
save_fig(p3, "luciferase_WT_predicted", 84, 80)

## ---- mutant versus wild-type, sign-concordant quadrants shaded ---------------------
n_by   <- mtwt %>% dplyr::count(construct)
lab_by <- setNames(sprintf("%s (n = %d)", n_by$construct, n_by$n), n_by$construct)
p4c <- mk_mtwt(mtwt, "delta_TE", TE_X, mt_pool, colour_by = "construct", labs_ = lab_by)
save_fig(p4c, "luciferase_mtWT_measured_combined", 86, 88)

## =============================================================================
## 5. the model object
## =============================================================================
# Only the model is written, and only locally: it is what export_figure_data.R
# reads to build the data behind Figures 7B, 7D and 7E. No table here is meant
# for the repository.
model <- list(
  description = sprintf("%s; 2:1 split %s; stability selection (%d half-samples, weighted cv.glmnet at %s); features above %s refitted by weighted linear regression",
                        TAG,
                        if (file.exists(SPLIT)) sprintf("supplied in %s", basename(SPLIT))
                        else sprintf("drawn with seed %d", SEED),
                        NUM_SS, LAMBDA_RULE, THR_TXT),
  intercept = intercept, coefficients_all = beta, stable_features = stable,
  table = cf, X = X, y = y, weights = w,
  modelling_table = d %>% dplyr::select(Group, seqName, y, Mbar),
  train_rows = tr_rows, test_rows = te_rows,
  selection_matrix = sel_mat, selection_probability = sel_freq,
  clusters = vs$membership,
  # the points behind the three panels, which export_Fig7_data.R turns into csv
  held_out = best, luciferase_WT = lw, luciferase_mtWT = mtwt,
  settings = list(MIN_COUNT = MIN_COUNT, MIN_BATCH = MIN_BATCH, WEIGHT_POW = WEIGHT_POW,
                  MAX_NA = MAX_NA, MIN_NONZERO = MIN_NONZERO, VIF_THRESH = VIF_THRESH,
                  NUM_SS = NUM_SS, NFOLD = NFOLD, LAMBDA_RULE = LAMBDA_RULE,
                  STABLE_THRESHOLD = THR, TEST_FRACTION = TEST_FRACTION, SEED = SEED))
saveRDS(model, file.path(OUT_DIR, "model.rds"))
cat(sprintf("\n[Written] %s\n", file.path(OUT_DIR, "model.rds")))
cat(sprintf("          %s\n", file.path(FIG_DIR, "coefficients(_with_probability)")))
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
cat("\n=== Done ===\n")
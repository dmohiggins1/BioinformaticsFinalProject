suppressPackageStartupMessages({
  library(edgeR)
})

# ---------- Settings ----------
BCV <- 0.2
DISP <- BCV^2  # 0.04

infile <- "gene_count.csv"   # must have column gene_id, then PU1 PT1 RU1 RT1
outdir <- "edgeR_pairwise_BCV02"
dir.create(outdir, showWarnings = FALSE)

# ---------- Load counts ----------
counts <- read.csv(infile, row.names = 1, check.names = FALSE)
needed <- c("PU1", "PT1", "RU1", "RT1")
stopifnot(all(needed %in% colnames(counts)))
counts <- counts[, needed]

# ---------- Helper: run a 2-sample edgeR GLM with fixed dispersion ----------
run_edger_2sample <- function(counts2, control_name, treat_name, disp, outpath) {
  # counts2: matrix/data.frame with exactly 2 columns: control, treatment
  stopifnot(ncol(counts2) == 2)
  
  group <- factor(c("Control", "Treatment"), levels = c("Control", "Treatment"))
  y <- DGEList(counts = counts2, group = group)
  
  # Filter (with only 2 samples, keep it simple and conservative)
  keep <- filterByExpr(y, group = group, min.count = 10)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~ group)
  
  fit <- glmFit(y, design, dispersion = disp)
  lrt <- glmLRT(fit, coef = "groupTreatment")
  
  res <- topTags(lrt, n = Inf)$table
  res$FDR <- p.adjust(res$PValue, method = "BH")
  
  # Add labels
  res$control <- control_name
  res$treatment <- treat_name
  
  write.table(res, file = outpath, sep = "\t", quote = FALSE)
}

# ---------- Pairwise comparisons ----------
comparisons <- list(
  c("PU1", "PT1"),
  c("PU1", "RU1"),
  c("PU1", "RT1"),
  c("RT1", "RU1")  # RU1 vs RT1 (treatment = RU1, control = RT1)
)

for (cmp in comparisons) {
  control <- cmp[1]
  treat   <- cmp[2]
  
  counts2 <- counts[, c(control, treat)]
  outpath <- file.path(outdir, paste0("edgeR_", treat, "_vs_", control, "_BCV02.tsv"))
  
  message("Running edgeR: ", treat, " vs ", control)
  run_edger_2sample(counts2, control_name = control, treat_name = treat, disp = DISP, outpath = outpath)
}

message("Done. Results in: ", outdir)

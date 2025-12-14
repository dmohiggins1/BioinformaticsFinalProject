suppressPackageStartupMessages({
  library(DESeq2)
})

# ---- Settings ----
infile <- "gene_count.csv"
outdir <- "DESeq2_fixedDisp_BCV02"
dir.create(outdir, showWarnings = FALSE)

BCV <- 0.2
disp_assumed <- BCV^2          # 0.04
alpha_assumed <- 0.05          # FDR threshold used in results()

# ---- Load counts ----
counts <- read.csv(infile, row.names = 1, check.names = FALSE)
needed <- c("PU1", "PT1", "RU1", "RT1")
stopifnot(all(needed %in% colnames(counts)))
counts <- counts[, needed]

# ---- Helper: run 1v1 with fixed dispersion ----
run_deseq2_fixed_disp_1v1 <- function(counts2, control_name, treat_name,
                                      disp_fixed, alpha, outpath) {
  stopifnot(ncol(counts2) == 2)
  
  coldata <- data.frame(
    row.names = colnames(counts2),
    condition = factor(c("Control", "Treatment"), levels = c("Control", "Treatment"))
  )
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts2)),
    colData = coldata,
    design = ~ condition
  )
  
  # Filter very low counts (optional but recommended)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  
  # Normalize
  dds <- estimateSizeFactors(dds)
  
  # Force a fixed dispersion for every gene
  dispersions(dds) <- rep(disp_fixed, nrow(dds))
  
  # Wald test using the fixed dispersions
  dds <- suppressWarnings(nbinomWaldTest(dds, betaPrior = FALSE))
  
  
  res <- results(
    dds,
    contrast = c("condition", "Treatment", "Control"),
    alpha = alpha,
    cooksCutoff = FALSE,          # <- prevents qf() NaN warnings in 1v1
    independentFiltering = FALSE  # <- optional; avoids filtering logic based on alpha
  )
  
  
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$control <- control_name
  res_df$treatment <- treat_name
  res_df$dispersion_assumed <- disp_fixed
  res_df$alpha_assumed <- alpha
  
  # Sort by padj then effect
  res_df <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange)), ]
  
  write.table(res_df, file = outpath, sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---- Comparisons you requested ----
comparisons <- list(
  c("PU1", "PT1"),
  c("PU1", "RU1"),
  c("PU1", "RT1"),
  c("RT1", "RU1")  # RU1 vs RT1 (control = RT1)
)

for (cmp in comparisons) {
  control <- cmp[1]
  treat <- cmp[2]
  counts2 <- counts[, c(control, treat)]
  
  outpath <- file.path(outdir, paste0("DESeq2_", treat, "_vs_", control,
                                      "_fixedDisp", disp_assumed, "_alpha", alpha_assumed, ".tsv"))
  
  message("Running DESeq2 fixed-disp: ", treat, " vs ", control)
  run_deseq2_fixed_disp_1v1(counts2, control, treat, disp_assumed, alpha_assumed, outpath)
}

message("Done. Results in: ", outdir)

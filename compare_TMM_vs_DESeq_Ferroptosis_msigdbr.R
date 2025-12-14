suppressPackageStartupMessages({
  library(edgeR)
  library(DESeq2)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(fgsea)
  library(msigdbr)
})

# ---------- Inputs ----------
counts <- read.csv("gene_count.csv", row.names = 1, check.names = FALSE)
colnames(counts) <- trimws(colnames(counts))
counts <- counts[, c("PU1","PT1","RU1","RT1")]

contrasts <- list(
  PT1_vs_PU1 = c("PU1","PT1"),
  RU1_vs_PU1 = c("PU1","RU1"),
  RT1_vs_PU1 = c("PU1","RT1"),
  RU1_vs_RT1 = c("RT1","RU1")
)

dir.create("NormCompare_Ferroptosis", showWarnings = FALSE)

# ---------- 1) Normalization ----------
# TMM -> CPM
y <- DGEList(counts = counts)
y <- calcNormFactors(y, method = "TMM")
tmm_cpm <- cpm(y, normalized.lib.sizes = TRUE, log = FALSE)

# DESeq -> size-factor normalized counts
dds0 <- DESeqDataSetFromMatrix(round(as.matrix(counts)),
                               colData = data.frame(row.names = colnames(counts)),
                               design = ~ 1)
dds0 <- estimateSizeFactors(dds0)
deseq_norm <- counts(dds0, normalized = TRUE)

write.csv(as.data.frame(tmm_cpm),
          "NormCompare_Ferroptosis/TMM_normalized_CPM.csv", quote = FALSE)
write.csv(as.data.frame(deseq_norm),
          "NormCompare_Ferroptosis/DESeq_normalized_counts.csv", quote = FALSE)

# ---------- 2) log2FC from normalized values ----------
log2fc_from_norm <- function(norm_mat, a, b, pseudocount) {
  log2((norm_mat[, b] + pseudocount) / (norm_mat[, a] + pseudocount))
}

# ---------- 3) Map Ensembl -> Entrez (mouse) ----------
ens <- rownames(counts)
map <- bitr(ens, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
map <- map[!duplicated(map$ENSEMBL), ]

to_entrez <- function(lfc_vec, map) {
  common <- intersect(names(lfc_vec), map$ENSEMBL)
  v <- lfc_vec[common]
  names(v) <- map$ENTREZID[match(common, map$ENSEMBL)]
  v
}

# ---------- 4) KEGG ferroptosis gene set from msigdbr ----------
ms <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
ferro <- ms[ms$gs_name == "KEGG_FERROPTOSIS", ]

if (nrow(ferro) == 0) {
  stop("msigdbr did not return KEGG_FERROPTOSIS. Check msigdbr install/version.")
}

ferro_list <- list(FERROPTOSIS = unique(ferro$entrez_gene))

# ---------- 5) fgsea per contrast, compare methods ----------
out_rows <- list()

for (nm in names(contrasts)) {
  a <- contrasts[[nm]][1]
  b <- contrasts[[nm]][2]
  
  lfc_tmm   <- log2fc_from_norm(tmm_cpm,   a, b, pseudocount = 0.5)
  lfc_deseq <- log2fc_from_norm(deseq_norm, a, b, pseudocount = 1)
  
  names(lfc_tmm) <- rownames(tmm_cpm)
  names(lfc_deseq) <- rownames(deseq_norm)
  
  r_tmm <- sort(to_entrez(lfc_tmm, map), decreasing = TRUE)
  r_ds  <- sort(to_entrez(lfc_deseq, map), decreasing = TRUE)
  
  fg_tmm <- fgsea(pathways = ferro_list, stats = r_tmm, minSize = 5, maxSize = 500)
  fg_ds  <- fgsea(pathways = ferro_list, stats = r_ds,  minSize = 5, maxSize = 500)
  
  # Gene-level table for ferroptosis genes
  ferro_ids <- ferro_list$FERROPTOSIS
  common_ids <- intersect(names(r_tmm), names(r_ds))
  ferro_common <- intersect(common_ids, ferro_ids)
  
  gene_tbl <- data.frame(
    ENTREZID = ferro_common,
    log2FC_TMM = r_tmm[ferro_common],
    log2FC_DESeq = r_ds[ferro_common]
  )
  gene_tbl <- gene_tbl[order(-abs(gene_tbl$log2FC_TMM)), ]
  
  write.csv(gene_tbl,
            file = paste0("NormCompare_Ferroptosis/FerroptosisGenes_", nm, ".csv"),
            row.names = FALSE, quote = FALSE)
  
  out_rows[[length(out_rows) + 1]] <- data.frame(
    contrast = nm,
    NES_TMM = fg_tmm$NES,
    padj_TMM = fg_tmm$padj,
    NES_DESeq = fg_ds$NES,
    padj_DESeq = fg_ds$padj,
    spearman_log2FC_ferro = if (nrow(gene_tbl) >= 3)
      cor(gene_tbl$log2FC_TMM, gene_tbl$log2FC_DESeq, method = "spearman") else NA_real_
  )
}

summary_tbl <- do.call(rbind, out_rows)
write.csv(summary_tbl,
          "NormCompare_Ferroptosis/Ferroptosis_TMM_vs_DESeq_summary.csv",
          row.names = FALSE, quote = FALSE)

summary_tbl

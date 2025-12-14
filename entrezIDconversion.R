suppressPackageStartupMessages({
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

folder <- "NormCompare_KEGG_mmu04216_Ferroptosis"
files <- list.files(folder, pattern = "^mmu04216_genes_.*\\.csv$", full.names = TRUE)

for (f in files) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  ann <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = as.character(df$ENTREZID),
    keytype = "ENTREZID",
    columns = c("SYMBOL", "GENENAME")
  )
  ann <- ann[!duplicated(ann$ENTREZID), ]
  
  df2 <- merge(df, ann, by = "ENTREZID", all.x = TRUE)
  df2 <- df2[, c("ENTREZID", "SYMBOL", "GENENAME",
                 setdiff(colnames(df2), c("ENTREZID","SYMBOL","GENENAME")))]
  
  out <- sub("\\.csv$", "_ANNOTATED.csv", f)
  write.csv(df2, out, row.names = FALSE, quote = FALSE)
}

message("Done. Wrote annotated files ending in _ANNOTATED.csv")


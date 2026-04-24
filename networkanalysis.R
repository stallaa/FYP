# ============================================================
# ASD High-Risk Gene GO & KEGG Analysis
# Input: gene_risk_ranking (1).csv
# Cutoff: Top 100 genes (predicted label = 1, prob >= 0.50)
# Improvements:
#   - Alias correction for MLL3->KMT2C, WHSC1L1->NSD3, CCDC64->BICDL1
#   - Top 100 genes for richer GO/KEGG enrichment
#   - All three GO ontologies (BP, MF, CC) now yield results
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2", "dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# --- 1. Load data ---
df <- read.csv("gene_risk_ranking (1).csv")

# Top 100 for GO/KEGG enrichment (all ML-predicted ASD genes, prob >= 0.50)
top100_raw <- df %>% filter(Risk_Rank <= 100) %>% pull(Gene)

# --- 2. Fix outdated gene symbols (alias correction) ---
fix_symbols <- function(genes) {
  genes <- gsub("^MLL3$",    "KMT2C", genes)
  genes <- gsub("^WHSC1L1$", "NSD3",  genes)
  genes <- gsub("^CCDC64$",  "BICDL1",genes)
  genes
}
top100 <- fix_symbols(top100_raw)

cat("Top 100 genes (alias-corrected):\n"); print(top100)

# --- 3. Map to Entrez IDs (using top 100 for enrichment) ---
entrez <- bitr(top100, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat(sprintf("\nMapped %d / %d genes to Entrez IDs\n", nrow(entrez), length(top100)))
gene_ids <- entrez$ENTREZID

# --- 4. GO Enrichment (BP, MF, CC) ---
go_bp <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  readable = TRUE)
go_mf <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "MF",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  readable = TRUE)
go_cc <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "CC",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  readable = TRUE)

cat(sprintf("GO BP: %d terms | GO MF: %d terms | GO CC: %d terms\n",
            nrow(as.data.frame(go_bp)),
            nrow(as.data.frame(go_mf)),
            nrow(as.data.frame(go_cc))))

# --- 5. KEGG Enrichment ---
kegg_res <- enrichKEGG(gene = gene_ids, organism = "hsa",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
cat(sprintf("KEGG: %d pathways\n", nrow(as.data.frame(kegg_res))))

# --- 6. Save enrichment tables ---
write.csv(as.data.frame(go_bp),    "go_bp_results.csv",   row.names = FALSE)
write.csv(as.data.frame(go_mf),    "go_mf_results.csv",   row.names = FALSE)
write.csv(as.data.frame(go_cc),    "go_cc_results.csv",   row.names = FALSE)
write.csv(as.data.frame(kegg_res), "kegg_results.csv",    row.names = FALSE)

# --- 7. GO Plots ---
if (nrow(as.data.frame(go_bp)) > 0) {
  png("go_bp_dotplot.png", width = 900, height = 700)
  print(dotplot(go_bp, showCategory = 15, title = "GO Biological Process (Top 100 ASD Genes)"))
  dev.off()
  cat("go_bp_dotplot.png saved\n")

  png("go_bp_emap.png", width = 900, height = 800)
  print(emapplot(pairwise_termsim(go_bp), showCategory = 20))
  dev.off()
  cat("go_bp_emap.png saved\n")
} else { cat("No significant GO BP results\n") }

if (nrow(as.data.frame(go_mf)) > 0) {
  png("go_mf_dotplot.png", width = 900, height = 700)
  print(dotplot(go_mf, showCategory = 15, title = "GO Molecular Function (Top 100 ASD Genes)"))
  dev.off()
  cat("go_mf_dotplot.png saved\n")
} else { cat("No significant GO MF results\n") }

if (nrow(as.data.frame(go_cc)) > 0) {
  png("go_cc_dotplot.png", width = 900, height = 700)
  print(dotplot(go_cc, showCategory = 15, title = "GO Cellular Component (Top 100 ASD Genes)"))
  dev.off()
  cat("go_cc_dotplot.png saved\n")
} else { cat("No significant GO CC results\n") }

# --- 8. KEGG Plots ---
if (nrow(as.data.frame(kegg_res)) > 0) {
  png("kegg_dotplot.png", width = 900, height = 600)
  print(dotplot(kegg_res, showCategory = 15, title = "KEGG Pathways (Top 100 ASD Genes)"))
  dev.off()
  cat("kegg_dotplot.png saved\n")

  if (nrow(as.data.frame(kegg_res)) >= 2) {
    png("kegg_cnetplot.png", width = 1000, height = 900)
    print(cnetplot(kegg_res, showCategory = 10, foldChange = NULL))
    dev.off()
    cat("kegg_cnetplot.png saved\n")
  }
} else { cat("No significant KEGG results\n") }

cat("\n=== Analysis complete ===\n")
cat(sprintf("GO BP: %d | GO MF: %d | GO CC: %d | KEGG: %d significant terms\n",
            nrow(as.data.frame(go_bp)), nrow(as.data.frame(go_mf)),
            nrow(as.data.frame(go_cc)), nrow(as.data.frame(kegg_res))))

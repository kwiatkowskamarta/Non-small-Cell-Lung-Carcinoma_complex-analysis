# ==============================================================================
# PROJECT: Transcriptomic RNA-Seq Analysis (GSE298476)
# GOAL: Identify expression changes (Knockdown vs Control) and biological roles
# ==============================================================================

# --- 0. Install and Load Libraries ---
required_packages <- c("limma", "edgeR", "clusterProfiler", "org.Hs.eg.db", 
                       "enrichplot", "pheatmap", "RColorBrewer", "GEOquery")

# Install missing packages (if necessary)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)

library(limma)
library(edgeR)           # For DGEList object and TMM normalization
library(clusterProfiler) # For GO enrichment analysis
library(org.Hs.eg.db)    # Human genome database
library(enrichplot)      # Enrichment visualizations
library(pheatmap)        # Heatmaps
library(RColorBrewer)
library(GEOquery)        # To read the .soft file

# --- 1. Load Data and Metadata ---
message("--- STEP 1: Loading Data ---")

# A. Load Count Matrix (Numerical Data)
raw_data <- read.delim("data/GSE298476_raw_count_all_sample.txt", stringsAsFactors = FALSE)

# Aggregate duplicates (summing counts for the same GeneID)
counts_agg <- aggregate(raw_data[, -c(1:3)], by = list(GeneID = raw_data$Gene.ID), FUN = sum)
rownames(counts_agg) <- counts_agg$GeneID
counts_matrix <- as.matrix(counts_agg[, -1]) # Remove ID column, keep numbers only

# B. Load Metadata (CORRECTED METHOD)
# Use getGEO on local file. Returns a GSE object, not ExpressionSet.
gse_soft <- getGEO(filename = "data/GSE298476_family.soft")

# Extract list of samples (GSM objects) from GSE object
gsm_list <- GSMList(gse_soft)

# Extract titles manually to bypass the phenoData error
sample_ids <- names(gsm_list)
sample_titles <- sapply(gsm_list, function(x) Meta(x)$title)

# Create a metadata dataframe
meta_data <- data.frame(
  geo_accession = sample_ids,
  title = sample_titles,
  stringsAsFactors = FALSE
)

# Sort metadata by sample ID to ensure alignment with columns
meta_data <- meta_data[order(meta_data$geo_accession), ]

# Display for verification
print("Column names in data file:")
print(colnames(counts_matrix))
print("Sample titles extracted from SOFT file:")
print(meta_data$title)

# --- 2. Experimental Group Definition (CORRECTED) ---
message("--- STEP 2: Defining Groups based on Column Names ---")

# Sprawdzamy nazwy kolumn w macierzy zliczeń
print("Current column names:")
print(colnames(counts_matrix))

# Tworzymy grupy na podstawie nazw kolumn, a nie metadanych (aby uniknąć błędu kolejności)
# Jeśli w nazwie kolumny jest "shCtrl" -> Control, w przeciwnym razie -> Knockdown (shIGF2BP3)
group <- factor(ifelse(grepl("shCtrl", colnames(counts_matrix)), "Control", "Knockdown"))

# Upewniamy się, że Control jest poziomem bazowym (referencyjnym)
group <- relevel(group, ref = "Control") 

print("Final Defined groups:")
print(group)

# --- 3. Normalization and Filtering (edgeR + voom) ---
message("--- STEP 3: Normalization and Filtering ---")

# Create DGEList object (Standard for RNA-Seq)
dge <- DGEList(counts = counts_matrix, group = group)

# Filter low expression genes (edgeR automatic method)
keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# TMM Normalization (Trimmed Mean of M-values)
dge <- calcNormFactors(dge)

print(paste("Genes before filtering:", nrow(counts_matrix)))
print(paste("Genes after filtering:", nrow(dge)))

# --- 4. Voom Transformation and Differential Analysis ---
message("--- STEP 4: Differential Analysis (Limma-Voom) ---")

# Define design matrix
design <- model.matrix(~ group)
colnames(design) <- c("Control", "KnockdownVsControl")

# Voom transformation: converts counts to precision weights
v <- voom(dge, design, plot = FALSE)

# Fit linear model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Get results table
deg_table <- topTable(fit, coef = "KnockdownVsControl", number = Inf, sort.by = "P")

# --- 5. Gene Selection (Fallback Strategy) ---
message("--- STEP 5: Gene Selection (Fallback Strategy) ---")

# Strategy 1: Gold standard (FDR < 0.05 & |logFC| > 1)
sig_genes <- deg_table[deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) >= 1, ]
filter_method <- "FDR < 0.05 & |logFC| >= 1"

# Strategy 2: If < 10 genes, relax to raw p-value < 0.01
if (nrow(sig_genes) < 10) {
  message("Insufficient results for FDR. Relaxing criteria...")
  sig_genes <- deg_table[deg_table$P.Value < 0.01 & abs(deg_table$logFC) >= 1, ]
  filter_method <- "P-value < 0.01 & |logFC| >= 1"
}

# Strategy 3: If still < 10, take top 50 genes by P-value (Guarantees plot)
if (nrow(sig_genes) < 10) {
  message("Still few genes. Taking TOP 50 by P-value.")
  sig_genes <- head(deg_table[order(deg_table$P.Value), ], 50)
  filter_method <- "Top 50 by P-value"
}

gene_ids <- rownames(sig_genes)
print(paste("Selected", length(gene_ids), "genes using criterion:", filter_method))

# --- 6. Functional Enrichment Analysis (GO) ---
message("--- STEP 6: Enrichment Analysis (GO) ---")

# Simple ID type detection (ENSEMBL vs ENTREZ)
id_type <- ifelse(grepl("ENS", gene_ids[1]), "ENSEMBL", "ENTREZID")

if(length(gene_ids) > 0) {
  # Biological Process (BP)
  ego_bp <- enrichGO(gene = gene_ids,
                     OrgDb = org.Hs.eg.db,
                     keyType = id_type,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.5, # Loose threshold for demonstration
                     qvalueCutoff = 0.5,
                     readable = TRUE) 
  
  if(!is.null(ego_bp) && nrow(ego_bp) > 0) {
    print(dotplot(ego_bp, showCategory=10, title="GO: Biological Process"))
  } else {
    message("No significant GO terms (BP) found.")
  }
}

# --- 7. Visualization (Heatmap) ---
message("--- STEP 7: Heatmap ---")

if(length(gene_ids) > 1) {
  # Get normalized expression data for selected genes
  heatmap_data <- v$E[gene_ids, ]
  
  # Map IDs to Gene Symbols
  gene_symbols <- mapIds(org.Hs.eg.db, keys = rownames(heatmap_data), 
                         column = "SYMBOL", keytype = id_type, multiVals = "first")
  
  # Use symbols as row names
  labels <- ifelse(is.na(gene_symbols), rownames(heatmap_data), gene_symbols)
  rownames(heatmap_data) <- labels
  
  # Draw heatmap
  pheatmap(heatmap_data,
           scale = "row",
           annotation_col = data.frame(Group = group, row.names = colnames(heatmap_data)),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           main = paste("Expression Heatmap (", filter_method, ")", sep=""),
           fontsize_row = 8,
           show_colnames = TRUE)
}

message("--- FINISHED ---")
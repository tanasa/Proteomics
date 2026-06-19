# https://satijalab.org/seurat/articles/multimodal_vignette.html

# ARTICLE CITE-seq : https://www.nature.com/articles/nmeth.4380

# In this vignette, we analyze a dataset of 8,617 cord blood mononuclear cells (CBMCs): 
# transcriptomic measurements are paired with abundance estimates for 11 surface proteins, 
# whose levels are quantified with DNA-barcoded antibodies.

# https://broadinstitute.github.io/2020_scWorkshop/cite-seq.html
# the code was re-written with Code

setwd("/home/tanasa/CITEseq")
list.files()
mydir = "/home/tanasa/CITEseq"

library(Seurat)
library(Matrix)

# ============================================================
# SECTION 1: Load Data
# ============================================================

# or MOUSE_ appended to the beginning of each gene
cbmc.rna <- as.sparse(read.csv(file = "GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
    sep = ",", header = TRUE, row.names = 1))

# To make life a bit easier going forward, we're going to discard all but the top 100 most
# highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# Load in the ADT UMI matrix
cbmc.adt <- as.sparse(read.csv(file = "GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz",
    sep = ",", header = TRUE, row.names = 1))

# Note that since measurements were made in the same cells, the two matrices have identical
# column names
all.equal(colnames(cbmc.rna), colnames(cbmc.adt))

# Remove poorly enriched antibodies (optional)
cbmc.adt <- cbmc.adt[setdiff(rownames(cbmc.adt), c("CCR5", "CCR7", "CD10")), ]
 
# Check overlap between RNA and ADT cells
length(intersect(colnames(cbmc.rna), colnames(cbmc.adt))) /
  length(union(colnames(cbmc.rna), colnames(cbmc.adt)))


head(cbmc.rna)

head(cbmc.adt)

dim(cbmc.rna)

dim(cbmc.adt)

# Set up a Seurat object with the RNA and ADT data


# ============================================================
# SECTION 2: RNA Clustering
# ============================================================

cbmc <- CreateSeuratObject(counts = cbmc.rna)

# Standard RNA preprocessing
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)

ElbowPlot(cbmc, ndims = 50)

# Cluster on first 25 PCs
cbmc <- FindNeighbors(cbmc, dims = 1:25)
cbmc <- FindClusters(cbmc, resolution = 0.8)
cbmc <- RunTSNE(cbmc, dims = 1:25, method = "FIt-SNE")

# Find RNA markers per cluster
cbmc.rna.markers <- FindAllMarkers(cbmc,
                                   max.cells.per.ident = 100,
                                   logfc.threshold = log(2),
                                   only.pos = TRUE,
                                   min.diff.pct = 0.3)

# Annotate clusters manually based on marker genes
new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", "CD14+ Mono",
                     "Mouse", "B", "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK",
                     "CD34+", "Multiplets", "Mouse", "Eryth", "Mk", "Mouse", "DC", "pDCs")

names(new.cluster.ids) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, new.cluster.ids)

# Visualize RNA clusters
DimPlot(cbmc, label = TRUE, reduction = "tsne") + NoLegend()

# save(cbmc, cbmc.rna.markers, cbmc.adt, file = Rda.RNA.path)

# ============================================================
# SECTION 3: Add ADT (Protein) Assay
# ============================================================

# Add ADT as a second assay
cbmc[["ADT"]] <- CreateAssayObject(counts = cbmc.adt)

# Check raw ADT counts
GetAssayData(cbmc, layer = "counts", assay = "ADT")[1:3, 1:3]  # v5: layer= not slot=

# Normalize ADT with CLR (recommended for protein data)
cbmc <- NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR")
cbmc <- ScaleData(cbmc, assay = "ADT")


options(repr.plot.width = 16,    # width in inches
        repr.plot.height = 9,    # height in inches
        repr.plot.res = 150)     # resolution


# ============================================================
# SECTION 4: Visualize Protein Levels on RNA Clusters
# ============================================================

DefaultAssay(cbmc) <- "RNA"

# Protein (top row) vs RNA (bottom row) side by side
FeaturePlot(cbmc,
            features = c("adt_CD3", "adt_CD11c", "adt_CD8", "adt_CD16",
                         "CD3E", "ITGAX", "CD8A", "FCGR3A"),
            min.cutoff = "q05", max.cutoff = "q95", ncol = 4)

options(repr.plot.width = 20,    # width in inches
        repr.plot.height = 9,    # height in inches
        repr.plot.res = 150)     # resolution


FeaturePlot(cbmc,
            features = c("adt_CD4", "adt_CD45RA", "adt_CD56", "adt_CD14",
                         "adt_CD19", "adt_CD34", "CD4", "PTPRC",
                         "NCAM1", "CD14", "CD19", "CD34"),
            min.cutoff = "q05", max.cutoff = "q95", ncol = 6)



# Ridge plots
RidgePlot(cbmc, features = c("adt_CD3", "adt_CD11c", "adt_CD8", "adt_CD16"), ncol = 2)

options(repr.plot.width = 8,    # width in inches
        repr.plot.height = 4,   # height in inches
        repr.plot.res = 150)     # resolution

# Biaxial scatter plots (like FACS)
FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")
FeatureScatter(cbmc, feature1 = "adt_CD3",  feature2 = "CD3E")  # protein vs RNA

# CD4 vs CD8 in T cells only
tcells <- subset(cbmc, idents = c("Naive CD4 T", "Memory CD4 T", "CD8 T"))
FeatureScatter(tcells, feature1 = "adt_CD4", feature2 = "adt_CD8")

options(repr.plot.width = 4,     # width in inches
        repr.plot.height = 4,    # height in inches
        repr.plot.res = 150)     # resolution

#  What fraction of T cells are double negative in gene expression? (CD4- and CD8-

ncol(subset(tcells, subset = CD4 == 0 & CD8A == 0)) / ncol(tcells)

# RNA double negatives (CD4- CD8-): high due to dropout

DefaultAssay(tcells) <- "RNA"  # work with ADT count matrix

FeatureScatter(tcells, feature1 = "CD4", feature2 = "CD8A")


options(repr.plot.width = 8,    # width in inches
        repr.plot.height = 4,   # height in inches
        repr.plot.res = 150)     # resolution

# What fraction of T cells are double negative in protein expression? (CD4- and CD8-)
# length(cells) / length(tcells@cell.names)

DefaultAssay(tcells) <- "ADT"  # work with ADT count matrix

FeatureScatter(tcells, feature1 = "adt_CD4", feature2 = "adt_CD8")

ncol(subset(tcells, subset = adt_CD4 < 1 & adt_CD8 < 1)) / ncol(tcells)


# ============================================================
# SECTION 5: Differential Protein Expression Between Clusters
# ============================================================

# Downsample to 300 cells per cluster for speed
cbmc.small <- subset(cbmc, downsample = 300)

# Find ADT markers per cluster
adt.markers <- FindAllMarkers(cbmc.small, assay = "ADT", only.pos = TRUE)

# Heatmap of protein markers (requires ScaleData first in Seurat v5)
cbmc.small <- ScaleData(cbmc.small, assay = "ADT")

DoHeatmap(cbmc.small,
          features = unique(adt.markers$gene),
          assay = "ADT",
          angle = 90) + NoLegend()


# You can see that our unknown cells co-express both myeloid and lymphoid markers (true at the
# RNA level as well). They are likely cell clumps (multiplets) that should be discarded. We'll
# remove the mouse cells now as well

# Remove multiplets and mouse cells
# cbmc <- subset(cbmc, idents = c("Multiplets", "Mouse"), invert = TRUE)

# ============================================================
# SECTION 6: Cluster Directly on Protein (ADT) Levels
# ============================================================

DefaultAssay(cbmc) <- "ADT"

# PCA on ADT (optional visualization — only 10 proteins so limited value)
cbmc <- ScaleData(cbmc)  # required in Seurat v5 before RunPCA
cbmc <- RunPCA(cbmc,
               features = rownames(cbmc),
               reduction.name = "pca_adt",
               reduction.key  = "pca_adt_",
               verbose = FALSE)

DimPlot(cbmc, reduction = "pca_adt")
ElbowPlot(cbmc)

# Use Euclidean distance matrix instead of PCA (better for small number of proteins)

adt.data <- GetAssayData(cbmc, assay = "ADT", layer = "data")  # v5: layer= not slot=
adt.dist  <- dist(t(adt.data))

# Stash RNA cluster IDs before reclustering
cbmc[["rnaClusterID"]] <- Idents(cbmc)

# tSNE and SNN graph from ADT distance matrix

cbmc[["tsne_adt"]] <- RunTSNE(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
cbmc[["adt_snn"]]  <- FindNeighbors(adt.dist)$snn

cbmc <- FindClusters(cbmc, resolution = 0.2, graph.name = "adt_snn")

# Compare RNA vs ADT clustering
clustering.table <- table(Idents(cbmc), cbmc$rnaClusterID)
clustering.table

# Annotate ADT-based clusters
new.cluster.ids <- c("CD4 T", "CD14+ Mono", "NK", "B", "CD8 T", "NK",
                     "CD34+", "T/Mono doublets", "CD16+ Mono", "pDCs", "B")

names(new.cluster.ids) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, new.cluster.ids)

library("ggplot2")
library(patchwork)

options(repr.plot.width = 12,    # width in inches
        repr.plot.height = 6,    # height in inches
        repr.plot.res = 150)     # resolution

# Side-by-side comparison: RNA clusters vs ADT clusters on ADT tSNE
tsne_rnaClusters <- DimPlot(cbmc, reduction = "tsne_adt",
                             group.by = "rnaClusterID", pt.size = 0.5) + NoLegend() +
  ggtitle("Clustering based on scRNA-seq") +
  theme(plot.title = element_text(hjust = 0.5))

tsne_rnaClusters <- LabelClusters(tsne_rnaClusters, id = "rnaClusterID", size = 4)


tsne_adtClusters <- DimPlot(cbmc, reduction = "tsne_adt", pt.size = 0.5) + NoLegend() +
  ggtitle("Clustering based on ADT signal") +
  theme(plot.title = element_text(hjust = 0.5))

tsne_adtClusters <- LabelClusters(tsne_adtClusters, id = "ident", size = 4)

wrap_plots(list(tsne_rnaClusters, tsne_adtClusters), ncol = 2)



cat("

The Problem with Single-Modality Clustering

RNA alone has limitations:

High dropout (many genes read as zero by chance)
CD4 and CD8 T cells look nearly identical transcriptomically
~83% of T cells appear <double negative> in RNA vs only ~1% in protein

Protein (ADT) alone has limitations:

Only 10-13 proteins measured — can't distinguish all cell types
Good for major lineages (T, B, NK, Mono) but loses resolution for rare populations 
like DCs, Mk, Eryth that lack distinguishing antibodies

What WNN Does : 

Weighted Nearest Neighbor combines both modalities, but learns how much to trust each one per cell:

Final cell similarity = w_RNA × (RNA distance) + w_ADT × (ADT distance)

The weights are per cell, not global — so:

T cells    → high ADT weight  (CD4/CD8 protein clearly separates them)
DCs        → high RNA weight  (no good DC antibody in the panel)
Monocytes  → balanced         (both modalities informative)

")



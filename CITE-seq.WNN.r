# https://hoohm.github.io/CITE-seq-Count/

setwd("/home/tanasa/CITEseq")

# We use the CITE-seq dataset from (Stuart*, Butler* et al, Cell 2019), 
# which consists of 30,672 scRNA-seq profiles measured alongside 
# a panel of 25 antibodies from bone marrow. 
# The object contains two assays, RNA and antibody-derived tags (ADT).

# In (Hao*, Hao* et al, Cell 2021 Login to Jenni), we introduce ‘weighted-nearest neighbor’ (WNN) analysis, 
# an unsupervised framework to learn the relative utility of each data type in each cell, 
# enabling an integrative analysis of multiple modalities.

# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz

library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)

# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis

cat("WNN analysis of CITE-seq, RNA + ADT")

InstallData("bmcite")
bm <- LoadData(ds = "bmcite")

DefaultAssay(bm) <- 'RNA'

bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(bm) <- 'ADT'

# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the default PCA

VariableFeatures(bm) <- rownames(bm[["ADT"]])

bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>%   ScaleData() %>% RunPCA(reduction.name = 'apca')

# str(bm)

bm[["RNA"]]

str(bm[["RNA"]])

# bm[["RNA"]]@cells
# bm[["RNA"]]@features

# bm[["ADT"]]@cells
bm[["ADT"]]@features

# For each cell, we calculate its closest neighbors in the dataset based on a weighted combination of RNA and protein similarities. 
# The cell-specific modality weights and multimodal neighbors are calculated in a single function, 
# which takes ~2 minutes to run on this datas

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]], 
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight

bm <- FindMultiModalNeighbors(
  bm, 
  reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), 
  modality.weight.name = "RNA.weight"
)

bm[['weighted.nn']]

bm[["wknn"]]

# a UMAP visualization of the data based on a weighted combination of RNA and protein data. 

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
bm <- FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

options(repr.plot.width = 10,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution

p1 <- DimPlot(bm, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2

# We can also compute UMAP visualization based on only the RNA and protein data and compare. 
# We find that the RNA analysis is more informative than the ADT analysis in identifying progenitor states 
# (the ADT panel contains markers for differentiated cells), 
# while the converse is true of T cell states (where the ADT analysis outperforms RNA).

bm <- RunUMAP(bm, reduction = 'pca', dims = 1:30, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

bm <- RunUMAP(bm, reduction = 'apca', dims = 1:18, assay = 'ADT', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

p3 <- DimPlot(bm, reduction = 'rna.umap', group.by = 'celltype.l2', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()

p4 <- DimPlot(bm, reduction = 'adt.umap', group.by = 'celltype.l2', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()

p3 + p4

p5 <- FeaturePlot(bm, features = c("adt_CD45RA","adt_CD16","adt_CD161"),
                  reduction = 'wnn.umap', max.cutoff = 2, 
                  cols = c("lightgrey","darkgreen"), ncol = 3)

p6 <- FeaturePlot(bm, features = c("rna_TRDC","rna_MPO","rna_AVP"), 
                  reduction = 'wnn.umap', max.cutoff = 3, ncol = 3)

p5 / p6

# Each of the populations with the highest RNA weights represent progenitor cells, 
# while the populations with the highest protein weights represent T cells.

str(bm)

bm@meta.data

colnames(bm@meta.data)

@ graphs      :List of 2
  .. ..$ wknn:Formal class 'Graph' [package "SeuratObject"] with 7 slots
  .. .. .. ..@ assay.used: Named chr "RNA"
  .. .. .. .. ..- attr(*, "names")= chr "pca"
  .. .. .. ..@ i         : int [1:1014222] 0 2247 3305 3637 7567 8633 12399 12915 15010 15447 ...
  .. .. .. ..@ p         : int [1:30673] 0 24 45 81 113 143 180 204 228 275 ...
  .. .. .. ..@ Dim       : int [1:2] 30672 30672
  .. .. .. ..@ Dimnames  :List of 2
  .. .. .. .. ..$ : chr [1:30672] "a_AAACCTGAGCTTATCG-1" "a_AAACCTGAGGTGGGTT-1" "a_AAACCTGAGTACATGA-1" "a_AAACCTGCAAACCTAC-1" ...
  .. .. .. .. ..$ : chr [1:30672] "a_AAACCTGAGCTTATCG-1" "a_AAACCTGAGGTGGGTT-1" "a_AAACCTGAGTACATGA-1" "a_AAACCTGCAAACCTAC-1" ...
  .. .. .. ..@ x         : num [1:1014222] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..@ factors   : list()
  .. ..$ wsnn:Formal class 'Graph' [package "SeuratObject"] with 7 slots
  .. .. .. ..@ assay.used: Named chr "RNA"
  .. .. .. .. ..- attr(*, "names")= chr "pca"
  .. .. .. ..@ i         : int [1:2195882] 0 161 340 1070 1185 2109 2247 2416 2540 2748 ...
  .. .. .. ..@ p         : int [1:30673] 0 109 129 177 247 331 438 476 529 620 ...
  .. .. .. ..@ Dim       : int [1:2] 30672 30672
  .. .. .. ..@ Dimnames  :List of 2
  .. .. .. .. ..$ : chr [1:30672] "a_AAACCTGAGCTTATCG-1" "a_AAACCTGAGGTGGGTT-1" "a_AAACCTGAGTACATGA-1" "a_AAACCTGCAAACCTAC-1" ...
  .. .. .. .. ..$ : chr [1:30672] "a_AAACCTGAGCTTATCG-1" "a_AAACCTGAGGTGGGTT-1" "a_AAACCTGAGTACATGA-1" "a_AAACCTGCAAACCTAC-1" ...
  .. .. .. ..@ x         : num [1:2195882] 1 0.2903 0.0811 0.25 0.0811 ...
  .. .. .. ..@ factors   : list()
  ..@ neighbors   :List of 1
  .. ..$ weighted.nn:Formal class 'Neighbor' [package "SeuratObject"] with 5 slots
  .. .. .. ..@ nn.idx    : num [1:30672, 1:20] 21408 26747 2354 14276 2344 ...
  .. .. .. ..@ nn.dist   : num [1:30672, 1:20] 0.146 0.204 0.295 0.336 0.296 ...
  .. .. .. ..@ alg.idx   : NULL
  .. .. .. ..@ alg.info  : list()
  .. .. .. ..@ cell.names: chr [1:30672] "a_AAACCTGAGCTTATCG-1" "a_AAACCTGAGGTGGGTT-1" "a_AAACCTGAGTACATGA-1" "a_AAACCTGCAAACCTAC-1" ...
  ..@ reductions  :List of 6

bm@reductions

bm@neighbors

bm@neighbors$weighted.nn

bm@graphs
str(bm@graphs)

bm@graphs$wknn

bm@graphs$wsnn

# Scatter plot of RNA vs ADT
FeatureScatter(bm, feature1 = "rna_CD3E", feature2 = "adt_CD3")

options(repr.plot.width = 6,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution

# RNA gene expression
FeaturePlot(bm, features = "CD3E", reduction = "wnn.umap")

# ADT protein levels (note the "adt_" prefix)
FeaturePlot(bm, features = "adt_CD3", reduction = "wnn.umap")

options(repr.plot.width = 12,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution

# Side by side comparison
FeaturePlot(bm, features = c("CD3E", "adt_CD3"), reduction = "wnn.umap")

options(repr.plot.width = 6,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution

# Visualize the weighted graph
DimPlot(bm, reduction = "wnn.umap", group.by = "seurat_clusters")

options(repr.plot.width = 12,   # width in inches
        repr.plot.height = 4,   # height in inches
        repr.plot.res = 150)    # 

# Check modality weights (RNA vs protein contribution per cell)
VlnPlot(bm, features = "RNA.weight", group.by = "seurat_clusters", pt.size = 0)
VlnPlot(bm, features = "ADT.weight", group.by = "seurat_clusters", pt.size = 0)

options(repr.plot.width = 26,   # width in inches
        repr.plot.height = 12,   # height in inches
        repr.plot.res = 150)    # resolution

# All ADT markers across clusters
DefaultAssay(bm) <- "ADT"
DotPlot(bm, features = rownames(bm[["ADT"]]), group.by = "seurat_clusters") +
  RotatedAxis()

# Scatter plot of RNA vs ADT
# FeatureScatter(bm, feature1 = "CD3E", feature2 = "adt_CD3")

cat("The differences in the protein expression between the clusters")

DefaultAssay(bm) <- "ADT"

adt_markers <- FindAllMarkers(
  bm,
  assay = "ADT",
  only.pos = FALSE,    # include both up and down
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# View top markers per cluster
library(dplyr)

adt_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Compare two specific clusters

# Cluster 1 vs Cluster 3

FindMarkers(
  bm,
  assay = "ADT",
  ident.1 = 1,
  ident.2 = 3,
  min.pct = 0.1
)

options(repr.plot.width = 18,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution

VlnPlot(bm, assay = "ADT",
        features = c("CD3", "CD19", "CD14"),
        group.by = "seurat_clusters",
        pt.size = 0)

options(repr.plot.width = 12,   # width in inches
        repr.plot.height = 12,   # height in inches
        repr.plot.res = 150)    # resolution

FeaturePlot(bm,
            features = c("CD3", "CD19", "CD14", "CD56"),
            reduction = "wnn.umap",
            ncol = 2)

cat('

Important Note on ADT Normalization : 

# CLR normalization is recommended for ADT before DE
bm <- NormalizeData(bm, assay = "ADT",
                    normalization.method = "CLR",
                    margin = 2)  # margin=2 normalizes across cells

')



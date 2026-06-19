# https://satijalab.org/seurat/articles/multimodal_vignette.html

# ARTICLE CITE-seq : https://www.nature.com/articles/nmeth.4380

# In this vignette, we analyze a dataset of 8,617 cord blood mononuclear cells (CBMCs): 
# transcriptomic measurements are paired with abundance estimates for 11 surface proteins, 
# whose levels are quantified with DNA-barcoded antibodies.

setwd("/home/tanasa/CITEseq")
list.files()

library(Seurat)
library(Matrix)

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

head(cbmc.rna)

head(cbmc.adt)

dim(cbmc.rna)

dim(cbmc.adt)

# Set up a Seurat object with the RNA and ADT data

# creates a Seurat object based on the scRNA-seq data
cbmc <- CreateSeuratObject(counts = cbmc.rna)

# We can see that by default, the cbmc object contains an assay storing RNA measurement
Assays(cbmc)

# Create a new assay to store ADT information
adt_assay <- CreateAssay5Object(counts = cbmc.adt)

# Add this assay to the previously created Seurat object
cbmc[["ADT"]] <- adt_assay

# Validate that the object now contains multiple assays
Assays(cbmc)

rownames(cbmc[["ADT"]])

cbmc[["RNA"]] 

cbmc[["ADT"]] 

# Cluster cells on the basis of their scRNA-seq profiles

# Note that all operations below are performed on the RNA assay Set and verify that the
# default assay is RNA

DefaultAssay(cbmc) <- "RNA"
DefaultAssay(cbmc)

# Perform visualization and clustering steps
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)
cbmc <- FindNeighbors(cbmc, dims = 1:30)
cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:30)

DimPlot(cbmc, label = TRUE)

# Normalize ADT data

library("ggplot2")

DefaultAssay(cbmc) <- "ADT"
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2)
DefaultAssay(cbmc) <- "RNA"

# Note that the following command is an alternative but returns the same result
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2, assay = "ADT")

# Now, we will visualize CD19 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(cbmc) <- "ADT"
p1 <- FeaturePlot(cbmc, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")

DefaultAssay(cbmc) <- "RNA"
p2 <- FeaturePlot(cbmc, "CD19") + ggtitle("CD19 RNA")


options(repr.plot.width = 12,    # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)     # resolution


# Place plots side-by-side
p1 | p2

# Alternately, we can use specific assay keys to specify a specific modality Identify the key
# for the RNA and protein assays
Key(cbmc[["RNA"]])

Key(cbmc[["ADT"]])

# Now, we can include the key in the feature name, which overrides the default assay
p1 <- FeaturePlot(cbmc, "adt_CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
p2 <- FeaturePlot(cbmc, "rna_CD19") + ggtitle("CD19 RNA")
p1 | p2

cat("Identify cell surface markers for scRNA-seq clusters")

# We can leverage our paired CITE-seq measurements to help annotate clusters derived from scRNA-seq, and to identify both protein and RNA markers.

# As we know that CD19 is a B cell marker, we can identify cluster 6 as expressing CD19 on the
# surface
VlnPlot(cbmc, "adt_CD19")

# We can also identify alternative protein and RNA markers for this cluster through
# differential expression

adt_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "ADT")
rna_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "RNA")


adt_markers

head(adt_markers)

rna_markers

head(rna_markers)

cat("Additional visualizations of multimodal data")

# Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if
# desired by using HoverLocator and FeatureLocator
FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")

# view relationship between protein and RNA
FeatureScatter(cbmc, feature1 = "adt_CD3", feature2 = "rna_CD3E")

FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8")

# Let's look at the raw (non-normalized) ADT counts. You can see the values are quite high,
# particularly in comparison to RNA values This is due to the significantly higher protein
# copy number in cells, which significantly reduces 'drop-out' in ADT data
FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8", slot = "counts")



cat("Identify differentially expressed proteins between clusters")

# https://broadinstitute.github.io/2020_scWorkshop/cite-seq.html

# Downsample the clusters to a maximum of 300 cells each (makes the heatmap easier to see for
# small clusters)
cbmc.small <- subset(cbmc, downsample = 300)

# Find protein markers for all clusters, and draw a heatmap
adt.markers <- FindAllMarkers(cbmc.small, assay = "ADT", only.pos = TRUE)

# Step 1: Set ADT as default assay
DefaultAssay(cbmc.small) <- "ADT"

# Step 2: Scale the ADT data first
cbmc.small <- ScaleData(cbmc.small, assay = "ADT")

# Step 3: Now run DoHeatmap
DoHeatmap(cbmc.small, 
          features = unique(adt.markers$gene), 
          assay = "ADT", 
          angle = 90) + NoLegend()





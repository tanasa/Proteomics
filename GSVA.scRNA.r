# https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html

# The GSVA package provides the implementation of four single-sample gene set enrichment methods, 
# concretely zscore, plage, ssGSEA and its own called GSVA

# We import the PBMC data using the TENxPBMCData package :

library(SingleCellExperiment)
library(TENxPBMCData)

sce <- TENxPBMCData(dataset="pbmc4k")
sce

head(rownames(sce))
tail(rownames(sce))

colData(sce)

rowData(sce)

assayNames(sce)

dim(assay(sce, "counts"))
assay(sce, "counts")

cat("Add log-normalized counts")

library(scuttle)
sce <- logNormCounts(sce)

assayNames(sce)

cat("Annotate cell types using GSVA")

library(GSEABase)
library(GSVA)

fname <- file.path(system.file("extdata", package="GSVAdata"),
                   "pbmc_cell_type_gene_set_signatures.gmt.gz")
gsets <- readGMT(fname)
gsets

# str(gsets)

names(gsets)

geneIds(gsets)

# Note that while gene identifers in the sce object correspond to Ensembl stable identifiers (ENSG...), 
# the gene identifiers in the gene sets are HGNC gene symbol

# str(gsets)

gsvaAnnotation(sce) <- ENSEMBLIdentifier("org.Hs.eg.db")
gsvaAnnotation(sce)

cat("Build parameter object")
# By default, the expression values in the logocounts assay will be selected for downstream analysis.

gsvapar <- gsvaParam(sce, gsets)
gsvapar

cat("Calculate GSVA scores")

gsvaranks <- gsvaRanks(gsvapar)
gsvaranks

#  we calculate the GSVA scores using the output of gsvaRanks() as input to the function gsvaScores().

es <- gsvaScores(gsvaranks)
es

class(es)
dim(es)

# gene sets x cells
# 20 4340

es[1:10, 1:10]

head(rownames(es))
head(colnames(es))

es["B_CELLS_MEMORY", ]

score.mat <- assay(es)

dim(score.mat)

score.mat[1:10,1:10]

head(rownames(score.mat))
head(colnames(score.mat))

colData(es)

colnames(colData(es))

dim(score.mat)

nrow(colData(es))

cat("Using GSVA scores to assign cell types")
# https://bioconductor.org/books/3.16/OSCA.basic/clustering.html

# we use GSVA scores to build a nearest-neighbor graph of the cells using the function buildSNNGraph() from the scran package

library(bluster)

g <- makeSNNGraph(t(assay(es)), k=20)

library(igraph)

colLabels(es) <- factor(cluster_walktrap(g)$membership)
table(colLabels(es))

library(RColorBrewer)

res <- prcomp(assay(es))
varexp <- res$sdev^2 / sum(res$sdev^2)
nclusters <- nlevels(colLabels(es))
hmcol <- colorRampPalette(brewer.pal(nclusters, "Set1"))(nclusters)
par(mar=c(4, 5, 1, 1))
plot(res$rotation[, 1], res$rotation[, 2], col=hmcol[colLabels(es)], pch=19,
     xlab=sprintf("PCA 1 (%.0f%%)", varexp[1]*100),
     ylab=sprintf("PCA 2 (%.0f%%)", varexp[2]*100),
     las=1, cex.axis=1.2, cex.lab=1.5)
legend("topright", gsub("_", " ", levels(colLabels(es))), fill=hmcol, inset=0.01)

# if we want to better understand why a specific cell type is annotated to a given cell, we can use the gsvaEnrichment() function, 
# which will show a GSEA enrichment plot. 

# firsteosinophilcell <- which(colLabels(es) == "EOSINOPHILS")[1]
# par(mar=c(4, 5, 1, 1))

# gsvaEnrichment(gsvaranks, column=firsteosinophilcell, geneSet="EOSINOPHILS",
#               cex.axis=1.2, cex.lab=1.5, plot="ggplot")






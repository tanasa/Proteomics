# https://www.bioconductor.org/packages/release/bioc/vignettes/scRepertoire/inst/doc/vignette.html

# https://www.borch.dev/uploads/screpertoire/
# https://zenodo.org/records/18187313

cat("scRepertoire")

library("scRepertoire")

setwd("/home/tanasa/scRepertoire/data")

# contig_list = load("contig_list.rda")
# obj = load("scRep_example.rda")

# https://pubmed.ncbi.nlm.nih.gov/33622974/

# data("contig_list")
# data("scRep_example")

str(contig_list)

str(contig_list, max.level=0)

sapply(contig_list, class)

length(contig_list)

str(contig_list, max.level=1, list.len=3)

head(contig_list[[1]])
tail(contig_list[[1]])

head(contig_list[[2]])
tail(contig_list[[2]])

head(contig_list[[8]])
tail(contig_list[[8]])

cat("Combining Contigs into Clones")

# The combineTCR() function processes a list of TCR sequencing results, consolidating them to the level of individual cell barcodes. 
# It handles potential issues with repeated barcodes by adding prefixes from samples and ID parameters. 
# The output includes combined reads into clone calls by nucleotide sequence (CTnt), amino acid sequence (CTaa), 
# VDJC gene sequence (CTgene), or a combination of nucleotide and gene sequence (CTstrict).

# Key Parameter(s) for combineTCR()

#    input.data: A list of filtered contig annotations (e.g., filtered_contig_annotations.csv from 10x Cell Ranger) 
#                or outputs from loadContigs().
#    samples: Labels for your samples (recommended).
#    ID: Additional sample labels (optional).
#    removeNA: If TRUE, removes any cell barcode with an NA value in at least one chain (default is FALSE).
#    removeMulti: If TRUE, removes any cell barcode with more than two immune receptor chains (default is FALSE).
#    filterMulti: If TRUE, isolates the top two expressed chains in cell barcodes with multiple chains (default is FALSE).
#    filterNonproductive: If TRUE, removes non-productive chains if the variable exists in the contig data (default is TRUE).

cat("combined.TCR")

combined.TCR <- combineTCR(contig_list, 
                           samples = c("P17B", "P17L", "P18B", "P18L", 
                                       "P19B","P19L", "P20B", "P20L"),
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

head(combined.TCR[[1]])

head(combined.TCR[[2]])

head(combined.TCR[[8]])

cat("combineBCR")

# The combineBCR() function is the primary tool for processing raw B cell receptor contig data into a format ready for analysis. 
# It is analogous to combineTCR() but includes specialized logic for handling the complexities of BCRs, 
# such as somatic hypermutation. 
# The function consolidates contigs into a single data frame per sample, with each row representing a unique cell.

# threshold: The similarity threshold for clonalCluster() (default: 0.85).
# dist_type: The metric for calculating difference: “levenshtein” (default), “hamming”, 
# “nw” (Needleman-Wunsch), or “sw” (Smith-Waterman).
# dist_mat: The substitution matrix used if dist_type is “nw” or “sw” (e.g., “BLOSUM62”, “PAM30”)
# normalize: How to handle the threshold: “none”, “length”, or “maxlen”.
# chain: The chain to use for clustering when call.related.clones = TRUE (default: both).
# sequence: The sequence type (nt or aa) for clustering (default: nt).
# use.V, use.J: If TRUE, sequences must share the same V/J gene to be clustered (default: TRUE)

BCR.contigs <- read.csv("https://www.borch.dev/uploads/contigs/b_contigs.csv")

dim(BCR.contigs)
BCR.contigs 

head(BCR.contigs)

# Combine using the default similarity clustering
combined.BCR.clustered <- combineBCR(BCR.contigs, 
                                     samples = "Patient1", 
                                     threshold = 0.85)


# The CTstrict column contains cluster IDs (e.g., "cluster.1")

head(combined.BCR.clustered[[1]][, c("barcode", "CTstrict", "IGH", "cdr3_aa1")])

# Advanced: Grouping by Alignment

# For more biological accuracy—specifically when analyzing amino acid sequences—you can use alignment metrics. 
# This allows the clustering to penalize conservative amino acid changes less than radical changes.

# Here we use Needleman-Wunsch alignment with the BLOSUM80 substitution matrix :

combined.BCR.aligned <- combineBCR(BCR.contigs, 
                                   samples = "Patient1",
                                   sequence = "aa",        
                                   dist_type = "nw",      
                                   dist_mat = "BLOSUM80",  
                                   threshold = 0.85)

head(combined.BCR.aligned[[1]][, c("barcode", "CTstrict", "IGH", "cdr3_aa1")])

# Filtering and Cleaning Data

cleaned.BCR <- combineBCR(BCR.contigs,
                          samples = "Patient1",
                          filterNonproductive = TRUE,
                          filterMulti = TRUE)

head(cleaned.BCR[[1]])

# combineBCR() is designed for processing B cell repertoire data, going beyond simple contig aggregation to incorporate advanced clustering 
# based on CDR3 sequence similarity. 

# this enables the identification of clonally related B cells, crucial for studying B cell development, affinity maturation, 
# and humoral immune responses

# addVariable: Adding Variables for Plotting

combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "Type", 
                            variables = rep(c("B", "L"), 4))

head(combined.TCR[[1]])

# subsetClones: Filter Out Clonal Information

subset1 <- subsetClones(combined.TCR, 
                        name = "sample", 
                        variables = c("P18L", "P18B"))

head(subset1[[1]][,1:4])


subset2 <- combined.TCR[c(3,4)]

head(subset2[[1]][,1:4])

# exportClones: Save Clonal Data

# exportClones(combined, 
#             write.file = TRUE,
#             dir = "~/"
#             file.name = "clones.csv")

immunarch <- exportClones(combined.TCR, 
                          format = "immunarch", 
                          write.file = FALSE)
head(immunarch[[1]][[1]])

packageVersion("scRepertoire")

# clonalBin: Bin Clones by Frequency or Proportion : deprecated function in 2.6.2 ??

# The clonalBin() function adds a clonal grouping variable (cloneSize) to the output of combineTCR(), combineBCR(), or combineExpression(). 
# This function calculates the clonal frequency and proportion, 
# then bins clones into categories based on customizable thresholds. 

# suppressMessages(library(scRepertoire))
# combined.TCR <- clonalBin(combined.TCR, clone.call = "strict")

grep("clonal", ls("package:scRepertoire"), value = TRUE)



cat("VISUALIZATIONS")

# cloneCall: How to call clones.

#    gene - use the VDJC genes comprising the TCR/Ig
#    nt - use the nucleotide sequence of the CDR3 region
#    aa - use the amino acid sequence of the CDR3 region
#    strict - use the VDJC genes comprising the TCR + the nucleotide sequence of the CDR3 region. 
# This is the proper definition of clonotype. For combineBCR() strict refers to the edit distance clusters + Vgene of the Ig.

cat("clonalQuant: Quantifying Unique Clones")

clonalQuant(combined.TCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = TRUE)  

clonalQuant(combined.TCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = FALSE)

cat("clonalAbundance: Distribution of Clones by Size")

clonalAbundance(combined.TCR, 
                cloneCall = "gene", 
                scale = FALSE)

clonalAbundance(combined.TCR, 
                cloneCall = "gene", 
                scale = TRUE)

cat("clonalLength: Distribution of Sequence Lengths")

clonalLength(combined.TCR, 
             cloneCall="aa", 
             chain = "both") 

clonalLength(combined.TCR, 
             cloneCall="aa", 
             chain = "TRA", 
             scale = TRUE) 

cat("clonalCompare: Clonal Dynamics Between Categorical Variables")

clonalCompare(combined.TCR, 
                  top.clones = 10, 
                  samples = c("P17B", "P17L"), 
                  cloneCall="aa", 
                  graph = "alluvial")

clonalCompare(combined.TCR, 
              top.clones = 10,
              highlight.clones = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF", "NA_CASSVRRERANTGELFF"),
              relabel.clones = TRUE,
              samples = c("P17B", "P17L"), 
              cloneCall="aa", 
              graph = "alluvial")

clonalCompare(combined.TCR, 
              clones = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF", "NA_CASSVRRERANTGELFF"),
              relabel.clones = TRUE,
              samples = c("P17B", "P17L"), 
              cloneCall="aa", 
              graph = "alluvial")

cat("clonalScatter: Scatterplot of Two Variables")

clonalScatter(combined.TCR, 
              cloneCall ="gene", 
              x.axis = "P18B", 
              y.axis = "P18L",
              dot.size = "total",
              graph = "proportion")

cat("Visualizing Clonal Dynamics")

cat("clonalHomeostasis: Examining Clonal Space")

# Set BEFORE plotting — controls the output cell size
options(repr.plot.width = 14,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution (dpi)

library(cowplot)

p1 <- clonalHomeostasis(combined.TCR, cloneCall = "gene")
p2 <- clonalHomeostasis(combined.TCR, 
                        cloneCall = "gene",
                        cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, 
                                      Large = 0.3, Hyperexpanded = 1))
plot_grid(p1, p2, 
          ncol = 2, 
          labels = c("Default bins", "Custom bins"),
          label_size = 10)

# p1 + p2

# clonalHomeostasis() provides an assessment of how different “sizes” of clones 
# (based on their proportional abundance) contribute to the overall repertoire.

cat("clonalProportion: Examining Space Occupied by Ranks of Clones")

# Like clonal space homeostasis, clonalProportion() also categorizes clones into separate bins. 
# The key difference is that instead of looking at the relative proportion of the clone to the total, 
# clonalProportion() ranks the clones by their total count or frequency of occurrence and then places them into predefined bins.

options(repr.plot.width = 14,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution 

# p3 = clonalProportion(combined.TCR, 
#                 cloneCall = "gene") 

# p4  = clonalProportion(combined.TCR, 
#                 cloneCall = "nt",
#                 clonalSplit = c(1, 5, 10, 100, 1000, 10000))

# p3 + p4

cat("Summarizing Repertoires")

options(repr.plot.width = 8,   # width in inches
        repr.plot.height = 12,   # height in inches
        repr.plot.res = 150)    # resolution 

vizGenes(combined.TCR,
         x.axis = "TRBV",
         y.axis = NULL, # No specific y-axis variable, will group all samples
         plot = "barplot",
         summary.fun = "proportion") 

# Peripheral Blood Samples
vizGenes(combined.TCR[c("P17B", "P18B", "P19B", "P20B")],
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap",
         summary.fun = "percent") # Display percentages

cat("The intensity reflects the percentage of each V-J pairing.")

# Lung Samples
vizGenes(combined.TCR[c("P17L", "P18L", "P19L", "P20L")],
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap",
         summary.fun = "percent") # Display percentages

cat("percentGenes: Quantifying Single Gene Usage")

cat("The percentGenes() function is a specialized alias for percentGeneUsage() 
designed to quantify the usage of a single V, D, or J gene locus for a specified immune receptor chain. ")

percentGenes(combined.TCR,
             chain = "TRB",
             gene = "Vgene",
             summary.fun = "percent")

options(repr.plot.width = 4,    # width in inches
        repr.plot.height = 4,   # height in inches
        repr.plot.res = 150)    # resolution 

df.genes <- percentGenes(combined.TCR,
                         chain = "TRB",
                         gene = "Vgene",
                         exportTable = TRUE,
                         summary.fun = "proportion") 

# Performing PCA on the gene usage matrix
pc <- prcomp(t(df.genes) )

# Getting data frame to plot from
df_plot <- as.data.frame(cbind(pc$x[,1:2], colnames(df.genes)))
colnames(df_plot) <- c("PC1", "PC2", "Sample")
df_plot$PC1 <- as.numeric(df_plot$PC1)
df_plot$PC2 <- as.numeric(df_plot$PC2)

ggplot(df_plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Sample), shape = 21, size = 5) +
  guides(fill=guide_legend(title="Samples")) +
  scale_fill_manual(values = hcl.colors(nrow(df_plot), "inferno")) +
  theme_classic() +
  labs(title = "PCA of TRBV Gene Usage")

cat("percentVJ: Quantifying V-J Gene Pairings")

options(repr.plot.width = 6,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution 

percentVJ(combined.TCR[1:2],
          chain = "TRB",
          summary.fun = "percent")

df.vj <- percentVJ(combined.TCR,
                   chain = "TRB",
                   exportTable = TRUE,
                   summary.fun = "proportion") # Export proportions for PCA

# Performing PCA on the V-J pairing matrix
pc.vj <- prcomp(t(df.vj))

# Getting data frame to plot from
df_plot_vj <- as.data.frame(cbind(pc.vj$x[,1:2], colnames(df.vj)))
colnames(df_plot_vj) <- c("PC1", "PC2", "Sample")
df_plot_vj$PC1 <- as.numeric(df_plot_vj$PC1)
df_plot_vj$PC2 <- as.numeric(df_plot_vj$PC2)

# Plotting the PCA results
ggplot(df_plot_vj, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Sample), shape = 21, size = 5) +
  guides(fill=guide_legend(title="Samples")) +
  scale_fill_manual(values = hcl.colors(nrow(df_plot_vj), "inferno")) +
  theme_classic() +
  labs(title = "PCA of TRBV-TRBJ Gene Pairings")

cat("percentAA: Amino Acid Composition by Residue")

options(repr.plot.width = 8,     # width in inches
        repr.plot.height = 16,   # height in inches
        repr.plot.res = 150)     # resolution 

percentAA(combined.TCR, 
          chain = "TRB", 
          aa.length = 20)

cat("positionalEntropy: Entropy across CDR3 Sequences")

cat("

PositionalEntropy() combines the quantification by residue of percentAA() with diversity calculations. 
Positions without variance will have a value reported as 0 for the purposes of comparison.

Key Parameter(s) for positionalEntropy()

method
       shannon - Shannon Index
       inv.simpson - Inverse Simpson Index
       gini.simpson - Gini-Simpson Index
       norm.entropy - Normalized Entropy
       pielou - Pielou’s Evenness
       hill1, hill2, hill3 - Hill Numbers

")

cat("
The plot generated by positionalEntropy() illustrates the diversity or entropy at each amino acid position within the CDR3 sequence. 
Higher entropy values indicate greater variability in amino acid usage at that position, 
suggesting less selective pressure or more promiscuous binding, 
while lower values suggest conserved positions critical for structural integrity or antigen recognition.
")

options(repr.plot.width = 8,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution 


positionalEntropy(combined.TCR, 
                  chain = "TRB", 
                  aa.length = 20)

cat("positionalProperty: Amino Acid Properties across CDR3 Sequence")

# to examine the Atchley factors of amino acids across the CDR3 sequence for the first two samples :

positionalProperty(combined.TCR[c(1,2)], 
                  chain = "TRB", 
                  aa.length = 20, 
                  method = "atchleyFactors") + 
  scale_color_manual(values = hcl.colors(5, "inferno")[c(2,4)])

cat("percentKmer: Motif Quantification")

library("cowplot")
library("patchwork")

p1 = percentKmer(combined.TCR, 
            cloneCall = "aa",
            chain = "TRB", 
            motif.length = 3, 
            top.motifs = 25)

p2 = percentKmer(combined.TCR, 
            cloneCall = "nt",
            chain = "TRB", 
            motif.length = 3, 
            top.motifs = 25)

p1+p2

cat("Comparing Clonal Diversity and Overlap")

cat("clonalDiversity: Clonal Diversity Quantification")

options(repr.plot.width = 6,   # width in inches
        repr.plot.height = 4,   # height in inches
        repr.plot.res = 150)    # resolution 

combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "Patient", 
                            variables = c("P17", "P17", "P18", "P18", 
                                           "P19","P19", "P20", "P20"))

p1 = clonalDiversity(combined.TCR, 
                cloneCall = "gene")
                                           
p2 = clonalDiversity(combined.TCR, 
                cloneCall = "gene", 
                group.by = "Patient", 
                metric = "inv.simpson")

p1 + p2


cat("clonalRarefaction: Sampling-based Extrapolation")

cat("

clonalRarefaction() uses Hill numbers to estimate rarefaction, or estimating species richness, 
based on the abundance of clones across groupings. 

The underlying rarefaction calculation uses the observed receptor abundance to compute diversity.

Hill Numbers and Their Interpretation

    0: Species richness
    1: Shannon Diversity
    2: Simpson Diversity

")

cat("

By visualizing rarefaction and extrapolation curves, researchers can assess the completeness of their sampling 
and make fair comparisons of diversity, which is essential for understanding immune complexity.

")

options(repr.plot.width = 15,   # width in inches
        repr.plot.height = 4,   # height in inches
        repr.plot.res = 150)    # resolution 

p1 = clonalRarefaction(combined.TCR,
                  plot.type = 1,
                  hill.numbers = 0,
                  n.boots = 2)

p2 = clonalRarefaction(combined.TCR,
                  plot.type = 2,
                  hill.numbers = 0,
                  n.boots = 2)

p3 = clonalRarefaction(combined.TCR,
                  plot.type = 3,
                  hill.numbers = 0,
                  n.boots = 2)


p1 + p2 + p3

cat("Rarefaction using Shannon Diversity (q = 1)")

options(repr.plot.width = 15,   # width in inches
        repr.plot.height = 4,   # height in inches
        repr.plot.res = 150)    # resolution 

p1 = clonalRarefaction(combined.TCR,
                  plot.type = 1,
                  hill.numbers = 1,
                  n.boots = 2) 

p2 = clonalRarefaction(combined.TCR,
                  plot.type = 2,
                  hill.numbers = 1,
                  n.boots = 2)

p3 = clonalRarefaction(combined.TCR,
                  plot.type = 3,
                  hill.numbers = 1,
                  n.boots = 2)

p1 + p2 + p3 

cat("clonalSizeDistribution: Modeling Clonal Composition")

# By applying a spliced statistical model, it provides a more accurate representation of the repertoire’s underlying clonal architecture

clonalSizeDistribution(combined.TCR, 
                       cloneCall = "aa", 
                       method= "ward.D2")

cat("clonalOverlap: Exploring Sequence Overlap")

cat("

method Parameters for clonalOverlap()

    overlap - Overlap coefficient
    morisita - Morisita’s overlap index
    jaccard - Jaccard index
    cosine - Cosine similarity
    raw - Exact number of overlapping clones

")

p1 = clonalOverlap(combined.TCR, 
              cloneCall = "strict", 
              method = "morisita")

p2 = clonalOverlap(combined.TCR, 
              cloneCall = "strict", 
              method = "raw")

p1 + p2



cat("Combining Clones and Single-Cell Objects")

# The data in the scRepertoire package is derived from a study of acute respiratory stress disorder 
# in the context of bacterial and COVID-19 infections.

scRep_example <- get(data("scRep_example"))

#Making a Single-Cell Experiment object
sce <- Seurat::as.SingleCellExperiment(scRep_example)

#Adding patient information
scRep_example$Patient <- substr(scRep_example$orig.ident, 1,3)

#Adding type information
scRep_example$Type <- substr(scRep_example$orig.ident, 4,4)

# Note on Dimensional Reduction

# In single-cell RNA sequencing workflows, dimensional reduction is typically performed by first identifying highly variable features. 
# These features are then used directly for UMAP/tSNE projection or as inputs for principal component analysis. 
# The same approach is commonly applied to clustering as well.

# However, in immune-focused datasets, VDJ genes from TCR and BCR are often among the most variable genes. 
# This variability arises naturally due to clonal expansion and diversity within lymphocytes. 
# As a result, UMAP projections and clustering outcomes may be influenced by clonal information rather 
# than broader transcriptional differences across cell types.

# To mitigate this issue, a common strategy is to exclude VDJ genes from the set of highly variable features 
# before proceeding with clustering and dimensional reduction. 

# We introduce a set of functions that facilitate this process by removing VDJ-related genes from either 
# a Seurat Object or a vector of gene names (useful for SCE-based workflows).

# quietVDJgenes() – Removes both TCR and BCR VDJ genes.
# quietTCRgenes() – Removes only TCR VDJ genes.
# quietBCRgenes() – Removes only BCR VDJ genes, but retains BCR VDJ pseudogenes in the variable features.

# Check the first 10 variable features before removal

# Remove TCR VDJ genes

scRep_example <- quietTCRgenes(scRep_example)

library(scRepertoire)
library(Seurat)
library(SingleCellExperiment)
library(scater)


# By applying these functions, you can ensure that clustering and dimensional reduction 
# are driven by broader transcriptomic differences across cell types rather than being 
# skewed by the inherent variability due to clonal expansion

# ── Step 1: Load built-in data ──────────────────────────────
data("contig_list")
data("scRep_example")

# ── Step 2: Create combined TCR ─────────────────────────────
combined.TCR <- combineTCR(contig_list,
                           samples = c("P17B", "P17L", "P18B", "P18L",
                                       "P19B", "P19L", "P20B", "P20L"))

# ── Step 3: Add metadata to Seurat object ───────────────────
scRep_example$Patient <- substr(scRep_example$orig.ident, 1, 3)
scRep_example$Type    <- substr(scRep_example$orig.ident, 4, 4)

# ── Step 4: Remove VDJ genes before reduction ───────────────
scRep_example <- quietTCRgenes(scRep_example)

# ── Step 5: Convert Seurat → SCE ────────────────────────────
sce <- Seurat::as.SingleCellExperiment(scRep_example)

# ── Step 6: Remove ident column (causes issues) ─────────────
sce$ident <- NULL


options(repr.plot.width = 15,   # width in inches
        repr.plot.height = 4,   # height in inches
        repr.plot.res = 150)    # resolution 

# ── Step 7a: SCE version with plotUMAP ──────────────────────
sce <- combineExpression(combined.TCR,
                         sce,
                         cloneCall = "gene",
                         group.by = "sample",
                         proportion = TRUE)

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

plotUMAP(sce, colour_by = "cloneSize") +
    scale_color_manual(values = rev(colorblind_vector[c(1,3,5,7)]))

# ── Step 7b: Seurat version with DimPlot (proportion) ───────
scRep_example <- combineExpression(combined.TCR,
                                   scRep_example,
                                   cloneCall = "gene",
                                   group.by = "sample",
                                   proportion = TRUE)

DimPlot(scRep_example, group.by = "cloneSize") +
    scale_color_manual(values = rev(colorblind_vector[c(1,3,5,7)]))

# ── Step 7c: Seurat version with DimPlot (frequency) ────────
scRep_example <- combineExpression(combined.TCR,
                                   scRep_example,
                                   cloneCall = "gene",
                                   group.by = "sample",
                                   proportion = FALSE,
                                   cloneSize = c(Single = 1,
                                                 Small = 5,
                                                 Medium = 20,
                                                 Large = 100,
                                                 Hyperexpanded = 500))

DimPlot(scRep_example, group.by = "cloneSize") +
    scale_color_manual(values = rev(colorblind_vector[c(1,3,5,7)]))


# Combining both TCR and BCR

# TCR <- combineTCR(...)
# BCR <- combineBCR(...)
# list.receptors <- c(TCR, BCR)

# seurat <- combineExpression(list.receptors, 
#                            seurat, 
#                            cloneCall="gene", 
#                            proportion = TRUE)

cat("clonalOverlay")

options(repr.plot.width = 10,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution 

clonalOverlay(scRep_example, 
              reduction = "umap", 
              cutpoint = 1, 
              bins = 10, 
              facet.by = "Patient") + 
              guides(color = "none")

cat("clonalNetwork")

# ggraph needs to be loaded due to issues with ggplot
library(ggraph)

options(repr.plot.width = 12,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution

# No Identity filter
a = clonalNetwork(scRep_example, 
              reduction = "umap", 
              group.by = "seurat_clusters",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "aa")

# Examining Cluster 3 only
b = clonalNetwork(scRep_example, 
              reduction = "umap", 
              group.by = "seurat_clusters",
              filter.identity = 3,
              cloneCall = "aa")

a + b

# shared clones and highlightClones

shared.clones <- clonalNetwork(scRep_example, 
                               reduction = "umap", 
                               group.by = "seurat_clusters",
                               cloneCall = "aa", 
                               exportClones = TRUE)
head(shared.clones)

scRep_example <- highlightClones(scRep_example, 
                    cloneCall= "aa", 
                    sequence = c("CAERGSGGSYIPTF_CASSDPSGRQGPRWDTQYF", 
                                 "CARKVRDSSYKLIF_CASSDSGYNEQFF"))

cat("clonalOccupy")

options(repr.plot.width = 12,   # width in inches
        repr.plot.height = 6,   # height in inches
        repr.plot.res = 150)    # resolution


a = clonalOccupy(scRep_example, 
              x.axis = "seurat_clusters") 

b = clonalOccupy(scRep_example, 
             x.axis = "ident", 
             proportion = TRUE, 
             label = FALSE)

a + b

cat("alluvialClones")


alluvialClones(scRep_example,
               cloneCall = "aa",          # changed from clone.call → cloneCall
               y.axes = c("Patient", "ident", "Type"),
               color = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF", 
                         "NA_CASSVRRERANTGELFF")) +
    scale_fill_manual(values = c("grey", colorblind_vector[3]))

# Fix: clone.call → cloneCall
alluvialClones(scRep_example,
               cloneCall = "gene",
               y.axes = c("Patient", "ident", "Type"),
               color = "ident")

# getCirclize and vizCirclize

# getCirclize() - Generates data for use with the circlize package
# vizCirclize() - A convenient wrapper that handles the circlize plotting 

# formals(alluvialClones)

# library(circlize)
# library(scales)

# vizCirclize(scRep_example,
#            group.by = "seurat_clusters",
#            clone.call = "aa")

library(scales)
library(circlize)

circles <- getCirclize(scRep_example,
                       group.by = "seurat_clusters")

#Just assigning the normal colors to each cluster
grid.cols <- hue_pal()(length(unique(scRep_example$seurat_clusters)))
names(grid.cols) <- unique(scRep_example$seurat_clusters)

#Graphing the chord diagram
chordDiagram(circles, self.link = 1, grid.col = grid.cols)

subset <- subset(scRep_example, Type == "L")

circles <- getCirclize(subset, group.by = "ident", proportion = TRUE)

grid.cols <- scales::hue_pal()(length(unique(subset@active.ident)))
names(grid.cols) <- levels(subset@active.ident)

chordDiagram(circles,
             self.link = 1,
             grid.col = grid.cols,
             directional = 1,
             direction.type = "arrows",
             link.arr.type = "big.arrow")

cat("Quantifying Clonal Bias")

cat("

StartracDiversity

From the excellent work by Zhang et al. (2018, Nature)
The authors introduced new methods for looking at clones by cellular origins and cluster identification. 
We strongly recommend you read and cite their publication when using this function.

")

# Calculate and plot all three STARTRAC indices

StartracDiversity(scRep_example, 
                  type = "Type", 
                  group.by = "Patient")

cat("Calculating a Single Index")

# Calculate and plot only the clonal expansion index
StartracDiversity(scRep_example, 
                  type = "Type", 
                  group.by = "Patient",
                  index = "expa")

cat("Pairwise Migration Analysis")

# # Calculate pairwise migration between tissues
# StartracDiversity(scRep_example, 
#                  type = "Type", 
#                  group.by = "Patient",
#                  index = "migr",
#                  pairwise = "Type")

cat("clonalBias")

cat("
clonalBias(), like STARTRAC, is a clonal metric that seeks to quantify how individual clones are skewed towards a specific cellular compartment or cluster. 

A clone bias of 1 indicates that a clone is composed of cells from a single compartment or cluster, 
while a clone bias of 0 matches the background subtype distribution. 

")

options(repr.plot.width = 6,   # width in inches
        repr.plot.height = 4,   # height in inches
        repr.plot.res = 150)    # resolution

clonalBias(scRep_example, 
           cloneCall = "aa", 
           split.by = "Patient", 
           group.by = "seurat_clusters",
           n.boots = 10, 
           min.expand =5)

cat("Clustering by Edit Distance")

cat("

clonalCluster: Cluster by Sequence Similarity

The clonalCluster() function provides a powerful method to group clonotypes based on sequence similarity. 
It calculates the edit distance between CDR3 sequences and uses this information to build a network, 
identifying closely related clusters of T or B cell receptors.

")

# Run clustering on the first two samples for the TRA chain
sub_combined <- clonalCluster(combined.TCR[c(1,2)], 
                              chain = "TRA", 
                              sequence = "aa", 
                              threshold = 0.85)

# View the new cluster column
head(sub_combined[[1]][, c("barcode", "TCR1", "TRA.Cluster")])

#Adding patient information
scRep_example$Patient <- substr(scRep_example$orig.ident, 1,3)

#Adding type information
scRep_example$Type <- substr(scRep_example$orig.ident, 4,4)

# Run clustering, but group calculations by "Patient"
scRep_example <- clonalCluster(scRep_example, 
                               chain = "TRA", 
                               sequence = "aa", 
                               threshold = 0.85, 
                               group.by = "Patient")

#Define color palette 
num_clusters <- length(unique(na.omit(scRep_example$TRA.Cluster)))
cluster_colors <- hcl.colors(n = num_clusters, palette = "inferno")

# Bug in Seurat 5.3.1 no longer handles NA for DimPlot
# Will Update in Future Commits

# DimPlot(scRep_example, group.by = "TRA.Cluster") +
#  scale_color_manual(values = cluster_colors, na.value = "grey") +

set.seed(42)
library(igraph)

# Clustering Patient 19 samples
igraph.object <- clonalCluster(combined.TCR[c(5,6)],
                               chain = "TRB",
                               sequence = "aa",
                               group.by = "sample",
                               threshold = 0.85, 
                               exportGraph = TRUE)

#Setting color scheme
col_legend <- factor(igraph::V(igraph.object)$group)
col_samples <- hcl.colors(2,"inferno")[as.numeric(col_legend)]
color.legend <- factor(unique(igraph::V(igraph.object)$group))
sample.vertices <- V(igraph.object)[sample(length(igraph.object), 500)]

subgraph.object <- induced_subgraph(igraph.object, vids = sample.vertices)
V(subgraph.object)$degrees <- igraph::degree(subgraph.object)
edge_alpha_color <- adjustcolor("gray", alpha.f = 0.3)

#Plotting
plot(subgraph.object,
     layout = layout_nicely(subgraph.object),
     vertex.label = NA,
     vertex.size = sqrt(igraph::V(subgraph.object)$degrees), 
     vertex.color = col_samples[sample.vertices],
     vertex.frame.color = "white", 
     edge.color = edge_alpha_color,
     edge.arrow.size = 0.05,
     edge.curved = 0.05, 
     margin = -0.1)
legend("topleft", 
       legend = levels(color.legend), 
       pch = 16, 
       col = unique(col_samples), 
       bty = "n")

# Generate the sparse matrix
adj.matrix <- clonalCluster(combined.TCR[c(1,2)],
                            chain = "TRB",
                            exportAdjMatrix = TRUE)

# View the dimensions and a snippet of the matrix
dim(adj.matrix)

print(adj.matrix[1:10, 1:10])



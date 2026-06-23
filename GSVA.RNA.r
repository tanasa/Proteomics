# https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html

# The GSVA package provides the implementation of four single-sample gene set enrichment methods, 
# concretely zscore, plage, ssGSEA and its own called GSVA

library(GSVA)

cat("

Gene set variation analysis (GSVA) provides an estimate of pathway activity by transforming an input gene-by-sample expression data matrix 
into a corresponding gene-set-by-sample expression data matrix. 

This resulting expression data matrix can be then used with classical analytical methods such as differential expression, classification, 
survival analysis, clustering or correlation analysis in a pathway-centric manner. 

One can also perform sample-wise comparisons between pathways and other molecular data types such as microRNA expression or binding data, 
copy-number variation (CNV) data or single nucleotide polymorphisms (SNPs).

")

# gsva()

# two input arguments:

# a normalized gene expression dataset
# a collection of gene sets

# kcdf: The first step of the GSVA algorithm brings gene expression profiles to a common scale by calculating 
# an expression statistic through the estimation of the CDF across samples. 

# maxDiff: The last step of the GSVA algorithm calculates the gene set enrichment score from 
# two Kolmogorov-Smirnov random walk statistics

library(org.Hs.eg.db)

# goannot <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns="GO")
# head(goannot)

library(GSEABase)
library(GSVAdata)

data(c2BroadSets)
class(c2BroadSets)
c2BroadSets

names(c2BroadSets)

length(c2BroadSets)

library(GSEABase)
library(GSVA)

URL <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c7.immunesigdb.v2024.1.Hs.symbols.gmt"
c7.genesets <- readGMT(URL)

c7.genesets 

gsvaAnnotation(c7.genesets) <- SymbolIdentifier("org.Hs.eg.db")

c7.genesets 



cat("Quantification of pathway activity in bulk microarray and RNA-seq data")

# gene expression data of lymphoblastoid cell lines (LCL) from HapMap individuals

library(Biobase)

data(commonPickrellHuang)

stopifnot(identical(featureNames(huangArrayRMAnoBatchCommon_eset),
                    featureNames(pickrellCountsArgonneCQNcommon_eset)))
stopifnot(identical(sampleNames(huangArrayRMAnoBatchCommon_eset),
                    sampleNames(pickrellCountsArgonneCQNcommon_eset)))

canonicalC2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                                      grep("^REACTOME", names(c2BroadSets)),
                                      grep("^BIOCARTA", names(c2BroadSets)))]
canonicalC2BroadSets

names(canonicalC2BroadSets)

gs.list <- geneIds(c2BroadSets)

str(gs.list)



huangPar <- gsvaParam(huangArrayRMAnoBatchCommon_eset, canonicalC2BroadSets,
                      minSize=5, maxSize=500)
esmicro <- gsva(huangPar)
esmicro

# esmicro is an ExpressionSet, not a SummarizedExperiment.

library(Biobase)
score.mat <- exprs(esmicro)

class(score.mat)
dim(score.mat)

score.mat[1:10, 1:10]

esmicro[1:5, 1:5]

pickrellPar <- gsvaParam(pickrellCountsArgonneCQNcommon_eset,
                         canonicalC2BroadSets, minSize=5, maxSize=500,
                         kcdf="Poisson")
esrnaseq <- gsva(pickrellPar)
esrnaseq 

library(Biobase)

score.mat <- exprs(esrnaseq)

class(score.mat)

dim(score.mat)

score.mat[1:10, 1:10]

head(rownames(score.mat))   # pathway names
head(colnames(score.mat))   # sample names



cat("Molecular signature identification")

grep("brain|Tx|Db", ls(), value = TRUE, ignore.case = TRUE)

data(brainTxDbSets)

brainTxDbSets

# In (Verhaak et al. 2010) four subtypes of glioblastoma multiforme (GBM) -proneural, classical, neural and mesenchymal- 
# were identified by the characterization of distinct gene-level expression patterns. 

data(gbm_VerhaakEtAl)
gbm_eset

gbmPar <- gsvaParam(gbm_eset, brainTxDbSets, maxDiff=FALSE)
gbm_es <- gsva(gbmPar)

library(RColorBrewer)

subtypeOrder <- c("Proneural", "Neural", "Classical", "Mesenchymal")

sampleOrderBySubtype <- sort(match(gbm_es$subtype, subtypeOrder),
                             index.return=TRUE)$ix
subtypeXtable <- table(gbm_es$subtype)
subtypeColorLegend <- c(Proneural="red", Neural="green",
                        Classical="blue", Mesenchymal="orange")

geneSetOrder <- c("astroglia_up", "astrocytic_up", "neuronal_up",
                  "oligodendrocytic_up")

geneSetLabels <- gsub("_", " ", geneSetOrder)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
hmcol <- hmcol[length(hmcol):1]

heatmap(exprs(gbm_es)[geneSetOrder, sampleOrderBySubtype], Rowv=NA,
        Colv=NA, scale="row", margins=c(3,5), col=hmcol,
        ColSideColors=rep(subtypeColorLegend[subtypeOrder],
                          times=subtypeXtable[subtypeOrder]),
        labCol="", gbm_es$subtype[sampleOrderBySubtype],
        labRow=paste(toupper(substring(geneSetLabels, 1,1)),
                     substring(geneSetLabels, 2), sep=""),
        cexRow=2, main=" \n ")
par(xpd=TRUE)

text(0.23,1.21, "Proneural", col="red", cex=1.2)
text(0.36,1.21, "Neural", col="green", cex=1.2)
text(0.47,1.21, "Classical", col="blue", cex=1.2)
text(0.62,1.21, "Mesenchymal", col="orange", cex=1.2)
mtext("Gene sets", side=4, line=0, cex=1.5)
mtext("Samples          ", side=1, line=4, cex=1.5)



cat("Differential expression at pathway level")

data(geneprotExpCostaEtAl2021)
se <- geneExpCostaEtAl2021
se

gsvaAnnotation(se) <- EntrezIdentifier("org.Hs.eg.db")

colData(se)

table(colData(se))

cat("Filtering of immunologic gene sets")

innatepat <- c("NKCELL_VS_.+_UP", "MAST_CELL_VS_.+_UP",
               "EOSINOPHIL_VS_.+_UP", "BASOPHIL_VS_.+_UP",
               "MACROPHAGE_VS_.+_UP", "NEUTROPHIL_VS_.+_UP")
innatepat <- paste(innatepat, collapse="|")
innategsets <- names(c7.genesets)[grep(innatepat, names(c7.genesets))]
length(innategsets)


adaptivepat <- c("CD4_TCELL_VS_.+_UP", "CD8_TCELL_VS_.+_UP", "BCELL_VS_.+_UP")
adaptivepat <- paste(adaptivepat, collapse="|")
adaptivegsets <- names(c7.genesets)[grep(adaptivepat, names(c7.genesets))]
excludepat <- c("NAIVE", "LUPUS", "MYELOID")
excludepat <- paste(excludepat, collapse="|")
adaptivegsets <- adaptivegsets[-grep(excludepat, adaptivegsets)]
length(adaptivegsets)


c7.genesets.filt <- c7.genesets[c(innategsets, adaptivegsets)]
length(c7.genesets.filt)


cat("Running GSVA")

gsvapar <- gsvaParam(se, 
                     c7.genesets.filt, 
                     assay="logCPM", 
                     minSize=5,
                     maxSize=300)

es <- gsva(gsvapar)
es

assayNames(se)

assay(es)[1:3, 1:3]

# head(lapply(geneSets(es), head))

head(geneSetSizes(es))

# MDS plot of GSVA enrichment scores.

library(limma)

plotMDS(assay(es))

es$FIR

cat("Differential expression")

library(sva)
library(limma)

## build design matrix of the model to which we fit the data
mod <- model.matrix(~ FIR, colData(es))
## build design matrix of the corresponding null model
mod0 <- model.matrix(~ 1, colData(es))
## estimate surrogate variables (SVs) with SVA
sv <- sva(assay(es), mod, mod0)

## add SVs to the design matrix of the model of interest
mod <- cbind(mod, sv$sv)
## fit linear models
fit <- lmFit(assay(es), mod)
## calculate moderated t-statistics using the robust regime
fit.eb <- eBayes(fit, robust=TRUE)
## summarize the extent of differential expression at 5% FDR
res <- decideTests(fit.eb)
summary(res)

gssizes <- geneSetSizes(es)
plot(sqrt(gssizes), sqrt(fit.eb$sigma), xlab="Sqrt(gene sets sizes)",
          ylab="Sqrt(standard deviation)", las=1, pch=".", cex=4)
lines(lowess(sqrt(gssizes), sqrt(fit.eb$sigma)), col="red", lwd=2)

fit.eb.trend <- eBayes(fit, robust=TRUE, trend=sqrt(gssizes))
res <- decideTests(fit.eb.trend)
summary(res)

tt <- topTable(fit.eb.trend, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.05]
DEpwys

library(limma)
library(RColorBrewer)

## get DE pathway GSVA enrichment scores, removing the covariates effect
DEpwys_es <- removeBatchEffect(
  assay(es[DEpwys, ]),
  covariates = mod[, 2:ncol(mod), drop = FALSE],
  design = mod[, 1:2, drop = FALSE]
)

## create FIR color map
table(es$FIR)

fir_levels <- unique(es$FIR)

fircolor <- setNames(
  brewer.pal(max(3, length(fir_levels)), "Set1")[seq_along(fir_levels)],
  fir_levels
)

## cluster samples
sam_col_map <- fircolor[as.character(es$FIR)]
names(sam_col_map) <- colnames(DEpwys_es)

sampleClust <- hclust(
  as.dist(1 - cor(DEpwys_es, method = "spearman")),
  method = "complete"
)

## cluster pathways
gsetClust <- hclust(
  as.dist(1 - cor(t(DEpwys_es), method = "pearson")),
  method = "complete"
)

## annotate pathways as innate or adaptive
labrow <- rownames(DEpwys_es)

mask <- rownames(DEpwys_es) %in% innategsets
labrow[mask] <- paste("(INNATE)", labrow[mask], sep = "_")

mask <- rownames(DEpwys_es) %in% adaptivegsets
labrow[mask] <- paste("(ADAPTIVE)", labrow[mask], sep = "_")

labrow <- gsub("_", " ", gsub("GSE[0-9]+_", "", labrow))

## pathway expression color scale from blue low to red high
pwyexpcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
pwyexpcol <- rev(pwyexpcol)

## generate heatmap
heatmap(
  DEpwys_es,
  ColSideColors = sam_col_map,
  xlab = "Samples",
  ylab = "Pathways",
  margins = c(2, 20),
  labCol = "",
  labRow = labrow,
  col = pwyexpcol,
  scale = "row",
  Colv = as.dendrogram(sampleClust),
  Rowv = as.dendrogram(gsetClust)
)

# res <- igsva()

# using proteomics data :
# https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA_proteomics.html






# load libraries

library(Seurat)
library(dplyr)
library(Matrix)
library(abind)
source(paste0(.libPaths(), "/MarcusFuncs/statsFunction.R"))

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read dge files into list and setup seurat object

dge <- list.files("/Users/basiri/LHA/objects/dge/individual/dge", pattern="*.txt", full.names=TRUE)
dge.list <- lapply(dge, read.table, sep="\t", header=TRUE, row.names=1)
names(dge.list) <- paste0(substr(dge, 123, nchar(dge)-38), ".data")
dge.list <- dge.list[order(names(dge.list))]

dge.seurat.list <- list()
for (i in seq_along(dge.list)){ 
dge.seurat.list[[i]] <- new("seurat", raw.data = dge.list[[i]])}
names(dge.seurat.list) <- substr(names(dge.list), 1, 4)

seurat.list <- list()
for (i in seq_along(dge.seurat.list)){  
seurat.list[[i]] <- Setup(dge.seurat.list[[i]], min.cells = 0, min.genes = 0, do.logNormalize = F, project = paste0(substr(names(dge.seurat.list[i]), 1, 4)), do.scale=F, do.center=F)}
names(seurat.list) <- substr(names(dge.list), 1, 4)

# merge seurat count object

counts <- MergeSeurat(seurat.list[[1]], seurat.list[[2]], project = "control.highfat.counts", min.cells = 10, min.genes = 0, do.logNormalize = F, names.field = 1, add.cell.id1 = paste0(substr(names(dge.seurat.list[1]), 1, 4)), add.cell.id2 = paste0(substr(names(dge.seurat.list[2]), 1, 4)), do.scale=F, do.center=F)
for (i in 3:length(seurat.list)){counts <- MergeSeurat(counts, seurat.list[[i]], project = "control.highfat.counts", min.cells = 10, min.genes = 0, do.logNormalize = F, names.field = 1, add.cell.id2 = paste0(substr(names(dge.seurat.list[i]), 1, 4)), do.scale=F, do.center=F)}

# add metadata-counts
# percentMito metadata
mitoGenes.counts <- grep("^mt-", rownames(counts@data), value = T)
percentMito.counts <- colSums(counts@data[mitoGenes.counts, ])/colSums(counts@data)
counts <- AddMetaData(counts, percentMito.counts, "percentMito")

# grep cells
cells.list <- list()
for (i in seq_along(dge.seurat.list)){ 
cells.list[[i]] <- grep(paste0(names(seurat.list[i])), rownames(counts@data.info), value = T)}
names(cells.list) <- paste0(substr(dge, 123, nchar(dge)-38), ".cells")

# group metadata
group.list <- list()
for (i in 1:length(seurat.list)){
    if (i %% 2 == 1) group.list[[i]] <- array(rep("control",length(cells.list[[i]])))
    if (i %% 2 != 1) group.list[[i]] <- array(rep("highfat",length(cells.list[[i]])))}

names(group.list) <- paste0(substr(dge, 123, nchar(dge)-38), ".group")
for (i in seq_along(seurat.list)){dimnames(group.list[[i]]) <- list(cells.list[[i]])}      
group <- do.call("abind", group.list)
counts <- AddMetaData(counts, group, "group")
        
# pool metadata
pool.list <- list()
MLB003 = c(1,2,3,5,6,8,13,14)
MLB004 = c(4,7,9,10,11,12)

for (i in 1:length(seurat.list)){
    if (i %in% MLB003) pool.list[[i]] <- array(rep("MLB003",length(cells.list[[i]])))
    if (i %in% MLB004) pool.list[[i]] <- array(rep("MLB004",length(cells.list[[i]])))}

names(pool.list) <- paste0(substr(dge, 123, nchar(dge)-38), ".pool")
for (i in seq_along(seurat.list)){dimnames(pool.list[[i]]) <- list(cells.list[[i]])}       
pool <- do.call("abind", pool.list)
counts <- AddMetaData(counts, pool, "pool")
        
# batch metadata
batch.list <- list()
batch.A <- list("1" = c(1,2), "2" = c(3,4), "3" = c(5,6), "4" = c(7,8), "5" = c(9,10), "6" = c(11,12), "7" = c(13,14))

for (i in 1:length(seurat.list)){for(j in 1:length(batch.A)){
    if (i %in% batch.A[[j]]) batch.list[[i]] <- array(rep(paste0(names(batch.A[j])),length(cells.list[[i]])))}}

names(batch.list) <- paste0(substr(dge, 123, nchar(dge)-38), ".batch.A")
for (i in seq_along(seurat.list)){dimnames(batch.list[[i]]) <- list(cells.list[[i]])}       
batch.A <- do.call("abind", batch.list)
counts <- AddMetaData(counts, batch.A, "batch.A")
                                 
head(counts@data.info,3)
stats(counts)

# subset data-counts
counts <- SubsetData(counts, subset.name = "nGene", accept.high = 5000, accept.low = 200)
counts <- SubsetData(counts, subset.name = "percentMito", accept.high = 0.10, accept.low = 0.001)
counts <- SubsetData(counts, subset.name = "nUMI", accept.high = 15000)

stats(counts)

# count histograms -- Figure S1B

pdf("/Users/basiri/LHA/figs/subset/hist/180109_MLB003MLB004_control.highfat.counts.hist.pdf")
par(mfrow=c(1,2), lwd=0.2)
hist(as.matrix(counts@data.info$nGene), breaks="FD", col="darkred", xlab = "total counts", ylab = "frequency", lty="blank", xlim=c(0,max(counts@data.info$nGene)))
abline(v = median(counts@data.info$nGene), col = "black", lwd = 1, lty=3)
hist(as.matrix(counts@data.info$nUMI), breaks="FD", col="steelblue", xlab = "total counts", ylab = "frequency", lty="blank", xlim=c(0,max(counts@data.info$nUMI)))
abline(v = median(counts@data.info$nUMI), col = "black", lwd = 1, lty=3)
dev.off()

# merge seurat-logMed

med = median(as.matrix(counts@data.info$nUMI))
logMed <- MergeSeurat(seurat.list[[1]], seurat.list[[2]], project = "logMed", min.cells = 10, min.genes = 0, do.logNormalize = T, total.expr = med, names.field = 1, add.cell.id1 = paste0(substr(names(dge.seurat.list[1]), 1, 4)), add.cell.id2 = paste0(substr(names(dge.seurat.list[2]), 1, 4)), do.scale=T, do.center=T)
for (i in 3:length(seurat.list)){logMed <- MergeSeurat(logMed, seurat.list[[i]], project = "logMed", min.cells = 10, min.genes = 0, do.logNormalize = T, total.expr = med, names.field = 1, add.cell.id2 = paste0(substr(names(dge.seurat.list[i]), 1, 4)), do.scale=T, do.center=T)}

# add metadata-logMed
# percentMito metadata
mitoGenes.logMed <- grep("^mt-", rownames(logMed@data), value = T)
percentMito.logMed <- colSums(expm1(logMed@data[mitoGenes.logMed, ]))/colSums(expm1(logMed@data))
logMed <- AddMetaData(logMed, percentMito.logMed, "percentMito")

# grep cells
cells.list <- list()
for (i in seq_along(dge.seurat.list)){ 
cells.list[[i]] <- grep(paste0(names(seurat.list[i])), rownames(logMed@data.info), value = T)}
names(cells.list) <- paste0(substr(dge, 123, nchar(dge)-38), ".cells")

# group metadata
group.list <- list()
for (i in 1:length(seurat.list)){
    if (i %% 2 == 1) group.list[[i]] <- array(rep("control",length(cells.list[[i]])))
    if (i %% 2 != 1) group.list[[i]] <- array(rep("highfat",length(cells.list[[i]])))}
names(group.list) <- paste0(substr(dge, 123, nchar(dge)-38), ".group")

        for (i in seq_along(seurat.list)){dimnames(group.list[[i]]) <- list(cells.list[[i]])}        
group <- do.call("abind", group.list)
logMed <- AddMetaData(logMed, group, "group")
 
# pool metadata
pool.list <- list()
MLB003 = c(1,2,3,5,6,8,13,14)
MLB004 = c(4,7,9,10,11,12)

for (i in 1:length(seurat.list)){
    if (i %in% MLB003) pool.list[[i]] <- array(rep("MLB003",length(cells.list[[i]])))
    if (i %in% MLB004) pool.list[[i]] <- array(rep("MLB004",length(cells.list[[i]])))}
names(pool.list) <- paste0(substr(dge, 123, nchar(dge)-38), ".pool")

        for (i in seq_along(seurat.list)){dimnames(pool.list[[i]]) <- list(cells.list[[i]])}       
pool <- do.call("abind", pool.list)
logMed <- AddMetaData(logMed, pool, "pool")
        
# batch metadata
batch.list <- list()
batch.A <- list("1" = c(1,2), "2" = c(3,4), "3" = c(5,6), "4" = c(7,8), "5" = c(9,10), "6" = c(11,12), "7" = c(13,14))

for (i in 1:length(seurat.list)){for(j in 1:length(batch.A)){
    if (i %in% batch.A[[j]]) batch.list[[i]] <- array(rep(paste0(names(batch.A[j])),length(cells.list[[i]])))}}
names(batch.list) <- paste0(substr(dge, 123, nchar(dge)-38), ".batch.A")
                                 
for (i in seq_along(seurat.list)){dimnames(batch.list[[i]]) <- list(cells.list[[i]])}       
batch.A <- do.call("abind", batch.list)
logMed <- AddMetaData(logMed, batch.A, "batch.A")
                                 
head(logMed@data.info,3)
stats(logMed)

# subset data and save-logMed
logMed <- SubsetData(logMed, subset.name = "nGene", accept.high = 5000, accept.low = 200)
logMed <- SubsetData(logMed, subset.name = "percentMito", accept.high = 0.10, accept.low = 0.001)
logMed <- SubsetData(logMed, subset.name = "nUMI", accept.high = 15000)
cat("total number of cells: ", nrow(logMed@data.info),"\n")

saveRDS(logMed, "/Users/basiri/LHA/objects/subset/180109_MLB003MLB004_control.highfat.logMed.rds")
saveRDS(counts, "/Users/basiri/LHA/objects/subset/180109_MLB003MLB004_control.highfat.counts.rds")

# load libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(abind)
library(reshape2)
library(sva)
library(scRNA.seq.funcs)
library(GGally)

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read files in
logMed <- readRDS("/Users/basiri/LHA/objects/subset/180109_MLB003MLB004_control.highfat.logMed.rds")

# remove zero variance genes from @data

data <- as.matrix(logMed@data)
logMed.var <- apply(data, 1, var)
zeroVarIndex <- which(logMed.var %in% 0)
zeroVarGenes <- rownames(data[c(zeroVarIndex),])

cat("number of zero variance genes: ", length(zeroVarGenes),"\n")

data <- data[!rownames(data) %in% zeroVarGenes, ]
logMed@data <- as(data, "dgCMatrix")

# ComBat

# ComBat Model Matrix
logMed.combatPar.mod <- model.matrix(~group + pool + as.numeric(nGene) + as.numeric(nUMI) + as.numeric(percentMito), data=logMed@data.info)

# ComBat parametric batch correction
logMed.combatPar <- ComBat(dat=logMed@data, batch=logMed@data.info$batch.A, mod=logMed.combatPar.mod, par.prior=T, prior.plots=F)
logMed.combatPar <- as.matrix(logMed.combatPar)

# save batch corrected matrix to @scale.data
logMed@scale.data <- as(logMed.combatPar, "matrix")
saveRDS(logMed, "/Users/basiri/LHA/objects/batchCorrection/combat/180112_MLB003MLB004_control.highfat.rds")

# assess batch correction -- Figure S1C and S1D

origIdent <- sort(unique(logMed@data.info$orig.ident))
j = seq_along(origIdent)

cells.list <- list()
for (i in j){ 
cells.list[[i]] <- grep(paste0(origIdent[i]), rownames(logMed@data.info), value = T)}
names(cells.list) <- paste0(origIdent, ".cells")

# raw.means
raw <- as.matrix(logMed@raw.data)
raw.means.list <- list()
for (i in j){raw.means.list[[i]] <- rowMeans(raw[, cells.list[[i]]])}
names(raw.means.list) <- paste0(origIdent, ".cells")

control.raw.means.list <- list()
highfat.raw.means.list <- list()
for (i in j){ 
    if (i %% 2 == 1) control.raw.means.list[[(i-1)/2+1]] <- raw.means.list[[i]]
    if (i %% 2 != 1) highfat.raw.means.list[[(i-2)/2+1]] <- raw.means.list[[i]]}
names(control.raw.means.list) <- origIdent[c(j %% 2 == 1)]
names(highfat.raw.means.list) <- origIdent[c(j %% 2 != 1)]
        
# log.means
log <- as.matrix(logMed@data)
log.means.list <- list()
for (i in j){log.means.list[[i]] <- rowMeans(log[, cells.list[[i]]])}
names(log.means.list) <- paste0(origIdent, ".cells")

control.log.means.list <- list()
highfat.log.means.list <- list()
for (i in j){ 
    if (i %% 2 == 1) control.log.means.list[[(i-1)/2+1]] <- log.means.list[[i]]
    if (i %% 2 != 1) highfat.log.means.list[[(i-2)/2+1]] <- log.means.list[[i]]}
names(control.log.means.list) <- origIdent[c(j %% 2 == 1)]
names(highfat.log.means.list) <- origIdent[c(j %% 2 != 1)]
        
# corrected.means
corrected <- as.matrix(logMed@scale.data)
corrected.means.list <- list()
for (i in j){corrected.means.list[[i]] <- rowMeans(corrected[, cells.list[[i]]])}
names(corrected.means.list) <- paste0(origIdent, ".cells")

control.corrected.means.list <- list()
highfat.corrected.means.list <- list()
for (i in j){ 
    if (i %% 2 == 1) control.corrected.means.list[[(i-1)/2+1]] <- corrected.means.list[[i]]
    if (i %% 2 != 1) highfat.corrected.means.list[[(i-2)/2+1]] <- corrected.means.list[[i]]}
names(control.corrected.means.list) <- origIdent[c(j %% 2 == 1)]
names(highfat.corrected.means.list) <- origIdent[c(j %% 2 != 1)]
        
# calculate and plot rle
log.rle <- calc_cell_RLE(expm1(log))
corrected.rle <- calc_cell_RLE(expm1(corrected))
rle <- melt(cbind(log = log.rle, corrected = corrected.rle))[,2:3]
names(rle) <- c("group", "rle")

ggplot(rle, aes(x=group, y=rle, fill=group)) + geom_violin(scale="width", width=1, trim=F, adjust = 3) + stat_summary(fun.y=median, geom="point", shape=21, size=2, stroke=1, fill="black", color="white") + theme(legend.position="none")
ggsave("/Users/basiri/LHA/figs/batchCorrection/rle/180112_MLB003MLB004_control.highfat.RLE.pdf")

# plot mean expression Pearsons R across animals
pdf("/Users/basiri/LHA/figs/batchCorrection/ggcorr/180112_MLB003MLB004_control.highfat.ggcorr.pdf")
ggcorr(control.raw.means.list, hjust = 0.75, size = 4, low = "lightgrey", high = "#115f70", layout.exp = 1, midpoint = 0.84, limits = c(0.68, 1), label = T, label_round = 2, label_size = 5, legend.position = "top", name="control.raw")
ggcorr(control.corrected.means.list, hjust = 0.75, size = 4, low = "lightgrey", high = "#115f70", layout.exp = 1, midpoint = 0.84, limits = c(0.68, 1), label = T, label_round = 2, label_size = 5, legend.position = "top", name="control.corrected")
ggcorr(highfat.raw.means.list, hjust = 0.75, size = 4, low = "lightgrey", high = "#115f70", layout.exp = 1, midpoint = 0.84, limits = c(0.68, 1), label = T, label_round = 2, label_size = 5, legend.position = "top", name="highfat.raw")
ggcorr(highfat.corrected.means.list, hjust = 0.75, size = 4, low = "lightgrey", high = "#115f70", layout.exp = 1, midpoint = 0.84, limits = c(0.68, 1), label = T, label_round = 2, label_size = 5, legend.position = "top", name="highfat.corrected")
dev.off()

# load libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(abind)
library(M3Drop)

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read files in
logMed <- readRDS("/Users/basiri/LHA/objects/batchCorrection/combat/180112_MLB003MLB004_control.highfat.rds")

# identify commGenes

dge <- list.files("/Users/basiri/LHA/objects/dge/individual/dge", pattern="*.txt", full.names=TRUE)
dge.list <- lapply(dge, read.table, sep="\t", header=TRUE, row.names=1)

# make vector of matching genes across animals in each group
genes.list <- list()
for (i in seq_along(dge.list)){genes.list[[i]] <- rownames(as.vector(dge.list[[i]]))}
names(genes.list) <- paste0(substr(dge, 123, nchar(dge)-38), ".genes")
list2env(genes.list ,.GlobalEnv)

control.commGenes <- Reduce(intersect, list(AJP1.genes, ALP1.genes, AOP1.genes, AQP1.genes, AZP1.genes, BBP1.genes, BDP1.genes))
highfat.commGenes <- Reduce(intersect, list(AKP1.genes, ANP1.genes, APP1.genes, ASP1.genes, BAP1.genes, BCP1.genes, BEP1.genes))
cat("number of control.commGenes: ", length(control.commGenes),"\n")
cat("number of highfat.commGenes: ", length(highfat.commGenes),"\n")

# combine unique commGenes by group
commGenes <- unique(abind(control.commGenes, highfat.commGenes))
commGenes <- Reduce(intersect, list(rownames(logMed@scale.data), commGenes))
cat("number of unique commGenes control.highfat: ", length(commGenes),"\n")

logMed@var.genes <- as(commGenes, "character")
saveRDS(logMed, "/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.rds")
saveRDS(commGenes, "/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.genes.rds")

# identify BrenneckeHVG

logMed.data <- as.matrix(logMed@data)

pdf("/Users/basiri/LHA/figs/genes/BrenneckeHVG/180112_MLB003MLB004_control.highfat.HVG.pdf")
logMed.BrenneckeHVG <- BrenneckeGetVariableGenes(logMed.data, fdr = 0.01, minBiolDisp = 0.5)
dev.off()

logMed.commHVG <- Reduce(intersect, list(commGenes, logMed.BrenneckeHVG))

cat("number of BrenneckeHVG: ", length(logMed.BrenneckeHVG),"\n")
cat("number of commHVG: ", length(logMed.commHVG),"\n")

logMed@var.genes <- as(logMed.commHVG, "character")
saveRDS(logMed, "/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.rds")

# run PCA

library(Seurat)
library(dplyr)
library(Matrix)

commHVG <- readRDS("/proj/stuberlb/Users/marcus/dropseq/data/MLB003-MLB004_170222/seurat/objects/genes/180112_MLB003MLB004_control.highfat.rds")
commHVG <- PCA(commHVG, pc.genes=commHVG@var.genes, pcs.store = 120, replace.pc = TRUE)
commHVG <- ProjectPCA(commHVG, pcs.store = 120, do.center=TRUE)

saveRDS(commHVG, "/proj/stuberlb/Users/marcus/dropseq/data/MLB003-MLB004_170222/seurat/objects/pca/180112_MLB003MLB004_control.highfat.pca.rds")

# load libraries
library(Seurat)

JonSnow <- colorRampPalette(c("#6be9ec", "black", "#FFD173"))(n = 100)

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read files in
pca <- readRDS("/Users/basiri/LHA/objects/pca/180112_MLB003MLB004_control.highfat.pca.rds")

# inspect plots, determine number of PCs for clustering

# pcaPlot
pdf("/Users/basiri/LHA/figs/pca/pcaPlot/180112_MLB003MLB004_control.highfat.pcaPlot.pdf")
for (i in 1:119){PCAPlot(pca, i, i+1, pt.size = 0.25, group.by = "orig.ident")}
dev.off()

# pcHeatmap
options(warn=-1)
pdf("/Users/basiri/LHA/figs/pca/pcHeatmap/180112_MLB003MLB004_control.highfat.pcHeatmap.pdf")
try(for (i in 1:22){PCHeatmap(pca, pc.use = 9*(i-1)+(1:9), cells.use = 500, do.balanced = T, use.scale = T, label.columns = F, col.use = JonSnow, use.full = T, disp.min = -2, disp.max = 2, cexRow = 0.45)}, silent=T)
dev.off()
options(warn=0)   

# pcElbow
pdf("/Users/basiri/LHA/figs/pca/pcElbow/180112_MLB003MLB004_control.highfat.pcElbow.pdf")
for (i in 1:10){pcElbow <- PCElbowPlot(pca, num.pc = i*50)
print(pcElbow + geom_point(aes(fill = "black"),show.legend = NA) + theme_minimal() + theme(legend.position="none"))}
dev.off()

# run tsne

library(Seurat)
library(dplyr)
library(Matrix)

dims = 100
res = 0.5

filename <- paste0("180112_MLB003MLB004_control.highfat.tsne",dims,".perp30.res",res,".rds")

tsne <- readRDS("/proj/stuberlb/Users/marcus/dropseq/data/MLB003-MLB004_170222/seurat/objects/pca/180112_MLB003MLB004_control.highfat.pca.rds")
tsne <- FindClusters(tsne, pc.use = 1:dims, resolution = res, print.output = T, save.SNN = T, algorithm = 2, k.param = 30)
tsne <- RunTSNE(tsne, dims.use = 1:dims, perplexity = 30)
saveRDS(tsne, paste0("/proj/stuberlb/Users/marcus/dropseq/data/MLB003-MLB004_170222/seurat/objects/tsne/",filename))

#load libraries
library(Seurat)
library(reshape2)

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read objects

tsne <- readRDS("/Users/basiri/LHA/objects/tsne/180112_MLB003MLB004_control.highfat.rds")
commGenes <- readRDS("/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.genes.rds")

# buildClusterTree

tsne <- BuildClusterTree(tsne, genes.use = commGenes, do.plot = T, do.reorder = T, reorder.numeric = T)
saveRDS(tsne, "/Users/basiri/LHA/objects/tree/180117_MLB003MLB004_control.highfat.tree.rds")


pdf("/Users/basiri/LHA/figs/tree/180117_MLB003MLB004_control.highfat.clusterTree.pdf")
PlotClusterTree(tsne)
dev.off()

# visualize tsnePlot and features -- Figure 1B and S1E

clust.colors <- c("#a42f2b", "#e4773a", "#ccaf47", "#97bd53", "#749d60", "#3d573c", "#70bb8e", "#3f877f", "gray", "#3480b0", "#2a2d6b", "#bea9cd", "#b874aa", "#ac5399", "#65315a")

pdf("/Users/basiri/LHA/figs/tsne/tsnePlot/180117_MLB003MLB004_control.highfat.tree.pdf")
TSNEPlot(tsne, do.label = F, pt.size = 0.6, color=clust.colors)
TSNEPlot(tsne, do.label = F, pt.size = 0.6, group.by='group')
TSNEPlot(tsne, do.label = F, pt.size = 0.6, group.by='orig.ident')
dev.off()

# markers FeaturePlots
markers <- c("Stmn2", "Scg2", "Camk2b", "Slc17a6", "Gad1", "Hcrt", "Pmch", "Mag", "Flt1", "C1qb", "Cacng4", "Agt", "Acta2")
for (i in seq_along(markers)){
svg(paste0("/Users/basiri/LHA/figs/tree/180117_MLB003MLB004_control.highfat.tree.",markers[i],".FeaturePlot.svg"))
FeaturePlot(tsne, c(markers[i]),cols.use = c("#bed3d3","#e04b2a"), pt.size = 1, nCol=2)
dev.off()}

# cluster num cells, mean nUMI, mean percentMito tsne

for(i in 1:max(tsne@ident)){cat('num cells cluster ',i,':', nrow(subset(tsne@data.info, tree.ident == i)), "\n")}
cat("\n")
for(i in 1:max(tsne@ident)){cat('mean nUMI cluster ',i,':', mean(subset(tsne@data.info, tree.ident == i)$nUMI), "\n")}
cat("\n")
for(i in 1:max(tsne@ident)){cat('mean nGene cluster ',i,':', mean(subset(tsne@data.info, tree.ident == i)$nGene), "\n")}
cat("\n")
for(i in 1:max(tsne@ident)){cat('mean percentMito cluster ',i,':', mean(subset(tsne@data.info, tree.ident == i)$percentMito), "\n")}

# plot nUMI and nGene  data

allCells.info <- as.data.frame(tsne@data.info)
allCells.info$tree.ident <- as.factor(allCells.info$tree.ident)

ggplot(allCells.info, aes(x=tree.ident, y=nGene, fill=tree.ident)) + geom_violin(width=1, scale="width") + stat_summary(fun.y=median, geom="point", shape=21, size=2, stroke=1, fill="black", color="white") + theme(legend.position="none") + scale_fill_manual(values=as.vector(clust.colors))
ggsave("/Users/basiri/LHA/figs/subset/vlnPlot/180117_MLB003MLB004_control.highfat.tree.nGene.vlnPlot.svg", width = 12, height = 4, dpi = 300)

ggplot(allCells.info, aes(x=tree.ident, y=nUMI, fill=tree.ident)) + geom_violin(width=1, scale="width") + stat_summary(fun.y=median, geom="point", shape=21, size=2, stroke=1, fill="black", color="white") + theme(legend.position="none") + scale_fill_manual(values=as.vector(clust.colors))
ggsave("/Users/basiri/LHA/figs/subset/vlnPlot/180117_MLB003MLB004_control.highfat.tree.nUMI.vlnPlot.svg", width = 12, height = 4, dpi = 300)

ggplot(allCells.info, aes(x=tree.ident, y=percentMito, fill=tree.ident)) + geom_violin(width=1, scale="width") + stat_summary(fun.y=median, geom="point", shape=21, size=2, stroke=1, fill="black", color="white") + theme(legend.position="none") + scale_fill_manual(values=as.vector(clust.colors))
ggsave("/Users/basiri/LHA/figs/subset/vlnPlot/180117_MLB003MLB004_control.highfat.tree.pctMito.vlnPlot.svg", width = 12, height = 4, dpi = 300)

point.colors <-  c(rep("#BD3786FF", 12), rep("grey50", 7)) 
ggplot(allCells.info, aes(x=nUMI, y=nGene, fill=tree.ident)) + geom_point(shape=21, size=0.75, stroke=0.2, color="black") + scale_fill_manual(values=as.vector(clust.colors)) + theme_minimal() + theme(legend.position="none") 
ggsave("/Users/basiri/LHA/figs/subset/genePlot/180117_MLB003MLB004_control.highfat.tree.nUMIxnGene.vlnPlot.svg", width = 5, height = 4, dpi = 300)

# percent of cells in each cluster and group

ident <- list()
for(i in 1:max(tsne@ident)){ ident[i] <-nrow(subset(tsne@data.info, tree.ident == i))}
ident <- melt(as.matrix(ident))
ident$value <- as.numeric(ident$value)
ident <- cbind(ident, pct = ident$value/sum(ident$value)*100)

# pct.highfat
control.highfat.c <- c(NA)
control.highfat.c.list <- list()
for (i in 1:max(tsne@ident)){
    control.highfat.c <- SubsetData(tsne, ident.use = c(i), do.center = F, do.scale = F)
        assign(paste0('control.highfat.c',i), control.highfat.c)
        control.highfat.c.list[[i]] <- control.highfat.c                         
        setNames(control.highfat.c.list, paste0('control.highfat.c',i))}

control.highfat.c.group <- c(NA)
control.highfat.c.group.list <- list()
for (i in 1:max(tsne@ident)){
    control.highfat.c.group <- SetAllIdent(control.highfat.c.list[[i]], "group")
        assign(paste0('control.highfat.c',i,'.group'), control.highfat.c.group)
        control.highfat.c.group.list[[i]] <- control.highfat.c.group                         
        setNames(control.highfat.c.group.list, paste0('control.highfat.c',i,'.group'))}

highfat <- list()
for (i in 1:max(tsne@ident)){(highfat[[i]] <- as.vector(control.highfat.c.group.list[[i]]@ident))}

pct.highfat <- list()
for (i in 1:max(tsne@ident)){pct.highfat[[i]] <- sum(highfat[[i]]=='highfat')/length(highfat[[i]])*100}

ident <- cbind(ident, pct.highfat = as.matrix(unlist(pct.highfat)))

ggplot(as.data.frame(ident), aes(Var1, pct, fill=Var1)) + geom_bar(stat="identity", col = "white", fill=clust.colors)
ggsave("/Users/basiri/LHA/figs/tsne/pct/180117_MLB003MLB004_control.highfat.oobeMerge.tree.clusterPct.pdf")

ggplot(as.data.frame(ident), aes(Var1, pct.highfat, fill=Var1)) + geom_bar(stat="identity", col = "white", fill=clust.colors)
ggsave("/Users/basiri/LHA/figs/tsne/pct/180117_MLB003MLB004_control.highfat.oobeMerge.tree.pctHighfat.pdf")

# load libraries

library(Seurat)
library(Matrix)
library(dplyr)

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read files in

tsne <- readRDS("/Users/basiri/LHA/objects/tree/180117_MLB003MLB004_control.highfat.tree.rds")
commGenes <- readRDS("/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.genes.rds")

# find markersClustV

markersClust.1v2 <- FindMarkers(tsne, 1, c(2), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.2v1 <- FindMarkers(tsne, 2, c(1), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.3vN25 <- FindMarkers(tsne, 3, c(4:6), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.4vN26 <- FindMarkers(tsne, 4, c(5:6), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.5v6 <- FindMarkers(tsne, 5, c(6), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.6v5 <- FindMarkers(tsne, 6, c(5), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.7v8 <- FindMarkers(tsne, 7, c(8), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.8v7 <- FindMarkers(tsne, 8, c(7), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.9vN29 <- FindMarkers(tsne, 9, c(10:11), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.10v11 <- FindMarkers(tsne, 10, c(11), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.11v10 <- FindMarkers(tsne, 11, c(10), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.12v13 <- FindMarkers(tsne, 12, c(13), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.13v12 <- FindMarkers(tsne, 13, c(12), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.14v15 <- FindMarkers(tsne, 14, c(15), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
markersClust.15v14 <- FindMarkers(tsne, 15, c(14), genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)

markersClust.1v2 <- cbind(markersClust.1v2, gene=rownames(markersClust.1v2))
markersClust.2v1 <- cbind(markersClust.2v1, gene=rownames(markersClust.2v1))
markersClust.3vN25 <- cbind(markersClust.3vN25, gene=rownames(markersClust.3vN25))
markersClust.4vN26 <- cbind(markersClust.4vN26, gene=rownames(markersClust.4vN26))
markersClust.5v6 <- cbind(markersClust.5v6, gene=rownames(markersClust.5v6))
markersClust.6v5 <- cbind(markersClust.6v5, gene=rownames(markersClust.6v5))
markersClust.7v8 <- cbind(markersClust.7v8, gene=rownames(markersClust.7v8))
markersClust.8v7 <- cbind(markersClust.8v7, gene=rownames(markersClust.8v7))
markersClust.9vN29 <- cbind(markersClust.9vN29, gene=rownames(markersClust.9vN29))
markersClust.10v11 <- cbind(markersClust.10v11, gene=rownames(markersClust.10v11))
markersClust.11v10 <- cbind(markersClust.11v10, gene=rownames(markersClust.11v10))
markersClust.12v13 <- cbind(markersClust.12v13, gene=rownames(markersClust.12v13))
markersClust.13v12 <- cbind(markersClust.13v12, gene=rownames(markersClust.13v12))
markersClust.14v15 <- cbind(markersClust.14v15, gene=rownames(markersClust.14v15))
markersClust.15v14 <- cbind(markersClust.15v14, gene=rownames(markersClust.15v14))

saveRDS(markersClust.1v2, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.1v2.rds")
saveRDS(markersClust.2v1, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.2v1.rds")
saveRDS(markersClust.3vN25, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.3vN25.rds")
saveRDS(markersClust.4vN26, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.4vN26.rds")
saveRDS(markersClust.5v6, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.5v6.rds")
saveRDS(markersClust.6v5, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.6v5.rds")
saveRDS(markersClust.7v8, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.7v8.rds")
saveRDS(markersClust.8v7, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.8v7.rds")
saveRDS(markersClust.9vN29, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.9vN29.rds")
saveRDS(markersClust.10v11, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.10v11.rds")
saveRDS(markersClust.11v10, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.11v10.rds")
saveRDS(markersClust.12v13, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.12v13.rds")
saveRDS(markersClust.13v12, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.13v12.rds")
saveRDS(markersClust.14v15, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.14v15.rds")
saveRDS(markersClust.15v14, "/Users/basiri/LHA/objects/clusterMarkers/allCells/180117_MLB003MLB004_control.highfat.tree.markersClust.15v14.rds")

# load libraries

library(Seurat)
library(Matrix)
library(dplyr)
library(abind)
library(reshape2)

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read objects and markersClustV

allCells.tsne <- readRDS("/Users/basiri/LHA/objects/tree/180117_MLB003MLB004_control.highfat.tree.rds")
commGenes <- readRDS("/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.genes.rds")

# read markersClustV into list

markersClustV <- list.files("/Users/basiri/LHA/objects/clusterMarkers/allCells/individual/", pattern="*.rds", full.names=TRUE)
markersClustV.list <- lapply(markersClustV, readRDS)
names(markersClustV.list) <- substr(markersClustV, 279, nchar(markersClustV)-4)
markersClustV.list <- markersClustV.list[c('markersClust.1v2', 'markersClust.2v1', 'markersClust.3vN25', 'markersClust.4vN26', 'markersClust.5v6', 'markersClust.6v5', 'markersClust.7v8', 'markersClust.8v7', 'markersClust.9vN29', 'markersClust.10v11', 'markersClust.11v10', 'markersClust.12v13', 'markersClust.13v12', 'markersClust.14v15', 'markersClust.15v14')]

# create matrix of pVals for all clusters

markersClustV.pVal.list <- lapply(markersClustV.list, function(x) {x[c("gene", "p_val")]})

for (i in seq_along(markersClustV.pVal.list)){
  colnames(markersClustV.pVal.list[[i]]) <- c("gene", paste0("c",i, ".pVal"))}

markersClustV.pVal <- Reduce(function(x,y)merge(x,y,by="gene", all=TRUE), markersClustV.pVal.list)

row.names(markersClustV.pVal) <- markersClustV.pVal[,1]
markersClustV.pVal <- markersClustV.pVal[,-1]

cat("number of genes markersClustV.pVal: ",nrow(markersClustV.pVal),"\n")

saveRDS(markersClustV.pVal, "/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pVal/180117_MLB003MLB004_control.highfat.tree.markersClustV.pVal.rds")
write.table(markersClustV.pVal, "/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pVal/180117_MLB003MLB004_control.highfat.tree.markersClustV.pVal.txt", sep="\t", col.names=NA)

# create cluster average expression matrix

clust.avg <- AverageExpression(allCells.tsne, genes.use = commGenes, return.seurat = F)
colnames(clust.avg) <- c(paste0("c",seq_along(colnames(clust.avg)), ".avg"))

cat("number of genes clust.avg:", nrow(clust.avg), "\n")

saveRDS(clust.avg, "/Users/basiri/LHA/objects/clusterMarkers/allCells/master/clust.avg/180117_MLB003MLB004_control.highfat.tree.clust.avg.rds")
write.table(clust.avg, "/Users/basiri/LHA/objects/clusterMarkers/allCells/master/clust.avg/180117_MLB003MLB004_control.highfat.tree.clust.avg.txt", sep="\t", col.names=NA)

# create matrix of percent positive cells per cluster -- Supplementary Data S1

rawData <- as.matrix(allCells.tsne@raw.data)

# create list of cell names for each cluster
allCells.ident <- as.matrix(allCells.tsne@ident)
cells <- c(NA)
cells.list <- list()
for(i in 1:max(allCells.tsne@ident)){cells <- names(allCells.ident[allCells.ident[,1]==i,])
                assign(paste0('c',i,'.cells'), cells)
        cells.list[[i]] <- cells                         
        setNames(cells.list, paste0('c',i,'.cells'))}

# slice rawData by cluster
cells.rawData <- c(NA)
cells.rawData.list <- list()
for(i in 1:max(allCells.tsne@ident)){cells.rawData <- rawData[,cells.list[[i]]]
                assign(paste0('c',i,'.cells.rawData'), cells.rawData)
        cells.rawData.list[[i]] <- cells.rawData                         
        setNames(cells.rawData.list, paste0('c',i,'.cells.rawData'))}

# sums by gene per cluster
sums <- c(NA)
sums.list <- list()
for(i in 1:max(allCells.tsne@ident)){sums <- rowSums(cells.rawData.list[[i]] != 0)
                assign(paste0('c',i,'.sums'), sums)
        sums.list[[i]] <- sums                         
        setNames(sums.list, paste0('c',i,'.sums'))}

# pct by gene per cluster
pct <- c(NA)
pct.list <- list()
for(i in 1:max(allCells.tsne@ident)){pct <- sums.list[[i]]/length(cells.list[[i]])
                assign(paste0('c',i,'.pct'), pct)
        pct.list[[i]] <- pct                         
        setNames(pct.list, paste0('c',i,'.pct'))}

pct <- data.frame(sapply(pct.list,c))

for(i in 1:max(allCells.tsne@ident)){names(pct)[i] <- c(paste0('c',i,'.pct'))}

saveRDS(pct, "/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pct/180117_MLB003MLB004_control.highfat.tree.pct.rds")
write.table(pct, "/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pct/180117_MLB003MLB004_control.highfat.tree.pct.txt", sep="\t", col.names=NA)

# load libraries

library(Seurat)
library(Matrix)
library(dplyr)
library(matrixStats)
library(made4)
library(reshape2)
library(viridis)
library(RColorBrewer)

RdBu <- brewer.pal(11, "RdBu")
RdBu <- c(colorRampPalette(c(RdBu[1], RdBu[6]))(51), colorRampPalette(c(RdBu[6], RdBu[11]))(51)[-1])

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read files in

allCells <- readRDS("/Users/basiri/LHA/objects/tree/180117_MLB003MLB004_control.highfat.tree.rds")
commGenes <- readRDS("/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.genes.rds")

clust.avg <- readRDS("/Users/basiri/LHA/objects/clusterMarkers/allCells/master/clust.avg/180117_MLB003MLB004_control.highfat.tree.clust.avg.rds")
pct <- readRDS("/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pct/180117_MLB003MLB004_control.highfat.tree.pct.rds")
markersClustV.pVal <- readRDS("/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pVal/180117_MLB003MLB004_control.highfat.tree.markersClustV.pVal.rds")

clust.avg <- as.matrix(clust.avg[,c(1:8,10:15)])
pct <- as.matrix(pct[,c(1:8,10:15)])
markersClustV.pVal <- as.matrix(markersClustV.pVal[,c(1:8,10:15)])

# subset genes by pct.min, pVal.max, avg.min
pct.min = 0.1
pVal.max = 10e-3
avg.min = 0
avg.mult.mean = 2

pct.list <- list()
for (i in seq_along(colnames(pct))){pct.list[[i]] <- names(subset(pct[,i], pct[,i] >= pct.min))}

avg.list <- list()
for (i in seq_along(colnames(clust.avg))){avg.list[[i]] <- names(subset(expm1(clust.avg[,i]), expm1(clust.avg[,i]) >= avg.min))}

sd.list <- list()
clust.avg.stats <- cbind(expm1(clust.avg), sd.lim = avg.mult.mean*rowSds(expm1(clust.avg))+rowMeans(expm1(clust.avg)))
for (i in seq_along(colnames(clust.avg))){sd.list[[i]] <- names(subset(clust.avg.stats[,i], clust.avg.stats[,i] >= clust.avg.stats[,ncol(clust.avg.stats)]))}

pVal.list <- list()
for (i in seq_along(colnames(markersClustV.pVal))){pVal.list[[i]] <- names(subset(markersClustV.pVal[,i], markersClustV.pVal[,i] <= pVal.max))}

intersect.list <- list()
for (i in seq_along(pVal.list)){intersect.list[[i]] <- Reduce(intersect, c(pct.list[i], avg.list[i], pVal.list[i], sd.list[i]))}

genes <- unique(unlist(intersect.list))
clust.avg.subs <- clust.avg[genes,]

cat("number of genes in subset:", length(genes))

# clusterMarkers heatmaps -- Figures 1D and 1E

labels <- c("Slc32a1", "Hcrt", "Pmch", "Opalin", "Flt1", "Cx3cr1", "Acta2", "Spp1", "S100b", "Lyz2", "Pdgfra", "Vtn", "Agt", "Slc17a6")
labels <- intersect(genes, labels)
labelsIndx <- match(labels, rownames(clust.avg.subs))
rowlabels <- rep("", nrow(clust.avg.subs))
rowlabels[labelsIndx] <- labels
rownames(clust.avg.subs) <- rowlabels

row.col <- rep(NA, nrow(clust.avg.subs))
row.col[labelsIndx] <- "#000000"

clust.colors <-  c(magma(ncol(clust.avg)), alpha = 1)[match(colnames(clust.avg), unique(colnames(clust.avg)))]

pdf("/Users/basiri/LHA/figs/clusterMarkers/allCells/heatmap/180202_MLB003MLB004_control.highfat.oobeMerge.tree.markers.heatmap.pdf")
heatmap.2(clust.avg.subs, Rowv=F, Colv=F, dendrogram="none", ColSideColors=clust.colors, RowSideColors=row.col, col=RdBu, trace="none", scale="row", key=T, cexRow=1, labCol=colnames(clust.avg))
dev.off()

genes <- c("Gad1","Gad2", "Slc32a1","Slc17a6", "Hcrt", "Pmch")
clust.avg.subs <- clust.avg[genes,9:12]
pdf("/Users/basiri/LHA/figs/clusterMarkers/allCells/heatmap/180202_MLB003MLB004_control.highfat.oobeMerge.tree.markers.heatmap.Neurons.pdf")
heatmap.2(clust.avg.subs, Rowv=F, Colv=F, dendrogram="none", col=rev(RdBu), trace="none", scale="row", key=T, cexRow=1, cexCol=1, labCol=colnames(clust.avg.subs))
dev.off()


# tsne plot by multiple features -- Figure 1C

features <- c("Slc17a6", "Slc32a1", "Hcrt", "Pmch", "Opalin", "Flt1", "Cx3cr1", "Acta2", "Spp1", "S100b", "Lyz2", "Pdgfra", "Vtn", "Agt")

raw.data <- as.matrix(allCells@raw.data)
feature.names <- list()
for (i in 1:length(features)){feature.names[[i]] <- names(subset(raw.data[unlist(features[i]),], raw.data[unlist(features[i]),] !=0))}

tsne.rot <- allCells@tsne.rot
feature.tsne.rot <- list()
for (i in 1:length(features)){feature.tsne.rot[[i]] <- na.omit(as.matrix(tsne.rot[feature.names[[i]],]))}

scale.data <- allCells@data
feature.data <- list()
for (i in 1:length(features)){feature.data[[i]] <- as.matrix(scale.data[features[i],])
                              feature.data[[i]] <- subset(feature.data[[i]], rownames(feature.data[[i]]) %in% feature.names[[i]])}

feature.data.bind <- list()
for (i in 1:length(features)){feature.data.bind[[i]] <- cbind(gene = rep(features[i], nrow(feature.data[[i]])), feature.tsne.rot[[i]], data = as.numeric(feature.data[[i]]), scale = as.numeric(feature.data[[i]]/max(feature.data[[i]])))}

data.bind <- as.data.frame(do.call("rbind", feature.data.bind))
data.bind[, 2:ncol(data.bind)] <- lapply(2:ncol(data.bind), function(x) as.numeric(as.character(data.bind[[x]])))
tsne.rot.diff <- tsne.rot[setdiff(rownames(tsne.rot), rownames(data.bind)),]
data.bind <- rbind(data.bind, cBind(gene = rep("other", nrow(tsne.rot.diff)), tsne.rot.diff, data = max(data.bind$data), scale = max(data.bind$scale)))

ggplot(data.bind, aes(x=tSNE_1,y=tSNE_2)) + geom_point(aes(color=gene, alpha=scale), shape = 20, size=2)
ggsave("/Users/basiri/LHA/figs/clusterMarkers/allCells/featurePlot/180202_MLB003MLB004_control.highfat.oobeMerge.tree.clusterMarkers.featurePlot.distinct.pdf")

# load libraries

library(Seurat)
library(Matrix)
library(dplyr)
library(abind)
library(matrixStats)
library(reshape2)

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read objects

allCells <- readRDS("/Users/basiri/LHA/objects/tree/180117_MLB003MLB004_control.highfat.tree.rds")
commGenes <- readRDS("/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.genes.rds")

pct <- readRDS("/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pct/180117_MLB003MLB004_control.highfat.tree.pct.rds")


# subset clusters and add group names

control.highfat.c <- c(NA)
control.highfat.c.list <- list()
for (i in 1:max(allCells@ident)){
    control.highfat.c <- SubsetData(allCells, ident.use = c(i), do.center = F, do.scale = F)
        assign(paste0('control.highfat.c',i), control.highfat.c)
        control.highfat.c.list[[i]] <- control.highfat.c                         
        setNames(control.highfat.c.list, paste0('control.highfat.c',i))}

control.highfat.c.group <- c(NA)
control.highfat.c.group.list <- list()
for (i in 1:max(allCells@ident)){
    control.highfat.c.group <- SetAllIdent(control.highfat.c.list[[i]], "group")
        assign(paste0('control.highfat.c',i,'.group'), control.highfat.c.group)
        control.highfat.c.group.list[[i]] <- control.highfat.c.group                         
        setNames(control.highfat.c.group.list, paste0('control.highfat.c',i,'.group'))}

# find groupMarkers by cluster -- Supplementary Data S1
c.groupMarkers <- c(NA)
c.groupMarkers.list <- list()
for (i in 1:max(allCells@ident)){
    c.groupMarkers <- FindMarkers(control.highfat.c.group.list[[i]], 'highfat', 'control', genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
        assign(paste0('c',i,'.groupMarkers'), c.groupMarkers)
        c.groupMarkers.list[[i]] <- c.groupMarkers                         
        setNames(c.groupMarkers.list, paste0('c.',i,'.groupMarkers'))}

# save groupMarkers
for (i in 1:max(allCells@ident)){
saveRDS(c.groupMarkers.list[[i]], paste0('/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/individual/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.c',i,'.rds'))}

for (i in 1:max(allCells@ident)){
write.table(c.groupMarkers.list[[i]], paste0('/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/individual/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.c',i,'.txt'), sep="\t", col.names=NA)}

# create and save groupMarkers pVal master file
c.groupMarkers.list.subs <- lapply(c.groupMarkers.list, function(x)cbind(gene = as.character(row.names(x)),x))
c.groupMarkers.list.subs <- lapply(c.groupMarkers.list.subs, function(x)x[,1:2])
groupMarkers.subs <- suppressWarnings(Reduce(function(x,y) merge(x,y, by = 'gene', all = TRUE), c.groupMarkers.list.subs))
groupMarkers.subs <- as.data.frame(groupMarkers.subs[,-1], row.names=as.character(groupMarkers.subs[,1]))

names <- c(NA)
names.list <- list()
for (i in 1:max(allCells@ident)){names <- paste0('c', rep(i),'.pVal')
    assign(paste0('c',i,'.pVal'), names)
        names.list[[i]] <- names                        
        setNames(names.list, paste0('c',i,'.pVal'))}

colnames(groupMarkers.subs)<- names.list[1:max(allCells@ident)]
                            
saveRDS(groupMarkers.subs, '/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.master.pVal.rds')
write.table(groupMarkers.subs, '/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.master.pVal.txt', sep="\t", col.names=NA)

# find asinhFC by cluster -- Supplementary Data S1

c.avg <- c(NA)
c.avg.list <- list()
for (i in 1:max(allCells@ident)){
    c.avg <- AverageExpression(control.highfat.c.group.list[[i]], genes.use = commGenes)
        assign(paste0('c',i,'.avg'), c.avg)
        c.avg.list[[i]] <- c.avg                         
        setNames(c.avg.list, paste0('c.',i,'.avg'))}

# add asinhFC to clustAvg
c.avg.asinhFC <- c(NA)
c.avg.asinhFC.list <- list()
for (i in 1:max(allCells@ident)){
    c.avg.asinhFC <- cbind(c.avg.list[[i]], asinhFC=(asinh(expm1(c.avg.list[[i]]$highfat))-asinh(expm1(c.avg.list[[i]]$control))))
        assign(paste0('c',i,'.avg.asinhFC'), c.avg.asinhFC)
        c.avg.asinhFC.list[[i]] <- c.avg.asinhFC                        
        setNames(c.avg.asinhFC.list, paste0('c.',i,'.avg.asinhFC'))}

# save clustAvg individual matrices
for (i in 1:max(allCells@ident)){
saveRDS(c.avg.asinhFC.list[[i]], paste0('/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/individual/rds/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.clustAvg.c',i,'.rds'))}

for (i in 1:max(allCells@ident)){
write.table(c.avg.asinhFC.list[[i]], paste0('/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/individual/txt/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.clustAvg.c',i,'.txt'), sep="\t", col.names=NA)}

# merge asinhFC for all clusters into a single data.frame and save
c.clustAvg.asinhFC <- do.call("cbind", lapply(c.avg.asinhFC.list, function(x) (x[,3])))
c.clustAvg.asinhFC <- cbind(rownames(c.avg.asinhFC.list[[5]]), c.clustAvg.asinhFC )
c.clustAvg.asinhFC <- as.data.frame(c.clustAvg.asinhFC[,-1], row.names=c.clustAvg.asinhFC[,1])
    
names <- c(NA)
names.list <- list()
for (i in 1:max(allCells@ident)){names <- paste0('c', rep(i),'.asinhFC')
    assign(paste0('c',i,'.asinhFC'), names)
        names.list[[i]] <- names                        
        setNames(names.list, paste0('c',i,'.asinhFC'))}
    
colnames(c.clustAvg.asinhFC)<- names.list[1:15]
c.clustAvg.asinhFC.indx <- sapply(c.clustAvg.asinhFC, is.factor)
c.clustAvg.asinhFC[c.clustAvg.asinhFC.indx] <- lapply(c.clustAvg.asinhFC[c.clustAvg.asinhFC.indx], function(x) as.numeric(as.character(x)))
    
saveRDS(c.clustAvg.asinhFC, '/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.clustAvg.master.asinhFC.rds') 
write.table(c.clustAvg.asinhFC, '/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.clustAvg.master.asinhFC.txt', sep="\t", col.names=NA)
                                                      

# load libraries

library(Seurat)
library(Matrix)
library(dplyr)
library(abind)
library(matrixStats)
library(reshape2)
library(viridis)
library(made4)
library(RColorBrewer)

RdBu <- brewer.pal(11, "RdBu")
RdBu <- c(colorRampPalette(c(RdBu[1], RdBu[6]))(51), colorRampPalette(c(RdBu[6], RdBu[11]))(51)[-1])

absMax <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
absMin <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read objects

allCells <- readRDS("/Users/basiri/LHA/objects/tree/180117_MLB003MLB004_control.highfat.tree.rds")
commGenes <- readRDS("/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.genes.rds")

groupMarkers <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.master.pVal.rds')
pct <- readRDS("/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pct/180117_MLB003MLB004_control.highfat.tree.pct.rds")
asinhFC <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.clustAvg.master.asinhFC.rds') 

asinhFC <- as.matrix(asinhFC[,c(1:8,10:15)])
pct <- as.matrix(pct[,c(1:8,10:15)])
groupMarkers <- as.matrix(groupMarkers[,c(1:8,10:15)])

# pct sig by cluster -- Figure 2B

pct.min = 0.5
pVal.max = 10e-5

groupMarkers <- as.matrix(groupMarkers)
pct <- as.matrix(pct)

groupMarkers.subs <- list()
for (i in 1:ncol(groupMarkers)) {groupMarkers.subs[[i]] <- names(subset(groupMarkers[,i], groupMarkers[,i] <= pVal.max))}

pct.subs <- list()
for (i in 1:ncol(pct)) {pct.subs[[i]] <- names(subset(pct[,i], pct[,i] >= pct.min))}

intersect <- list()
for (i in 1:ncol(pct)) {intersect[[i]] <- intersect(unlist(pct.subs[[i]]), unlist(groupMarkers.subs[[i]]))}

pctSig <- list()
for (i in 1:ncol(groupMarkers)) {pctSig[[i]] <- (((length(na.omit(intersect[[i]])))/(length(na.omit(groupMarkers[,i]))))*100)}
pctSig <- as.data.frame(cbind(cluster = 1:ncol(groupMarkers), pct = unlist(pctSig)))

write.table(pctSig, paste0('/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/pct/190130_MLB003MLB004_control.highfat.tree.bimod.clustAvg.pVal',pVal.max,'.pctSig',pct.min,'.subs.txt'), sep="\t", col.names=NA)

clust.colors <- c("#a42f2b", "#e4773a", "#ccaf47", "#97bd53", "#749d60", "#3d573c", "#70bb8e", "#3f877f", "#3480b0", "#2a2d6b", "#bea9cd", "#b874aa", "#ac5399", "#65315a")

ggplot(pctSig, aes(x=cluster, y=pct, fill=factor(cluster))) + geom_bar(stat="identity", color="black", size=0.3) + scale_fill_manual(values=as.character(clust.colors)) + theme(legend.position="none") + labs(title=paste0("pct.min = ",pct.min, "\n pVal.max = ",pVal.max))
ggsave(paste0("/Users/basiri/LHA/figs/groupMarkers/allCells/pVals/pct/190130_MLB003MLB004_control.highfat.tree.bimod.clustAvg.pVal",pVal.max,".pctSig",pct.min,".subs.pdf"), width=6, height=4)

# calculate snr of each gene by cell diff per clust -- Figure 2A and Supplementary Data S1

pct.min = 0
pVal.max = 10e-4
asinhFC.min = 0

pct.list <- list()
for (i in seq_along(colnames(pct))){pct.list[[i]] <- names(subset(pct[,i], pct[,i] >= pct.min))}

asinhFC.list <- list()
for (i in seq_along(colnames(asinhFC))){asinhFC.list[[i]] <- names(subset(asinhFC[,i], asinhFC[,i] >= asinhFC.min))}

pVal.list <- list()
for (i in seq_along(colnames(groupMarkers))){pVal.list[[i]] <- names(subset(groupMarkers[,i], groupMarkers[,i] <= pVal.max))}

intersect.list <- list()
for (i in seq_along(pVal.list)){intersect.list[[i]] <- Reduce(intersect, c(pct.list[i], asinhFC.list[i], pVal.list[i]))}

genes <- unique(unlist(intersect.list))
asinhFC.subs <- asinhFC[genes,]

data <- as.matrix(expm1(allCells@scale.data))
control.cells <- row.names(subset(allCells@data.info, allCells@data.info$group == "control"))
highfat.cells <- row.names(subset(allCells@data.info, allCells@data.info$group == "highfat"))

clust.cells.list <- list()
for (i in c(1:8,10:15)){clust.cells.list[[i]] <- row.names(subset(allCells@data.info, allCells@data.info$tree.ident == i))}
clust.cells.list <- Filter(Negate(is.null), clust.cells.list)

clust.control.data.list <- list()
for (i in seq_along(clust.cells.list)){clust.control.data.list[[i]] <- asinh(data[intersect.list[[i]],intersect(clust.cells.list[[i]], control.cells)])}

clust.highfat.data.list <- list()
for (i in seq_along(clust.cells.list)){clust.highfat.data.list[[i]] <- asinh(data[intersect.list[[i]],intersect(clust.cells.list[[i]], highfat.cells)])}

clust.diff.list <- list()
for (i in seq_along(clust.cells.list)){clust.diff.list[[i]] <- clust.highfat.data.list[[i]]-rowMedians(clust.control.data.list[[i]])}

clust.snr.list <- list()
for (i in seq_along(clust.cells.list)){clust.snr.list[[i]] <- cbind(cluster=rep(i, nrow(clust.diff.list[[i]])), snr=(rowMeans(clust.control.data.list[[i]]))/(rowSds(clust.diff.list[[i]])))}  

clust.snr <- as.data.frame(do.call("rbind", clust.snr.list))
clust.snr <- subset(clust.snr, abs(clust.snr$snr) < 2)

ggplot(clust.snr, aes(x=cluster, y=snr, color=factor(cluster))) + geom_jitter(size = 0.1, width=0.2) + scale_color_manual(values=as.character(clust.colors)) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", color = "#000000", width = 0.65, size=0.25) + theme(legend.position="none") + labs(subtitle = paste0("pct.min = ",pct.min, "\n pVal.max = ",pVal.max))
ggsave("/Users/basiri/LHA/figs/groupMarkers/allCells/asinhFC/snr/190202_MLB003MLB004_control.highfat.tree.bimod.asinhFC.snr.pdf", width=6, height=4)


# create matrix with TRUE/FALSE by significance per cluster

pVal.subs <- groupMarkers[genes,]
pVal.subs[pVal.subs <= pVal.max] <- "TRUE"
pVal.subs[pVal.subs != "TRUE"] <- "FALSE"

colnames(pVal.subs) <- c("MG1", "MG2", "EOC", "End", "Peri", "VSM", "Olig1", "Olig2", "Vgat", "Vglut2", "Orx", "Mch", "Astro", "OPC")

write.table(pVal.subs, paste0('/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/summary/190423_MLB003MLB004_control.highfat.pVal.summary.txt'), sep="\t", col.names=NA)


# make master asinhFC.subs heatmap -- Figure S4A

pct.min = 0.1
pVal.max = 10e-4
delta_asinhFC = 0.05

pct.list <- list()
for (i in seq_along(colnames(pct))){pct.list[[i]] <- names(subset(pct[,i], pct[,i] >= pct.min))}

asinhFC.list <- list()
for (i in seq_along(colnames(asinhFC))){asinhFC.list[[i]] <- names(subset(asinhFC[,i], abs(asinhFC[,i]) >= delta_asinhFC))}

pVal.list <- list()
for (i in seq_along(colnames(groupMarkers))){pVal.list[[i]] <- names(subset(groupMarkers[,i], groupMarkers[,i] <= pVal.max))}

intersect.list <- list()
for (i in seq_along(pVal.list)){intersect.list[[i]] <- Reduce(intersect, c(pct.list[i], asinhFC.list[i], pVal.list[i]))}

genes <- unique(unlist(intersect.list))
asinhFC.subs <- asinhFC[genes,]

cat("number of genes in subset:", length(genes))

pdf("/Users/basiri/LHA/figs/groupMarkers/allCells/asinhFC/heatmap/190117_MLB003MLB004_control.highfat.tree.bimod.pVal.asinhFC.subs.heatmap.pdf")
heatmap.2(as.matrix(asinhFC.subs), Rowv=F, Colv=F, dendrogram="none", ColSideColors=clust.colors, col=RdBu, trace="none", scale="col", key=F, labRow = F, cexCol=0.5, labCol=colnames(asinhFC.subs), breaks=seq(-2, 2, length.out=102), lhei = c(1,8), lwid = c(1,1))
dev.off()

pdf("/Users/basiri/LHA/figs/groupMarkers/allCells/asinhFC/heatmap/190117_MLB003MLB004_control.highfat.tree.bimod.pVal.asinhFC.subs.heatmapKey.pdf")
heatmap.2(as.matrix(asinhFC.subs), Rowv=F, Colv=F, dendrogram="none", ColSideColors=clust.colors, col=RdBu, trace="none", scale="col", key=T, labRow = F, cexCol=0.5, labCol=colnames(asinhFC.subs), breaks=seq(-2, 2, length.out=102))
dev.off()

# asinhFC correlation matrix

cor.asinhFC <- cor(asinhFC)
cor.asinhFC[cor.asinhFC == 1] <- NA

orangeBlueDark100 <- colorRampPalette(c("#215872", "#FFFFFF", "#995c2c"))(100)

pdf("/Users/basiri/LHA/figs/groupMarkers/allCells/asinhFC/corr/190117_MLB003MLB004_control.highfat.tree.bimod.pVal.asinhFC.corr.pdf")
suppressWarnings(heatmap.2(cor.asinhFC, Rowv=F, Colv=F, dendrogram="none", col=rev(orangeBlueDark100), trace="none", na.color = "#102b38", scale="none", cexRow=0.5, cexCol=0.5, key=T))
dev.off()

# CDF plot of pVals in neurons -- Figure 2C

pct.min = 0.5
pVal.max = 10e-2

groupMarkers <- as.matrix(groupMarkers)
pct <- as.matrix(pct)

groupMarkers.subs <- list()
for (i in 1:ncol(groupMarkers)) {groupMarkers.subs[[i]] <- names(subset(groupMarkers[,i], groupMarkers[,i] <= pVal.max))}

pct.subs <- list()
for (i in 1:ncol(pct)) {pct.subs[[i]] <- names(subset(pct[,i], pct[,i] >= pct.min))}

intersect <- list()
for (i in 1:ncol(pct)) {intersect[[i]] <- intersect(unlist(pct.subs[[i]]), unlist(groupMarkers.subs[[i]]))}

neurons = c(9:12)

groupMarkers.subs.melt <- c(NA)
groupMarkers.subs.melt.list <- list()
for (i in neurons){groupMarkers.subs.melt.list[[i]] <- cbind(melt(na.omit(groupMarkers[intersect[[i]],i])), cluster = rep(i))}
groupMarkers.subs.melt.list <- Reduce(rbind, groupMarkers.subs.melt.list)

ecdf.data <- ggplot(groupMarkers.subs.melt.list, aes(value, color = as.factor(cluster))) + stat_ecdf(geom = "step") + scale_color_manual(values=as.character(clust.colors[neurons]))
ecdf.data <- as.matrix(ecdf.data$data)

ecdf.data.gaba <- as.numeric(subset(ecdf.data, ecdf.data[,"cluster"] == 9)[,"value"])
ecdf.data.glut <- subset(ecdf.data, ecdf.data[,"cluster"] == 10)[,"value"]
ecdf.data.orx <- subset(ecdf.data, ecdf.data[,"cluster"] == 11)[,"value"]
ecdf.data.pmch <- subset(ecdf.data, ecdf.data[,"cluster"] == 12)[,"value"]

# ks test vglut cluster vs other neuronal clusters
ecdf.data.gaba <- ks.test(ecdf.data.glut, ecdf.data.gaba, alternative = c("two.sided"), tol = 1e-20, exact = T, B = 10000)
ecdf.data.orx <- ks.test(ecdf.data.glut, ecdf.data.orx, alternative = c("two.sided"), tol = 1e-20, exact = T, B = 10000)
ecdf.data.pmch <- ks.test(ecdf.data.glut, ecdf.data.pmch, alternative = c("two.sided"), tol = 1e-20, exact = T, B = 10000)

ecdf.data.gaba
ecdf.data.orx
ecdf.data.pmch

groupMarkers.subs.melt.list[groupMarkers.subs.melt.list == 9] <- "GABA"
groupMarkers.subs.melt.list[groupMarkers.subs.melt.list == 10] <- "Glut"
groupMarkers.subs.melt.list[groupMarkers.subs.melt.list == 11] <- "Orx"
groupMarkers.subs.melt.list[groupMarkers.subs.melt.list == 12] <- "Mch"

write.table(groupMarkers.subs.melt.list, paste0('/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/cdf/190130_MLB003MLB004_control.highfat.tree.bimod.clustAvg.pVal",pVal.max,".pctSig",pct.min,".subs.CDF.txt'), sep="\t", col.names=NA)

ggplot(groupMarkers.subs.melt.list, aes(value, color = as.factor(cluster))) + stat_ecdf(geom = "step", size = 0.75) + scale_color_manual(values=as.character(clust.colors[neurons])) + labs(title="Empirical Cumulative Density Function", subtitle = paste0("pct.min = ",pct.min, "\n pVal.max = ",pVal.max)) + xlab("p-value") + ylab("CDF")
ggsave(paste0("/Users/basiri/LHA/figs/groupMarkers/allCells/pVals/cdf/190130_MLB003MLB004_control.highfat.tree.bimod.clustAvg.pVal",pVal.max,".pctSig",pct.min,".subs.CDF.pdf"), width=6, height=4)

# plot cluster 11 asinhFC Figure 2D

asinhFC.c11 <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/individual/rds/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.clustAvg.c11.rds')
pVal.c11 <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/individual/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.c11.rds')

intersect <- intersect(rownames(asinhFC.c11), rownames(pVal.c11))
asinhFC.c11.subs <- asinhFC.c11[intersect,]
pVal.c11.subs <- pVal.c11[intersect,]
c11.subs <- cbind(asinhFC.c11.subs, pVal= pVal.c11.subs[,1])
c11.subs <- subset(c11.subs, abs(c11.subs$asinhFC) <=0.1)
c11.subs <- subset(c11.subs, abs(c11.subs$pVal) >=10e-20)
c11.subs <- cbind(c11.subs, color=sign(c11.subs$asinhFC))

ggplot(as.data.frame(c11.subs), aes(x=(asinhFC), y=-log10(pVal))) + geom_point(aes(color=abs(asinhFC)),size=0.5, color = "#000000") + geom_hline(yintercept=-log10(0.01), color = "red")
ggsave("/Users/basiri/LHA/figs/groupMarkers/allCells/asinhFC/scatter/180121_MLB003MLB004_control.highfat.tree.bimod.pVal.asinhFC.scatter.pdf", width = 12, height = 8, dpi = 300)

# load libraries

library(Seurat)
library(Matrix)
library(dplyr)
library(abind)
library(RColorBrewer)
library(monocle)
library(tibble)
library(scico)

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read objects

allCells <- readRDS("/Users/basiri/LHA/objects/tree/180117_MLB003MLB004_control.highfat.tree.rds")
commGenes <- readRDS("/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.genes.rds")

groupMarkers <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.master.pVal.rds')
pct <- readRDS("/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pct/180117_MLB003MLB004_control.highfat.tree.pct.rds")
asinhFC <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.clustAvg.master.asinhFC.rds') 

asinhFC <- as.matrix(asinhFC[,c(1:8,10:15)])
pct <- as.matrix(pct[,c(1:8,10:15)])
groupMarkers <- as.matrix(groupMarkers[,c(1:8,10:15)])

# run monocle on all vglut cells

pct.min = 0.10
qval.min = 0.01

clust11 <- SubsetData(allCells, ident.use = c(11), do.center = T, do.scale = T)
clust11 <- SubsetData(clust11, subset.name = "nGene", accept.high = 5000, accept.low = 500)
clust11 <- SubsetData(clust11, subset.name = "percentMito", accept.high = 0.05, accept.low = 0.001)
clust11 <- SubsetData(clust11, subset.name = "nUMI", accept.high = 15000)

clust11@data.info$tree.ident <- as.character(clust11@data.info$tree.ident)

clust11.monocle <- newCellDataSet(as.matrix(clust11@raw.data[commGenes, colnames(clust11@data)]), phenoData = AnnotatedDataFrame(clust11@data.info), expressionFamily = negbinomial.size())
clust11.monocle <- estimateSizeFactors(clust11.monocle)
clust11.monocle <- estimateDispersions(clust11.monocle)

diff_test_res <- differentialGeneTest(clust11.monocle, fullModelFormulaStr = "~group")
ordering_genes <- intersect(rownames(subset(pct, pct[,"c11.pct"] >= pct.min)), row.names(subset(diff_test_res, qval <= qval.min)))
cat("number of genes in ordering subset:", length(ordering_genes), "\n")

clust11.monocle <- setOrderingFilter(clust11.monocle, ordering_genes)

clust11.monocle <- reduceDimension(clust11.monocle, max_components = 2, method = 'DDRTree')
clust11.monocle <- orderCells(clust11.monocle)

saveRDS(clust11.monocle, '/Users/basiri/LHA/objects/monocle/190211_MLB003MLB004_control.highfat.clust11Subs.tree.monocle.rds')

# plot trajectory and pseudotimes -- Figure 2F and 2G

colorBy <- c("group", "tree.ident", "Pseudotime", "State")
for (i in seq_along(colorBy)){
plot_cell_trajectory(clust11.monocle, color_by = colorBy[i])
ggsave(paste0("/Users/basiri/LHA/figs/monocle/trajectory/190211_MLB003MLB004_control.highfat.clust11Subs.tree.monocle.",colorBy[i],".trajectory.pdf"), dpi = 300, limitsize = F)}

plot_cell_trajectory(clust11.monocle, color_by = "group", cell_size = 0.75) + facet_wrap(~group, nrow = 1)
ggsave(paste0("/Users/basiri/LHA/figs/monocle/trajectory/190211_MLB003MLB004_control.highfat.clust11Subs.tree.monocle.groupFacet.trajectory.pdf"), dpi = 300, height = 8, width = 12)


# plot pesudotime mst

sample_state <- pData(clust11.monocle)$State
lib_info_with_pseudo <- pData(clust11.monocle)

reduced_dim_coords <- reducedDimK(clust11.monocle)

  ica_space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() %>%
    select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
    mutate(sample_name = rownames(.), sample_state = rownames(.))

 edges <- minSpanningTree(clust11.monocle) %>%
    igraph::as_data_frame() %>%
    select_(source = "from", target = "to") %>%
    left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
    left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")

 data_df <- t(monocle::reducedDimS(clust11.monocle)) %>%
    as.data.frame() %>%
    select_(data_dim_1 = 1, data_dim_2 = 2) %>%
    rownames_to_column("sample_name") %>%
    mutate(sample_state) %>%
    left_join(lib_info_with_pseudo %>% rownames_to_column("sample_name"), by = "sample_name")

theta = 0
return_rotation_mat <- function(theta) {theta <- theta / 180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)}

data_df[, c("data_dim_1", "data_dim_2")] <- as.matrix(data_df[, c("data_dim_1", "data_dim_2")]) %*% t(return_rotation_mat(theta))
edges[, c("source_prin_graph_dim_1", "source_prin_graph_dim_2")] <- as.matrix(edges[, c("source_prin_graph_dim_1", "source_prin_graph_dim_2")]) %*% t(return_rotation_mat(theta))
edges[, c("target_prin_graph_dim_1", "target_prin_graph_dim_2")] <- as.matrix(edges[, c("target_prin_graph_dim_1", "target_prin_graph_dim_2")]) %*% t(return_rotation_mat(theta))

ptime <- as.data.frame(pData(clust11.monocle))
ptime <- ptime[match(edges$target, rownames(ptime)),]

edges <- cbind(edges, ptime = ptime$Pseudotime)

ggplot(edges) + geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2", color = "ptime"), size=2, data=edges) + scale_colour_gradientn(colors = rev(scico(1000, begin = 0.05, end = 1, direction = 1, palette = "grayC")))
ggsave(paste0("/Users/basiri/LHA/figs/monocle/mst/190211_MLB003MLB004_control.highfat.clust11Subs.tree.monocle.trajectory.mst.pdf"), dpi = 300)

# rug plots
colorBy <- c("group", "tree.ident", "State")

for (i in seq_along(colorBy)){
pdata <- pData(clust11.monocle)
pdata <- pdata[order(pdata$Pseudotime),]
pdata <- cbind(pdata[c(paste0(colorBy[i]), "Pseudotime")], order = rep(seq_along(1:nrow(pdata))))
ggplot(pdata, aes_string("order", colour=paste0(colorBy[i]), group=paste0(colorBy[i]))) + geom_line(stat="density") + geom_rug()
ggsave(paste0("/Users/basiri/LHA/figs/monocle/density/190211_MLB003MLB004_control.highfat.clust11Subs.tree.monocle.",colorBy[i],".rugs.pdf"), dpi = 300, limitsize = F)}

# load libraries

library(Seurat)
library(Matrix)
library(dplyr)
library(reshape2)
library(monocle)

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

# read objects

allCells <- readRDS("/Users/basiri/LHA/objects/tree/180117_MLB003MLB004_control.highfat.tree.rds")
commGenes <- readRDS("/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.genes.rds")

groupMarkers <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.master.pVal.rds')
pct <- readRDS("/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pct/180117_MLB003MLB004_control.highfat.tree.pct.rds")
asinhFC <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.clustAvg.master.asinhFC.rds') 

asinhFC <- as.matrix(asinhFC[,c(1:8,10:15)])
pct <- as.matrix(pct[,c(1:8,10:15)])
groupMarkers <- as.matrix(groupMarkers[,c(1:8,10:15)])

clust11.monocle <- readRDS('/Users/basiri/LHA/objects/monocle/190211_MLB003MLB004_control.highfat.clust11Subs.tree.monocle.rds')

# find pseudotimeMarkers between top and bottom quintile of highfat cells ordered by pseudotime -- Supplementary Data S2

clust11.monocle.pData <- pData(clust11.monocle)
clust11.monocle.pData <- clust11.monocle.pData[order(clust11.monocle.pData$Pseudotime),]
clust11 <- SubsetData(allCells, rownames(clust11.monocle.pData))
clust11@data.info <- clust11.monocle.pData

clust11.highfat <- SubsetData(clust11, rownames(subset(clust11@data.info, clust11@data.info$group == "highfat")), scale = T, center = T)

clust11.monocle.early <- rownames(head(clust11.highfat@data.info, nrow(clust11.highfat@data.info)/5))
clust11.monocle.late <- rownames(tail(clust11.highfat@data.info, nrow(clust11.highfat@data.info)/5))

clust11.subs <- SubsetData(clust11.highfat, c(clust11.monocle.early, clust11.monocle.late), scale = T, center = T)

clust11.subs@data.info <- cbind(clust11.subs@data.info, pseudotimeBin = c(rep("early", length(clust11.monocle.early)), rep("late", length(clust11.monocle.late))))
clust11.subs <- SetAllIdent(clust11.subs, "pseudotimeBin")

pseudotimeMarkers <- FindMarkers(clust11.subs, 'late', 'early', genes.use = commGenes, thresh.use = 0, test.use = "bimod", min.pct = 0, min.diff.pct = 0, only.pos = F)
pseudotimeMarkers <- pseudotimeMarkers[,1, drop=F]

write.table(pseudotimeMarkers, "/Users/basiri/LHA/objects/monocle/pseudotimeMarkers/190430_MLB003MLB004_control.highfat.pseudotimeMarkers.txt", sep="\t", col.names=NA)

# load libraries

library(Seurat)
library(Matrix)
library(dplyr)
library(abind)
library(viridis)
library(gplots)
library(made4)
library(matrixStats)
library(reshape2)
library(corrplot)
library(plotrix)
library(spatstat)

viridis100 <- viridis(100, alpha = 1)
magma100 <- magma(100, alpha = 1)
plasma100 <- plasma(100, alpha = 1)
inferno100 <- inferno(100, alpha = 1)

brownwhiteblue <- colorRampPalette(c("#9e968c", "white", "#6e98a3"))(40)

gray100 <- gray.colors(100, start = 0, end = 1, gamma = 0.5, alpha = NULL)

library(RColorBrewer)
brbg <- brewer.pal(11, "BrBG")
brbg <- c(colorRampPalette(c(brbg[1], brbg[6]))(51), colorRampPalette(c(brbg[6], brbg[11]))(51)[-1])

PuOr <- brewer.pal(11, "PuOr")
PuOr <- c(colorRampPalette(c(PuOr[1], PuOr[6]))(51), colorRampPalette(c(PuOr[6], PuOr[11]))(51)[-1])

PiYG <- brewer.pal(11, "PiYG")
PiYG <- c(colorRampPalette(c(PiYG[1], PiYG[6]))(51), colorRampPalette(c(PiYG[6], PiYG[11]))(51)[-1])

library("RColorBrewer")

orangeBlueDark <- colorRampPalette(c("#4D6B7F", "#d8d8d8", "#CC8411"))(80)

RdBu <- brewer.pal(11, "RdBu")
RdBu <- c(colorRampPalette(c(RdBu[1], RdBu[6]))(51), colorRampPalette(c(RdBu[6], RdBu[11]))(51)[-1])


IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
}

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

# read objects

groupMarkers <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.master.pVal.rds')
pct <- readRDS("/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pct/180117_MLB003MLB004_control.highfat.tree.pct.rds")
asinhFC <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.clustAvg.master.asinhFC.rds') 

commGenes <- readRDS("/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.genes.rds")

asinhFC <- as.matrix(asinhFC[,c(1:8,10:15)])
pct <- as.matrix(pct[,c(1:8,10:15)])
groupMarkers <- as.matrix(groupMarkers[,c(1:8,10:15)])

# subset genes by pct.min, pVal.max, delta_asinhFC

pct.min = 0
pVal.max = 10e-4
delta_asinhFC = 0

pct.list <- list()
for (i in seq_along(colnames(pct))){pct.list[[i]] <- names(subset(pct[,i], pct[,i] >= pct.min))}

asinhFC.list <- list()
for (i in seq_along(colnames(asinhFC))){asinhFC.list[[i]] <- names(subset(asinhFC[,i], abs(asinhFC[,i]) >= delta_asinhFC))}

pVal.list <- list()
for (i in seq_along(colnames(groupMarkers))){pVal.list[[i]] <- names(subset(groupMarkers[,i], groupMarkers[,i] <= pVal.max))}

intersect.list <- list()
for (i in seq_along(pVal.list)){intersect.list[[i]] <- Reduce(intersect, c(pct.list[i], asinhFC.list[i], pVal.list[i]))}

genes <- unique(unlist(intersect.list))
asinhFC.subs <- asinhFC[genes,]

cat("number of genes in subset:", length(genes))

genes <- do.call(cbind, intersect.list)

for(i in 1:ncol(genes)){genes[,i][duplicated(genes[,i])] <- NA}
genes <- as.data.frame(genes)
colnames(genes) <- c("c1","c2","c3","c4","c5","c6","c7","c8","c10","c11","c12","c13","c14", "c15")

write.table(genes, "/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/GOterms/ENRICHr/genes/190209_MLB003MLB004_control.highfat.ENRICHr.genes.txt", sep="\t", row.names = F) # replace NA values in excel with blank

for (i in seq_along(intersect.list)){cat(paste0(colnames(genes[i]), " : ", length(intersect.list[[i]]),"\n"))}

# read ENRICHr files into list

ENRICHr <- list.files("/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/GOterms/ENRICHr/out/0pct_10e-4", pattern="*.txt", full.names=T, recursive=T)
ENRICHr.list <- lapply(ENRICHr, read.table,sep="\t", quote = "", header = T, stringsAsFactors=F)
ENRICHr.summary <- lapply(ENRICHr.list, summary)

c.names <- substr(ENRICHr, 156, 157)
ENRICHr.list <- lapply(seq_along(ENRICHr.list), function(i) cbind(ENRICHr.list[[i]][,c(1,4)], cluster = c.names[i])) 
names(ENRICHr.summary) <- substr(ENRICHr, 156, nchar(ENRICHr)-10)
names(ENRICHr.list) <- substr(ENRICHr, 156, nchar(ENRICHr)-10)
                   
GO_Biological_Process_2018 <- do.call(rbind, ENRICHr.list[grepl("GO_Biological_Process_2018", names(ENRICHr.list))])
GO_Cellular_Component_2018 <- do.call(rbind, ENRICHr.list[grepl("GO_Cellular_Component_2018", names(ENRICHr.list))])
GO_Molecular_Function_2018 <- do.call(rbind, ENRICHr.list[grepl("GO_Molecular_Function_2018", names(ENRICHr.list))])
KEGG_2016 <- do.call(rbind, ENRICHr.list[grepl("KEGG_2016", names(ENRICHr.list))])
Reactome_2016 <- do.call(rbind, ENRICHr.list[grepl("Reactome_2016", names(ENRICHr.list))])
Panther_2016 <- do.call(rbind, ENRICHr.list[grepl("Panther_2016", names(ENRICHr.list))])

rownames(GO_Biological_Process_2018) <- NULL
rownames(GO_Cellular_Component_2018) <- NULL
rownames(GO_Molecular_Function_2018) <- NULL
rownames(KEGG_2016) <- NULL
rownames(Reactome_2016) <- NULL
rownames(Panther_2016) <- NULL
                       
clusters <- c("07","08","10","11")

# GO_Biological_Process_2018

database = GO_Biological_Process_2018
db.terms <- c("0048143", "0030497", "1902950", "0006665", "0046467", "0042761", "0006006", "0032869", "0005513", "1901021", "0000165", "0036465", "0099504", "0048268", "0006811")

database.subset.list <-list()
for (i in seq_along(db.terms)){
database.subset.list[[i]] <- database[grepl(db.terms[i],database$Term, ignore.case = TRUE),]}
    
missing.list <- list()
for (i in seq_along(db.terms)){
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) == 0) 
        missing.list[[i]] <- NULL
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) > 0) 
        missing.list[[i]] <- cbind(Term = unique(database.subset.list[[i]]$Term), Adjusted.P.value = 1, cluster = setdiff(clusters, database.subset.list[[i]]$cluster))}


missing.list <- do.call(rbind, missing.list)
database.subset.list <- do.call(rbind, database.subset.list)
database.subset.list <- rbind(database.subset.list, missing.list)
database.subset.list$Adjusted.P.value <- as.numeric(database.subset.list$Adjusted.P.value)

database.subset.list <- distinct(database.subset.list, Term, Adjusted.P.value, cluster, .keep_all = TRUE)
db.subs <- cbind(dcast(melt(database.subset.list), Term~cluster), database = rep("GO_Biological_Process_2018"))
rownames(db.subs) <- db.subs[,1]
GO_BP.subs <- db.subs[-1]

# GO_Cellular_Component_2018

database = GO_Cellular_Component_2018
db.terms <- c("0030425", "0030424")

database.subset.list <-list()
for (i in seq_along(db.terms)){
database.subset.list[[i]] <- database[grepl(db.terms[i],database$Term, ignore.case = TRUE),]}
    
missing.list <- list()
for (i in seq_along(db.terms)){
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) == 0) 
        missing.list[[i]] <- NULL
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) > 0) 
        missing.list[[i]] <- cbind(Term = unique(database.subset.list[[i]]$Term), Adjusted.P.value = 1, cluster = setdiff(clusters, database.subset.list[[i]]$cluster))}


missing.list <- do.call(rbind, missing.list)
database.subset.list <- do.call(rbind, database.subset.list)
database.subset.list <- rbind(database.subset.list, missing.list)
database.subset.list$Adjusted.P.value <- as.numeric(database.subset.list$Adjusted.P.value)

database.subset.list <- distinct(database.subset.list, Term, Adjusted.P.value, cluster, .keep_all = TRUE)
db.subs <- cbind(dcast(melt(database.subset.list), Term~cluster), database = rep("GO_Cellular_Component_2018"))
rownames(db.subs) <- db.subs[,1]
GO_CC.subs <- db.subs[-1]

# GO_Molecular_Function_2018

database = GO_Molecular_Function_2018
db.terms <- c("0005246", "0019855", "0019855")

database.subset.list <-list()
for (i in seq_along(db.terms)){
database.subset.list[[i]] <- database[grepl(db.terms[i],database$Term, ignore.case = TRUE),]}
    
missing.list <- list()
for (i in seq_along(db.terms)){
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) == 0) 
        missing.list[[i]] <- NULL
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) > 0) 
        missing.list[[i]] <- cbind(Term = unique(database.subset.list[[i]]$Term), Adjusted.P.value = 1, cluster = setdiff(clusters, database.subset.list[[i]]$cluster))}


missing.list <- do.call(rbind, missing.list)
database.subset.list <- do.call(rbind, database.subset.list)
database.subset.list <- rbind(database.subset.list, missing.list)
database.subset.list$Adjusted.P.value <- as.numeric(database.subset.list$Adjusted.P.value)

database.subset.list <- distinct(database.subset.list, Term, Adjusted.P.value, cluster, .keep_all = TRUE)
db.subs <- cbind(dcast(melt(database.subset.list), Term~cluster), database = rep("GO_Molecular_Function_2018"))
rownames(db.subs) <- db.subs[,1]
GO_MF.subs <- db.subs[-1]

# KEGG_2016

database = KEGG_2016
db.terms <- c("hsa04915", "hsa05032", "hsa04727", "hsa04723", "hsa04720", "hsa04915", "hsa04720")

database.subset.list <-list()
for (i in seq_along(db.terms)){
database.subset.list[[i]] <- database[grepl(db.terms[i],database$Term, ignore.case = TRUE),]}
    
missing.list <- list()
for (i in seq_along(db.terms)){
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) == 0) 
        missing.list[[i]] <- NULL
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) > 0) 
        missing.list[[i]] <- cbind(Term = unique(database.subset.list[[i]]$Term), Adjusted.P.value = 1, cluster = setdiff(clusters, database.subset.list[[i]]$cluster))}


missing.list <- do.call(rbind, missing.list)
database.subset.list <- do.call(rbind, database.subset.list)
database.subset.list <- rbind(database.subset.list, missing.list)
database.subset.list$Adjusted.P.value <- as.numeric(database.subset.list$Adjusted.P.value)

database.subset.list <- distinct(database.subset.list, Term, Adjusted.P.value, cluster, .keep_all = TRUE)
db.subs <- cbind(dcast(melt(database.subset.list), Term~cluster), database = rep("KEGG_2016"))
rownames(db.subs) <- db.subs[,1]
KEGG.subs <- db.subs[-1]

# Panther_2016

database = Panther_2016
db.terms <- c("P05734", "P05731", "P00040", "P06959")

database.subset.list <-list()
for (i in seq_along(db.terms)){
database.subset.list[[i]] <- database[grepl(db.terms[i],database$Term, ignore.case = TRUE),]}
    
missing.list <- list()
for (i in seq_along(db.terms)){
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) == 0) 
        missing.list[[i]] <- NULL
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) > 0) 
        missing.list[[i]] <- cbind(Term = unique(database.subset.list[[i]]$Term), Adjusted.P.value = 1, cluster = setdiff(clusters, database.subset.list[[i]]$cluster))}


missing.list <- do.call(rbind, missing.list)
database.subset.list <- do.call(rbind, database.subset.list)
database.subset.list <- rbind(database.subset.list, missing.list)
database.subset.list$Adjusted.P.value <- as.numeric(database.subset.list$Adjusted.P.value)

database.subset.list <- distinct(database.subset.list, Term, Adjusted.P.value, cluster, .keep_all = TRUE)
db.subs <- cbind(dcast(melt(database.subset.list), Term~cluster), database = rep("Panther_2016"))
rownames(db.subs) <- db.subs[,1]
Panther.subs <- db.subs[-1]

# Reactome_2016

database = Reactome_2016
db.terms <- c("R-HSA-622312", "R-HSA-983712", "R-HSA-112310", "R-HSA-112315", "R-HSA-5653656", "R-HSA-166520", "R-HSA-2586552", "R-HSA-112399", "R-HSA-2428924")

database.subset.list <-list()
for (i in seq_along(db.terms)){
database.subset.list[[i]] <- database[grepl(db.terms[i],database$Term, ignore.case = TRUE),]}
    
missing.list <- list()
for (i in seq_along(db.terms)){
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) == 0) 
        missing.list[[i]] <- NULL
    if(length(setdiff(clusters, database.subset.list[[i]]$cluster)) > 0) 
        missing.list[[i]] <- cbind(Term = unique(database.subset.list[[i]]$Term), Adjusted.P.value = 1, cluster = setdiff(clusters, database.subset.list[[i]]$cluster))}


missing.list <- do.call(rbind, missing.list)
database.subset.list <- do.call(rbind, database.subset.list)
database.subset.list <- rbind(database.subset.list, missing.list)
database.subset.list$Adjusted.P.value <- as.numeric(database.subset.list$Adjusted.P.value)

database.subset.list <- distinct(database.subset.list, Term, Adjusted.P.value, cluster, .keep_all = TRUE)
db.subs <- cbind(dcast(melt(database.subset.list), Term~cluster), database = rep("Reactome_2016"))
rownames(db.subs) <- db.subs[,1]
Reactome.subs <- db.subs[-1]

# aggregate all into final matrix and plot -- Figure 2H and Supplementary Data S2

db.master <- rbind(GO_BP.subs, GO_CC.subs, GO_MF.subs, KEGG.subs, Panther.subs, Reactome.subs)
db.master$names <- rownames(db.master)

clust.colors <- c("#a42f2b", "#e4773a", "#ccaf47", "#97bd53", "#749d60", "#3d573c", "#70bb8e", "#3f877f", "gray", "#3480b0", "#2a2d6b", "#bea9cd", "#b874aa", "#ac5399", "#65315a")
clust.colors <- clust.colors[as.numeric(clusters)]

row.colors <- c(NA)
row.colors.list <- list()
for (i in seq_along(unique(db.master$database))){
row.colors <- rep(c(magma(length(unique(db.master$database)), alpha = 1))[i], table(db.master$database)[i])
        assign(paste0('row.colors.',i), row.colors)
        row.colors.list[[i]] <- row.colors                       
        setNames(row.colors.list, paste0('row.colors.',i))}
row.colors <- unlist(row.colors.list)

db.master <- cbind(db.master, db.colors = as.character(row.colors))

hm <- heatmap.2(as.matrix(-log(db.master[,1:4])), Rowv=T, Colv=F, dendrogram="none", col=(orangeBlueDark), ColSideColors=clust.colors, RowSideColors=as.character(db.master$db.colors), trace="none",  scale="row", key=T, cexRow=0.5, cexCol=0.5)
db.master <- db.master[rev(hm$rowInd),]
db.master[db.master == 1] <- NA
rownames(db.master) <- as.character(seq_along(1:nrow(db.master)))

db.master.reorder <- db.master[rev(c('26', '27', '30','21', '33', '31', '32', '34', '8', '6', '7', '3', '4', '2', '1', '9', '5', '23', '25', '13', '36', '11', '16', '17', '18', '19', '22', '20', '37', '35', '12', '14', '15', '24','10', '28', '29')),]
rownames(db.master.reorder) <- as.character(db.master.reorder$names)

pdf("/Users/basiri/LHA/figs/groupMarkers/allCells/asinhFC/GOterms/ENRICHr/190209_MLB003MLB004_control.highfat.ENRICHr.colorblind.pdf")
heatmap.2(as.matrix(-log(db.master.reorder[,1:4])), Rowv=F, Colv=F, dendrogram="none", col=(PuOr), na.color="white", ColSideColors=clust.colors, RowSideColors=as.character(db.master.reorder$db.colors), trace="none",  scale="row", key=F, cexRow=0.5, cexCol=0.5, lwid=c(0.1,10), lhei=c(0.1,4),margin=c(4,36))
dev.off()

pdf("/Users/basiri/LHA/figs/groupMarkers/allCells/asinhFC/GOterms/ENRICHr/190209_MLB003MLB004_control.highfat.ENRICHr.colorblind.KEY.pdf")
heatmap.2(as.matrix(-log(db.master.reorder[,1:4])), Rowv=F, Colv=F, dendrogram="none", col=(PuOr), na.color="white", ColSideColors=clust.colors, RowSideColors=as.character(db.master.reorder$db.colors), trace="none",  scale="row", key=T, cexRow=0.1, cexCol=0.1)
dev.off()

row.names(db.master.reorder) <- db.master.reorder$names
db.master <- db.master.reorder[,1:5]
names(db.master.reorder) <- c("cluster07.pVal","cluster08.pVal","cluster10.pVal","cluster11.pVal", "database")

db.master.reorder <- cbind(db.master.reorder, type = c(rep("other", 3), rep("synapse", 8), rep("ion", 6), rep("signaling", 9), "other", "synapse", "signaling", rep("lipid",4), "other", "signaling", "synapse", "other"))
# db.master.reorder[,1:4] <- -log(db.master.reorder[,1:4])
db.master.reorder <- db.master.reorder[,c(1:5,8)]

saveRDS(db.master.reorder, "/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/GOterms/ENRICHr/master/190209_MLB003MLB004_control.highfat.ENRICHr.rds")
write.table(db.master.reorder, "/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/GOterms/ENRICHr/master/190209_MLB003MLB004_control.highfat.ENRICHr.txt", sep="\t", col.names=NA)

# Create matrix of all databases and all terms -- Supplementary Data S2

db.list <- list(GO_Biological_Process_2018, GO_Cellular_Component_2018, GO_Molecular_Function_2018, KEGG_2016, Panther_2016, Reactome_2016)
names(db.list) <- c("GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018", "KEGG_2016", "Panther_2016", "Reactome_2016")

db.list.subs <- list()
for (i in seq_along(db.list)){ db.list.subs[[i]] <- cbind(dcast(melt(db.list[[i]]), Term~cluster), database = rep(names(db.list[i])))}

db.list.subs <- as.matrix(do.call(rbind, db.list.subs))

saveRDS(db.list.subs, "/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/GOterms/ENRICHr/master/190425_MLB003MLB004_control.highfat.ENRICHr.ALL.rds")
write.table(db.list.subs, "/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/GOterms/ENRICHr/master/190425_MLB003MLB004_control.highfat.ENRICHr.ALL.txt", sep='\t', col.names=NA)

# load libraries

library(Seurat)
library(Matrix)
library(dplyr)
library(abind)
library(matrixStats)
library(reshape2)
library(made4)
library(viridis)
library(RColorBrewer)
library(monocle)
source(paste0(.libPaths(), "/MarcusFuncs/multiplot.R"))

absMax <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
absMin <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)

PiYG <- brewer.pal(11, "PiYG")
PiYG <- c(colorRampPalette(c(PiYG[1], PiYG[6]))(51), colorRampPalette(c(PiYG[6], PiYG[11]))(51)[-1])
 
RdBu <- brewer.pal(11, "RdBu")
RdBu <- c(colorRampPalette(c(RdBu[1], RdBu[6]))(51), colorRampPalette(c(RdBu[6], RdBu[11]))(51)[-1])
    
PuOr <- brewer.pal(11, "PuOr")
PuOr <- c(colorRampPalette(c(PuOr[1], PuOr[6]))(51), colorRampPalette(c(PuOr[6], PuOr[11]))(51)[-1])

magma100 <- magma(100, alpha = 1)

IRdisplay::display_html("<style> .container { width:95% !important; } </style>")

 firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# read objects

allCells <- readRDS("/Users/basiri/LHA/objects/tree/180117_MLB003MLB004_control.highfat.tree.rds")
commGenes <- readRDS("/Users/basiri/LHA/objects/genes/180112_MLB003MLB004_control.highfat.logMed.combatPar.commGenes.genes.rds")

groupMarkers <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/pVal/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.master.pVal.rds')
pct <- readRDS("/Users/basiri/LHA/objects/clusterMarkers/allCells/master/pct/180117_MLB003MLB004_control.highfat.tree.pct.rds")
asinhFC <- readRDS('/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/master/180117_MLB003MLB004_control.highfat.tree.groupMarkers.bimod.clustAvg.master.asinhFC.rds') 

asinhFC <- as.matrix(asinhFC[,c(1:8,10:15)])
pct <- as.matrix(pct[,c(1:8,10:15)])
groupMarkers <- as.matrix(groupMarkers[,c(1:8,10:15)])

clust11.monocle <- readRDS('/Users/basiri/LHA/objects/monocle/190211_MLB003MLB004_control.highfat.clust11Subs.tree.monocle.rds')
db.master <- readRDS("/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/GOterms/ENRICHr/master/190209_MLB003MLB004_control.highfat.ENRICHr.rds")

# read ENRICHr files into list and extract term.genes

ENRICHr <- list.files("/Users/basiri/LHA/objects/groupMarkers/allCells/asinhFC/GOterms/ENRICHr/out/0pct_10e-4", pattern="*.txt", full.names=T, recursive=T)
ENRICHr.list <- lapply(ENRICHr, read.table,sep="\t", quote = "", header = T, stringsAsFactors=F)
ENRICHr.summary <- lapply(ENRICHr.list, summary)

c.names <- substr(ENRICHr, 156, 157)
ENRICHr.list <- lapply(seq_along(ENRICHr.list), function(i) cbind(ENRICHr.list[[i]], cluster = c.names[i])) 
names(ENRICHr.summary) <- substr(ENRICHr, 156, nchar(ENRICHr)-10)
names(ENRICHr.list) <- substr(ENRICHr, 156, nchar(ENRICHr)-10)
                   
GO_Biological_Process_2018 <- do.call(rbind, ENRICHr.list[grepl("GO_Biological_Process_2018", names(ENRICHr.list))])
GO_Cellular_Component_2018 <- do.call(rbind, ENRICHr.list[grepl("GO_Cellular_Component_2018", names(ENRICHr.list))])
GO_Molecular_Function_2018 <- do.call(rbind, ENRICHr.list[grepl("GO_Molecular_Function_2018", names(ENRICHr.list))])
KEGG_2016 <- do.call(rbind, ENRICHr.list[grepl("KEGG_2016", names(ENRICHr.list))])
Reactome_2016 <- do.call(rbind, ENRICHr.list[grepl("Reactome_2016", names(ENRICHr.list))])
Panther_2016 <- do.call(rbind, ENRICHr.list[grepl("Panther_2016", names(ENRICHr.list))])

ENRICHr.list.bind <- do.call(rbind, ENRICHr.list)
  

# plot pseudotime heatmap for term.genes -- Figure S4B

pct.min = 0.5
pVal.max = 1
delta_asinhFC = 0

pct.list <- list()
for (i in seq_along(colnames(pct))){pct.list[[i]] <- names(subset(pct[,i], pct[,i] >= pct.min))}

asinhFC.list <- list()
for (i in seq_along(colnames(asinhFC))){asinhFC.list[[i]] <- names(subset(asinhFC[,i], abs(asinhFC[,i]) >= delta_asinhFC))}

pVal.list <- list()
for (i in seq_along(colnames(groupMarkers))){pVal.list[[i]] <- names(subset(groupMarkers[,i], groupMarkers[,i] <= pVal.max))}

intersect.list <- list()
for (i in seq_along(pVal.list)){intersect.list[[i]] <- Reduce(intersect, c(pct.list[i], asinhFC.list[i], pVal.list[i]))}

genes <- unique(unlist(intersect.list))

db.master.terms <-row.names(subset(db.master, db.master$type == "signaling" | db.master$type == "ion" | db.master$type == "synapse"))

ENRICHr.list.bind.subs <- ENRICHr.list.bind[ENRICHr.list.bind$Term %in% db.master.terms, ]

term.genes <- strsplit(ENRICHr.list.bind.subs$Genes,';')
term.genes <- tolower(unlist(term.genes))                       
term.genes <- intersect(unique(firstup(term.genes)), commGenes)
term.genes <- intersect(term.genes, genes)

len = nrow(pData(clust11.monocle[term.genes,]))

newdata <- data.frame(Pseudotime = seq(min(pData(clust11.monocle[term.genes,])$Pseudotime), max(pData(clust11.monocle[term.genes,])$Pseudotime),length.out = len)) 
  
m <- genSmoothCurves(clust11.monocle[term.genes,], cores=1, trend_formula = '~sm.ns(Pseudotime, df=3)', relative_expr = T, new_data = newdata)

#remove genes with no expression in any condition
m=m[!apply(m,1,sum)==0,]
m = vstExprs(clust11.monocle[term.genes,], expr_matrix=m)

# Row-center the data
m=m[!apply(m,1,sd)==0,]
m=Matrix::t(scale(Matrix::t(m),center=TRUE))
m=m[is.na(row.names(m)) == FALSE,]
m[is.nan(m)] = 0
scale_max=3
scale_min=-3

m[m>scale_max] = scale_max
m[m<scale_min] = scale_min
heatmap_matrix <- m

heatmap_matrix.out <- heatmap.2(heatmap_matrix, Rowv=T, Colv=F, dendrogram="none", col=magma100, trace="none",  scale="none", key=T, cexRow=0.5, labCol = F, distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="average"))
row.order <- rownames(heatmap_matrix)[heatmap_matrix.out$rowInd]

pdf("/Users/basiri/LHA/figs/monocle/heatmap/190212_MLB003MLB004_control.highfat.clust11Subs.tree.monocle.heatmap.colorblind.pdf")                            
heatmap.2(heatmap_matrix, Rowv=T, Colv=F, dendrogram="none", col=magma100, trace="none",  scale="none", key=T, cexRow=0.5, labCol = F, distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="average"))
dev.off()

# plot dots for type
          
genes.terms <- ENRICHr.list.bind.subs[c("Term","Genes")]
genes.terms <- cbind(genes.terms, type = (db.master[match(genes.terms$Term, rownames(db.master)),][6]))
rownames(genes.terms) <- NULL

genes.terms$Genes <- strsplit(genes.terms$Genes,';')
genes.terms <- genes.terms[c("Genes", "type")]
genes.terms <- cbind(genes = firstup(tolower(unlist(with(genes.terms, rep(genes.terms$Genes, vapply(genes.terms$type, length, 1L)))))), type = as.character(with(genes.terms, rep(genes.terms$type, vapply(genes.terms$Genes, length, 1L)))))

genes.terms <- genes.terms[!duplicated(genes.terms), ]
genes.terms.counts <- as.matrix(cbind(genes = unique(genes.terms[,"genes"]), synapse = rep(as.numeric(0), length(unique(genes.terms[,"genes"]))), ion = rep(as.numeric(0), length(unique(genes.terms[,"genes"]))), signaling = rep(as.numeric(0), length(unique(genes.terms[,"genes"])))))

genes.terms.synapse <- subset(genes.terms, genes.terms[,"type"] == "synapse")[,"genes"]
genes.terms.ion <- subset(genes.terms, genes.terms[,"type"] == "ion")[,"genes"]
genes.terms.signaling <- subset(genes.terms, genes.terms[,"type"] == "signaling")[,"genes"]

synapse.indx <- with(data.frame(genes.terms.counts), replace(genes.terms.counts[,"synapse"], genes.terms.counts[,"genes"] %in% genes.terms.synapse, as.numeric(1)))
ion.indx <- with(data.frame(genes.terms.counts), replace(genes.terms.counts[,"ion"], genes.terms.counts[,"genes"] %in% genes.terms.ion, as.numeric(1)))
signaling.indx <- with(data.frame(genes.terms.counts), replace(genes.terms.counts[,"signaling"], genes.terms.counts[,"genes"] %in% genes.terms.signaling, as.numeric(1)))

genes.terms.counts[,"synapse"] <- as.numeric(synapse.indx)
genes.terms.counts[,"ion"] <- as.numeric(ion.indx)
genes.terms.counts[,"signaling"] <- as.numeric(signaling.indx)

genes.terms.counts <- as.data.frame(genes.terms.counts)
          
genes.terms.counts[,2] <- as.numeric(genes.terms.counts[,2])-1
genes.terms.counts[,3] <- as.numeric(genes.terms.counts[,3])-1
genes.terms.counts[,4] <- as.numeric(genes.terms.counts[,4])-1

genes.terms.counts.melt <- melt(genes.terms.counts)

ggplot(genes.terms.counts.melt, aes(y = variable, x = factor(genes, level = rev(row.order)))) + geom_point(aes(colour = value))+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
ggsave("/Users/basiri/LHA/figs/monocle/heatmap/190212_MLB003MLB004_control.highfat.clust11Subs.tree.monocle.heatmap.types.dots.pdf", height = 1.14, width = 10)     

# monocle featurePlot -- Figure S4B

pct.min = 0.50
pVal.max = 10e-4
delta_asinhFC = 0.05

pct.list <- list()
for (i in seq_along(colnames(pct))){pct.list[[i]] <- names(subset(pct[,i], pct[,i] >= pct.min))}

asinhFC.list <- list()
for (i in seq_along(colnames(asinhFC))){asinhFC.list[[i]] <- names(subset(asinhFC[,i], abs(asinhFC[,i]) >= delta_asinhFC))}

pVal.list <- list()
for (i in seq_along(colnames(groupMarkers))){pVal.list[[i]] <- names(subset(groupMarkers[,i], groupMarkers[,i] <= pVal.max))}

intersect.list <- list()
for (i in seq_along(pVal.list)){intersect.list[[i]] <- Reduce(intersect, c(pct.list[i], asinhFC.list[i], pVal.list[i]))}

genes <- unique(unlist(intersect.list))
genes <- unique(intersect(genes, rownames(heatmap_matrix)))

clust11.monocle.dims <- plot_cell_trajectory(clust11.monocle)
clust11.monocle.dims <- clust11.monocle.dims$data
rownames(clust11.monocle.dims) <- clust11.monocle.dims[,1]
clust11.monocle.dims <- clust11.monocle.dims[order(clust11.monocle.dims$Pseudotime),]
row.names(clust11.monocle.dims) <- as.character(seq_along(1:nrow(clust11.monocle.dims)))

clust11.monocle.data.dims.list <- list()

for (i in seq_along(genes)){
clust11.monocle.data.dims.list[[i]] <- merge(clust11.monocle.dims,heatmap_matrix[genes[i],],by="row.names",all.x=TRUE)

rownames(clust11.monocle.data.dims.list[[i]]) <- clust11.monocle.data.dims.list[[i]][,1]
clust11.monocle.data.dims.list[[i]] <- clust11.monocle.data.dims.list[[i]][-1]
rownames(clust11.monocle.data.dims.list[[i]]) <- NULL
names(clust11.monocle.data.dims.list[[i]])[ncol(clust11.monocle.data.dims.list[[i]])] = "feature"}

for (i in seq_along(clust11.monocle.data.dims.list)){
svg(paste0("/Users/basiri/LHA/figs/monocle/featurePlot/190212_MLB003MLB004_control.highfat.monocle.featureplot.",genes[i],".svg"))
print(ggplot(as.data.frame(clust11.monocle.data.dims.list[[i]]), aes(x=data_dim_1,y=data_dim_2)) + geom_point(aes(fill=feature), shape = 21, size = 3, stroke = 0) + scale_fill_gradientn(limits = c(min(heatmap_matrix.out$breaks),max(heatmap_matrix.out$breaks)), colours=magma100) + ggtitle(genes[i]))
dev.off()}

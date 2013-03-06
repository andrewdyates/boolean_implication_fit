### Testing suite on synthetic data
library("RColorBrewer")
library("gplots")
library("energy")
source("compute_mine.R")

all.pairs.dcor <- function(M) {
  n <- nrow(M)
  outer(1:n,1:n, FUN = Vectorize(function(i,j) dcor(M[i,],M[j,])))
}

add.names <- function(M, D) {
  rownames(M) <- rownames(D)
  colnames(M) <- rownames(D)
  M
}
make.hist <- function(M, name, xlim=c(0,1)) {
  pdf(paste0("hist_",name,"_faked.pdf"), width=8, height=8)
  hist(M[lower.tri(M, diag=FALSE)], main=name, xlim=xlim)
  dev.off()
}

# Heatmap Colors
breaks <- seq(-1,1,0.05) # 21
hmcol <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(40)
hmcol[20:21] <- "#FFFFFF"
# --------------------------------------------------


# 1. Read Fake Data
# ==============================
D <- as.matrix(read.table("example/fake_data.tab", sep="\t", header=TRUE, row.names=1))
# Classes of genes
RowSideColors <- c(rep("#FF7F00", 4), rep("#984EA3", 12))

# There is some dependence with the enumeration
x<-D[5,]
y<-D[1,]
reg1 <- lm(y~x)
pdf("RBFOX1_enumeration_scatterplot.pdf", 8, 8)
plot(x,y, xlab="Enumeration", ylab="RBFOX1", main="RBFOX1 PCC:-0.18 (pv: 3e-4), dCOR:0.18")
par(cex=.8)
abline(reg1, col="red", lwd=2, lty=1)
dev.off()
y<-D[3,]
reg1 <- lm(y~x)
pdf("A2M_enumeration_scatterplot.pdf", 8, 8)
plot(x,y, xlab="Enumeration", ylab="A2M", main="RBFOX1 PCC:+0.14 (pv: 4e-3), dCOR:0.14")
par(cex=.8)
abline(reg1, col="red", lwd=2, lty=1)
dev.off()



# 2. Get Dependency Matrices
# ==============================
# MIC Statistical significance n=380
# from http://exploredata.net
# 0.27211 p-value: 0.000001295
# Convert result table into matrices
MINE.T <- read.table("example/fake_data.csv,allpairs,cv=0.0,B=n^0.6,Results.csv", header=TRUE, sep=",")

MIC <- mat.or.vec(16,16); MIC[,] <- 1
MAS  <- mat.or.vec(16,16)
NONLIN <- mat.or.vec(16,16)
MEV <- mat.or.vec(16,16)
MCN <- mat.or.vec(16,16)
  
MINE.T$x <- match(MINE.T$X.var, rownames(D))
MINE.T$y <- match(MINE.T$Y.var, rownames(D))

for (i in 1:nrow(MINE.T)) {
  r <- MINE.T[i,]
  MIC[r$x,r$y] <- r$MIC..strength.; MIC[r$y,r$x] <- MIC[r$x,r$y]
  MAS[r$x,r$y] <- r$MAS..non.monotonicity.; MAS[r$y,r$x] <- MAS[r$x,r$y]
  MEV[r$x,r$y] <- r$MEV..functionality.; MEV[r$y,r$x] <- MEV[r$x,r$y]
  MCN[r$x,r$y] <- r$MCN..complexity.; MCN[r$y,r$x] <- MCN[r$x,r$y]
}

MIC <- sqrt(MIC)
# Compute all-pairs DCOR
DCOR <- all.pairs.dcor(as.matrix(D))
# Compute all-pairs PCC
PCC <- cor(t(D))
# Calculate Spearman's Rho
SPEAR <- cor(apply(D, 1, rank))


# Add names to all matrices
PCC <- add.names(PCC, D)
SPEAR <- add.names(SPEAR, D)
DCOR <- add.names(DCOR, D)
MIC <- add.names(MIC, D)
MAS <- add.names(MAS, D)
MEV <- add.names(MEV, D)
MCN <- add.names(MCN , D)

make.hist(abs(PCC), "abs(PCC)")
make.hist(PCC, "PCC", xlim=c(-1,1))
make.hist(abs(SPEAR), "abs(Spearman)")
make.hist(SPEAR, "Spearman", xlim=c(-1,1))
make.hist(DCOR, "DCOR")
make.hist(MIC, "MIC")
make.hist(MAS, "MAS")
make.hist(MEV, "MEV")
make.hist(MCN, "MCN", xlim=c(0,max(MCN)))
make.hist(DCOR-abs(PCC), "DCOR-abs(PCC)", xlim=c(-0.1,1))
make.hist(SPEAR-PCC, "SPEAR-PCC", xlim=c(-0.5,0.5))
make.hist(abs(SPEAR)-abs(PCC), "abs(SPEAR)-abs(PCC)")
make.hist(MIC, "sqrt(MIC)-abs(PCC)")
           
pdf("splom_fake_data_original_order.pdf", width=30, height=30)
par(mar = rep(0, 4))
pairs(t(D), pch=1)
dev.off()

# Try to pick out clusters based on different dependency measures
# ----------------------------------------
cols <- RowSideColors
# PCC
h <- hclust(dist(abs(PCC)), method="complete")
h <- as.dendrogram(h)
h <- reorder(h, 1:16)
pdf("pcc_clust_faked.pdf", width=10, height=10)
par(mar = rep(0, 4))
heatmap.2(abs(PCC), Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="abs(PCC)")
heatmap.2(PCC, Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="PCC")
heatmap.2(abs(PCC), Rowv=h, Colv=h, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="row", main="abs(PCC) clustered")
heatmap.2(PCC, Rowv=h, Colv=h, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="row", main="PCC; clustered on abs(PCC)")
dev.off()
# ----------------------------------------

# SPEARMAN
h <- hclust(dist(abs(SPEAR)), method="complete")
h <- as.dendrogram(h)
h <- reorder(h, 1:16)
pdf("spearman_clust_faked.pdf", width=10, height=10)
par(mar = rep(0, 4))
par(mar = rep(0, 4))
heatmap.2(abs(SPEAR), Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="abs(SPEARMAN)")
heatmap.2(SPEAR, Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="SPEARMAN")
heatmap.2(SPEAR-PCC, Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="SPEARMAN - PCC")
heatmap.2(abs(SPEAR)-abs(PCC), Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="abs(SPEARMAN) - abs(PCC)")
heatmap.2(abs(SPEAR), Rowv=h, Colv=h, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="row", main="abs(SPEARMAN) clustered")
heatmap.2(SPEAR, Rowv=h, Colv=h, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="row", main="SPEARMAN; clustered on abs(SPEARMAN)")
dev.off()
# ----------------------------------------

# DCOR
h <- hclust(dist(DCOR), method="complete")
h <- as.dendrogram(h)
h <- reorder(h, 1:16)
pdf("dcor_cluster_faked.pdf", width=10, height=10)
par(mar = rep(0, 4))
heatmap.2(DCOR, Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="DCOR")
heatmap.2(DCOR-abs(PCC), Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="DCOR - abs(PCC)")
heatmap.2(DCOR-abs(PCC), Rowv=FALSE, Colv=FALSE, RowSideColors=cols, col=hmcol[18:40], trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="DCOR - abs(PCC) (high intensity)")
heatmap.2(DCOR-abs(SPEAR), Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="DCOR - abs(SPEAR)")
heatmap.2(DCOR-MIC, Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="DCOR - MIC")
heatmap.2(DCOR, Rowv=h, Colv=h, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="row", main="DCOR Clustered")
dev.off()
# ----------------------------------------

# MIC
h <- hclust(dist(MIC), method="complete")
pdf("mic_cluster_faked.pdf", width=10, height=10)
par(mar = rep(0, 4))
heatmap.2(MIC, Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="sqrt(MIC)")
heatmap.2(MIC-abs(PCC), Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="sqrt(MIC)-abs(PCC)")
heatmap.2(MIC-abs(SPEAR), Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="sqrt(MIC)-abs(SPEAR)")
heatmap.2(MIC, Rowv=as.dendrogram(h), Colv=as.dendrogram(h), RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="row", main="sqrt(MIC)")
dev.off()
# ----------------------------------------

# MAS
pdf("MAS_MEV_heatmap_faked.pdf", width=10, height=10)
par(mar = rep(0, 4))
heatmap.2(MAS, Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="MAS")
heatmap.2(MEV, Rowv=FALSE, Colv=FALSE, RowSideColors=cols, breaks=breaks, col=hmcol, trace="none", cexRow=0.7, cexCol=0.7, symm=TRUE, dendrogram="none", main="MEV")
dev.off()

save(D, DCOR, PCC, SPEAR, MIC, MAS, MEV, file="simulated.data.RData")

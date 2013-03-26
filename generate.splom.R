load("example/simulated.data.RData")
load("example/cls.faked.RData")
# Obviously, this path would have to be corrected in actual usage.
source("/Users/z/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")
source("bool.R")

dim(DCOR)
dim(CLS)
dim(MAS)
dim(PCC)

CLS.text <- CLS
## CLS[CLS=="NA"]<-0
## CLS[CLS=="HIH"]<-1
## CLS[CLS=="PC"]<-2
## CLS[CLS=="LIL"]<-3
## CLS[CLS=="UNL"]<-4
## CLS[CLS=="HIL"]<-5
## CLS[CLS=="NC"]<-6
## CLS[CLS=="LIH"]<-7
## CLS <- matrix(as.numeric(CLS),nrow(CLS),ncol(CLS))
CLS <- apply(CLS, c(1,2), cls.to.enum)

pdf("glyph_splom_preserve_order_faked.pdf", width=10, height=10)
#par(mar = rep(0, 4))
splom(CLS, DCOR, asGlyphs=TRUE, lwd=3, grid.col="white", reorder=FALSE, high.sat=FALSE, MIN=0.18)
splom(CLS, DCOR, asGlyphs=TRUE, lwd=3, grid.col="white", reorder=FALSE, high.sat=FALSE, MAX=1, MIN=0.18)
splom(CLS, DCOR, asGlyphs=T, lwd=3, grid.col="white", MIN=0.18, MAX=0.8)
SPLOM.R <- splom(CLS, DCOR, asGlyphs=T, lwd=3, grid.col="white", MIN=0.18, MAX=1)
dev.off()
pdf("synthetic.dendrogram.pdf", width=15, height=8)
plot(SPLOM.R$Rhclust)
dev.off()

pdf("synthetic.data.splom.pdf", width=20, height=20)
pairs(t(D), main="original order")
x <- order.dendrogram(SPLOM.R$Rhclust)
pairs(t(D[x,]), main="GSPLOM order")
dev.off()

# Order D in GSPLOM order
D.original <- D
D <- D[x,]
# add gaps in variable names to handle rows and columns
varnames <- sapply(1:32, function(i) if(i%%2==1) rownames(D)[(i/2)+1] else "")
MAS.CLS <- expand.cls(CLS[x,x], MAS[x,x], bg=NA)
rownames(MAS.CLS) <- varnames
colnames(MAS.CLS) <- varnames
#Img <- t(MAS.CLS)[,seq(nrow(MAS.CLS),1,-1)]
w<-ncol(MAS.CLS); h<-nrow(MAS.CLS)
#image(1:w, 1:h, Img, col=c("black"), breaks=c(-1,1), axes=FALSE, xlab="", ylab="")

library("gplots")
breaks <- seq(-1,1,0.05) # 21
hmcol <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(40)
hmcol[20:21] <- "#FFFFFF"
pdf("MAS_glyph_splom_overlay_faked.pdf", width=10, height=10)
heatmap.2(MAS.CLS, Rowv=FALSE, Colv=FALSE, breaks=breaks, col=hmcol, trace="none", cexRow=0.6, cexCol=0.6, symm=TRUE, dendrogram="none", main="MAS Overlay on SPLOM")
dev.off()


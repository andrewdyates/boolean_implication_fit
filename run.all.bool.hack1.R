## WARNING: THIS IS SUPER HACKY
## Not optimized: boolean classification
## Sample:
## Rscript run.all.bool.R $HOME/c.elegans/clean_data/GSE2180.clean.RData

library(Biobase)
#library(lumi)
source("step.up.R")
source("bool.R")

print.eval = TRUE
argv <- commandArgs(trailingOnly = TRUE)

## 1. Load matrices (WARNING: SUPER HACKY!)
load(argv[1])
ROW <- exprs(GSE2180)
dim(ROW)
if(length(argv) >= 2) {
  #load(argv[2])
  COL <- M
  outname <- paste0(argv[1],".vs.",argv[2],".cls.RData")
} else {
  COL <- NULL
  outname <- paste0(argv[1],".cls.RData")
}
## 2. Step-fit matrix 1 (rows)
RS.ROW <- list()
for(i in 1:nrow(ROW)) {
  name <- rownames(ROW)[i]
  RS.ROW[[name]] <- fit.upstep(ROW[i,])
}
row.stds <- apply(ROW,1,sd)
b.row <- quantile(row.stds,0.03)*2
print("Row Quantile")
print(b.row)

## 3. Step-fit matrix 2 (cols, if it exists) 
if (!is.null(COL)) {
  RS.COL <- list()
  for(i in 1:nrow(COL)) {
    name <- rownames(COL)[i]
    RS.COL[[name]] <- fit.upstep(COL[i,])
  }
  col.stds <- apply(COL,1,sd)
  b.col <- quantile(col.stds,0.03)*2
  print("Col Quantile")
  print(b.col)
} else {
  COL <- ROW
  RS.COL <- RS.ROW
  b.col <- b.row
}
## 4. For every pair, classify
m <- length(RS.ROW)
n <- length(RS.COL)
CLS <- mat.or.vec(m,n)

for(i in 1:n) {
  print(i)
  for(j in 1:n) {
    y <- ROW[i,]
    y.th <- RS.ROW[[i]]$th
    x <- COL[j,]
    x.th <- RS.COL[[j]]$th
    RR <- cls.pair(x,y,x.th,y.th,b.col,b.row, do.plot=F)
    CLS[i,j] <- cls.to.enum(RR$CLS)
  }
}
print(outname)
save(CLS, RS.ROW, RS.COL, b.row, b.col, file=outname)
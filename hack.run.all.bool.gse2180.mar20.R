## This is a horrible hack that only exists to compute all pairs bools for one dataset once.
## Once I learn how to manage packages in R, this type of script will be irrelevant.
source("bool.R")
source("step.up.R")
load("/nfs/01/osu6683/c.elegans/GSE2180.SCAN.N2.ms.max_0.5.RData")

all.steps <- function(M, do.plot=FALSE) {
  STEPS <- apply(M, 1, fit.upstep)
  if(do.plot) {
    for (i in 1:dim(M)[1]) {
      title=rownames(M)[i]
      plot.stepfit(STEPS[[i]], v=M[i,], add.mean.median=T, main=paste(title, "Stepfit"))
      plot.sse(STEPS[[i]], add.mean.median=T, main=paste(title, "SSE"))
    }
  }
  STEPS
}

all.pairs.cls <- function(M, steps, b) {
  n <- dim(M)[1]
  CLS <- mat.or.vec(n,n)
  for(i in 1:n) { # row
    for(j in 1:n) { # col
      y <- M[i,]
      y.th <- steps[[i]]$th
      x <- M[j,]
      x.th <- steps[[j]]$th
      x.title=rownames(M)[j]
      y.title=rownames(M)[i]
      RR <- cls.pair(x,y,x.th,y.th, b.x=b, b.y=b, do.plot=F, xlab=x.title, ylab=y.title)
      CLS[i,j] <- cls.to.enum(RR$CLS)
    }
  }
  rownames(CLS) <- rownames(M)
  colnames(CLS) <- rownames(M)
  CLS
}

b <- 0.1726519 # manually
steps <- all.steps(M)
CLS <- all.pairs.cls(M, steps, b=b)
save(steps, CLS, b, file="/nfs/01/osu6683/c.elegans/mar20.gse2180.steps.cls.RData")
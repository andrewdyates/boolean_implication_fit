# requires bool.R and step.up.R

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

all.pairs.cls <- function(M, steps, b, do.plot=FALSE, r.th=2/3) {
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
      RR <- cls.pair(x,y,x.th,y.th, b.x=b, b.y=b, do.plot=do.plot, xlab=x.title, ylab=y.title, r.th=r.th)
      CLS[i,j] <- cls.to.enum(RR$CLS)
    }
  }
  rownames(CLS) <- rownames(M)
  colnames(CLS) <- rownames(M)
  CLS
}

# compute all pairs given a particular row
single.pairs.cls <- function(M, steps, b, i, do.plot=FALSE) {
  n <- dim(M)[1]
  CLS <- mat.or.vec(n,1)
  names(CLS) <- rownames(M)
  for(j in 1:n) { # col
    y <- M[i,]
    y.th <- steps[[i]]$th
    x <- M[j,]
    x.th <- steps[[j]]$th
    x.title=rownames(M)[j]
    y.title=rownames(M)[i]
    RR <- cls.pair(x,y,x.th,y.th, b.x=b, b.y=b, do.plot=do.plot, xlab=x.title, ylab=y.title)
    CLS[j] <- cls.to.enum(RR$CLS)
  }
  CLS
}


#STEPS <- all.steps(MyData, do.plot=F)
#all.pairs.cls <- (MyData, STEPS, b, do.plot=F)

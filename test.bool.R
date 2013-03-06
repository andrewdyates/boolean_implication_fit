load("simulated.data.RData")
source("step.up.R")
source("bool.R")


RS = list()
for(i in 1:nrow(D)) {
  name <- rownames(D)[i]
  R <- fit.upstep(D[i,])
  pdf(paste0("sse.",name,".pdf"), width=8, height=8)
  plot.sse(R, add.mean.median=TRUE, main=paste("SSE", name))
  dev.off()
  pdf(paste0("stepfit.",name,".pdf"), width=8, height=8)
  plot.stepfit(R, D[i,], add.mean.median=TRUE, main=paste("Stepfit", name))
  dev.off()
  RS[[name]] <- R
}



i<-6; j<-9
x <- D[i,]
y <- D[j,]
x.th <- RS[[i]]$th
y.th <- RS[[j]]$th
R <- cls.pair(x,y,x.th,y.th, do.plot=TRUE)

## Classify all pairs using boolean implication
n <- nrow(D)
b <- 0.3
CLS <- mat.or.vec(n,n)
CLS[,] <- "NA"
for(i in 1:n) {
  for(j in 1:n) {
    y <- D[i,]
    y.th <- RS[[i]]$th
    x <- D[j,]
    x.th <- RS[[j]]$th
    name <- paste0(rownames(D)[j],".vs.",rownames(D)[i])
    png(paste0(name,".scatter.png"), width=500, height=500)
    RR <- cls.pair(x,y,x.th,y.th,b, do.plot=TRUE)
    dev.off()
    CLS[i,j] <- RR$CLS
  }
}

save(RS, CLS, file="example/cls.faked.RData")

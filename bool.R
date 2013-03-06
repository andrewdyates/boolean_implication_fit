b <- 0.5

in.th <- function(v, th, b=0.5) v>=th-b & v<=th+b

cls.pair <- function(x, y, x.th, y.th, b=0.5, do.plot=FALSE) {
discard <- in.th(x,x.th,b) | in.th(y,y.th,b)
pct.discard <- sum(discard)/length(x)*100

x0.y0 <- sum(x<x.th-b & y<y.th-b)
x0.y1 <- sum(x<x.th-b & y>y.th+b)
x1.y0 <- sum(x>x.th+b & y<y.th-b)
x1.y1 <- sum(x>x.th+b & y>y.th+b)
x0 <- x0.y0 + x0.y1
x1 <- x1.y0 + x1.y1
y0 <- x0.y0 + x1.y0
y1 <- x0.y1 + x1.y1

# quad sparse tests
total <- x0.y0 + x0.y1 + x1.y0 + x1.y1

# x0.y0
q00 <- list()
q00$expect  <- x0*y0/total
q00$stat <- (q00$exp-x0.y0)/sqrt(q00$expect)
q00$err  <- (x0.y0/x0 + x0.y0/y0)/2
q00$sparse <- q00$stat > 3 && q00$err < 0.1
# x1.y0
q10 <- list()
q10$expect  <- x1*y0/total
q10$stat <- (q10$exp-x1.y0)/sqrt(q10$expect)
q10$err  <- (x1.y0/x1 + x1.y0/y0)/2
q10$sparse <- q10$stat > 3 && q10$err < 0.1
# x0.y1
q01 <- list()
q01$expect  <- x0*y1/total
q01$stat <- (q01$exp-x0.y1)/sqrt(q01$expect)
q01$err  <- (x0.y1/x0 + x0.y1/y1)/2
q01$sparse <- q01$stat > 3 && q01$err < 0.1
# x1.y1
q11 <- list()
q11$expect  <- x1*y1/total
q11$stat <- (q11$exp-x1.y1)/sqrt(q11$expect)
q11$err  <- (x1.y1/x1 + x1.y1/y1)/2
q11$sparse <- q11$stat > 3 && q11$err < 0.1

R <- list(q00=q00,q01=q01,q10=q10,q11=q11,total=total,all=length(x),pct.discard=pct.discard,
          x0.y0=x0.y0, x0.y1=x0.y1, x1.y0=x1.y0, x1.y1=x1.y1, discard=discard)

# Choose class.
if(pct.discard > (2/3*100)) {
  R$CLS <- "NA"
} else {
  R$CLS <- "UNL"
  if(q00$sparse && !q01$sparse && !q10$sparse && !q11$sparse)
    R$CLS <- "LIH"
  if(!q00$sparse && q01$sparse && !q10$sparse && !q11$sparse)
    R$CLS <- "LIL"
  if(!q00$sparse && !q01$sparse && q10$sparse && !q11$sparse)
    R$CLS <- "HIH"
  if(!q00$sparse && !q01$sparse && !q10$sparse && q11$sparse)
    R$CLS <- "HIL"
  if(!q00$sparse && q01$sparse && q10$sparse && !q11$sparse)
    R$CLS <- "PC"
  if(q00$sparse && !q01$sparse && !q10$sparse && q11$sparse)
    R$CLS <- "NC"
}
if(do.plot) {
  col<-rep("black",length(x))
  col[discard]<-"#ffcccc"
  plot(x,y, col=col, xlab=rownames(D)[i], ylab=rownames(D)[j], main=R$CLS)
  abline(v=x.th, col="red", lwd=2)
  abline(v=x.th+b, col="red", lty=3)
  abline(v=x.th-b, col="red", lty=3)
  abline(h=y.th, col="blue", lwd=2)
  abline(h=y.th+b, col="blue",lty=3)
  abline(h=y.th-b, col="blue",lty=3)
}
R
} # END cls.pair()



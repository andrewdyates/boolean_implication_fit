

q2txt <- function(q, name, Z=3, Err=0.1, cnt=1)
  paste(paste0(name, " n/E: ", q$count, "/", formatC(q$expect,digits=2,format="f")),
        paste0("(Z>",Z,",Err<",Err,",n>,",cnt,"): "),
        paste(formatC(q$stat,digits=2), formatC(q$err,digits=2), q$count>cnt, sep=", "),
        paste0("sparse?", q$sparse), sep=" ")

in.th <- function(v, th, b=0.5) v>=th-b & v<=th+b

cls.to.enum <- function(cls) {
  if(cls=="NA") r<-0
  if(cls=="HIH") r<-1
  if(cls=="PC") r<-2
  if(cls=="LIL") r<-3
  if(cls=="UNL") r<-4
  if(cls=="HIL") r<-5
  if(cls=="NC") r<-6
  if(cls=="LIH") r<-7  
  r
}

## MAIN FUNCTION
# x is ROW. y is COLUMN
cls.pair <- function(x, y, x.th, y.th, b.x=0.5, b.y=0.5, do.plot=FALSE, xlab="", ylab="", quad.min=1, Z=3, E=0.1, r.th=2/3) {
discard <- in.th(x,x.th,b.x) | in.th(y,y.th,b.y)
pct.discard <- sum(discard)/length(x)*100

x0.y0 <- sum(x<x.th-b.x & y<y.th-b.y)
x0.y1 <- sum(x<x.th-b.x & y>y.th+b.y)
x1.y0 <- sum(x>x.th+b.x & y<y.th-b.y)
x1.y1 <- sum(x>x.th+b.x & y>y.th+b.y)
x0 <- x0.y0 + x0.y1
x1 <- x1.y0 + x1.y1
y0 <- x0.y0 + x1.y0
y1 <- x0.y1 + x1.y1

# quad sparse tests
total <- x0.y0 + x0.y1 + x1.y0 + x1.y1

# x0.y0
q00 <- list()
q00$count <- x0.y0
q00$expect  <- x0*y0/total
q00$stat <- (q00$exp-x0.y0)/sqrt(q00$expect)
q00$err  <- (x0.y0/x0 + x0.y0/y0)/2
q00$sparse <- q00$stat > Z && q00$err < E || x0.y0 <= quad.min
if (x0*y0==0 && is.na(q00$sparse)) q00$sparse <- T

# x1.y0
q10 <- list()
q10$count <- x1.y0
q10$expect  <- x1*y0/total
q10$stat <- (q10$exp-x1.y0)/sqrt(q10$expect)
q10$err  <- (x1.y0/x1 + x1.y0/y0)/2
q10$sparse <- q10$stat > Z && q10$err < E || x1.y0 <= quad.min
if (x1*y0==0 && is.na(q10$sparse)) q10$sparse <- T

# x0.y1
q01 <- list()
q01$count <- x0.y1
q01$expect  <- x0*y1/total
q01$stat <- (q01$exp-x0.y1)/sqrt(q01$expect)
q01$err  <- (x0.y1/x0 + x0.y1/y1)/2
q01$sparse <- q01$stat > Z && q01$err < E || x0.y1 <= quad.min
if (x0*y1==0 && is.na(q01$sparse)) q01$sparse <- T

# x1.y1
q11 <- list()
q11$count <- x1.y1
q11$expect  <- x1*y1/total
q11$stat <- (q11$exp-x1.y1)/sqrt(q11$expect)
q11$err  <- (x1.y1/x1 + x1.y1/y1)/2
q11$sparse <- q11$stat > Z && q11$err < E || x1.y1 <= quad.min
if (x1*y1==0 && is.na(q11$sparse)) q11$sparse <- T

R <- list(q00=q00,q01=q01,q10=q10,q11=q11,total=total,all=length(x),pct.discard=pct.discard,
          x0.y0=x0.y0, x0.y1=x0.y1, x1.y0=x1.y0, x1.y1=x1.y1, discard=discard)

# Choose class.
if(pct.discard > (r.th*100)) {
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
  if(xlab != "" && ylab != "")
    title <- paste(R$CLS, paste0("(",cls.to.enum(R$CLS),")"), ":", xlab, "vs", ylab)
  else
    title <- paste(R$CLS, paste0("(",cls.to.enum(R$CLS),")"))
  par(mar=c(9.1,4.1,4.1,2.1))
  plot(x,y, col=col, main=title, xlab=xlab, ylab=ylab)

  if (!R$q00$sparse)
    rect(xleft=min(x), ybottom=min(y), xright=x.th-b.x, ytop=y.th-b.y, col=rgb(0,0,0,0.3), density=NA)
  if (!R$q01$sparse)
    rect(xleft=min(x), ybottom=y.th+b.y, xright=x.th-b.x, ytop=max(y), col=rgb(0,0,0,0.3), density=NA)
  if (!R$q10$sparse)
    rect(xleft=x.th+b.x, ybottom=min(y), xright=max(x), ytop=y.th-b.y, col=rgb(0,0,0,0.3), density=NA)
  if (!R$q11$sparse)
    rect(xleft=x.th+b.x, ybottom=y.th+b.y, xright=max(x), ytop=max(y), col=rgb(0,0,0,0.3), density=NA)
  
  abline(v=x.th, col="red", lwd=2)
  abline(v=x.th+b.x, col="red", lty=3)
  abline(v=x.th-b.x, col="red", lty=3)
  abline(h=y.th, col="blue", lwd=2)
  abline(h=y.th+b.y, col="blue",lty=3)
  abline(h=y.th-b.y, col="blue",lty=3)

  mtext(paste(q2txt(q00,"q00"), q2txt(q01,"q01"), q2txt(q10,"q10"), q2txt(q11,"q11"), sep="\n"), side=1,line=8)

}
R
} # END cls.pair()

## Example use of MINE.R
## ------------------------------
#install.packages("rJava")
source("MINE.R")

D <- read.table("example/fake_data.tab", sep="\t", header=TRUE, row.names=1)
write.table(as.matrix(t(D)), sep=",", row.names=FALSE, file="example/fake_data.csv")
MINE("example/fake_data.csv","all.pairs")
#MINE.T <- read.table("fake_data.csv,allpairs,cv=0.0,B=n^0.6,Results.csv", header=TRUE, sep=",")

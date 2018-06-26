
# Setup: load libraries & datasets

library(VGAM)
library(tm)

dir <- DirSource(directory = "E:/Thesis/test/")$filelist

for (i in 1:16) {
  nam <- paste("PG", i, sep = "")
  assign(nam, read.table(dir[i], quote="\"", comment.char=""))
}

# Estimation of pdf of text length L (fig 1)

for (i in 1:16) {
  nam <- paste("fit", i, sep = "")
  PGname <- paste("PG", i, sep = "")
  assign(nam, vglm(V1~., data = noquote(PGname), family = 'diffzeta'))
}

######## testing and cleaing

df <- cbind(nrow(PG1):1, PG1)
plot(df)

qzeta(0.05, shape = 0.74)


######


# Histograms of p-values (fit 2)

# cdf of p-values


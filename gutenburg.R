### Setup

library(VGAM)
library(ggplot2)
library(gridExtra)

# Load 200 datasets from the raw frequency data: 
# into book1, book2, book3... and a list of all.

book.dir <- dir('gutenburgSubset', full.names = T) # load directory
book.list <- list()

for (i in 1:200) { 
  nam <- paste('book', i, sep = '')
  sortFreq <- data.frame('frequency' = 
                           sort(read.table(book.dir[i])$V1, decreasing = T),
                         'rank' = 
                           1:nrow(read.table(book.dir[i])))
  assign(nam, sortFreq)
  book.list[[i]] <- sortFreq
  names(book.list)[i] <- as.character(i)
}


### Quick log-log plot of 12 books with lm fits

group.plot <- function(x) {
  
  qplot(x$rank, x$frequency, log = 'xy', geom = 'line',
        main = expression(which(book.list[1:12] %in% x))) + # to fix
    geom_smooth(method = lm) +
    theme(axis.title = element_blank())
  
}

book.plots <- lapply(book.list[1:12], group.plot)
do.call('grid.arrange', c(book.plots, ncol = 3, nrow = 4,
                          left = 'Frequency', bottom = 'Rank'))



### Quick fit of 1 dataset: book1, with zeta family function
# (Q:) Wanted zetaff(lshape = loglog) for shape > 1 but errors

fit.zeta <- vglm(rank ~ 1, zetaff, data = book1, w = frequency,
                 trace = TRUE, criterion = "coef", maxit = 30)
summary(fit.zeta)



### Quick evaluations and insights

Coef(fit.zeta)  # coeficient of shape parameter = 0.2395966
coef(summary(fit.zeta)) # estimate, stdE, z value & p-value
head(resid(fit.zeta)) # residuals
vcov(fit.zeta) # variance-covariance matrix
df.residual(fit.zeta) # 1371 (number of observation nobs(fit.zeta) - 1)
logLik(fit.zeta) # log likelihood of the fit
npred(fit.zeta) # 1 parameter distribution, hence 1 predictor

has.intercept(fit.zeta) # There is intercept (formula = y ~ 1)
is.parallel(fit.zeta) # Hat matrix = Identity matrix

head(hatvalues(fit.zeta) * sum(book1$frequency), n = 10) # Would be identical
head(book1$frequency, n = 10)                            # to this.
head(weights(fit.zeta), n = 10) # Where frequency is the weights of the case.

head(QR.Q(fit.zeta)) # The Q matrix and
QR.R(fit.zeta)       # the R matrix of
head(model.matrix(fit.zeta)) # the model matrix used to solve
# the linear least square problems in IRLS.



### Model comparisons

# diffzeta (the second variant) with loglog link function for shape > 1
fit.diffzeta <- vglm(rank ~ 1, diffzeta(lshape = 'loglog'),
                     data = book1, w = frequency,
                     trace = TRUE, criterion = 'coef', maxit = 30)

# pareto family: paretoff(), paretoII(), paretoIII() with default loge link.
fit.pareto <- vglm(rank ~ 1, paretoff, data = book1, w = frequency,
                   trace = TRUE, criterion = 'coef')
fit.paretoII <- vglm(rank ~ 1, paretoII, data = book1, w = frequency,
                     trace = TRUE, criterion = 'coef')
fit.paretoIII <- vglm(rank ~ 1, paretoIII, data = book1, w = frequency,
                      trace = TRUE, criterion = 'coef')

AIC(fit.zeta) 
AIC(fit.diffzeta) # somehow not working well
AIC(fit.pareto)
AIC(fit.paretoII)
AIC(fit.paretoIII) # best

AICc(fit.zeta)
AICc(fit.diffzeta)
AICc(fit.pareto)
AICc(fit.paretoII)
AICc(fit.paretoIII) # best

BIC(fit.zeta)
BIC(fit.diffzeta) 
BIC(fit.pareto)
BIC(fit.paretoII)
BIC(fit.paretoIII) # best

TIC(fit.zeta)
TIC(fit.diffzeta) 
TIC(fit.pareto)
TIC(fit.paretoII)
TIC(fit.paretoIII) # best

lrtest(fit.zeta, fit.pareto) # seems like pareto model fits better



### Analysis in a batch of 200 books (out of 31075 books for now)
### With 4 models: Zeta, Pareto, ParetoII, ParetoIII

# 1. Coeficients

vglm.list <- function(x, y) vglm(rank ~ 1, family = y, data = x,
                                 w = frequency, criterion = 'coef')

## Zeta
coefs.zeta <- lapply(book.list,
                     function(x, y) Coef(vglm.list(x, y = zetaff)))
coefs.zeta <- as.data.frame(unlist(coefs.zeta, use.names = F))
colnames(coefs.zeta) <- 'coeficient'
head(coefs.zeta) # quick check

## Pareto
coefs.pareto <- lapply(book.list,
                       function(x, y) Coef(vglm.list(x, y = paretoff)))
coefs.pareto <- as.data.frame(unlist(coefs.pareto, use.names = F))
colnames(coefs.pareto) <- 'coeficient'
head(coefs.pareto) # quick check

## ParetoII (Problems: doesn't converge here... 2 maximums?)
coefs.paretoII <- lapply(book.list,
                         function(x, y) Coef(vglm.list(x, y = paretoII)))
coefs.paretoII <- as.data.frame(unlist(coefs.paretoII, use.names = F))
colnames(coefs.paretoII) <- 'coeficient'
head(coefs.paretoII) # quick check

## ParetoIII
coefs.paretoIII <- lapply(book.list,
                          function(x, y) Coef(vglm.list(x, y = paretoIII)))
coefs.paretoIII <- as.data.frame(unlist(coefs.paretoIII, use.names = F))
colnames(coefs.paretoIII) <- 'coeficient'
head(coefs.paretoIII) # quick check


# Density plot of 200 coeficients

cutoff <- function(coef) { # maximum of density curve
  peak <- which.max(density(coef$coeficient)$y)
  ans <- density(coef$coeficient)$x[peak]
  ans
}

## zeta (cutoff = 0.210116)

coef.plot <- ggplot(coefs.zeta, aes(x = coeficient)) +
  geom_density(colour = 'black',
               fill = 'gray', size = 1, alpha = 0.5) +
  geom_vline(xintercept = cutoff(coefs.zeta), size = 1,
             color = 'blue', linetype = 'dashed') +
  xlim(c(0, 0.3))

coef.plot

## pareto (cutoff = 0.2364423)

coef.plot

## paretoII (cutoff = 0.5580624) problems
## paretoIII (cutoff = 1.373475) problems
  
  
  

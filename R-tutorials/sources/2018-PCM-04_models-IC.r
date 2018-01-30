
aic <- function(x, ...) {
    dev <- -2 * logLik(x)
    K = attr(logLik(x), 'df')
    n = attr(logLik(x), 'nall')
    aic <- dev + 2*K
    aic.c <- aic + (2 * K * (K + 1)) / (n - K - 1)
    out <- list(aic = aic, aic.c = aic.c, n = n, K = K, lnL = logLik(x), param = as.numeric(x$modelStruct))
    class(out) <- 'aic'
    return(out)
}

aic.w <- function(x, ...) {
    aic.lnL <- exp(-0.5 * (x - min(x)))
    aic.w <- aic.lnL / sum(aic.lnL)
    return(aic.w)
}

aic.w.modelSet <- function(x, which.use = 'aic.c', ...) {
    aic.set <- lapply(x, aic)
    aic.set <- sapply(aic.set, function(x) x[[which.use]])
    out <- aic.w(aic.set, ...)
    return(out)
}

### TRYING IT OUT

library(phytools)
library(geiger)
library(nlme)
                      
tr <- pbtree(n = 100)
dat <- as.data.frame(sim.char(tr, matrix(c(1,0.8,
                                           0.8,1),2,2, byrow = T))[, , 1])
names(dat) <- c('y','x')
                      
models <- list(y.x.brown = gls(y ~ x, dat, correlation = corBrownian(1, tr)),
               y.x.pagel = gls(y ~ x, dat, correlation = corPagel(1, tr)), 
               y.brown = gls(y ~ 1, dat, correlation = corBrownian(1, tr)),
               y.x.star = gls(y ~ x, dat, correlation = corPagel(0, tr, fixed = TRUE)))

sapply(models, aic)    
round(aic.w.modelSet(models), 4)

K = c(2, 3, 5, 6, 17)
lnL = c(17.33, 15.69, 24.82, 26.69, 44.17)
aic.boett <- -2*lnL + 2*K
aicc.boett <- -2*lnL + 2*(K*(K + 1)) / (13 - K - 1) 
names(aic.boett) <- names(aicc.boett) <- c('BM', 'OU.1', 'OU.3', 'OU.4', 'OU.15')
message('AIC weights corresponding to Boettiger et al. 2012, Table 1')
cbind(aic.w = round(aic.w(aic.boett), 8),
      aic_c.w = round(aic.w(aicc.boett), 8)
     )


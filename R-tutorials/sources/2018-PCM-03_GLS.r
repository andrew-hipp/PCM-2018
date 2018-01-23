
library(phytools)
library(geiger)

a=tr <- pbtree(n = 100)
dat <- sim.char(tr, matrix(c(1,0,0,1), 2, 2))[,,1]
colnames(dat) <- c('y', 'x')
#layout(matrix(c(1:2, 1, 2)))
par(mfrow = c(1,2))
plot(tr)
plot(dat)

tr.v <- vcv.phylo(tr)
round(tr.v[1:10, 1:10], 2) # shows just the upper lefthand corner of your matrix

round(solve(tr.v), 2)[1:10, 1:10]
round(solve(tr.v) %*% tr.v, 2)[1:10, 1:10]

one <- matrix(1, length(tr$tip.label), 1) # a matrix of 1s, 1 column
tr.w <- t(one) %*% solve(tr.v)
a=plot(tr, show.tip.label = F, x.lim = nodeheight(tr,1) * 3) # sets the plot window width to be thrice the tree depth
phydataplot(tr.w[1,] / max(tr.w[1,]) * nodeheight(tr, 1), tr)

b.gls <- function(tr, dat, y = 'y', x = 'x') {
    #print(dat)
    v.solved <- solve(vcv.phylo(tr))
    one <- matrix(1, length(tr$tip.label))
    X <- cbind(origin = one, x = dat[, x])
    out <- solve(t(X) %*% v.solved %*% X) %*% (t(X) %*% v.solved %*% dat[, y])
    out
}

b.ols <- function(dat, y = 'y', x = 'x') {
    one <- matrix(1, dim(dat)[1])
    X <- cbind(origin = one, x = dat[, x])
    out <- solve(t(X) %*% X) %*% (t(X) %*% dat[, y])
    out
}
message('GLS coefficients (intercept and slope)')
b.gls(tr, dat)
message('OLS coefficients (intercept and slope)')
b.ols(dat)

library(nlme)
out <- gls(y ~ x, as.data.frame(dat), corPagel(1, tr))
summary(out)

gls.r.squared <- function(x) {
# based on Judge et al. 1985, eq. 2.3.16
  e = x$resid
  V <- corMatrix(x$modelStruct$corStruct)
  Y = x$resid + x$fitted
  one <- matrix(1, length(Y), 1)
  a <- as.numeric(solve(t(one) %*% solve(V) %*% one) %*% (t(one) %*% solve(V) %*% Y))
  r.squared= 1-(t(e) %*% solve(V) %*% e) / (t(Y-a) %*% solve(V) %*% (Y-a))
  return(r.squared[1,1])
  }

## example:
gls.r.squared(out) # probably will be pretty low, as we simulated data with no correlation

tr.set <- lapply(rep(1, 500), pbtree, n = 50)
dat.set <- sapply(tr.set, sim.char, par = 1, root = 0)
phylo.mean <- function(tr, dat) {
    #print(dat)
    v.solved = solve(vcv.phylo(tr))
    one = matrix(1, length(tr$tip.label))
    out = solve(t(one) %*% v.solved %*% one) %*% (t(one) %*% v.solved %*% dat)
    out
}
dat.phyloMeans <- sapply(1:length(tr.set), function(x) phylo.mean(tr.set[[x]], dat.set[, x]))
par(mfrow = c(1, 2))
hist(apply(dat.set, 1, mean), 20, main = 'ordinary means', xlim = c(-3, 3))
abline(v = mean(apply(dat.set, 1, mean)), col= 'red', lty = 'dashed')
hist(dat.phyloMeans, 20, main = 'phylogenetic means', xlim = c(-3, 3))
abline(v = mean(dat.phyloMeans), col = 'red', lty = 'dashed')

dat.set <- lapply(tr.set, sim.char, par = matrix(c(1, 0, 0, 1), 2, 2,), root = 0)
dat.set <- lapply(dat.set, function(x) x[,,1])

out <- list(ols = t(sapply(1:length(tr.set), function(i) b.ols(dat.set[[i]], y = 1, x = 2))),
             gls = t(sapply(1:length(tr.set), function(i) b.gls(tr.set[[i]], dat.set[[i]], y = 1, x = 2)))
                 )
colnames(out$ols) <- colnames(out$gls) <- c('intercept', 'slope')
                 
par(mfrow = c(3, 2)) # set up plotting area
                 
for(j in c('intercept', 'slope')) {
    for(i in c('ols', 'gls')) {
        hist(out[[i]][, j], 20, main = paste(i, j), xlim = c(-3, 3))
        abline(v = mean(out[[i]][, j]), col= 'red', lty = 'dashed')
    } # close i
    plot(out$ols[, j], out$gls[, j], xlab = paste('ols', j), ylab = paste('gls', j))
    abline(0, 1, col = 'red', lty = 'dashed')
  } # close j

summary(lm(ols.out[, 2] ~ gls.out[, 2]))
summary(lm(ols.out[, 1] ~ gls.out[, 1]))

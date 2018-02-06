
library(geiger)
tr <- sim.bdtree(n = 20)
dat <- sim.char(tr, 1, 1)

tr2 <- ex.ratesimulator(tr, min = 5)
dat2 <- sim.char(tr2,1,1)

layout(matrix(1:3, 1))
plot(tr, main = 'original tree')
plot(tr2, main = 'rescaled tree')
plot(dat[,,1], dat2[,,1], xlab = 'data, raw tree', ylab = 'data, rescaled tree', main = '')

models <- list(
    origOnOrig = fitContinuous(tr, dat[,,1]),
    origOnRescaled = fitContinuous(tr2, dat[,,1]),
    rescaledOnRescaled = fitContinuous(tr2, dat2[,,1]),
    rescaledOnOrig = fitContinuous(tr,dat2[,,1])
)
models.mat <- cbind(orig = models$origOnOrig$opt,
                    origOnRescaled = models$origOnRescaled$opt,
                    rescaled = models$rescaledOnRescaled$opt,
                    rescaledOnOrig = models$rescaledOnOrig$opt)
print(models.mat)

sims <- list(
    s12 = sim.char(tr2, par = as.numeric(models.mat['sigsq', 'origOnRescaled']), root = as.numeric(models.mat['z0', 'origOnRescaled']), nsim = 20),
    s11 = sim.char(tr, par = as.numeric(models.mat['sigsq', 'orig']), root = as.numeric(models.mat['z0', 'orig']), nsim = 20),
    s22 = sim.char(tr2, par = as.numeric(models.mat['sigsq', 'rescaled']), root = as.numeric(models.mat['z0', 'rescaled']), nsim = 20),
    s21 = sim.char(tr, par = as.numeric(models.mat['sigsq', 'rescaledOnOrig']), root = as.numeric(models.mat['z0', 'rescaledOnOrig']), nsim = 20)
)
sims.out <- sapply(names(sims), function(x) {
    if(substr(x,3,3) == 2) simTree <- tr2
        else simTree <- tr
    out <- apply(sims[[x]], 3, function(y) fitContinuous(simTree, y)$opt$lnL)
    out
        })
sims.out

library(ape)

rootState <- function(tr, X) {
  ## the maximum likelihood estimate of root state based on O'Meara, p. 925
  vcvMat <- vcv(tr)
  one <- matrix(1, length(X), 1) # a matrix of 1s, as tall as the length of X
  B0 <- as.numeric(solve(t(one) %*% solve(vcvMat) %*% one) %*% (t(one) %*% solve(vcvMat) %*% X))
  return(B0)
  }

sigmaSq <- function(tr, X) {
  ## following O'Meara 2006, eq. 2
  vcvMat <- vcv(tr)
  N <- length(X)
  E.X <- matrix(rootState(tr, X), length(X), 1) # expected value of X is just the rootstate of X
  sigmaSq <- as.numeric((t(X - E.X) %*% solve(vcvMat) %*% (X - E.X)) / N)
  out <- list(V = vcvMat * sigmaSq, sigmaSq = sigmaSq, X = X, E.X. = E.X, N = N)
  class(out) <- "phylogSigmaSq"
  return(out)
  }

logLik.phylogSigmaSq <- function(object, ...) {
    # following O'Meara 2006, eq. 3
    # takes output from sigmaSq
numer <- exp(-0.5 * t(object$X - object$E.X) %*% solve(object$V) %*% (object$X - object$E.X))
  denom <- sqrt((2*pi)^object$N * det(object$V))
  out <- log(numer / denom)
  attr(out, "nobs") <- object$N
  attr(out, "df") <- 2
  class(out) <- "logLik"
  return(out)
  }

library(geiger)
library(phytools)

tr <- pbtree(n = 70)
x.raw <- sim.char(tr, 1)[,,]
x.lambda.5 <- sim.char(rescale(tr, 'lambda', 0.5), 1)[,,]
par(mfrow=c(1,2))
x.vals <- seq(from = 0, to = 1, by = 0.02)
y.raw <- sapply(x.vals, function(w) {
    tr.temp <- rescale(tr, 'lambda', w)
    logLik.phylogSigmaSq(sigmaSq(tr.temp, x.raw))
    }
    )
y.lambda <- sapply(x.vals, function(w) {
    tr.temp <- rescale(tr, 'lambda', w)
    logLik.phylogSigmaSq(sigmaSq(tr.temp, x.lambda.5))
    }
    )

plot(x.vals, y.raw, main = 'original tree', xlab = 'lambda estimate', ylab = 'lnL', type = 'l')
abline(h = max(y.raw)-2, lty = 'dashed')
plot(x.vals, y.lambda, main = 'lambda tree (0.5)', xlab = 'lambda estimate', ylab = 'lnL', type = 'l')
abline(h = max(y.lambda)-2, lty = 'dashed')

censor <- function(tr, taxa = 2) {
  ## tr is a tree of class "phylo"
  ## taxa is a list of taxa defining subtrees, or a single number for the number of subtrees you'd like to define somewhat interactively
  if(class(taxa) == "numeric" && length(taxa) == 1) {
    taxa <- vector("list", taxa)
    for(i in seq(length(taxa))) taxa[[i]] <- select.list(tr$tip.label, multiple = TRUE, title = "Select taxa comprising one subtree")
    }
  if(!identical(sort(unlist(taxa)), sort(tr$tip.label))) warning("Taxa do not represent non-overlapping sets of taxa on the tree")
  treeList <- vector("list",length(taxa))
  for(i in seq(length(taxa))) {
    tipsToDrop <- tr$tip.label[-which(tr$tip.label %in% taxa[[i]])]
    treeList[[i]] <- drop.tip(tr, tipsToDrop)
    }
  return(treeList)
  }

testModels <- function(tr, X, multiRateModels = 1, modelList = NULL, ...) {
  ## tests a single-rate model against a single multiple-rate models
    # -- should be generalized to compare the single-rate model against as many multiple-rate models as you like
  ## currently only the censored option is implemented
  ## Arguments:
  ##  tr = tree of class "phylo"
  ##  X = vector of trait values
  ##  multiRateModels = number of multiple-rate models to test
  ##  modelList = list of list of trees; if NULL, trees are created interactively using censor function
   
  ## 1. Check the data out a bit to catch obvious problems
  vcvMat <- vcv(tr)
  if(length(X) != dim(vcvMat)[1]) stop("This function currently only takes vectors of trait values equal in length to the number of taxa in the tree")
  if(identical(names(X), NULL)) {
    warning("X has no labels; assumed to be ordered the same as VCV matrix")
    names(X) <- dimnames(vcvMat)[[1]]
    }
  else {
    if(!all(names(X) %in% dimnames(vcvMat)[[1]])) warning("X labels not identical to tip labels; traits assumed to be ordered the same as VCV matrix")
    vcvMat <- vcvMat[names(X), names(X)]
    }
  
  ## 2. Set up the models
  treeSets <- vector("list", 1 + multiRateModels)
  treeSets[[1]] <- list(tr) ## necessary to make this a list of length=1 so that the nested loop below works
  maxTrees <- 1
  for(i in 2:(1 + multiRateModels)) {
    if(class(modelList) == 'list') {
      taxa <- length(modelList[[i]])
      treeSets[[i]] <- modelList[[i]] # use trees if provided, or ...
      } ## close if
    else {
      taxa <- as.numeric(select.list(c('2','3','4','5','6','7','8','9','10'), preselect = "2", title = "Select the number of subtrees you would like"))
      treeSets[[i]] <- censor(tr, taxa) # ... do interactive tree creation if no trees are handed in via modelList
      } ## close else
    maxTrees <- max(taxa, maxTrees)
    } ## close i
  
  ## 3. Set up matrix to capture output
  paramsBase <- c("sigma.sq", "root")
  params <- paste(sort(rep(paramsBase, maxTrees)), rep(1:maxTrees,length(paramsBase)), sep = ".")
  results <- matrix(NA, nrow = length(treeSets), ncol = 2 + length(paramsBase)*maxTrees, dimnames = list(1:length(treeSets), c('lnL', 'df', params)))
  for(i in 1:length(treeSets)) {
    results[i, c('lnL', 'df')] <- c(0,0)
    for(j in 1:length(treeSets[[i]])) {
      tr.temp <- treeSets[[i]][[j]]
      sigma.sq.temp <- sigmaSq(vcv(tr.temp), X[tr.temp$tip.label])
      lnL.temp <- logLik(sigma.sq.temp)
      results[i, c('lnL', 'df')] <- results[i, c('lnL', 'df')] + c(lnL.temp, attr(lnL.temp, 'df'))
      results[i, paste(paramsBase, j, sep = ".")] <- c(sigma.sq.temp$sigmaSq,                         ## sigma.sq
                                                       rootState(vcv(tr.temp), X[tr.temp$tip.label])  ## root
                                                       )                                              ## and close that concatenation
    } ## close j
  } ## close i
  results <- list(treeSets = treeSets, modelSummary = results)
  return(results)
  }

k.blomberg <- function(tr, X) {
    N = length(X)
    ones = matrix(1, N)
    C = vcv(tr)
    C.inv = solve(C)
    a.hat = solve(t(ones) %*% C.inv %*% ones) %*% (t(ones) %*% C.inv %*% X)
    a.hat = rep(a.hat, N)
    MSE0 = (t(X - a.hat) %*% (X - a.hat)) / (N - 1)
    MSE = (t(X - a.hat) %*% C.inv %*% (X - a.hat)) / (N - 1)
    MSE.ratio.brown = (1 / (N - 1)) * (sum(diag(C)) - N/sum(C.inv))
    out = list(K = (MSE0 / MSE)/MSE.ratio.brown,
               MSE0 = MSE0,
               MSE = MSE,
               MSE.ratio.brown = MSE.ratio.brown)
    out
    }


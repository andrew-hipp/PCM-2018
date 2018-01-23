---
title       : Generalized Least Squares
subtitle    : PCM Week 3
author      : Andrew Hipp (ahipp@mortonarb.org)
job         :
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      #
widgets     : [mathjax]            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
knit        : slidify::knit2slides
---

<style>
{
  background-color: #FFFFFF
}
</style>

# Remember that in phylogenetic independent contrasts (PIC), our goal was to render the tip states:
1. Independent
2. Equal variance (normalized by variance)

![Demonstration phylogeny, birth-death tree](assets/fig/unnamed-chunk-1-1.png)

---
# We achieve these goals in two ways:
1. <b>Independence:</b> contrasts between the tips are independent of one another.
$$S1 - S2$$

2. <b>Equal variance:</b> dividing contrasts by the standard deviation (the square-root of the branch length separating tips) normalizes contrasts to unit variance
$${S1 - S2} \over \sqrt{V1 + V2}$$

![plot of chunk unnamed-chunk-2](assets/fig/unnamed-chunk-2-1.png)

```
## Error in BOTHlabels(text, sel, XX, YY, adj, frame, pch, thermo, pie, piecol, : object 'layouts' not found
```

---
# We have to do this because in ordinary least squares (OLS), data points are assumed constant and independent of one another. Thus, when we calculate $\beta$ using OLS, we assume that the covariance matrix is <b>I</b>, the identity matrix:

$$
  \begin{align}
    \beta
    & = (X'X)^{-1} X'y \\
    & = (X'IX)^{-1} X'Iy \\
  \end{align}  
$$

where, for our five-taxon tree:

$$
I =
\begin{bmatrix}
  1 & 0 & 0 & 0 & 0 \\
  0 & 1 & 0 & 0 & 0 \\
  0 & 0 & 1 & 0 & 0 \\
  0 & 0 & 0 & 1 & 0 \\
  0 & 0 & 0 & 0 & 1 \\
\end{bmatrix}
$$

This covariance matrix simply indicates that variances are equal for all five tips, and covariances are all zero.

---

# Now let's back up a minute and get our covariance matrix, <b>C</b>. Recall our tree... this time, we'll plot it with branch lengths instead of labels:

![plot of chunk unnamed-chunk-3](assets/fig/unnamed-chunk-3-1.png)

---
# The covariance matrix takes the distance from the tips to the root as the expected variance; these go on the diagonals. The off-diagonals are then the expected covariance, estimated as the amount of shared evolutionary history between each pair of taxa

![plot of chunk unnamed-chunk-4](assets/fig/unnamed-chunk-4-1.png)

$$
  C =
      \begin{bmatrix}
        1.332 & 0.877 & 0.654 & 0.654 & 0 \\
        0.877 & 1.332 & 0.654 & 0.654 & 0 \\
        0.654 & 0.654 & 1.332 & 0.874 & 0 \\
        0.654 & 0.654 & 0.874 & 1.332 & 0 \\
        0 & 0 & 0 & 0 & 1.332 \\
      \end{bmatrix}
$$

---
# Inverting the covariance matrix (using `solve`) yields



$$
    C^{-1} =
    \begin{bmatrix}
      1.447 & -0.749 & -0.207 & -0.207 & 0 \\
      -0.749 & 1.447 & -0.207 & -0.207 & 0 \\
      -0.207 & -0.207 & 1.441 & -0.742 & 0 \\
      -0.207 & -0.207 & -0.742 & 1.441 & 0 \\
      0 & 0 & 0 & 0 & 0.751 \\
    \end{bmatrix}
$$

$C^{-1}$ has an interesting property. Where the column sums of $C$ estimate overall relatedness of each taxon to all others, the column sums of $C^{-1}$ estimate phylogenetic distinctiveness.

---
# Phylogenetic distinctiveness estimated as $C^{-1}$


```r
library(scales)
tr2 <- ladderize(sim.bdtree(n = 20)) # generate a 20 taxon birth-death tree
tr.w <- rescale(colSums(solve(vcv(tr2))), to = c(1, 4))
plot.phylo(tr2, show.tip.label = FALSE)
tiplabels(pch = 19, cex = tr.w, offset = 0.02)
```

![plot of chunk unnamed-chunk-5](assets/fig/unnamed-chunk-5-1.png)

---
# The covariance matrix $C$ is also used to obtain the phylogenetic mean, the ancestral character state under a Brownian motion model:

$$
\hat{a} = (1' C^{-1} 1)^{-1} (1' C^{-1} X)
$$

which is the same as the PIC estimate of the root of the tree. Note that the term $1' C^{-1}$ yields the column sums of $C^{-1}$; multiplying by the row vector $1$ yields the sum of the column sums. Placing this whole thing in the denominator normalizes the phylogenetic mean by $C$, intuitively satisfying given that $\hat{a}$ is really just a weighted mean.

# $C$ is also used to calculate the phylogenetic variance $\sigma^2$:

$$
  \sigma^2 = {{(x - \hat{a})' C^{-1} (x - \hat{a})} \over N}
$$

$$
  V = \sigma^2 C
$$

# and the likelihood, here represented in `R` notation

```
lnL <- dmvnorm(x, rep(a.hat, N), V, log = TRUE) # from mvtnorm package
```

---
# But we're getting a bit ahead of ourselves! These are all the special case of the no-predictor model, which is what we'll use when we are looking at stretching the tree. For purposes of this lecture, it's important to understand the special property of an inverted matrix: $C^{-1}$ has the desirable property of yielding <b>I</b> when it is multiplied by <b>C</b> ($C^{-1}C = I$).


```r
C = vcv(tr)
print(cbind(round(C, 3),
            round(solve(C), 3)
            )
      )
```

```
##       s1    s2    s3    s4    s5     s1     s2     s3     s4    s5
## s1 1.332 0.877 0.654 0.654 0.000  1.447 -0.749 -0.207 -0.207 0.000
## s2 0.877 1.332 0.654 0.654 0.000 -0.749  1.447 -0.207 -0.207 0.000
## s3 0.654 0.654 1.332 0.874 0.000 -0.207 -0.207  1.441 -0.742 0.000
## s4 0.654 0.654 0.874 1.332 0.000 -0.207 -0.207 -0.742  1.441 0.000
## s5 0.000 0.000 0.000 0.000 1.332  0.000  0.000  0.000  0.000 0.751
```

---
# ... and when we multiply these together:


```r
print(round(C %*% solve(C), 3))
```

```
##    s1 s2 s3 s4 s5
## s1  1  0  0  0  0
## s2  0  1  0  0  0
## s3  0  0  1  0  0
## s4  0  0  0  1  0
## s5  0  0  0  0  1
```

---
# Incorporating $C$ into both the numerator and denominator of the least squares estimator generalizes that estimator to cases where $X$ and $Y$ are nonindependent with unequal variances. Thus, where the OLS estimator of $\beta$ was:

$$
  \begin{align}
    \beta & = (X'X)^{-1} X'y \\
    & = (X'IX)^{-1} X'Iy \\
  \end{align}  
$$

the GLS estimator of $\beta$ is:

$$
\beta = (X'CX)^{-1} X'Cy
$$

---
## OKAY!! Enough equations!! How do you actually do this stuff?!

There are a variety of generalized least squares (GLS) implementations floating around, but I'll present the one I've used the most, which comes from the `nlme` package.

Let's start by simulating some data and getting the package opened.


```r
library(geiger)
library(magrittr) # for formatting data
tr2 <- sim.bdtree(n = 50) # simulate a 50-tip birth-death bdtree...
# ... and two characters with r = 0.3, mean = 10
dat <- sim.char(tr2, par = matrix(c(1, 0.3, 0.3, 1), 2, 2), root = 10) %>%
  as.data.frame
names(dat) <- c('x', 'y')
```

---
# ... and let's plot it, just to get a sense of what the data look like.


```r
library(ggplot2, gridExtra)
dat.pic <- data.frame(x = pic(dat$x, tr2), y = pic(dat$y, tr2))
p1 <- qplot(x, y, data = dat.pic) + geom_smooth(method = 'lm') + ggtitle('Independent contrasts')
p2 <- qplot(x, y, data = dat) + geom_smooth(method = 'lm') + ggtitle('Raw data')
grid.arrange(p1, p2, nrow = 1)
```

![Biplot of raw data and independent contrasts](assets/fig/unnamed-chunk-9-1.png)

---
# If we are willing to assume that the model of evolution is correct, that these traits actually evolved according to a Brownian motion model, we can just use GLS without adjusting the branch lengths or modifying $C$ in any way. This is equivalent to PIC on an untransformed tree:


```r
library(nlme, ape)
dat.fits <- list(
  gls = gls(y ~ x, data = dat, correlation = corBrownian(value = 1, phy = tr2)),
  pic = lm(y ~ x + 0, data = dat.pic),
  ols = lm(y ~ x, data = dat)
  )
dat.fits$pic$coefficients = c(NA, dat.fits$pic$coefficients)
sapply(dat.fits, coef)
```

```
##                   gls       pic       ols
## (Intercept) 7.7272827        NA 8.9455430
## x           0.2323257 0.2703705 0.1182874
```

---
# We can get more information from each of these analyses using the `summary` function, which has a separate method for an `lm` object than for a `gls` object.


```r
summary(dat.fits$gls)
```

```
## Generalized least squares fit by REML
##   Model: y ~ x 
##   Data: dat 
##        AIC      BIC    logLik
##   138.8839 144.4975 -66.44196
## 
## Correlation Structure: corBrownian
##  Formula: ~1 
##  Parameter estimate(s):
## numeric(0)
## 
## Coefficients:
##                Value Std.Error  t-value p-value
## (Intercept) 7.727283  1.555181 4.968736   0.000
## x           0.232326  0.134704 1.724713   0.091
## 
##  Correlation: 
##   (Intr)
## x -0.919
## 
## Standardized residuals:
##         Min          Q1         Med          Q3         Max 
## -1.99594091 -0.39455711 -0.07507029  0.38367479  1.75896130 
## 
## Residual standard error: 1.63941 
## Degrees of freedom: 50 total; 48 residual
```

---
# But it may well be that a Brownian motion model doesn't fit our data. One of the lessons of Revell 2010 and Rohlf 2006 is that GLS is BLUE (best linear unbiased estimator), but only if we specify the correlation structure correctly. Let's plot the log-likelihood for our original data, which were evolved on the tree, and one for which the data were not evolved on the tree.


```r
tr3 <- tr2
tr3$tip.label <- sample(tr3$tip.label) # randomizes tip label order
lambdaVals <- seq(from = 0, to = 1, by = 0.02)
## create rescaled trees from 0 to 1, where 0 is no phylogenetic structure, 1 is original:
tr2.set <- lapply(lambdaVals,
                 function(x) geiger:::rescale.phylo(tr2, 'lambda', x))
tr3.set <- lapply(lambdaVals,
                 function(x) geiger:::rescale.phylo(tr3, 'lambda', x))
names(tr2.set) <- names(tr3.set) <- as.character(lambdaVals)
## and now fit all the models:
dat.fits2 <- list(
                   gls = lapply(tr2.set, function(trInd)
                     gls(y ~ x, data = dat, correlation = corBrownian(value = 1, phy = trInd))),
                   gls.rnd = lapply(tr3.set, function(trInd)
                     gls(y ~ x, data = dat, correlation = corBrownian(value = 1, phy = trInd)))
                   )
```

---
# So you have some perspective on what we've just done, here are three example trees, scaled by Pagel's $\lambda$


```r
layout(matrix(1:3, 1))
for(i in c("0", "0.5", "1")) plot(tr2.set[[i]], show.tip.label = F, main = paste('lambda =', i))
```

![plot of chunk unnamed-chunk-13](assets/fig/unnamed-chunk-13-1.png)

---
# Now let's put the results into a couple of data frames, so we can plot lnL against Pagel's $\lambda$. Here we have the log-likelihood plot with a dashed line at $lnL - 2$ for the raw data... this is a commonly used value that approximates the 95% confidence interval when sample sizes are large enough.


```r
tr2.lnL <- data.frame(lambda = lambdaVals, lnL = sapply(dat.fits2$gls, logLik))
tr3.lnL <- data.frame(lambda = lambdaVals, lnL = sapply(dat.fits2$gls.rnd, logLik))
layout(matrix(1:2, 1))
plot(tr2.lnL,  type = 'l', main = "Original data, real lambda = 1")
abline(h = max(tr2.lnL$lnL) - 2, lty = 'dashed')
plot(tr3.lnL,  type = 'l', main = "Shuffled data, real lambda << 1")
abline(h = max(tr3.lnL$lnL) - 2, lty = 'dashed')
```

![plot of chunk unnamed-chunk-14](assets/fig/unnamed-chunk-14-1.png)

---
# When you are doing GLS regression in real time, you can simultaneously fit the model using generalized least squares, and Pagel's $\lambda$ or any other tree scalar using maximum likelihood. Likelihood plots are helpful for seeing what's really going on with your data, but `R` will help you fit these models more smoothly. Let's look at tr3, because that's the one where assuming Brownian motion would be wrong.


```r
# first fit the model with just a Brownian assumption
trs.glsAlt <- list(
  brown.rnd = gls(y ~ x, data = dat, correlation = corBrownian(phy = tr3)),
  pagel.rnd = gls(y ~ x, data = dat, correlation = corPagel(value = 1, phy = tr3)),
  brown.raw = gls(y ~ x, data = dat, correlation = corBrownian(phy = tr2)),
  pagel.raw = gls(y ~ x, data = dat, correlation = corPagel(value = 1, phy = tr2))
  )
```
Now, let's look at the results.

---
# first, the coefficients. Notice that for the raw data, both the Brownian assumption (traditional GLS / PIC) and the model fitting Pagel's $\lambda$ give approximately the same results:


```r
sapply(trs.glsAlt, coef) %>% t
```

```
##           (Intercept)          x
## brown.rnd    8.124274 0.20963551
## pagel.rnd   10.180093 0.03557905
## brown.raw    7.727283 0.23232574
## pagel.raw    7.729468 0.23211213
```

# When we look at the value of Pagel's $\lambda$ for each model, we can see why:


```r
paste('shuffled data lambda =', trs.glsAlt$pagel.rnd$modelStruct,
      '\nraw data lambda =', trs.glsAlt$pagel.raw$modelStruct) %>%
message
```

```
## shuffled data lambda = -0.307742864461561 
## raw data lambda = 1.00056334390611
```

---
## Now you get to try it yourselves! Let's discuss the articles, then you can work on the tutorial.

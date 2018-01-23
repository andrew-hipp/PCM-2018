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

which is the same as the PIC estimate of the root of the tree.

---
# More importantly for our purposes, $C^{-1}$ has the desirable property of yielding <b>I</b> when it is multiplied by <b>C</b> ($C^{-1}C = I$).


```r
C = vcv(tr)
print(cbind(round(C, 3), rep(NA, 5), round(solve(C), 3)))
```

```
##       s1    s2    s3    s4    s5        s1     s2     s3     s4    s5
## s1 1.332 0.877 0.654 0.654 0.000 NA  1.447 -0.749 -0.207 -0.207 0.000
## s2 0.877 1.332 0.654 0.654 0.000 NA -0.749  1.447 -0.207 -0.207 0.000
## s3 0.654 0.654 1.332 0.874 0.000 NA -0.207 -0.207  1.441 -0.742 0.000
## s4 0.654 0.654 0.874 1.332 0.000 NA -0.207 -0.207 -0.742  1.441 0.000
## s5 0.000 0.000 0.000 0.000 1.332 NA  0.000  0.000  0.000  0.000 0.751
```

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


library(ape)
options(repr.plot.width=10, repr.plot.height=4)
par(mar = c(0,0,2,0))
tr <- list(edge = matrix(c(7L, 9L,
                          9L, 10L,
                          10L, 11L,
                          11L, 4L,
                          11L, 5L,
                          10L, 3L, 
                          9L, 6L,
                          7L, 8L,
                          8L, 1L,
                          8L, 2L), ncol = 2, byrow = T, dimnames = list(NULL, c('anc', 'desc'))), 
          edge.length = c(0.527646570653822, 0.0921654720483081, 
                          0.0570646922113782, 0.323123264983397, 
                          0.323123264983397, 0.380187957176582, 
                          0.472353429285534, 0.673568669210027, 
                          0.326431330789973, 0.326431330789973), 
          Nnode = 5L, 
          tip.label = c("s1", "s2", "s3", "s4", "s5", "s6")
          )
class(tr) <- 'phylo'
plot(tr)
edgelabels(text = round(tr$edge.length, 3))

C = round(vcv(tr), 3) # rounding just so it is easier to read
C

V <- function(tr, alpha = 0.1, rescale = T) {
     Tm = max(node.depth.edgelength(tr))
     C = vcv(tr)
     out = (1 / 2*alpha) * exp(-2 * alpha * (Tm-C)) * (1 - exp(-2 * alpha * C))
     if(rescale) out <- out / max(out)
     return(out)
     }
round(V(tr, 0.00001), 3)

library(geiger)
b = 'black'; r = 'red'; y = 'yellow'
regimes = c(b, r, y, y, y, r, r, b, b, b)
par(mar = c(0,0,2,0))
plot(tr, edge.color = regimes, main = 'Original tree')

plot(geiger:::rescale.phylo(tr, 'OU', 3), edge.color = regimes, main = 'O-U tree, alpha = 3')

w.branch = function(tr, alpha = 0.1) {
    t_i.mat <- apply(tr$edge, 1:2, function(x) node.depth.edgelength(tr)[x])
    exp(alpha * t_i.mat[, 2]) - exp(alpha * t_i.mat[, 1]) 
}

print('example weights')
round(cbind(edge.t = node.depth.edgelength(tr)[tr$edge[, 1]],
            edge.l = tr$edge.length,
            w_0.1= w.branch(tr),
            w_1.0 = w.branch(tr, alpha = 1),
            w_5.0 = w.branch(tr, alpha = 5)),
      3)

library(phangorn)
node.list <- lapply(1:length(tr$tip.label), function(i) c(i, Ancestors(tr, i)))
edge.list <- lapply(node.list, function(x) which(tr$edge[, 2] %in% x))
names(edge.list) <- tr$tip.label
print(edge.list)

W <- function(tr, regimes, alpha = 1) {
    Tm = max(node.depth.edgelength(tr))
    out <- matrix(NA, 
                  nrow = length(tr$tip.label), 
                  ncol = length(unique(regimes)) + 1,
                  dimnames = list(sort(tr$tip.label), c(0,sort(unique(regimes))))
                  )
    w.temp <- w.branch(tr, alpha)
    for(i in tr$tip.label) {
    for(k in unique(regimes)) {
        out[i, as.character(k)] <- exp(-alpha * Tm) * sum(w.temp[intersect(edge.list[[i]], which(regimes == k))])
    }
  }
  out[, '0'] <- exp(-alpha * Tm) 
  out
}
W(tr, regimes)

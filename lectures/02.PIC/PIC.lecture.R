## make slides for PIC lecture

slideCounter = 1 ## not currently used
plotAsYouGo = TRUE
w = 1800
h = 1350
un = 'px'
pt = 48
q = 100

## slide1: phylogeny, 2 clades
library(ape)
ntips = 20
tr.base = rcoal(2)
tr.base$edge.length = tr.base$edge.length * 10
tr.base$tip.label = c('NULL1', 'NULL2')
tr <- bind.tree(tr.base, rcoal(ntips), 1)
tr <- bind.tree(tr, rcoal(ntips), 1)
tr$tip.label <- make.unique(tr$tip)
jpeg('slide01.tree.jpg', w, h, un, pt, q)
plot(tr, show.tip.label = F, main = 'Phylogenetic tree with high relative within-clade covariance')
dev.off()

## Slide 2: data on this tree
library(geiger)
sqrtN = 4
xy <- sim.char(tr, matrix(c(2,0,0,2), 2, byrow = T), sqrtN ^ 2)

jpeg('slide02.charsims.jpg', w, h, un, pt, q)
par(mfrow = c(sqrtN, sqrtN),
    oma = c(0,0,3,0) + 0.5,
    mar = c(0,0,0,0) + 0.1)

for(i in 1:(sqrtN ^ 2)) {
  plot(xy[,,i],
       xlab = '', ylab = '',
       xaxt='n', yaxt = 'n',
       pch = 19,
       cex = 1,
       col = c(rep('red', ntips), rep('blue', ntips))
       )
  text(min(xy[,1,i]), max(xy[,2,i]), paste('r =', round(cor(xy[,,i])[1,2], 4),
                                           '| p =', round(cor.test(xy[,1,i], xy[,2,i])$p.val,4)),
                                           cex = 1,
                                           pos = 4,
                                           offset = 0,
                                           #adj = c(-1,1)
                                           )
                                         }
  title('Correlations of normally distributed data,\nno correlation, phylogenetically structured',
       outer = TRUE, cex.sub = 3)
dev.off()

jpeg('slide03.charsims.normal.jpg', w, h, un, pt, q)
par(mfrow = c(sqrtN, sqrtN),
    oma = c(0,0,3,0) + 0.5,
    mar = c(0,0,0,0) + 0.1)

for(i in 1:(sqrtN ^ 2)) {
  xy.temp <- cbind(x = rnorm(ntips * 2, sd = 2),
                   y = rnorm(ntips * 2, sd = 2))
  plot(xy.temp,
       xlab = '', ylab = '',
       xaxt='n', yaxt = 'n',
       pch = 19,
       cex = 1,
       col = c(rep('red', ntips), rep('blue', ntips))
       )
  text(min(xy.temp[, 'x']), max(xy.temp[, 'y']), paste('r =', round(cor(xy.temp)[1,2], 4),
                                           '| p =', round(cor.test(xy.temp[, 'x'], xy.temp[, 'y'])$p.val,4)),
                                           cex = 1,
                                           pos = 4,
                                           offset = 0,
                                           #adj = c(-1,1)
                                           )
                                         }
  title('Correlations of normally distributed data, no tree', outer = TRUE, cex.sub = 3)
   dev.off()

tr.demo1 <-  tr

## slide 4 -- demo tree

ntips2 = 5

tr <- sim.bdtree(n = ntips2)
tr <- reorder(tr, 'postorder')
tr.orig <- tr
tr$node.label = sort(unique(tr$edge[, 1]))
row.names(tr$edge) = paste('V',1:dim(tr$edge)[1], sep = '')
tip.dat <- sample(1:20, ntips2)
names(tip.dat) <- tr$tip.label


plot(tr)
pp.xlim <- get("last_plot.phylo", envir=.PlotPhyloEnv)$x.lim
jpeg('slide04.exampleTree.jpg', w, h, un, pt, q)
doPlotTree <- function(extra.scale = 0.3, node.cex = 0.7, nodesBold = NA, ...) {
  plot(tr,
       x.lim = c(pp.xlim[1], pp.xlim[2] + extra.scale * abs(diff(pp.xlim))),
       ...)
  pp <- get("last_plot.phylo", envir=.PlotPhyloEnv)
  edgelabels(row.names(tr$edge), adj = c(0,1.2), frame = 'n', cex = node.cex)
  nodesFont = rep(1, tr$Nnode)
  nodesFont[nodesBold-ntips2] = 2
  nodelabels(cex = node.cex, font = nodesFont, adj = -1.2, frame = 'n', text = tr$node.label)
  text(pp$x.lim[2], pp$yy[1:ntips2], tip.dat[tr$tip.label], cex = node.cex + 0.1)
}
doPlotTree(main = "Demonstration phylogeny, birth-death tree, with tip data")
dev.off()

## START PIC BEHIND THE SCENES

nodesDone = numeric(0)
nodesToDo <- unique(tr$edge[, 1])
ic <- ic.sd <- rep(NA, length(nodesToDo))
names(ic) <- names(ic.sd) <- sort(as.character(nodesToDo))
means <- c(tip.dat, ic)
names(means)[which(names(means) %in% tr$tip.label)] <-
  match(names(means)[which(names(means) %in% tr$tip.label)], tr$tip.label)

nodeStep = 0

if(plotAsYouGo) {
  jpeg(paste('slide05.a.pic.step', nodeStep, 'jpg', sep='.'), w, h, un, pt, q)
  layout(matrix(1:2, 1))
  doPlotTree(extra.scale = 0.9, show.tip.label = T, node.cex = 0.3, cex = 0.5)
  plot(1, type = 'n', xlim = c(0, 10), ylim = c(0, 10),
       axes = F, xlab = '', ylab = '')
  text(0,10, "Starting tree",
        cex = 0.8,
        pos = 4)
  dev.off()
}
for(working.node in nodesToDo) {
    nodesDone <- c(nodesDone, working.node)
    # do a little book-keeping
    nodeStep = nodeStep + 1
    edges <- which(tr$edge[, 1] == working.node)
    desc.values <- means[as.character(tr$edge[edges, 2])]
    desc.lengths <- tr$edge.length[edges]
    branch.to.rescale <- which(tr$edge[, 2] == working.node)
#    cumulativeEdges = c(cumulativeEdges, edges)

    # 1. make and store the weighted average
    means[as.character(working.node)] <- weighted.mean(desc.values, 1/desc.lengths)

    # 2. do the contrast and its sd
    ic[as.character(working.node)] <- diff(desc.values)
    ic.sd[as.character(working.node)] <- sqrt(sum(desc.lengths))

    # 3. rescale the remaining branch
    tr$edge.length[branch.to.rescale] <-
        tr$edge.length[branch.to.rescale] +
        (desc.lengths[1] * desc.lengths[2]) / sum(desc.lengths)

#    tr.orig$edge <- tr$edge[-cumulativeEdges, ]
#    tr.orig$edge.length <- tr$edge.length[-cumulativeEdges]
#    tr.orig$tip.label %in%
tr$node.label[working.node - ntips2] <-
  round(weighted.mean(desc.values, 1/desc.lengths), 3)

if(plotAsYouGo) {
  jpeg(paste('slide05', letters[nodeStep + 1], 'pic.step', nodeStep, 'jpg', sep='.'), w, h, un, pt, q)
  layout(matrix(1:2, 1))
  doPlotTree(extra.scale = 0.9, show.tip.label = T,
             node.cex = 0.5, cex = 0.5,
             nodesBold = nodesDone)
  plot(1, type = 'n', xlim = c(0, 10), ylim = c(0, 10),
       axes = F, xlab = '', ylab = '')
  text(0, 10, paste('Working tips:', paste(c(tr$tip.label, tr$node.label)[tr$edge[edges, 2]], collapse = ", ")), cex = 0.8, pos = 4)
  text(0,9, paste('New value =', round(weighted.mean(desc.values, 1/desc.lengths), 3)),
        cex = .8, pos = 4)
  text(0,8, paste('Contrast =', round(diff(desc.values), 3)), cex = .8, pos = 4)
  text(0,7, paste('sd =', round(sqrt(sum(desc.lengths)), 3)), cex = .8, pos = 4)
  text(0,6, paste('Normalized contrast =', round(diff(desc.values) / sum(desc.lengths), 3)),
                  cex = 0.8, pos = 4)
  if(working.node != tail(nodesToDo,1)) {
    text(0,4, paste('added to branch =',
                    round((desc.lengths[1] * desc.lengths[2]) / sum(desc.lengths), 3)),
          cex = .8,
          pos = 4)
    text(0,3, paste(row.names(tr$edge)[branch.to.rescale],
                    '=',
                    row.names(tr$edge)[branch.to.rescale],
                    '+',
                    row.names(tr$edge)[edges[1]],
                    '*',
                    row.names(tr$edge)[edges[2]],
                    '/ (',
                    row.names(tr$edge)[edges[1]],
                    '+',
                    row.names(tr$edge)[edges[2]],
                    ')'),
                    cex = 0.5,
                    pos = 4)
    }

  dev.off()
}

} # close working.node


## slide 6 -- demo tree

ntips2 = 8

tr <- sim.bdtree(n = ntips2)
tr <- reorder(tr, 'postorder')
tr.orig <- tr
tr$node.label = sort(unique(tr$edge[, 1]))
row.names(tr$edge) = paste('V',1:dim(tr$edge)[1], sep = '')
tip.dat <- sample(1:20, ntips2)
names(tip.dat) <- tr$tip.label


plot(tr)
pp.xlim <- get("last_plot.phylo", envir=.PlotPhyloEnv)$x.lim
jpeg('slide06.exampleTree.jpg', w, h, un, pt, q)
doPlotTree <- function(extra.scale = 0.3, node.cex = 0.7, nodesBold = NA, ...) {
  plot(tr,
       x.lim = c(pp.xlim[1], pp.xlim[2] + extra.scale * abs(diff(pp.xlim))),
       ...)
  pp <- get("last_plot.phylo", envir=.PlotPhyloEnv)
  edgelabels(row.names(tr$edge), adj = c(0,1.2), frame = 'n', cex = node.cex)
  nodesFont = rep(1, tr$Nnode)
  nodesFont[nodesBold-ntips2] = 2
  nodelabels(cex = node.cex, font = nodesFont, adj = -1.2, frame = 'n', text = tr$node.label)
  text(pp$x.lim[2], pp$yy[1:ntips2], tip.dat[tr$tip.label], cex = node.cex + 0.1)
}
doPlotTree(main = "Demonstration phylogeny, birth-death tree, just a little bigger")
dev.off()

## START PIC BEHIND THE SCENES

nodesDone = numeric(0)
nodesToDo <- unique(tr$edge[, 1])
ic <- ic.sd <- rep(NA, length(nodesToDo))
names(ic) <- names(ic.sd) <- sort(as.character(nodesToDo))
means <- c(tip.dat, ic)
names(means)[which(names(means) %in% tr$tip.label)] <-
  match(names(means)[which(names(means) %in% tr$tip.label)], tr$tip.label)

nodeStep = 0

if(plotAsYouGo) {
  jpeg(paste('slide07.a.pic.step', nodeStep, 'jpg', sep='.'), w, h, un, pt, q)
  layout(matrix(1:2, 1))
  doPlotTree(extra.scale = 0.9, show.tip.label = T, node.cex = 0.3, cex = 0.5)
  plot(1, type = 'n', xlim = c(0, 10), ylim = c(0, 10),
       axes = F, xlab = '', ylab = '')
  text(0,10, "Starting tree",
        cex = 0.8,
        pos = 4)
  dev.off()
}
for(working.node in nodesToDo) {
    nodesDone <- c(nodesDone, working.node)
    # do a little book-keeping
    nodeStep = nodeStep + 1
    edges <- which(tr$edge[, 1] == working.node)
    desc.values <- means[as.character(tr$edge[edges, 2])]
    desc.lengths <- tr$edge.length[edges]
    branch.to.rescale <- which(tr$edge[, 2] == working.node)
#    cumulativeEdges = c(cumulativeEdges, edges)

    # 1. make and store the weighted average
    means[as.character(working.node)] <- weighted.mean(desc.values, 1/desc.lengths)

    # 2. do the contrast and its sd
    ic[as.character(working.node)] <- diff(desc.values)
    ic.sd[as.character(working.node)] <- sqrt(sum(desc.lengths))

    # 3. rescale the remaining branch
    tr$edge.length[branch.to.rescale] <-
        tr$edge.length[branch.to.rescale] +
        (desc.lengths[1] * desc.lengths[2]) / sum(desc.lengths)

#    tr.orig$edge <- tr$edge[-cumulativeEdges, ]
#    tr.orig$edge.length <- tr$edge.length[-cumulativeEdges]
#    tr.orig$tip.label %in%
tr$node.label[working.node - ntips2] <-
  round(weighted.mean(desc.values, 1/desc.lengths), 3)

if(plotAsYouGo) {
  jpeg(paste('slide07', letters[nodeStep + 1], 'pic.step', nodeStep, 'jpg', sep='.'), w, h, un, pt, q)
  layout(matrix(1:2, 1))
  doPlotTree(extra.scale = 0.9, show.tip.label = T,
             node.cex = 0.5, cex = 0.5,
             nodesBold = nodesDone)
  plot(1, type = 'n', xlim = c(0, 10), ylim = c(0, 10),
       axes = F, xlab = '', ylab = '')
  text(0, 10, paste('Working tips:', paste(c(tr$tip.label, tr$node.label)[tr$edge[edges, 2]], collapse = ", ")), cex = 0.8, pos = 4)
  text(0,9, paste('New value =', round(weighted.mean(desc.values, 1/desc.lengths), 3)),
        cex = .8, pos = 4)
  text(0,8, paste('Contrast =', round(diff(desc.values), 3)), cex = .8, pos = 4)
  text(0,7, paste('sd =', round(sqrt(sum(desc.lengths)), 3)), cex = .8, pos = 4)
  text(0,6, paste('Normalized contrast =', round(diff(desc.values) / sum(desc.lengths), 3)),
                  cex = 0.8, pos = 4)
  if(working.node != tail(nodesToDo,1)) {
    text(0,4, paste('added to branch =',
                    round((desc.lengths[1] * desc.lengths[2]) / sum(desc.lengths), 3)),
          cex = .8,
          pos = 4)
    text(0,3, paste(row.names(tr$edge)[branch.to.rescale],
                    '=',
                    row.names(tr$edge)[branch.to.rescale],
                    '+',
                    row.names(tr$edge)[edges[1]],
                    '*',
                    row.names(tr$edge)[edges[2]],
                    '/ (',
                    row.names(tr$edge)[edges[1]],
                    '+',
                    row.names(tr$edge)[edges[2]],
                    ')'),
                    cex = 0.5,
                    pos = 4)
    }

  dev.off()
}

} # close working.node


## slide 8 -- PIC for simulated data
jpeg('slide08.charsims.revisited.jpg', w, h, un, pt, q)
par(mfrow = c(sqrtN, sqrtN),
    oma = c(0,0,3,0) + 0.5,
    mar = c(0,0,0,0) + 0.1)

for(i in 1:(sqrtN ^ 2)) {
  plot(xy[,,i],
       xlab = '', ylab = '',
       xaxt='n', yaxt = 'n',
       pch = 19,
       cex = 1,
       col = c(rep('red', ntips), rep('blue', ntips))
       )
  text(min(xy[,1,i]), max(xy[,2,i]), paste('r =', round(cor(xy[,,i])[1,2], 4),
                                           '| p =', round(cor.test(xy[,1,i], xy[,2,i])$p.val,4)),
                                           cex = 1,
                                           pos = 4,
                                           offset = 0,
                                           #adj = c(-1,1)
                                           )
                                         }
  title('Recall the first simulated data:\nno correlation, phylogenetically structured',
       outer = TRUE, cex.sub = 3)
dev.off()

## slide 9
jpeg('slide09.charsims-PIC.jpg', w, h, un, pt, q)
par(mfrow = c(sqrtN, sqrtN),
    oma = c(0,0,3,0) + 0.5,
    mar = c(0,0,0,0) + 0.1)


for(i in 1:(sqrtN ^ 2)) {
  xy.pic <- cbind(x = pic(xy[tr.demo1$tip.label,1,i], tr.demo1),
                  y = pic(xy[tr.demo1$tip.label,2,i], tr.demo1)
                  )
  plot(xy.pic,
       xlab = '', ylab = '',
       xaxt='n', yaxt = 'n',
       pch = 19,
       cex = 1,
       col = c(rep('red', ntips), rep('blue', ntips))
       )
  text(min(xy.pic[, 'x']), max(xy.pic[, 'y']), paste('r =', round(cor(xy.pic)[1,2], 4),
                                           '| p =', round(cor.test(xy.pic[, 'x'], xy.pic[, 'y'])$p.val,4)),
                                           cex = 1,
                                           pos = 4,
                                           offset = 0,
                                           #adj = c(-1,1)
                                           )
                                         }
  title('Correlations of PICs,\nno correlation, phylogenetically structured',
       outer = TRUE, cex.sub = 3)
dev.off()

## Slide 10: data on this tree
library(geiger)
sqrtN = 4
xy <- sim.char(tr.demo1, matrix(c(2,1.5,1.5,2), 2, byrow = T), sqrtN ^ 2)

jpeg('slide10.charsims-correlated.jpg', w, h, un, pt, q)
par(mfrow = c(sqrtN, sqrtN),
    oma = c(0,0,3,0) + 0.5,
    mar = c(0,0,0,0) + 0.1)

for(i in 1:(sqrtN ^ 2)) {
  plot(xy[,,i],
       xlab = '', ylab = '',
       xaxt='n', yaxt = 'n',
       pch = 19,
       cex = 1,
       col = c(rep('red', ntips), rep('blue', ntips))
       )
  text(min(xy[,1,i]), max(xy[,2,i]), paste('r =', round(cor(xy[,,i])[1,2], 4),
                                           '| p =', round(cor.test(xy[,1,i], xy[,2,i])$p.val,4)),
                                           cex = 1,
                                           pos = 4,
                                           offset = 0,
                                           #adj = c(-1,1)
                                           )
                                         }
  title('Now we simulate phylogenetic data,\nwith correlation. Here is the tips correlation...',
       outer = TRUE, cex.sub = 3)
dev.off()

## slide 11
jpeg('slide11.charsims-PIC-withCorrelation.jpg', w, h, un, pt, q)
par(mfrow = c(sqrtN, sqrtN),
    oma = c(0,0,3,0) + 0.5,
    mar = c(0,0,0,0) + 0.1)


for(i in 1:(sqrtN ^ 2)) {
  xy.pic <- cbind(x = pic(xy[tr.demo1$tip.label,1,i], tr.demo1),
                  y = pic(xy[tr.demo1$tip.label,2,i], tr.demo1)
                  )
  plot(xy.pic,
       xlab = '', ylab = '',
       xaxt='n', yaxt = 'n',
       pch = 19,
       cex = 1,
       col = c(rep('red', ntips), rep('blue', ntips))
       )
  text(min(xy.pic[, 'x']), max(xy.pic[, 'y']), paste('r =', round(cor(xy.pic)[1,2], 4),
                                           '| p =', round(cor.test(xy.pic[, 'x'], xy.pic[, 'y'])$p.val,4)),
                                           cex = 1,
                                           pos = 4,
                                           offset = 0,
                                           #adj = c(-1,1)
                                           )
                                         }
  title('...and here is the correlation of the PICs,\nsame phylogenetically structured data (r = 0.75)',
       outer = TRUE, cex.sub = 3)
dev.off()

## slide 12 -- mapping PIC back into original trait space

jpeg('slide12.compare-PIC-withCorrelation.jpg', w, h, un, pt / 2, q)
layout(matrix(1:2, 1))
par(oma = c(0,0,4,0) + 0.2)
plot(xy[,,sqrtN ^ 2],
  xlab = 'x', ylab = 'y',
  pch = 19,
  cex = 2,
  col = c(rep('red', ntips), rep('blue', ntips)),
  main = 'Raw data with fitted line'
  )

abline(lm(xy[,2,sqrtN ^ 2] ~ xy[,1,sqrtN ^ 2]))

  plot(xy.pic,
    xlab = 'pic(x)', ylab = 'pic(y)',
    pch = 19,
    cex = 2,
    col = c(rep('red', ntips), rep('blue', ntips)),
    main = 'Independent contrasts with fitted line'
    )

abline(lm(xy.pic[,2] ~ xy.pic[,1] - 1), lty = 'dashed', col = 'red')

dev.off()

## slide 13
jpeg('slide13.a.project-PIC-intoOriginalCoordSpace.jpg', w, h, un, pt / 2, q)
layout(matrix(1:2, 1))
par(oma = c(0,0,4,0) + 0.2)


      plot(xy[,,sqrtN ^ 2],
        xlab = 'x', ylab = 'y',
        pch = 19,
        cex = 2,
        col = c(rep('red', ntips), rep('blue', ntips)),
        main = 'Raw data with fitted line'
        )

      abline(lm(xy[,2,sqrtN ^ 2] ~ xy[,1,sqrtN ^ 2]))
      legend(x = 'topleft', legend = 'fitted line, OLS', lty = 'solid', col = 'black')

      rootNode <- length(tr.demo1$tip.label) + 1
      rootAnc <- c(x = as.numeric(ace(xy[,1,sqrtN ^ 2], tr.demo1)$ace[as.character(rootNode)]),
                   y = as.numeric(ace(xy[,2,sqrtN ^ 2], tr.demo1)$ace[as.character(rootNode)]))

                   points(rootAnc['x'], rootAnc['y'], pch = 19, cex = 3)
                   text(rootAnc['x'], rootAnc['y'], 'Root state', pos = 4)

plot(tr.demo1, show.tip.label = F)
rootNode <- length(tr.demo1$tip.label) + 1
rootAnc <- c(x = as.numeric(ace(xy[,1,sqrtN ^ 2], tr.demo1)$ace[as.character(rootNode)]),
             y = as.numeric(ace(xy[,2,sqrtN ^ 2], tr.demo1)$ace[as.character(rootNode)]))
nodelabels(node = rootNode)
pp <- get("last_plot.phylo", envir=.PlotPhyloEnv)
text(pp$xx[as.numeric(rootNode)] + 0.1, pp$yy[as.numeric(rootNode)], paste(
  "Root state from PIC:\n",
  "x =", round(rootAnc['x'], 3),
  '\ny =', round(rootAnc['y'], 3)
  ), cex = 1.6, pos = 4)

  title("To project the PIC regression into original coordinate space,\nregress through the phylogenetic mean (root state)",
        outer = TRUE)


      dev.off()

## slide 13.b -- adding OLS and PIC regression line to same plot
jpeg('slide13.b.project-PIC-intoOriginalCoordSpace.jpg', w, h, un, pt / 2, q)
layout(matrix(1:2, 1))
par(oma = c(0,0,4,0) + 0.2)

plot(xy[,,sqrtN ^ 2],
  xlab = 'x', ylab = 'y',
  pch = 19,
  cex = 2,
  col = c(rep('red', ntips), rep('blue', ntips)),
  main = 'Raw data with fitted line'
  )

  points(rootAnc['x'], rootAnc['y'], pch = 19, cex = 3)
  text(rootAnc['x'], rootAnc['y'], 'Root state', pos = 4)

abline(lm(xy[,2,sqrtN ^ 2] ~ xy[,1,sqrtN ^ 2]))
slope = lm(xy.pic[, 2] ~ xy.pic[, 1] -1)$coef
intercept = rootAnc['y'] - slope * rootAnc['x'] # from y = mx + b... thus b = y - mx

abline(a = intercept, b = slope, lty = 'dashed', col = 'red')
legend('topleft', legend = c('fitted line, OLS', 'fitted line, PIC'),
        lty = c('solid', 'dashed'),
        col = c('black', 'red'))


plot(tr.demo1, show.tip.label = F)
pp <- get("last_plot.phylo", envir=.PlotPhyloEnv)
text(pp$xx[as.numeric(rootNode)] + 0.1, pp$yy[as.numeric(rootNode)], paste(
"Root state from PIC:\n",
"x =", round(rootAnc['x'], 3),
'\ny =', round(rootAnc['y'], 3)
), cex = 1.6, pos = 4)

title("To project the PIC regression into original coordinate space,\nregress through the phylogenetic mean (root state)",
  outer = TRUE)


dev.off()

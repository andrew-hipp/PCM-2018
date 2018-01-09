
setwd('PCM2018-tutorials')
dir.create('tutorial.02') # create a directory for the next tutorial
setwd('tutorial.02')
dir.create('workspace') # create our new workspace
setwd('workspace') # move into the new workspace
load('../../PCM.session01a.Rdata') # load data from last tutorial
print(ls()) # show what objects are here to confirm that we have the objects we need

library(ape)
library(ggtree)
library(phytools)

babyTree <- rtree(5) # creates a random tree with 5 tips
par(mfcol = c(1,3))
#layout(matrix(1:3, 1))
for(trType in c('fan', 'phylo', 'unrooted')) {
    plot(babyTree, trType, main = trType)
    nodelabels()
}

print(babyTree$edge)

print(babyTree$tip.label)

print(round(babyTree$edge.length, 4))

print(unique(dat.sections$Section))
dat.sections$color <- NA
dat.sections$color[grep('Quercus', dat.sections$Section)] <- 'black'
dat.sections$color[dat.sections$Section == 'Virentes'] <- 'green'
dat.sections$color[dat.sections$Section == 'Sadlerianae'] <- 'blue'
dat.sections$color[dat.sections$Section == 'Protobalanus'] <- 'yellow'
dat.sections$color[dat.sections$Section == 'Lobatae'] <- 'red'

## plot the tree first:
plot(tr, show.tip.label = T, cex = 0.4, label.offset = 4)

## now, get the plotting parameters, stick them into an object pp:
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

## now we'll plot the points at just the tips
points(pp$xx[1:length(tr$tip.label)] + 2, 
       pp$yy[1:length(tr$tip.label)], 
       pch = 15, # makes the plotting character a square
       cex = 0.5,
       col = dat.sections[tr$tip.label, 'color']) # orders the colors by the actual tip labels

## and we can add a legend, using the minimum and maximum values of x and y respectively:
legend(x = pp$x.lim[1], y = pp$y.lim[2], 
       legend = c('sect Quercus',
                 'sect Virentes',
                 'sect Sadlerianae',
                 'sect Protobalanus',
                 'sect Lobatae'), 
       pch = c(15),
       col = c('black', 'green', 'blue', 'yellow', 'red'), 
       cex = 0.75, pt.cex = 0.75, 
       title = "Quercus sections"
      )

## first, let's plot the tree in white to just set the plotting parameters
plot(tr, show.tip.label = F, edge.color = 'white')

## then we'll grab those plotting parameters, as we did above
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)

## then we'll replot, setting the rightmost limit of the plot inward so that there is room for annotation on the right side
par(new = T) # this command overplots on the existing plot window rather than creating a new one
plot(tr, show.tip.label = F, x.lim = c(pp$x.lim[1], pp$x.lim[2] + 0.6 * abs(diff(pp$x.lim))))

## now let's plot three characters at different relative distances on the x scale, and in different colors
points(pp$xx[1:length(tr$tip.label)] + diff(pp$x.lim) * 0.1 + scale(dat.bio.means[tr$tip.label, 'bio1']),
       pp$yy[1:length(tr$tip.label)],
       pch = 19,
       cex = 0.4,
       col = 'blue'
       )
points(pp$xx[1:length(tr$tip.label)] + diff(pp$x.lim) * 0.2 + scale(dat.bio.means[tr$tip.label, 'bio4']),
       pp$yy[1:length(tr$tip.label)],
       pch = 19,
       cex = 0.4,
       col = 'red'
       )
points(pp$xx[1:length(tr$tip.label)] + diff(pp$x.lim) * 0.4 + scale(dat.bio.means[tr$tip.label, 'bio12']),
       pp$yy[1:length(tr$tip.label)],
       pch = 19,
       cex = 0.4,
       col = 'black'
       )

## and let's add the SEM to each of our bio12 points;
##   this is actually just a relative SEM, because it is rescaled to unit variance for visualization purposes
segments(x0 = pp$xx[1:length(tr$tip.label)] + diff(pp$x.lim) * 0.4 + 
         scale(dat.bio.means[tr$tip.label, 'bio12']) + 
         scale(dat.bio.sem[tr$tip.label, 'bio12']),
         y0 = pp$yy[1:length(tr$tip.label)],
         x1 = pp$xx[1:length(tr$tip.label)] + diff(pp$x.lim) * 0.4 + 
         scale(dat.bio.means[tr$tip.label, 'bio12']) - 
         scale(dat.bio.sem[tr$tip.label, 'bio12']),
         y1 = pp$yy[1:length(tr$tip.label)]
         )



## when you add the legend, points and lines and other symbols can be intermixed; 
##   just set NA for any that you don't want to show
legend("topleft", legend = c('BIO1 -- mean annual temp', 
                             'BIO4 -- temperature seasonality',
                            'BIO12 -- mean annual precip\n+/- 1 standard error'), 
       pch = c(19,19,19), 
       lty = c(NA,NA,'solid'),
       col = c('blue', 'red', 'black'), 
       cex = 0.6, pt.cex = 0.6, 
       title = "Bioclim variables, rescaled units",
       seg.len = 1.5
      )

dat.sect.vector <- unique(dat.sections$Section)
print(dat.sect.vector)

dat.sect.vector <- dat.sect.vector[-c(4)]

dat.sect.colors <- c('gray75', 'gray95', 'green', 'yellow', 'red')
names(dat.sect.colors) <- dat.sect.vector
print(dat.sect.colors)

dat.sect.mrca <- sapply(dat.sect.vector, function(x) {
    getMRCA(tr, row.names(dat.sections)[dat.sections$Section == x])
})
print(dat.sect.mrca)

dat.sect.desc <- lapply(dat.sect.mrca, getDescendants, tree = tr)
print(dat.sect.desc)

dat.sect.edges <- lapply(dat.sect.desc, function(x) which(tr$edge[, 2] %in% x))
print(dat.sect.edges)

dat.sect.colVect <- rep('black', dim(tr$edge)[1]) # a character as long as the first dimension of the edge matrix
for(i in names(dat.sect.edges)) { # set up a loop to run through the section names
    dat.sect.colVect[dat.sect.edges[[i]] ] <- # index dat.sect.colVect by the edges in element i of the dat.sect.edges list...
       dat.sect.colors[i] # ... and fill those elements with the appropriate color
} # close i
print(dat.sect.colVect)

plot.phylo(tr, 'fan', show.tip.label = TRUE, edge.color = dat.sect.colVect, cex = 0.5)
nodelabels(text = names(dat.sect.mrca), node = dat.sect.mrca, cex = 0.6, bg = dat.sect.colors)

p <- ggtree(tr, 
            layout = 'fan', # makes a circular layout
            open.angle = 5, # leaves room for the labels
            colour = 'gray') # and makes the branches gray so they fade into the backgroun

dat.geog.cat <- dat.geog
for(i in 1:4) dat.geog.cat[[i]] <- factor(dat.geog.cat[[i]])
    
p <- gheatmap(p,
              dat.geog.cat,
             colnames_position = 'bottom',
             font.size = 2, 
             width = 0.15)

p <- p + scale_fill_manual(name = "Geographic region",
                            values = c('gray90', 'black'))

p <- p + theme(legend.position = 'none')

print(p)

dat.pca <- prcomp(dat.bio.means[, grep('bio', dimnames(dat.bio.means)[[2]])], scale = T)
pc.dat <- as.data.frame(cbind(dat.pca$x, dat.sections[row.names(dat.pca$x), ]))

layout(matrix(1:2, 1))
p <- ggplot(pc.dat, aes(PC1, PC2))
#p1 <- p + scale_fill_manual(c('yellow', 'green', 'black', 'orange', 'gray10', 'gray30', 'gray50'))
p1 <- p + geom_point(aes(color = subclade, shape = subclade), size = 3)
p1 <- p1 + ggtitle('Climatic PCA, colored by subclade (geographic region)')

p2 <- p + geom_point(aes(color = Section, shape = Section), size = 3)
p2 <- p2 + ggtitle('Climatic PCA, colored by taxonomic section')

multiplot(p2, p1, ncol = 1)

phylomorphospace(tr, dat.pca$x[, 1:2], label = "off")

save.image(file = '../../PCM.session01b.Rdata') # places this in the PCM tutorials folder

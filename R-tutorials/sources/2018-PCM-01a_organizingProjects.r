
## if you haven't already installed the needed packages, take a moment to do so. Use:
#install.packages('ape', repos='http://cran.us.r-project.org')
#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree")

## Take a moment and do this now. Use the commands:
getwd() # see where we are starting 
dir.create('PCM2018-tutorials')
setwd('PCM2018-tutorials')
dir.create('tutorial.01')
setwd('tutorial.01')
dir.create('scripts')
dir.create('data')
dir.create('workspace')
dir()
setwd('workspace')

## take a moment and create dummy files for this worksession. You'll fill them in as we go along.
writeLines('', '../scripts/01.readData.R')
writeLines('', '../scripts/02.plotData.R')
dir('../scripts/')

path = 'https://raw.githubusercontent.com/andrew-hipp/PCM-2018/master/R-tutorials/DATA/oaks/'
dat.bio <- read.csv(paste(path, 'oak.dat.eco.trimmed.csv', sep = ''), as.is = T)

dim(dat.bio)

head(dat.bio)

dat.bio.bySp <- split(dat.bio, dat.bio$species) # produces a list, one data frame each

dir.create('maps')

library(maps)

for(i in names(dat.bio.bySp)[1:10]) { # this creates an index, i, of the unique species names; 
                                      # we'll just do the first 10
    message(paste('doing map', i)) # gives us a status message so we know something is going on...
    jpeg(paste('maps/', i, '.jpg', sep = ''), 800, 600) # opens a jpg file named by species, 800 x 600 pt
    map() # plots the base map; you can limit the range by lat and long using ylim and xlim
    title(i) # adds a title at the top of our map
    points(dat.bio.bySp[[i]]$longitude, 
           dat.bio.bySp[[i]]$latitude,
           pch = 21, col = 'black', bg = 'red', 
           cex = 2) # plot points, red with black outline, 2x normal size
    dev.off() # close pdf file
    } # close for loop

layout(matrix(1:4, 2, 2)) # sets up a plotting layout, for multiple plots on a single panel 
for(i in names(dat.bio.bySp)[11:14]) {
    map(
        xlim = extendrange(r = range(dat.bio$longitude), f = 0.2),
        ylim = extendrange(r = range(dat.bio$latitude), f = 0.2)
        )
    title(i)
    points(dat.bio.bySp[[i]]$longitude, 
           dat.bio.bySp[[i]]$latitude,
           pch = 21, col = 'black', bg = 'red', 
           cex = 2) # plot points, red with black outline, 2x normal size
} # close for loop

numCols <- grep('lat|long|bio', names(dat.bio), value = T)
dat.bio.means <- t(sapply(dat.bio.bySp, function(x) apply(x[, numCols], 2, mean, na.rm = T)))
head(dat.bio.means)

dat.bio.sem <- t(sapply(dat.bio.bySp, function(x) {
    apply(x[, numCols], 2, sd, na.rm = T) / sqrt(dim(x)[1])
    } # close function
                        ) # close sapply
                 ) # close t
head(dat.bio.sem)

temp.na = names(which(table(dat.bio$species) == 1))
print(temp.na)
dat.bio.sem[temp.na, ]

temp.sd <- apply(t(sapply(dat.bio.bySp, function(x) apply(x[, numCols], 2, sd, na.rm = T))),
    2,
    mean, 
    na.rm = T) # the mean of the standard deviation for all these variables
dat.bio.sem[temp.na, ] <- matrix(temp.sd, length(temp.na), length(temp.sd), byrow = T)
dat.bio.sem[temp.na, ]

# sections for each species:
dat.sections <- read.csv(paste(path, 'sect.species.translate.csv', sep = ''), as.is = T, row.names = 1)

# geography for each species:
dat.geog <- read.csv(paste(path, 'spp.geog.csv', sep = ''), as.is = T, row.names = 1)

# and leaf traits:
dat.lf <- read.delim(paste(path, 'lfPhenology.2016-03-09.jcb.tsv', sep = ''), row.names = 1, as.is = T)

head(dat.sections)

head(dat.geog)

head(dat.lf[, c('lfPhenology', 'References.used'), drop = FALSE])

library(ape)

tr <- read.nexus(paste(path, 'trs.calib.jackknife.4.annotated.tre', sep = ''))
print(tr)
plot(tr, cex = 0.4)

if(identical(row.names(dat.bio.means), row.names(dat.bio.sem)))
    row.names(dat.bio.means) <- row.names(dat.bio.sem) <- paste('Quercus', row.names(dat.bio.means), sep = '_')
print(head(row.names(dat.bio.means)))


spp.intersect <- intersect(tr$tip.label, row.names(dat.bio.means))
spp.intersect <- intersect(spp.intersect, row.names(dat.geog))
spp.intersect <- intersect(spp.intersect, row.names(dat.lf))
spp.intersect <- intersect(spp.intersect, row.names(dat.sections))
length(spp.intersect)

spp.intersect <- Reduce(intersect, list(tr$tip.label,
                                      row.names(dat.bio.means),
                                      row.names(dat.geog),
                                      row.names(dat.lf),
                                      row.names(dat.sections)))
length(spp.intersect)

message('Extra tips : Tree')
print(setdiff(tr$tip.label, spp.intersect))
message('Extra tips : bioclim')
print(setdiff(row.names(dat.bio.means), spp.intersect))
message('Extra tips : geog')
print(setdiff(row.names(dat.geog), spp.intersect))
message('Extra tips : leaf data')
print(setdiff(row.names(dat.lf), spp.intersect))
message('Extra tips : sections')
print(setdiff(row.names(dat.sections), spp.intersect))

tr <- drop.tip(tr, which(!tr$tip.label %in% spp.intersect))

dat.bio.means <- dat.bio.means[spp.intersect, ]
dat.bio.sem <- dat.bio.sem[spp.intersect, ]
dat.geog <- dat.geog[spp.intersect, ]
dat.lf <- dat.lf[spp.intersect, ]
dat.sections <- dat.sections[spp.intersect, ]

print(head(
    cbind(bio = row.names(dat.bio.means), geog = row.names(dat.geog), 
            lf = row.names(dat.lf), sect = row.names(dat.sections)),
    10))


rm(i, path, temp.na, temp.sd)
save.image(file = '../../PCM.session01a.Rdata') # places this in the PCM tutorials folder


library(ape)

prairiePath = ('https://raw.githubusercontent.com/andrew-hipp/PCM-2018/master/R-tutorials/DATA/prairie/')
pr <- function(x, path = prairiePath) paste(path, x, sep = '')

tr.prairie <- read.tree(pr('tree.pruned.tre'))

dat <- list(
    blocks = read.csv(pr('dat.blocksSoilCover.csv'), row.names = 1, as.is = T),
    composition = read.csv(pr('dat.composition.2017.csv'), row.names = 1, as.is = T),
    plotMeta = read.csv(pr('dat.cover.meta.2017.csv'), row.names = 1, as.is = T)
    )
## a little cleanup on two names
names(dat$composition) <- gsub('[.-]', '', names(dat$composition))
tr.prairie$tip.label <- gsub('[.-]', '', tr.prairie$tip.label)
dat$bin <- dat$composition
dat$bin[!is.na(dat$bin)] <- 1
dat$bin[is.na(dat$bin)] <- 0
tr.prairie <- drop.tip(tr.prairie, which(!tr.prairie$tip.label %in% names(dat$bin)))


phyD <- function(phy, roundDigits = 2) {
    out <- c(
        pd = sum(phy$edge.length),
        mpd = mean(as.dist(cophenetic(phy))),
        mntd = mean(apply(cophenetic(phy), 1, function(x) min(x[x > 0])))
        )
    out <- round(out, roundDigits)
    out
}

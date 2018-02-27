
library(geiger)
tr <- sim.bdtree(n = 100)
options(repr.plot.width=10, repr.plot.height=4)
par(mar = c(0,0,0,0))
plot(tr, show.tip.label = F)

library(ggplot2)
library(gridExtra)

dat <- sim.char(tr, matrix(c(1,0.6,0.6,1), 2, 2), nsim = 500)
a = as.data.frame(dat[,,1])
names(a) <- c('x','y')

p1 <- ggplot(a, aes(x=x, y=y)) + geom_point() + geom_smooth(method = 'lm')

a.pic <- data.frame(x_pic = pic(a$x, tr),
                     y_pic = pic(a$y, tr))
p2 <- ggplot(a.pic, aes(x=x_pic, y=y_pic)) + geom_point() + geom_smooth(method = 'lm')
grid.arrange(grobs = list(p1, p2), ncol = 2, main = "Raw data and independent contrasts")

dat.cor <- apply(dat, 3, cor)
dat.pic.cor <- apply(dat, 3, function(x) cor(pic(x[, 1], tr), pic(x[, 2], tr)))
layout(matrix(1:2, 1))
hist(dat.cor, 20, xlab = 'correlation', main = 'Raw data, 500 sims',
     xlim = c(-0.3, 1))
abline(v = c(quantile(dat.cor, c(0.025, 0.975)), mean(dat.cor)), 
       lty = c('dashed', 'dashed', 'solid'), lwd = 2, col = 'red')
legend(-0.3, 1000, c('95% quantile', '', 'Mean'), 
       lty = c('dashed', NA, 'solid'), lwd = 2, col = 'red', 
       bty = 'n', cex = 0.7)
hist(dat.pic.cor, 20, xlab = 'correlation', main = 'Independent contrasts, 500 sims')
abline(v = c(quantile(dat.pic.cor, c(0.025, 0.975)), mean(dat.pic.cor)), 
       lty = c('dashed', 'dashed', 'solid'), lwd = 2, col = 'red')

sigma = 0.58
n = c(1, 2, 5, 10, 20, 50, 100, 200, 500)
dat.sampled <- lapply(n, function(nSamples) {
    apply(dat, 1:3, function(x) mean(rnorm(mean = x, sd = sigma, n = nSamples)))
        })
names(dat.sampled) <- paste('n =', n)

layout(matrix(c(1:9), 3, 3, byrow = T))
options(repr.plot.width=10, repr.plot.height=8)

for(i in names(dat.sampled)) {
    dat.sampled.pic <- apply(dat.sampled[[i]], 3, function(x) cor(pic(x[, 1], tr), pic(x[, 2], tr)))
    a=hist(dat.sampled.pic, 20, xlab = 'correlation, pic', xlim = c(-0.5, 1), main = i)
    abline(v = c(quantile(dat.sampled.pic, c(0.025, 0.975)), mean(dat.sampled.pic)),
           lty = c('dashed', 'dashed', 'solid'), lwd = 2, col = 'red')
    text(-0.5, 0.9*max(a$counts), paste('r =', round(mean(dat.sampled.pic), 3)), pos = 4)
}

options(warn=-1)
dat.noME.lambda <- apply(dat[,,1:20], 3, function(x) fitContinuous(tr, x[,1], model = 'OU'))
a=round(sapply(dat.noME.lambda, function(x) x$opt$alpha), 3)
print(paste('Data known without error; OU alpha:', round(mean(a), 3), '+/-', round(sd(a), 3)))
    
dat.ME5.lambda <- apply(dat.sampled$'n = 5'[,,1:20], 3, function(x) fitContinuous(tr, x[,1], model = 'OU'))
a=round(sapply(dat.ME5.lambda, function(x) x$opt$alpha), 3)
print(paste('Data based on 5 individuals per species; OU alpha:', round(mean(a), 3), '+/-', round(sd(a), 3)))
    
dat.ME1.lambda <- apply(dat.sampled$'n = 2'[,,1:20], 3, function(x) fitContinuous(tr, x[,1], model = 'OU'))
a=round(sapply(dat.ME1.lambda, function(x) x$opt$alpha), 3)
print(paste('Data based on 2 individuals per species; OU alpha:', round(mean(a), 3), '+/-', round(sd(a), 3)))

options(warn=-1)

a=round(sapply(dat.noME.lambda, function(x) x$opt$alpha), 3)
print(paste('Data known without error; OU alpha:', round(mean(a), 3), '+/-', round(sd(a), 3)))
    
dat.ME5.lambda <- apply(dat.sampled$'n = 5'[,,1:20], 3, function(x) fitContinuous(tr, x[,1], model = 'OU', SE = sigma / sqrt(5)))
a=round(sapply(dat.ME5.lambda, function(x) x$opt$alpha), 3)
print(paste('Data based on 5 individuals per species with measurement error; OU alpha:', round(mean(a), 3), '+/-', round(sd(a), 3)))
    
dat.ME1.lambda <- apply(dat.sampled$'n = 2'[,,1:20], 3, function(x) fitContinuous(tr, x[,1], model = 'OU', SE = sigma / sqrt(2)))
a=round(sapply(dat.ME1.lambda, function(x) x$opt$alpha), 3)
print(paste('Data based on 2 individuals per species with measurement error; OU alpha:', round(mean(a), 3), '+/-', round(sd(a), 3)))

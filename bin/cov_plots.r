library(RColorBrewer)
# mean coverge for target and non-target region (and total for chr19)
bedgraphs <- list.files(pattern=".bedgraph")

bgrph <- list()
for (i in 1:length(bedgraphs)) {
    bgrph[[i]] <- read.table(bedgraphs[i])
}

# get summary stats & sd 

avgCov <- lapply(bgrph, '[[', 4)
names(avgCov) <- c("total", "non-target", "target")

avgCovMu <- do.call(c, lapply(avgCov, mean))
avgCovSd <- do.call(c, lapply(avgCov, sd))

# check to see what the mean & sd for total would be, if we took out the 0-depth positions):

summary(avgCov$total[avgCov$total > 0])

# I want the order to be reversed for the plot
muc <- rev(avgCovMu)
sdc <- rev(avgCovSd)
# make the barplot of average coverage values 
png("figs/bar_cov.png", h=1000, w=1000, pointsize=20)
mid <- barplot(muc)
barplot(muc, col = NA, border = NA, axes = FALSE, ylim= c(0, 270))
abline(h = 100, col = "gray60")
abline(h = 150, col = "gray60")
abline(h = 200, col = "gray60")
abline(h = 250, col = "gray60")
barplot(muc, col= t(cols), ylim= c(0, 270), add= T)
arrows(x0=mid, y0=muc-sdc, x1=mid, y1=muc+sdc, code=3, angle=90, length=0.1)
box()
dev.off()


# read hist_all files (i.e. the cumulative portion of the 'coverage -hist' files):
files <- list.files(pattern="hist_all.cov")

labs <- gsub("_hist_all\\.cov", "", files, perl= T)

cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}


# create figure 
cols <- brewer.pal(length(cov), "Dark2")
png("figs/hist_cov.png", h=1000, w=1000, pointsize=20)
plot(cov[[1]][2:1067, 2], cov_cumul[[1]][1:1066], type='n', xlab="Depth", ylab="Fraction of captured bases \u2265 depth", ylim=c(0,1.0), main="Region Coverage")
abline(v = 20, col = "gray60")
abline(v = 50, col = "gray60")
abline(v = 80, col = "gray60")
abline(v = 100, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80))
axis(2, at=c(0.90), labels=c(0.90))
axis(2, at=c(0.50), labels=c(0.50))

# plot the actual data
for (i in 1:length(cov)) points(cov[[i]][2:1067, 2], cov_cumul[[i]][1:1066], type='l', lwd=4, col=rev(cols)[i])
legend("topright", legend=c('target', 'non-target', 'total'), col=cols, lty=1, lwd=5)
dev.off()


# nucleotide composition 

nuc <- list.files(pattern='_tmp.nuc')

nucls <- list()
for(i in 1:length(nuc)){
	nucls[[i]] <- read.table(nuc[i], header= T)
}
names(nucls) <- gsub("_tmp\\.nuc", "", nuc, perl= T)

gcmu <- do.call(c, lapply(lapply(nucls, '[[', 5), mean))
gcsd <- do.call(c, lapply(lapply(nucls, '[[', 5), sd))


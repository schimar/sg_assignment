
library(RColorBrewer)


# read derivates of confusion matrices

derivats <- read.table('confMatDerivat.txt', header= T)

dstrict <- derivats[,c(2,4)]
drlxd <- derivats[,c(3,5)]


cols <- brewer.pal(3, "Set2")

png("../figs/bar_rlxd.png", h=1000, w=1000, pointsize=20)
barplot(as.matrix(drlxd), beside= T, ylim= c(0,100), col=cols, names.arg=c("bbtools", "freebayes"))
abline(h = 20, col = "gray60")
abline(h = 40, col = "gray60")
abline(h = 60, col = "gray60")
abline(h = 80, col = "gray60")
legend('topleft', legend= c("Sensitivity", "Precision", "Specificity"), fill= cols)
box()
dev.off()


# strict variant calling
png("../figs/bar_strict.png", h=1000, w=1000, pointsize=20)
barplot(as.matrix(dstrict), beside= T, ylim= c(0,100), col=cols, names.arg=c("bbtools", "freebayes"))
abline(h = 20, col = "gray60")
abline(h = 40, col = "gray60")
abline(h = 60, col = "gray60")
abline(h = 80, col = "gray60")
legend('topleft', legend= c("Sensitivity", "Precision", "Specificity"), fill= cols)
box()
dev.off()



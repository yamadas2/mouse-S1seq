# Copyright 2021 by Soonjoung Kim

# Author: Soonjoung Kim

# Project: Mouse Resection
# Date: 2021/09/01
# Aim: Scatter plot comparing SPO11-oligo and S1-seq heat


# Source and Library
library(here)
source(here("scripts", "MouseResection", "lib", "MouseResection.R.library.20192018.R"))


# Function


# Environmental setting
options(scipen = 10)

workdir <- here("R", "scatter.plot.comparing.SPO11-oligo.and.S1-seq.heat")
dir.create(workdir)
setwd(workdir)


# Main
########## Load ##########
# SPO11-oligo hotspots
dir.peak <- here("peaks", "mm10", "primary")
peak <- "hotspot.center.SPO11.B6.Cell.mmc2.bed"
hs <- read.table(file.path(dir.peak, peak))
hs.sex.chr <- is.element(hs[, 1], c("chrX", "chrY"))

# SPO11 oligo
indir <- here("matrix.Rdata", "keeneylab", "mm10")
fname <- "hotspot.center.SPO11.B6.Cell.mmc2.bed.SPO11-oligo.unique.RPM.halfwin.5600.matrix.Rdata"
load(file.path(indir, fname))
data.oligo <- list(f = data, r = data)

indir <- here("matrix.Rdata", "wt")
sample.names <- c("wt1", "wt2")

data.s1 <- mean.matrixRdata(indir.path = indir,
                         peak.name = peak, 
                         sample.name.list = sample.names)

########## Scatter plots of hotspot signal ##########
# follow a method of Dr. Julian Lange to make scatter plots
# in Lange et al., 2016, PMID: 27745971

dataset.name <- basename(indir)

# SPO11 oligo signal
data.oligo.sum <- apply(data.oligo[["f"]][, -1000:1000 + 5601], 1, sum)

# resection signal
x.region <- list(f = c(250, 2000), r = c(-2000, -250))
data.s1.sum <- cbind(apply(data.s1[["f"]][, x.region[["f"]][1]:x.region[["f"]][2] + 5601], 1, sum)
                     , apply(data.s1[["r"]][, x.region[["r"]][1]:x.region[["r"]][2] + 5601], 1, sum))
data.s1.sum <- apply(data.s1.sum, 1, sum)
hs.nonzero <- (data.s1.sum > 0)
data.s1.sum.plot.min <- 0.01
sum(data.s1.sum < data.s1.sum.plot.min)
data.s1.sum[data.s1.sum < data.s1.sum.plot.min] <- data.s1.sum.plot.min

pdf(paste("Scatter.plot.oligo.vs.S1-seq", dataset.name, "pdf", sep = "."), height=1.48, width=1.48)
par(mar=c(0.5,0.5,0.5,0.5))
plot(log2(data.oligo.sum), log2(data.s1.sum), xlim=c(0,log2(1200)), ylim=c(log2(data.s1.sum.plot.min), log2(100)), xlab="", ylab="", xaxt="n", yaxt="n", cex=0.1, pch=16, col=rgb(51/255, 51/255, 51/255, 0.3), bty="n")
axis(side=1, at=log2(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)), labels=F, cex.axis=0.5, tck=-0.02, mgp=c(1,0.2,0))
axis(side=2, at=log2(c(as.numeric(1:9 %*% t(10^(-2:1))), 100)), labels=F, cex.axis=0.5, tck=-0.02, mgp=c(1,0.5,0))
dev.off()

cor.test(log2(data.oligo.sum[hs.nonzero]), log2(data.s1.sum[hs.nonzero]))


# central signal
x.cen <- list(f = c(-300, 100), r = c(-100, 300))
data.s1.cen.sum <- cbind(apply(data.s1[["f"]][, x.cen[["f"]][1]:x.cen[["f"]][2] + 5601], 1, sum)
                     , apply(data.s1[["r"]][, x.cen[["r"]][1]:x.cen[["r"]][2] + 5601], 1, sum))
data.s1.cen.sum <- apply(data.s1.cen.sum, 1, sum)
hs.nonzero <- (data.s1.cen.sum > 0)
data.s1.cen.sum.plot.min <- 0.01
sum(data.s1.cen.sum < data.s1.cen.sum.plot.min)
data.s1.cen.sum[data.s1.cen.sum < data.s1.cen.sum.plot.min] <- data.s1.cen.sum.plot.min

pdf(paste("Scatter.plot.oligo.vs.S1-seq.center", dataset.name, "pdf", sep = "."), height=1.48, width=1.48)
par(mar=c(0.5,0.5,0.5,0.5))
plot(log2(data.oligo.sum), log2(data.s1.cen.sum), xlim=c(0,log2(1200+1)), ylim=c(log2(data.s1.cen.sum.plot.min), log2(100)), xlab="", ylab="", xaxt="n", yaxt="n", cex=0.1, pch=16, col=rgb(51/255, 51/255, 51/255, 0.3), bty="n")
axis(side=1, at=log2(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)), labels=F, cex.axis=0.5, tck=-0.02, mgp=c(1,0.2,0))
axis(side=2, at=log2(c(as.numeric(1:9 %*% t(10^(-2:1))), 100)), labels=F, cex.axis=0.5, tck=-0.02, mgp=c(1,0.5,0))
dev.off()

cor.test(log2(data.oligo.sum[hs.nonzero]), log2(data.s1.cen.sum[hs.nonzero]))

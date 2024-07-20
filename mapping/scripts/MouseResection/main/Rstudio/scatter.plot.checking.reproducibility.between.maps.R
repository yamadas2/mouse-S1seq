# Copyright 2021 by Soonjoung Kim

# Author: Soonjoung Kim

# Project: Mouse Resection
# Date: 2024/07/11
# Aim: Scatter plot checking reproducibility between maps


# Source and Library
library(here)
source(here("scripts", "MouseResection", "lib", "MouseResection.R.library.20192018.R"))


# Function


# Environmental setting
options(scipen = 10)

workdir <- here("R", "scatter.plot.checking.reproducibility.between.maps")
dir.create(workdir)
setwd(workdir)


# Main
########## Load ##########
# SPO11-oligo hotspots
dir.peak <- here("peaks", "mm10", "primary")
peak <- "hotspot.center.SPO11.B6.Cell.mmc2.bed"
hs <- read.table(file.path(dir.peak, peak))
hs.sex.chr <- is.element(hs[, 1], c("chrX", "chrY"))

indir <- here("matrix.Rdata", "wt")
sample.names <- c("wt1")

data.s1.1 <- mean.matrixRdata(indir.path = indir,
                         peak.name = peak, 
                         sample.name.list = sample.names)


indir <- here("matrix.Rdata", "wt")
sample.names <- c("wt2")

data.s1.2 <- mean.matrixRdata(indir.path = indir,
                              peak.name = peak, 
                              sample.name.list = sample.names)


########## Scatter plots of hotspot signal ##########
# follow a method of Dr. Julian Lange to make scatter plots
# in Lange et al., 2016, PMID: 27745971

dataset.name <- basename(indir)

# resection signal
x.region <- list(f = c(250, 2000), r = c(-2000, -250))
data.s1.1.sum <- cbind(apply(data.s1.1[["f"]][, x.region[["f"]][1]:x.region[["f"]][2] + 5601], 1, sum)
                     , apply(data.s1.1[["r"]][, x.region[["r"]][1]:x.region[["r"]][2] + 5601], 1, sum))
data.s1.1.sum <- apply(data.s1.1.sum, 1, sum)
data.s1.1.sum.plot.min <- 0.01
sum(data.s1.1.sum < data.s1.1.sum.plot.min)
data.s1.1.sum[data.s1.1.sum < data.s1.1.sum.plot.min] <- data.s1.1.sum.plot.min

data.s1.2.sum <- cbind(apply(data.s1.2[["f"]][, x.region[["f"]][1]:x.region[["f"]][2] + 5601], 1, sum)
                     , apply(data.s1.2[["r"]][, x.region[["r"]][1]:x.region[["r"]][2] + 5601], 1, sum))
data.s1.2.sum <- apply(data.s1.2.sum, 1, sum)
data.s1.2.sum.plot.min <- 0.01
sum(data.s1.2.sum < data.s1.2.sum.plot.min)
data.s1.2.sum[data.s1.2.sum < data.s1.2.sum.plot.min] <- data.s1.2.sum.plot.min

hs.nonzero <- ((data.s1.1.sum > 0) & (data.s1.2.sum > 0))

pdf(paste("Scatter.plot.comparison.between.maps", dataset.name, "pdf", sep = "."), height=1.48, width=1.48)
par(mar=c(0.5,0.5,0.5,0.5))
plot(log2(data.s1.1.sum[hs.nonzero]), log2(data.s1.2.sum[hs.nonzero]), xlim=c(log2(data.s1.2.sum.plot.min), log2(100)), ylim=c(log2(data.s1.2.sum.plot.min), log2(100)), xlab="", ylab="", xaxt="n", yaxt="n", cex=0.1, pch=16, col=rgb(51/255, 51/255, 51/255, 0.3), bty="n")
axis(side=1, at=log2(c(as.numeric(1:9 %*% t(10^(-2:1))), 100)), labels=F, cex.axis=0.5, tck=-0.02, mgp=c(1,0.2,0))
axis(side=2, at=log2(c(as.numeric(1:9 %*% t(10^(-2:1))), 100)), labels=F, cex.axis=0.5, tck=-0.02, mgp=c(1,0.5,0))
dev.off()

cor.test(log2(data.s1.1.sum[hs.nonzero]), log2(data.s1.2.sum[hs.nonzero]))

# central signal
x.cen <- list(f = c(-300, 100), r = c(-100, 300))
data.s1.1.cen.sum <- cbind(apply(data.s1.1[["f"]][, x.cen[["f"]][1]:x.cen[["f"]][2] + 5601], 1, sum)
                     , apply(data.s1.1[["r"]][, x.cen[["r"]][1]:x.cen[["r"]][2] + 5601], 1, sum))
data.s1.1.cen.sum <- apply(data.s1.1.cen.sum, 1, sum)
data.s1.1.cen.sum.plot.min <- 0.01
sum(data.s1.1.cen.sum < data.s1.1.cen.sum.plot.min)
data.s1.1.cen.sum[data.s1.1.cen.sum < data.s1.1.cen.sum.plot.min] <- data.s1.1.cen.sum.plot.min

data.s1.2.cen.sum <- cbind(apply(data.s1.2[["f"]][, x.cen[["f"]][1]:x.cen[["f"]][2] + 5601], 1, sum)
                         , apply(data.s1.2[["r"]][, x.cen[["r"]][1]:x.cen[["r"]][2] + 5601], 1, sum))
data.s1.2.cen.sum <- apply(data.s1.2.cen.sum, 1, sum)
data.s1.2.cen.sum.plot.min <- 0.01
sum(data.s1.2.cen.sum < data.s1.2.cen.sum.plot.min)
data.s1.2.cen.sum[data.s1.2.cen.sum < data.s1.2.cen.sum.plot.min] <- data.s1.2.cen.sum.plot.min

hs.nonzero <- ((data.s1.1.cen.sum > 0) & (data.s1.2.cen.sum > 0))

pdf(paste("Scatter.plot.comparison.between.maps.center", dataset.name, "pdf", sep = "."), height=1.48, width=1.48)
par(mar=c(0.5,0.5,0.5,0.5))
plot(log2(data.s1.1.sum[hs.nonzero]), log2(data.s1.2.cen.sum[hs.nonzero]), xlim=c(log2(data.s1.1.cen.sum.plot.min), log2(100)), ylim=c(log2(data.s1.2.cen.sum.plot.min), log2(100)), xlab="", ylab="", xaxt="n", yaxt="n", cex=0.1, pch=16, col=rgb(51/255, 51/255, 51/255, 0.3), bty="n")
axis(side=1, at=log2(c(as.numeric(1:9 %*% t(10^(-2:1))), 100)), labels=F, cex.axis=0.5, tck=-0.02, mgp=c(1,0.2,0))
axis(side=2, at=log2(c(as.numeric(1:9 %*% t(10^(-2:1))), 100)), labels=F, cex.axis=0.5, tck=-0.02, mgp=c(1,0.5,0))
dev.off()

cor.test(log2(data.s1.1.cen.sum[hs.nonzero]), log2(data.s1.2.cen.sum[hs.nonzero]))

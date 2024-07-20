# Copyright 2024 by Soonjoung Kim

# Author: Soonjoung Kim

# Project: Mouse Resection
# Date: 2024/01/26
# Aim: Plot a resection tract histogram, calculate the mean resection length and plot a mean profile normalized to the resection peak height


# Source and Library
library(here)
source(here("scripts", "MouseResection", "lib", "MouseResection.R.library.20192018.R"))


# Function


# Environmental setting
options(scipen = 10)

workdir <- here("R", "resection.tract.histogram.mean.resection.length.and.normalized.mean.profiles")
dir.create(workdir)
setwd(workdir)


# Main
########## Load ##########
# SPO11-oligo hotspots
dir.peak <- here("peaks", "mm10", "primary")
peak <- "hotspot.center.SPO11.B6.Cell.mmc2.bed"

indir <- here("matrix.Rdata", "wt")
sample.names <- c("wt1", "wt2")

data <- mean.matrixRdata(indir.path = indir,
                         peak.name = peak,
                         sample.name.list = sample.names)


########## Resection histogram and mean plots ##########
dataset.name <- basename(indir)

x <- -5600:5600

# Left and right end positions of the resection peak
offset_l <- 100
offset_r <- 2500

# Left end position of the central peak
offset_s <- -1500

# Mean profile with 100 and 2500 bp lines; original scale; binned in 11 bp
mean.dat <- colMeans(data[["f"]] + data[["r"]][, rev(1:ncol(data[["r"]]))]) / 2
flt <- rep(1/11, 11)
plot(seq(-5.6, 5.6, 0.01), filter(mean.dat, flt)[seq(1, 11201, 10)], xlim = c(-2, 3), type = "l", bty = "n", xlab = "Position relative to hotspot (kb)", ylab = "S1-seq signals")
abline(v = 0.1, col = 2, lty = 2, lwd = 2); abline(v = 2.5, col = 2, lty = 2, lwd = 2)
legend("topright", legend = c("100 bp", "2500 bp"), col = 2, lty = 2, lwd = 2)

# Mean profile with the fixed y-axis; original scale; smoothed with a Hann 151-bp filter
x.plot <- -1000:3000
mean.dat.plot <- filter(mean.dat, hanning.flt(151))[x.plot + 5601]
pdf(paste("plot.offset", dataset.name, "pdf", sep = "."), height=0.74, width=1.7)
par(mar=c(1,1,0,1))
plot(x.plot, mean.dat.plot, type = "l",
     xlim = range(x.plot), ylim = c(0, 0.002),
     bty="n", axes = F, ann = F)
abline(h = 0, lwd = 0.5, lty = 2)
abline(v = offset_s, col = 2, lty = 2, lwd = 0.5)
abline(v = offset_l, col = 2, lty = 2, lwd = 0.5)
abline(v = offset_r, col = 2, lty = 2, lwd = 0.5)
abline(h = mean.dat.plot[offset_r - x.plot[1] + 1], col = 2, lty = 2, lwd = 0.5)
axis(side=1, lwd = 0.67, cex.axis=0.6, tck=-0.02, mgp=c(1,0.2,0))
axis(side=2, lwd = 0.67, cex.axis=0.6, tck=-0.02, mgp=c(1,0.5,0), las = 1)
dev.off()

# Resection tract histogram
mean.dat.offset <- mean.dat[offset_l < x & x <= offset_r]
x.offset <- x[offset_l < x & x <= offset_r]
x.offset.bin <- round(x.offset / 100 + 0.491) * 100
mean.dat.offset.bin <- tapply(mean.dat.offset, x.offset.bin, mean)
bg <- mean.dat.offset.bin[length(mean.dat.offset.bin)]
mean.dat.offset.bin <- mean.dat.offset.bin - bg
mean.dat.offset.bin[mean.dat.offset.bin < 0] <- 0
mean.dat.offset.bin <- mean.dat.offset.bin / sum(mean.dat.offset.bin) * 100
pdf(paste("histogram", dataset.name, "pdf", sep = "."), height=0.74, width=1.7)
par(mar=c(1,1,0,1))
barplot(mean.dat.offset.bin, space = 0, col = "blue", axes = F, ann = F)
abline(h = 0, lwd = 0.67)
axis(side=2, lwd = 0.67, cex.axis=0.6, tck=-0.02, mgp=c(1,0.5,0), las = 1)
dev.off()
(rlm <- sum(as.numeric(names(mean.dat.offset.bin)) * mean.dat.offset.bin) / 100)

# Map parameters, mean resection length and other statistics
flt.n <- 151
mean.dat.smoothed <- filter(mean.dat, hanning.flt(flt.n))
bg <- mean.dat.smoothed[x == offset_r]
mean.dat.smoothed.bg.subtracted <- mean.dat.smoothed - bg 
mean.dat.smoothed.bg.subtracted[mean.dat.smoothed.bg.subtracted < 0] <- 0
res <- mean.dat.smoothed.bg.subtracted[offset_l < x & x <= offset_r]
cen <- mean.dat.smoothed.bg.subtracted[offset_s < x & x <= offset_l]
ratio.dimnames <- list(dataset.name,
                       c("Boundary.left", "Boundary.middle", "Boundary.right",
                         "Mean.resection.length",
                         "Backgound.rpm.level", "Hanning.filter.window.size",
                         "Central.peak.position", "Resection.peak.position",
                         "Central.to.resection.peak.hight.ratio", 
                         "Central.to.resection.peak.area.ratio",
                         "Central.peak.height", "Resection.peak.height",
                         "Central.peak.area", "Resection.peak.area"))
(ratio <- matrix(c(offset_s, offset_l, offset_r, rlm, bg, flt.n,
                   which.max(cen) + offset_s - 1, which.max(res) + offset_l - 1,
                   max(cen) / max(res), sum(cen) / sum(res),
                   max(cen)+bg, max(res)+bg, sum(cen), sum(res)),
                 1, 14, dimnames = ratio.dimnames))
out <- paste("Map.parameters.mean.resectoin.length.and.other.statistics", dataset.name, "txt", sep = ".")
write.table(ratio, out, quote = F, sep = "\t")

# Mean profile smoothed with a Hann 151-bp filter, background subtracted
# and normalized to the resection peak height
scale.factor <- max(res)
mean.dat.smoothed.bg.subtracted.scaled <- mean.dat.smoothed.bg.subtracted / scale.factor

x.plot <- -1000:3000
mean.dat.smoothed.bg.subtracted.scaled.plot <- mean.dat.smoothed.bg.subtracted.scaled[x.plot + 5601]
pdf(paste("plot.offset.height.normalized", dataset.name, "pdf", sep = "."),  height=0.74, width=1.7)
par(mar=c(1,1,0,1))
plot(x.plot, mean.dat.smoothed.bg.subtracted.scaled.plot, type = "l",
     xlim = range(x.plot), ylim = c(0, 1.2),
     bty="n", axes = F, ann = F)
abline(h = 0, lwd = 0.5, lty = 2)
abline(v = offset_s, col = 2, lty = 2, lwd = 0.5)
abline(v = offset_l, col = 2, lty = 2, lwd = 0.5)
abline(v = offset_r, col = 2, lty = 2, lwd = 0.5)
abline(h = mean.dat.smoothed.bg.subtracted.scaled.plot[offset_r - x.plot[1] + 1], col = 2, lty = 2, lwd = 0.5)
axis(side=1, lwd = 0.67, cex.axis=0.6, tck=-0.02, mgp=c(1,0.2,0))
axis(side=2, lwd = 0.67, cex.axis=0.6, tck=-0.02, mgp=c(1,0.5,0), las = 1)
dev.off()

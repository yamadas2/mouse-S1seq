# Copyright 2021 by Soonjoung Kim

# Author: Soonjoung Kim

# Project: Mouse Resection
# Date: 2021/08/18
# Aim: Plot heatmaps


# Source and Library
library(here)
source(here("scripts", "MouseResection", "lib", "MouseResection.R.library.20192018.R"))


# Function


# Environmental setting
options(scipen = 10)

workdir <- here("R", "heatmaps.around.SPO11-oligo.hotspot.centers")
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
sample.names <- c("wt1", "wt2")

data <- mean.matrixRdata(indir.path = indir,
                         peak.name = peak, 
                         sample.name.list = sample.names)

########## Heatmap ##########
# follow a method of Dr. Julian Lange to plot heatmaps
# showing H3K4 asymmetry in Lange et al., 2016, PMID: 27745971

# require:
## signals normalized with 40-bp bins

dataset.name <- basename(indir)
x.image <- list(f = (-2000:2000) + 5601, r = (-2000:2000) + 5601)
strands <- c("f", "r")

data.image <- list(f = t(data[["f"]][, x.image[["f"]]]), r = t(data[["r"]][, x.image[["r"]]]))
scale.factor <- colSums(data.image[["f"]]) + colSums(data.image[["r"]])
data.image.scaled <- list(f = scale(data.image[["f"]], center = F, scale = scale.factor),
                          r = scale(data.image[["r"]], center = F, scale = scale.factor))
data.image.scaled[["f"]][is.na(data.image.scaled[["f"]])] <- 0
data.image.scaled[["r"]][is.na(data.image.scaled[["r"]])] <- 0

bins <- list(f = apply(data.image.scaled[["f"]], 2, function(n) sapply(seq(1, 4001, 40), function(i) mean(n[i:(i + 39)], na.rm = T))),
             r = apply(data.image.scaled[["r"]], 2, function(n) sapply(seq(1, 4001, 40), function(i) mean(n[i:(i + 39)], na.rm = T))))
bins[["f"]][is.na(bins[["f"]])] <- 0
bins[["r"]][is.na(bins[["r"]])] <- 0

# Sort order: hotspot heat
sortorder <- order(hs[, 4])

plotcols <- colorRampPalette(c("white", "blue"))
ibreaks <- quantile(cbind(bins[["f"]], bins[["r"]]), prob = seq(0, 1, length = 11), type = 5)
pdf(paste(dataset.name, "top.strand.ordered.by.hotspot.heat.pdf", sep = "."), height=1.48, width=1)
par(mar=c(0.5, 0.5, 0.5, 0.5))
image(-50:50, 1:ncol(bins[["f"]]), bins[["f"]][, sortorder], ann = F, col = (plotcols(10)), breaks = ibreaks, axes = F)
dev.off()

plotcols <- colorRampPalette(c("white", "red"))
ibreaks <- quantile(cbind(bins[["f"]], bins[["r"]]), prob = seq(0, 1, length = 11), type = 5)
pdf(paste(dataset.name, "bottom.strand.ordered.by.hotspot.heat.pdf", sep = "."), height=1.48, width=1)
par(mar=c(0.5, 0.5, 0.5, 0.5))
image(-50:50, 1:ncol(bins[["r"]]), bins[["r"]][, sortorder], ann = F, col = (plotcols(10)), breaks = ibreaks, axes = F)
dev.off()

# after saving pdf:
# 1) open pdf in illustrator
# 2) export pdf as jpg
# 3) open jpg in illustrator and copy to figure
# 4) use x-axis from other part of figure and resize width of heatmap to fit x-axis

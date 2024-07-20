# Copyright 2019 by Shintaro Yamada

# Author: Shintaro Yamada

# Project: Mouse Resection
# Date: 2019/12/28
# Aim: Convert deepTools-generated matrix.txt.gz to Rdata
#      for visualizing S1-seq profiles at hotspots


# Source and Library


# Function


# Environmental setting
options(scipen = 10)


# Main
args <- commandArgs(trailingOnly = T)


if(length(args) < 2) {
  write("R --vanilla --slave --args <PEAK> <MAT> < file.R", stderr())
  q()
}


PEAK <- args[1]
MAT <- args[2]


# Get peak positions
pk <- read.table(PEAK)
pk.pos <- paste(pk$V1, ":", pk$V2, "-", pk$V3, sep = "")


# Get a tag matrix
mt <- read.table(MAT, skip = 1)
mt.pos <- paste(mt$V1, ":", mt$V2, "-", mt$V3, sep = "")
mt.dat <- as.matrix(mt[, -1 * 1:6])


# Re-order a tag matrix
data <- matrix(0, nrow = nrow(pk), ncol = ncol(mt.dat))
flg <- sapply(mt.pos, function (x) which(x == pk.pos))
data[flg, ] <- mt.dat
data[is.na(data)] = 0
save(data, file = sub("txt.gz", "Rdata", MAT), compress = T)

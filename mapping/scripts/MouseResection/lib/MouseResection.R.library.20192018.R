# Copyright 2019 by Shintaro Yamada

# Author: Shintaro Yamada

# Project: Mouse Resection
# Date: 2019/12/28
# Aim: Functions for mouse resection analysis


# Source and Library


# Function

hanning.flt <- function(n = 51) {
  # Hanning filter function from Sam
  #
  # Arguments
  #   n: hanning filter window size (bp)

  if (n == 1)
    c <- 1
  else {
    n <- n - 1
    c <- 0.5 - 0.5 * cos(2 * pi * (0:n)/n)
  }
  return(c/sum(c))
}


load.matrixRdata <- function(sample.num = NULL
                           , allele = NULL
                           , strand.num = c(1, 2)
                           , half.win = 5600
                           , rpm.norm = T
                           , indir = "./"
                           , indir.rpm = NULL
                           , peak = ""
                           , sample.name.prefix ="DSB"
                           , sample.name = NULL
                           , Rdata.name = NULL
                           , strand.info = T) {
  # Load matrix Rdata
  #
  # Arguments
  #   sample.num: sample number
  #   allele: allele name
  #   strand.num: 1 (forward) or 2 (reverse)
  #   half.win: half window size around peaks
  #   rpm.norm: if true, signals will be normalized to RPM
  #   indir: input file directory
  #   indir.rpm: input file directory of sample.rpm.txt
  #              If not specified, indir.rpm is assigned using indir
  #   peak: peak file name
  #   sample.name.prefix: sample name prefix
  #   sample.name: sample name.
  #                If specified, names are assigned using sample.name
  #   Rdata.name: input Rdata name
  #   strand.info: if true, add strand information to the name

  strand <- c("F", "R")
  if (!is.null(Rdata.name)) {
    if (length(Rdata.name) == 1) {
      Rdata.name <- c(Rdata.name, sub(".F.halfwin", ".R.halfwin", Rdata.name))
    }
    name.p <- paste(peak, "", sep = ".")
    name.s <- paste("", strand[strand.num], "halfwin", half.win, "matrix"
                  , "Rdata", sep = ".")
    name <- sub(name.p, "", sub(name.s, "", basename(Rdata.name[strand.num])))
    file.name <- Rdata.name[strand.num]
  } else {
    allele.name <- ifelse(is.null(allele), "", paste(".", allele, sep = ""))
    name <- ifelse(!is.null(sample.name), sample.name, paste(sample.name.prefix
                 , sprintf("%02d", sample.num), allele.name, sep = ""))
    name.strand <- ifelse(!strand.info, name
                        , paste(name, strand[strand.num], sep = "."))
    file.name <- paste(peak, name.strand, "halfwin", half.win, "matrix"
                     , "Rdata", sep = ".")
  }
  load(file.path(indir, file.name))
  if (rpm.norm) {
    file.name.rpm <- paste(name, "rpm.txt", sep = ".")
    indir.rpm <- ifelse(is.null(indir.rpm), indir, indir.rpm)
    rpm <- read.table(file.path(indir.rpm, file.name.rpm))[, 2] / 1e6
    return(data / rpm)
  } else {
    return(data)
  }
}

mean.matrixRdata <- function(sample.num.list = c()
                           , allele.name = NULL
                           , mat.ini = matrix(0, nrow = 13960, ncol = 11201)
                           , half.win.size = 5600
                           , rpm.norm.each = T
                           , indir.path = "./"
                           , indir.rpm.path = NULL
                           , peak.name = ""
                           , sample.name.prefix.list = "DSB"
                           , sample.name.list = NULL
                           , Rdata.name.list = NULL
                           , strand.specific = T) {
  # A wrapper function of load.matixRdata
  # Load multiple matrix Rdata and return a matrix of their mean
  #
  # Arguments
  #   sample.num.list: a vector of sample numbers
  #   allele.name: an F1 hybrid genome name
  #   mat.ini: initial matrix whose dimention is the same as that of inputs
  #   half.win.size: half window size around peaks
  #   rpm.norm.each: each data will be normalized to RPM before averaging them
  #   indir.path: a path to input file directory
  #   indir.rpm.path: a path to input file directory of sample.rpm.txt
  #                   If not specified, indir.rpm is assigned using indir
  #   peak.name: a peak file name
  #   sample.name.prefix.list: sample name prefix list
  #   sample.name.list: sample name list.
  #                     If specified, names are assigned using sample.name.list
  #   Rdata.name.list: input Rdata name list
  #   strand.specific: if false, return list(f = data, r = data)

  mat.f <- mat.ini
  mat.r <- mat.ini
  if (is.null(sample.num.list)) {
    if (!is.null(sample.name.list)) {
      sample.num.list <- seq(along = sample.name.list)
    } else {
      sample.num.list <- seq(along = Rdata.name.list)
    }
  }
  for (i in seq(along = sample.num.list)) {
    mat.f <- mat.f + load.matrixRdata(sample.num.list[i]
      , allele = allele.name, 1, half.win.size, rpm.norm = rpm.norm.each
      , indir = indir.path, indir.rpm = indir.rpm.path, peak = peak.name
      , sample.name.prefix = ifelse(length(sample.name.prefix.list) > 1
        , sample.name.prefix.list[i], sample.name.prefix.list)
      , sample.name = sample.name.list[i], Rdata.name = Rdata.name.list[[i]]
      , strand.info = strand.specific)
    if (strand.specific) {
      mat.r <- mat.r + load.matrixRdata(sample.num.list[i]
        , allele = allele.name, 2, half.win.size, rpm.norm = rpm.norm.each
        , indir = indir.path, indir.rpm = indir.rpm.path, peak = peak.name
        , sample.name.prefix = ifelse(length(sample.name.prefix.list) > 1
          , sample.name.prefix.list[i], sample.name.prefix.list)
        , sample.name = sample.name.list[i], Rdata.name = Rdata.name.list[[i]]
        , strand.info = strand.specific)
    }
  }

  if (!strand.specific) { mat.r <- t(apply(mat.f, 1, rev)) }
  s.len <- length(sample.num.list)
  data.tmp <- list(f = mat.f / s.len, r = mat.r / s.len)
  return(data.tmp)
}


plot.ind <- function(data.s1
                   , hs.num, y.max
                   , x = (-3000:3000), x.core = (-3000:3000) + 5601
                   , flt = hanning.flt(151)) {
  # plot S1-seq at a specified hotspot
  #
  # Arguments
  #   data.s1: S1-seq matrix
  #   hs.num: hotspot number
  #   x: x position
  #   x.core: sub-region of an input matrix
  #   y.max: max of y-axis
  #   flt: filter argument of a filter function

  fy <- filter(data.s1[["f"]][hs.num, ], flt)[x.core]
  ry <- filter(data.s1[["r"]][hs.num, ], flt)[x.core]
  plot(NA, type = "n", xlim = range(x), ylim = c(0, y.max)
     , bty="n", axes = F, ann = F)
  polygon(c(x, rev(x)), c(numeric(length(x)), rev(ry))
        , col = rgb(1, 0, 0, alpha = 0.5), border = F)
  polygon(c(x, rev(x)), c(numeric(length(x)), rev(fy))
        , col = rgb(0, 0, 1, alpha = 0.5), border = F)
  lines(x, ry, col = rgb(1, 0, 0))
  lines(x, fy, col = rgb(0, 0, 1))
  axis(side=1, labels = F, lwd = 0.5, tck = -0.04, mgp = c(3, 0.2, 0))
  axis(side=2, labels = F, lwd = 0.5, tck = -0.04, mgp = c(3, 0.5, 0))
}


plot.ind.oligo <- function(data.oligo
                         , x = (-3000:3000), x.core = (-3000:3000) + 5601
                         , flt = hanning.flt(151)) {
  # plot SPO11 oligos at a specified hotspot
  #
  # Arguments
  #   data.oligo: SPO11-oligo matrix
  #   x: x position
  #   x.core: sub-region of an input matrix
  #   flt: filter argument of a filter function

  y <- filter(data.oligo, flt)[x.core]
  plot(NA, type = "n", xlim = range(x), ylim = c(0, max(y))
     , bty="n", axes = F, ann = F)
  polygon(c(x, rev(x)), c(numeric(length(x)), rev(y))
        , col = rgb(0.5, 0, 0.5, alpha = 0.5), border = F)
  lines(x, y, col = rgb(0.5, 0, 0.5))
  axis(side=1, labels = F, lwd = 0.5, tck = -0.04, mgp = c(3, 0.2, 0))
  axis(side=2, labels = F, lwd = 0.5, tck = -0.04, mgp = c(3, 0.5, 0))
}


plot.mean.t.and.b <- function(data.s1, data.ctrl
                            , col.s1 = rgb(1, 0, 0), col.ctrl = rgb(0, 0, 1)
                            , x = (-3000:3000), x.core = (-3000:3000) + 5601
                            , flt = hanning.flt(401)) {
  # Plot mean S1-seq on both strands
  #
  # Argument
  #   data.s1: S1-seq matrix
  #   data.ctrl: control S1-seq matrix
  #   col.s1: color for S1-seq
  #   col.ctrl: color for control S1-seq
  #   x: x position
  #   x.core: sub-region of an input matrix
  #   flt: filter argument of a filter function

  fy <- filter(colMeans(data.s1[["f"]]), flt)[x.core]
  ry <- filter(colMeans(data.s1[["r"]]), flt)[x.core]
  fcy <- filter(colMeans(data.ctrl[["f"]]), flt)[x.core]
  rcy <- filter(colMeans(data.ctrl[["r"]]), flt)[x.core]
  y.max <- max(fy, ry, fcy, rcy, na.rm = T)
  plot(NA, type = "n", xlim = range(x), ylim = c(-y.max, y.max)
     , bty="n", axes = F, ann = F)
  lines(x, -rcy, col = col.ctrl)
  lines(x, fcy, col = col.ctrl)
  lines(x, -ry, col = col.s1)
  lines(x, fy, col = col.s1)
  abline(h = 0, lwd = 0.5, lty = 2)
  axis(side=1, lwd = 0.67, cex.axis=0.6, tck=-0.02, mgp=c(1,0.2,0))
  axis(side=2, lwd = 0.67, cex.axis=0.6, tck=-0.02, mgp=c(1,0.5,0), las = 1)
}


pdf.plot.mean.t.and.b <- function(data.s1, data.ctrl
                                , col.s1 = rgb(1, 0, 0), col.ctrl = rgb(0, 0, 1)
                                , x = (-3500:3500), x.core = (-3500:3500) + 5601
                                , flt = hanning.flt(401), dataset.name = "") {
  # Plot mean S1-seq on both strands and create a pdf
  #
  # Argument
  #   data.s1: S1-seq matrix
  #   data.ctrl: control S1-seq matrix
  #   col.s1: color for S1-seq
  #   col.ctrl: color for control S1-seq
  #   x: x position
  #   x.core: sub-region of an input matrix
  #   flt: filter argument of a filter function
  #   dataset.name: characters that will be a part of an output file name

  pdf(paste("plot.mean", dataset.name, "pdf", sep = ".")
    , height=1.48, width=1.7)
  par(mar=c(1,1,0,1))
  plot.mean.t.and.b(data.s1, data.ctrl, col.s1, col.ctrl, x, x.core, flt)
  dev.off()
}


plot.mean.t.and.b.averaged <- function(data.s1, data.ctrl
                                , col.s1 = rgb(1, 0, 0), col.ctrl = rgb(0, 0, 1)
                                , x = (-3000:3000), x.core = (-3000:3000) + 5601
                                , flt = hanning.flt(401)) {
  # Plot mean S1-seq with signals on both strands averaged
  #
  # Argument
  #   data.s1: S1-seq matrix
  #   data.ctrl: control S1-seq matrix
  #   col.s1: color for S1-seq
  #   col.ctrl: color for control S1-seq
  #   x: x position
  #   x.core: sub-region of an input matrix
  #   flt: filter argument of a filter function

  y  <- filter(colMeans(data.s1[["f"]]
        + data.s1[["r"]][, rev(1:ncol(data.s1[["r"]]))]) / 2, flt)[x.core]
  cy <- filter(colMeans(data.ctrl[["f"]]
        + data.ctrl[["r"]][, rev(1:ncol(data.ctrl[["r"]]))]) / 2, flt)[x.core]

  ylim <- range(y, cy, na.rm = T)
  plot(NA, type = "n", xlim = range(x), ylim = ylim, bty="n", axes = F, ann = F)
  lines(x, cy, col = col.ctrl)
  lines(x, y, col = col.s1)
  abline(h = 0, lwd = 0.5, lty = 2)
  axis(side=1, lwd = 0.67, cex.axis=0.6, tck=-0.02, mgp=c(1,0.2,0))
  axis(side=2, lwd = 0.67, cex.axis=0.6, tck=-0.02, mgp=c(1,0.5,0), las = 1)
}


pdf.plot.mean.t.and.b.averaged <- function(data.s1, data.ctrl
                                , col.s1 = rgb(1, 0, 0), col.ctrl = rgb(0, 0, 1)
                                , x = (-3500:3500), x.core = (-3500:3500) + 5601
                                , flt = hanning.flt(401), dataset.name = "") {
  # Plot mean S1-seq with signals on both strands averaged, and create a pdf
  #
  # Argument
  #   data.s1: S1-seq matrix
  #   data.ctrl: control S1-seq matrix
  #   col.s1: color for S1-seq
  #   col.ctrl: color for control S1-seq
  #   x: x position
  #   x.core: sub-region of an input matrix
  #   flt: filter argument of a filter function
  #   dataset.name: characters that will be a part of an output file name

  pdf(paste("plot.mean.strand.averaged", dataset.name, "pdf", sep = ".")
    , height=0.74, width=1.7)
  par(mar=c(1,1,0,1))
  plot.mean.t.and.b.averaged(data.s1, data.ctrl
                           , col.s1, col.ctrl, x, x.core, flt)
  dev.off()
}


plot.mean.t.and.b.averaged.grouped <- function(data.s1, group
                                , cols.s1 = seq(along = group)
                                , x = (-3000:3000), x.core = (-3000:3000) + 5601
                                , flt = hanning.flt(401)) {
  # Plot mean S1-seq for subgroups with signals on both strands averaged
  #
  # Argument
  #   data.s1: S1-seq matrix
  #   group: a vector of factors to specify sub-groups for each row of a matrix
  #   cols.s1: a vector of colors for each sub-group
  #   x: x position
  #   x.core: sub-region of an input matrix
  #   flt: filter argument of a filter function

  group.levels <- levels(group)
  y <- list()
  for (i in seq(along = group.levels)) {
    y[[i]] <- filter(colMeans((data.s1[["f"]] + data.s1[["r"]][
      , rev(1:ncol(data.s1[["r"]]))])[group == group.levels[i], ]) / 2
      , flt)[x.core]
  }
  ylim <- range(y, na.rm = T)
  plot(NA, type = "n", xlim = range(x), ylim = ylim, bty="n", axes = F, ann = F)
  for (i in rev(seq(along = group.levels))) {
    lines(x, y[[i]], col = cols.s1[i])
  }
  axis(side=1, lwd = 0.67, cex.axis=0.6, tck=-0.02, mgp=c(1,0.2,0))
  axis(side=2, lwd = 0.67, cex.axis=0.6, tck=-0.02, mgp=c(1,0.5,0), las = 1)
}


pdf.plot.mean.t.and.b.averaged.grouped <- function(data.s1, group
                                , cols.s1 = seq(along = levels(group))
                                , x = (-3500:3500), x.core = (-3500:3500) + 5601
                                , flt = hanning.flt(401), dataset.name = "") {
  # Plot mean S1-seq with signals on both strands averaged, and create a pdf
  #
  # Argument
  #   data.s1: S1-seq matrix
  #   group: a vector of factors to specify a sub group for each row of a matrix
  #   cols.s1: a vector of colors for each sub-group
  #   x: x position
  #   x.core: sub-region of an input matrix
  #   flt: filter argument of a filter function
  #   dataset.name: characters that will be a part of an output file name

  pdf(paste("plot.mean.strand.averaged.grouped", dataset.name, "pdf", sep = ".")
    , height=1.11, width=1.7)
  par(mar=c(1,1,0,1))
  plot.mean.t.and.b.averaged.grouped(data.s1, group, cols.s1, x, x.core, flt)
  dev.off()
}


filter.s1.mat <- function(data.s1, hs.num) {
  # Make a subset of S1-seq matrix with selected hotspots
  #
  # Argument
  #   data.s1: S1-seq matrix
  #   hs.num: a vector of hotspot numbers (rows) selected

  list(f = data.s1[["f"]][hs.num, ], r = data.s1[["r"]][hs.num, ])
}


subtract.s1.mat <- function(data.s1, data.s1.ctrl) {
  # Subtract data.s1.ctrl from data.s1
  #
  # Argument
  #   data.s1: S1-seq matrix
  #   data.s1.ctrl: S1-seq matrix to subtract

  list(f = data.s1[["f"]] - data.s1.ctrl[["f"]]
     , r = data.s1[["r"]] - data.s1.ctrl[["r"]])
}


local.average.s1 <- function(data.s1, norm.factor) {
  # locally average S1-seq data
  #
  # Argument
  #   data.s1: S1-seq matrix
  #   norm.factor: a vector of normalization values

  list(f = local.average(data.s1[["f"]], factor = norm.factor)
     , r = local.average(data.s1[["r"]], factor = norm.factor))
}


local.average <- function(data, transpose = T, factor = NA) {
  # locally average data
  #
  # Argument
  #   data: data matrix (position x locus)
  #   transpose: if true, data will be treated as locus x position
  #   factor: scale factor

  if (transpose) {
    data <- t(data)
  }

  if (is.na(factor[1])) {
    data <- scale(data, center = F, scale = colMeans(data))
  } else {
    data <- scale(data, center = F, scale = factor)
  }
  data[is.na(data)] <- 0

  if (transpose) {
    return(t(data))
  } else {
    return(data)
  }
}


trim.pos.s1 <- function(data.s1, trim.pos, offset.for.each.row = NA) {
  # trim positions of S1-seq data
  #
  # Argument
  #   data.s1: S1-seq matrix
  #   trim.pos: positions (columns) to output
  #   offset.for.each.row: a vector of values by which positions for each row
  #                        will be shifted

  list(f = trim.pos(data.s1[["f"]], trim.pos, offset.for.each.row)
     , r = trim.pos(data.s1[["r"]], trim.pos, offset.for.each.row))
}


trim.pos <- function(data, trim.pos, offset.for.each.row = NA) {
  # trim positions of data
  #
  # Argument
  #   data: data matrix (position x locus)
  #   shift.center: a vector of values by which positions for each row will be
  #                 shifted

  if (is.na(offset.for.each.row[1])) {
    return(data[, trim.pos])
  } else {
    if (length(offset.for.each.row) != nrow(data)) {
      stop("The length of offset.for.each.row is not the same as the number of \
           data rows\n")
    }
    data.t <- matrix(NA, nrow = nrow(data), ncol = length(trim.pos))
    for (i in seq(along = offset.for.each.row)) {
      data.t[i, ] <- data[i, trim.pos - offset.for.each.row[i]]
    }
    return(data.t)
  }
}


Bin <- function(data, bin.size, fun = sum) {
  # Bin signal
  #
  # Args:
  #   data: Signal matrix
  #         Column represents hotspot numbers and row represents positions
  #   bin.size: Interval length by which to bin signals (bp)
  #   fun: Operation to apply to each bin (default: sum)
  #
  # Return:
  #   Binned signal matrix
  pos <- data[, 1]
  pos.bin.index <- round(pos / bin.size + 0.01) * bin.size
  signal.bin <- tapply(data[, 2], INDEX = pos.bin.index, fun)
  pos.bin <- as.numeric(names(signal.bin))
  data.binned <- cbind(pos.bin, signal.bin)
  return(data.binned)
}

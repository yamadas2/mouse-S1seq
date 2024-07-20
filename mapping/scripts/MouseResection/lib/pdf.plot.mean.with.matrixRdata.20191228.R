# Copyright 2019 by Shintaro Yamada

# Author: Shintaro Yamada

# Project: Mouse Resection
# Date: 2019/12/28
# Aim: Load matrix Rdata and plot mean profiles around peaks in a PDF


# Arguments
args <- commandArgs(trailingOnly = T)


if(length(args) < 2) {
  write("Rscript --vanilla --slave file.R <PEAK> <MATRIX_RDATA_FORWARD> \
    [DIR_MATRIX_RDATA_CONTROL]", stderr())
  q()
}

cargs <- commandArgs()
SDIR <- dirname(dirname(gsub("^--file=", "", cargs[grep("^--file=", cargs)])))
pk <- args[1]
sample.f <- args[2]
dir.ctrl <- args[3]


# Source and Library
source(file.path(SDIR, "lib", "MouseResection.R.library.20192018.R"))


# Function
get.random.data <- function(dat
                          , mat.ini = matrix(NA, nrow = 13960, ncol = 11201)) {
  dat.r <- list(f = mat.ini, r = mat.ini)
  dat.r[["f"]][, -1000:4000 + 5601] <- dat[["f"]][, -5600:-600 + 5601]
  dat.r[["r"]][, -4000:1000 + 5601] <- dat[["r"]][, 600:5600 + 5601]
  return(dat.r)
}


# Environmental setting
options(scipen = 10)


# Main
############### Load matrixRdata files ###############
# Sample
peak            <- read.table(pk)
mat.ini         <- matrix(0, nrow = nrow(peak), ncol = 11201)
indir.path      <- dirname(sample.f)
peak.name       <- basename(pk)
Rdata.name.list <- basename(sample.f)
dat.s <- mean.matrixRdata(mat.ini       = mat.ini
                      , peak.name       = peak.name
                      , indir.path      = indir.path
                      , Rdata.name.list = Rdata.name.list)
name.prefix <- paste(peak.name, "", sep = ".")
name.sample <- gsub(name.prefix, "", gsub(".F.halfwin.5600.matrix.Rdata$", ""
                                          , basename(sample.f)))


# Control
name.posctrl <- "none"
name.negctrl <- "none"
if (is.na(dir.ctrl)) {
  dat.posctrl <- list(f = mat.ini, r = mat.ini)
  dat.negctrl <- list(f = mat.ini, r = mat.ini)
} else {
  name.prefix <- paste(peak.name, "posctrl", "", sep = ".")
  file.ctrl <- list.files(dir.ctrl, pattern = name.prefix, full.names = T)
  if (length(file.ctrl) == 0) {
    dat.posctrl <- list(f = mat.ini, r = mat.ini)
  } else {
    load(file.ctrl)
    dat.posctrl <- data
    name.posctrl <- gsub(name.prefix, "", gsub(".M.halfwin.5600.matrix.Rdata$"
                      , "", basename(file.ctrl)))
  }
  name.prefix <- paste(peak.name, "negctrl", "", sep = ".")
  file.ctrl <- list.files(dir.ctrl, pattern = name.prefix, full.names = T)
  if (length(file.ctrl) == 0) {
    dat.negctrl <- list(f = mat.ini, r = mat.ini)
  } else {
    load(file.ctrl)
    dat.negctrl <- data
    name.negctrl <- gsub(name.prefix, "", gsub(".M.halfwin.5600.matrix.Rdata$"
                      , "", basename(file.ctrl)))
  }
}


# Control subtraction
dat.s.c <- subtract.s1.mat(dat.s, dat.negctrl)
dat.p.c <- subtract.s1.mat(dat.posctrl, dat.negctrl)

# Random sites with control subtraction
dat.r.c <- get.random.data(dat.s.c, mat.ini)

# Garbage collection
rm(mat.ini, data, dat.negctrl)
gc()


############### Plot average profiles ###############
col.list <- c("orange")

x <- (-3500:3500)
x.core <- (-3500:3500) + 5601

dataset.name <- paste(peak.name, name.sample, "and", name.posctrl, sep =".")
pdf.plot.mean.t.and.b(dat.s, dat.posctrl, dataset.name = dataset.name
                    , col.s1 = rgb(1, 0, 0), col.ctrl = rgb(0, 0, 1)
                    , x = x, x.core = x.core, flt = hanning.flt(401))

dataset.name <- paste(peak.name, name.negctrl, "subtracted", name.sample, "and"
                      , name.posctrl, sep =".")
pdf.plot.mean.t.and.b(dat.s.c, dat.p.c, dataset.name = dataset.name
                      , col.s1 = rgb(1, 0, 0), col.ctrl = rgb(0, 0, 1)
                      , x = x, x.core = x.core, flt = hanning.flt(151))

# Top and bottom strand averaged
x <- (-1000:4000)
x.core <- (-1000:4000) + 5601
x.cen <- (-500:500)
x.cen.core <- (-500:500) + 5601

dataset.name <- paste(peak.name, name.sample, "and", name.posctrl, sep =".")
pdf.plot.mean.t.and.b.averaged(dat.s, dat.posctrl, dataset.name = dataset.name
                             , col.s1 = col.list[1], col.ctrl = 1
                             , x = x, x.core = x.core
                             , flt = hanning.flt(151))

dataset.name <- paste(peak.name, name.sample, "and", name.posctrl, "center"
                    , sep =".")
pdf.plot.mean.t.and.b.averaged(dat.s, dat.posctrl, dataset.name = dataset.name
                             , col.s1 = col.list[1], col.ctrl = 1
                             , x = x.cen, x.core = x.cen.core
                             , flt = hanning.flt(51))

dataset.name <- paste(peak.name, name.negctrl, "subtracted", name.sample
                    , sep =".")
pdf.plot.mean.t.and.b.averaged(dat.s.c, dat.r.c, dataset.name = dataset.name
                             , col.s1 = col.list[1], col.ctrl = 1
                             , x = x, x.core = x.core
                             , flt = hanning.flt(151))

dataset.name <- paste(peak.name, name.negctrl, "subtracted", name.sample
                    , "center", sep =".")
pdf.plot.mean.t.and.b.averaged(dat.s.c, dat.r.c, dataset.name = dataset.name
                             , col.s1 = col.list[1], col.ctrl = 1
                             , x = x.cen, x.core = x.cen.core
                             , flt = hanning.flt(51))

dataset.name <- paste(peak.name, name.negctrl, "subtracted", name.sample, "and"
                    , name.posctrl, sep =".")
pdf.plot.mean.t.and.b.averaged(dat.s.c, dat.p.c, dataset.name = dataset.name
                             , col.s1 = col.list[1], col.ctrl = 1
                             , x = x, x.core = x.core
                             , flt = hanning.flt(151))

dataset.name <- paste(peak.name, name.negctrl, "subtracted", name.sample, "and"
                    , name.posctrl, "center", sep =".")
pdf.plot.mean.t.and.b.averaged(dat.s.c, dat.p.c, dataset.name = dataset.name
                             , col.s1 = col.list[1], col.ctrl = 1
                             , x = x.cen, x.core = x.cen.core
                             , flt = hanning.flt(51))

# Sex chromosomes
peak.sex.chr <- is.element(peak[, 1], c("chrX", "chrY"))


dataset.name <- paste(peak.name, name.negctrl, "subtracted", "sexchr"
                    , name.sample, sep =".")
pdf.plot.mean.t.and.b.averaged(filter.s1.mat(dat.s.c, peak.sex.chr)
                             , filter.s1.mat(dat.s.c, !peak.sex.chr)
                             , dataset.name = dataset.name
                             , col.s1 = col.list[1], col.ctrl = 1
                             , x = x, x.core = x.core
                             , flt = hanning.flt(151))

dataset.name <- paste(peak.name, name.negctrl, "subtracted", "sexchr"
                    , name.sample, "center", sep =".")
pdf.plot.mean.t.and.b.averaged(filter.s1.mat(dat.s.c, peak.sex.chr)
                             , filter.s1.mat(dat.s.c, !peak.sex.chr)
                             , dataset.name = dataset.name
                             , col.s1 = col.list[1], col.ctrl = 1
                             , x = x.cen, x.core = x.cen.core
                             , flt = hanning.flt(51))

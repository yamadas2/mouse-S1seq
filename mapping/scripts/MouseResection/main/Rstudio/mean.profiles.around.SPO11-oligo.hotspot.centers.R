# Copyright 2024 by Shintaro Yamada

# Author: Shintaro Yamada

# Project: Mouse Resection
# Date: 2024/07/11
# Aim: plot mean profiles


# Source and Library
library(here)
source(here("scripts", "MouseResection", "lib", "MouseResection.R.library.20192018.R"))


# Function


# Environmental setting
options(scipen = 10)

workdir <- here("R", "mean.profiles.around.SPO11-oligo.hotspot.centers")
dir.create(workdir)
setwd(workdir)


# Main
########## Load ##########
# SPO11-oligo hotspots
dir.peak <- here("peaks", "mm10", "primary")
peak <- "hotspot.center.SPO11.B6.Cell.mmc2.bed"

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


############### Plot average profiles ###############
dataset.name <- basename(indir)

col.list <- c("orange")

# Top and bottom strand averaged
x <- (-1000:4000)
x.core <- (-1000:4000) + 5601
x.cen <- (-500:500)
x.cen.core <- (-500:500) + 5601

plot.name <- paste(peak.name, "comparison.between.maps", dataset.name, sep =".")
pdf.plot.mean.t.and.b.averaged(data.s1.1,
                               data.s1.2,
                               dataset.name = plot.name,
                               col.s1 = col.list[1], col.ctrl = 1,
                               x = x, x.core = x.core,
                               flt = hanning.flt(151))


plot.name <- paste(peak.name, "comparison.between.maps", dataset.name, "center", sep =".")
pdf.plot.mean.t.and.b.averaged(data.s1.1,
                               data.s1.2,
                               dataset.name = plot.name,
                               col.s1 = col.list[1], col.ctrl = 1,
                               x = x.cen, x.core = x.cen.core,
                               flt = hanning.flt(51))

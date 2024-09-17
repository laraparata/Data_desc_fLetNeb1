### OG14 circos plot


#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("GenomicRanges")


## Load packages and files
library(circlize)
library(RColorBrewer)
library(GenomicRanges)
library(data.table)
library(grDevices)

chr <- read.table("OG14_chr.sizes", header = TRUE, colClasses = c("character", "numeric", "numeric"))
gaps <- read.table("OG14_gaps.track", header = TRUE, colClasses = c("character", "numeric", "numeric"))
gc <- read.table("OG14_10000bps.gc", header = TRUE, colClasses = c("character", "numeric", "numeric", "numeric"))
repeats <- read.table("repeat_density.txt", header = TRUE, colClasses = c("character", "numeric", "numeric", "numeric"))

#only keep rows that are chromosomes (e.g., remove entries that contain unloc or are scaffolds)

filtered_chr <- chr[grepl("^SUPER_(?!.*_unloc_)", chr$chr, perl = TRUE), ]
filtered_gaps <- gaps[grepl("^SUPER_(?!.*_unloc_)", gaps$chr, perl = TRUE), ]
filtered_gc <- gc[grepl("^SUPER_(?!.*_unloc_)", gc$chr, perl = TRUE), ]
filtered_repeats <- repeats[grepl("^SUPER_(?!.*_unloc_)", repeats$chr, perl = TRUE), ]

#Make sure that the column idenifying the order of the chromosomes can be numeric, to do this Split the strings in "chr" column by "_"
split_chr <- strsplit(filtered_chr$chr, "_")

# Extract the part after "_" and make this the new chr column. #this needs to be done to all files to ensure they have matching labels to link data.

#chromosomes
filtered_chr$chr <- sapply(split_chr, function(x) ifelse(length(x) > 1, x[[2]], NA))

#gaps
split_gaps <- strsplit(filtered_gaps$chr, "_")
filtered_gaps$chr <- sapply(split_gaps, function(x) ifelse(length(x) > 1, x[[2]], NA))

#GC
split_gc <- strsplit(filtered_gc$chr, "_")
filtered_gc$chr <- sapply(split_gc, function(x) ifelse(length(x) > 1, x[[2]], NA))

#repeats
split_repeats <- strsplit(filtered_repeats$chr, "_")
filtered_repeats$chr <- sapply(split_repeats, function(x) ifelse(length(x) > 1, x[[2]], NA))

#see how many chromosomes you have so you can select the correct number of colours 
str(filtered_chr$chr)

#update this number based off how many chromosomes there are (i.e., this fish has 24)

pal <- colorspace::rainbow_hcl(24)
pal80 <- add_transparency("#96BBDA66", transparency = 0.8)
pal60 <- add_transparency("#96BBDA66", transparency = 0.6)
pal40 <- add_transparency("#96BBDA66", transparency = 0.4)
pal20 <- add_transparency("#96BBDA66", transparency = 0.2)


## Plot
circos.clear()
par(mar = c(1, 1, 1, 1), lwd = 2, cex = 3)
circos.par("start.degree" = 90)
# Track 1 Chromosomes
circos.initializeWithIdeogram(filtered_chr, plotType = NULL)
circos.track(ylim = c(0, 1), bg.col = pal20, panel.fun = function(x, y) {
  filtered_chr2 = CELL_META$sector.numeric.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), filtered_chr2, cex = 1, col = "black",  # increase cex for larger image
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.1, bg.border = NA)

# Track 2 Gaps
circos.genomicTrack(ylim = c(0, 1), filtered_gaps, bg.col = pal40, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, ...)  # setting col does not work
}, track.height = 0.1, bg.border = NA)

# Track 3 GC
circos.genomicTrack(filtered_gc, bg.col = pal60, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, ..., col = "black")
}, track.height = 0.1, bg.border = NA)

# Track 4 rEPEATS
circos.genomicTrack(filtered_repeats, bg.col = pal80, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, ..., col = "black")
}, track.height = 0.1, bg.border = NA)

#rename to match OG
tiff(filename = "OG14_Hap1_circos_blue.tiff", units = "cm", bg = "transparent", width = 80, height = 80, res = 600)
# run plot code
dev.off()

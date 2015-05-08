
# OPEN LIBRARY
library(MASS, quietly=TRUE)
library(car, quietly=TRUE)
library(lsr, quietly=TRUE)
library(lme4, quietly=TRUE)
library(nlme)
library(ape)
library(boot)
library(Hmisc)
library(plyr)

library(metafor)
library(ade4, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(gridExtra, quietly=TRUE)
library(MAc,quietly=TRUE)
library(maptools)
library(downloader)

library(bibtex)
library(RefManageR)

for( f in list.files("R", full.names=TRUE)) {
  source(f)
}

CompileTable <- clean_raw_data()
RawData <- standardise_data(CompileTable)
CompleteData <- build_complete_data(RawData)
IdealData  <- build_ideal_data(CompleteData)
CoordTable <- build_map_data(RawData)

## split datasets by trait - makes a list with named elements
## given subset of data for each trait
list_by_trait <- function(df) {
  split(df, df$trait)
}

GC <- list_by_trait(CompleteData)
GI <- list_by_trait(IdealData)

## Create restricted Dataset: coefficient of correlation are averaged by study and by trait
RC <- list_by_trait(EffectSizeSum(CompleteData))
RI <- list_by_trait(EffectSizeSum(IdealData))

## output directory

dir.create("output", showWarnings =FALSE)

## Export compiled dataset
write.csv(RawData, "output/data.csv", row.names=FALSE)

## Make output
pdf("output/Fig1.pdf",height=3, width=4)
  figure_1(RawData)
dev.off()

pdf("output/Fig2.pdf",height=6,width=5)
  figure_2(GC, GI)
dev.off()

pdf("output/FigA1.pdf", height=6)
  figure_A1(CompleteData)
dev.off()


pdf("output/FigA2.pdf", height=6, width=6)
  figure_allometry()
dev.off()

pdf("output/FigA3.pdf", height=3, width=6)
  figure_map(CoordTable)
dev.off()

pdf("output/FigA4.pdf",height=6,width=6)
  figure_A4(GC)
dev.off()


pdf("output/FigA5.pdf",height=6,width=5)
  figure_A5(GI)
dev.off()


pdf("output/FigA6.pdf",  height=6, width=5)
  figure_A6(RI, RC)
dev.off()


pdf("output/FigA7a.pdf")
  figure_A7(GI, GC, "SLA", c("a","b"))
dev.off()


pdf("output/FigA7b.pdf")
  figure_A7(GI, GC, "WD", c("c","d"))
dev.off()


pdf("output/FigA7c.pdf")
  figure_A7(GI, GC, "Hmax",  c("e","f"))
dev.off()

pdf("output/FigA7d.pdf")
  figure_A7(GI, GC, "Seedmass", c("g","h"))
dev.off()

pdf("output/FigA7e.pdf")
  figure_A7(GI, GC, "Aarea", c("i","j"))
dev.off()


pdf("output/FigA8.pdf")
  figure_A8(GC)
dev.off()


pdf("output/FigA9.pdf")
  figure_A9(GC)
dev.off()

# Reference list

paper <- read.bib("references/paper.bib")
meta <- read.bib("references/metaanalyses.bib")

combined <- c(paper[setdiff(names(paper), names(meta))], meta)
write.bib(combined[[sort(names(combined))]], file = "output/refs.bib")


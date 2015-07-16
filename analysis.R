
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
library(stringr) # Trim whitespace from start and end of string

library(metafor)
library(ade4, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(gridExtra, quietly=TRUE)
library(MAc,quietly=TRUE)
library(maptools)
library(maps)
library(downloader)

library(bibtex)
library(RefManageR)

for( f in list.files("R", full.names=TRUE)) {
  source(f)
}

CompileTable <- clean_raw_data()
RawData <- standardise_data(CompileTable)
CompleteData <- build_complete_data(RawData)
CompleteData_inter <- Build_intersp_complete_data (RawData)
IdealData_inter  <- build_ideal_data(CompleteData)
CoordTable <- build_map_data(RawData)

## split datasets by trait - makes a list with named elements
## given subset of data for each trait
list_by_trait <- function(df) {
  split(df, df$trait)
}

GC <- list_by_trait(CompleteData)
GCi <- list_by_trait(CompleteData_inter)
GIi <- list_by_trait(IdealData_inter)


## Create restricted Dataset: coefficient of correlation are averaged by study and by trait
RC <- list_by_trait(EffectSizeSum(CompleteData))
RCi <- list_by_trait(EffectSizeSum(CompleteData_inter))
RIi <- list_by_trait(EffectSizeSum(IdealData_inter))

## output directory

dir.create("output", showWarnings =FALSE)

## Export compiled dataset
write.csv(RawData, "output/data.csv", row.names=FALSE)

## Make output
pdf("output/Fig1.pdf",height=3, width=4)
  figure_1(RawData)
dev.off()

pdf("output/Fig2.pdf",height=6,width=5)
  figure_2(GCi, GIi)
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
  figure_A5(GIi)
dev.off()


pdf("output/FigA6.pdf",height=6,width=5)
  figure_A6(GC)
dev.off()


pdf("output/FigA7.pdf",  height=6, width=5)
  figure_A7(RIi, RCi)
dev.off()


pdf("output/FigA8a.pdf")
  figure_A8(GIi, GC, "SLA", c("a","b"))
dev.off()


pdf("output/FigA8b.pdf")
  figure_A8(GIi, GC, "WD", c("c","d"))
dev.off()


pdf("output/FigA8c.pdf")
  figure_A8(GIi, GC, "Hmax",  c("e","f"))
dev.off()

pdf("output/FigA8d.pdf")
  figure_A8(GIi, GC, "Seedmass", c("g","h"))
dev.off()

pdf("output/FigA8e.pdf")
  figure_A8(GIi, GC, "Aarea", c("i","j"))
dev.off()


pdf("output/FigA9.pdf")
  figure_A9(GC)
dev.off()


pdf("output/FigA10.pdf")
  figure_A10(GC)
dev.off()


pdf("output/FigA11a.pdf")
figure_A11(GIi, GC, "SLA", c("a","b"))
dev.off()

pdf("output/FigA11b.pdf")
figure_A11(GIi, GC, "WD", c("c","d"))
dev.off()

pdf("output/FigA11c.pdf")
figure_A11(GIi, GC, "Hmax",  c("e","f"))
dev.off()

pdf("output/FigA11d.pdf")
figure_A11(GIi, GC, "Seedmass", c("g","h"))
dev.off()

pdf("output/FigA11e.pdf")
figure_A11(GIi, GC, "Aarea", c("i","j"))
dev.off()

# Reference list

paper <- read.bib("references/paper.bib")
meta <- read.bib("references/metaanalyses.bib")

combined <- c(meta[setdiff(names(meta), names(paper))], paper)
write.bib(combined[[sort(names(combined))]], file = "output/refs.bib")


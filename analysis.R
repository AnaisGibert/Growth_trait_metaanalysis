
library(bibtex)
library(metafor)
library(grid)
library(gridExtra)
library(lme4)
library(ggplot2)

for (f in list.files("R", full.names = TRUE)) {
  source(f)
}

x <- readLines("data/CompileData.csv")
writeLines(x, "data/CompileData.csv")

CompileTable <- clean_raw_data()
RawData <- standardise_data(CompileTable)

CompleteData <- build_complete_data(RawData)
CompleteData_inter <- Build_intersp_complete_data(RawData)

IdealData_inter <- build_ideal_data(CompleteData)

CoordTable <- build_map_data(RawData)

IdealData_rgr <- subset(IdealData_inter, IdealData_inter$growth == "RGR")
IdealData_agr <- subset(IdealData_inter, IdealData_inter$growth == "AbGR")

CompleteData_rgr <- subset(CompleteData_inter, CompleteData_inter$growth == "RGR")
CompleteData_agr <- subset(CompleteData_inter, CompleteData_inter$growth == "AbGR")

## split datasets by trait - makes a list with named elements given subset of
## data for each trait
list_by_trait <- function(df) {
  split(df, df$trait)
}

GCi <- list_by_trait(CompleteData_inter)
GIi <- list_by_trait(IdealData_inter)
GIrgr <- list_by_trait(IdealData_rgr)
GIagr <- list_by_trait(IdealData_agr)

GCrgr <- list_by_trait(CompleteData_rgr)
GCagr <- list_by_trait(CompleteData_agr)

## Create restricted Dataset: coefficient of correlation are averaged by study
## and by trait
RCi <- list_by_trait(EffectSizeSum(CompleteData_inter))
RIi <- list_by_trait(EffectSizeSum(IdealData_rgr))

## output directory

dir.create("output", showWarnings = FALSE)

## Export compiled dataset
write.csv(RawData, "output/data.csv", row.names = FALSE)

## Make output
pdf("output/Fig1.pdf", height = 3, width = 4)
figure_1(RawData)
dev.off()

pdf("output/Fig2.pdf", height = 8, width = 6)
figure_2(CompleteData_inter)
dev.off()


pdf("output/Fig3.pdf", height = 6, width = 5)
figure_3(GIrgr)
dev.off()


pdf("output/Fig4.pdf", height = 6, width = 5)
figure_4(GCi)
dev.off()


## Figure appendix
pdf("output/FigA1.pdf", height = 6, width = 6)  # allometry
figure_A1()
dev.off()

pdf("output/FigA2.pdf", height = 3, width = 6)  #map
figure_A2(CoordTable)
dev.off()

pdf("output/FigA3.pdf", height = 6, width = 6)
figure_A3(GCi)
dev.off()

pdf("output/FigA4.pdf", height = 7, width = 5)
figure_A4(GIi, GIrgr, GIagr)
dev.off()

pdf("output/FigA4.1.pdf", height = 6, width = 5)
figure_3.3(GCi, GIi)
dev.off()

pdf("output/FigA5.pdf", height = 6, width = 5)
figure_A5(RIi, RCi)
dev.off()

pdf("output/FigA6a.pdf")
figure_A6(GCi, "SLA", "WD", c("a", "b"))
dev.off()

pdf("output/FigA6b.pdf")
figure_A6(GCi, "Aarea", "Seedmass", c("c", "d"))
dev.off()

pdf("output/FigA6c.pdf")
figure_A6.2(GCi, "Hmax", "e")
dev.off()

pdf("output/FigA7.pdf")
figure_A7(GIi)
dev.off()

pdf("output/FigA8.pdf")
figure_A8(GCi)
dev.off()

pdf("output/FigA9.pdf")
figure_A9(GIi)
dev.off()


# Reference list
paper <- read.bib("references/paper.bib")
meta <- read.bib("references/metaanalyses.bib")

combined <- c(meta[setdiff(names(meta), names(paper))], paper)
write.bib(combined[[sort(names(combined))]], file = "output/refs.bib")

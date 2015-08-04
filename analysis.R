
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
IdealData_rgr <- subset(IdealData_inter,IdealData_inter$growth=="RGR")
IdealData_agr <- subset(IdealData_inter,IdealData_inter$growth=="AbGR")
## split datasets by trait - makes a list with named elements
## given subset of data for each trait
list_by_trait <- function(df) {
  split(df, df$trait)
}

GC <- list_by_trait(CompleteData)
GCi <- list_by_trait(CompleteData_inter)
GIi <- list_by_trait(IdealData_inter)
GIrgr <- list_by_trait(IdealData_rgr)
GIagr <- list_by_trait(IdealData_agr)



## Create restricted Dataset: coefficient of correlation are averaged by study and by trait
RC <- list_by_trait(EffectSizeSum(CompleteData))
RCi <- list_by_trait(EffectSizeSum(CompleteData_inter))
RIi <- list_by_trait(EffectSizeSum(IdealData_inter))

## output directory

dir.create("output", showWarnings =FALSE)

## Export compiled dataset
write.csv(RawData, "output/data.csv", row.names=FALSE)

## Make output
pdf("output/FigA1.pdf", height=6) # In the main text, need to be changed
figure_A1(CompleteData_inter)
dev.off()

pdf("output/Fig1.pdf",height=3, width=4) # ok
  figure_1(RawData)
dev.off()


pdf("output/Fig2.pdf",height=6,width=5) # ok
  figure_2(GCi, GIi)
dev.off()


pdf("output/FigA4.pdf",height=6,width=6) # maybe in the main text
figure_A4(GCi)
dev.off()


pdf("output/FigA2.pdf", height=6, width=6) # ok , maybe need to put only the line address in the text
  figure_allometry()
dev.off()

pdf("output/FigA3.pdf", height=3, width=6)# ok
  figure_map(CoordTable)
dev.off()


pdf("output/FigA5.pdf",height=6,width=5) # need to be change, both graph on the same page
figure_A5(GIi)
figure_A5b(GIrgr, GIagr)
dev.off()


# pdf("output/FigA6.pdf",height=6,width=5)
#   figure_A6(GCi)
# dev.off()


pdf("output/FigA7.pdf",  height=6, width=5) # ok, heterogenity
  figure_A7(RIi, RCi)
dev.off()


pdf("output/FigA8.pdf")  # ok, test other class
 figure_8(GCi, GIi)
dev.off()
# 
pdf("output/FigA8.b.pdf")
figure_A8(GCi, "Aarea","Seedmass", c("c","d"))
dev.off()
# 
# pdf("output/FigA8.c.pdf")
# figure_A82(GCi, "Hmax", "e")
# dev.off()


pdf("output/FigA9.pdf") # ok funnel plot
  figure_A9(GC)
dev.off()


pdf("output/FigA10a.pdf") # a mixer 
  figure_A10a(GC)
  figure_A10b(GC)
dev.off()

# pdf("output/FigA11.pdf")
# figure_A12(GIi)
# 
# dev.off()

pdf("output/FigA11.pdf")
figure_A12(GIi)
dev.off()


# pdf("output/FigA11a.pdf")
# figure_A11(GIi, GC, "SLA", c("a","b"))
# dev.off()
# 
# pdf("output/FigA11b.pdf")
# figure_A11(GIi, GC, "WD", c("c","d"))
# dev.off()
# 
# pdf("output/FigA11c.pdf")
# figure_A11(GIi, GC, "Hmax",  c("e","f"))
# dev.off()
# 
# pdf("output/FigA11d.pdf")
# figure_A11(GIi, GC, "Seedmass", c("g","h"))
# dev.off()
# 
# pdf("output/FigA11e.pdf")
# figure_A11(GIi, GC, "Aarea", c("i","j"))
# dev.off()

# Heterogeneity
fun_HeterogenityI2(GIi[["SLA"]], plant.stage="NA") 
fun_HeterogenityI2(GIi[["SLA"]], plant.stage="juvenile") 
fun_HeterogenityI2(GIi[["SLA"]], plant.stage="sapling") 
fun_HeterogenityI2(GIi[["SLA"]], plant.stage="adult") 

fun_HeterogenityI2(GIi[["WD"]], plant.stage="NA") 
# fun_HeterogenityI2(GIi[["WD"]], plant.stage="juvenile") 
# fun_HeterogenityI2(GIi[["WD"]], plant.stage="sapling") 
fun_HeterogenityI2(GIi[["WD"]], plant.stage="adult")

fun_HeterogenityI2(GIi[["Hmax"]], plant.stage="NA") 
# fun_HeterogenityI2(GIi[["Hmax"]], plant.stage="juvenile") 
# fun_HeterogenityI2(GIi[["Hmax"]], plant.stage="sapling") 
fun_HeterogenityI2(GIi[["Hmax"]], plant.stage="adult")
# 
fun_HeterogenityI2(GIi[["Seedmass"]], plant.stage="NA") 
fun_HeterogenityI2(GIi[["Seedmass"]], plant.stage="juvenile") 
# fun_HeterogenityI2(GIi[["Seedmass"]], plant.stage="sapling") 
# fun_HeterogenityI2(GIi[["Seedmass"]], plant.stage="adult")
fun_HeterogenityI2(GIi[["Aarea"]], plant.stage="NA") 

# Reference list
paper <- read.bib("references/paper.bib")
meta <- read.bib("references/metaanalyses.bib")

combined <- c(meta[setdiff(names(meta), names(paper))], paper)
write.bib(combined[[sort(names(combined))]], file = "output/refs.bib")


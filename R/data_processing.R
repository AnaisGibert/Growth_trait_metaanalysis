## split datasets by trait - makes a list with named elements given subset of
## data for each trait
list_by_trait <- function(df) {
  split(df, df$trait)
}

subset_growth <- function(data, measure) {
  subset(data, data$growth == measure)
}


clean_raw_data <- function(filename = "data/CompileData.csv") {

  Table <- read.csv(filename, sep = ";", dec = ",")

  for (f in c("id", "ref", "idcor", "authors", "stress", "RGR", "trait", "coef.type",
    "experiment.type", "veg.type", "measure.size", "trait.stage", "stage",
    "RGR.stage", "relation.sign", "stage.simi", "life.form")) {
    Table[[f]] <- as.factor(Table[[f]])
  }

  for (f in c("year", "size.min", "size.max", "size.mean.range", "sample.size",
    "nb.sp.reported", "coef", "RGR.min", "RGR.max", "trait.min", "trait.max")) {
    Table[[f]] <- as.numeric(Table[[f]])
  }

  ### Transformation of raw data ####

  ## Transformation LMA en SLA
  fun_LMA <- function(x) (x * -1)

  for (i in 1:length(Table$coef)) {
    if (Table$trait[i] == "LMA") {
      Table$coef[i] <- fun_LMA(x = Table$coef[i])
    }
  }
  Table$trait[Table$trait == "LMA"] <- "SLA"

  # Transformation seedvolume in seedmass
  Table$trait[Table$trait == "SeedVolume"] <- "Seedmass"

  # Transformation of correlation sign in correlation significance
  Table["sign"] <- "NA"
  Table$sign[Table$relation.sign == "P"] <- "S"
  Table$sign[Table$relation.sign == "N"] <- "S"
  Table$sign[Table$relation.sign == "O"] <- "NS"

  # Add variable for growth and growth measurmemnt
  levels(Table$RGR)

  Table["growth"] <- "NA"
  Table$growth[Table$RGR %in% c("RGR(?)", "RGR(CSAi)", "RGR(Di)", "RGR(Hi)",
    "RGR(Mi)", "RGR(Vi)")] <- "RGR"
  Table$growth[Table$RGR %in% c("GR(Di)", "GR(Hi)", "GR(Mi)", "GR(Shoot)")] <- "AbGR"

  Table$growth <- factor(Table$growth, levels = c("RGR", "AbGR"))
  Table$growth <- as.factor(Table$growth)

  Table["measurement"] <- "NA"
  Table$measurement[Table$RGR %in% c("RGR(Di)", "GR(Di)")] <- "Diameter"
  Table$measurement[Table$RGR %in% c("RGR(Hi)", "GR(Hi)")] <- "Height"
  Table$measurement[Table$RGR %in% c("RGR(Mi)", "GR(Mi)", "GR(Shoot)")] <- "Mass"
  Table$measurement[Table$RGR %in% c("RGR(?)", "RGR(Vi)", "RGR(CSAi)")] <- "Other"

  Table["growth.form"] <- "NA"
  Table$growth.form[Table$life.form %in% c("tree", "tree ", "NA")] <- "tree"
  Table$growth.form[Table$life.form %in% c("tree, shrub", "woody", "woody and semi-woody",
    "woody vines")] <- "woody"
  Table$growth.form[Table$life.form %in% c("liana and tree", "shrub and tree",
    "tree and herbs", "tree, herbs", "tree, shrub, herbs", "woody and herbs",
    "woody, forbs, herbs", "woody, shrub, perennial herbs")] <- "across growth form"

  Table$growth.form <- as.factor(Table$growth.form)
  # Nb of species used to performed the correlation
  names(Table)[names(Table) == "nb.sp.reported"] <- "nb.sp"

  i <- Table$nb.sp > Table$sample.size
  Table$nb.sp[i] <- Table$sample.size[i]

  i <- Table$bio.scale == "intrasp"
  Table$nb.sp[i] <- 1

  ## Clean the data trait stage
  Table["stageTrait"] <- NA
  Table$stageTrait <- Table$trait.stage
  Table$stageTrait[Table$trait.stage == "adult?"] <- "adult"
  Table$stageTrait[Table$trait.stage == "sapling?"] <- "sapling"
  Table$stageTrait[Table$trait.stage == "juvenile?"] <- "juvenile"
  Table$stageTrait[Table$trait.stage == "seedling?"] <- "seedling"
  Table$stageTrait[Table$trait.stage == "mix?"] <- "mix"
  Table$stageTrait[Table$trait.stage == "juveniladult"] <- "mix"

  ## RGR stage
  Table["stageRGR"] <- NA
  Table$stageRGR <- Table$RGR.stage
  Table$stageRGR[Table$RGR.stage == "adult?"] <- "adult"
  Table$stageRGR[Table$RGR.stage == "sapling?"] <- "sapling"
  Table$stageRGR[Table$RGR.stage == "juvenile?"] <- "juvenile"
  Table$stageRGR[Table$RGR.stage == "seedling?"] <- "seedling"
  Table$stageRGR[Table$RGR.stage == "mix?"] <- "mix"
  Table$stageRGR[Table$RGR.stage == "juveniladult"] <- "mix"


  # Table$stage[Table$stageRGR == "juvenile" & Table$ ] <- "seedling"
  Table$stage <- as.factor(Table$stage)
  Table$stage <- factor(Table$stage, levels = c("seedling", "juvenile", "sapling", "adult","mix"))
  Table$stage <- replace(Table$stage, Table$stage == "juvenile", "seedling")
  Table$stage <- factor(Table$stage, levels = c("seedling", "sapling", "adult"))


  ##### less restrictive, all experiments using data from database in parallel of
  ##### field experiment are considered as filed experiment or natural for map
  Table["experiment"] <- NA
  # levels(Table$experiment.type)
  Table$experiment[Table$experiment.type == "common garden"] <- "field"
  Table$experiment[Table$experiment.type == "field_forest experiment"] <- "field"
  Table$experiment[Table$experiment.type == "field_forest plantation"] <- "field"
  Table$experiment[Table$experiment.type == "shadehouse"] <- "field"

  Table$experiment[Table$experiment.type == "field_savanna"] <- "nature"
  Table$experiment[Table$experiment.type == "field_forest understory"] <- "nature"
  Table$experiment[Table$experiment.type == "field_gradient"] <- "nature"


  Table$experiment[Table$experiment.type == "glasshouse"] <- "control"
  Table$experiment[Table$experiment.type == "growth chamber"] <- "control"
  Table$experiment[Table$experiment.type == "hydroponic culture"] <- "control"

  Table$experiment[Table$experiment.type == "database"] <- "database"
  Table$experiment[Table$experiment.type == "mix.field.gchamber.database"] <- "database"  # j affecte ces etudes mixtes aux database pour montrer que ces donnes sont moins interessantes car basee sur d autres etudes peut etre deja dans la manip
  Table$experiment[Table$experiment.type == "mix.field.database"] <- "database"
  Table$experiment[Table$experiment.type == "mix.glasshouse.database"] <- "database"

  Table$experiment <- as.factor(Table$experiment)

  ##### more restrictive, all experiments using data from database are considered as
  ##### database data in model multiple
  Table["experiment.coord"] <- NA
  Table$experiment.coord[Table$experiment.type == "common garden"] <- "field"
  Table$experiment.coord[Table$experiment.type == "field_forest experiment"] <- "field"
  Table$experiment.coord[Table$experiment.type == "field_forest plantation"] <- "field"
  Table$experiment.coord[Table$experiment.type == "shadehouse"] <- "field"

  Table$experiment.coord[Table$experiment.type == "field_savanna"] <- "nature"
  Table$experiment.coord[Table$experiment.type == "field_forest understory"] <- "nature"
  Table$experiment.coord[Table$experiment.type == "field_gradient"] <- "nature"


  Table$experiment.coord[Table$experiment.type == "glasshouse"] <- "control"
  Table$experiment.coord[Table$experiment.type == "growth chamber"] <- "control"
  Table$experiment.coord[Table$experiment.type == "hydroponic culture"] <- "control"

  Table$experiment.coord[Table$experiment.type == "database"] <- "database"
  Table$experiment.coord[Table$experiment.type == "mix.field.gchamber.database"] <- "nature"  # j affecte ces etudes mixtes aux database pour montrer que ces donnes sont moins interessantes car basee sur d autres etudes peut etre deja dans la manip
  Table$experiment.coord[Table$experiment.type == "mix.field.database"] <- "nature"
  Table$experiment.coord[Table$experiment.type == "mix.glasshouse.database"] <- "control"

  Table$life.form <- as.factor(Table$life.form)

  Table$experiment.coord <- as.factor(Table$experiment.coord)

  ## Rename colomns
  Table <- rename(Table, c(stage.simi = "similarity"))

  Table <- subset(Table, select = c(id, idcor, authors, year, ref, doi, experiment, experiment.coord,
    stress, growth.form, veg.type, bio.scale, nb.sp, location, nb.site, names.s1,
    names.s2, names.s3, names.s4, names.s5, trace, lat.dms.s1, lat.dir.s1,
    long.dms.s1, long.dir.s1, lat.dms.s2, lat.dir.s2, long.dms.s2, long.dir.s2,
    lat.dms.s3, lat.dir.s3, long.dms.s3, long.dir.s3, lat.dms.s4, lat.dir.s4,
    long.dms.s4, long.dir.s4, lat.dms.s5, lat.dir.s5, long.dms.s5, long.dir.s5,
    RGR, growth, measurement, RGR.min, RGR.max, RGR.unit, stageRGR, similarity,
    stageTrait, stage, measure.size, size.min, size.max, size.mean.range, trait,
    trait.min, trait.max, sample.size, sign, coef, coef.type))

  Table <- subset(Table, Table$trait %in% c("SLA", "WD", "Amass",
    "Aarea", "Nmass", "Narea", "Hmax", "LA", "Pmass", "Seedmass", "LMR", "LAR",
    "NARmass", "NARarea", "Ks", "Vessel density", "Vessel size", "Thickness",
    "SA/LA", "Parea"))
}

standardise_data <- function(AllData) {

  RawData <- AllData

  ## correlation coefficient transformation in corr. r ####
  RawData["corr.r"] <- NA

  fun_Kendall <- function(x) (sin(0.5 * pi * x))
  fun_Spearman <- function(x) (2 * sin(pi * x/6))

  for (i in 1:length(RawData$coef.type)) {
    if (RawData$coef.type[i] == "Pearson") {
      RawData$corr.r[i] <- RawData$coef[i]
    }
    if (RawData$coef.type[i] == "CoefRegR") {
      RawData$corr.r[i] <- RawData$coef[i]
    }
    if (RawData$coef.type[i] == "CoefDetRsq") {
      RawData$corr.r[i] <- RawData$coef[i]
    }
    if (RawData$coef.type[i] == "CoefCorR") {
      RawData$corr.r[i] <- RawData$coef[i]
    }
    if (RawData$coef.type[i] == "no") {
      RawData$corr.r[i] <- "NA"
    }
    if (RawData$coef.type[i] == "slope") {
      RawData$corr.r[i] <- "NA"
    }
    if (RawData$coef.type[i] == "PartialCoefReg") {
      RawData$corr.r[i] <- "NA"
    }
    if (RawData$coef.type[i] == "Spearman" && RawData$sample.size[i] < 90) {
      RawData$corr.r[i] <- fun_Spearman(x = RawData$coef[i])
    }
    if (RawData$coef.type[i] == "Spearman" && RawData$sample.size[i] >= 90) {
      RawData$corr.r[i] <- RawData$coef[i]
    }
    if (RawData$coef.type[i] == "KendallRankCorr") {
      RawData$corr.r[i] <- fun_Kendall(x = RawData$coef[i])
    }
  }

  RawData$corr.r <- as.numeric(RawData$corr.r)

  #### Z-transformation: to normalise corr.r ####

  fn_z <- function(x) 0.5 * (log((1 + x)/(1 - x)))
  RawData$corr.z <- fn_z(RawData$corr.r)

  fn_vz <- function(x) (1/(x - 3))
  RawData$vr.z <- fn_vz(RawData$nb.sp)

  se.z <- NA
  fn_sez <- function(x) {
    ret <- rep(NA_real_, length(x))
    ii <- x >=0
    ret[ii] <- sqrt(x[ii])
    ret
  }
  se.z <- fn_sez(RawData$vr.z)

  fn_wiz <- function(x) 1/x
  RawData$wi.z <- fn_wiz(se.z)

  RawData
}

build_complete_data <- function(RawData) {
  #### Complete data set: focus on 5 traits:

  CompleteData <- subset(RawData, RawData$trait %in% c("SLA", "WD", "Aarea",
    "Hmax", "Seedmass"))
  CompleteData$trait <- factor(CompleteData$trait, levels = c("SLA", "WD", "Aarea",
    "Hmax", "Seedmass"))

  # select the column of interest
  CompleteData <- subset(CompleteData, select = c(id, idcor, authors, year, ref,
    doi, experiment, stress, growth.form, veg.type, bio.scale, nb.sp, RGR,
    growth, measurement, RGR.min, RGR.max, RGR.unit, stageRGR, similarity,
    stageTrait, stage, trait, trait.min, trait.max, measure.size, size.min,
    size.max, size.mean.range, coef, sample.size, corr.r, corr.z, vr.z, wi.z))
  CompleteData$id <- as.factor(CompleteData$id)
  CompleteData
}

Build_intersp_complete_data <- function(RawData) {
  #### Complete data set: focus on 5 traits:

  CompleteData <- subset(RawData, RawData$trait %in% c("SLA", "WD", "Aarea",
    "Hmax", "Seedmass") & RawData$bio.scale == "intersp")
  CompleteData$trait <- factor(CompleteData$trait, levels = c("SLA", "WD", "Aarea",
    "Hmax", "Seedmass"))

  # select the column of interest
  CompleteData <- subset(CompleteData, select = c(id, idcor, authors, year, ref,
    doi, experiment, stress, growth.form, veg.type, bio.scale, nb.sp, RGR,
    growth, measurement, RGR.min, RGR.max, RGR.unit, stageRGR, similarity,
    stageTrait, stage, trait, trait.min, trait.max, measure.size, size.min,
    size.max, size.mean.range, coef, sample.size, corr.r, corr.z, vr.z, wi.z))
  CompleteData$id <- as.factor(CompleteData$id)
  CompleteData
}

# Ideal data set: only the correlation coefficient measured at interspecific
# level, under unstressed conditions and for measurement of growth and trait
# perfomed at the same plant stage
build_ideal_data <- function(CompleteData) {

  IdealData <- subset(CompleteData, CompleteData$stress == "unstressed" & CompleteData$nb.sp > 10
     & CompleteData$similarity == "Sim")
  IdealData$stage <- factor(IdealData$stage, levels = c("seedling", "sapling",
    "adult"))
  IdealData$stageTrait <- factor(IdealData$stageTrait, levels = c("seedling",
    "juvenile", "sapling", "adult", "mix"))
  IdealData$stageRGR <- factor(IdealData$stageRGR, levels = c("seedling", "juvenile",
    "sapling", "adult", "mix"))
  levels(IdealData$RGR) <- factor(IdealData$RGR, levels = c("GR(Di)", "GR(Hi)",
    "GR(Mi)", "GR(Shoot)", "RGR(?)", "RGR(CSAi)", "RGR(Di)", "RGR(Hi)", "RGR(Mi)",
    "RGR(Vi)"))
  IdealData
}

EffectSizeSum <- function(data) {
  Fun_var <- function(y) {
    n <- length(y)
    z <- prod(y)
    var <- (1/n)^2 * (sum(y)) + 2 * 1 * sqrt(z)
    var
  }
  
  data2 <- data[, c("corr.z", "id", "stage","trait", "vr.z", "nb.sp")]

  res <- aggregate(corr.z ~ id + stage + trait, mean, data = na.omit(data2))
  res2 <- aggregate(vr.z ~ id + stage + trait, Fun_var, data = na.omit(data2))
  res3 <- aggregate(nb.sp ~ id + stage + trait, mean, data = na.omit(data2))
  table <- cbind(res, res2[, "vr.z"], res3[, "nb.sp"])
  colnames(table) <- c("id", "stage", "trait", "corr.z", "vr.z", "nb.sp")
  table
}

build_map_data <- function(RawData) {

  ######### File Coordinate ######## corresponds to the coordinate of the experiments
  ######### performed

  RawDataF <- subset(RawData, RawData$experiment.coord %in% c("field", "nature"))

  site1 <- subset(RawDataF, select = c("id", "ref", "doi", "experiment", "location",
    "names.s1", "lat.dms.s1", "lat.dir.s1", "long.dms.s1", "long.dir.s1"))
  site2 <- subset(RawDataF, select = c("id", "ref", "doi", "experiment", "location",
    "names.s2", "lat.dms.s2", "lat.dir.s2", "long.dms.s2", "long.dir.s2"))
  site3 <- subset(RawDataF, select = c("id", "ref", "doi", "experiment", "location",
    "names.s3", "lat.dms.s3", "lat.dir.s3", "long.dms.s3", "long.dir.s3"))
  site4 <- subset(RawDataF, select = c("id", "ref", "doi", "experiment", "location",
    "names.s4", "lat.dms.s4", "lat.dir.s4", "long.dms.s4", "long.dir.s4"))
  site5 <- subset(RawDataF, select = c("id", "ref", "doi", "experiment", "location",
    "names.s5", "lat.dms.s5", "lat.dir.s5", "long.dms.s5", "long.dir.s5"))

  site1 <- rename(site1, c(names.s1 = "name", lat.dms.s1 = "lat.dms", lat.dir.s1 = "lat.dir",
    long.dms.s1 = "long.dms", long.dir.s1 = "long.dir"))
  site2 <- rename(site2, c(names.s2 = "name", lat.dms.s2 = "lat.dms", lat.dir.s2 = "lat.dir",
    long.dms.s2 = "long.dms", long.dir.s2 = "long.dir"))
  site3 <- rename(site3, c(names.s3 = "name", lat.dms.s3 = "lat.dms", lat.dir.s3 = "lat.dir",
    long.dms.s3 = "long.dms", long.dir.s3 = "long.dir"))
  site4 <- rename(site4, c(names.s4 = "name", lat.dms.s4 = "lat.dms", lat.dir.s4 = "lat.dir",
    long.dms.s4 = "long.dms", long.dir.s4 = "long.dir"))
  site5 <- rename(site5, c(names.s5 = "name", lat.dms.s5 = "lat.dms", lat.dir.s5 = "lat.dir",
    long.dms.s5 = "long.dms", long.dir.s5 = "long.dir"))

  CoordTable <- rbind(site1, site2, site3, site4, site5)
  CoordTable <- unique(CoordTable)
  CoordTable <- na.omit(CoordTable)

  ## Transformation of coordinates from degree minute second to degree decimal
  ## Ensure Input non-negative and define Direction Prefix (south and west are
  ## with a minus sign)
  CoordTable["lat.dd"] <- "NA"
  z <- sapply(strsplit(as.character(CoordTable$lat.dms), "[:;/]"), as.character)
  for (i in 1:nrow(CoordTable)) {
    CoordTable$lat.dd[[i]] <- (as.numeric(z[1, i])) + (as.numeric(z[2, i])/60) +
      (as.numeric(z[3, i])/3600)
  }

  for (i in 1:nrow(CoordTable)) {
    if (CoordTable$lat.dir[[i]] %in% c("N", "NA")) {
      CoordTable$lat.dd[[i]] <- CoordTable$lat.dd[[i]]
    } else {
      CoordTable$lat.dd[[i]] <- (0 - (as.numeric(CoordTable$lat.dd[[i]])))
    }
  }

  CoordTable["long.dd"] <- "NA"
  z <- sapply(strsplit(as.character(CoordTable$long.dms), "[:;/]"), as.character)
  for (i in 1:nrow(CoordTable)) {
    CoordTable$long.dd[[i]] <- (as.numeric(z[1, i])) + (as.numeric(z[2, i])/60) +
      (as.numeric(z[3, i])/3600)
  }

  for (i in 1:nrow(CoordTable)) {
    if (CoordTable$long.dir[[i]] %in% c("E", "NA")) {
      CoordTable$long.dd[[i]] <- CoordTable$long.dd[[i]]
    } else {
      CoordTable$long.dd[[i]] <- (0 - (as.numeric(CoordTable$long.dd[[i]])))
    }
  }

  CoordTable$lat.dd <- as.numeric(CoordTable$lat.dd)
  CoordTable$long.dd <- as.numeric(CoordTable$long.dd)

  CoordTable
}


merge_bib_files <- function(meta, paper) {
  combined <- c(meta[setdiff(names(meta), names(paper))], paper)
  combined[[sort(names(combined))]]
}

snapshot_websci <- function(filename = "data/ref.traits/WebSci_all.csv"){
    Table <- read.csv(filename, sep = ";", dec = ",")
}
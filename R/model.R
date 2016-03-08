
# function for calculating effect size and intercept of the model with and
# without the explanatory variable : stage

fun_model <- function(x, y) {
  # models for x ('ideal' dataset) and y ('complete' datatset)
  null <- lmer(corr.z ~ 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # modele null cad sans la variable que je veux analyser
  m <- lmer(corr.z ~ stage + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # model avec la variable
  m1 <- lmer(corr.z ~ stage - 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # model avec les intercepts qui sont a 1

  null.s <- lmer(corr.z ~ 1 + (1 | id), data = y, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # idem mais pour le jeux y
  m.s <- lmer(corr.z ~ stage + (1 | id), data = y, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))
  m1.s <- lmer(corr.z ~ stage - 1 + (1 | id), data = y, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))
  summ.m <- summary(m)

  # extraction intercept of fixed parameters
  FP <- as.data.frame(fixef(m1))
  FP["stage"] <- "NA"
  FP["stage"] <- names(fixef(m1))

  FixSEdiag <- vcov(m1, useScale = FALSE)
  FixSE <- as.data.frame(sqrt(diag(FixSEdiag)))
  FixSE["stage"] <- "NA"
  FixSE["stage"] <- names(fixef(m1))

  FP.s <- as.data.frame(fixef(m1.s))
  FP.s["stage"] <- "NA"
  FP.s["stage"] <- names(fixef(m1.s))

  FixSEdiag.s <- vcov(m1.s, useScale = FALSE)
  FixSE.s <- as.data.frame(sqrt(diag(FixSEdiag.s)))
  FixSE.s["stage"] <- "NA"
  FixSE.s["stage"] <- names(fixef(m1.s))

  # new table with extracted data from models
  l <- (length(FP$stage) + length(FP.s$stage))
  l1 <- (length(FP$stage))
  l2 <- (length(FP.s$stage))
  a <- (l1 + 1)

  name <- list(1:l, c("stage", "stress", "Inte", "SE", "CIupper", "CIlower"))
  CoefModel <- as.data.frame(matrix(nrow = length(name[[1]]), ncol = length(name[[2]]),
    dimnames = name))
  CoefModel[1:l1, 1] <- c((FP$stage))
  CoefModel[a:l, 1] <- c((FP.s$stage))
  CoefModel[1:l1, 2] <- "ideal"
  CoefModel[a:l, 2] <- "complete"

  for (i in 1:l1) {
    CoefModel[i, 3] <- FP[i, 1]
    CoefModel[i, 4] <- FixSE[i, 1]
    CoefModel[i, 5] <- CoefModel[i, 4] * 1.96 + CoefModel[i, 3]
    CoefModel[i, 6] <- CoefModel[i, 4] * 1.96 - CoefModel[i, 3]
  }

  for (j in a:l) {
    CoefModel[j, 3] <- FP.s[j - l1, 1]
    CoefModel[j, 4] <- FixSE.s[j - l1, 1]
    CoefModel[j, 5] <- CoefModel[j, 4] * 1.96 + CoefModel[j, 3]
    CoefModel[j, 6] <- CoefModel[j, 4] * 1.96 - CoefModel[j, 3]
  }

  CoefModel$stage[CoefModel$stage == "stageseedling"] <- "seedling"
  CoefModel$stage[CoefModel$stage == "stagesapling"] <- "sapling"
  CoefModel$stage[CoefModel$stage == "stageadult"] <- "adult"


  fun1 <- function(z) length(na.omit(z$corr.z))
  fun2 <- function(z) length(na.omit(z$corr.z[z$stage == "seedling"]))
  fun5 <- function(z) length(na.omit(z$corr.z[z$stage == "sapling"]))
  fun3 <- function(z) length(na.omit(z$corr.z[z$stage == "adult"]))


  N <- data.frame(value = c(fun2(x), fun5(x), fun3(x), fun2(y), fun5(y), fun3(y)))

  CoefModel["N"] <- "NA"
  # CoefModel['N'] <- N[!(apply(N, 1, function(s) any(s == 0))),]
  CoefModel$N[CoefModel$stage == "seedling" & CoefModel$stress == "ideal"] <- fun2(x)
  CoefModel$N[CoefModel$stage == "seedling" & CoefModel$stress == "complete"] <- fun2(y)

  CoefModel$N[CoefModel$stage == "sapling" & CoefModel$stress == "ideal"] <- fun5(x)
  CoefModel$N[CoefModel$stage == "sapling" & CoefModel$stress == "complete"] <- fun5(y)

  CoefModel$N[CoefModel$stage == "adult" & CoefModel$stress == "ideal"] <- fun3(x)
  CoefModel$N[CoefModel$stage == "adult" & CoefModel$stress == "complete"] <- fun3(y)

  CoefModel
}

fun_model1 <- function(x, y) {
#   y <- data2
#   x <- data1
  # models for x ('ideal' dataset) and y ('complete' datatset)
  null <- lmer(corr.z ~ 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
               control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
                                     check.nobs.vs.nRE = "ignore"))  # modele null cad sans la variable que je veux analyser
  m <- lmer(corr.z ~ stage + (1 | id), data = x, weights = nb.sp, REML = TRUE,
            control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
                                  check.nobs.vs.nRE = "ignore"))  # model avec la variable
  m1 <- lmer(corr.z ~ stage - 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
             control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
                                   check.nobs.vs.nRE = "ignore"))  # model avec les intercepts qui sont a 1
  
  null.s <- lmer(corr.z ~ 1 + (1 | id), data = y, weights = nb.sp, REML = TRUE,
                 control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
                                       check.nobs.vs.nRE = "ignore"))  # idem mais pour le jeux y
  m.s <- lmer(corr.z ~ stage + (1 | id), data = y, weights = nb.sp, REML = TRUE,
              control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nRE = "ignore"))
  m1.s <- lmer(corr.z ~ stage - 1 + (1 | id), data = y, weights = nb.sp, REML = TRUE,
               control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
                                     check.nobs.vs.nRE = "ignore"))
  summ.m <- summary(m)
  
  # extraction intercept of fixed parameters
  FP <- as.data.frame(fixef(m1))
  FP["stage"] <- "NA"
  FP["stage"] <- names(fixef(m1))
  
  FixSEdiag <- vcov(m1, useScale = FALSE)
  FixSE <- as.data.frame(sqrt(diag(FixSEdiag)))
  FixSE["stage"] <- "NA"
  FixSE["stage"] <- names(fixef(m1))
  
  FP.s <- as.data.frame(fixef(m1.s))
  FP.s["stage"] <- "NA"
  FP.s["stage"] <- names(fixef(m1.s))
  
  FixSEdiag.s <- vcov(m1.s, useScale = FALSE)
  FixSE.s <- as.data.frame(sqrt(diag(FixSEdiag.s)))
  FixSE.s["stage"] <- "NA"
  FixSE.s["stage"] <- names(fixef(m1.s))
  
  # new table with extracted data from models
  l <- (length(FP$stage) + length(FP.s$stage))
  l1 <- (length(FP$stage))
  l2 <- (length(FP.s$stage))
  a <- (l1 + 1)
  
  name <- list(1:l, c("stage", "growth", "Inte", "SE", "CIupper", "CIlower"))
  CoefModel <- as.data.frame(matrix(nrow = length(name[[1]]), ncol = length(name[[2]]),
                                    dimnames = name))
  CoefModel[1:l1, 1] <- c((FP$stage))
  CoefModel[a:l, 1] <- c((FP.s$stage))
  CoefModel[1:l1, 2] <- "RGR"
  CoefModel[a:l, 2] <- "AGR"
  
  for (i in 1:l1) {
    CoefModel[i, 3] <- FP[i, 1]
    CoefModel[i, 4] <- FixSE[i, 1]
    CoefModel[i, 5] <- CoefModel[i, 3] + CoefModel[i, 4] * 1.96 
    CoefModel[i, 6] <- CoefModel[i, 3] - CoefModel[i, 4] * 1.96
  }
  
  for (j in a:l) {
    CoefModel[j, 3] <- FP.s[j - l1, 1]
    CoefModel[j, 4] <- FixSE.s[j - l1, 1]
    CoefModel[j, 5] <- CoefModel[j, 3] + CoefModel[j, 4] * 1.96 
    CoefModel[j, 6] <- CoefModel[j, 3] - CoefModel[j, 4] * 1.96
  }
  
  CoefModel$stage[CoefModel$stage == "stageseedling"] <- "seedling"
  CoefModel$stage[CoefModel$stage == "stagesapling"] <- "sapling"
  CoefModel$stage[CoefModel$stage == "stageadult"] <- "adult"
  
  
  fun1 <- function(z) length(na.omit(z$corr.z))
  fun2 <- function(z) length(na.omit(z$corr.z[z$stage == "seedling"]))
  fun5 <- function(z) length(na.omit(z$corr.z[z$stage == "sapling"]))
  fun3 <- function(z) length(na.omit(z$corr.z[z$stage == "adult"]))
  
  
  N <- data.frame(value = c(fun2(x), fun5(x), fun3(x), fun2(y), fun5(y), fun3(y)))
  
  CoefModel["N"] <- "NA"
  # CoefModel['N'] <- N[!(apply(N, 1, function(s) any(s == 0))),]
  CoefModel$N[CoefModel$stage == "seedling" & CoefModel$growth == "RGR"] <- fun2(x)
  CoefModel$N[CoefModel$stage == "seedling" & CoefModel$growth == "AGR"] <- fun2(y)
  
  CoefModel$N[CoefModel$stage == "sapling" & CoefModel$growth == "RGR"] <- fun5(x)
  CoefModel$N[CoefModel$stage == "sapling" & CoefModel$growth == "AGR"] <- fun5(y)
  
  CoefModel$N[CoefModel$stage == "adult" & CoefModel$growth == "RGR"] <- fun3(x)
  CoefModel$N[CoefModel$stage == "adult" & CoefModel$growth == "AGR"] <- fun3(y)
  
#   CoefModel$growth <- as.factor(CoefModel$growth)
#   CoefModel$growth <- factor(CoefModel$growth, levels = rev(levels(CoefModel$growth)))

  CoefModel
}

fun_model2 <- function(x) {
  # models for x ('ideal' dataset) and y ('complete' datatset)
  null <- lmer(corr.z ~ 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # modele null cad sans la variable que je veux analyser
  m <- lmer(corr.z ~ stage + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # model avec la variable
  m1 <- lmer(corr.z ~ stage - 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # model avec les intercepts qui sont a 1

  summ.m <- summary(m)

  # extraction intercept of fixed parameters
  FP <- as.data.frame(fixef(m1))
  FP["stage"] <- "NA"
  FP["stage"] <- names(fixef(m1))

  FixSEdiag <- vcov(m1, useScale = FALSE)
  FixSE <- as.data.frame(sqrt(diag(FixSEdiag)))
  FixSE["stage"] <- "NA"
  FixSE["stage"] <- names(fixef(m1))


  # new table with extracted data from models
  l <- (length(FP$stage))
  l1 <- (length(FP$stage))

  a <- (l1 + 1)

  name <- list(1:l, c("stage", "growth", "Inte", "SE", "CIupper", "CIlower"))
  CoefModel <- as.data.frame(matrix(nrow = length(name[[1]]), ncol = length(name[[2]]),
    dimnames = name))
  CoefModel[1:l1, 1] <- c((FP$stage))
  CoefModel[1:l1, 2] <- "all"


  for (i in 1:l1) {
    CoefModel[i, 3] <- FP[i, 1]
    CoefModel[i, 4] <- FixSE[i, 1]
    CoefModel[i, 5] <- CoefModel[i, 4] * 1.96 + CoefModel[i, 3]
    CoefModel[i, 6] <- CoefModel[i, 4] * 1.96 - CoefModel[i, 3]
  }

  CoefModel$stage[CoefModel$stage == "stageseedling"] <- "seedling"
  CoefModel$stage[CoefModel$stage == "stagesapling"] <- "sapling"
  CoefModel$stage[CoefModel$stage == "stageadult"] <- "adult"


  fun1 <- function(z) length(na.omit(z$corr.z))
  fun2 <- function(z) length(na.omit(z$corr.z[z$stage == "seedling"]))
  fun5 <- function(z) length(na.omit(z$corr.z[z$stage == "sapling"]))
  fun3 <- function(z) length(na.omit(z$corr.z[z$stage == "adult"]))


  N <- data.frame(value = c(fun2(x), fun5(x), fun3(x)))

  CoefModel["N"] <- "NA"
  # CoefModel['N'] <- N[!(apply(N, 1, function(s) any(s == 0))),]
  CoefModel$N[CoefModel$stage == "seedling" & CoefModel$growth == "all"] <- fun2(x)

  CoefModel$N[CoefModel$stage == "sapling" & CoefModel$growth == "all"] <- fun5(x)

  CoefModel$N[CoefModel$stage == "adult" & CoefModel$growth == "all"] <- fun3(x)

  CoefModel
}

fun_model3 <- function(x, y) {
  # models for x ('ideal' dataset) and y ('complete' datatset)
  null <- lmer(corr.z ~ 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # modele null cad sans la variable que je veux analyser
  m <- lmer(corr.z ~ stageRGR + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # model avec la variable
  m1 <- lmer(corr.z ~ stageRGR - 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # model avec les intercepts qui sont a 1

  null.s <- lmer(corr.z ~ 1 + (1 | id), data = y, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # idem mais pour le jeux y
  m.s <- lmer(corr.z ~ stageRGR + (1 | id), data = y, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))
  m1.s <- lmer(corr.z ~ stageRGR - 1 + (1 | id), data = y, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))

  summ.m <- summary(m)

  # extraction intercept of fixed parameters
  FP <- as.data.frame(fixef(m1))
  FP["stage"] <- "NA"
  FP["stage"] <- names(fixef(m1))

  FixSEdiag <- vcov(m1, useScale = FALSE)
  FixSE <- as.data.frame(sqrt(diag(FixSEdiag)))
  FixSE["stage"] <- "NA"
  FixSE["stage"] <- names(fixef(m1))

  FP.s <- as.data.frame(fixef(m1.s))
  FP.s["stage"] <- "NA"
  FP.s["stage"] <- names(fixef(m1.s))

  FixSEdiag.s <- vcov(m1.s, useScale = FALSE)
  FixSE.s <- as.data.frame(sqrt(diag(FixSEdiag.s)))
  FixSE.s["stage"] <- "NA"
  FixSE.s["stage"] <- names(fixef(m1.s))

  # new table with extracted data from models
  l <- (length(FP$stage) + length(FP.s$stage))
  l1 <- (length(FP$stage))
  l2 <- (length(FP.s$stage))
  a <- (l1 + 1)

  name <- list(1:l, c("stage", "stress", "Inte", "SE", "CIupper", "CIlower"))
  CoefModel <- as.data.frame(matrix(nrow = length(name[[1]]), ncol = length(name[[2]]),
    dimnames = name))
  CoefModel[1:l1, 1] <- c((FP$stage))
  CoefModel[a:l, 1] <- c((FP.s$stage))
  CoefModel[1:l1, 2] <- "ideal"
  CoefModel[a:l, 2] <- "complete"

  for (i in 1:l1) {
    CoefModel[i, 3] <- FP[i, 1]
    CoefModel[i, 4] <- FixSE[i, 1]
    CoefModel[i, 5] <- CoefModel[i, 4] * 1.96 + CoefModel[i, 3]
    CoefModel[i, 6] <- CoefModel[i, 4] * 1.96 - CoefModel[i, 3]
  }

  for (j in a:l) {
    CoefModel[j, 3] <- FP.s[j - l1, 1]
    CoefModel[j, 4] <- FixSE.s[j - l1, 1]
    CoefModel[j, 5] <- CoefModel[j, 4] * 1.96 + CoefModel[j, 3]
    CoefModel[j, 6] <- CoefModel[j, 4] * 1.96 - CoefModel[j, 3]
  }

  CoefModel$stage[CoefModel$stage == "stageRGRseedling"] <- "seedling"
  CoefModel$stage[CoefModel$stage == "stageRGRjuvenile"] <- "juvenile"
  CoefModel$stage[CoefModel$stage == "stageRGRsapling"] <- "sapling"
  CoefModel$stage[CoefModel$stage == "stageRGRadult"] <- "adult"
  CoefModel$stage[CoefModel$stage == "stageRGRmix"] <- "mix"

  fun1 <- function(z) length(na.omit(z$corr.z))
  fun2 <- function(z) length(na.omit(z$corr.z[z$stageRGR == "juvenile"]))
  fun5 <- function(z) length(na.omit(z$corr.z[z$stageRGR == "sapling"]))
  fun3 <- function(z) length(na.omit(z$corr.z[z$stageRGR == "adult"]))
  fun4 <- function(z) length(na.omit(z$corr.z[z$stageRGR == "seedling"]))
  fun6 <- function(z) length(na.omit(z$corr.z[z$stageRGR == "mix"]))

  N <- data.frame(value = c(fun2(x), fun5(x), fun3(x), fun2(y), fun5(y), fun3(y)))

  CoefModel["N"] <- "NA"
  # CoefModel['N'] <- N[!(apply(N, 1, function(s) any(s == 0))),]
  CoefModel$N[CoefModel$stage == "seedling" & CoefModel$stress == "ideal"] <- fun4(x)
  CoefModel$N[CoefModel$stage == "seedling" & CoefModel$stress == "complete"] <- fun4(y)


  CoefModel$N[CoefModel$stage == "juvenile" & CoefModel$stress == "ideal"] <- fun2(x)
  CoefModel$N[CoefModel$stage == "juvenile" & CoefModel$stress == "complete"] <- fun2(y)

  CoefModel$N[CoefModel$stage == "sapling" & CoefModel$stress == "ideal"] <- fun5(x)
  CoefModel$N[CoefModel$stage == "sapling" & CoefModel$stress == "complete"] <- fun5(y)

  CoefModel$N[CoefModel$stage == "adult" & CoefModel$stress == "ideal"] <- fun3(x)
  CoefModel$N[CoefModel$stage == "adult" & CoefModel$stress == "complete"] <- fun3(y)

  CoefModel$N[CoefModel$stage == "mix" & CoefModel$stress == "ideal"] <- fun6(x)
  CoefModel$N[CoefModel$stage == "mix" & CoefModel$stress == "complete"] <- fun6(y)


  CoefModel
}

fun_model_multiple3 <- function(x) {
  # x data en opt
  null <- lmer(corr.z ~ 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
        tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
        tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))  # modele null cad sans la variable que je veux analyser
  m1 <- lmer(corr.z ~ stageRGR - 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
        tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
        tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))  # model avec les intercepts qui sont a 1
  gr1 <- lmer(corr.z ~ RGR - 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
        tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
        tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))
  veg1 <- lmer(corr.z ~ growth.form - 1 + (1 | id), data = x, weights = nb.sp,
    REML = TRUE, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
        tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
        tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))
  exp1 <- lmer(corr.z ~ experiment - 1 + (1 | id), data = x, weights = nb.sp,
    REML = TRUE, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
        tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
        tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))

  a <- list(m1, gr1, veg1, exp1)

  ## extraction des valeurs du modeles
  fun_FE <- function(z) {
    fixef(z)
  }  #fixef permet de extraire les intercepts des effet fixes des modeles
  fun_name <- function(z) {
    names(fun_FE(z))
  }  #fixef permet de extraire les SE des modeles
  fun_FixFE <- function(z) {
    sqrt(diag(vcov(z, useScale = FALSE)))
  }
  l <- (length(unlist(lapply(a, fun_FE))))

  ## creation du tableau de donnes
  name <- list(1:l, c("model", "params", "mean", "SE", "upper", "lower"))
  dat <- as.data.frame(matrix(nrow = length(name[[1]]), ncol = length(name[[2]]),
    dimnames = name))

  dat["mean"] <- data.frame(cbind(unlist(lapply(a, fun_FE))))
  dat["params"] <- data.frame(cbind(unlist(lapply(a, fun_name))))
  dat["SE"] <- as.data.frame(cbind(unlist(lapply(a, fun_FixFE))))
  dat["upper"] <- dat$mean + (dat$SE * 1.96)
  dat["lower"] <- dat$mean - (dat$SE * 1.96)

  s <- length(fun_FE(m1))
  g <- length(fun_FE(gr1))
  v <- length(fun_FE(veg1))
  e <- length(fun_FE(exp1))

  dat[1:s, 1] <- "M1"
  dat[(s + 1):(s + g), 1] <- "G1"
  dat[(s + g + 1):(s + g + v), 1] <- "V1"
  dat[(s + g + v + 1):l, 1] <- "E1"

  dat
}

fun_model_multiple3.1 <- function(x) {
  # x data en opt
  null <- lmer(corr.z ~ 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
        tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
        tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))  # modele null cad sans la variable que je veux analyser
  m1 <- lmer(corr.z ~ stageRGR - 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
        tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
        tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))  # model avec les intercepts qui sont a 1
  gr1 <- lmer(corr.z ~ RGR - 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
        tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
        tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))
  exp1 <- lmer(corr.z ~ experiment - 1 + (1 | id), data = x, weights = nb.sp,
    REML = TRUE, control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
        tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
        tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))

  a <- list(m1, gr1, exp1)

  ## extraction des valeurs du modeles
  fun_FE <- function(z) {
    fixef(z)
  }  #fixef permet de extraire les intercepts des effet fixes des modeles
  fun_name <- function(z) {
    names(fun_FE(z))
  }  #fixef permet de extraire les SE des modeles
  fun_FixFE <- function(z) {
    sqrt(diag(vcov(z, useScale = FALSE)))
  }
  l <- (length(unlist(lapply(a, fun_FE))))

  ## creation du tableau de donnes
  name <- list(1:l, c("model", "params", "mean", "SE", "upper", "lower"))
  dat <- as.data.frame(matrix(nrow = length(name[[1]]), ncol = length(name[[2]]),
    dimnames = name))

  dat["mean"] <- data.frame(cbind(unlist(lapply(a, fun_FE))))
  dat["params"] <- data.frame(cbind(unlist(lapply(a, fun_name))))
  dat["SE"] <- as.data.frame(cbind(unlist(lapply(a, fun_FixFE))))
  dat["upper"] <- dat$mean + (dat$SE * 1.96)
  dat["lower"] <- dat$mean - (dat$SE * 1.96)

  s <- length(fun_FE(m1))
  g <- length(fun_FE(gr1))
  e <- length(fun_FE(exp1))

  dat[1:s, 1] <- "M1"
  dat[(s + 1):(s + g), 1] <- "G1"
  dat[(s + g + 1):l, 1] <- "E1"

  dat
}


fun_model_growth <- function(x) {
  x[["stagegrowth"]] <- NA
  x$stagegrowth <- do.call(paste, c(x[c("stage", "growth")], sep = "")) 

  # models for x ('ideal' dataset) and y ('complete' datatset)
  null <- lmer(corr.z ~ stage + (1 | id), data = x, weights = nb.sp, REML = FALSE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # modele null cad sans la variable que je veux analyser
  m <- lmer(corr.z ~ stagegrowth + (1 | id), data = x, weights = nb.sp, REML = FALSE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # model avec la variable
  m1 <- lmer(corr.z ~ stagegrowth  - 1 + (1 | id), data = x, weights = nb.sp, REML = TRUE,
    control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
      check.nobs.vs.nRE = "ignore"))  # model avec les intercepts qui sont a 1
  
  summ.m <- summary(m)
  
  # extraction intercept of fixed parameters
  FP <- as.data.frame(fixef(m1))
  FP["stage"] <- "NA"
  FP["stage"] <- names(fixef(m1))
  
  FixSEdiag <- vcov(m1, useScale = FALSE)
  FixSE <- as.data.frame(sqrt(diag(FixSEdiag)))
  FixSE["stage"] <- "NA"
  FixSE["stage"] <- names(fixef(m1))
  
  
  # new table with extracted data from models
  l <- (length(FP$stage))
  l1 <- (length(FP$stage))
  
  a <- (l1 + 1)
  
  name <- list(1:l, c("comp","stage", "growth", "Inte", "SE", "CIupper", "CIlower"))
  CoefModel <- as.data.frame(matrix(nrow = length(name[[1]]), ncol = length(name[[2]]),
    dimnames = name))
  CoefModel[1:l1, 1] <- c((FP$stage))
  CoefModel[1:l1, 2] <- c((FP$stage))
  
  
  for (i in 1:l1) {
    CoefModel[i, 4] <- FP[i, 1]
    CoefModel[i, 5] <- FixSE[i, 1]
    CoefModel[i, 6] <- CoefModel[i, 4] + CoefModel[i, 5] * 1.96 
    CoefModel[i, 7] <- CoefModel[i, 4] - CoefModel[i, 5] * 1.96
  }
  
  CoefModel$stage[CoefModel$comp %in% c("stagegrowthseedlingAbGR", "stagegrowthseedlingRGR")] <- "seedling"
  CoefModel$stage[CoefModel$comp %in% c("stagegrowthsaplingAbGR", "stagegrowthsaplingRGR")] <- "sapling"
  CoefModel$stage[CoefModel$comp %in% c("stagegrowthadultAbGR", "stagegrowthadultRGR")] <- "adult"
  
  CoefModel$growth[CoefModel$comp %in% c("stagegrowthseedlingAbGR","stagegrowthsaplingAbGR","stagegrowthadultAbGR")] <- "AGR"
  CoefModel$growth[CoefModel$comp %in% c("stagegrowthsaplingRGR", "stagegrowthadultRGR","stagegrowthseedlingRGR")] <- "RGR"
  
  
  fun1 <- function(z) length(na.omit(z$corr.z))
  fun2 <- function(z) length(na.omit(z$corr.z[z$stagegrowth == "seedlingAbGR"]))
  fun5 <- function(z) length(na.omit(z$corr.z[z$stagegrowth == "saplingAbGR"]))
  fun3 <- function(z) length(na.omit(z$corr.z[z$stagegrowth == "adultAbGR"]))
  fun4 <- function(z) length(na.omit(z$corr.z[z$stagegrowth == "seedlingRGR"]))
  fun6 <- function(z) length(na.omit(z$corr.z[z$stagegrowth == "saplingRGR"]))
  fun7 <- function(z) length(na.omit(z$corr.z[z$stagegrowth == "adultRGR"]))
  
  N <- data.frame(value = c(fun2(x), fun5(x), fun3(x),fun4(x), fun6(x), fun7(x)))
  
  CoefModel["N"] <- "NA"
  # CoefModel['N'] <- N[!(apply(N, 1, function(s) any(s == 0))),]
  CoefModel$N[CoefModel$stage == "seedling" & CoefModel$growth == "AGR"] <- fun2(x)
  
  CoefModel$N[CoefModel$stage == "sapling" & CoefModel$growth == "AGR"] <- fun5(x)
  
  CoefModel$N[CoefModel$stage == "adult" & CoefModel$growth == "AGR"] <- fun3(x)
  
  CoefModel$N[CoefModel$stage == "seedling" & CoefModel$growth == "RGR"] <- fun4(x)
  
  CoefModel$N[CoefModel$stage == "sapling" & CoefModel$growth == "RGR"] <- fun6(x)
  
  CoefModel$N[CoefModel$stage == "adult" & CoefModel$growth == "RGR"] <- fun7(x)
  
  model.table <- list(CoefModel=CoefModel, LRT=anova(null, m)$Chisq[2], PVAL=anova(null, m)$Pr[2])
  model.table
}

# Weighted correlation for FigA4
table_overall.stage <- function(y) {
  name <- list(1:3, c("stage", "corr.r", "SD", "CI"))
  MeanTab <- as.data.frame(matrix(nrow = length(name[[1]]), ncol = length(name[[2]]),
    dimnames = name))
  y <- y[!is.na(y[, "corr.r"]), ]
  MeanTab$stage <- c("seedling", "sapling", "adult")
  MeanTab$corr.r <- sapply(split(y, na.omit(y$stage)), function(x) wtd.mean(x$corr.r,
    x$nb.sp))
  MeanTab$SD <- sapply(split(y, na.omit(y$stage)), function(x) wtd.var(x$corr.r,
    x$nb.sp))
  MeanTab$CI <- sapply(split(y, na.omit(y$stage)), function(x) {
    1.96 * (sqrt(wtd.var(x$corr.r, x$nb.sp))/sqrt(length(x)))
  })

  n <- count(y, "stage")

  MeanTab <- merge(MeanTab, n, by = "stage")
  MeanTab$freq <- as.character(MeanTab$freq)

  as.data.frame(MeanTab)

}

table_overall <- function(y) {
  name <- list(1:1, c("stage", "corr.r", "CI", "SD", "freq"))
  y <- y[!is.na(y[, "corr.r"]), ]
  MeanTab <- as.data.frame(matrix(nrow = length(name[[1]]), ncol = length(name[[2]]),
    dimnames = name))
  MeanTab[, 1] <- "total"
  MeanTab$corr.r <- wtd.mean(y$corr.r, y$nb.sp)
  MeanTab$SD <- sqrt(wtd.var(y$corr.r, y$nb.sp))
  MeanTab$CI <- 1.96 * (sqrt(wtd.var(y$corr.r, y$nb.sp))/sqrt(length(y)))
  n <- count(y, "stage")
  MeanTab$freq <- sum(n$freq)
  MeanTab$freq <- as.character(MeanTab$freq)
  as.data.frame(MeanTab)

}



# Fit and summarise model
fit_lmer <- function(data, effect) {
  default_lmer_control <- function() {
    # Ignore warnings when we have fewer observations than levels
    lmerControl(check.nobs.vs.nlev = "ignore",
                check.nobs.vs.rankZ = "ignore",
                check.nobs.vs.nRE = "ignore")
  }
  # Throws messages when there are less than 3 observations for each stage,
  # but we're fine with that
  suppressWarnings({
    null <- lmer(corr.z ~ 1 + (1 | id), data = data, weights = nb.sp, REML = FALSE,
                 control = default_lmer_control())
    m <- lmer(as.formula(paste("corr.z ~", effect, "+ (1 | id)")), data = data, weights = nb.sp, REML = FALSE,
              control = default_lmer_control())
    ## log likelihood ratio test, anova makes sure to use the ML criterion, even
    ## though the models were fit with REML
    av <- anova(null, m)
    list(LRT = av$Chisq[2], PVAL = av$Pr[2])
  })
}

fun_OneLR <- function(x) {
  fit_lmer(x, "stage")[["LRT"]]
}

fun_Onepvalue <- function(x) {
  fit_lmer(x, "stage")[["PVAL"]]
}

fun_OneLR_year <- function(x) {
  fit_lmer(x, "year")[["LRT"]]
}

fun_Onepvalue_year <- function(x) {
  fit_lmer(x, "year")[["PVAL"]]
}


fun_Heterogeneity.CI <- function(x) {
  x <- droplevels(x)
  
  res_null <- rma(corr.z, vr.z, 
                  data = x[!is.na(x$corr.z) & !is.na(x$nb.sp),])
  confint(res_null)
}

fun_Heterogeneity.H <- function(x, mods) {
  
  x <- droplevels(x)
  res_null <- rma(corr.z, vr.z, 
                  data = x[!is.na(x$corr.z) & !is.na(x$nb.sp),])
  res_stage <- rma(corr.z, vr.z, mods = mods, 
                  data = x[!is.na(x$corr.z) & !is.na(x$vr.z),])
  Hmodel <- (res_null$tau2 - res_stage$tau2)/res_null$tau2 * 100  #% of the total amount of heterogeneity can be accounted for by including the moderators of the model
  H <- res_null$I2 
  l <- list(Hlevel=H, Hexplained = Hmodel, pvalue = res_stage$QEp, Qdf = res_stage$QE,
    effectsize_moderator = res_stage$zval, pvalueeffectsize = res_stage$pval)
  l
}


fun_trim.and.fill_number <- function(x) {
  x <- droplevels(x)
  
  res <- rma(corr.z, vr.z, data = x[!is.na(x$corr.z) & !is.na(x$vr.z), ])
  ### carry out trim-and-fill analysis
  taf <- trimfill(res)
  l <- list(N= taf$k0, estimate= taf$b, se =taf$se, CImax= taf$ci.ub, CImin=taf$ci.lb)
  l
}


fun_fsn <- function(x){
  a <- fsn(corr.z, vr.z, data = x[!is.na(x$corr.z) & !is.na(x$vr.z), ], type = "Rosenthal", alpha = 0.05, digits = 4)
  a
}


fun_year <- function(a) {
  table <- a

  res <- rma(corr.z, vr.z, mod = ~year, 
                  data = table[!is.na(table$corr.z) & !is.na(table$nb.sp),])

  x <- subset(table, table$stage == "seedling")
  res_juv <- rma(corr.z, vr.z, mod = ~year, 
                  data = x[!is.na(x$corr.z) & !is.na(x$nb.sp),])

  y <- subset(table, table$stage == "sapling")
  res_sap <- rma(corr.z, vr.z, mod = ~year, 
                  data = y[!is.na(y$corr.z) & !is.na(y$nb.sp),])

  z <- subset(table, table$stage == "adult")
  browser()
  res_adu <- rma(corr.z, vr.z, mod = ~year, 
                  data = z[!is.na(z$corr.z) & !is.na(z$nb.sp),])

  list(all = res, seedling = res_juv, sapling = res_sap, adult = res_adu)
}


fun_snapshot <- function(Snapshot){
  
  Snapshot$ref <- with(Snapshot, paste(AU, TI,PY, sep = "_"))
    nb_ref <- length(unique(Snapshot$ref))
  
  x_2 <- subset(Snapshot, Snapshot$trait %in% c("SLA", "WD"))
   nb_ref_uniq_2  <- length(unique(x_2$ref))
  
  x_5 <- subset(Snapshot, Snapshot$trait %in% c("SLA", "WD", "Aarea", "Hmax", "Seedmass"))
   nb_ref_uniq_5 <- length(unique(x_5$ref))
 
  perc_19traits = (nb_ref_uniq_2/nb_ref) *100
  perc_5traits = (nb_ref_uniq_2/nb_ref_uniq_5) *100
  
  list(main.traits =perc_5traits, all.traits =perc_19traits)
}


fun_gap <- function(CleanData){
  library(plyr)
  
  table <- ddply(CleanData, c("trait", "stage"), function(x) length(x$corr.r))
  table2 <- ddply(table, c("trait"), function(x) max(x$V1)/ sum(x$V1))
  
  table3 <- subset(table2, table2$V1 > 0.75)
  length(table3$V1)
} 

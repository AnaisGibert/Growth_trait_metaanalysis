## Figure main text
figure_1 <- function(RawData) {

  p2 <- my_plot_1("", ggplot(RawData, aes(x = reorder(factor(trait), factor(trait),
    function(x) length(x) * 1), fill = stage, order = stage)))
  p2
}

figure_2 <- function(CompleteData_inter) {

  s <- unique(CompleteData_inter[c("ref", "id", "stageRGR", "measure.size", "size.min",
    "size.max", "growth.form", "stage")])
  s <- s[order(s$size.max, s$size.min, s$stageRGR), ]
  s <- na.omit(s)

  s$ref <- gsub(".", " ", as.character(s$ref), fixed=TRUE)
  s$ref <- gsub(" and ", " & ", s$ref, fixed=TRUE)

  # sort by trait et stage
  s$stageRGR <- factor(s$stageRGR, levels = c("seedling", "juvenile", "sapling",
    "adult", "mix"))
  s$stage <- factor(s$stage, levels = c("seedling", "juvenile", "sapling", "adult",
    "mix"))
  s$stage <- replace(s$stage, s$stage == "juvenile", "seedling")

  s$id <- as.factor(s$id)
  sh <- subset(s, s$measure.size == "height")
  sa <- subset(s, s$measure.size == "age")
  sd <- subset(s, s$measure.size == "diameter")
  sd <- sd[order(sd$size.min, sd$size.max), ]

  heights <- c(nrow(sa), nrow(sh), nrow(sd)) + 10
  heights <- heights/sum(heights)

  layout(matrix(c(1, 2, 3), 3, 1, byrow = TRUE), widths = c(3, 3, 3), heights = heights)
  par(mar = c(5, 8, 0.1, 0.5), oma = c(0, 2, 3, 1))

  plotCI2 <- function(data, cuts, xlab, lab) {

    rescale <- function(x, cuts) {
      x[x > max(cuts)] <- max(cuts)
      x2 <- x
      for (i in seq_along(cuts)[-1]) {
        low <- cuts[i - 1]
        high <- cuts[i]
        ii <- (low < x & x <= high)
        x2[ii] <- i - 2 + (x[ii] - low)/(high - low)
      }
      x2
    }

    data$size.min <- rescale(data$size.min, cuts)
    data$size.max <- rescale(data$size.max, cuts)

    data <- data[order(data$size.max, data$size.min), ]
    data["growthf"] <- -0.02

    n <- nrow(data)
    y <- rev(seq_len(n))
    cols <- c("#13519E", "#73005C", "#F57C34", "#E6224C", "grey")[data$stageRGR]
    plot(NA, xlim = c(0, 3), ylim = c(0, n+1), yaxt = "n", xaxt = "n", xlab = "",
      ylab = "", yaxs="i")
    mtext(xlab, 1, line = 2, cex = 0.75)
    segments(data$size.min, y, data$size.max, col = cols)
    i <- data$size.max == data$size.min
    points(data$size.min[i], y[i], col = cols[i], pch = "-")

    j <- data$stageRGR != data$stage
    points(data$size.max[j] + 0.05, y[j], col = "black", pch = "*", cex=1.5)

    axis(1, at = 0:3, labels = cuts, las = 1)
    axis(2, at = y, labels = data$ref, las = 1, cex.axis = 0.6)
    abline(v = 1, col = "black", lty = 1, lwd = 0.5)
    abline(v = 2, col = "black", lty = 1, lwd = 0.5)

    text(-0.5, n+2, lab, xpd=NA)
  }

  plotCI2(sa, c(0, 1, 5, 10), "Age (yrs)", "a)")
  legend(2.3, 42, c("seedling", "juvenile", "sapling", "adult", "mix"), lwd = 1,  title = "Original stage:",
    col = c("#13519E", "#73005C", "#F57C34", "#E6224C", "grey"), pch = NA, bty = "n", cex = 0.75,
   x.intersp = 0.3, y.intersp = 1, seg.len = 0.5,  title.adj = 0)
  legend(2.3, 32, c("reassigned"), lwd = NA, col = "black",
    pch = c( "*"), title = "Flag:", bty = "n", cex = 0.75,
    x.intersp = 0.3, y.intersp = 1, seg.len = 0.5,  title.adj = 0)

  text(0.5, 45, "Seedlings", xpd=NA, col = "black")
  text(1.5, 45, "Saplings", xpd=NA, col = "black")
  text(2.5, 45, "Adults", xpd=NA, col = "black")

  plotCI2(sh, c(0, 0.5, 2, 20), "Height (m)", "b)")

  plotCI2(sd, c(0, 1, 10, 80), "Diameter (cm)", "c)")

}

figure_panels_traits_model <- function(fits, ...) {

  coeff.plot <- function(x, title = "",
    limit.x.min = -1, limit.x.max = 1.5,
    vjust.value = 2, width = 0.4,
    colors = c("#9E5F3A", "black"),
    category_variable = "growth",
    labels = c("RGR, conservative dataset", "RGR, entire dataset"),
    breaks = c("RGR", "AGR")) {

    data <- x[["model"]]
    LRT <- x[["LRT"]]
    PVAL <- x[["PVAL"]]

    data[["by"]] <- as.factor(data[[category_variable]])

    label.n.x <- limit.x.min + 0.85*(limit.x.max - limit.x.min)

    dodge <- position_dodge(width = width)

    significativite <- ""
    round.value <- 2
    if(!is.na(PVAL)) {
      significativite <- "ns"
      if(PVAL <= 0.05) {
        significativite <- "*"
      } else if(PVAL <= 0.01) {
        significativite <- "***"
        round.value <- 4
      }
    }

    p <- ggplot(data, aes(stage, Inte, ymin = Inte - 1.96 * SE, ymax = Inte +1.96 * SE, color = by)) +
      labs(title = title) + xlab("Effect size") + ylab("Effect size") +
      geom_vline(xintercept = 0,color = "white") + geom_hline(yintercept = 0, size = 0.3, linetype = "dashed") +
      scale_x_discrete("Stage", limit = c("seedling", "sapling", "adult")) +
      coord_flip() + mytheme() + theme(legend.position = "none") +
      geom_point(aes(x = stage, y = Inte), position = dodge) +
      geom_errorbar(aes(x = stage, y = Inte), size = 0.4, width = 0.2, position = dodge) +
      scale_y_continuous("Effect size",limits = c(limit.x.min, limit.x.max)) +
      geom_text(aes(stage, label.n.x, label = paste("n =", N), color = by),size = 2, data = data, parse = FALSE, position = dodge, hjust = 0) +
      scale_color_manual(name = "", values = colors, breaks = breaks, labels = labels) +
      annotate("text", x = 0.6, y = limit.x.max, size = 1.75, hjust = 1,
        label = paste0("LRT: ", round(LRT, 2),", p: ", round(PVAL, round.value)," (", significativite, ")"))
    p
  }

  p1 <- coeff.plot(fits[["SLA"]], title = "a) SLA",
                limit.x.min = -1, limit.x.max = 1.5,  ...) +
        theme(plot.margin  = unit(c(2.5, 2, 1.5, 0), "mm"), axis.title.x = element_blank(),
              axis.title.y = element_text(colour = "white"))

  p2 <- coeff.plot(fits[["WD"]], title = "b) Wood density",
               limit.x.min = -1, limit.x.max = 1,  ...) +
        theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
              plot.margin  = unit(c(2.5, 0, 1.5, 1), "mm"), axis.title.x = element_blank())

  p3 <- coeff.plot(fits[["Hmax"]], title = "c) Hmax",
                limit.x.min = -1, limit.x.max = 1.5, ...) +
        theme(plot.margin  = unit(c(2.5, 2, 1.5, 0), "mm"), axis.title.x = element_blank())

  p4 <- coeff.plot(fits[["Seedmass"]], title = "d) Seed mass",
                 limit.x.min = -1.5, limit.x.max = 1.5, ...) +
        theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
              plot.margin  = unit(c(2.5, 0, -3.5, 1), "mm"))

  p5 <- coeff.plot(fits[["Aarea"]], title = "e) Aarea",
                limit.x.min = -2.4, limit.x.max =3.5, ...)  +
        theme(legend.title = element_blank(), legend.justification = c(0, 0), legend.position = c(1, 0.3),
          legend.key = element_blank(), plot.margin  = unit(c(2.5, 2, 1.5, 0), "mm"),
          axis.title.y = element_text(colour = "white"))

  grid.arrange(p1, p2, p3, p4, p5, ncol = 2, nrow = 3, widths = c(1.3, 1),heights = c(1, 1, 1.1))
}

figure_3 <- function(GIrgr , GCrgr) {

  fit_model <- function(trait, GIrgr , GCrgr) {
    ret <- list()
    ret[["model"]] <- fun_model1(GIrgr[[trait]],GCrgr[[trait]])
    model <- fit_lmer(GIrgr[[trait]], "stage")
    ret[["LRT"]] <- model[["LRT"]]
    ret[["PVAL"]] <- model[["PVAL"]]
    ret
  }

  traits <- c("SLA", "WD", "Hmax", "Seedmass", "Aarea")
  fits <- lapply(traits, fit_model, GIrgr , GCrgr)
  names(fits) <- traits

  figure_panels_traits_model(fits)
}


figure_4 <- function(GCi) {

  fit_model <- function(trait, GCi) {
    ret <- list()
    fit <- fun_model_growth(GCi[[trait]])
    ret[["model"]] <- fit[["CoefModel"]]
    ret[["LRT"]] <- fit[["LRT"]]
    ret[["PVAL"]] <- fit[["PVAL"]]
    ret
  }

  traits <- c("SLA", "WD", "Hmax", "Seedmass", "Aarea")

  fits <- lapply(traits, fit_model, GCi)
  names(fits) <- traits

  figure_panels_traits_model(fits, labels = c("RGR", "AGR"))
}

download_baad <- function(destination_filename) {
  url <- "https://github.com/dfalster/baad/releases/download/v1.0.0/baad.rds"
  download(url, destination_filename, mode = "wb")
}

## Figure appendix
figure_A1 <- function(baad) {

  single_plot <- function(px, py, data) {

    plot(data[[px$var]], data[[py$var]], log = "xy", xlim = px$lim, ylim = py$lim,
      xlab = px$lab, ylab = py$lab, col = "grey", pch = 16)
    abline(v = px$sap, col = "#13519E", lty = 5)
    abline(h = py$sap, col = "#13519E", lty = 5)
    abline(v = px$adult, col = "#E6224C", lty = 5)
    abline(h = py$adult, col = "#E6224C", lty = 5)
  }

  data <- baad[["data"]]

  pars <- list()
  pars[["dia"]] <- list(var = "d.ba", lab = "Basal diameter (m)", lim = c(1e-04,
    1), sap = 0.0025, adult = 0.1)

  pars[["ht"]] <- list(var = "h.t", lab = "Height (m)", lim = c(0.01, 100), sap = 0.5,
    adult = 5)

  pars[["age"]] <- list(var = "age", lab = "Age (yr)", lim = c(0.1, 200), sap = 1,
    adult = 30)


  par(mfrow = c(1, 3))
  single_plot(pars[["dia"]], pars[["age"]], data)
  single_plot(pars[["ht"]], pars[["age"]], data)
  single_plot(pars[["dia"]], pars[["ht"]], data)
}

figure_A2 <- function(CoordTable) {

  mapWorld <- borders("world", colour = "#FFCC33", fill = "#FFCC33", lty = 0)  # create a layer of borders
  mp <- ggplot() + mapWorld + theme(text = element_text(size = 9, colour = "black"),
    title = element_text(size = 9), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "#333333", fill = "#6699FF"))

  x <- CoordTable$lat.dd
  y <- CoordTable$long.dd
  coordinate.map <- mp + geom_point(aes(x = y, y = x), color = "red", alpha = I(7/10))
  coordinate.map
}

figure_A3 <- function(GC) {
  GC[["SLA"]] <- GC[["SLA"]][!is.na(GC[["SLA"]][, "corr.r"]), ]
  GC[["WD"]] <- GC[["WD"]][!is.na(GC[["WD"]][, "corr.r"]), ]
  GC[["Hmax"]] <- GC[["Hmax"]][!is.na(GC[["Hmax"]][, "corr.r"]), ]
  GC[["Seedmass"]] <- GC[["Seedmass"]][!is.na(GC[["Seedmass"]][, "corr.r"]),
    ]
  GC[["Aarea"]] <- GC[["Aarea"]][!is.na(GC[["Aarea"]][, "corr.r"]), ]

  p1 <- my_plot_corr.r(GC[["SLA"]], title = "a) SLA") + theme(plot.margin = unit(c(0,
    0, 0, 0), "mm"), legend.position = "none")
  p2 <- my_plot_corr.r(GC[["WD"]], title = "b)WD") + theme(axis.text.y = element_blank(),
    axis.title.y = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "mm"),
    legend.position = "none")
  p3 <- my_plot_corr.r(GC[["Hmax"]], title = "c) Hmax") + theme(plot.margin = unit(c(0,
    0, 0, 0), "mm"), legend.position = "none")
  p4 <- my_plot_corr.r(GC[["Seedmass"]], title = "d) Seed mass", xlab = "Case studies ranked by coefficient of correlation r") +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0,
      0, 0, 0), "mm"), legend.position = "none")
  p5 <- my_plot_corr.r(GC[["Aarea"]], title = "e) Aarea", xlab = "Case studies ranked by coefficient of correlation r") +
    theme(legend.title = element_blank(), legend.justification = c(0, 0), legend.position = c(1.2,
      0.5), legend.key = element_blank(), plot.margin = unit(c(0, 0, 0, 0),
      "mm"))
  grid.arrange(p1, p2, p3, p4, p5, ncol = 2, nrow = 3, widths = c(1.1, 1))
}

figure_A4 <- function(GIi, GIrgr, GIagr) {

  table_trait <- function(trait) {
    rbind(table_overall(GIi[[trait]]), table_overall.stage(GIi[[trait]]))
  }

  SLA <- table_trait("SLA")
  WD <- table_trait("WD")
  Hmax <- table_trait("Hmax")
  Seedmass <- table_trait("Seedmass")
  Aarea <- table_trait("Aarea")

  p1 <- my_plot_overall(SLA, title = "a) SLA") 
  p2 <- my_plot_overall(WD, title = "b) Wood density") 
  p3 <- my_plot_overall(Hmax, title = "c) Hmax") 
  p4 <- my_plot_overall(Seedmass, title = "d) Seed mass") 
  p5 <- my_plot_overall(Aarea, title = "e) Aarea") 

  p1 <- p1 + theme(plot.margin  = unit(c(2.5, 2, 1.5, 0), "mm"), axis.title.x = element_blank(),axis.title.y = element_text(colour = "white"))
  p2 <- p2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                   plot.margin  = unit(c(2.5, 0, 1.5, 1), "mm"), axis.title.x = element_blank())
  p3 <- p3 + theme(plot.margin  = unit(c(2.5, 2, 1.5, 0), "mm"), axis.title.x = element_blank())
  p4 <- p4 + theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
                   plot.margin  = unit(c(2.5, 0, -3.5, 1), "mm"))
  p5 <- p5 + theme(plot.margin  = unit(c(2.5, 2, 1.5, 0), "mm"), axis.title.y = element_text(colour = "white"))
  
  # Missing data for some stages, that's fine but causes a warning message
  suppressWarnings({
    grid.arrange(p1, p2, p3, p4, p5, ncol = 2, nrow = 3, widths = c(1.3, 1),heights = c(1, 1, 1.1))
  })
}

figure_A5 <- function(RIi, RCi) {

  fit_model <- function(trait, RIi , RCi) {
    ret <- list()
    ret[["model"]] <- fun_model(RIi[[trait]], RCi[[trait]])
    model <- fit_lmer(RIi[[trait]], "stage")
    ret[["LRT"]] <- model[["LRT"]]
    ret[["PVAL"]] <- model[["PVAL"]]
    ret
  }

  traits <- c("SLA", "WD", "Hmax", "Seedmass", "Aarea")
  fits <- lapply(traits, fit_model, RIi , RCi)
  names(fits) <- traits

  figure_panels_traits_model(fits,
      colors = c("grey", "black"), width = 0.6,
      labels = c("conservative dataset", "entire dataset"),
      breaks =   c("ideal", "complete"),
      category_variable = "stress")
}

figure_A6 <- function(GC, trait1, trait2, titles) {
  par(mfcol = c(1, 2))
  par(mar = c(2, 5, 2, 0))

  coeff.plot.multiple3(GC[[trait1]], params = rev(c("stageRGRseedling", "stageRGRjuvenile",
    "stageRGRsapling", "stageRGRadult", "stageRGRmix", "RGRGR(Di)", "RGRGR(Hi)",
    "RGRGR(Mi)", "RGRRGR(CSAi)", "RGRRGR(Di)", "RGRRGR(Hi)", "RGRRGR(Mi)",
    "RGRRGR(Vi)", "growth.formtree", "growth.formwoody", "growth.formacross growth form",
    "experimentcontrol", "experimentdatabase", "experimentfield", "experimentnature")),
    labels = rev(c("seedling", "juvenile", "sapling", "adult", "all stage",
      "GR(D)", "GR(H)", "GR(M)", "RGR(CSA)", "RGR(D)", "RGR(H)", "RGR(M)",
      "RGR(V)", "tree", "woody", "across GF", "greenhouse", "database", "field exp",
      "forest")), title = paste0(titles[1], ") ", trait1))

  mtext("mod4", side = 2, line = 4.2, cex = 0.8, at = 2.2)
  mtext("mod3", side = 2, line = 4.2, cex = 0.8, at = 6)
  mtext("mod2", side = 2, line = 4.2, cex = 0.8, at = 12)
  mtext("mod1", side = 2, line = 4.2, cex = 0.8, at = 18)

  par(mar = c(2, 1.5, 2, 3.5))
  coeff.plot.multiple3(GC[[trait2]], params = rev(c("stageRGRseedling", "stageRGRjuvenile",
    "stageRGRsapling", "stageRGRadult", "stageRGRmix", "RGRGR(Di)", "RGRGR(Hi)",
    "RGRGR(Mi)", "RGRRGR(CSAi)", "RGRRGR(Di)", "RGRRGR(Hi)", "RGRRGR(Mi)",
    "RGRRGR(Vi)", "growth.formtree", "growth.formwoody", "growth.formacross growth form",
    "experimentcontrol", "experimentdatabase", "experimentfield", "experimentnature")),
    title = paste0(titles[2], ") ", trait2))
}

figure_A6.2 <- function(GC, trait1, titles) {
  par(mfcol = c(1, 2))
  par(mar = c(2, 5, 2, 0))

  coeff.plot.multiple3.1(GC[[trait1]], params = rev(c("stageRGRseedling", "stageRGRjuvenile",
    "stageRGRsapling", "stageRGRadult", "stageRGRmix", "RGRGR(Di)", "RGRGR(Hi)",
    "RGRGR(Mi)", "RGRRGR(CSAi)", "RGRRGR(Di)", "RGRRGR(Hi)", "RGRRGR(Mi)",
    "RGRRGR(Vi)", "experimentcontrol", "experimentdatabase", "experimentfield",
    "experimentnature")), labels = rev(c("seedling", "juvenile", "sapling",
    "adult", "all stage", "GR(D)", "GR(H)", "GR(M)", "RGR(CSA)", "RGR(D)",
    "RGR(H)", "RGR(M)", "RGR(V)", "greenhouse", "database", "field exp", "forest")),
    title = paste0(titles[1], ") ", trait1))

  mtext("mod4", side = 2, line = 4.2, cex = 0.8, at = 2.2)
  mtext("mod2", side = 2, line = 4.2, cex = 0.8, at = 8.5)
  mtext("mod1", side = 2, line = 4.2, cex = 0.8, at = 15.2)
}

figure_A7 <- function(GIi) {
  par(mfrow = c(3, 2))
  par(mar = c(5, 4, 1, 2))

  p1 <- figure_trim.and.fill(GIi[["SLA"]], title = "a) SLA")
  p2 <- figure_trim.and.fill(GIi[["WD"]], title = "b) Wood density")
  p3 <- figure_trim.and.fill(GIi[["Hmax"]], title = "c) Hmax")
  p4 <- figure_trim.and.fill(GIi[["Seedmass"]], title = "d) Seedmass")
  p5 <- figure_trim.and.fill(GIi[["Aarea"]], title = "e) Aarea")
}

figure_A8 <- function(GCi) {

  plotgrowth <- function(data, title) {

    data <- data[, c("growth", "stage", "measurement")]
    data$measurement <- factor(data$measurement, levels = c("Mass", "Diameter",
      "Height", "Other"))

    rgr <- subset(data, data$growth == "RGR")
    abs <- subset(data, data$growth != "RGR")

    counts <- table(rgr$measurement, rgr$stage)
    counts1 <- table(abs$measurement, abs$stage)

    a <- table(rgr$stage)
    a1 <- table(abs$stage)

    n <- max(a, a1)

    par(mar = c(4, 4, 2.5, 0))
    barplot(counts1, xlim = c(0, n + 5), xlab = "", col = c("#9E5F3A", "#E5B38F",
      "#F4395B", "#FA8C3D"), horiz = TRUE, border = NA, yaxt = "n", cex.axis = 0.8)
    axis(2, at = 1:3, labels = c("seedling", "sapling", "adult"), las = 1,
      cex.axis = 1)
    mtext("AGR", side = 3, line = 0, cex = 0.8)
    mtext(title, side = 3, line = 1, at = max(n+5))

    par(mar = c(4, 0.5, 2.5, 3.5))
    barplot(counts, xlim = c(0, n + 5), xlab = "", col = c("#9E5F3A", "#E5B38F",
      "#F4395B", "#FA8C3D"), horiz = TRUE, border = NA, yaxt = "n", cex.axis = 0.8)
    abline(v = 0, col = "black")
    mtext("RGR", side = 3, line = 0, cex = 0.8)
  }

  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 3, 4, byrow = TRUE))
  plotgrowth(GCi[["SLA"]], "a) SLA")
  plotgrowth(GCi[["WD"]], "b) Wood density")
  plotgrowth(GCi[["Hmax"]], "c) Hmax")
  plotgrowth(GCi[["Seedmass"]], "d) Seed mass")
  mtext("Number of correlation extracted", side = 1, line = 2.5, at = 0,
        cex = 0.8)
  plotgrowth(GCi[["Aarea"]], "e) Aarea")
  mtext("Number of correlation extracted", side = 1, line = 2.5, at = 0,
        cex = 0.8)

  plot(1, type = "n", axes = F, xlab = "", ylab = "", bty = "n", xaxt = "n",
    yaxt = "n")
  legend("topright", c("mass", "diameter", "height", "other (volume, etc.)"),
    fill = c("#9E5F3A", "#E5B38F", "#F4395B", "#FA8C3D"), border = c("#9E5F3A",
      "#E5B38F", "#F4395B", "#FA8C3D"), col = c("#9E5F3A", "#E5B38F", "#F4395B",
      "#FA8C3D"), bty = "n")
}

figure_A9 <- function(GIi) {

  LRT_SLA <- fun_OneLR_year(GIi[["SLA"]])
  PVAL_SLA <- fun_Onepvalue_year(GIi[["SLA"]])

  LRT_WD <- fun_OneLR_year(GIi[["WD"]])
  PVAL_WD <- fun_Onepvalue_year(GIi[["WD"]])

  LRT_Hmax <- fun_OneLR_year(GIi[["Hmax"]])
  PVAL_Hmax <- fun_Onepvalue_year(GIi[["Hmax"]])

  LRT_Seedmass <- fun_OneLR_year(GIi[["Seedmass"]])
  PVAL_Seedmass <- fun_Onepvalue_year(GIi[["Seedmass"]])

  LRT_Aarea <- fun_OneLR_year(GIi[["Aarea"]])
  PVAL_Aarea <- fun_Onepvalue_year(GIi[["Aarea"]])

  funnel_SLA_year <- my_funnelplot("a) SLA", ggplot(GIi[["SLA"]], aes(x = year,
    y = corr.r, colour = factor(stage), size = 2, alpha = 0.6))) + geom_point() +
    scale_y_continuous("Correlation coefficient  r", limits = c(-1, 1)) + theme(legend.position = "none") +
    scale_x_continuous("", limits = c(1990, 2015)) + annotate("text",
    x = 1995, y = -0.9, label = paste("LRT:", round(LRT_SLA, 0)), size = 2) +
    annotate("text", x = 1995, y = -1, label = paste("p.value =", round(PVAL_SLA,
      3), "***"), size = 2)


  funnel_WD_year <- my_funnelplot("b) Wood density", ggplot(GIi[["WD"]], aes(x = year,
    y = corr.r, colour = factor(stage), size = 2, alpha = 0.6))) + geom_point() +
    scale_y_continuous("", limits = c(-1, 1)) + scale_x_continuous("", limits = c(1990,
    2015)) + theme(legend.position = "none") + annotate("text", x = 1995, y = -0.9,
    label = paste("LRT:", round(LRT_WD, 0)), size = 2) + annotate("text", x = 1995,
    y = -1, label = paste("p.value =", round(PVAL_WD, 3), "ns"), size = 2)


  funnel_Hmax_year <- my_funnelplot("c) Hmax", ggplot(GIi[["Hmax"]], aes(x = year,
    y = corr.r, colour = factor(stage), size = 2, alpha = 0.6))) + geom_point() +
    scale_y_continuous("Correlation coefficient r", limits = c(-1, 1)) + scale_x_continuous("",
    limits = c(1990, 2015)) + theme(legend.position = "none") + annotate("text",
    x = 1995, y = -0.9, label = paste("LRT:", round(LRT_Hmax, 0)), size = 2) +
    annotate("text", x = 1995, y = -1, label = paste("p.value =", round(PVAL_Hmax,
      3), "ns"), size = 2)


  funnel_Seedmass_year <- my_funnelplot("d) Seed mass", ggplot(GIi[["Seedmass"]],
    aes(x = year, y = corr.r, colour = factor(stage), size = 2, alpha = 0.6))) +
    geom_point() + scale_x_continuous("Year of publication", limits = c(1990,
    2015)) + scale_y_continuous("", limits = c(-1, 1)) + theme(legend.position = "none") +
    annotate("text", x = 1995, y = -0.9, label = paste("LRT:", round(LRT_Seedmass,
      0)), size = 2) + annotate("text", x = 1995, y = -1, label = paste("p.value =",
    round(PVAL_Seedmass, 3), "**"), size = 2)


  funnel_Aarea_year <- my_funnelplot("e) Aarea", ggplot(GIi[["Aarea"]], aes(x = year,
    y = corr.r, colour = factor(stage), size = 2, alpha = 0.6))) + geom_point() +
    scale_alpha(guide = "none") + scale_size(guide = "none") + scale_y_continuous("Correlation coefficient  r",
    limits = c(-1, 1)) + scale_x_continuous("Year of publication", limits = c(1990,
    2015)) + theme(legend.title = element_blank(), legend.justification = c(0,
    0), legend.position = c(1.2, 0.5), legend.key = element_blank(), plot.margin = unit(c(0,
    0, 0, 0), "mm")) + annotate("text", x = 1995, y = -0.9, label = paste("LRT:",
    round(LRT_Aarea, 0)), size = 2) + annotate("text", x = 1995, y = -1, label = paste("p.value =",
    round(PVAL_Aarea, 3), "**"), size = 2)

  # Missing data for some stages, that's fine but causes a warning message
  suppressWarnings({
    grid.arrange(funnel_SLA_year, funnel_WD_year, funnel_Hmax_year, funnel_Seedmass_year,
      funnel_Aarea_year, nrow = 3, ncol = 2)
  })
}


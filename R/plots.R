
mytheme <- function() {
  theme_bw() +
  theme(text = element_text(size = 9, colour = "black"),
        title = element_text(size = 9, hjust = 0),
        axis.title = element_text(size = 9, hjust = 0.5),
        axis.text = element_text(size = 8),
    # required for ggplot2_2.0.0
        axis.line = element_line(colour = 'black', size=0.5, linetype='solid'),
    # required for ggplot2_2.1.0
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        legend.justification = c(1, 0), legend.position = c(1, 0),
        legend.key = element_rect(colour = "white"),
        axis.text.x = element_text(margin=margin(9,5,5,5,"pt")),
        axis.text.y = element_text(margin=margin(9,5,5,5,"pt")),
        axis.ticks = element_line (colour = "black", size = 0.5), 
        axis.ticks.length = unit(-0.08 , "cm"))
  }


## Figure 3 and 4
figure_panels_traits_model <- function(fits, ...) {
  
  coeff.plot <- function(x, title = "",
    limit.x.min = -1, limit.x.max = 1.5,
    vjust.value = 2, width = 0.4,
    colors = c("grey", "black"),
    category_variable = "growth",
    labels = c("RGR, conservative dataset", "RGR, entire dataset"),
    breaks = c("RGR", "AGR")) {
    
    data <- x[["model"]]
    LRT <- x[["LRT"]]
    PVAL <- x[["PVAL"]]
    
    data[["by"]] <- as.factor(data[[category_variable]])
    
    label.n.x <- limit.x.min + 0.88*(limit.x.max - limit.x.min)
    
    dodge <- position_dodge(width = width)
    
    significativite <- ""
    round.value <- 2
    if(!is.na(PVAL)) {
      significativite <- "ns"
      if(PVAL <= 0.01) {
        significativite <- "***"
        round.value <- 4
      } else if(PVAL <= 0.05) {
        significativite <- "*"
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
    limit.x.min = -1, limit.x.max = 1.6,  ...) +
    theme(plot.margin  = unit(c(2.5, 2, 1.5, 0), "mm"), axis.title.x = element_blank(),
      axis.title.y = element_text(colour = "white"))
  
  p2 <- coeff.plot(fits[["WD"]], title = "b) Wood density",
    limit.x.min = -1, limit.x.max = 1.05,  ...) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
      plot.margin  = unit(c(2.5, 0, 1.5, 1), "mm"), axis.title.x = element_blank())
  
  p3 <- coeff.plot(fits[["Hmax"]], title = expression("c) H"[max]),
    limit.x.min = -1, limit.x.max = 1.6, ...) +
    theme(plot.margin  = unit(c(2.5, 2, 1.5, 0), "mm"), axis.title.x = element_blank())
  
  p4 <- coeff.plot(fits[["Seedmass"]], title = "d) Seed mass",
    limit.x.min = -1.5, limit.x.max = 1.6, ...) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
      plot.margin  = unit(c(2.5, 0, -3.5, 1), "mm"))
  
  p5 <- coeff.plot(fits[["Aarea"]], title =  expression("e) A"[area]),
    limit.x.min = -2.4, limit.x.max =3.5, ...)  +
    theme(legend.title = element_blank(), legend.justification = c(0, 0), legend.position = c(1, 0.3),
      legend.key = element_blank(), plot.margin  = unit(c(2.5, 2, 1.5, 0), "mm"),
      axis.title.y = element_text(colour = "white"))
  
  grid.arrange(p1, p2, p3, p4, p5, ncol = 2, nrow = 3, widths = c(1.3, 1),heights = c(1, 1, 1.1))
}


my_plot_corr.r <- function(data1, title = "", xlab = "Case studies ranked by coefficient of correlation r") {

  my_funnelplot2 <- function(title, ggobj, xlab = "", ylab = "Stage") {

    p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + guides(fill = guide_legend(reverse = TRUE)) +
      mytheme() + scale_y_continuous("Correlation coefficient r", limits = c(-1,
      1)) + scale_colour_manual(limits = c("seedling", "sapling", "adult"),
      values = c(seedling = "#D5C9B1", sapling = "#805A3B", adult = "#C60000"),
      breaks = c("seedling", "sapling", "adult")) + geom_hline(yintercept = 0,
      size = 0.3, linetype = "dashed") + theme(text = element_text(size = 9,
      colour = "black"), title = element_text(size = 9, hjust = 0), axis.title = element_text(size = 9,
      hjust = 0.5), axis.text.x = element_blank(), axis.text.y = element_text(size = 9),
      axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(0,
        0, 0, 0), "mm"))
    p
  }

  my_funnelplot2(title, ggplot(data1, aes(x = reorder(id, corr.r), y = corr.r,
    colour = factor(stage)))) + geom_point(aes(size = wi.z), alpha = 0.7) +
    scale_size(guide = "none")

}


my_plot_overall <- function(data, title = "", name.xlab ="Coefficient of correlation r (+SD)", size=2) {

  plot1 <- function(title, ggobj,  xlab = "Effect size", ylab = "Stage") {

    p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + geom_vline(xintercept = 0,
      color = "white") + geom_hline(yintercept = 0, size = 0.3, linetype = "dashed") +
      geom_vline(xintercept = 3.5, size = 0.3, linetype = "dashed") + scale_x_discrete(ylab,
      limit = c("seedling", "sapling", "adult", "total")) + coord_flip() +
      mytheme()
    p
  }

  data["x.coord"] <- max(data$corr.r + data$SD)

  plot1(title,
    ggplot(data, aes(stage, corr.r, ymin = corr.r - SD, ymax = corr.r + SD, colour=stage))) +
    geom_point(aes(x = stage, y = corr.r), size=size) +
    geom_errorbar(aes(x = stage,y = corr.r), width = 0) +
    geom_text(aes(stage, 0.8, label = paste("n=",freq)), size = 2, data = data, parse = F, position = "identity", vjust = 0.2, hjust = 0)+
    scale_y_continuous(name.xlab, limits = c(-1, 1)) +
    scale_colour_manual(values = c(seedling = "#D5C9B1", sapling = "#805A3B", adult = "#C60000", total ="black"), breaks = c("seedling", "sapling",
      "adult", "total"), labels = c("seedling", "sapling", "adult", "total")
    )
}


## Plot Appendix:
coeff.plot.multiple3 <- function(data, params, labels = NA, xlab = "Effect size (z) +CI 95%",
  title = "") {
  fun_List_N <- function(data) {
    c(length(na.omit(data$coef[data$stageRGR == "seedling"])), length(na.omit(data$coef[data$stageRGR ==
      "juvenile"])), length(na.omit(data$coef[data$stageRGR == "sapling"])),
      length(na.omit(data$coef[data$stageRGR == "adult"])), length(na.omit(data$coef[data$stageRGR ==
        "mix"])), length(na.omit(data$coef[data$RGR == "GR(Di)"])), length(na.omit(data$coef[data$RGR ==
        "GR(Hi)"])), length(na.omit(data$coef[data$RGR == "GR(Mi)"])),
      length(na.omit(data$coef[data$RGR == "RGR(CSAi)"])), length(na.omit(data$coef[data$RGR ==
        "RGR(Di)"])), length(na.omit(data$coef[data$RGR == "RGR(Hi)"])),
      length(na.omit(data$coef[data$RGR == "RGR(Mi)"])), length(na.omit(data$coef[data$RGR ==
        "RGR(Vi)"])), length(na.omit(data$coef[data$growth.form == "tree"])),
      length(na.omit(data$coef[data$growth.form == "woody"])), length(na.omit(data$coef[data$growth.form ==
        "across growth form"])), length(na.omit(data$coef[data$experiment ==
        "control"])), length(na.omit(data$coef[data$experiment == "database"])),
      length(na.omit(data$coef[data$experiment == "field"])), length(na.omit(data$coef[data$experiment ==
        "nature"])))
  }

  fun_AICm <- function(data) {
    m <- lmer(corr.z ~ stage + (1 | id), data = data, weights = nb.sp, REML = FALSE,
      control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
        check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
          tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
          tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))  # model avec les intercepts qui sont a 1
    gr <- lmer(corr.z ~ RGR + (1 | id), data = data, weights = nb.sp, REML = FALSE,
      control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
        check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
          tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
          tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))
    veg <- lmer(corr.z ~ growth.form + (1 | id), data = data, weights = nb.sp,
      REML = FALSE, control = lmerControl(check.nobs.vs.nlev = "ignore",
        check.nobs.vs.rankZ = "ignore", check.nobs.vs.nRE = "ignore", optimizer = "bobyqa",
        check.conv.grad = .makeCC("ignore", tol = 0.002, relTol = NULL),
        check.conv.singular = .makeCC(action = "ignore", tol = 1e-04),
        check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))
    exp <- lmer(corr.z ~ experiment + (1 | id), data = data, weights = nb.sp,
      REML = FALSE, control = lmerControl(check.nobs.vs.nlev = "ignore",
        check.nobs.vs.rankZ = "ignore", check.nobs.vs.nRE = "ignore", optimizer = "bobyqa",
        check.conv.grad = .makeCC("ignore", tol = 0.002, relTol = NULL),
        check.conv.singular = .makeCC(action = "ignore", tol = 1e-04),
        check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))
    c(AIC(exp), AIC(gr), AIC(veg), AIC(m))
  }

  name <- list(1:4, c("value", "x", "y"))
  aic <- as.data.frame(matrix(nrow = length(name[[1]]), ncol = length(name[[2]]),
    dimnames = name))
  aic["value"] <- data.frame(fun_AICm(data))
  aic["x"] <- c(2.3, 2.3, 2.3, 2.3)
  aic["y"] <- c(1, 5, 8, 16)
  aic$value <- round(aic$value, digits = 0)

  # nom de la liste
  List_N <- fun_List_N(data)

  # creation du plot
  dat <- fun_model_multiple3(data)
  dat <- dat[match(params, dat$params), ]  # je mets dans dat la liste des parametres qui matche les parametres indique danq ma fonction coefplot
  if (is.null(labels))
    labels <- dat$params  #si labels=NULL dans ma fonction coefplot alors les labels sont ceux du tableau de donnee

  plot(0, xaxt = "n", yaxt = "n", bty = "n", pch = "", ylab = "", xlab = "",
    xlim = c(-3, 3), ylim = c(1, 20))
  # plot(0seq_len(nrow(dat)), xlab='', ylab='',xlim=c(-2,3))
  rect(-4, 0, 5, 4.5, col = "grey", border = "transparent", density = 70, xpd = FALSE)
  rect(-4, 7.5, 5, 15.5, col = "grey", border = "transparent", density = 70,
    xpd = FALSE)
  points(dat[, "mean"], seq_len(nrow(dat)), xlab = "", ylab = "", xlim = c(-3,
    3), yaxt = "n", pch = 21, bg = "black", cex.lab = 0.7, cex.axis = 0.7,
    panel.first = {
      abline(v = 0, lty = 3)
      segments(dat$lower, seq_len(nrow(dat)), dat$upper, seq_len(nrow(dat)),
        lend = 1)
    })
  axis(2, at = seq_len(nrow(dat)), labels = labels, las = 1, cex.axis = 0.7,tck=0.02)
  box()
  text(paste("n=", rev(List_N)), x = dat$lower - 0.5, y = seq_len(nrow(dat)),
    cex = 0.5)
  text(paste("AIC=", aic$value), x = aic$x, y = aic$y, cex = 0.7)
  mtext(xlab, side = 1, line = 0.5, cex = 0.8)
  mtext(title, side = 3, line = 0.5, cex = 0.9, at = 0.5)
}

coeff.plot.multiple3.1 <- function(data, params, labels = NA, xlab = "Effect size (z) +CI 95%",
  title = "") {
  
  fun_List_N <- function(data) {
    c(length(na.omit(data$coef[data$stageRGR == "seedling"])), length(na.omit(data$coef[data$stageRGR ==
      "juvenile"])), length(na.omit(data$coef[data$stageRGR == "sapling"])),
      length(na.omit(data$coef[data$stageRGR == "adult"])), length(na.omit(data$coef[data$stageRGR ==
        "mix"])), length(na.omit(data$coef[data$RGR == "GR(Di)"])), length(na.omit(data$coef[data$RGR ==
        "GR(Hi)"])), length(na.omit(data$coef[data$RGR == "GR(Mi)"])),
      length(na.omit(data$coef[data$RGR == "RGR(CSAi)"])), length(na.omit(data$coef[data$RGR ==
        "RGR(Di)"])), length(na.omit(data$coef[data$RGR == "RGR(Hi)"])),
      length(na.omit(data$coef[data$RGR == "RGR(Mi)"])), length(na.omit(data$coef[data$RGR ==
        "RGR(Vi)"])), length(na.omit(data$coef[data$experiment == "control"])),
      length(na.omit(data$coef[data$experiment == "database"])), length(na.omit(data$coef[data$experiment ==
        "field"])), length(na.omit(data$coef[data$experiment == "nature"])))
  }

  fun_AICm <- function(data) {
    m <- lmer(corr.z ~ stage + (1 | id), data = data, weights = nb.sp, REML = FALSE,
      control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
        check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
          tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
          tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))  # model avec les intercepts qui sont a 1
    gr <- lmer(corr.z ~ RGR + (1 | id), data = data, weights = nb.sp, REML = FALSE,
      control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.rankZ = "ignore",
        check.nobs.vs.nRE = "ignore", optimizer = "bobyqa", check.conv.grad = .makeCC("ignore",
          tol = 0.002, relTol = NULL), check.conv.singular = .makeCC(action = "ignore",
          tol = 1e-04), check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))
    exp <- lmer(corr.z ~ experiment + (1 | id), data = data, weights = nb.sp,
      REML = FALSE, control = lmerControl(check.nobs.vs.nlev = "ignore",
        check.nobs.vs.rankZ = "ignore", check.nobs.vs.nRE = "ignore", optimizer = "bobyqa",
        check.conv.grad = .makeCC("ignore", tol = 0.002, relTol = NULL),
        check.conv.singular = .makeCC(action = "ignore", tol = 1e-04),
        check.conv.hess = .makeCC(action = "ignore", tol = 1e-06)))
    c(AIC(exp), AIC(gr), AIC(m))
  }

  name <- list(1:3, c("value", "x", "y"))
  aic <- as.data.frame(matrix(nrow = length(name[[1]]), ncol = length(name[[2]]),
    dimnames = name))
  aic["value"] <- data.frame(fun_AICm(data))
  aic["x"] <- c(2.3, 2.3, 2.3)
  aic["y"] <- c(1, 5, 13)
  aic$value <- round(aic$value, digits = 0)

  # nom de la liste
  List_N <- fun_List_N(data)

  # creation du plot
  dat <- fun_model_multiple3.1(data)
  dat <- dat[match(params, dat$params), ]
  # je mets dans dat la liste des parametres qui matche les parametres indique
  # danq ma fonction coefplot
  if (is.null(labels))
    labels <- dat$params
  # si labels=NULL dans ma fonction coefplot alors les labels sont ceux du
  # tableau de donnee

  plot(0, xaxt = "n", yaxt = "n", bty = "n", pch = "", ylab = "", xlab = "",
    xlim = c(-3, 3), ylim = c(1, 17))
  # plot(0seq_len(nrow(dat)), xlab='', ylab='',xlim=c(-2,3))
  rect(-4, 0, 5, 4.5, col = "grey", border = "transparent", density = 70, xpd = FALSE)
  rect(-4, 12.5, 5, 19, col = "grey", border = "transparent", density = 70, xpd = FALSE)
  points(dat[, "mean"], seq_len(nrow(dat)), xlab = "", ylab = "", xlim = c(-3,
    3), yaxt = "n", pch = 21, bg = "black", cex.lab = 0.7, cex.axis = 0.7,
    panel.first = {
      abline(v = 0, lty = 3)
      segments(dat$lower, seq_len(nrow(dat)), dat$upper, seq_len(nrow(dat)),
        lend = 1)
    })
  axis(2, at = seq_len(nrow(dat)), labels = labels, las = 1, cex.axis = 0.7,tck=0.02)
  box()
  text(paste("n=", rev(List_N)), x = dat$lower - 0.5, y = seq_len(nrow(dat)),
    cex = 0.5)
  text(paste("AIC=", aic$value), x = aic$x, y = aic$y, cex = 0.7)
  mtext(xlab, side = 1, line = 0.5, cex = 0.8)
  mtext(title, side = 3, line = 0.5, cex = 0.9, at = 0.5)
}


my_funnelplot <- function(title, ggobj, xlab = "", ylab = "") {
  p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + mytheme() + scale_color_manual(name = "Stage",
    limits = c("seedling", "sapling", "adult"), breaks = c("seedling", "sapling",
      "adult"), values = c(seedling = "#D5C9B1",sapling = "#805A3B", adult = "#C60000")) +
    geom_hline(yintercept = 0, color = "grey")
  p
}


figure_trim.and.fill <- function(x, title) {
  res <- rma(corr.z, vr.z, data = x[!is.na(x$corr.z) & !is.na(x$vr.z), ])
  ### carry out trim-and-fill analysis
  taf <- trimfill(res)
  funnel(taf, xlab = "Correlation", main = title)
}

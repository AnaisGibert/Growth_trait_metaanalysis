
mytheme <- theme_bw() + theme(text = element_text(size = 9, colour = "black"),
  title = element_text(size = 9, hjust = 0), axis.title = element_text(size = 9,
    hjust = 0.5), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
  panel.background = element_blank(), legend.justification = c(1, 0), legend.position = c(1,
    0))

my_plot_1 <- function(title, ggobj, xlab = expression(paste("")), ylab = "Number of correlations extracted") {
  p <- ggobj + geom_bar(alpha = 0.95) + coord_cartesian(ylim = ylim) + coord_flip() +
    labs(title = title) + xlab(xlab) + ylab(ylab) + scale_x_discrete() + theme(text = element_text(size = 9),
    axis.text.x = element_text(size = 9, angle = 0, vjust = 1)) + mytheme +
    theme(axis.title = element_text(size = 10, hjust = 0.5)) + scale_fill_manual(values = c(juvenile = "#538936",
    sapling = "#F57C34", adult = "#E6224C"), breaks = c("juvenile", "sapling",
    "adult"), labels = c("juvenile", "sapling", "adult"))
  p

}


coeff.plot.2 <- function(data, data.complete, data.ideal, LRT, PVAL, title = "",
  significativite = "", round.value, limit.x.min, limit.x.max, limit.x.text,
  limit.x.n, limit.y.text.l1, limit.y.text.l2, vjust.value, color1, color2) {
  data.complete["x.coord"] <- limit.x.n

  data.ideal["x.coord"] <- limit.x.n

  my_size_effect_plot1 <- function(title, ggobj, xlab = "Effect size + CI (95%) ",
    ylab = "Stages") {
    p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + geom_vline(yintercept = 0,
      color = "white") + geom_hline(yintercept = 0, size = 0.3, linetype = "dashed") +
      scale_x_discrete("stage", limit = c("seedling", "juvenile", "sapling",
        "adult", "mix")) + coord_flip() + mytheme + theme(legend.position = "none")
    p
  }

  dodge <- position_dodge(width = 0.4)

  my_size_effect_plot1(title, ggplot(data, aes(stage, Inte, ymin = Inte - 1.96 *
    SE, ymax = Inte + 1.96 * SE, color = factor(stress)))) + geom_point(aes(x = stage,
    y = Inte), position = dodge) + geom_errorbar(aes(x = stage, y = Inte),
    size = 0.4, width = 0.2, position = dodge) + scale_y_continuous("Effect size (z) + CI 95%",
    limits = c(limit.x.min, limit.x.max)) + geom_text(aes(stage, x.coord, label = paste("(n=",
    N, ")"), color = factor(stress), group = "all"), size = 2, data = data.complete,
    parse = F, position = "identity", vjust = +vjust.value, hjust = 0) + geom_text(aes(stage,
    x.coord, label = paste("(n=", N, ")"), color = factor(stress), group = "ideal"),
    size = 2, data = data.ideal, parse = F, position = "identity", vjust = -vjust.value,
    hjust = 0) + scale_color_manual(values = c(color1, color2)) + annotate("text",
    x = limit.y.text.l1, y = limit.x.text, label = paste("LRT:", round(LRT,
      2)), size = 2) + annotate("text", x = limit.y.text.l2, y = limit.x.text,
    label = paste("p.value =", round(PVAL, round.value), significativite),
    size = 2)
}



coeff.plot <- function(data, data.complete, data.ideal, LRT, PVAL, title = "",
  significativite = "", round.value, limit.x.min, limit.x.max, limit.x.text,
  limit.x.n, limit.y.text.l1, limit.y.text.l2, vjust.value, color1, color2) {

  data.complete["x.coord"] <- limit.x.n
  data.ideal["x.coord"] <- limit.x.n

  plot1 <- function(title, ggobj, xlab = "Effect size + CI (95%) ", ylab = "Stages") {

    p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + geom_vline(yintercept = 0,
      color = "white") + geom_hline(yintercept = 0, size = 0.3, linetype = "dashed") +
      scale_x_discrete("stage", limit = c("juvenile", "sapling", "adult")) +
      coord_flip() + mytheme + theme(legend.position = "none")
    p
  }

  dodge <- position_dodge(width = 0.4)

  plot1(title, ggplot(data, aes(stage, Inte, ymin = Inte - 1.96 * SE, ymax = Inte +
    1.96 * SE, color = factor(stress)))) + geom_point(aes(x = stage, y = Inte),
    position = dodge) + geom_errorbar(aes(x = stage, y = Inte), size = 0.4,
    width = 0.2, position = dodge) + scale_y_continuous("Effect size (z) + CI 95%",
    limits = c(limit.x.min, limit.x.max)) + geom_text(aes(stage, x.coord, label = paste("(n=",
    N, ")"), color = factor(stress), group = "all"), size = 2, data = data.complete,
    parse = F, position = "identity", vjust = +vjust.value, hjust = 0) + geom_text(aes(stage,
    x.coord, label = paste("(n=", N, ")"), color = factor(stress), group = "ideal"),
    size = 2, data = data.ideal, parse = F, position = "identity", vjust = -vjust.value,
    hjust = 0) + scale_color_manual(values = c(color1, color2)) + annotate("text",
    x = limit.y.text.l1, y = limit.x.text, label = paste("LRT:", round(LRT,
      2)), size = 2) + annotate("text", x = limit.y.text.l2, y = limit.x.text,
    label = paste("p.value =", round(PVAL, round.value), significativite),
    size = 2)
}


coeff.plot.ideal.2 <- function(data.ideal, LRT, PVAL, title = "", significativite = "",
  round.value, limit.x.min, limit.x.max, limit.x.text, limit.x.n, limit.y.text.l1,
  limit.y.text.l2, vjust.value, color1) {

  data.ideal["x.coord"] <- limit.x.n

  my_size_effect_plot2 <- function(title, ggobj, xlab = "Effect size + CI (95%) ",
    ylab = "Stages") {
    p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + geom_vline(yintercept = 0,
      color = "white") + geom_hline(yintercept = 0, size = 0.3, linetype = "dashed") +
      scale_x_discrete("stage", limit = c("seedling", "juvenile", "sapling",
        "adult", "mix")) + coord_flip() + mytheme + theme(legend.position = "none")
    p
  }

  my_size_effect_plot2(title, ggplot(data.ideal, aes(stage, Inte, ymin = Inte -
    1.96 * SE, ymax = Inte + 1.96 * SE, color = factor(stress)))) + geom_point(aes(x = stage,
    y = Inte)) + geom_errorbar(aes(x = stage, y = Inte), size = 0.4, width = 0.2) +
    scale_y_continuous("Effect size (z) + CI 95%", limits = c(limit.x.min,
      limit.x.max)) + geom_text(aes(stage, x.coord, label = paste("(n=",
    N, ")"), color = factor(stress), group = "ideal"), size = 2, data = data.ideal,
    parse = F, position = "identity", vjust = -vjust.value, hjust = 0) + scale_color_manual(values = c(color1)) +
    annotate("text", x = limit.y.text.l1, y = limit.x.text, label = paste("LRT:",
      round(LRT, 2)), size = 2) + annotate("text", x = limit.y.text.l2, y = limit.x.text,
    label = paste("p.value =", round(PVAL, round.value), significativite),
    size = 2)
}


coeff.plot.ideal <- function(data.ideal, LRT, PVAL, title = "", significativite = "",
  round.value, limit.x.min, limit.x.max, limit.x.text, limit.x.n, limit.y.text.l1,
  limit.y.text.l2, vjust.value, color1) {

  data.ideal["x.coord"] <- limit.x.n

  my_size_effect_plot2 <- function(title, ggobj, xlab = "Effect size + CI (95%) ",
    ylab = "Stages") {
    p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + geom_vline(yintercept = 0,
      color = "white") + geom_hline(yintercept = 0, size = 0.3, linetype = "dashed") +
      scale_x_discrete("stage", limit = c("juvenile", "sapling", "adult")) +
      coord_flip() + mytheme + theme(legend.position = "none")
    p
  }

  my_size_effect_plot2(title, ggplot(data.ideal, aes(stage, Inte, ymin = Inte -
    1.96 * SE, ymax = Inte + 1.96 * SE, color = factor(stress)))) + geom_point(aes(x = stage,
    y = Inte)) + geom_errorbar(aes(x = stage, y = Inte), size = 0.4, width = 0.2) +
    scale_y_continuous("Effect size (z) + CI 95%", limits = c(limit.x.min,
      limit.x.max)) + geom_text(aes(stage, x.coord, label = paste("(n=",
    N, ")"), color = factor(stress), group = "ideal"), size = 2, data = data.ideal,
    parse = F, position = "identity", vjust = -vjust.value, hjust = 0) + scale_color_manual(values = c(color1)) +
    annotate("text", x = limit.y.text.l1, y = limit.x.text, label = paste("LRT:",
      round(LRT, 2)), size = 2) + annotate("text", x = limit.y.text.l2, y = limit.x.text,
    label = paste("p.value =", round(PVAL, round.value), significativite),
    size = 2)
}

plot_stage <- function(title, xlab, ggobj, ylab = "ID") {
  p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + geom_errorbarh(aes(xmin = size.min,
    xmax = size.max, y = id), size = 0.7, width = 0.05, height = 0.2) + mytheme +
    scale_colour_manual(limits = c("seedling", "juvenile", "sapling", "adult",
      "mix"), values = c(seedling = "green", juvenile = "blue", sapling = "orange",
      adult = "red", mix = "grey"), name = "stage", breaks = c("seedling",
      "juvenile", "sapling", "adult", "mix")) + scale_linetype_manual(limits = c("tree",
    "woody", "across growth form"), values = c("solid", "dashed", "dotted"),
    name = "growth.form", breaks = c("tree", "woody", "across growth form")) +
    theme(axis.text = element_text(size = 6), legend.position = "none")
  p
}

my_plot_corr.r <- function(data1, title = "", xlab = "Case studies ranked by coefficient of correlation r") {

  my_funnelplot2 <- function(title, ggobj, xlab = "", ylab = "") {

    p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + guides(fill = guide_legend(reverse = TRUE)) +
      mytheme + scale_y_continuous("Correlation coefficient r", limits = c(-1,
      1)) + scale_colour_manual(limits = c("juvenile", "sapling", "adult"),
      values = c(juvenile = "#538936", sapling = "#F57C34", adult = "#E6224C"),
      breaks = c("juvenile", "sapling", "adult")) + geom_hline(yintercept = 0,
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


my_plot_overall <- function(data, title = "") {

  plot1 <- function(title, ggobj, xlab = "Effect size + CI (95%) ", ylab = "") {

    p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + geom_vline(yintercept = 0,
      color = "white") + geom_hline(yintercept = 0, size = 0.3, linetype = "dashed") +
      geom_vline(xintercept = 3.5, size = 0.3, linetype = "dashed") + scale_x_discrete("stage",
      limit = c("juvenile", "sapling", "adult", "total")) + coord_flip() +
      mytheme
    p
  }

  data["x.coord"] <- max(data$corr.r + data$SD)

  plot1(title, ggplot(data, aes(stage, corr.r, ymin = corr.r - SD, ymax = corr.r +
    SD))) + geom_point(aes(x = stage, y = corr.r)) + geom_errorbar(aes(x = stage,
    y = corr.r), size = 0.4, width = 0.2) + scale_y_continuous("Coefficient of correlation r (+SD) ",
    limits = c(-1, 1))
}

my_plot_overall_gr <- function(data, title = "") {

  plot1 <- function(title, ggobj, xlab = "Effect size + CI (95%) ", ylab = "") {
    p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + geom_vline(yintercept = 0,
      color = "white") + geom_hline(yintercept = 0, size = 0.3, linetype = "dashed") +
      geom_vline(xintercept = 3.5, size = 0.3, linetype = "dashed") + scale_x_discrete("stage",
      limit = c("juvenile", "sapling", "adult", "total")) + coord_flip() +
      mytheme
    p
  }

  dodge <- position_dodge(width = 0.4)

  data["v.coord"] <- NA




  data1 <- subset(data, data$growth == "RGR")
  data2 <- subset(data, data$growth == "AbsGR")

  p <- plot1(title, ggplot(data, aes(stage, corr.r, ymin = corr.r - SD, ymax = corr.r +
    SD, color = factor(growth)))) + geom_point(aes(x = stage, y = corr.r),
    position = dodge) + geom_errorbar(aes(x = stage, y = corr.r), position = dodge,
    size = 0.4, width = 0.2) + scale_y_continuous("Coefficient of correlation r (+SD) ",
    limits = c(-1, 1)) + coord_flip() + geom_text(aes(stage, 0.8, label = paste("(n=",
    freq, ")"), color = factor(growth), group = "RGR"), size = 2, data = data1,
    parse = F, position = "identity", vjust = -0.5, hjust = 0) + geom_text(aes(stage,
    0.8, label = paste("(n=", freq, ")"), color = factor(growth), group = "AbsGR"),
    size = 2, data = data2, parse = F, position = "identity", vjust = +0.5,
    hjust = 0) + scale_colour_manual(values = c(RGR = "#9E5F3A", AbsGR = "#F4395B"),
    breaks = c("RGR", "AbsGR")) + theme(legend.position = "none")

  p
}

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
  axis(2, at = seq_len(nrow(dat)), labels = labels, las = 1, cex.axis = 0.7)
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
  axis(2, at = seq_len(nrow(dat)), labels = labels, las = 1, cex.axis = 0.7)
  box()
  text(paste("n=", rev(List_N)), x = dat$lower - 0.5, y = seq_len(nrow(dat)),
    cex = 0.5)
  text(paste("AIC=", aic$value), x = aic$x, y = aic$y, cex = 0.7)
  mtext(xlab, side = 1, line = 0.5, cex = 0.8)
  mtext(title, side = 3, line = 0.5, cex = 0.9, at = 0.5)
}


my_plot_3 <- function(title, ggobj, xlab = expression(paste("")), ylab = "growth measurements") {
  p <- ggobj + geom_bar() + coord_cartesian(ylim = ylim) + coord_flip() + labs(title = title) +
    xlab(xlab) + ylab(ylab) + scale_x_discrete("growth measurements", limits = c("GR(Di)",
    "GR(Hi)", "GR(Mi)", "GR(Shoot)", "RGR(?)", "RGR(CSAi)", "RGR(Di)", "RGR(Hi)",
    "RGR(Mi)", "RGR(Vi)"), labels = c(`GR(Di)` = "GR(D)", `GR(Hi)` = "GR(H)",
    `GR(Mi)` = "GR(M)", `GR(Shoot)` = "GR(Shoot)", `RGR(?)` = "RGR", `RGR(CSAi)` = "RGR(CSA)",
    `RGR(Di)` = "RGR(D)", `RGR(Hi)` = "RGR(H)", `RGR(Mi)` = "RGR(M)", `RGR(Vi)` = "RGR(V)")) +
    theme(text = element_text(size = 9), axis.text.x = element_text(size = 9,
      angle = 0, vjust = 1)) + mytheme + theme(axis.title = element_text(size = 10,
    hjust = 0.5)) + scale_fill_manual("", values = c(juvenile = "#538936",
    sapling = "#F57C34", adult = "#E6224C"), breaks = c("juvenile", "sapling",
    "adult"), labels = c("juvenile", "sapling", "adult"))

  p
}


my_plot_3b <- function(title, ggobj, xlab = expression(paste("")), ylab = "growth measurements") {
  p <- ggobj + geom_bar() + coord_cartesian(ylim = ylim) + coord_flip() + labs(title = title) +
    xlab(xlab) + ylab(ylab) + scale_x_discrete("stage", limit = c("juvenile",
    "sapling", "adult")) + theme(text = element_text(size = 9), axis.text.x = element_text(size = 9),
    axis.text.x = element_text(size = 9, angle = 0, vjust = 1)) + mytheme +
    theme(axis.title = element_text(size = 10, hjust = 0.5)) + scale_fill_manual("growth",
    limits = c("AbGR", "RGR"), values = c(AbGR = "grey", RGR = "blue"))
  p
}

my_funnelplot <- function(title, ggobj, xlab = "", ylab = "") {
  p <- ggobj + labs(title = title) + xlab(xlab) + ylab(ylab) + mytheme + scale_color_manual(name = "stage",
    limits = c("juvenile", "sapling", "adult"), breaks = c("juvenile", "sapling",
      "adult"), values = c(juvenile = "#538936", sapling = "#F57C34", adult = "#E6224C")) +
    geom_hline(xintercept = 0, color = "grey")
  p
}


figure_trim.and.fill <- function(x, title) {
  res <- rma(corr.z, vr.z, data = x[!is.na(x$corr.z) & !is.na(x$vr.z), ])
  ### carry out trim-and-fill analysis
  taf <- trimfill(res)
  funnel(taf, axis = "corr.z", main = title)
}

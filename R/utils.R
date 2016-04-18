
## for table_overall.stage, for FigA4
wtd.mean <- function(x, weights = NULL, normwt = "ignored", na.rm = TRUE) {
  if (!length(weights))
    return(mean(x, na.rm = na.rm))
  if (na.rm) {
    s <- !is.na(x + weights)
    x <- x[s]
    weights <- weights[s]
  }
  sum(weights * x)/sum(weights)
}

## for table_overall.stage, for FigA4
wtd.var <- function(x, weights = NULL, normwt = FALSE, na.rm = TRUE, method = c("unbiased",
  "ML")) {
  method <- match.arg(method)
  if (!length(weights)) {
    if (na.rm)
      x <- x[!is.na(x)]
    return(var(x))
  }
  if (na.rm) {
    s <- !is.na(x + weights)
    x <- x[s]
    weights <- weights[s]
  }
  if (normwt)
    weights <- weights * length(x)/sum(weights)
  if (method == "ML")
    return(as.numeric(stats::cov.wt(cbind(x), weights, method = "ML")$cov))
  sw <- sum(weights)
  xbar <- sum(weights * x)/sw
  sum(weights * ((x - xbar)^2))/(sw - (if (normwt)
    sum(weights^2)/sw else 1))
}


## for table_overall.stage, for FigA4
count <- function(data, var) {

  ret <- as.data.frame(table(data[[var]]))
  names(ret) <- c(var, "freq")
  ret
}

# Used in Data processing
rename <- function(data, replace, warn_missing = TRUE, warn_duplicated = TRUE) {

  from <- names(replace)
  to <- replace

  x <- names(data)
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  if (!is.atomic(x)) {
    stop("`x` must be an atomic vector.")
  }

  mapidx <- match(x, from)
  mapidxNA <- is.na(mapidx)
  from_found <- sort(unique(mapidx))
  if (warn_missing && length(from_found) != length(from)) {
    message("The following `from` values were not present in `x`: ", paste(from[!(1:length(from) %in%
      from_found)], collapse = ", "))
  }
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]

  names(data) <- x

  duplicated_names <- names(data)[duplicated(names(data))]
  if (warn_duplicated && (length(duplicated_names) > 0L)) {
    duplicated_names_message <- paste0("`", duplicated_names, "`", collapse = ", ")
    warning("The rename operation has created duplicates for the ", "following name(s): (",
      duplicated_names_message, ")", call. = FALSE)
  }
  data
}

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


default_lmer_control <- function() {
  # Ignore warnings when we have fewer observations than levels
  lmerControl(check.nobs.vs.nlev = "ignore",
    check.nobs.vs.rankZ = "ignore",
    check.nobs.vs.nRE = "ignore")
}

default_lmer_control.b <- function() {
  # Ignore warnings when we have fewer observations than levels
  lmerControl(check.nobs.vs.nlev = "ignore", 
    check.nobs.vs.rankZ = "ignore",
    check.nobs.vs.nRE = "ignore", 
    optimizer = "bobyqa", 
    check.conv.grad = .makeCC("ignore",
    tol = 0.002, relTol = NULL), 
    check.conv.singular = .makeCC(action = "ignore", tol = 1e-04), 
    check.conv.hess = .makeCC(action = "ignore", tol = 1e-06))
}

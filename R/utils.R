
# Imported from Hmisc to save loading that package ‘wtd.mean’, ‘wtd.var’
# compute weighted means, and variances For full info see package
# ?Hmisc::wtd.mean
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

# Count records in data set Behaves similarly to plyr::count
count <- function(data, var) {

  ret <- as.data.frame(table(data[[var]]))
  names(ret) <- c(var, "freq")
  ret
}

# Renames columns of data frame imported and modified from plyr::rename
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




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


list_by_sym <- function(df) {
  split(df, df$SymArd)
}

list_by_fert <- function(df) {
  split(df, df$Fertilization)
}





simplefreq <- function (x,y=NULL,margin=1,html=NULL) {
  if (is.null(y)) {simplefreq.tab<-cbind(data.frame(addmargins(
    prop.table(table(x)))), 
    data.frame(addmargins(table(x)))[2]) 
  colnames(simplefreq.tab)<-c(deparse(substitute(x)),"Percentage","Count") 
  RESTABLE<-simplefreq.tab 
  
  MINCHI2=NULL 
  MINOR=NULL 
  CHI2=NULL 
  } 
  if (!is.null(y)) { 
    simplefreq.tab0<-table(x,y) 
    simplefreq.chi<-chisq.test(simplefreq.tab0,
                               correct=FALSE) 
    if (margin==1) 
    { 
      simplefreq.tab<-cbind(addmargins(prop.table(
        addmargins(simplefreq.tab0,1),1),2), 
        rowSums(addmargins(simplefreq.tab0,1))) 
      colnames(simplefreq.tab)<-c(colnames(simplefreq.tab0),"Total","Count") 
      rownames(simplefreq.tab)<-c(rownames(simplefreq.tab0),"All") 
    } 
    if (margin==2) 
    { 
      simplefreq.tab<-rbind(addmargins(prop.table(
        addmargins(simplefreq.tab0,2),2),1), 
        t(colSums(addmargins(simplefreq.tab0,2)))) 
      rownames(simplefreq.tab)<-c(rownames(simplefreq.tab0),"Total","Count") 
      colnames(simplefreq.tab)<-c(colnames(simplefreq.tab0),"All") 
    } 
    RESTABLE<-simplefreq.tab 
    CHI2<-simplefreq.chi 
    
    WARNING1 <- "Mining chi-square : what would be the chi-square if we were analyzing a 2*2 table. " 
    WARNING2 <- "Sign is used for over (+) or under(-) representation." 
    
    MINCHI2<-sign(simplefreq.chi$observed-simplefreq.chi$expected)* 
      (simplefreq.chi$observed-simplefreq.chi$expected)**2/
      simplefreq.chi$expected/ 
      ((1-prop.table(simplefreq.tab0)/prop.table(
        simplefreq.tab0,2))* 
         (1-prop.table(simplefreq.tab0)/prop.table(
           simplefreq.tab0,1))) 
    
    
    WARNING1 <- "Mining odds ratio : what would be the odds ratio if we were analyzing a 2*2 table. " 
    n11<-simplefreq.tab0 
    n12<-(n11/prop.table(n11,1))-n11 
    n21<-(n11/prop.table(n11,2))-n11 
    n22<-sum(n11)-n12-n21-n11 
    MINOR<-(n11/n12/(n21/n22)) 
  } 
  if (!is.null(html)) 
  { 
    library(R2HTML) 
    HTML(simplefreq.tab,html) 
  } 
  structure(list(simplefreq=RESTABLE,chi2=CHI2,
                 minor=MINOR,minnchi2=MINCHI2)) 
}

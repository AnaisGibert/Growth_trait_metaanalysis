to.pdf <- function(expr, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

figure_1 <- function(RawData) {
  p2 <- my_plot_1("",
          ggplot(RawData, aes(x=reorder(factor(trait), factor(trait), function(x) length(x)*1), fill=stage, order=stage)))
  p2
}

figure_allometry <- function() {

  download_baad <- function(destination_filename) {
    url <-
      "https://github.com/dfalster/baad/releases/download/v0.9.0/baad.rds"
    download(url, destination_filename, mode="wb")
  }

  single_plot <- function(px, py, data) {

    plot(data[[px$var]], data[[py$var]], log="xy",
      xlim = px$lim, ylim = py$lim,
      xlab = px$lab, ylab = py$lab, col=adjustcolor("#00000033", alpha=0.5))
    abline(v=px$sap, col="green", lty=5)
    abline(h=py$sap, col="green", lty=5)
    abline(v=px$adult, col="red", lty=5)
    abline(h=py$adult, col="red", lty=5)
  }

  filename <- "downloads/baad.rds"
  if(!file.exists(filename)) {
   dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
   download_baad(filename)
  }

  data <- readRDS(filename)[["data"]]

  pars <- list()
  pars[["dia"]] <- list(
  						var = "d.ba",
  						lab = "basal diameter (m)",
  						lim = c(0.0001, 1),
  						sap = 0.0025,
  						adult = 0.10)

  pars[["ht"]] <- list(
  						var = "h.t",
  						lab = "height (m)",
  						lim = c(0.01, 100),
  						sap = 0.5,
  						adult = 5)

  pars[["age"]] <- list(
  						var = "age",
  						lab = "age (yr)",
  						lim = c(0.1, 200),
  						sap = 1,
  						adult = 30)


  par(mfrow=c(1,3))
  single_plot(pars[["dia"]], pars[["age"]], data)
  single_plot(pars[["ht"]], pars[["age"]], data)
  single_plot(pars[["dia"]], pars[["ht"]], data)
}

####################################################################
#### Geographical distribution of site locations from studies  #####
###  performed in the filed used in the meta-analysis. (Fig A3)  ###
####################################################################

figure_map <- function(CoordTable) {
  mapWorld <- borders("world", colour="#FFCC33", fill="#FFCC33",lty=0) # create a layer of borders
  mp <- ggplot() + mapWorld + theme(text=element_text(size = 9, colour = "black"),
                                      title=element_text(size=9),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.background = element_rect(colour="#333333", fill="#6699FF"))

  coordinate.map <- mp  +
    geom_point(aes(x=CoordTable$long.dd,y=CoordTable$lat.dd),color='red',alpha=I(5/10))
  coordinate.map
}

figure_2 <- function(GC, GI) {
  CoefModel.SLA<-fun_model(GI[["SLA"]],GC[["SLA"]])
    CoefModel.SLA["trait"] <- "SLA"
    CoefModel.SLA.s <- subset(CoefModel.SLA,stress=='complete')
    CoefModel.SLA.opt <- subset(CoefModel.SLA,stress=='ideal')
    LRT.sla <- fun_OneLR(GI[["SLA"]])
    LRT.sla.2 <- fun_OneLR(GC[["SLA"]])
    PVAL.sla <- fun_Onepvalue(GI[["SLA"]])
    PVAL.sla.2 <- fun_Onepvalue(GC[["SLA"]])

  CoefModel.WD<-fun_model(GI[["WD"]],GC[["WD"]])
    CoefModel.WD["trait"] <- "WD"
    CoefModel.WD.s <- subset(CoefModel.WD,stress=='complete')
    CoefModel.WD.opt <- subset(CoefModel.WD,stress=='ideal')
    LRT.wd <- fun_OneLR(GI[["WD"]])
    LRT.wd.2 <- fun_OneLR(GC[["WD"]])
    PVAL.wd <- fun_Onepvalue(GI[["WD"]])
    PVAL.wd.2 <- fun_Onepvalue(GC[["WD"]])

  CoefModel.Hmax<-fun_model(GI[["Hmax"]],GC[["Hmax"]])
    CoefModel.Hmax["trait"] <- "Hmax"
    CoefModel.Hmax.s <- subset(CoefModel.Hmax,stress=='complete')
    CoefModel.Hmax.opt <- subset(CoefModel.Hmax,stress=='ideal')
    LRT.h <- fun_OneLR(GI[["Hmax"]])
    LRT.h.2 <- fun_OneLR(GC[["Hmax"]])
    PVAL.h <- fun_Onepvalue(GI[["Hmax"]])
    PVAL.h.2 <- fun_Onepvalue(GC[["Hmax"]])

  CoefModel.Seedmass<-fun_model(GI[["Seedmass"]],GC[["Seedmass"]])
    CoefModel.Seedmass["trait"] <- "Seedmass"
    CoefModel.Seedmass.s <- subset(CoefModel.Seedmass,stress=='complete')
    CoefModel.Seedmass.opt <- subset(CoefModel.Seedmass,stress=='ideal')
    LRT.sm <- fun_OneLR(GI[["Seedmass"]])
    LRT.sm.2 <- fun_OneLR(GC[["Seedmass"]])
    PVAL.sm <- fun_Onepvalue(GI[["Seedmass"]])
    PVAL.sm.2 <- fun_Onepvalue(GC[["Seedmass"]])

  CoefModel.Aarea<-fun_model(GI[["Aarea"]],GC[["Aarea"]])
    CoefModel.Aarea["trait"] <- "Aarea"
    CoefModel.Aarea.s <- subset(CoefModel.Aarea,stress=='complete')
    CoefModel.Aarea.opt <- subset(CoefModel.Aarea,stress=='ideal')
    LRT.a <- fun_OneLR(GI[["Aarea"]])
    LRT.a.2 <- fun_OneLR(GC[["Aarea"]])
    PVAL.a <- fun_Onepvalue(GI[["Aarea"]])
    PVAL.a.2 <- fun_Onepvalue(GC[["Aarea"]])

  p1 <- coeff.plot(data=CoefModel.SLA, data.complete=CoefModel.SLA.s, data.ideal=CoefModel.SLA.opt,
                   LRT=LRT.sla, PVAL= PVAL.sla, title="a) SLA",
                   significativite="***",round.value=4,
                   limit.x.min=-0.5,limit.x.max=1.5,
                   limit.x.n=1, vjust.value=1,
                   limit.x.text=-0.1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   color1="grey",color2="black")

  p2 <- coeff.plot(data=CoefModel.WD, data.complete=CoefModel.WD.s, data.ideal=CoefModel.WD.opt,
                   LRT=LRT.wd, PVAL= PVAL.wd, title="b) WD",
                   significativite="ns",round.value=2,
                   limit.x.min=-1,limit.x.max=0.5,
                   limit.x.n=0.3, vjust.value=1,
                   limit.x.text=-0.75,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   color1="grey",color2="black")

  p3 <- coeff.plot.ideal(data.ideal=CoefModel.Hmax.opt,
                   LRT=LRT.h, PVAL= PVAL.h, title="c) Hmax",
                   significativite="ns",round.value=2,
                   limit.x.min=-0.5,limit.x.max=1.5,
                   limit.x.text=-0.1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   limit.x.n=1.2, vjust.value=0,color1="black")

  p4 <- coeff.plot(data=CoefModel.Seedmass, data.complete=CoefModel.Seedmass.s, data.ideal=CoefModel.Seedmass.opt,
                   LRT=LRT.sm, PVAL= PVAL.sm, title="d) Seed mass",
                   significativite="***",round.value=3,
                   limit.x.min=-1.5,limit.x.max=1,
                   limit.x.text=-1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   limit.x.n=0.7, vjust.value=1,
                   color1="grey",color2="black")

  p5 <- coeff.plot(data=CoefModel.Aarea, data.complete=CoefModel.Aarea.s, data.ideal=CoefModel.Aarea.opt,
                   LRT=LRT.a, PVAL= PVAL.a, title="e) Aarea",
                   significativite="ns",round.value=3,
                   limit.x.min=-2,limit.x.max=3,
                   limit.x.text=-1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   limit.x.n=2.5, vjust.value=1,
                   color1="grey",color2="black") +
    scale_fill_discrete(name="",breaks=c("ideal", "complete"),
                        labels=c("ideal"="ideal", "complete"="complete"))+
                  theme (legend.title=element_blank(),
                         legend.justification=c(0,0),
                         legend.position=c(1.2,0.5),
                         legend.key = element_blank())

  p1 <- p1 + theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank())
  p2 <- p2 + theme(axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank())
  p3 <- p3 + theme(plot.margin=unit(c(0,0,1.5,0),"mm"),axis.title.x=element_blank())
  p4 <- p4 + theme(axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"))
  p5 <- p5 + theme(plot.margin=unit(c(0,0,0,0),"mm"))

  grid.arrange(p1,p2,p3,p4,p5,ncol=2, nrow=3,widths=c(1.2,1))
}

figure_8 <- function(GCi, GIi) {
  CoefModel.SLA<-fun_model3(GIi[["SLA"]],GCi[["SLA"]])
  CoefModel.SLA["trait"] <- "SLA"
  CoefModel.SLA.s <- subset(CoefModel.SLA,stress=='complete')
  CoefModel.SLA.opt <- subset(CoefModel.SLA,stress=='ideal')
  LRT.sla <- fun_OneLR2(GIi[["SLA"]])
  LRT.sla.2 <- fun_OneLR2(GCi[["SLA"]])
  PVAL.sla <- fun_Onepvalue2(GIi[["SLA"]])
  PVAL.sla.2 <- fun_Onepvalue2(GCi[["SLA"]])
  
  CoefModel.WD<-fun_model3(GIi[["WD"]],GCi[["WD"]])
  CoefModel.WD["trait"] <- "WD"
  CoefModel.WD.s <- subset(CoefModel.WD,stress=='complete')
  CoefModel.WD.opt <- subset(CoefModel.WD,stress=='ideal')
  LRT.wd <- fun_OneLR2(GIi[["WD"]])
  LRT.wd.2 <- fun_OneLR2(GCi[["WD"]])
  PVAL.wd <- fun_Onepvalue2(GIi[["WD"]])
  PVAL.wd.2 <- fun_Onepvalue2(GCi[["WD"]])
  
  CoefModel.Hmax<-fun_model3(GIi[["Hmax"]],GCi[["Hmax"]])
  CoefModel.Hmax["trait"] <- "Hmax"
  CoefModel.Hmax.s <- subset(CoefModel.Hmax,stress=='complete')
  CoefModel.Hmax.opt <- subset(CoefModel.Hmax,stress=='ideal')
  LRT.h <- fun_OneLR2(GIi[["Hmax"]])
  LRT.h.2 <- fun_OneLR2(GCi[["Hmax"]])
  PVAL.h <- fun_Onepvalue2(GIi[["Hmax"]])
  PVAL.h.2 <- fun_Onepvalue2(GCi[["Hmax"]])
  
  CoefModel.Seedmass<-fun_model3(GIi[["Seedmass"]],GCi[["Seedmass"]])
  CoefModel.Seedmass["trait"] <- "Seedmass"
  CoefModel.Seedmass.s <- subset(CoefModel.Seedmass,stress=='complete')
  CoefModel.Seedmass.opt <- subset(CoefModel.Seedmass,stress=='ideal')
  LRT.sm <- fun_OneLR2(GIi[["Seedmass"]])
  LRT.sm.2 <- fun_OneLR2(GCi[["Seedmass"]])
  PVAL.sm <- fun_Onepvalue2(GIi[["Seedmass"]])
  PVAL.sm.2 <- fun_Onepvalue2(GCi[["Seedmass"]])
  
  CoefModel.Aarea<-fun_model3(GIi[["Aarea"]],GCi[["Aarea"]])
  CoefModel.Aarea["trait"] <- "Aarea"
  CoefModel.Aarea.s <- subset(CoefModel.Aarea,stress=='complete')
  CoefModel.Aarea.opt <- subset(CoefModel.Aarea,stress=='ideal')
  LRT.a <- fun_OneLR2(GIi[["Aarea"]])
  LRT.a.2 <- fun_OneLR2(GCi[["Aarea"]])
  PVAL.a <- fun_Onepvalue2(GIi[["Aarea"]])
  PVAL.a.2 <- fun_Onepvalue2(GCi[["Aarea"]])
  
  p1 <- coeff.plot.2(data=CoefModel.SLA, data.complete=CoefModel.SLA.s, data.ideal=CoefModel.SLA.opt,
                   LRT=LRT.sla, PVAL= PVAL.sla, title="a) SLA",
                   significativite="***",round.value=4,
                   limit.x.min=-0.5,limit.x.max=1.5,
                   limit.x.n=1, vjust.value=1,
                   limit.x.text=-0.1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   color1="grey",color2="black")
  
  p2 <- coeff.plot.2(data=CoefModel.WD, data.complete=CoefModel.WD.s, data.ideal=CoefModel.WD.opt,
                   LRT=LRT.wd, PVAL= PVAL.wd, title="b) WD",
                   significativite="ns",round.value=2,
                   limit.x.min=-1,limit.x.max=0.5,
                   limit.x.n=0.3, vjust.value=1,
                   limit.x.text=-0.75,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   color1="grey",color2="black")
  
  p3 <- coeff.plot.ideal.2(data.ideal=CoefModel.Hmax.opt,
                         LRT=LRT.h, PVAL= PVAL.h, title="c) Hmax",
                         significativite="ns",round.value=2,
                         limit.x.min=-0.5,limit.x.max=1.5,
                         limit.x.text=-0.1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                         limit.x.n=1.2, vjust.value=0,color1="black")
  
  p4 <- coeff.plot.2(data=CoefModel.Seedmass, data.complete=CoefModel.Seedmass.s, data.ideal=CoefModel.Seedmass.opt,
                   LRT=LRT.sm, PVAL= PVAL.sm, title="d) Seed mass",
                   significativite="***",round.value=3,
                   limit.x.min=-1.5,limit.x.max=1,
                   limit.x.text=-1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   limit.x.n=0.7, vjust.value=1,
                   color1="grey",color2="black")
  
  p5 <- coeff.plot.2(data=CoefModel.Aarea, data.complete=CoefModel.Aarea.s, data.ideal=CoefModel.Aarea.opt,
                   LRT=LRT.a, PVAL= PVAL.a, title="e) Aarea",
                   significativite="ns",round.value=3,
                   limit.x.min=-2,limit.x.max=3,
                   limit.x.text=-1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   limit.x.n=2.5, vjust.value=1,
                   color1="grey",color2="black") +
    scale_fill_discrete(name="",breaks=c("ideal", "complete"),
                        labels=c("ideal"="ideal", "complete"="complete"))+
    theme (legend.title=element_blank(),
           legend.justification=c(0,0),
           legend.position=c(1.2,0.5),
           legend.key = element_blank())
  
  p1 <- p1 + theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank())
  p2 <- p2 + theme(axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank())
  p3 <- p3 + theme(plot.margin=unit(c(0,0,1.5,0),"mm"),axis.title.x=element_blank())
  p4 <- p4 + theme(axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"))
  p5 <- p5 + theme(plot.margin=unit(c(0,0,0,0),"mm"))
  
  grid.arrange(p1,p2,p3,p4,p5,ncol=2, nrow=3,widths=c(1.2,1))
}

figure_A1 <- function(CompleteData_inter) {

    s <- unique(CompleteData_inter[c("ref","id", "stageRGR","measure.size","size.mean.range","size.min","size.max","growth.form")])
    s <- s[order(s$measure.size,s$size.max,s$size.min,s$stageRGR),]
    s <- na.omit(s)

  # sort by trait et stage
  s$stageRGR <- factor(s$stageRGR, levels=c("seedling","juvenile","sapling","adult","mix"))
  s$id <- as.factor(s$id)
  sh <- subset(s,s$measure.size=="height")
  sa <- subset(s,s$measure.size=="age")
  sd <- subset(s,s$measure.size=="diameter")
  sd <- sd[order(sd$size.min,sd$size.max),]

  heights <- c(nrow(sa), nrow(sh), nrow(sd))
  heights <- heights/sum(heights) +0.1

 layout(matrix(c(1,2,3),3,1,byrow=TRUE), widths=c(3,3,3),heights= heights)
  # layout.show(3)
  par(mar=c(4,4,0,0.5), oma = c(1,2,1,1))

  plotCI2 <- function(data,  cuts, xlab) {

    rescale <- function(x, cuts) {
      x[x > max(cuts)] <- max(cuts)
      x2 <- x
      for( i in seq_along(cuts)[-1]) {
        low <- cuts[i-1]
        high <- cuts[i]
        ii <- (low < x & x <= high)
        x2[ii] <- i-2 + (x[ii] - low)/(high-low)
      }
      x2
    }
    
    data$size.min <- rescale( data$size.min, cuts)
    data$size.max <- rescale( data$size.max, cuts)

    data <- data[order(data$size.max, data$size.min),]

    n <- nrow(data)
    y <- rev(seq_len(n))
    cols <- c("green", "blue", "orange",
               "red", "grey")[data$stageRGR]

    plot(NA, xlim=c(0,3), ylim= c(0,n), yaxt="n",xaxt="n",
         xlab="",ylab="")
    mtext(xlab, 1, line=2, cex=0.75)
    segments(data$size.min, y, data$size.max, col=cols)
    i <- data$size.max == data$size.min
    points(data$size.min[i], y[i], col=cols[i], pch='-')
    axis(1, at=0:3, labels = cuts, las=1)
    axis(2, at=y, labels = data$id, las=1, cex.axis=0.5)
    abline(v=1, col="grey")
    abline(v=2, col="grey")
  }

  plotCI2(sa, c(0,1,5,10), "age (yrs)")
  abline(v = 1, col = "green")

  plotCI2(sh, c(0,0.5,2,20), "height (m)")
  abline(v = 1, col = "green")

  plotCI2(sd, c(0,1,10,80), "diameter (cm)")
  abline(v = 2, col = "red")

  mtext("Study ID", 2, line=0, outer=TRUE,  cex=0.75)
}


figure_A4 <- function(GC) {
  GC[["SLA"]]  <- GC[["SLA"]] [!is.na(GC[["SLA"]] [,"corr.r"]),]
  GC[["WD"]] <- GC[["WD"]] [!is.na(GC[["WD"]] [,"corr.r"]),]
  GC[["Hmax"]] <- GC[["Hmax"]][!is.na(GC[["Hmax"]][,"corr.r"]),]
  GC[["Seedmass"]] <- GC[["Seedmass"]][!is.na(GC[["Seedmass"]][,"corr.r"]),]
  GC[["Aarea"]] <- GC[["Aarea"]][!is.na(GC[["Aarea"]][,"corr.r"]),]

  p1 <- my_plot_corr.r(GC[["SLA"]],title="a) SLA")+ theme(plot.margin=unit(c(0,0,0,0),"mm"), legend.position="none")
  p2 <- my_plot_corr.r(GC[["WD"]], title="b)WD")+ theme(axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"), legend.position="none")
  p3 <- my_plot_corr.r(GC[["Hmax"]],title="c) Hmax")+ theme(plot.margin=unit(c(0,0,0,0),"mm"), legend.position="none")
  p4 <- my_plot_corr.r(GC[["Seedmass"]], title="d) Seed mass",xlab='Case studies ranked by coefficient of correlation r')+ theme(axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"), legend.position="none")
  p5 <- my_plot_corr.r(GC[["Aarea"]], title="e) Aarea",xlab='Case studies ranked by coefficient of correlation r') + theme(legend.title=element_blank(),
                               legend.justification=c(0,0),
                               legend.position=c(1.2,0.5),
                               legend.key = element_blank(),
                               plot.margin=unit(c(0,0,0,0),"mm"))
  grid.arrange(p1,p2,p3,p4,p5,ncol=2, nrow=3,widths=c(1.1,1))
}

figure_A5 <- function(GI) {

    table_trait <- function(trait) {
       rbind(
        table_overall(GI[[trait]]),
        table_overall.stage(GI[[trait]]))
    }
    SLA <- table_trait("SLA")
    WD <- table_trait("WD")
    Hmax <- table_trait("Hmax")
    Seedmass <- table_trait("Seedmass")
    Aarea <- table_trait("Aarea")

    p1 <- my_plot_overall(SLA, title='a) SLA') + theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank())
    p2 <- my_plot_overall(WD, title='b) WD') + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank())
    p3 <- my_plot_overall(Hmax, title='c) Hmax') + theme(plot.margin=unit(c(0,0,1.5,0),"mm"),axis.title.x=element_blank())
    p4 <- my_plot_overall(Seedmass, title='d) Seed mass') + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"))
    p5 <- my_plot_overall(Aarea, title='e) Aarea') + theme (legend.title=element_blank(), legend.justification=c(0,0),legend.position=c(1.2,0.5), legend.key = element_blank(), plot.margin=unit(c(0,0,0,0),"mm"))

    grid.arrange(p1,p2,p3,p4,p5,ncol=2, nrow=3,widths=c(1.2,1))
}

figure_A5b <- function(GIrgr, GIagr) {
  
  table_trait <- function(GIrgr, GIagr, trait) {
      x <-  rbind(
          table_overall(GIrgr[[trait]]),
          table_overall.stage(GIrgr[[trait]]),
          table_overall(GIagr[[trait]]),
          table_overall.stage(GIagr[[trait]]))
      x["growth"] <- NA
      x$growth <- c("RGR","RGR","RGR","RGR","AbsGR","AbsGR", "AbsGR","AbsGR")
      for(f in c("growth","stage")) {
        x[[f]] <- as.factor(x[[f]])
      }
      x
   }
  
  SLA <- table_trait(GIrgr, GIagr,"SLA")
  WD <- table_trait(GIrgr, GIagr,"WD")
  Hmax <- table_trait(GIrgr, GIagr,"Hmax")
  Seedmass <- table_trait(GIrgr, GIagr,"Seedmass")
  Aarea <- table_trait(GIrgr, GIagr,"Aarea")
  
  
  p1 <- my_plot_overall.b(SLA, title='a) SLA') + theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank())
  p2 <- my_plot_overall.b(WD, title='b) WD') + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank())
  p3 <- my_plot_overall.b(Hmax, title='c) Hmax') + theme(plot.margin=unit(c(0,0,1.5,0),"mm"),axis.title.x=element_blank())
  p4 <- my_plot_overall.b(Seedmass, title='d) Seed mass') + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"))
  p5 <- my_plot_overall.b(Aarea, title='e) Aarea') + theme (legend.title=element_blank(), legend.justification=c(0,0),legend.position=c(1.2,0.5), legend.key = element_blank(), plot.margin=unit(c(0,0,0,0),"mm"))
  
  grid.arrange(p1,p2,p3,p4,p5,ncol=2, nrow=3,widths=c(1.2,1))
}

figure_A6 <- function(GC) {
  CoefModel.SLA<-fun_model2(GC[["SLA"]])
  CoefModel.SLA["trait"] <- "SLA"
  LRT.sla <- fun_OneLR(GC[["SLA"]])
  PVAL.sla <- fun_Onepvalue(GC[["SLA"]])

  CoefModel.WD<-fun_model2(GC[["WD"]])
  CoefModel.WD["trait"] <- "WD"
  LRT.wd <- fun_OneLR(GC[["WD"]])
  PVAL.wd <- fun_Onepvalue(GC[["WD"]])

  
  CoefModel.Hmax<-fun_model2(GC[["Hmax"]])
  CoefModel.Hmax["trait"] <- "Hmax"
  LRT.h <- fun_OneLR(GC[["Hmax"]])
  PVAL.h <- fun_Onepvalue(GC[["Hmax"]])

  
  CoefModel.Seedmass<-fun_model2(GC[["Seedmass"]])
  CoefModel.Seedmass["trait"] <- "Seedmass"
  LRT.sm <- fun_OneLR(GC[["Seedmass"]])
  PVAL.sm <- fun_Onepvalue(GC[["Seedmass"]])

  
  CoefModel.Aarea<-fun_model2(GC[["Aarea"]])
  CoefModel.Aarea["trait"] <- "Aarea"
  LRT.a <- fun_OneLR(GC[["Aarea"]])
  PVAL.a <- fun_Onepvalue(GC[["Aarea"]])
  
  p1 <- coeff.plot.ideal(data=CoefModel.SLA,
                   LRT=LRT.sla, PVAL= PVAL.sla, title="a) SLA",
                   significativite="***",round.value=4,
                   limit.x.min=-0.5,limit.x.max=1.5,
                   limit.x.text=-0.1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   limit.x.n=1.2, vjust.value=0,color1="black")
  

  p2 <- coeff.plot.ideal(data=CoefModel.WD, 
                   LRT=LRT.wd, PVAL= PVAL.wd, title="b) WD",
                   significativite="ns",round.value=2,
                   limit.x.min=-1,limit.x.max=0.5,
                   limit.x.text=-0.75,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   limit.x.n=0.3, vjust.value=0,color1="black")
 
  p3 <- coeff.plot.ideal(data.ideal=CoefModel.Hmax,
                         LRT=LRT.h, PVAL= PVAL.h, title="c) Hmax",
                         significativite="ns",round.value=2,
                         limit.x.min=-0.5,limit.x.max=1.5,
                         limit.x.text=-0.1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                         limit.x.n=1.2, vjust.value=0,color1="black")
  
  p4 <- coeff.plot.ideal(data=CoefModel.Seedmass,
                   LRT=LRT.sm, PVAL= PVAL.sm, title="d) Seed mass",
                   significativite="***",round.value=4,
                   limit.x.min=-1.5,limit.x.max=1,
                   limit.x.text=-1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   limit.x.n=0.7, vjust.value=0,color1="black")
  
  
  p5 <- coeff.plot.ideal(data=CoefModel.Aarea,
                   LRT=LRT.a, PVAL= PVAL.a, title="e) Aarea",
                   significativite="ns",round.value=3,
                   limit.x.min=-2,limit.x.max=3,
                   limit.x.text=-1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   limit.x.n=2.5, vjust.value=0,color1="black")+
    scale_fill_discrete(name="",breaks=c("raw"),
                        labels=c("raw"))+
    theme (legend.title=element_blank(),
           legend.justification=c(0,0),
           legend.position=c(1.2,0.5),
           legend.key = element_blank())
    
  
  p1 <- p1 + theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank())
  p2 <- p2 + theme(axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank())
  p3 <- p3 + theme(plot.margin=unit(c(0,0,1.5,0),"mm"),axis.title.x=element_blank())
  p4 <- p4 + theme(axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"))
  p5 <- p5 + theme(plot.margin=unit(c(0,0,0,0),"mm"))
  
  grid.arrange(p1,p2,p3,p4,p5,ncol=2, nrow=3,widths=c(1.2,1))
}

figure_A7 <- function(RIi, RCi) {
  CoefModel.SLA<-fun_model(RIi[["SLA"]],RCi[["SLA"]])
  CoefModel.SLA.s <- subset(CoefModel.SLA,stress=='complete')
  CoefModel.SLA.opt <- subset(CoefModel.SLA,stress=='ideal')

  LRT.sla <- fun_OneLR(RIi[["SLA"]])
  LRT.sla.2 <- fun_OneLR(RCi[["SLA"]])
  PVAL.sla <- fun_Onepvalue(RIi[["SLA"]])
  PVAL.sla.2 <- fun_Onepvalue(RCi[["SLA"]])

  p1 <- coeff.plot(data=CoefModel.SLA, data.complete=CoefModel.SLA.s,data.ideal=CoefModel.SLA.opt,
                   LRT=LRT.sla, PVAL= PVAL.sla, title="a) SLA",
                   significativite="***", round.value=3,
                   limit.x.min=-0.5,limit.x.max=1.5,
                   limit.x.n=1.2, vjust.value=0.8,
                   limit.x.text=-0.15,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   color1="grey",color2="black")

  CoefModel.WD<-fun_model(RIi[["WD"]],RCi[["WD"]])
  CoefModel.WD.s <- subset(CoefModel.WD,stress=='complete')
  CoefModel.WD.opt <- subset(CoefModel.WD,stress=='ideal')

  LRT.wd <- fun_OneLR(RIi[["WD"]])
  LRT.wd.2 <- fun_OneLR(RCi[["WD"]])
  PVAL.wd <- fun_Onepvalue(RIi[["WD"]])
  PVAL.wd.2 <- fun_Onepvalue(RCi[["WD"]])

  p2 <- coeff.plot(data=CoefModel.WD, data.complete=CoefModel.WD.s, data.ideal=CoefModel.WD.opt,
                   LRT=LRT.wd, PVAL= PVAL.wd, title="b) WD",
                   significativite="ns",round.value=3,
                   limit.x.min=-1,limit.x.max=0.5,
                   limit.x.n=0.3, vjust.value=0.8,
                   limit.x.text=-0.75,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   color1="grey",color2="black")


  CoefModel.Hmax<-fun_model(RIi[["Hmax"]],RCi[["Hmax"]])
  CoefModel.Hmax.s <- subset(CoefModel.Hmax,stress=='complete')
  CoefModel.Hmax.opt <- subset(CoefModel.Hmax,stress=='ideal')

  LRT.h <- fun_OneLR(RIi[["Hmax"]])
  LRT.h.2 <- fun_OneLR(RCi[["Hmax"]])
  PVAL.h <- fun_Onepvalue(RIi[["Hmax"]])
  PVAL.h.2 <- fun_Onepvalue(RCi[["Hmax"]])

  p3 <- coeff.plot.ideal(data.ideal=CoefModel.Hmax.opt,
                     LRT=0.0, PVAL= 0.98, title="c) Hmax",
                     significativite="ns",round.value=3,
                     limit.x.min=-0.5,limit.x.max=1.5,
                     limit.x.text=-0.25,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                     limit.x.n=1.3, vjust.value=0,color1="black")

  CoefModel.Seedmass<-fun_model(RIi[["Seedmass"]],RCi[["Seedmass"]])
  CoefModel.Seedmass.s <- subset(CoefModel.Seedmass,stress=='complete')
  CoefModel.Seedmass.opt <- subset(CoefModel.Seedmass,stress=='ideal')

  LRT.sm <- fun_OneLR(RIi[["Seedmass"]])
  LRT.sm.2 <- fun_OneLR(RCi[["Seedmass"]])
  PVAL.sm <- fun_Onepvalue(RIi[["Seedmass"]])
  PVAL.sm.2 <- fun_Onepvalue(RCi[["Seedmass"]])

  p4 <- coeff.plot(data=CoefModel.Seedmass, data.complete=CoefModel.Seedmass.s, data.ideal=CoefModel.Seedmass.opt,
                   LRT=LRT.sm, PVAL= PVAL.sm, title="d) Seed mass",
                   significativite="***",round.value=3,
                   limit.x.min=-1.5,limit.x.max=1,
                   limit.x.text=-1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   limit.x.n=0.7, vjust.value=0.8,
                   color1="grey",color2="black")


  CoefModel.Aarea<-fun_model(RIi[["Aarea"]],RCi[["Aarea"]])
  CoefModel.Aarea.s <- subset(CoefModel.Aarea,stress=='complete')
  CoefModel.Aarea.opt <- subset(CoefModel.Aarea,stress=='ideal')

  LRT.a <- fun_OneLR(RIi[["Aarea"]])
  LRT.a.2 <- fun_OneLR(RCi[["Aarea"]])
  PVAL.a <- fun_Onepvalue(RIi[["Aarea"]])
  PVAL.a.2 <- fun_Onepvalue(RCi[["Aarea"]])

  p5 <- coeff.plot(data=CoefModel.Aarea, data.complete=CoefModel.Aarea.s, data.ideal=CoefModel.Aarea.opt,
                   LRT=LRT.a, PVAL= PVAL.a, title="e) Aarea",
                   significativite="ns",round.value=3,
                   limit.x.min=-2,limit.x.max=3,
                   limit.x.text=-1,limit.y.text.l1=0.5, limit.y.text.l2=0.25,
                   limit.x.n=2.5, vjust.value=1,
                   color1="grey",color2="black") +
                    scale_fill_discrete(name="",breaks=c("complete", "ideal"),
                                        labels=c("complete"="complete", "ideal"="ideal"))+
                    theme (legend.title=element_blank(),
                           legend.justification=c(0,0),
                           legend.position=c(1.2,0.5),
                           legend.key = element_blank())


  p1 <- p1 + theme(plot.margin=unit(c(0,0,0,0),"mm"), axis.title.x=element_blank())
  p2 <- p2 + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"))
  p3 <- p3 + theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank())
  p4 <- p4 + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"))
  p5 <- p5 + theme(plot.margin=unit(c(0,0,0,0),"mm"))

 grid.arrange(p1,p2,p3,p4,p5,ncol=2, nrow=3,widths=c(1.2,1))
}

# # figure_A8 <- function(GI, GC, trait, titles) {
#   par(mfcol=c(1,2))
#   par(mar=c(2,5,2,0))
# 
#   coeff.plot.multiple(GI[[trait]], params=rev(c("stageRGRseedling","stageRGRjuvenile","stageRGRsapling","stageRGRadult","stageRGRmix",
#                                                 "veg.typeboreal and temperate deciduous forest","veg.typeboreal forest","veg.typemediteranean","veg.typemix" ,
#                                                 "veg.typetemperate deciduous and mediteranean forest","veg.typetemperate deciduous forest" , "veg.typetemperate rain forest" ,
#                                                 "veg.typetropical rain forest"  , "veg.typetropical seasonal forest"   ,"RGRGR(Di)", "RGRGR(Hi)", "RGRGR(Mi)"   , "RGRRGR(CSAi)"  ,
#                                                 "RGRRGR(Di)" ,"RGRRGR(Hi)" ,  "RGRRGR(Mi)", "RGRRGR(Vi)" ,"experimentcontrol" , "experimentdatabase", "experimentfield" ,"experimentnature")),
#                       labels=rev(c('seedling','juvenile','sapling','adult','all stage','boreal&temp','boreal','med','across veg','temp&med','temp','temp rain','trop rain','trop seas',
#                                    'GR(D)','GR(H)','GR(M)','RGR(CSA)','RGR(D)','RGR(H)','RGR(M)','RGR(V)','greenhouse','database','field exp','forest')),
#                       title=paste0(titles[1], ") ", trait, "- ideal dataset"))
# 
#   mtext("mod4", side=2, line=4.2, cex=0.8, at=2.2)
#   mtext("mod3", side=2, line=4.2, cex=0.8, at=8.8)
#   mtext("mod2", side=2, line=4.2, cex=0.8, at=17)
#   mtext("mod1", side=2, line=4.2, cex=0.8, at=24)
# 
#   par(mar=c(2,1.5,2,3.5))
#   coeff.plot.multiple(GC[[trait]], params=rev(c("stageRGRseedling","stageRGRjuvenile","stageRGRsapling","stageRGRadult","stageRGRmix",
#                                                 "veg.typeboreal and temperate deciduous forest","veg.typeboreal forest","veg.typemediteranean","veg.typemix" ,
#                                                 "veg.typetemperate deciduous and mediteranean forest","veg.typetemperate deciduous forest" , "veg.typetemperate rain forest" ,
#                                                 "veg.typetropical rain forest"  , "veg.typetropical seasonal forest"   ,"RGRGR(Di)", "RGRGR(Hi)", "RGRGR(Mi)"   , "RGRRGR(CSAi)"  ,
#                                                 "RGRRGR(Di)" ,"RGRRGR(Hi)" ,  "RGRRGR(Mi)", "RGRRGR(Vi)" ,"experimentcontrol" , "experimentdatabase", "experimentfield" ,
#                                                 "experimentnature")),
#                       title=paste0(titles[2], ") ", trait, "- raw dataset"))
# }

figure_A8 <- function(GC, trait1, trait2, titles) {
  par(mfcol=c(1,2))
  par(mar=c(2,5,2,0))
  
  coeff.plot.multiple3(GC[[trait1]], params=rev(c("stageRGRseedling","stageRGRjuvenile","stageRGRsapling","stageRGRadult","stageRGRmix",
                                                "RGRGR(Di)", "RGRGR(Hi)", "RGRGR(Mi)"   , "RGRRGR(CSAi)"  ,
                                                "RGRRGR(Di)" ,"RGRRGR(Hi)" ,  "RGRRGR(Mi)", "RGRRGR(Vi)" ,"growth.formtree","growth.formwoody","growth.formacross growth form","experimentcontrol" , "experimentdatabase", "experimentfield" ,"experimentnature")),
                      labels=rev(c('seedling','juvenile','sapling','adult','all stage',
                                   'GR(D)','GR(H)','GR(M)','RGR(CSA)','RGR(D)','RGR(H)','RGR(M)','RGR(V)',"tree","woody","across GF",'greenhouse','database','field exp','forest')),
                      title=paste0(titles[1], ") ", trait1, "- ideal dataset"))
  
  mtext("mod4", side=2, line=4.2, cex=0.8, at=2.2)
  mtext("mod3", side=2, line=4.2, cex=0.8, at=6)
  mtext("mod2", side=2, line=4.2, cex=0.8, at=12)
  mtext("mod1", side=2, line=4.2, cex=0.8, at=18)
  
  par(mar=c(2,1.5,2,3.5))
  coeff.plot.multiple3(GC[[trait2]], params=rev(c("stageRGRseedling","stageRGRjuvenile","stageRGRsapling","stageRGRadult","stageRGRmix",
                                                  "RGRGR(Di)", "RGRGR(Hi)", "RGRGR(Mi)"   , "RGRRGR(CSAi)"  ,
                                                  "RGRRGR(Di)" ,"RGRRGR(Hi)" ,  "RGRRGR(Mi)", "RGRRGR(Vi)" ,"growth.formtree","growth.formwoody","growth.formacross growth form","experimentcontrol" , "experimentdatabase", "experimentfield" ,"experimentnature")),
                      title=paste0(titles[2], ") ", trait2, "- ideal dataset"))
}

figure_A82 <- function(GC, trait1, titles) {
  par(mfcol=c(1,2))
  par(mar=c(2,5,2,0))
  
  coeff.plot.multiple3.1(GC[[trait1]], params=rev(c("stageRGRseedling","stageRGRjuvenile","stageRGRsapling","stageRGRadult","stageRGRmix",
                                                  "RGRGR(Di)", "RGRGR(Hi)", "RGRGR(Mi)"   , "RGRRGR(CSAi)"  ,
                                                  "RGRRGR(Di)" ,"RGRRGR(Hi)" ,  "RGRRGR(Mi)", "RGRRGR(Vi)" ,"experimentcontrol" , "experimentdatabase", "experimentfield" ,"experimentnature")),
                       labels=rev(c('seedling','juvenile','sapling','adult','all stage',
                                    'GR(D)','GR(H)','GR(M)','RGR(CSA)','RGR(D)','RGR(H)','RGR(M)','RGR(V)','greenhouse','database','field exp','forest')),
                       title=paste0(titles[1], ") ", trait1, "- ideal dataset"))
  
  mtext("mod4", side=2, line=4.2, cex=0.8, at=2.2)
  mtext("mod2", side=2, line=4.2, cex=0.8, at=8.5)
  mtext("mod1", side=2, line=4.2, cex=0.8, at=15.2)
}


figure_A9 <- function(GCi){
  par(mfrow=c(3,2))
  par(mar=c(5,4,1,2))
  
  p1 <- figure_trim.and.fill(GIi[["SLA"]], title= 'a) SLA') 
  p2 <- figure_trim.and.fill(GIi[["WD"]], title= 'b) WD')
  p3 <- figure_trim.and.fill(GIi[["Hmax"]],title= 'c) Hmax') 
  p4 <- figure_trim.and.fill(GIi[["Seedmass"]],title= 'd) Seedmass') 
  p5 <- figure_trim.and.fill(GIi[["Aarea"]],title= 'e) Aarea') 
}


# figure_A8b <- function(GC, trait1, trait2, titles) {
#   par(mfcol=c(1,2))
#   par(mar=c(2,5,2,0))
#   
#   coeff.plot.multiple(GC[[trait1]], params=rev(c("stageRGRseedling","stageRGRjuvenile","stageRGRsapling","stageRGRadult","stageRGRmix",
#                                                 "veg.typeboreal and temperate deciduous forest","veg.typeboreal forest","veg.typemediteranean","veg.typemix" ,
#                                                 "veg.typetemperate deciduous and mediteranean forest","veg.typetemperate deciduous forest" , "veg.typetemperate rain forest" ,
#                                                 "veg.typetropical rain forest"  , "veg.typetropical seasonal forest"   ,"RGRGR(Di)", "RGRGR(Hi)", "RGRGR(Mi)"   , "RGRRGR(CSAi)"  ,
#                                                 "RGRRGR(Di)" ,"RGRRGR(Hi)" ,  "RGRRGR(Mi)", "RGRRGR(Vi)" ,"experimentcontrol" , "experimentdatabase", "experimentfield" ,"experimentnature")),
#                       labels=rev(c('seedling','juvenile','sapling','adult','all stage','boreal&temp','boreal','med','across veg','temp&med','temp','temp rain','trop rain','trop seas',
#                                    'GR(D)','GR(H)','GR(M)','RGR(CSA)','RGR(D)','RGR(H)','RGR(M)','RGR(V)','control','database','field exp','forest')),
#                       title=paste0(titles[1], ") ", trait1))
#   
#   mtext("mod4", side=2, line=4.2, cex=0.8, at=2.2)
#   mtext("mod3", side=2, line=4.2, cex=0.8, at=8.8)
#   mtext("mod2", side=2, line=4.2, cex=0.8, at=17)
#   mtext("mod1", side=2, line=4.2, cex=0.8, at=24)
#   
#   par(mar=c(2,1.5,2,3.5))
#   coeff.plot.multiple(GC[[trait2]], params=rev(c("stageRGRseedling","stageRGRjuvenile","stageRGRsapling","stageRGRadult","stageRGRmix",
#                                                 "veg.typeboreal and temperate deciduous forest","veg.typeboreal forest","veg.typemediteranean","veg.typemix" ,
#                                                 "veg.typetemperate deciduous and mediteranean forest","veg.typetemperate deciduous forest" , "veg.typetemperate rain forest" ,
#                                                 "veg.typetropical rain forest"  , "veg.typetropical seasonal forest"   ,"RGRGR(Di)", "RGRGR(Hi)", "RGRGR(Mi)"   , "RGRRGR(CSAi)"  ,
#                                                 "RGRRGR(Di)" ,"RGRRGR(Hi)" ,  "RGRRGR(Mi)", "RGRRGR(Vi)" ,"experimentcontrol" , "experimentdatabase", "experimentfield" ,
#                                                 "experimentnature")),
#                       title=paste0(titles[2], ") ", trait2))
# }
# 
# figure_A8c <- function(GC, trait1, titles) {
#   par(mfcol=c(1,2))
#   par(mar=c(2,5,2,0))
#   
#   coeff.plot.multiple(GC[[trait1]], params=rev(c("stageRGRseedling","stageRGRjuvenile","stageRGRsapling","stageRGRadult","stageRGRmix",
#                                                  "veg.typeboreal and temperate deciduous forest","veg.typeboreal forest","veg.typemediteranean","veg.typemix" ,
#                                                  "veg.typetemperate deciduous and mediteranean forest","veg.typetemperate deciduous forest" , "veg.typetemperate rain forest" ,
#                                                  "veg.typetropical rain forest"  , "veg.typetropical seasonal forest"   ,"RGRGR(Di)", "RGRGR(Hi)", "RGRGR(Mi)"   , "RGRRGR(CSAi)"  ,
#                                                  "RGRRGR(Di)" ,"RGRRGR(Hi)" ,  "RGRRGR(Mi)", "RGRRGR(Vi)" ,"experimentcontrol" , "experimentdatabase", "experimentfield" ,"experimentnature")),
#                       labels=rev(c('seedling','juvenile','sapling','adult','all stage','boreal&temp','boreal','med','across veg','temp&med','temp','temp rain','trop rain','trop seas',
#                                    'GR(D)','GR(H)','GR(M)','RGR(CSA)','RGR(D)','RGR(H)','RGR(M)','RGR(V)','control','database','field exp','forest')),
#                       title=paste0(titles[1], ") ", trait1))
#   
#   mtext("mod4", side=2, line=4.2, cex=0.8, at=2.2)
#   mtext("mod3", side=2, line=4.2, cex=0.8, at=8.8)
#   mtext("mod2", side=2, line=4.2, cex=0.8, at=17)
#   mtext("mod1", side=2, line=4.2, cex=0.8, at=24)
#   
# }

# figure_A9 <- function(GC) {
# 
#   funnel_SLA_nbsp <- my_funnelplot("a) SLA",
#     ggplot(GC[["SLA"]],aes(x=nb.sp,y=corr.r, colour=factor(stage), size=2, alpha=0.6))) +
#     geom_point() + scale_y_continuous("Correlation coefficient  r",limits=c(-1,1)) +theme( legend.position="none")
# 
# 
#   funnel_WD_nbsp <- my_funnelplot("b) WD",
#     ggplot(GC[["WD"]],aes(x=nb.sp,y=corr.r,colour=factor(stage), size=2, alpha=0.6))) +
#     geom_point()+
#     scale_y_continuous("",limits=c(-1,1)) +
#     scale_x_continuous("",limits=c(0,150)) +
#     theme( legend.position="none")
# 
#   funnel_Hmax_nbsp <- my_funnelplot("c) Hmax",
#     ggplot(GC[["Hmax"]],aes(x=nb.sp,y=corr.r,colour=factor(stage), size=2, alpha=0.6))) +
#     geom_point() + scale_y_continuous("Correlation coefficient r",limits=c(-1,1))+
#     scale_x_continuous("",limits=c(0,150)) +
#     theme( legend.position="none")
# 
# 
#   funnel_Seedmass_nbsp <- my_funnelplot("d) Seed mass",
#     ggplot(GC[["Seedmass"]],aes(x=nb.sp,y=corr.r,colour=factor(stage), size=2, alpha=0.6))) +
#     geom_point() + scale_x_continuous("Number of species",limits=c(0,150))+
#     scale_y_continuous("",limits=c(-1,1)) +
#     theme( legend.position="none")
# 
# 
#   funnel_Aarea_nbsp <- my_funnelplot("e) Aarea",
#     ggplot(GC[["Aarea"]],aes(x=nb.sp,y=corr.r,colour=factor(stage),size=2, alpha=0.6))) +
#     geom_point() + scale_alpha(guide = 'none')+ scale_size(guide = 'none') + scale_y_continuous("Correlation coefficient  r",limits=c(-1,1)) + scale_x_continuous("Number of species",limits=c(0,150))+
#     theme (legend.title=element_blank(), legend.justification=c(0,0),legend.position=c(1.2,0.5), legend.key = element_blank(), plot.margin=unit(c(0,0,0,0),"mm"))
# 
# 
#   grid.arrange(funnel_SLA_nbsp,funnel_WD_nbsp,funnel_Hmax_nbsp, funnel_Seedmass_nbsp ,funnel_Aarea_nbsp, nrow=3, ncol=2)
# }

figure_A10a <- function(GC) {
  p1 <- my_plot_3("a) SLA",
          ggplot(GC[["SLA"]],aes(x=reorder(factor(RGR),factor(stage),function(x) length(x)*1),fill=stage,order=stage)))

  p2 <- my_plot_3("b) WD",
          ggplot(GC[["WD"]],aes(x=reorder(factor(RGR),factor(stage),function(x) length(x)*1),fill=stage,order=stage)))

  p3 <- my_plot_3("c) Hmax",
          ggplot(GC[["Hmax"]],aes(x=reorder(factor(RGR),factor(stage),function(x) length(x)*1),fill=stage,order=stage)))

  p4 <- my_plot_3("d) Seed mass",
          ggplot(GC[["Seedmass"]],aes(x=reorder(factor(RGR),factor(stage),function(x) length(x)*1),fill=stage,order=stage)))

  p5 <- my_plot_3("e) Aarea",
          ggplot(GC[["Aarea"]],aes(x=reorder(factor(RGR),factor(stage),function(x) length(x)*1),fill=stage,order=stage)))

  p1 <- p1 + theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),
                   legend.position="none")
  p2 <- p2 + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),
                    legend.position="none")
  p3 <- p3 + theme(plot.margin=unit(c(0,0,1.5,0),"mm"),axis.title.x=element_blank(),
                   legend.position="none")
  p4 <- p4 + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"),
                    legend.position="none")
  p5 <- p5 + theme(plot.margin=unit(c(0,0,0,0),"mm"), legend.justification=c(0,0),legend.position=c(1.2,0.2), legend.key = element_blank())

  grid.arrange(p1,p2,p3,p4,p5,ncol=2, nrow=3,widths=c(1.2,1))
}

figure_A10b <- function(GC) {
  p1 <- my_plot_3b("a) SLA",
                  ggplot(GC[["SLA"]],aes(x=reorder(factor(stage),factor(growth),function(x) length(x)*1), fill=growth)))
  
  p2 <- my_plot_3b("b) WD",
                  ggplot(GC[["WD"]],aes(x=reorder(factor(stage),factor(growth),function(x) length(x)*1),fill=growth, order=stage)))
  
  p3 <- my_plot_3b("c) Hmax",
                  ggplot(GC[["Hmax"]],aes(x=reorder(factor(stage),factor(growth),function(x) length(x)*1),fill=growth, order=stage)))
  
  p4 <- my_plot_3b("d) Seed mass",
                  ggplot(GC[["Seedmass"]],aes(x=reorder(factor(stage),factor(growth),function(x) length(x)*1),fill=growth, order=stage)))
  
  p5 <- my_plot_3b("e) Aarea",
                  ggplot(GC[["Aarea"]],aes(x=reorder(factor(stage),factor(growth),function(x) length(x)*1),fill=growth, order=stage)))
  
  p1 <- p1 + theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),
                   legend.position="none")
  p2 <- p2 + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),
                    legend.position="none")
  p3 <- p3 + theme(plot.margin=unit(c(0,0,1.5,0),"mm"),axis.title.x=element_blank(),
                   legend.position="none")
  p4 <- p4 + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"),
                    legend.position="none")
  p5 <- p5 + theme(plot.margin=unit(c(0,0,0,0),"mm"), legend.justification=c(0,0),legend.position=c(1.2,0.2), legend.key = element_blank())
  
  grid.arrange(p1,p2,p3,p4,p5,ncol=2, nrow=3,widths=c(1.2,1))
}

figure_A11 <- function(GI, GC, trait, titles) {
  par(mfcol=c(1,2))
  par(mar=c(2,5,2,0))
  
  coeff.plot.multiple2(GI[[trait]], params=rev(c("stagejuvenile","stagesapling","stageadult",
                                                "growthAbGR","growthRGR",
                                                "measurementDiameter","measurementHeight" , "measurementMass" ,"measurementOther")),
                      labels=rev(c('juvenile','sapling','adult','AbGR','RGR','diameter','height','mass','other')),
                      title=paste0(titles[1], ") ", trait, "- ideal dataset"))
  
  mtext("mod3", side=2, line=4.2, cex=0.8, at=2.5)
  mtext("mod2", side=2, line=4.2, cex=0.8, at=6.5)
  mtext("mod1", side=2, line=4.2, cex=0.8, at=8)
  
  par(mar=c(2,1.5,2,3.5))
  coeff.plot.multiple2(GC[[trait]], params=rev(c("stagejuvenile","stagesapling","stageadult",
                                                 "growthAbGR","growthRGR",
                                                 "measurementDiameter","measurementHeight" , "measurementMass" ,"measurementOther")),
                      title=paste0(titles[2], ") ", trait, "- raw dataset"))
}


# figure_A12 <- function(GC) {
#   p1 <- my_plot_3c("a) SLA",
#                   ggplot(GC[["SLA"]],aes(x=(year,factor(stage),function(x) length(x)*1),fill=stage,order=stage)))
#   
#   p2 <- my_plot_3c("b) WD",
#                   ggplot(GC[["WD"]],aes(x=reorder(year, factor(stage),function(x) length(x)*1),fill=stage,order=stage)))
#   
#   p3 <- my_plot_3c("c) Hmax",
#                   ggplot(GC[["Hmax"]],aes(x=reorder(year,factor(stage),function(x) length(x)*1),fill=stage,order=stage)))
#   
#   p4 <- my_plot_3c("d) Seed mass",
#                   ggplot(GC[["Seedmass"]],aes(x=reorder(year,factor(stage),function(x) length(x)*1),fill=stage,order=stage)))
#   
#   p5 <- my_plot_3c("e) Aarea",
#                   ggplot(GC[["Aarea"]],aes(x=reorder(year,factor(stage),function(x) length(x)*1),fill=stage,order=stage)))
#   
#   p1 <- p1 + theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),
#                    legend.position="none")
#   p2 <- p2 + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),
#                     legend.position="none")
#   p3 <- p3 + theme(plot.margin=unit(c(0,0,1.5,0),"mm"),axis.title.x=element_blank(),
#                    legend.position="none")
#   p4 <- p4 + theme( axis.text.y = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,0,0,0),"mm"),
#                     legend.position="none")
#   p5 <- p5 + theme(plot.margin=unit(c(0,0,0,0),"mm"), legend.justification=c(0,0),legend.position=c(1.2,0.2), legend.key = element_blank())
#   
#   grid.arrange(p1,p2,p3,p4,p5,ncol=2, nrow=3,widths=c(1.2,1))
# }

figure_A12 <- function(GIi) {
  
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
  
  funnel_SLA_year <- my_funnelplot("a) SLA",
                                    ggplot(GIi[["SLA"]],aes(x=year,y=corr.r, colour=factor(stage), size=2, alpha=0.6))) +
                                    geom_point() +
                                    scale_y_continuous("Correlation coefficient  r",limits=c(-1,1)) +theme( legend.position="none")+ 
                                    scale_x_continuous("Years of publication",limits=c(1990,2015)) +
                                    annotate("text", x = 1995, y = -0.9, label = paste("LRT:", round(LRT_SLA,0)), size = 2) +
                                    annotate("text", x = 1995, y = -1, label = paste("p.value =", round(PVAL_SLA, 3), "***"), size = 2)
                                  
  
  funnel_WD_year  <- my_funnelplot("b) WD",
                                  ggplot(GIi[["WD"]],aes(x=year,y=corr.r,colour=factor(stage), size=2, alpha=0.6))) +
                                  geom_point()+
                                  scale_y_continuous("",limits=c(-1,1)) +
                                  scale_x_continuous("",limits=c(1990,2015)) +
                                  theme( legend.position="none")+
                                  annotate("text", x = 1995, y = -0.9, label = paste("LRT:", round(LRT_WD,0)), size = 2) +
                                  annotate("text", x = 1995, y = -1, label = paste("p.value =", round(PVAL_WD, 3), "ns"), size = 2)
  
  
  funnel_Hmax_year  <- my_funnelplot("c) Hmax",
                                    ggplot(GIi[["Hmax"]],aes(x=year,y=corr.r,colour=factor(stage), size=2, alpha=0.6))) +
                                    geom_point() + 
                                    scale_y_continuous("Correlation coefficient r",limits=c(-1,1))+
                                    scale_x_continuous("",limits=c(1990,2015)) +
                                    theme( legend.position="none")+
                                    annotate("text", x = 1995, y = -0.9, label = paste("LRT:", round(LRT_Hmax,0)), size = 2) +
                                    annotate("text", x = 1995, y = -1, label = paste("p.value =", round(PVAL_Hmax, 3), "ns"), size = 2)
                                  
  
  funnel_Seedmass_year  <- my_funnelplot("d) Seed mass",
                                        ggplot(GIi[["Seedmass"]],aes(x=year,y=corr.r,colour=factor(stage), size=2, alpha=0.6))) +
                                        geom_point() + 
                                        scale_x_continuous("Years of publication",limits=c(1990,2015))+
                                        scale_y_continuous("",limits=c(-1,1)) +
                                        theme( legend.position="none")+
                                        annotate("text", x = 1995, y = -0.9, label = paste("LRT:", round(LRT_Seedmass,0)), size = 2) +
                                        annotate("text", x = 1995, y = -1, label = paste("p.value =", round(PVAL_Seedmass, 3), "**"), size = 2)
                                      
  
  funnel_Aarea_year <- my_funnelplot("e) Aarea",
                                     ggplot(GIi[["Aarea"]],aes(x=year,y=corr.r,colour=factor(stage),size=2, alpha=0.6))) +
                                      geom_point() + 
                                      scale_alpha(guide = 'none')+
                                      scale_size(guide = 'none') + 
                                      scale_y_continuous("Correlation coefficient  r",limits=c(-1,1)) +
                                      scale_x_continuous("Years of publication",limits=c(1990,2015))+
                                      theme (legend.title=element_blank(), legend.justification=c(0,0),legend.position=c(1.2,0.5), legend.key = element_blank(), plot.margin=unit(c(0,0,0,0),"mm"))+
                                      annotate("text", x = 1995, y = -0.9, label = paste("LRT:", round(LRT_Aarea,0)), size = 2) +
                                      annotate("text", x = 1995, y = -1, label = paste("p.value =", round(PVAL_Aarea, 3), "**"), size = 2)
                                    
  
  grid.arrange(funnel_SLA_year,funnel_WD_year,funnel_Hmax_year, funnel_Seedmass_year ,funnel_Aarea_year, nrow=3, ncol=2)
}




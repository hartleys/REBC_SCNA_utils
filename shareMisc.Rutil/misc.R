
#THIS SCRIPT ASSUMES THE EXISTANCE OF A VARIABLE CALLED "SHAREMISC.DIR"!!!
#  SHAREMISC.DIR <- ifelse(dir.exists("T:/"),"T:/shared/hartleys/software/shareMisc/","/mnt/nfs/gigantor/ifs/Shared/hartleys/software/shareMisc/");
#  source( paste0( SHAREMISC.DIR , "shareMisc.Rutil/misc.R" ) )
SOFTWARE.DIR <- paste0(  SHAREMISC.DIR , "shareMisc.Rutil/" )
JS.dir <- paste0(SOFTWARE.DIR,"JunctionSeq/")
QT.dir <- paste0(SOFTWARE.DIR,"QoRTs/")

source(paste0(SOFTWARE.DIR,"misc.p2.R"));
source(paste0(JS.dir,"00.minor.utils.R"));
source(paste0(JS.dir,"plotting.helpers.R"));
source(paste0(QT.dir,"minor.utils.R"));
source(paste0(QT.dir,"internal.plotting.func.R"));
source(paste0(SOFTWARE.DIR,"qqPlot.R"));

z.na <- function(x,z){
  ifelse(is.na(x),z,x)
}
z.inf <- function(x,z){
  ifelse(is.finite(x),x,z)
}
z.low <- function(x,k,z){
  ifelse(x > k,x,z)
}

z.nan.inf <- function(x,z){
  ifelse(is.nan(x) | is.infinite(x),z,x);
}

nonNumsToDots <- function(ddd,NA.MASK=".",NAN.MASK=NULL,verbose=FALSE){
  looks.like.number <- sapply(seq_along(names(ddd)),function(i){
    any(! is.na( as.numeric( ddd[[i]] ) ) )
  })
  ddd.fmt <- ddd;
  for( i in which(looks.like.number) ){
    if(verbose){
      examples.1 <- ddd[[i]][1:min(c(nrow(ddd),5))];
      examples.2 <- ddd[[i]][! is.na( as.numeric( ddd[[i]] ) )];
      examples.2 <- examples.2[1:min(c(length(examples.2),5))]
      examples.3 <- ddd[[i]][ is.na( as.numeric( ddd[[i]] ) )];
      examples.3 <- examples.3[1:min(c(length(examples.3),5))]
      message("Masking To Numeric: ",names(ddd)[[i]],"(examples: ",p(examples.1,collapse=","),")(examples: ",p(examples.2,collapse=","),")(examples: ",p(examples.3,collapse=","),")");
    }
    ddd.fmt[[i]] <- ddd[[ii]]
    if(! is.null(NAN.MASK)){
      nan.ct <- sum(ddd[[i]] == "NaN" | is.nan(ddd[[i]]),na.rm=T)
      message("     NAN.CT = ",nan.ct);
      ddd.fmt[[i]] <- ifelse(ddd[[i]] == "NaN" | is.nan(ddd[[i]]) ,NAN.MASK,  ddd[[i]])
      ddd.fmt[[i]] <- as.numeric(ddd.fmt[[i]])
    } else {
      ddd.fmt[[i]] <- as.numeric(ddd[[i]])
    }
    ddd.fmt[[i]] <- z.na( ddd.fmt[[i]], NA.MASK);
  }
  return(ddd.fmt);
}

getFirstDigit <- function(x){
  ifelse(x == 0,0, x / 10^floor(log10(x)))
}

isFALSE <- function(x){
  identical(x,FALSE);
}
#ease-of-use aliases:
p <- paste0

# R has row.names and rownames, and only colnames but no col.names.
# Why? Why would you do this?
rowNames <- row.names
colNames <- colnames
col.names <- colnames


simplePlotWrapper <- function(plotfunc,file,PNG, height=7,width=7,units="in",pointsize=9,res=600){
   if(PNG){
     png(file,height=height,width=width,units=units,pointsize=pointsize,res=res);
   } else {
     x11(height*2,width*2)
   }
   plotfunc();
   if(PNG){
     dev.off();
   }
}


makeListWithDefaults <- function(x,dx){
  out <- dx;
  if(isTRUE(x)){
    out;
  } else if(length(x) == 0){
    out;
  } else if( ! is.list(x)){
    stop("error in makeListWithDefaults, input must be a list!");
  } else {
    for(i in seq_along(x)){
      out[[ names(x)[[i]] ]] <- x[[i]];
    }
  }
  message("LIST:")
  print(out)
  return(out);
}

pretty.mini <- function(axis.at){
      axis.diff <- diff(axis.at)[[1]]
      axis.diff.factor <- getFirstDigit(axis.diff);
      if(axis.diff.factor == 1){
        minitick.factor <- 5;
      } else if(axis.diff.factor == 5){
        minitick.factor <- 5;
      } else if(axis.diff.factor == 2){
        minitick.factor <- 2;
      } else if(axis.diff.factor == 4){
        minitick.factor <- 4;
      } else {
        minitick.factor <- 4;
      }
      miniticks <- unlist(lapply(c(min(axis.at)-axis.diff ,axis.at),function(z){
        z + axis.diff * (1:(minitick.factor-1) / (minitick.factor))
      }))
      return(miniticks);
}

setupPlot <- function(xlim,ylim,parList = list(),
   x.guidelines = FALSE,
   y.guidelines = "white",
   y.axis.at = NULL, x.axis.at=NULL, x.axis = T, y.axis =T,
   x.axis.miniticks = FALSE, y.axis.miniticks = FALSE,
   bg = "lightgray",
   borders = c(T,T,F,F),
   adjust.to.fit.yaxis.label = TRUE,
   log.y = FALSE){

  plot.new();
  do.call(par,parList);

  if(isTRUE(log.y)){
    if(is.null(y.axis.at)){
      y.axis.at <- floor(min(log10(ylim))):ceiling(max(log10(ylim)));
    }
  } else {
    if(is.null(y.axis.at)){
      y.axis.at <- pretty(ylim,10)
    }
  }
  if(is.null(x.axis.at)){
      x.axis.at <- pretty(xlim,10)
  }

  if(adjust.to.fit.yaxis.label){
    max.yax.width <- max( strwidth(y.axis.at,cex = par("cex.axis")) );
    yax.margin.size <- plotting.limits()[[1]] - device.limits()[[1]]
    yax.margin.leftover <- yax.margin.size - max.yax.width - strwidth("X") * 2 - strwidth("X",cex=par("cex.lab"))*2
    if(yax.margin.leftover < 0){
      message("ADJUSTING Y axis margin!")
      mar <- par("mar");
      addme <- - yax.margin.leftover  ##/ (strwidth("X",cex=par("cex.lab")) )
      mar[[2]] <- mar[[2]] + addme
      par(mar=mar);
      plot.window(ylim=ylim,xlim=xlim)
    }
  }
  ylim.actual <- ylim
  if(isTRUE(log.y)){
    ylim.actual <- log10(ylim);
  }

  plot.window(ylim=ylim.actual,xlim=xlim)
  pu <- par("usr")
  rect(pu[[1]],pu[[3]],pu[[2]],pu[[4]],col=bg,border=FALSE)


  if(y.axis){
    if(isTRUE(log.y)){
            decades_log10 <- (floor((ylim.actual[1]))-1):(ceiling((ylim.actual[2]))+1);
            decades_unscaled <- 10 ^ decades_log10;
            decade.labels <- c();
            if(any(decades_log10 == 0)){
              decade.labels <- c(expression(1));
            }
            if(any(decades_log10 == 1)){
              decade.labels <- c(decade.labels,expression(10));
            }
            if(decades_log10[1] < 0){
              decade.labels <- c(sapply(decades_log10[1]:(-1), function(X){
                                          substitute(10 ^ x, list(x = X));
                                       }),
                                decade.labels);
            }
            if(max(decades_log10) > 1){
              decade.labels <- c(decade.labels,
                                sapply(2:max(decades_log10), function(X){
                                          substitute(10 ^ x, list(x = X));
                                       })
                                );
            }
            axis(2,at=decades_log10,labels=decade.labels,las=1);
    } else {
      axis(2,at=y.axis.at ,las=1);
    }
  }
  if(! isFALSE(y.guidelines)){
      if(!is.list( y.guidelines )){
         y.guidelines <- list(col=y.guidelines)
      }
      y.guidelines[["h"]] <- y.axis.at;
      do.call(abline,y.guidelines);
      #abline(v = x.axis.at,col=x.guidelines,lty=2);
  }
  if(x.axis){
    axis(1,at=x.axis.at ,las=1);
  }
  if(! isFALSE(x.guidelines)){
      if(!is.list( x.guidelines )){
         x.guidelines <- list(col=x.guidelines)
      }
      x.guidelines[["v"]] <- x.axis.at;
      print(x.guidelines)
      do.call(abline,x.guidelines);
      #abline(v = x.axis.at,col=x.guidelines,lty=2);
  }

  if(! isFALSE(x.axis.miniticks)){
    miniticks <- pretty.mini(x.axis.at);
    axis.call <- makeListWithDefaults(x.axis.miniticks,list(1,at=miniticks,labels=FALSE,tcl=-0.25));
    do.call(axis,axis.call)
  }
  if(! isFALSE(y.axis.miniticks)){
    if(isTRUE(log.y)){
      axis.at <- y.axis.at;
      minitick.offsets <- log10( 1:9 )
      miniticks <- unlist(lapply(c(min(axis.at)-1 ,axis.at),function(z){
        z + minitick.offsets
      }))
    } else {
      miniticks <- pretty.mini(y.axis.at);
    }
    axis.call <- makeListWithDefaults(y.axis.miniticks,list(2,at=miniticks,labels=FALSE,tcl=-0.25));
    do.call(axis,axis.call)
  }

  if(borders[[1]]){ lines(c(pu[[1]],pu[[2]]),rep(pu[[3]],2), xpd=NA) }
  if(borders[[2]]){ lines(rep(pu[[1]],2),c(pu[[3]],pu[[4]]), xpd=NA) }
  if(borders[[3]]){ lines(c(pu[[1]],pu[[2]]),rep(pu[[4]],2), xpd=NA) }
  if(borders[[4]]){ lines(rep(pu[[2]],2),c(pu[[3]],pu[[4]]), xpd=NA) }

}






#
plotText <- function(x = c("center","left","right","CENTER","LEFT","RIGHT"),
                     y=c("center","top","bottom","CENTER","TOP","BOTTOM"),
                     labels="",
                     plotlim=TRUE,outer=FALSE,
                     srt=0,xpd=NA, spacerChar = "X", spacerCex = 0.5,debug.alignment.line = FALSE,...){
  devlim <- device.limits();
  poslim <- plotting.limits();
  posX <- tolower(match.arg(x));
  posY <- tolower(match.arg(y));
  if(length(plotlim) == 1){
    plotlim <- rep(plotlim,2);
  } else if(length(plotlim) != 2){
    stop("Error: plotlimmust be of length 1 or 2!");
  }
  if(length(outer) == 1){
    outer<- rep(outer,2);
  } else if(length(outer) != 2){
    stop("Error: outermust be of length 1 or 2!");
  }
  if( length(spacerChar) == 1){
    spacerChar <- rep(spacerChar,2);
  } else if(length(spacerChar) != 2){
    stop("Error: spacerChar must be of length 1 or 2!");
  }
  if( length(spacerCex) == 1){
    spacerCex <- rep(spacerCex,2);
  } else if(length(spacerCex) != 2){
    stop("Error: spacerCex must be of length 1 or 2!");
  }
  ss <- c(
      strwidth(spacerChar[[1]],cex=spacerCex[[1]]),
      strheight(spacerChar[[2]],cex=spacerCex[[2]])
  )

  if(( ! plotlim[[1]]) && outer[[1]] && posX != "center"){
    warning("WARNING: text is being placed outside the device limits!")
  }
  if(( ! plotlim[[2]]) && outer[[2]] && posY != "center"){
    warning("WARNING: text is being placed outside the device limits!")
  }
  xlims <- poslim;
  ylims <- poslim;
  if(! plotlim[[1]]){
    xlims <- devlim;
  }
  if(! plotlim[[2]]){
    ylims <- devlim;
  }
  x <- mean(xlims[1:2]);
  y <- mean(ylims[3:4])
  adj <- c(0.5,0.5);
  spc <- c(0,0)
  if(posX == "left"){
    x <- xlims[[1]];
    if(outer[[1]]){
      adj[[1]] <- 1; spc[[1]] <- -ss[[1]]
    } else {
      adj[[1]] <- 0; spc[[1]] <-  ss[[1]]
    }
  } else if(posX == "right"){
    x <- xlims[[2]];
    if(outer[[1]]){
      adj[[1]] <- 0; spc[[1]] <- ss[[1]]
    } else {
      adj[[1]] <- 1; spc[[1]] <- -ss[[1]]
    }
  }
  if(posY == "bottom"){
    y <- ylims[[3]];
    if(outer[[2]]){
      adj[[2]] <- 1;  spc[[2]] <- -ss[[2]]
    } else {
      adj[[2]] <- 0;  spc[[2]] <- ss[[2]]
    }
  } else if(posY == "top"){
    y <- ylims[[4]];
    if(outer[[2]]){
      adj[[2]] <- 0;  spc[[2]] <- ss[[2]]
    } else {
      adj[[2]] <- 1;  spc[[2]] <- -ss[[2]]
    }
  }
  if(debug.alignment.line){
    lines(devlim[1:2],c(y,y),lty=3,col="red",xpd=T,lwd=1);
    lines(c(x,x),devlim[3:4],lty=3,col="red",xpd=T,lwd=1);
    points(x + ss[[1]],y + ss[[2]],pch=4,xpd=NA)
  }

  text(x + spc[[1]],y + spc[[2]],labels,adj=adj,srt=srt,xpd=xpd,...);

}


darken <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}
lighten <- function(color, factor=1.4){
    col <- col2rgb(color) / 255
    colInv <- 1 - col;
    colInv <- colInv / factor;
    col <- (1 - colInv) * 255
    col <- rgb(t(col), maxColorValue=255)
    col
}


rainbowSpectrum <- function(N, alpha = 255, darkenYellows = TRUE, darkenRanges = list(c(0.7,0.8),c(0.86,0.8)), maxDarkenFactor = 1.2){
  cols <- rev(rainbow(n = ceiling(N * 1.1))[1:N])
  pct <- 1:N / N;
  spanWds <- sapply(darkenRanges,function(ab){ ab[2] - ab[1] })
  ctrDarken <- maxDarkenFactor - 1
  if(darkenYellows ){
   for(ppidx in 1:N){
    pp <- pct[ppidx];
    for(i in 1:length(spanWds)){
      dr <- darkenRanges[[i]];
      if( pp >= min(dr) & pp < max(dr)){
         ppr <- ((pp - dr[[1]])/spanWds[[i]])
         darkFactor <- (ppr * ctrDarken) + 1;
         darkCol <- darken(cols[[ppidx]],darkFactor);
         message("ppidx: ",ppidx," - ",darkFactor,":",cols[ppidx], " => ",darkCol );
         cols[[ppidx]] <- darkCol;
      }# else {
      #   message("pp(",pp,") not in [",min(dr),",",max(dr),"]");
      #}
    }
   }
  }
  cols <- sapply(cols,color2transparent,alpha);
  return(cols);
  #0.775-0.84
}



redgreenSpectrum <- function(N,alpha=1){
  oddn <- if(N %% 2 == 1){
     1
  } else { 0 }
  cols <- c( terrain.colors(N,alpha=alpha)[1:(N/2)],rev(heat.colors(N*0.66,alpha=alpha)[1:((N/2)+oddn)]))
  blend0ix <- N/2;
  blend1ix <- ceiling(N/2 + N * 0.1)
  if(blend0ix + 1 < blend1ix){
  blendZone <- blend0ix:blend1ix;
  blend0 <- col2rgb(cols[blend0ix])[,1];
  blend1 <- col2rgb(cols[blend1ix])[,1];
  blendFactor <- (length(blendZone)-1):0 / (length(blendZone)-1)
  #blendRGB <- col2rgb(cols[blendZone])
  cols[blendZone] <- sapply(1:length(blendZone),function(i){
   #r <- col2rgb(c)[,1]
   r <- blend0 * blendFactor[[i]] + blend1 * (1 - blendFactor[[i]]);
   rgb(red = r[[1]]/255, green = r[[2]]/255, blue = r[[3]]/255, alpha = alpha);
  })
  }

  return(cols[1:N]);
}


rainbowSpectrum <- function(N, alpha = 255, darkenYellows = TRUE, darkenRanges = list(c(0.7,0.8),c(0.86,0.8)), maxDarkenFactor = 1.2){
  cols <- rev(rainbow(n = ceiling(N * 1.1))[1:N])
  pct <- 1:N / N;
  spanWds <- sapply(darkenRanges,function(ab){ ab[2] - ab[1] })
  ctrDarken <- maxDarkenFactor - 1
  if(darkenYellows ){
   for(ppidx in 1:N){
    pp <- pct[ppidx];
    for(i in 1:length(spanWds)){
      dr <- darkenRanges[[i]];
      if( pp >= min(dr) & pp < max(dr)){
         ppr <- ((pp - dr[[1]])/spanWds[[i]])
         darkFactor <- (ppr * ctrDarken) + 1;
         darkCol <- darken(cols[[ppidx]],darkFactor);
         message("ppidx: ",ppidx," - ",darkFactor,":",cols[ppidx], " => ",darkCol );
         cols[[ppidx]] <- darkCol;
      }# else {
      #   message("pp(",pp,") not in [",min(dr),",",max(dr),"]");
      #}
    }
   }
  }
  cols <- sapply(cols,color2transparent,alpha);
  return(cols);
  #0.775-0.84
}


cheat.sheet <- function(lim = 122,cex=1,cex.text=cex){
  if(lim <= 25){
    pch.set <- 1:lim;
  } else {
    pch.set <- c(1:25,32:lim);
  }
  lim <- length(pch.set);

  plot.new();
  plot.window(xlim=c(0,1),ylim=c(0,1));
  charht <- strheight("A",cex=cex.text) * 1.3;
  charwd <- strwidth("A",cex=cex.text)
  ht <- abs(par("usr")[4] - par("usr")[3])
  col.len <- ceiling(ht / charht)
  nc <- ceiling(lim/ col.len);

  chunks <- split(pch.set, ceiling(nc * seq_along(1:lim)/ (lim) ))
  col.len <- length(chunks[[1]]);
  col.Y <- rev(seq(from=par("usr")[3] + charht,to=par("usr")[4] - charht,length.out=col.len));
  col.X <- seq(from=par("usr")[1]+charwd,to=par("usr")[2]-charwd,length.out=nc + 1);

  for(j in seq_len(nc)){
    chunk <- chunks[[j]];
    points(rep(col.X[j],length(chunk)),col.Y[1:length(chunk)],pch=chunk);
    text(rep(col.X[j],length(chunk)) + charwd * 1.5,col.Y[1:length(chunk)],paste0(" ",chunk));
  }
  points(mean(par("usr")[1:2]),par("usr")[4]-charht,pch=3,col="gray",xpd=NA);
  text(mean(par("usr")[1:2]),par("usr")[4]-charht,"(0,0)",adj=c(0,0),xpd=NA);
  text(mean(par("usr")[1:2]),par("usr")[4]-charht,"(0,1)",adj=c(0,1),xpd=NA);
  text(mean(par("usr")[1:2]),par("usr")[4]-charht,"(1,0)",adj=c(1,0),xpd=NA);
  text(mean(par("usr")[1:2]),par("usr")[4]-charht,"(1,1)",adj=c(1,1),xpd=NA);
  text(mean(par("usr")[1:2]) - charwd * 4,par("usr")[3]-charht,"adj: ",xpd=NA);
}


spaceItemsToFit <- function(x,xscale,minWidth, verbose = TRUE){
    xscale <- sort(xscale)
    coor <- sort(x)
    xspan <- diff(xscale)
    x <- (coor - xscale[1]) / xspan
    xt <- function(xx){ xx * xspan + xscale[1] }
    w <- minWidth/ xspan
    xd <- diff(c(0,x,1))
    xd.raw <- xd;
    smalldiff <- xd < w
    rr <- function(a){round(a / w,2)}
    rrd <- function(){
      data.frame(raw = rr(xd.raw),xd = rr(xd), sum=rr(cumsum(xd)))
    }
    if(verbose){message("minWidth=",minWidth)}
    if(verbose){print(data.frame(xd=xd,smalldiff=smalldiff))}
    if(verbose){print(rrd())}

    if(length(coor) * w > xspan){
        x.jit <- 1:length(coor) / length(coor) - (0.5 / length(coor))
        if(verbose){ message("steveJit label jitter fcn: Jittering labels overpacked. Distributing evenly instead.")}
        return(xt(x.jit))
    } else if( all(xd > w)){
        if(verbose){message("steveJit label jitter fcn: Jittering unnecessary. Returning raw coords")}
        return(coor)
    } else {
        if(verbose){message("steveJit label jitter fcn: Jittering labels to fit")}

        g.rle <- rle(smalldiff)
        g.rle <- data.frame(lengths = g.rle$lengths, values = g.rle$values)
        g.rle$i0 <- -1;
        g.rle$i1 <- cumsum(g.rle$lengths)
        g.rle$i0 <- c(1,g.rle$i1[-nrow(g.rle)] +1)
        s.rle <- g.rle[g.rle$values,]
        if(verbose){message("Element Packing RLE:")}
        if(verbose){print(g.rle)}
        if(verbose){message("Found ",nrow(s.rle)," groups of overpacked elements.")}
        if(verbose){print(s.rle)}
        if(verbose){print(rrd())}

        for( i in 1:nrow(s.rle)){
            out.wd <- s.rle$lengths[i] * w;
            i0 <- s.rle$i0[i]
            i1 <- s.rle$i1[i]
            raw.wd <- sum(xd[(i0):(i1)])
            xd[i0:i1] <- w
            amt <- out.wd - raw.wd;
            dist <- 1;
            if(verbose){message("Starting distribution: raw.wd=",rr(raw.wd),", out.wd=",rr(out.wd)," amt=",rr(amt),", dist=",dist,", i0=",i0,", i1=",i1)}
            while(amt > 0){
                ltAvail <- 0;
                rtAvail <- 0;
                lt <- i0 - dist;
                rt <- i1 + dist
                if(verbose){message("  iter:  xdlt=",rr(xd[lt]),", xdrt=",rr(xd[rt]))}
                if(lt <= 0 && rt > length(xd)){
                  warning("JITTER FAILED!")
                  amt <- 0;
                }
                if( lt > 0 && xd[lt] > w){
                  ltAvail <- xd[lt] - w
                }
                if( rt <= length(xd) && xd[rt] > w){
                  rtAvail <- xd[rt] - w;
                }
                ltRemove <- 0;
                rtRemove <- 0;
                if(ltAvail > amt / 2 && rtAvail > amt / 2){
                  ltRemove <- amt / 2
                  rtRemove <- amt / 2
                } else if(ltAvail < rtAvail){
                  ltRemove <- ltAvail
                  rtRemove <- min(rtAvail,amt - ltRemove)
                } else {
                  rtRemove <- rtAvail
                  ltRemove <- min(ltAvail,amt - rtRemove)
                }
                if(verbose){message("         dist = ",dist,", lt=",lt,", rt=",rt,", ltAvail=",rr(ltAvail),", rtAvail=",rr(rtAvail),", amt=",rr(amt))}
                amt <- amt - ltRemove - rtRemove
                xd[lt] <- xd[lt] - ltRemove
                xd[rt] <- xd[rt] - rtRemove
                dist <- dist + 1
                if(verbose){message("         dist = ",dist,", lt=",lt,", rt=",rt,", ltAvail=",rr(ltAvail),", rtAvail=",rr(rtAvail),", amt=",rr(amt))}
                if(verbose){message("         ltRemove=",rr(ltRemove),", rtRemove=",rr(rtRemove))}
                if(verbose){message("         xdlt=",rr(xd[lt]),", xdrt=",rr(xd[rt]))}
                if(verbose){print(rrd())}
                if(verbose){message("------------------------")}

            }
        }
        x.jit <- cumsum(xd[-length(xd)])
        return(xt(x.jit))

    }

}


pack.spans.into.rows <- function(spans = NULL,start = spans[[1]],end = spans[[2]], xlim = c(min(start),max(end)),buffer = 0){
  if(is.null(spans)){
    spans <- data.frame(start = start, end = end)
  }
  if(length(buffer) == 1){
    buffer <- rep(buffer,2)
  }
  #span.is.placed <- rep(FALSE,nrow(spans))
  #message("TEST: ",xlim[[1]],"-",xlim[[2]])
  r <- list( data.frame(index = numeric(),start = numeric(), end = numeric()) )
  for(i in 1:nrow(spans)){
    s <- spans[[1]][[i]]
    e <- spans[[2]][[i]]
    currRow <- 1;
    while(nrow(r[[currRow]]) != 0 && any(r[[currRow]]$start < e + buffer[[1]] & r[[currRow]]$end + buffer[[2]] > s) ){
      if(currRow == length(r)){
        r[[currRow + 1]] <- data.frame(index = numeric(),start = numeric(), end = numeric())
      }
      currRow <- currRow + 1;
    }
    xx <- r[[currRow]]
    r[[currRow]] <- rbind.data.frame(xx, data.frame(index = i, start = s, end = e))
  }
  print(r)
  out <- rep(-1,nrow(spans));
  for(i in 1:length(r)){
    out[ r[[i]]$index ] <- i;
  }
  if(any(out == -1)){
    stop("ERROR: pack.spans.into.rows() impossible state: element not assigned a row!")
  }
  out;
}

plot.miniticks <- function(ax = 2,log.y = TRUE,log.x=TRUE, lwd = -1, lwd.ticks = 0.5, tcl = -0.25,lim = NULL, invert=FALSE,...){
  if(invert){
  if(log.y & ax == 2){
    if(is.null(lim)){
      lim <- c(par("usr")[3],par("usr")[4]);
    }
    miniticks <- -log10(getAllOrdersOfTen(10 ^ (-lim[2]),10^(-lim[1]), skipBase=FALSE));
    axis(ax,at=miniticks,labels=FALSE,lwd = lwd, lwd.ticks = lwd.ticks,tcl=tcl, ...);
  } else if(log.x & ax == 1){
    if(is.null(lim)){
      lim <- c(par("usr")[1],par("usr")[2]);
    }
    miniticks <- -log10(getAllOrdersOfTen(10 ^ (-lim[2]),10^(-lim[1]), skipBase=FALSE));
    axis(ax,at=miniticks,labels=FALSE,lwd = lwd, lwd.ticks = lwd.ticks,tcl=tcl, ...);
  } else {
    #do nothing!
  }
  } else {


  if(log.y & ax == 2){
    if(is.null(lim)){
      lim <- c(par("usr")[3],par("usr")[4]);
    }
    miniticks <- log10(getAllOrdersOfTen(10 ^ lim[1],10^lim[2], skipBase=FALSE));
    axis(ax,at=miniticks,labels=FALSE,lwd = lwd, lwd.ticks = lwd.ticks,tcl=tcl, ...);
  } else if(log.x & ax == 1){
    if(is.null(lim)){
      lim <- c(par("usr")[1],par("usr")[2]);
    }
    miniticks <- log10(getAllOrdersOfTen(10 ^ lim[1],10^lim[2], skipBase=FALSE));
    axis(ax,at=miniticks,labels=FALSE,lwd = lwd, lwd.ticks = lwd.ticks,tcl=tcl, ...);
  } else {
    #do nothing!
  }
  }
}

f.na <- function(x){
  ifelse(is.na(x),FALSE,x);
}

JS.axis <- function(side, at, labels = at, tcl = -0.5, lwd = 1, lwd.ticks = 1, cex.axis = 1, srt = 0, line = NA, mgp.axis = 1, pos = NA, adj = c(0.5,1.1), font = 1, ...){

  if(side == 1){
    axis.height <- par("usr")[3]
    axis.tick.floor <- axis.height + (tcl * (2 * par("cxy")[2]))
    segments(at, axis.tick.floor,at, axis.height, lwd = lwd.ticks, cex = cex.axis, xpd = NA, ...)
    devlim <- device.limits()

    #if(is.na(line)){
    #  axis.label.height <- abs((axis.tick.floor - devlim[3])/2) + devlim[3];
    #} else {
      axis.label.height <- axis.height - par("cxy")[2] * mgp.axis
    #}

    text(x = at, y = axis.label.height, labels = labels, cex = cex.axis, adj = adj, srt = srt, xpd = NA, font = font,...)
    axis.floor <- axis.label.height - max(strheight(labels, cex = cex.axis))

    return(axis.floor)
  } else {
    mgp <- par("mgp");
    mgp[2] <- mgp.axis;
    axis(side,at=at,labels=labels,tcl=tcl,lwd=lwd,lwd.ticks=lwd.ticks,cex.axis=cex.axis,srt=srt,line=line,pos=pos,font=font,mgp =mgp,...);
  }
}


draw.axis.adv <- function(pos,at,labels = c(), mini.at = c(), side = 2,lim = c(min(c(at,mini.at),na.rm=T),max(c(at,mini.at),na.rm=T)), tcl = 0.25, mini.tcl = tcl / 2, axis.lty = 1, tick.lty = 1, axis.col = "black", tick.col = "black", label.spacing = 0.25,charwd = strwidth("-"),charht = strheight("|"), label.cex = 1,label.adj = NULL,...){
  if(side == 2){
    #charwd <- charwd;
    if(is.null(label.adj)){
      label.adj <- c(1.25,0.5)
    }
    tick.len <- tcl * charwd;
    mini.len <- mini.tcl * charwd;
    label.spacing.wd <- label.spacing * charwd;
    segments(pos,lim[1],pos,lim[2],lty=axis.lty,col = axis.col,...);
    segments(pos,at,pos - tick.len,at,lty=tick.lty,col = tick.col,...);
    if(! is.null(mini.at)){
      segments(pos,mini.at,pos - mini.len,mini.at,lty=tick.lty,col = tick.col,...);
    }
    if(! is.null(labels)){
      text(rep(pos - tick.len - label.spacing.wd,length(at)),at,labels,cex=label.cex,adj = label.adj,...)
    }
  } else if(side == 1){
    if(is.null(label.adj)){
      label.adj <- c(0.5,1.25)
    }
    #charht <- charht;
    tick.len <- tcl * charht * 0.75;
    mini.len <- mini.tcl * charht * 0.75;
    label.spacing.ht <- label.spacing * charwd;
    segments(lim[1],pos,lim[2],pos,lty=axis.lty,col = axis.col,...);
    segments(at,pos,at,pos - tick.len,lty=tick.lty,col = tick.col,...);
    if(! is.null(mini.at)){
      segments(mini.at,pos,mini.at,pos - mini.len,lty=tick.lty,col = tick.col,...);
    }
    if(! is.null(labels)){
      text(at,rep(pos - tick.len - label.spacing.ht,length(at)),labels,cex=label.cex,adj = label.adj,...)
    }
  } else if( side == 3){
    if(is.null(label.adj)){
      label.adj <- c(0.5,-0.25)
    }
    #charht <- charht;
    tick.len <- - tcl * charht * 0.75;
    mini.len <- - mini.tcl * charht * 0.75;
    label.spacing.ht <- - label.spacing * charwd;
    segments(lim[1],pos,lim[2],pos,lty=axis.lty,col = axis.col,...);
    message("lim = ",lim[1],"-",lim[2])
    segments(at,pos,at,pos - tick.len,lty=tick.lty,col = tick.col,...);
    if(! is.null(mini.at)){
      segments(mini.at,pos,mini.at,pos - mini.len,lty=tick.lty,col = tick.col,...);
    }
    if(! is.null(labels)){
      text(at,rep(pos - tick.len - label.spacing.ht,length(at)),labels,cex=label.cex,adj = label.adj,...)
    }
  }
}


calcTwoWayHistogram <- function(X,Y,X.breaks,Y.breaks, includeNA = TRUE){
  cat("[");
  counts <- sapply(1:length(Y.breaks),function(i){
    if(i == length(Y.breaks)){
      yvector <- is.na(Y);
    } else {
      lim.Y <- c(Y.breaks[i],Y.breaks[i+1])
      yvector <- f.na(Y >= lim.Y[1] & Y < lim.Y[2]);
    }
    cat(".");
    if(i %% 5 == 0){
      message("] Finished ",i," out of ",length(Y.breaks)-1," rows.");
      cat("[");
    }
    sapply(1:length(X.breaks),function(j){
      if(j == length(X.breaks)){
        xvector <- is.na(X);
      } else {
        lim.X <- c(X.breaks[j],X.breaks[j+1])
        xvector <- f.na(X >= lim.X[1] & X < lim.X[2]);
      }
      sum( xvector & yvector );
    });
  })

  X.titles <- c(paste0( X.breaks[-length(X.breaks)],"-",X.breaks[-1]),"NA");
  Y.titles <- c(paste0( Y.breaks[-length(Y.breaks)],"-",Y.breaks[-1]),"NA");

  if(includeNA){
    list(counts = counts, X = X, Y = Y,
     X.breaks = X.breaks,Y.breaks = Y.breaks,
     X.titles = X.titles, Y.titles = Y.titles,
     includesNA = TRUE);
   } else {
    list(counts = counts[-nrow(counts),-ncol(counts)], X = X, Y = Y,
     X.breaks = X.breaks,Y.breaks = Y.breaks,
     X.titles = X.titles[-length(X.titles)], Y.titles = Y.titles[-length(Y.titles)],
     includesNA = FALSE);
   }
}

plot.matrix.heatmap <- function(
        counts,
        X.axis.titles=NULL,
        Y.axis.titles = NULL,
        colorCt = 100,
        totals = TRUE,
        plot.fracs=FALSE,
        count.cex = 1,
        display.counts = TRUE,
        cex.axis = 1,
        color.spectrum.func = heat.colors,...){
  N <- sum(counts);
  if(plot.fracs){
    counts <- apply(counts,MAR=c(1,2),function(a){ round(100 * a / sum(counts),4) })
  }

  plot.new();
  par(xaxs="i",yaxs="i",mar=c(5.5,5.5,5.5,5.5)+0.1);
  plot.window(xlim=c(-0.5,nrow(counts) - 0.5 + totals),ylim=c(-0.5,ncol(counts) - 0.5 + totals));
  if(! is.null(X.axis.titles) ){
    JS.axis(1,at=0:nrow(counts),labels=X.axis.titles,tcl=0,las=1, cex = cex.axis);
    axis(1,at=0:nrow(counts) - 0.5, labels=FALSE, cex = cex.axis);
  }
  if(! is.null(Y.axis.titles) ){
    JS.axis(2,at=0:ncol(counts),labels=Y.axis.titles,tcl=0,las=1, cex = cex.axis);
    axis(2,at=0:ncol(counts) - 0.5, labels=FALSE, cex = cex.axis);
  }

  maxVal <- max(counts)
  color.spectrum <- rev(heat.colors(colorCt + 1));
  for(i in 1:nrow(counts)){
    #normVal <- log(GQ.vs.AD.counts.2[i,,drop=TRUE]) / log(maxVal)
    #normVal <- ifelse(is.finite(normVal),normVal,0);
    normVal <- counts[i,,drop=TRUE] / maxVal
    color.idx <- floor( normVal * colorCt ) + 1;
    rect(xleft=i-1.5, ybottom=1:length(normVal) - 1.5, xright=i-0.5, ytop=1:length(normVal) - 0.5,col=color.spectrum[color.idx],border="transparent");
    if(display.counts){
      text(i-1,1:length(normVal) - 1,labels=counts[i,,drop=TRUE],cex=count.cex)
    }
  }
  if(totals){
    ytotals <- rowSums(counts)
    xtotals <- colSums(counts)
    for(i in 1:length(ytotals)){
      normVal <- ytotals / max(ytotals)
      color.idx <- floor( normVal * colorCt ) + 1;
      rect(xleft=1:length(ytotals) - 1.5, ybottom=ncol(counts)-0.5, xright=1:length(ytotals) - 0.5, ytop=ncol(counts)+0.5,col=color.spectrum[color.idx],border="transparent");
      if(display.counts){
        text(1:length(ytotals) - 1,ncol(counts),labels=ytotals ,cex=count.cex)
      }
    }
    for(i in 1:length(xtotals)){
      normVal <- xtotals / max(xtotals)
      color.idx <- floor( normVal * colorCt ) + 1;
      rect(xleft=nrow(counts)-0.5, ybottom=1:length(xtotals)-1.5, xright=nrow(counts)+0.5, ytop=1:length(xtotals)-0.5,col=color.spectrum[color.idx],border="transparent");
      if(display.counts){
        text(nrow(counts),1:length(xtotals)-1,labels=xtotals ,cex=count.cex)
      }
    }
    abline(v=nrow(counts)-0.5,col="black",lty=3);
    abline(h=ncol(counts)-0.5,col="black",lty=3);
  }
  text(par("usr")[2],par("usr")[4],paste0("N=",N),adj=c(1.1,-1.25),xpd=NA)
  box();
}


draw.logyaxis.stdScalePlot <- function(ylim, ylim.truncate, lwd, lwd.mini, tcl = -0.5, side = 2, label.style = c("base10","raw","none"),...){
  label.style <- match.arg(label.style);
  if(missing(ylim)){
    ylim <- c(par("usr")[3], par("usr")[4]);
  }
  if(missing(lwd)){
    lwd <- par("lwd");
  }
  if(missing(lwd.mini)){
    lwd.mini <- lwd / 2;
  }
  decades_log10 <- (floor((ylim[1]))-1):(ceiling((ylim[2]))+1);
  decades_unscaled <- 10 ^ decades_log10;

  ticks_unscaled <- c();
  for(i in 1:length(decades_unscaled)){
    ticks_unscaled <- c(ticks_unscaled,
                        decades_unscaled[i]*2,
                        decades_unscaled[i]*3,
                        decades_unscaled[i]*4,
                        decades_unscaled[i]*5,
                        decades_unscaled[i]*6,
                        decades_unscaled[i]*7,
                        decades_unscaled[i]*8,
                        decades_unscaled[i]*9
                        );
  }
  ticks_log10 <- log10(ticks_unscaled);

  if(label.style == "base10"){
      decade.labels <- c();
      if(any(decades_log10 == 0)){
        decade.labels <- c(expression(1));
      }
      if(any(decades_log10 == 1)){
        decade.labels <- c(decade.labels,expression(10));
      }
      if(decades_log10[1] < 0){
        decade.labels <- c(sapply(decades_log10[1]:(-1), function(X){
                                    substitute(10 ^ x, list(x = X));
                                 }),
                          decade.labels);
      }
      if(max(decades_log10) > 1){
        decade.labels <- c(decade.labels,
                          sapply(2:max(decades_log10), function(X){
                                    substitute(10 ^ x, list(x = X));
                                 })
                          );
      }
  } else if(label.style == "raw"){
    decade.labels <- decades_unscaled;
  } else if(label.style == "none"){
    decade.labels <- rep("",length(decades_log10));
  }

  if(! missing(ylim.truncate)){
    keep_decades <- decades_log10 >= ylim.truncate[1] & decades_log10 <= ylim.truncate[2];
    decades_unscaled <- decades_unscaled[keep_decades];
    decade.labels <- decade.labels[keep_decades];
    decades_log10 <- decades_log10[keep_decades];
    ticks_unscaled <- ticks_unscaled[ticks_log10 >= ylim.truncate[1] & ticks_log10 <= ylim.truncate[2]];
    ticks_log10 <-  ticks_log10[ticks_log10 >= ylim.truncate[1] & ticks_log10 <= ylim.truncate[2]];
  }

  axis(side, at = decades_log10, labels=as.expression(decade.labels), lwd = -1, lwd.ticks = lwd, tcl = tcl, las = 1, ...);
  axis(side, at = ticks_log10, labels=FALSE, lwd = -1, lwd.ticks = lwd.mini, tcl = tcl / 2, ...);
}


plotKD.violin <- function(v, type = c("violin","left","right","up","down"),
                             col="magenta",border="black",
                             lty=1,
                             lwd=1,
                             rectCol="black",
                             colMed="white",
                             pchMed=19,
                             at=NULL,
                             add=FALSE,
                             wex=1,
                             drawRect=TRUE,
                             horizontal=FALSE,
                             wd = 0.95,
                             rectWd = wd / 16,
                             truncateToPlot = TRUE,
                             whiskerQuantile = NULL,
                             densityOptions = list(),
                             ...){
  if(type == "up" || type == "down") horizontal <- TRUE;
  if(type == "left" || type == "right") horizontal <- FALSE;

  if(!is.null(whiskerQuantile)){
    if(length(whiskerQuantile) == 1){
      whiskerQuantile <- c(whiskerQuantile,1-whiskerQuantile);
      whiskerQuantile <- c(min(whiskerQuantile),max(whiskerQuantile));
    }
  }

  #wd <- 0.95;
  #rectWd <- wd / 16;
  h <- do.call(density,c(list(x=v),densityOptions))
  #h <- density(v);
  qt <- quantile(v, probs = c(0.25,0.5,0.75));
  if(! is.null(whiskerQuantile)){
    whiskers <- quantile(v, probs = whiskerQuantile);
  }

  type <- match.arg(type);
  if(is.null(at)){
    at <- 1
  }

  if(! add){
    plot.new();
    if(horizontal){
      plot.window(ylim=c(at-wd,at+wd),xlim=c(min(h$x),max(h$x)));
      axis(1)
    } else {
      plot.window(xlim=c(at-wd,at+wd),ylim=c(min(h$x),max(h$x)));
      axis(2)
    }
    box();
  }


  if(truncateToPlot){
    if(horizontal){
      hkeep <- h$x > par("usr")[1] & h$x < par("usr")[2];
      h$x <- c(par("usr")[1],h$x[hkeep],par("usr")[2]);
      h$y <- c(0,h$y[hkeep],0);

    } else {
      hkeep <- h$x > par("usr")[3] & h$x < par("usr")[4];
      h$x <- c(par("usr")[3],h$x[hkeep],par("usr")[4]);
      h$y <- c(0,h$y[hkeep],0);
    }
  }

  if(! horizontal){

    if(type=="violin"){
      y.left <- h$x;
      y.rght <- h$x;
      rawXSpan <- abs(max(h$y) - min(h$y));
      x.left <- at - ((h$y / rawXSpan) * (wd/2))
      x.rght <- ((h$y / rawXSpan) * (wd/2)) + at;

      polygon(c(x.left,rev(x.rght)),
              c(y.left,rev(y.rght)),border=border,col=col,lty=lty,lwd=lwd,...);

      if(drawRect) rect(at-rectWd/2,qt[1],at+rectWd/2,qt[3],border=border,col=rectCol,...);
      lines(c(at-rectWd/2,at+rectWd/2),c(qt[2],qt[2]),col=colMed,...);
      if(!is.null(whiskerQuantile)) segments(x0=at-rectWd/4,x1=at+rectWd/4,y0=c(qt[1],qt[2]),y1=c(whiskers[1],whiskers[2]),col=border,...)
    } else if(type == "left"){
      y <- h$x;
      rawXSpan <- abs(max(h$y) - min(h$y));
      x <- (h$y / rawXSpan) * (wd) + at;

      polygon(c(x),
              c(y),border=border,col=col,lty=lty,lwd=lwd,...);

      if(drawRect) rect(at,qt[1],at+rectWd,qt[3],border=border,col=rectCol,...);
      #points(at+rectWd,qt[2],pch="<",col=colMed,...);
      lines(c(at,at+rectWd),c(qt[2],qt[2]),col=colMed,...);
      if(!is.null(whiskerQuantile)) segments(x0=at+rectWd/2,x1=at+rectWd/2,y0=c(qt[1],qt[2]),y1=c(whiskers[1],whiskers[2]),col=border,...)
    } else if(type == "right"){
      y <- h$x;
      rawXSpan <- abs(max(h$y) - min(h$y));
      x <- at - (h$y / rawXSpan) * (wd)

      polygon(c(x),
              c(y),border=border,col=col,lty=lty,lwd=lwd,...);

      if(drawRect) rect(at,qt[1],at-rectWd,qt[3],border=border,col=rectCol,...);
      #points(at-rectWd,qt[2],pch=">",col=colMed,...);
      lines(c(at,at-rectWd),c(qt[2],qt[2]),col=colMed,...);
      if(!is.null(whiskerQuantile)) segments(x0=at-rectWd/2,x1=at-rectWd/2,y0=c(qt[1],qt[2]),y1=c(whiskers[1],whiskers[2]),col=border,...)

    }
  } else {
    if(type=="violin"){


      y.left <- h$x;
      y.rght <- h$x;
      rawXSpan <- abs(max(h$y) - min(h$y));
      x.left <- at - ((h$y / rawXSpan) * (wd/2))
      x.rght <- ((h$y / rawXSpan) * (wd/2)) + at;

      polygon(c(y.left,rev(y.rght)),
              c(x.left,rev(x.rght)),border=border,col=col,lty=lty,lwd=lwd,...);

      if(drawRect) rect(qt[1],at-rectWd/2,qt[3],at+rectWd/2,border=border,col=rectCol,...);
      lines(c(qt[2],qt[2]),c(at-rectWd/2,at+rectWd/2),col=colMed,...);
      if(!is.null(whiskerQuantile)) segments(y0=at-rectWd/4,y1=at+rectWd/4,x0=c(qt[1],qt[2]),x1=c(whiskers[1],whiskers[2]),col=border,...)

    } else if(type == "up"){
      y <- h$x;
      rawXSpan <- abs(max(h$y) - min(h$y));
      x <- (h$y / rawXSpan) * (wd) + at;

      polygon(c(y),c(x),border=border,col=col,lty=lty,lwd=lwd,...);

      if(drawRect) rect(qt[1],at,qt[3],at+rectWd,border=border,col=rectCol,...);
      #points(qt[2],at+rectWd,pch="v",col=colMed,...);
      lines(c(qt[2],qt[2]),c(at,at+rectWd),col=colMed,...);
      if(!is.null(whiskerQuantile)) segments(y0=at+rectWd/2,y1=at+rectWd/2,x0=c(qt[1],qt[2]),x1=c(whiskers[1],whiskers[2]),col=border,...)

    } else if(type == "down"){
      y <- h$x;
      rawXSpan <- abs(max(h$y) - min(h$y));
      x <- at - (h$y / rawXSpan) * (wd)

      polygon(c(y),c(x),border=border,col=col,lty=lty,lwd=lwd,...);

      if(drawRect) rect(qt[1],at,qt[3],at-rectWd,border=border,col=rectCol,...);
      #points(qt[2],at-rectWd,pch="^",col=colMed,...);
      lines(c(qt[2],qt[2]),c(at,at-rectWd),col=colMed,...);
      if(!is.null(whiskerQuantile)) segments(y0=at-rectWd/2,y1=at-rectWd/2,x0=c(qt[1],qt[2]),x1=c(whiskers[1],whiskers[2]),col=border,...)

    }
  }
}



#############################################
#############################################


functionFactory <- function(fun){
    function(...) {
        warn <- err <- NULL
        res <- withCallingHandlers(
            tryCatch(fun(...), error=function(e) {
                err <<- paste0("Error thrown: ",e$call,":",e$message)
                NULL
            }), warning=function(w) {
                warn <<- append(warn, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        list(res, warn=warn, err=err)
    }
}

test <- function(i){
    switch(i, "1"={stop("oops"); i}, "2"={ warning("hmm"); i }, i)
}


catchWarning <- function(expr,warn,fullMatch=FALSE){
  withCallingHandlers(
    tryCatch({
      eval(expr);
    }), warning=function(w) {
      if((fullMatch & w$message != warn) | (! startsWith(w$message,warn))){
        message(paste0("Warning thrown in: ",w$call,":\n   ",w$message));
      } else {
        #suppress matched warning...
      }
      invokeRestart("muffleWarning")
    }
  )
}

warnFunct <- function(){
  warning("TEST WARNING");
  return(10);
}

rho.test <- function(X,Y,alpha=0.05){
  test <- catchWarning({
    cor.test(Y,X,method="spearman")
  },warn="Cannot compute exact p-value with ties");
  #Z <- qnorm(c(alpha/2,1-(alpha/2)))
  #CI <- tanh(atanh(test$estimate) + Z/sqrt(sum(keep)-3))
  c(
       rho = as.numeric(test$estimate),
       rho.Lo  = NA,#as.numeric(CI[1]),
       rho.Hi  = NA,#as.numeric(CI[2]),
       p.value = as.numeric(test$p.value)
  );
}


rho.test.ci <- function(X,Y,alpha=0.05){
  test <- catchWarning({
    cor.test(Y,X,method="spearman")
  },warn="Cannot compute exact p-value with ties");
  #Z <- qnorm(c(alpha/2,1-(alpha/2)))
  #CI <- tanh(atanh(test$estimate) + Z/sqrt(sum(keep)-3);
  rho <- as.numeric(test$estimate)
  N <- length(X);
  Z <- qnorm(1-(alpha/2))
  ZR <- (1/2) * log((1 + rho) / (1-rho));
  ZDIFF <- Z * sqrt((1 + (rho^2/2)) / (N-3));
  Z.LO <- ZR - ZDIFF;
  Z.HI <- ZR + ZDIFF;
  R.LO <- (exp(2 * Z.LO) - 1) / (exp(2 * Z.LO) + 1)
  R.HI <- (exp(2 * Z.HI) - 1) / (exp(2 * Z.HI) + 1)

  pval <- as.numeric(test$p.value);

  c(
       rho = as.numeric(test$estimate),
       rho.Lo  = Z.LO,
       rho.Hi  = Z.HI,
       p.value = as.numeric(test$p.value)
  );
}

rho.CI <- function(rho,N,alpha){
  Z <- qnorm(1-(alpha/2))
  ZR <- (1/2) * log((1 + rho) / (1-rho));
  ZDIFF <- Z * sqrt((1 + (rho^2/2)) / (N-3));
  Z.LO <- ZR - ZDIFF;
  Z.HI <- ZR + ZDIFF;
  R.LO <- (exp(2 * Z.LO) - 1) / (exp(2 * Z.LO) + 1)
  R.HI <- (exp(2 * Z.HI) - 1) / (exp(2 * Z.HI) + 1)
  c(R.LO,R.HI);
}




#singlePoint.AOC <- function(TPR,FPR){
#  left.triangle.area <- TPR * FPR / 2
#  top.triangle.area <- (1-TPR)*(1-FPR)/2
#  square.area <- TPR * (1-FPR)
#  return( left.triangle.area + top.triangle.area + square.area );
#}

#Calculates the AUC of the ROC curve given a single TPR/FPR point.
#   This AUC will be equal for any point on the isocost line, which is
#   parallel to the unity line.
#For example, any point between TPR/FPR: 0.5/0 and 1/0.5 will have an AUC of 0.75

singlePoint.AOC <- function(TPR,FPR){
  return( (1 + TPR - FPR) / 2 );
}


#singlePoint.AOC(0.5,   0)
#singlePoint.AOC(1,   0.5)
#singlePoint.AOC(1,   0.5)
isocost.line.fcn <- function(fpr, tpr.intersect){
  tpr <- tpr.intersect + fpr;
  if(tpr > 1 | tpr < 0){
     return(NA);
  } else {
     return(tpr);
  }
}




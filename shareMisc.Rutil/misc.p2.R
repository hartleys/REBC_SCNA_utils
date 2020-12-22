f.na <- function(x){
  ifelse(is.na(x),FALSE,x);
}
library("Biostrings");






randRound = function(x){
  amt = x - floor(x);
  return( floor(x) + rbinom(x,size=1,prob = amt) );
}


library(RColorBrewer)

addSpacesToWidth <- function(A,B,wd,cex = 1){
    i <- 1;
    while( strwidth(paste0(A,paste0(rep(" ",i),collapse=""),B),cex=cex) < wd ){
      i <- i + 1;
    }
    return( paste0(A,paste0(rep(" ",i),collapse=""),B) )
}
addSpacesToWidth.multiLine <- function(A,B,wd,cex = 1){
    As <- strsplit(A,"\n")[[1]]
    Bs <- strsplit(B,"\n")[[1]]
    if(length(As) > length(Bs)){
      Bs <- c(Bs,rep("",length(As)-length(Bs)))
    } else if(length(As) < length(Bs)){
      As <- c(As,rep("",length(Bs)-length(As)))
    }
    ABs <- sapply(1:length(As),function(i){
      addSpacesToWidth(As[[i]],Bs[[i]],wd=wd,cex=cex)
    })
    return( paste0(ABs,collapse="\n") )
}

grepv <- function(pattern,x,...){
   out <- x;
   for(i in 1:length(pattern)){
     out <- grep(pattern[[i]],out,...,value=TRUE)
   }
   return(out);
}

getMode <- function(v){
   uniqv <- unique(v);
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

mfloor <- function(x,base){
   base * floor(x / base);
}
mround <- function(x,base){
   base * round(x / base);
}

z.nan <- function(x,z){
  ifelse(is.nan(x),z,x);
}














find.pretty.lims <- function(x,n=10){
  px <- pretty(range(x),n)
  hi <- px[length(px)]
  if(any(x > hi)){
    hi <- hi + (px[length(px)] - px[length(px) - 1])
  }
  lo <- px[1];
  if(any(x < lo)){
    lo <- lo - (px[2] - px[1])
  }
  c(lo,hi)
}



convert.txCoord.to.genoCoord <- function(txCoord,refData){
  exonNum <- sapply(1:length(txCoord),function(i){
    if(is.na(txCoord[i])){
      NA;
    } else {
      which(
        refData$tx.starts <= txCoord[i] &
        refData$tx.ends   >= txCoord[i]
      );
    }
  });
  genoPos <- sapply(1:length(txCoord),function(i){
    if(is.na(txCoord[i])){
      NA;
    } else {
      refData$ref.starts[exonNum[i]] + (txCoord[i] - refData$tx.starts[exonNum[i]]);
    }
  });
  genoPos;
}



convert.HGVS.variant.to.coord <- function(variantPosition, refData){
  buf <- substr(variantPosition,3,nchar(variantPosition));
  buf  <- ifelse(startsWith(variantPosition,"c."),buf,NA);
  buf2 <- strsplit(buf,"[-\\+]");
  offset.direction <- grepl("-",buf);
  offset.direction <- ifelse(offset.direction,-1,1);

  exonBase <- as.numeric(sapply(buf2, function(x){
    if(length(x) == 1){
      x[1]
    } else if (length(x) == 2){
      x[1]
    } else {
      stop("!!!!");
    }
  }))
  offset <- offset.direction * as.numeric(sapply(buf2, function(x){
    if(length(x) == 1){
      0
    } else if (length(x) == 2){
      x[2]
    } else {
      stop("!!!!");
    }
  }))

  exonNum <- sapply(1:length(variantPosition),function(i){
    if(is.na(exonBase[i])){
      NA;
    } else {
      which(
        refData$tx.starts <= exonBase[i] &
        refData$tx.ends   >= exonBase[i]
      );
    }
  });

  exonicStartPos <- sapply(1:length(variantPosition),function(i){
    if(is.na(exonBase[i])){
      NA;
    } else {
      refData$ref.starts[exonNum[i]] + (exonBase[i] - refData$tx.starts[exonNum[i]]);
    }
  });

  finalRefPos <- exonicStartPos + offset;

  out <- data.frame(
    txPos = variantPosition,
    exonBase = exonBase,
    offset = offset,
    exonNum = exonNum,
    exonicStartPos = exonicStartPos,
    refPos = finalRefPos
  )
  return(out);
}

convert.exon.del.to.coord <- function(varStart,varEnd, refData){
  buf  <- substr(varStart,3,nchar(varStart));
  buf  <- ifelse(startsWith(varStart,"e."),buf,NA);
  startExonNum <- as.numeric(buf);
  startPos <- sapply(startExonNum ,function(x){ if(is.na(x)){ NA } else { refData$ref.starts[max(x,1)] } })
  buf  <- substr(varEnd,3,nchar(varEnd));
  buf  <- ifelse(startsWith(varEnd,"e."),buf,NA);
  endExonNum <- as.numeric(buf);
  endPos <- sapply(endExonNum,function(x){ if(is.na(x)){ NA } else { refData$ref.ends[max(x,1)] }  })
  
  startExonBase <- sapply(startExonNum ,function(x){ if(is.na(x)){ NA } else { refData$tx.starts[max(x,1)] } })
  endExonBase   <- sapply(endExonNum,   function(x){ if(is.na(x)){ NA } else { refData$tx.ends[max(x,1)] }  })

  out <- data.frame(
    varStart = varStart ,
    varEnd = varEnd,
    startExonBase = startExonBase ,
    endExonBase = endExonBase   ,
    startExonNum = startExonNum ,
    endExonNum = endExonNum ,
    startPos = startPos,
    endPos = endPos
  )
  return(out);
}


color2transparentVector <- function(c,t){
   sapply(c, FUN = color2transparent, t = t)
}
color2transparent <- function(c,t){
   r <- col2rgb(c,alpha=TRUE)
   return(rgb(r[1],r[2],r[3],t,maxColorValue = 255))
}

generate.interval.scale <- function(gene.features, exon.fraction, rescaleFunction = c("sqrt","log","linear","34root"), debug.mode = FALSE){
  rescaleFunction <- match.arg(rescaleFunction)
  if(rescaleFunction == "34root"){
    resFunc <- function(x){ x ^ 0.75 }
  } else if(rescaleFunction == "log"){
    resFunc <- function(x){ ifelse(x == 0, 0, log(x)) }
  } else if(rescaleFunction == "sqrt"){
    resFunc <- function(x){ sqrt(x) }
  } else {
    resFunc <- function(x){x}
  }

  if(! "mod" %in% colnames(gene.features)){
    gene.features$mod <- rep(1,nrow(gene.features));
  }
  
  MAX <- max(c(gene.features$start, gene.features$end))
  MIN <- min(c(gene.features$start, gene.features$end))
  SPAN <- MAX - MIN
  
  if(is.na(exon.fraction) || exon.fraction <= 0 || exon.fraction >= 1){ #if rescaling is off, do not rescale!
    iv <- data.frame(start = MIN, end = MAX, span = SPAN,
                       rescale.start = 0, rescale.end = 1, normSpan = 1)
    return(iv)
  }
  
  exon.iv <- data.frame(start = gene.features$start[gene.features$is.exon], 
                          end = gene.features$end[gene.features$is.exon],
                          mod = gene.features$mod[gene.features$is.exon]);
  exon.iv$span <- exon.iv$end - exon.iv$start
  intr.iv <- data.frame(start = exon.iv$end[-nrow(exon.iv)], end = exon.iv$start[-1], mod = 1)
  if(max(exon.iv$end) < MAX){
    intr.iv <- rbind.data.frame(intr.iv,data.frame(start = max(exon.iv$end), end = MAX, mod = 1))
  }
  if(min(exon.iv$start) > MIN){
    intr.iv <- rbind.data.frame(data.frame(start = MIN, end = min(exon.iv$start), mod = 1),intr.iv)
  }
  exon.iv$span <- exon.iv$end - exon.iv$start
  intr.iv$span <- intr.iv$end - intr.iv$start
  
  
  if((nrow(intr.iv) > 0) && sum(intr.iv$span) > 0 ){
    exon.iv$type <- "E"
    intr.iv$type <- "I"
    exon.iv$normSpan <- exon.fraction * resFunc(exon.iv$span) * exon.iv$mod / sum(resFunc(exon.iv$span) * exon.iv$mod)
    intr.iv$normSpan <- (1 - exon.fraction) * resFunc(intr.iv$span) / sum(resFunc(intr.iv$span))
    
    #alternate intron and exon:
    
    iv <- rbind(exon.iv, intr.iv)
    iv <- iv[order(iv$start, iv$end),,drop=FALSE]
    cumulativeSum <- cumsum(iv$normSpan)
    
    iv$rescale.start <- c(0,cumulativeSum[-length(cumulativeSum)])
    iv$rescale.end <- cumulativeSum
    iv$normSpan <- iv$rescale.end - iv$rescale.start
    iv$simple.normSpan <- iv$span / sum(iv$span)
    iv$simple.rescale.start <- (iv$start - MIN) / SPAN
    iv$simple.rescale.end   <- (iv$end - MIN) / SPAN
    return(iv)
  } else { #Catch special case: the entire gene is exonic:
    exon.iv$type <- "E"
    exon.iv$normSpan <- 1 * resFunc(exon.iv$span) * exon.iv$mod / sum(resFunc(exon.iv$span)* exon.iv$mod)
    
    iv <- exon.iv
    cumulativeSum <- cumsum(iv$normSpan)
    iv$rescale.start <- c(0,cumulativeSum[-length(cumulativeSum)])
    iv$rescale.end <- cumulativeSum
    iv$normSpan <- iv$rescale.end - iv$rescale.start
    iv$simple.normSpan <- iv$span / sum(iv$span)
    iv$simple.rescale.start <- (iv$start - MIN) / SPAN
    iv$simple.rescale.end   <- (iv$end - MIN) / SPAN
    return(iv)
  }
}

rescale.coords <- function(x, rescale.iv){
  sapply(x, FUN=function(y){
    if( y < rescale.iv$start[[1]] ){
      rescale.iv$rescale.start[[1]]
    } else if( y > rescale.iv$end[ nrow(rescale.iv) ]){
      rescale.iv$rescale.end[[ nrow(rescale.iv)]]
    } else {
      idy <- which(y >= rescale.iv$start & y <= rescale.iv$end)[1]
      pct <- (y - rescale.iv$start[idy]) / rescale.iv$span[idy]
      rescale.iv$rescale.start[idy] + rescale.iv$normSpan[idy] * pct
    }
  })
}

rect.adv <- function(xleft,ybottom,xright,ytop,
           draw.border=rep(TRUE,4),
           border=rep("black",4),lty=c(1,1,1,1),lwd=c(1,1,1,1),
           ...){
  if(length(border) == 1) border = rep(border,4);
  if(length(draw.border) == 1) draw.border = rep(draw.border,4);
  if(length(lty) == 1) lty= rep(lty,4);
  if(length(lwd) == 1) lwd= rep(lwd,4);


  rect(xleft,ybottom,xright,ytop,border=NA,...);

  if(draw.border[1]) lines(c(xleft,xleft),c(ybottom,ytop),lty=lty[1],lwd=lwd[1]);
  if(draw.border[2]) lines(c(xleft,xright),c(ybottom,ybottom),lty=lty[1],lwd=lwd[1]);
  if(draw.border[3]) lines(c(xright,xright),c(ybottom,ytop),lty=lty[1],lwd=lwd[1]);
  if(draw.border[4]) lines(c(xleft,xright),c(ytop,ytop),lty=lty[1],lwd=lwd[1]);

}
beveled.rect <- function(xleft,xmid,xright,
                         ybottomleft,ybottomright,
                         ytopleft,ytopright,
                         rect.col = "gray",rect.border = "black",
                         ...){
  rect.adv(xleft,ybottomleft,xmid,ytopleft,   draw.border=c(TRUE,TRUE,FALSE,TRUE),col = rect.col,border=rect.border,...);
  rect.adv(xmid,ybottomright,xright,ytopright,draw.border=c(FALSE,TRUE,TRUE,TRUE),col = rect.col,border=rect.border,...);
  lines(c(xmid,xmid),c(ybottomright,ybottomleft),col=rect.border,...);
  lines(c(xmid,xmid),c(ytopright,ytopleft),col=rect.border,...);
}

draw.TX <- function(ybottom,ytop,exon.iv,cds.start,cds.end, rescale.iv,intron.col = NA, intron.lty = 1,
           rect.col = "gray",rect.border = "black", cex.text = 0.6,label.exons=FALSE,text.y = mean(c(ybottom[1],ytop[1])), ...){
  iv <- exon.iv[,c("start","end")];
  iv$start <-  rescale.coords(iv$start,rescale.iv);
  iv$end   <-  rescale.coords(iv$end,rescale.iv);
  cds <- c(rescale.coords(cds.start,rescale.iv), rescale.coords(cds.end,  rescale.iv))
  
  if((! is.na( intron.col)) && (nrow(iv) > 1)){
    #message("drawing introns!")
    ymid <- ybottom + (ytop - ybottom)/2
    segments( iv$end[-1],ymid,iv$end[-length(iv$end)],ymid,col = intron.col, lty=  intron.lty )
  }

  curr <- 1;
  while(cds[1] > iv$end[curr]){
    rect(iv$start[curr],ybottom[1],iv$end[curr],ytop[1],col=rect.col,border=rect.border,...);
    if(label.exons) text(mean(c(iv$start[curr],iv$end[curr])),mean(c(ybottom[1],ytop[1])),"5' UTR",cex=cex.text,...)
    curr <- curr+1;
  }
  promoterExons <- curr - 1;
  beveled.rect(iv$start[curr],cds[1],iv$end[curr],
               ybottom[1],ybottom[2],ytop[1],ytop[2],rect.col=rect.col,rect.border=rect.border,...);
    if(label.exons) text(mean(c(cds[1],iv$end[curr])),  mean(c(ybottom[1],ytop[1])),paste0("E",curr-promoterExons),cex=cex.text,...)
    if(label.exons) text(mean(c(iv$start[curr],cds[1])),mean(c(ybottom[1],ytop[1])),"5' UTR",cex=cex.text,...)

  curr <- curr + 1;
  while(cds[2] > iv$end[curr] & curr <= nrow(iv)){
    rect(iv$start[curr],ybottom[2],iv$end[curr],ytop[2],col=rect.col,border=rect.border,...);
    if(label.exons) text(mean(c(iv$start[curr],iv$end[curr])),mean(c(ybottom[1],ytop[1])),paste0("E",curr-promoterExons),cex=cex.text,...)

    curr <- curr + 1;
  }

  beveled.rect(iv$start[curr],cds[2],iv$end[curr],
               ybottom[2],ybottom[1],ytop[2],ytop[1],rect.col=rect.col,rect.border=rect.border,...);
    if(label.exons) text(mean(c(iv$start[curr],cds[2])),  mean(c(ybottom[1],ytop[1])),paste0("E",curr-promoterExons),cex=cex.text,...)
    if(label.exons) text(mean(c(cds[2],iv$end[curr])),mean(c(ybottom[1],ytop[1])),"3' UTR",cex=cex.text,...)

  curr <- curr + 1;
  while(curr <= nrow(iv)){
    rect(iv$start[curr],ybottom[1],iv$end[curr],ytop[1],col=rect.col,border=rect.border,...);
    if(label.exons) text(mean(c(iv$start[curr],iv$end[curr])),mean(c(ybottom[1],ytop[1])),"3' UTR",cex=cex.text,...)

    curr <- curr+1;
  }



}





JS.axis <- function(side, at, labels = at, tcl = -0.5, lwd = 1, lwd.ticks = 1, cex.axis = 1, srt = 0, line = NA, pos = NA, adj = c(0.5,1.1), font = 1, ...){
  
  if(side == 1){
    axis.height <- par("usr")[3]
    axis.tick.floor <- axis.height + (tcl * (2 * par("cxy")[2]))
    segments(at, axis.tick.floor,at, axis.height, lwd = lwd.ticks, cex = cex.axis, xpd = NA, ...)
    devlim <- device.limits()

    axis.label.height <- abs((axis.tick.floor - devlim[3])/2) + devlim[3]
    
    text(x = at, y = axis.label.height, labels = labels, cex = cex.axis, adj = adj, srt = srt, xpd = NA, font = font,...)
    axis.floor <- axis.label.height - max(strheight(labels, cex = cex.axis))
    
    
    return(axis.floor)
  }
}

###################################################################
###################################################################
###################################################################
################################################################


replaceCharAt <- function(s,pos,newChar){
  maxLen <- max(length(s),length(pos),length(newChar));
  if(maxLen > 1){
    if(length(s) == 1)       s <- rep(s,maxLen);
    if(length(pos) == 1)     s <- rep(pos,maxLen);
    if(length(newChar) == 1) s <- rep(newChar,maxLen);
  }
  if(length(s) != maxLen || length(pos) != maxLen || length(newChar) != maxLen){
    stop("Error: replaceCharAt(s,p,n): all inputs must be of the same length or of length 1");
  }
  anyNA <- is.na(s) | is.na(pos) | is.na(newChar);
  sapply(1:length(s),function(i){
    if(anyNA[i]){
      NA;
    } else if(pos[i] == 1){
      paste0(newChar[i],substring(s[i],2));
    } else if(pos[i] == nchar(s[i])){
      paste0(substring(s[i],1,nchar(s[i])-1),newChar[i]);
    } else {
      paste0(substring(s[i],1,pos[i]-1),newChar[i],substring(s[i],pos[i]+1));
    }
  })
}
codonSpacing <- function(fa){
  c = 1:(nchar(fa)/3)
  s = 1:(nchar(fa)/3) * 3 - 2
  e = 1:(nchar(fa)/3) * 3
  seq <- substring(rep(fa,length(c)),s,e);
  return(paste0(seq,collapse=" "));
}

AAtoAA3 <- function(x){
  ifelse(is.na(x),NA,
    ifelse(x == "*","*",AMINO_ACID_CODE[x])
  )
}

###################################################################

testing <- function(){

plot.new();
plot.window(xlim=c(0,1),ylim=c(0,1));
n <- 4;

adv.point(0.5,0.5,
  bg=c("red","blue","green","purple")[1:n],
  density = NA,
  angle = 45,
  pch = "left", cex=5,lwd=1,lty=1,border="black"
)

}


adv.point <- function(x,y,
              bg="transparent",density=NA,angle=45,
              pch="circle",cex=1,lwd=1,lty=1,border="black",
              innerLabel=NA,
              sideLabel=NA,sideLabel.pos = 4,sideLabel.cex=cex/3,
              numVert = 100){

  maxLen <- max(sapply(list(bg,density,angle),length));
  if(length(bg) < maxLen) bg <- rep(bg,length.out=maxLen);

  if(length(density) < maxLen) density<- rep(density,length.out=maxLen);
  if(length(angle) < maxLen) angle <- rep(angle,length.out=maxLen);

  #message("   maxLen :",paste0(maxLen ,collapse=","));
  #message("   bg:",paste0(bg,collapse=","));
  #message("   density:",paste0(density,collapse=","));
  #message("   angle:",paste0(angle,collapse=","));

  #Equivalent:
  #usr <- par("usr")
  #pin <- par("pin")
  ## usr coords per inch
  #upi <- c(usr[2L] - usr[1L], usr[4L] - usr[3L]) / pin
  #userScaleAspectRatio <- upi[1]/upi[2]

  Xwide    <- strwidth("o",units="user",cex=1);
  Xwide.in <- strwidth("o",units="in",cex=1);
  Xtall    <- strheight("o",units="user",cex=1);
  Xtall.in <- strheight("o",units="in",cex=1);
  userScaleAspectRatio <- (Xwide/Xwide.in) / (Xtall / Xtall.in);
  
  point.tall <- Xtall * cex;
  point.wide <- point.tall * userScaleAspectRatio
  vertices.x <- c(x- ((0:numVert)/numVert) * (point.wide/2),       x+ ((0:numVert)/numVert) * (point.wide/2));
  vertices.y <- c();

  x1 <- x + point.wide/2;
  x0 <- x - point.wide/2;
  y1 <- y + point.tall/2;
  y0 <- y - point.tall/2; 

  xm <- mean(c(x0,x1));
  ym <- mean(c(y0,y1));

  if(pch == "circle"){
      theta <- (1:(numVert*4) / (numVert*4)) * 2 * pi
      A <- abs(x1-x0)/2
      B <- abs(y1-y0)/2
      X <- A * cos(theta) + abs(x1-x0)/2 + x0
      Y <- B * sin(theta) + mean(c(y0,y1))
    if(maxLen == 1){
      polygon(X,Y,lwd=lwd,lty=lty,col=bg,density=density,angle=angle,border=border);
    } else {
      if(maxLen == 2){
        vert <- which(X <= x)
        polygon(X[vert],Y[vert],col=bg[1],density=density[1],angle=angle[1],border="transparent");
        vert <- which(X >= x)
        polygon(X[vert],Y[vert],col=bg[2],density=density[2],angle=angle[2],border="transparent");
      } else if(maxLen == 3){
        vert <- which(X <= x)
        polygon(X[vert],Y[vert],col=bg[1],density=density[1],angle=angle[1],border="transparent");
        vert <- which(X >= x & Y >= y)
        polygon(c(x,X[vert],x),c(y,Y[vert],y),col=bg[2],density=density[2],angle=angle[2],border="transparent");
        vert <- which(X >= x & Y <= y)
        polygon(c(x,X[vert],x),c(y,Y[vert],y),col=bg[3],density=density[3],angle=angle[3],border="transparent");
      } else if(maxLen == 4){
        vert <- which(X <= x & Y >= y)
        polygon(X[vert],Y[vert],col=bg[1],density=density[1],angle=angle[1],border="transparent");
        vert <- which(X <= x & Y <= y)
        polygon(c(x,X[vert],x),c(y,Y[vert],y),col=bg[2],density=density[2],angle=angle[2],border="transparent");
        vert <- which(X >= x & Y >= y)
        polygon(c(x,X[vert],x),c(y,Y[vert],y),col=bg[3],density=density[3],angle=angle[3],border="transparent");
        vert <- which(X >= x & Y <= y)
        polygon(c(x,X[vert],x),c(y,Y[vert],y),col=bg[4],density=density[4],angle=angle[4],border="transparent");
      }
      polygon(X,Y,lwd=lwd,lty=lty,col="transparent",density=NA,border=border);
    }
  } else if(pch == "up"){
    splitFrac <- 0.6;
    ym <- abs(y0-y1) * (1-splitFrac)  + y0;
    x.midWidth <- abs(x1-x0) * ((splitFrac)/2);
    xL <- xm - x.midWidth;
    xR <- xm + x.midWidth;

    if(maxLen == 1){
      polygon(c(x0,x1,xm,x0),
              c(y0,y0,y1,y0),
              lwd=lwd,lty=lty,col=bg,density=density,angle=angle,border=border);
    } else if(maxLen == 2){
      polygon(c(x0,xm,xm,x0),
              c(y0,y0,y1,y0),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(xm,x1,xm,xm),
              c(y0,y0,y1,y0),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(x0,x1,xm,x0),
              c(y0,y0,y1,y0),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 3){
      polygon(c(xL,xR,xm,xL),
              c(ym,ym,y1,ym),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(x0,xm,xm,xL,x0),
              c(y0,y0,ym,ym,y0),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(xm,x1,xR,xm,xm),
              c(y0,y0,ym,ym,y0),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(x0,x1,xm,x0),
              c(y0,y0,y1,y0),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 4){
      polygon(c(xL,xm,xm,xL),
              c(ym,ym,y1,ym),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(x0,xm,xm,xL,x0),
              c(y0,y0,ym,ym,y0),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(xm,x1,xR,xm,xm),
              c(y0,y0,ym,ym,y0),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(xm,xR,xm,xm),
              c(ym,ym,y1,ym),
              lwd=lwd,lty=lty,col=bg[4],density=density[4],angle=angle[4],border="transparent");
      polygon(c(x0,x1,xm,x0),
              c(y0,y0,y1,y0),
              lwd=lwd,lty=lty,col="transparent",border=border);
    }
  } else if(pch == "down"){
    splitFrac <- 0.6;
    ym <- abs(y0-y1) *splitFrac  + y0;
    x.midWidth <- abs(x1-x0) * ((splitFrac)/2);
    xL <- xm - x.midWidth;
    xR <- xm + x.midWidth;

    if(maxLen == 1){
      polygon(c(x0,x1,xm,x0),
              c(y1,y1,y0,y1),
              lwd=lwd,lty=lty,col=bg,density=density,angle=angle,border=border);
    } else if(maxLen == 2){
      polygon(c(x0,xm,xm,x0),
              c(y1,y1,y0,y1),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(xm,x1,xm,xm),
              c(y1,y1,y0,y1),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(x0,x1,xm,x0),
              c(y1,y1,y0,y1),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 3){
      polygon(c(xL,xR,xm,xL),
              c(ym,ym,y0,ym),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(x0,xm,xm,xL,x0),
              c(y1,y1,ym,ym,y1),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(xm,x1,xR,xm,xm),
              c(y1,y1,ym,ym,y1),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(x0,x1,xm,x0),
              c(y1,y1,y0,y1),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 4){
      polygon(c(xL,xm,xm,xL),
              c(ym,ym,y0,ym),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(x0,xm,xm,xL,x0),
              c(y1,y1,ym,ym,y1),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(xm,x1,xR,xm,xm),
              c(y1,y1,ym,ym,y1),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(xm,xR,xm,xm),
              c(ym,ym,y0,ym),
              lwd=lwd,lty=lty,col=bg[4],density=density[4],angle=angle[4],border="transparent");
      polygon(c(x0,x1,xm,x0),
              c(y1,y1,y0,y1),
              lwd=lwd,lty=lty,col="transparent",border=border);
    }
  } else if(pch == "left"){
    splitFrac <- 0.6;
    xm <- abs(x0-x1) * (splitFrac)  + x0;
    y.midHT <- abs(y1-y0) * ((splitFrac)/2);
    yD <- ym - y.midHT ;
    yU <- ym + y.midHT ;

    if(maxLen == 1){
      polygon(c(x0,x1,x1,x0),
              c(ym,y1,y0,ym),
              lwd=lwd,lty=lty,col=bg,density=density,angle=angle,border=border);
    } else if(maxLen == 2){
      polygon(c(x0,xm,xm),
              c(ym,yU,yD),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(xm,x1,x1,xm),
              c(yU,y1,y0,yD),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(x0,x1,x1,x0),
              c(ym,y1,y0,ym),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 3){
      polygon(c(x0,xm,xm),
              c(ym,yU,yD),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(xm,x1,x1,xm),
              c(yU,y1,ym,ym),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(xm,xm,x1,x1),
              c(yD,ym,ym,y0),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(x0,x1,x1,x0),
              c(ym,y1,y0,ym),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 4){
      polygon(c(x0,xm,xm),
              c(ym,yU,ym),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(x0,xm,xm),
              c(ym,yD,ym),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(xm,x1,x1,xm),
              c(yU,y1,ym,ym),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(xm,xm,x1,x1),
              c(yD,ym,ym,y0),
              lwd=lwd,lty=lty,col=bg[4],density=density[4],angle=angle[4],border="transparent");
      polygon(c(x0,x1,x1,x0),
              c(ym,y1,y0,ym),
              lwd=lwd,lty=lty,col="transparent",border=border);
    }
  } else if(pch == "right"){
    ##stop("UNKNOWN/UNSUPPORTED SHAPE!")

    splitFrac <- 0.6;
    xm <- abs(x0-x1) * (splitFrac)  + x0;
    y.midHT <- abs(y1-y0) * ((splitFrac)/2);
    yD <- ym - y.midHT ;
    yU <- ym + y.midHT ;

    if(maxLen == 1){
      polygon(c(x1,x0,x0,x1),
              c(ym,y1,y0,ym),
              lwd=lwd,lty=lty,col=bg,density=density,angle=angle,border=border);
    } else if(maxLen == 2){
      polygon(c(x1,xm,xm),
              c(ym,yU,yD),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(xm,x0,x0,xm),
              c(yU,y1,y0,yD),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(x1,x0,x0,x1),
              c(ym,y1,y0,ym),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 3){
      ##stop("UNKNOWN/UNSUPPORTED SHAPE!")
      polygon(c(x1,xm,xm),
              c(ym,yU,yD),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(xm,x0,x0,xm),
              c(yU,y1,ym,ym),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(xm,xm,x0,x0),
              c(yD,ym,ym,y0),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(x1,x0,x0,x1),
              c(ym,y1,y0,ym),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 4){
      ##stop("UNKNOWN/UNSUPPORTED SHAPE!")      
      polygon(c(x1,xm,xm),
              c(ym,yU,ym),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(x1,xm,xm),
              c(ym,yD,ym),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(xm,x0,x0,xm),
              c(yU,y1,ym,ym),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(xm,xm,x0,x0),
              c(yD,ym,ym,y0),
              lwd=lwd,lty=lty,col=bg[4],density=density[4],angle=angle[4],border="transparent");
      polygon(c(x1,x0,x0,x1),
              c(ym,y1,y0,ym),
              lwd=lwd,lty=lty,col="transparent",border=border);
    }



  } else if(pch == "dia"){
    if(maxLen == 1){
      polygon(c(x0,xm,x1,xm),
              c(ym,y1,ym,y0),
              lwd=lwd,lty=lty,col=bg,density=density,angle=angle,border=border);
    } else if(maxLen == 2){
      polygon(c(x0,xm,xm),
              c(ym,y1,y0),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(x1,xm,xm),
              c(ym,y0,y1),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(x0,xm,x1,xm),
              c(ym,y1,ym,y0),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 3){
      polygon(c(x0,xm,xm),
              c(ym,y1,y0),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(xm,x1,xm),
              c(y1,ym,ym),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(x1,xm,xm),
              c(ym,y0,ym),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(x0,xm,x1,xm),
              c(ym,y1,ym,y0),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 4){
      polygon(c(x0,xm,xm),
              c(ym,y1,ym),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(xm,x1,xm),
              c(y1,ym,ym),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(x1,xm,xm),
              c(ym,y0,ym),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(x0,xm,xm),
              c(ym,ym,y0),
              lwd=lwd,lty=lty,col=bg[4],density=density[4],angle=angle[4],border="transparent");
      polygon(c(x0,xm,x1,xm),
              c(ym,y1,ym,y0),
              lwd=lwd,lty=lty,col="transparent",border=border);
    }
  } else if(pch == "rect"){

    if(maxLen == 1){
      polygon(c(x0,x1,x1,x0),
              c(y0,y0,y1,y1),
              lwd=lwd,lty=lty,col=bg,density=density,angle=angle,border=border);
    } else if(maxLen == 2){
      polygon(c(x0,xm,xm,x0),
              c(y0,y0,y1,y1),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(xm,x1,x1,xm),
              c(y0,y0,y1,y1),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(x0,x1,x1,x0),
              c(y0,y0,y1,y1),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 3){
      polygon(c(x0,xm,xm,x0),
              c(y0,y0,y1,y1),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(xm,x1,x1,xm),
              c(ym,ym,y1,y1),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(xm,x1,x1,xm),
              c(y0,y0,ym,ym),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(x0,x1,x1,x0),
              c(y0,y0,y1,y1),
              lwd=lwd,lty=lty,col="transparent",border=border);
    } else if(maxLen == 4){
      polygon(c(x0,xm,xm,x0),
              c(y0,y0,ym,ym),
              lwd=lwd,lty=lty,col=bg[1],density=density[1],angle=angle[1],border="transparent");
      polygon(c(xm,x1,x1,xm),
              c(ym,ym,y1,y1),
              lwd=lwd,lty=lty,col=bg[2],density=density[2],angle=angle[2],border="transparent");
      polygon(c(xm,x1,x1,xm),
              c(y0,y0,ym,ym),
              lwd=lwd,lty=lty,col=bg[3],density=density[3],angle=angle[3],border="transparent");
      polygon(c(x0,xm,xm,x0),
              c(ym,ym,y1,y1),
              lwd=lwd,lty=lty,col=bg[4],density=density[4],angle=angle[4],border="transparent");
      polygon(c(x0,x1,x1,x0),
              c(y0,y0,y1,y1),
              lwd=lwd,lty=lty,col="transparent",border=border);
    }
  } else {
    stop("UNKNOWN SHAPE!")
  }
  if(!is.na(sideLabel)){
    if(sideLabel.pos == 4){
      text(x1,ym,sideLabel,cex=sideLabel.cex,adj=c(0,0.5));
    } else if(sideLabel.pos == 2){
      text(x0,ym,sideLabel,cex=sideLabel.cex,adj=c(1,0.5));
    } else{
      stop("UNSUPPORTED SIDELABEL POSITION!")
    }
  }
}


  getParam <- function(c,class.param){
    if(length(c) == 1 && c == ""){ class.param[["default"]];
    } else if(sum(c %in% names(class.param)) == 0){ class.param[["default"]];
    } else { 
       matchP <- c %in% names(class.param)
       class.param[c[matchP]];
    }
  }


advancedPoints <- function(
  x,y,pointClass,
  class.pch = c(default="circle"),
  class.bg = c(default="transparent"),
  class.cex = c(default=1),
  class.lwd = c(default=1),
  class.lty = c(default=1),
  class.density = c(default=NA),
  class.angle = c(default=45),
  class.innerLabel = c(default=NA),
  class.sideLabel = c(default=NA),
  sideLabels = NULL, innerLabels = NULL,
  sideLabel.cex = 0.5,  sideLabel.pos = 4
){

  
  if(missing(x))     stop("Parameter x is mandatory.")
  if(missing(y))     stop("Parameter y is mandatory.")
  if(missing(pointClass)) stop("Parameter pointClass is mandatory.")
  
  if(length(x) < length(y)) x <- rep(x,length.out = length(y));
  if(length(y) < length(x)) y <- rep(y,length.out = length(x));
  if(length(x) != length(pointClass)) stop("pointClass length must equal the number of points given")
  
  if( ! all( sapply(pointClass,function(pc){ length(pc) > 0 }) ) ){
    stop("All points must have a class (or \"\" for the default)");
  }  

    for(i in 1:length(pointClass)){
      pc <- pointClass[[i]];

     # message("classes:",paste0(pc ,collapse=","));

      pch <- getParam(pc,class.pch );
      bg  <- getParam(pc,class.bg  );
      cex <- getParam(pc,class.cex );
      lwd <- getParam(pc,class.lwd );
      lty <- getParam(pc,class.lty );
      density <- getParam(pc,class.density );
      angle <- getParam(pc,class.angle );
      if(is.null(innerLabels)){
        innerLabel <- getParam(pc,class.innerLabel );
      } else {
        innerLabel <- innerLabels[i];
      }
      if( is.null(sideLabels)){
        sideLabel <- getParam(pc,class.sideLabel);
      } else {
        sideLabel <- sideLabels[i];
      }
      adv.point(x[i],y[i],pch = pch, bg = bg, cex=cex,lwd=lwd,lty=lty,density=density,angle=angle,innerLabel=innerLabel,sideLabel=sideLabel,sideLabel.pos=sideLabel.pos,sideLabel.cex=sideLabel.cex);
    }
}

advancedPointsLegend.helper <- function(
    x,y,
    pch,
    labels,label.cex = 0.5,
    lineSpacing = 2,
    indent = 1,
    y.bot = NULL
  ) {
    sw = strwidth("X") 
    sh = strheight("X")
    pt.x <- x + sw * indent
    pt.y.max <- y - sh
    pt.y.min <- y - lineSpacing * sh * (length(pch) + 1)
    pt.y.span <- pt.y.max - pt.y.min
    pt.y <- pt.y.max - (((0:(length(pch)-1))/(length(pch))) * pt.y.span) 
    
    if(! is.null(y.bot)){
      pt.y.bot <- min(pt.y)
      message("attempting center :",y.bot," vs ",pt.y.bot)
      if(pt.y.bot > y.bot){
        message("attempting center...")
        pt.y.top <- max(pt.y)
        pt.y.top.space <- y - pt.y.top;
        pt.y.bot.space <- pt.y.bot - y.bot;
        space.avg <- mean(c(pt.y.bot.space,pt.y.top.space));
        if(pt.y.top.space < pt.y.bot.space){
          message("centering...")
          pt.diff <- pt.y.bot.space - space.avg;
          pt.y <- pt.y - pt.diff
        }
      }
      
    }
    
    classes <- paste0("p",1:length(pch))    
    class.pch <- as.list(pch)
    names(class.pch) <- classes
    message("pt.y = ",paste0(pt.y,collapse=","))
    message("classes = ",paste0(classes,collapse=","))
    message("class.pch = ",paste0(class.pch,collapse=","))

    advancedPoints(
      rep(pt.x,length(pch)),
      pt.y,
      pointClass = classes,
      class.pch = class.pch,
      sideLabels = paste0(" ",labels),
      sideLabel.cex = label.cex,  sideLabel.pos = 4
    )
  }


TEST.AREA <- function(){

###################################################################
###################################################################
################################################################

classList <- list(
c("S1","sarc.OOF","other.IF","del"),
c("S2","sarc.OOF","other.IF","sarc.IF","ins"),
c("S3","sarc.OOF","other.IF","sarc.IF","mel","ins"),
c("S4","sarc.OOF","ins"),
c("S5","other.IF","ins"),
c("S6","other.IF","del"),
c("S7","other.IF","sub"),
c("S8","UNKNOWN","sub"),
c("S9","UNKNOWN"),
c("")
)

plot.new();
par(mar=c(0.25,0.25,0.25,0.25))
plot.window(xlim=c(0,1),ylim=c(0,1));



advancedPoints(
  seq(from=0.25,to=0.75,length.out = length(classList)),
  rep(c(0.25,0.5,0.75),length.out=length(classList)),
  pointClass = classList,
  class.pch = c(default="circle",del="down",ins="up",sub="circle",indel="dia"),
  class.bg = c(default="gray", sarc.OOF="darkred",sarc.IF="red",other.IF="darkorchid",mel="deepskyblue"),
  class.cex = c(default=5),
  class.lwd = c(default=1),
  class.lty = c(default=1),
  class.density = c(default=20),
  class.angle = c(default=45),
  class.innerLabel = c(default=NA),
  class.sideLabel = c(default=NA,S1="1",S2="2",S3=3,S4=4,S5=5,S6=6,S7=7,S8=8,S9=9)
)


x = seq(from=0.25,to=0.75,length.out = length(classList))
y =   rep(c(0.25,0.5,0.75),length.out=length(classList))
  pointClass = classList
  class.pch = c(default="circle",del="down",ins="up",sub="circle",indel="dia")
  class.bg = c(default="gray", sarc.OOF="darkred",sarc.IF="red",other.IF="darkorchid",mel="deepskyblue")
  class.cex = c(default=5)
  class.lwd = c(default=1)
  class.lty = c(default=1)
  class.density = c(default=NA)
  class.angle = c(default=45)
  class.innerLabel = c(default=NA)
  class.sideLabel = c(default=NA)

  class.bg = c(default="gray", sarc.OOF="darkred",sarc.IF="red",other.IF="darkorchid",mel="deepskyblue")

getParam("S1",class.sideLabel)
getParam("S2",class.sideLabel)


###################################################################
###################################################################
################################################################

}
###################################################################
###################################################################
################################################################



#Modified code originally from:
#  File src/library/graphics/R/polygon.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright 1995-2016 The R Core Team
#  In part (C) 2001 Kevin Buhr
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

### polyhatch -  a pure R implementation of polygon hatching
### Copyright (C) 2001 Kevin Buhr
### Provided to the R project for release under GPL.
### Original nice clean structure destroyed by Ross Ihaka
### Further modifications by Stephen Hartley, National Cancer Institute

#mod notes:
#Added parameters to mix-and-match overlaid hatching.
#Removed support for multiple separate polygons.

polygon.mod <-
  function(x, y = NULL, density = NULL, angle = 45,
           border = NULL, col = NA, lty = par("lty"), lwd = par("lwd"),
           hatching.lty = lty, hatching.lwd = lwd, hatching.col = col, hatching.bg = "transparent",
           ..., fillOddEven=FALSE)
{
    ## FIXME: remove this eventually
    ..debug.hatch <- FALSE
    ##-- FIXME: what if `log' is active, for x or y?
    xy <- xy.coords(x, y)

    if (is.numeric(density) && all(is.na(density) | density < 0))
        density <- NULL
    if (!is.null(angle) && !is.null(density)) {

        ## hatch helper functions

        polygon.onehatch <-
            function(x, y, x0, y0, xd, yd, ..debug.hatch = FALSE, ...)
        {
            ## draw the intersection of one line with polygon
            ##
            ##  x,y - points of polygon (MUST have first and last points equal)
            ##  x0,y0 - origin of line
            ##  xd,yd - vector giving direction of line
            ##  ... - other parameters to pass to "segments"

            if (..debug.hatch) {
                points(x0, y0)
                arrows(x0, y0, x0 + xd, y0 + yd)
            }

            ## halfplane[i] is 0 or 1 as (x[i], y[i]) lies in left or right
            ##   half-plane of the line

            halfplane <- as.integer(xd * (y - y0) - yd * (x - x0) <= 0)

            ## cross[i] is -1,0, or 1 as segment (x[i], y[i]) -- (x[i+1], y[i+1])
            ##   crosses right-to-left, doesn't cross, or crosses left-to-right

            cross <- halfplane[-1L] - halfplane[-length(halfplane)]
            does.cross <- cross != 0
            if (!any(does.cross)) return() # nothing to draw?

            ## calculate where crossings occur

            x1 <- x[-length(x)][does.cross]; y1 <- y[-length(y)][does.cross]
            x2 <- x[-1L][does.cross]; y2 <- y[-1L][does.cross]

            ## t[i] is "timepoint" on line at which segment (x1, y1)--(x2, y2)
            ##   crosses such that (x0,y0) + t*(xd,yd) is point of intersection

            t <- (((x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1))/
                  (xd * (y2 - y1) - yd * (x2 - x1)))

            ## sort timepoints along line

            o <- order(t)
            tsort <- t[o]

            ## we draw the part of line from t[i] to t[i+1] whenever it lies
            ##   "inside" the polygon --- the definition of this depends on
            ##   fillOddEven:  if FALSE, we crossed
            ##   unequal numbers of left-to-right and right-to-left polygon
            ##   segments to get there.  if TRUE, an odd number of crossings.
            ##

	    crossings <- cumsum(cross[does.cross][o])
	    if (fillOddEven) crossings <- crossings %% 2
            drawline <- crossings != 0

            ## draw those segments

            lx <- x0 + xd * tsort
            ly <- y0 + yd * tsort
            lx1 <- lx[-length(lx)][drawline]; ly1 <- ly[-length(ly)][drawline]
            lx2 <- lx[-1L][drawline]; ly2 <- ly[-1L][drawline]
            segments(lx1, ly1, lx2, ly2, ...)
        }

        polygon.fullhatch <-
            function(x, y, density, angle, lty, lwd, col, ..debug.hatch = FALSE, ...)
        {
            ## draw the hatching for a given polygon
            ##
            ##  x,y - points of polygon (need not have first and last points
            ##        equal, but no NAs are allowed)
            ##  density,angle - of hatching
            ##  ... - other parameters to pass to "segments"

            x <- c(x, x[1L])
            y <- c(y, y[1L])
            angle <- angle %% 180

            if (par("xlog") || par("ylog")) {
                warning("cannot hatch with logarithmic scale active")
                return()
            }
            usr <- par("usr"); pin <- par("pin")

            ## usr coords per inch

            upi <- c(usr[2L] - usr[1L], usr[4L] - usr[3L]) / pin

            ## handle "flipped" usr coords

            if (upi[1L] < 0) angle <- 180 - angle
            if (upi[2L] < 0) angle <- 180 - angle
            upi <- abs(upi)

            ## usr-coords direction vector for hatching

            xd <- cos(angle / 180 * pi) * upi[1L]
            yd <- sin(angle / 180 * pi) * upi[2L]

            ## to generate candidate hatching lines for polygon.onehatch,
            ##   we generate those lines necessary to cover the rectangle
            ##   (min(x),min(y)) to (max(x),max(y)) depending on the
            ##   hatching angle

            ## (Note:  We choose hatch line origins such that the hatching,
            ##   if extended outside polygon, would pass through usr-coordinate
            ##   origin.  This ensures that all hatching with same density,
            ##   angle in figure will be aligned.)

            if (angle < 45 || angle > 135) {

                ## first.x and last.x are x-coords of first and last points
                ##  of rectangle to hit, as y-coord moves from bottom up

                if (angle < 45) {
                    first.x <- max(x)
                    last.x <- min(x)
                }
                else {
                    first.x <- min(x)
                    last.x <- max(x)
                }

                ## y.shift is vertical shift between parallel hatching lines

                y.shift <- upi[2L] / density / abs(cos(angle / 180 * pi))

                ## choose line origin (of first line) to align hatching
                ##   with usr origin

                x0 <- 0
                y0 <- floor((min(y) - first.x * yd / xd) / y.shift) * y.shift

                ## line origins above y.end won't hit figure

                y.end <- max(y) - last.x * yd / xd

                ## hatch against all candidate lines
                i <- 1;
                while (y0 < y.end) {
                    polygon.onehatch(x, y, x0, y0, xd, yd,
                                     ..debug.hatch=..debug.hatch,lty=lty[(i%%length(lty))+1],lwd=lwd[(i%%length(lwd))+1],col=col[(i%%length(col))+1],...)
                    y0 <- y0 + y.shift
                    i <- i+1;
                }
            }
            else {
                ## first.y, last.y are y-coords of first and last points
                ##   of rectangle to hit, as x-coord moves from left to right

                if (angle < 90) {
                    first.y <- max(y)
                    last.y <- min(y)
                }
                else {
                    first.y <- min(y)
                    last.y <- max(y)
                }

                ## x.shift is horizontal shift between parallel hatching lines

                x.shift <- upi[1L] / density / abs(sin(angle / 180 * pi))

                ## choose line origin to align with usr origin

                x0 <- floor((min(x) - first.y * xd / yd) / x.shift) * x.shift
                y0 <- 0

                ## line origins to right of x.end won't hit figure

                x.end <- max(x) - last.y * xd / yd

                ## hatch!
                i <- 1;
                while (x0 < x.end) {
                    polygon.onehatch(x, y, x0, y0, xd, yd,
                                     ..debug.hatch=..debug.hatch,lty=lty[(i%%length(lty))+1],lwd=lwd[(i%%length(lwd))+1],col=col[(i%%length(col))+1],...)
                    x0 <- x0 + x.shift;
                    i <- i + 1;
                }
            }
        }

        ## end of hatch helper functions


        if (missing(col) || is.null(col) || is.na(col)) col <- par("fg")
        if (is.null(border)) border <- col
        if (is.logical(border)) {
            if (!is.na(border) && border) border <- col
            else border <- NA
        }

        ## process multiple polygons separated by NAs

        #start <- 1
        #ends <- c(seq_along(xy$x)[is.na(xy$x) | is.na(xy$y)], length(xy$x) + 1)

        #num.polygons <- length(ends)
        #col <- rep_len(col, num.polygons)
        #if(length(border))
        #    border <- rep_len(border, num.polygons)
        #if(length(lty))
        #    lty <- rep_len(lty, num.polygons)
        #if(length(density))
        #    density <- rep_len(density, num.polygons)
        #angle <- rep_len(angle, num.polygons)

        #i <- 1L
        #.External.graphics(C_polygon, xy$x, xy$y, hatching.bg, border, lty, ...)
        polygon(xy$x,xy$y,col=hatching.bg,border=border,lty=lty,density=NA,lwd=lwd,...);

        polygon.fullhatch(xy$x,
                          xy$y,
                          col = hatching.col, lty = hatching.lty,lwd = hatching.lwd,
                          density = density,
                          angle = angle,
                          ..debug.hatch = ..debug.hatch, ...)



       # for (end in ends) {
       #     if (end > start) {
       #         if(is.null(density) || is.na(density[i]) || density[i] < 0){ 
       #             #.External(C_polygon, xy$x, xy$y, col, border, lty, ...,PACKAGE="graphics")
       #             #.External(C_polygon, xy$x[start:(end - 1)],
       #             #                   xy$y[start:(end - 1)],
       #             #                   col[i], NA, lty[i], ...,
       #             #                   PACKAGE="graphics")
       #             polygon(xy$x[start:(end - 1)],
       #                                xy$y[start:(end - 1)],
       #                                col=col[i],border=NA,lty=lty[i],lwd=lwd,...);#
#
 #                   #.External.graphics(C_polygon, xy$x[start:(end - 1)],
  #                  #                   xy$y[start:(end - 1)],
   #                 #                   col[i], NA, lty[i], ...)
    #            } else if (density[i] > 0) {
#
 #                       ## note: if col[i]==NA, "segments" will fill with par("fg")
#
 #                       polygon.fullhatch(xy$x[start:(end - 1)],
  #                                        xy$y[start:(end - 1)],
   #                                       col = hatching.col, lty = hatching.lty,lwd = hatching.lwd,
    #                                      density = density[i],
     #                                     angle = angle[i],
      #                                    ..debug.hatch = ..debug.hatch, ...)
       #             }
#
 #               ## compatible with C_polygon:
  #              ## only cycle through col, lty, etc. on non-empty polygons
   #             i <- i + 1
    #        }
    #        start <- end + 1
    #    }

    }
    else {
        if (is.logical(border)) {
            if (!is.na(border) && border) border <- par("fg")
            else border <- NA
        }
        #.External.graphics(C_polygon, xy$x, xy$y, col, border, lty, ...)
        #.External("C_polygon", xy$x, xy$y, col, border, lty, ...,PACKAGE="graphics")
        polygon(xy$x,xy$y,col=col,border=border,lty=lty,density=NA,lwd=lwd,...);
    }
    invisible()
}

###################################################################
###################################################################
###################################################################
################################################################

###################################################################
###################################################################
###################################################################
################################################################

###################################################################
###################################################################
###################################################################
################################################################

###################################################################
###################################################################
###################################################################
################################################################

###################################################################
###################################################################
###################################################################
################################################################










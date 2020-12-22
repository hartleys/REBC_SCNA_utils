setwd("C:\\Users\\hartleys\\Desktop\\BACKUP\\SCNA\\scna.github.repo\\")

#####################
# Initial Parameters:
#####################

SHAREMISC.DIR <- "./"
source( paste0( SHAREMISC.DIR , "shareMisc.Rutil/misc.R" ) )
options(stringsAsFactors=FALSE);


#####################
# Initial Parameters:
#####################

#SET THIS:
#RUNVER="v211"

input.data.file <- "../raw/REBC_primary_pair_393.alleliccapseg_pp_ccf_fit_v3_collapsed_seg.txt"

CLONALITY.THRESHOLD <- 0.75
DEBUG.MODE <- FALSE

#This file should have the following columns:
#   sample: unique sample ID
#   Chromosome
#   Start_bp
#   End_bp
#   NA
#   NB
#   NA_CCF
#   NB_CCF
#   

used.columns <- c("sample","Chromosome","Start_bp","End_bp","NA.","NB","NA_CCF","NB_CCF","IS_SCNA")



################################################################################
################################################################################
# A few more simple helper functions

is.equal <- function(x,y){
   ifelse( is.na(x), is.na(y), 
   ifelse( is.na(y), FALSE, x == y))
}

z.nan.inf <- function(x,z){
  ifelse(is.nan(x) | is.infinite(x),z,x);
}

renameDF <- function(d,x){
   for( i in seq_along(x)){
      names(d)[ names(d) == x[[i]][[1]] ] <- x[[i]][[2]]
   }
   return(d);
}

copyDF <- function(d,x){
   for( i in seq_along(x)){
      d[[ x[[i]][[2]] ]] <- d[[ x[[i]][[1]] ]]
   }
   return(d);
}

rescaleToQuantile <- function(x,p = c(0.25,0.75), ctrFunc = median){
  qt <- quantile(x , p)
  iqr <- abs(qt[[1]] - qt[[2]])
  (x - ctrFunc( x )) / iqr
}

rescaleToPhi <- function(x){
  (x - mean( x )) /  sd(x)
}
library(DESeq2)

library(data.table)
#setDT(df, keep.rownames = TRUE)[]

pd <- function(x=NULL,...,sep=",",collapse=NULL){
   if(is.null(x) || f.na(x == "") || is.na(x)){
     paste(...,sep=sep,collapse=collapse);
   } else {
     paste(x,...,sep=sep,collapse=collapse);
   }
}

################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################

BADCHROMS <- c(
 "chr13","chr14", "chr15", "chr21", "chr22"
);
#Entire SCNA has to be inside centromere
#In above 5, any intersect leads to removal

acrocentric.chroms <- c(
 "13","14", "15", "21", "22"
);
#Entire SCNA has to be inside centromere
#In above 5, any intersect leads to removal

centromere.data <- read.table("hg19.centromere.pos.txt",header=T,sep='\t');
centromere.data$chr <- gsub("chr","",centromere.data$chrom)



################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################


#####################
#####################
#####################

if(FALSE){
   RELOAD.BASE.DATA=TRUE
   SINK=FALSE
}

#####################
#####################
#####################
################################################################################


if(RELOAD.BASE.DATA){

	dda <- read.table(input.data.file,header=T,sep='\t',quote="",comment.char="");

       dda <- dda[,used.columns]

	samplist <- unique(dda$sample);
	chromlist <- sort(unique(dda$Chromosome));

	dda$event.type <- 
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 0 )   & f.na( dda$NA_CCF > 0.2 ), "deletion",
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 2 )   & f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 ), "CopyNeutralLOH",
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB >= 2 )   & f.na( dda$NB_CCF > 0.2 ), "gain",
	   ifelse(f.na(dda$NA. >= 2) & f.na( dda$NB >= 2 )   & f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 ), "gain",
	   ifelse(is.na(dda$NA.) | is.na(dda$NB) | is.na(dda$NA_CCF) | is.na(dda$NB_CCF),"NA",
	   "unknown"
	)))))

dda$event.type.v5 <- 
   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 0 )   & f.na( dda$NA_CCF > 0.2 ), "deletion",
   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 2 )   & f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 ), "CopyNeutralLOH",
   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB >= 2 )   & f.na( dda$NB_CCF > 0.2 ), "gain",
   ifelse(f.na(dda$NA. >= 2) & f.na( dda$NB >= 2 )   & f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 ), "gain",
   ifelse(is.na(dda$NA.) | is.na(dda$NB) | is.na(dda$NA_CCF) | is.na(dda$NB_CCF),"unknown",
   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 0 )   & f.na( dda$NA_CCF < 0.2 ) & f.na( dda$NB_CCF < 0.2 ), "nonEvent",
   "unknown"
))))))

dda$event.review.flag <- 
   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 0 )   & f.na( dda$NA_CCF > 0.2 ), F,
   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 2 )   & f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 ), F,
   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB >= 2 )   & f.na( dda$NB_CCF > 0.2 ), F,
   ifelse(f.na(dda$NA. >= 2) & f.na( dda$NB >= 2 )   & f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 ), F,
   ifelse(is.na(dda$NA.) | is.na(dda$NB) | is.na(dda$NA_CCF) | is.na(dda$NB_CCF),F,
   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 0 )   & f.na( dda$NA_CCF < 0.2 ) & f.na( dda$NB_CCF < 0.2 ),T,
   T
))))))


	dda$event.clonality.09 <- 
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 0 )   & f.na( dda$NA_CCF > 0.9 ), "clonal",
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 2 )   & (f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 )) & (f.na( dda$NA_CCF > 0.9 ) | f.na( dda$NB_CCF > 0.9 )), "clonal",
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB >= 2 )   & f.na( dda$NB_CCF > CLONALITY.THRESHOLD ), "clonal",
	   ifelse(f.na(dda$NA. >= 2) & f.na( dda$NB >= 2 )   & (f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 )) & (f.na( dda$NA_CCF > 0.9 ) | f.na( dda$NB_CCF > 0.9 )), "clonal",
	   ""
	))))

	dda$event.clonality <- 
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 0 )   & f.na( dda$NA_CCF > CLONALITY.THRESHOLD ), "clonal",
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 2 )   & (f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 )) & (f.na( dda$NA_CCF > CLONALITY.THRESHOLD ) | f.na( dda$NB_CCF > CLONALITY.THRESHOLD )), "clonal",
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB >= 2 )   & f.na( dda$NB_CCF > CLONALITY.THRESHOLD ), "clonal",
	   ifelse(f.na(dda$NA. >= 2) & f.na( dda$NB >= 2 )   & (f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 )) & (f.na( dda$NA_CCF > CLONALITY.THRESHOLD ) | f.na( dda$NB_CCF > CLONALITY.THRESHOLD )), "clonal",
	   ""
	))))

	dda$NXCCF <- 
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 0 )   & f.na( dda$NA_CCF > 0.2 ), dda$NA_CCF,
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB == 2 )   & (f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 )), pmax(dda$NA_CCF,dda$NB_CCF ,na.rm=T),
	   ifelse(f.na(dda$NA. == 0) & f.na( dda$NB >= 2 )   & f.na( dda$NB_CCF > 0.2 ), dda$NB_CCF,
	   ifelse(f.na(dda$NA. >= 2) & f.na( dda$NB >= 2 )   & (f.na( dda$NA_CCF > 0.2 ) & f.na( dda$NB_CCF > 0.2 )), pmax(dda$NA_CCF,dda$NB_CCF ,na.rm=T),
	   pmax(dda$NA_CCF,dda$NB_CCF ,na.rm=T)
	))))

	dda$event.info      <- ifelse(dda$event.clonality=="",dda$event.type,p(dda$event.clonality,".",dda$event.type))

}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

dda$Nstat <- p(dda$NA.,"/",dda$NB,"/",dda$NA_CCF,"/",dda$NB_CCF)
dda$rownum <- 1:nrow(dda);
dda$SU <- p(dda$sample,"/",dda$Chromosome);
dda$notes <- ""

stats <- c("rawCT","chrCT","mrgCT");
eventlist <- c("deletion","gain","CopyNeutralLOH")
fulllist <- c("deletion","gain","CopyNeutralLOH","clonal.deletion","clonal.gain","clonal.CopyNeutralLOH")
BUFFER.WINDOW <- 1000000000;

anytype.multi.on.chrom.ct = 0;
multi.on.chrom.ct = 0;

multidata.info <- list();

sampCounts <- list()
for( ss in samplist ){
  message("-",ss);
  sampCounts[[ss]] <- list();
  kk <- dda$sample == ss & dda$IS_SCNA == 1;
  ddax <- dda[kk,]
        kkn <-  dda$sample == ss & (dda$IS_SCNA == 1 | dda$event.type == "NA" | dda$event.type == "unknown");
        ddaxn <- dda[kkn,]
  multidata.info[[ p(ss) ]] <- list()
  multidata.info[[ p(ss) ]][["data"]] <- ddaxn;
  multidata.info[[ p(ss) ]][["events"]] <- character();
  multidata.info[[ p(ss) ]][["events.clonal"]] <- character();
  multidata.info[[ p(ss) ]][["events.info"]] <- list();
  multidata.info[[ p(ss) ]][["events.idx"]]  <- list();
  multidata.info[[ p(ss) ]][["events.data"]] <- list();
  multidata.info[[ p(ss) ]][["fullchrom"]] <- list();
  multidata.info[[ p(ss) ]][["flags"]] <- c();

    sampCounts[[ss]][[p("rawCT.ALL")]] <- sum(kk);
    sampCounts[[ss]][[p("chrCT.ALL")]] <- 0;
    sampCounts[[ss]][[p("mrgCT.ALL")]] <- 0;
    for(cc in chromlist){
      kkx  <- kk  & dda$Chromosome == cc
      kkxc <- kkx & dda$event.clonality == "clonal"
      if(sum(kkx) > 0){
        sampCounts[[ss]][[p("chrCT.ALL")]] <- sampCounts[[ss]][[p("chrCT.ALL")]]+1;
      }
    }

  for(ee in eventlist ){
    kk <- dda$sample == ss & dda$event.type == ee & dda$IS_SCNA == 1;

    sampCounts[[ss]][[p("rawCT.",ee)]] <- sum(kk);
    sampCounts[[ss]][[p("chrCT.",ee)]] <- 0;
    sampCounts[[ss]][[p("mrgCT.",ee)]] <- 0;
    sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- 0;
    for(cc in chromlist){
      kkx  <- kk  & dda$Chromosome == cc
      kkxc <- kkx & dda$event.clonality == "clonal"
      if(sum(kkx) > 0){
        sampCounts[[ss]][[p("chrCT.",ee)]] <- sampCounts[[ss]][[p("chrCT.",ee)]]+1;
        idx <- which(kkx | (dda$sample == ss & (dda$event.info == "NA" | dda$event.type == "unknown") & dda$Chromosome == cc))
        while( length(idx) > 0 && ("NA" == ( dda$event.type[idx[1]] ) | "unknown" == ( dda$event.type[idx[1]] ) )     ){
          idx <- idx[-1];
        }
        while( length(idx) > 0 && ("NA" == ( dda$event.type[[ idx[[length(idx)]] ]] ) | "unknown" == ( dda$event.type[[ idx[[length(idx)]] ]] )  )){
          idx <- idx[-length(idx)];
        }
        if(length(idx) == 1){
          sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]+1;
          if( sum(kkxc) > 0){
            sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- sampCounts[[ss]][[p("mrgCT.clonal.",ee)]]+1;
            multidata.info[[ p(ss) ]][["events.clonal"]] <- c(multidata.info[[ p(ss) ]][["events.clonal"]],p(ee,":chr",cc))

          }
          is.clonal <- sum(kkxc) > 0
          clonality.string <- ifelse(sum(kkxc) > 0,"clonal","nonclonal")

            ddax <- multidata.info[[ p(ss) ]][["data"]]
            zzix <- which( ddax$rownum == dda$rownum[ idx[1] ])
            ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("SINGLETON[",ee,"]"))
            multidata.info[[ p(ss) ]][["data"]] <- ddax;
            multidata.info[[ p(ss) ]][["events"]] <- c(multidata.info[[ p(ss) ]][["events"]],p(ee,":chr",cc))
                  multidata.info[[ p(ss) ]][["events.info"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- list(
                      event = ee,chrom=cc,
                      start = ddax$Start_bp[[ zzix ]],
                      end   = ddax$End_bp[[zzix]],
                      clonal = is.clonal,
                      idx    = p(zzix,collapse=","),
                      N = 1
                  )
                  #message("class: ",class( multidata.info[[ p(ss) ]][["events.idx"]] ))
                  multidata.info[[ p(ss) ]][["events.idx"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- c(zzix);
                  multidata.info[[ p(ss) ]][["events.data"]][[ length(multidata.info[[ p(ss) ]][["events.data"]])+1  ]] <- ddax[zzix,,drop=F]
        } else if(length(idx) > 1){
          ddax <- multidata.info[[ p(ss) ]][["data"]]
          message("Multi On Chrom: ",ss,"/",cc);
          multi.on.chrom.ct = multi.on.chrom.ct  + 1;
          sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]+1;
          zzix <- which( ddax$rownum == dda$rownum[ idx[1] ])
          ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("START[",ee,"]"))
          if(DEBUG.MODE){ message("notes[i=",1,"]/[idx=",idx[1],"]/[ddax=",(zzix),"]=END") }
          multidata.info[[ p(ss) ]][["events"]] <- c(multidata.info[[ p(ss) ]][["events"]],p(ee,":chr",cc))
          #multidata.info[[ p(ss) ]][["events.info"]] <- c(multidata.info[[ p(ss) ]][["events.info"]],p(ee,":chr",cc,":", ddax$Start_bp[[zzix]],":"  ))

          for(i in seq_along(idx)[-1]){
            if(DEBUG.MODE){ message("i=",i) }
            if(DEBUG.MODE){ print(dda[idx[(i-1):(i)],c("sample","Chromosome","Start_bp","End_bp","event.type","rownum","notes")] )}
            if( dda$End_bp[idx[i-1]] + BUFFER.WINDOW < dda$Start_bp[idx[i]] || dda$rownum[idx[i-1]] + 1 < dda$rownum[idx[i]] ){
              if(DEBUG.MODE){ message("BREAK-step 1") }
              zzix <- which( ddax$rownum == dda$rownum[ idx[i-1] ])
              ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("END[",ee,"]"))
              if(DEBUG.MODE){ message("notes[i=",i-1,"]/[idx=",idx[i-1],"]/[ddax=",(zzix),"]=END") }

              idx.prev   <- 1:(i-1);
              zzix.prev  <- which( ddax$rownum %in% dda$rownum[ idx[idx.prev] ] )
              notes.prev <- ddax$notes[ zzix.prev ];
              ixx.start  <- max( which(grepl(p("START[",ee,"]"),notes.prev,fixed=T) ))
              idx.curr   <- ixx.start:(i-1);
              zzix.curr  <- which( ddax$rownum %in% dda$rownum[ idx[idx.curr] ])
              zzix.last  <- zzix.curr[[length(zzix.curr)]]
              is.clonal  <- any( f.na( ddax$event.clonality[zzix.curr] == "clonal" ))
              is.event   <- any( f.na( ddax$event.type[zzix.curr] == ee ))
              clonality.string <- ifelse(is.clonal,"clonal","nonclonal")
              if(is.clonal & is.event){
                sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- sampCounts[[ss]][[p("mrgCT.clonal.",ee)]]+1;
                multidata.info[[ p(ss) ]][["events.clonal"]] <- c(multidata.info[[ p(ss) ]][["events.clonal"]],p(ee,":chr",cc))
              }
              if(! is.event ){
                message("DELETING EVENT!")
                sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]-1;
                multidata.info[[ p(ss) ]][["events"]] <- multidata.info[[ p(ss) ]][["events"]][- length(multidata.info[[ p(ss) ]][["events"]])]
              } else {
                  multidata.info[[ p(ss) ]][["events.info"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- list(
                      event = ee,chrom=cc,
                      start = ddax$Start_bp[[ zzix.curr[[1]] ]],
                      end   = ddax$End_bp[[zzix.last]],
                      clonal = is.clonal,
                      idx    = p(zzix.curr,collapse=","),
                      N = length(zzix.curr)
                  )
                #multidata.info[[ p(ss) ]][["events.info"]][[length(multidata.info[[ p(ss) ]][["events.info"]])]] <- 
                #            p(multidata.info[[ p(ss) ]][["events.info"]][[length(multidata.info[[ p(ss) ]][["events.info"]])]],ddax$End_bp[[zzix.last]],":",clonality.string)
                  multidata.info[[ p(ss) ]][["events.idx"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- zzix.curr;
                  multidata.info[[ p(ss) ]][["events.data"]][[ length(multidata.info[[ p(ss) ]][["events.data"]])+1  ]] <- ddax[zzix.curr,,drop=F]
              }

              if(any( dda$event.type[idx[i:length(i)]] == ee )){
                if(DEBUG.MODE){ message("BREAK-step 2") }
                sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]+1;
                zzix <- which( ddax$rownum == dda$rownum[ idx[i] ])
                ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("START[",ee,"]"))
                if(DEBUG.MODE){ message("notes[i=",i,"]/[idx=",idx[i],"]/[ddax=",(zzix),"]=START") }
                multidata.info[[ p(ss) ]][["events"]] <- c(multidata.info[[ p(ss) ]][["events"]],p(ee,":chr",cc))
              }
            }  
          }
          zzix <- which( ddax$rownum == dda$rownum[ idx[length(idx)] ])
          ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("END[",ee,"]"))

              idx.prev   <- 1:length(idx);
              zzix.prev  <- which( ddax$rownum %in% dda$rownum[ idx[idx.prev] ] )
              notes.prev <- ddax$notes[ zzix.prev ];
              ixx.start  <- max( which(grepl(p("START[",ee,"]"),notes.prev,fixed=T) ))
              idx.curr   <- ixx.start:length(idx);
              zzix.curr  <- which( ddax$rownum %in% dda$rownum[ idx[idx.curr] ])
              zzix.last  <- zzix.curr[[length(zzix.curr)]]
              is.clonal  <- any( f.na( ddax$event.clonality[zzix.curr] == "clonal" ))
              is.event   <- any( f.na( ddax$event.type[zzix.curr] == ee ))
              clonality.string <- ifelse(is.clonal,"clonal","nonclonal")
              if(is.clonal){
                sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- sampCounts[[ss]][[p("mrgCT.clonal.",ee)]]+1;
                multidata.info[[ p(ss) ]][["events.clonal"]] <- c(multidata.info[[ p(ss) ]][["events.clonal"]],p(ee,":chr",cc))
              }
              if(is.event){
                  multidata.info[[ p(ss) ]][["events.info"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- list(
                      event = ee,chrom=cc,
                      start = ddax$Start_bp[[ zzix.curr[[1]] ]],
                      end   = ddax$End_bp[[zzix.last]],
                      clonal = is.clonal,
                      idx    = p(zzix.curr,collapse=","),
                      N = length(zzix.curr)
                  )
                  multidata.info[[ p(ss) ]][["events.idx"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- zzix.curr;
                  multidata.info[[ p(ss) ]][["events.data"]][[ length(multidata.info[[ p(ss) ]][["events.data"]])+1  ]] <- ddax[zzix.curr,,drop=F]
              }
          if(DEBUG.MODE){ message("notes[i=",length(idx),"]/[idx=",idx[length(idx)],"]/[ddax=",(zzix),"]=END") }
          multidata.info[[ p(ss) ]][["data"]] <- ddax;
        }
      }
    }
  }
  sampCounts[[ss]][[p("mrgCT.ALL")]] <- 0;
  for(ee in eventlist ){
    sampCounts[[ss]][[p("mrgCT.ALL")]] <- sampCounts[[ss]][[p("mrgCT.ALL")]] + sampCounts[[ss]][[p("mrgCT.",ee)]]
  }
}

sampCountData <- do.call(rbind.data.frame,sampCounts);
setDT(sampCountData, keep.rownames = TRUE)[]
names(sampCountData)[[1]] <- "sample";
sampCountData.LABELED_REBC_ID <- data.frame( LABELED_REBC_ID = sapply(strsplit(sampCountData$sample,"-"),function(ss){
  p(ss[[1]],"-",ss[[2]]);
}))
sampCountData <- cbind.data.frame(sampCountData.LABELED_REBC_ID,sampCountData);





zzlist <- grepv("mrgCT.",names(sampCounts[[1]]))
zzlist <- zzlist[zzlist != "mrgCT.ALL"]
for( ss in names(multidata.info) ){
   ddax <- multidata.info[[ p(ss) ]][["data"]]
   ddad <- data.frame(SCNA.type = character(), CT = numeric());
   ddsc <- sampCounts[[ss]];
   #message("ss = ",ss);
   if(! is.null(ddax)){
      for(zz in zzlist){
        if( ddsc[[zz]] > 0 ){
          ddad <- rbind(ddad,data.frame(SCNA.type = zz, CT = ddsc[[zz]]));
        }
      }
   }
   multidata.info[[ p(ss) ]][["result"]] <- ddad
}


for( ss in samplist ){

   if(is.null(multidata.info[[ p(ss) ]])){
     ddad <- data.frame(SCNA.type = character(), CT = numeric());
     ddsc <- sampCounts[[ss]];
     multidata.info[[ss]] <- list();
     multidata.info[[ss]][["data"]] <- 
     for(zz in zzlist){
        if( ddsc[[zz]] > 0 ){
          ddad <- rbind(ddad,data.frame(SCNA.type = zz, CT = ddsc[[zz]]));
        }
     }
     kkn <-  dda$sample == ss & (dda$IS_SCNA == 1 | dda$event.type == "NA" | dda$event.type == "unknown" )
     ddaxn <- dda[kkn,]
     multidata.info[[ p(ss) ]][["data"]] <- ddaxn;
     multidata.info[[ p(ss) ]][["result"]] <- ddad
   }

}

multidata.info.B <- multidata.info
sampCounts.B <- sampCounts
sampCountData.B <- sampCountData 

multidata.info.B.events.len <- sapply(sapply(multidata.info.B,'[[',"events"),length);
multidata.info.B.eventsINFO.len <- sapply(sapply(multidata.info.B,'[[',"events.info"),length);
multidata.info.B.eventsIDX.len <- sapply(sapply(multidata.info.B,'[[',"events.idx"),length);
multidata.info.B.eventsDATA.len <- sapply(sapply(multidata.info.B,'[[',"events.data"),length);


table(multidata.info.B.events.len,multidata.info.B.eventsINFO.len)
table(multidata.info.B.events.len,multidata.info.B.eventsINFO.len)
table(multidata.info.B.events.len,multidata.info.B.eventsIDX.len)
table(multidata.info.B.events.len,multidata.info.B.eventsDATA.len)

which(multidata.info.B.events.len !=  multidata.info.B.eventsINFO.len)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


#RUNVER="v070"

dda$rownum <- 1:nrow(dda);
dda$SU <- p(dda$sample,"/",dda$Chromosome);
dda$notes <- ""

stats <- c("rawCT","chrCT","mrgCT");
eventlist <- c("deletion","gain","CopyNeutralLOH")
fulllist <- c("deletion","gain","CopyNeutralLOH",
                     "clonal.deletion","clonal.gain","clonal.CopyNeutralLOH")
BUFFER.WINDOW <- 1000000000;

anytype.multi.on.chrom.ct = 0;
multi.on.chrom.ct = 0;

multidata.info <- list();



sampCounts <- list()
for( ss in samplist ){
  message("-",ss);
  sampCounts[[ss]] <- list();
  kk <- dda$sample == ss & dda$IS_SCNA == 1;
  kkn <-  dda$sample == ss & (dda$IS_SCNA == 1 | dda$event.type == "NA");
  ddax <- dda[kk,]
  ddaxn <- dda[kkn,]
  multidata.info[[ p(ss) ]] <- list()
  multidata.info[[ p(ss) ]][["data"]] <- ddaxn;
  multidata.info[[ p(ss) ]][["events"]] <- character();
  multidata.info[[ p(ss) ]][["fullchrom"]] <- list();
  multidata.info[[ p(ss) ]][["fullchrom.FULL"]] <- list();
  multidata.info[[ p(ss) ]][["flags"]] <- c();
  multidata.info[[ p(ss) ]][["events.clonal"]] <- character();
  multidata.info[[ p(ss) ]][["events.info"]] <- list();
  multidata.info[[ p(ss) ]][["events.idx"]]  <- list();
  multidata.info[[ p(ss) ]][["events.data"]] <- list();



    sampCounts[[ss]][[p("rawCT.ALL")]] <- sum(kk);
    sampCounts[[ss]][[p("chrCT.ALL")]] <- 0;
    #sampCounts[[ss]][[p("mrgCT.ALL")]] <- 0;
    for(cc in chromlist){
      kkx  <- kk  & dda$Chromosome == cc
      kkxc <- kkx & dda$event.clonality == "clonal"
      if(sum(kkx) > 0){
        sampCounts[[ss]][[p("chrCT.ALL")]] <- sampCounts[[ss]][[p("chrCT.ALL")]]+1;
      }
    }

    for(cc in chromlist){
      kk <- dda$sample == ss & dda$Chromosome == cc
      if(cc %in% acrocentric.chroms){
        centromere.end <- centromere.data$chromEnd[ centromere.data$chr == cc ]
        kk <- kk & dda$End_bp >= centromere.end
      }
      event.listing <- unique(dda$event.type.v5[kk])
      event.listing <- event.listing[ ! event.listing %in% c("NA") ]
      simple.ee <- event.listing;
      if(length(event.listing) == 1 && event.listing %in% c("deletion","gain","CopyNeutralLOH")){
        if( any(dda$event.clonality[kk] == "clonal")){
           event.listing <- p("clonal.",event.listing);
        }
        multidata.info[[ p(ss) ]][["fullchrom"]][[ length(multidata.info[[ p(ss) ]][["fullchrom"]]) + 1 ]] <-
            list(chrom=cc,ee=event.listing,event=simple.ee,flag=FALSE)
      } else {
        event.listing <- event.listing[ ! event.listing %in% c("unknown") ]
        if(length(event.listing) == 1 && event.listing %in% c("deletion","gain","CopyNeutralLOH")){
          if( any(dda$event.clonality[kk] == "clonal")){
             event.listing <- p("clonal.",event.listing);
          }
          multidata.info[[ p(ss) ]][["fullchrom"]][[ length(multidata.info[[ p(ss) ]][["fullchrom"]]) + 1 ]] <-
              list(chrom=cc,ee=event.listing,event=simple.ee,flag=TRUE) 
        }
      }
    }
    for(cc in chromlist){
      kk <- dda$sample == ss & dda$Chromosome == cc
      event.listing <- unique(dda$event.type.v5[kk])
      event.listing <- event.listing[ ! event.listing %in% c("NA") ]
      simple.ee <- event.listing;
      if(length(event.listing) == 1 && event.listing %in% c("deletion","gain","CopyNeutralLOH")){
        if( any(dda$event.clonality[kk] == "clonal")){
           event.listing <- p("clonal.",event.listing);
        }
        multidata.info[[ p(ss) ]][["fullchrom.FULL"]][[ length(multidata.info[[ p(ss) ]][["fullchrom.FULL"]]) + 1 ]] <-
            list(chrom=cc,ee=event.listing,event=simple.ee,flag=FALSE)
      } else {
        event.listing <- event.listing[ ! event.listing %in% c("unknown") ]
        if(length(event.listing) == 1 && event.listing %in% c("deletion","gain","CopyNeutralLOH")){
          if( any(dda$event.clonality[kk] == "clonal")){
             event.listing <- p("clonal.",event.listing);
          }
          multidata.info[[ p(ss) ]][["fullchrom.FULL"]][[ length(multidata.info[[ p(ss) ]][["fullchrom.FULL"]]) + 1 ]] <-
              list(chrom=cc,ee=event.listing,event=simple.ee,flag=TRUE) 
        }
      }
    }


  for(ee in eventlist ){
    kk <- dda$sample == ss & dda$event.type == ee & dda$IS_SCNA == 1;

    sampCounts[[ss]][[p("rawCT.",ee)]] <- sum(kk);
    sampCounts[[ss]][[p("chrCT.",ee)]] <- 0;
    sampCounts[[ss]][[p("mrgCT.",ee)]] <- 0;
    #sampCounts[[ss]][[p("rawCT.clonal.",ee)]] <- sum(kk);
    #sampCounts[[ss]][[p("chrCT.clonal.",ee)]] <- 0;
    sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- 0;
    for(cc in chromlist){
      kkx  <- kk  & dda$Chromosome == cc
      kkxc <- kkx & dda$event.clonality == "clonal"
      if(sum(kkx) > 0){
        sampCounts[[ss]][[p("chrCT.",ee)]] <- sampCounts[[ss]][[p("chrCT.",ee)]]+1;
        idx <- which(kkx | (dda$sample == ss & dda$event.info == "NA" & dda$Chromosome == cc))
        while( length(idx) > 0 && "NA" == ( dda$event.type[[ idx[[1]] ]] )){
          idx <- idx[-1];
        }
        while( length(idx) > 0 && "NA" == ( dda$event.type[[ idx[[length(idx)]] ]] )){
          idx <- idx[-length(idx)];
        }


        if(length(idx) == 1){
          sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]+1;
          if( sum(kkxc) > 0){
            sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- sampCounts[[ss]][[p("mrgCT.clonal.",ee)]]+1;
            multidata.info[[ p(ss) ]][["events.clonal"]] <- c(multidata.info[[ p(ss) ]][["events.clonal"]],p(ee,":chr",cc))
          }
          is.clonal <- sum(kkxc) > 0
          clonality.string <- ifelse(sum(kkxc) > 0,"clonal","nonclonal")

          ddax <- multidata.info[[ p(ss) ]][["data"]]
          zzix <- which( ddax$rownum == dda$rownum[ idx[1] ])
          ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("SINGLETON[",ee,"]"))
          multidata.info[[ p(ss) ]][["data"]] <- ddax;
          multidata.info[[ p(ss) ]][["events"]] <- c(multidata.info[[ p(ss) ]][["events"]],p(ee,":chr",cc))
                  multidata.info[[ p(ss) ]][["events.info"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- list(
                      event = ee,chrom=cc,
                      start = ddax$Start_bp[[ zzix ]],
                      end   = ddax$End_bp[[zzix]],
                      clonal = is.clonal,
                      idx    = p(zzix,collapse=","),
                      N = 1
                  )
                  #message("class: ",class( multidata.info[[ p(ss) ]][["events.idx"]] ))
                  multidata.info[[ p(ss) ]][["events.idx"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- c(zzix);
                  multidata.info[[ p(ss) ]][["events.data"]][[ length(multidata.info[[ p(ss) ]][["events.data"]])+1  ]] <- ddax[zzix,,drop=F]


        } else if(length(idx) > 1){
          ddax <- multidata.info[[ p(ss) ]][["data"]]
          message("Multi On Chrom: ",ss,"/",cc);
          multi.on.chrom.ct = multi.on.chrom.ct  + 1;
          sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]+1;
          multidata.info[[ p(ss) ]][["events"]] <- c(multidata.info[[ p(ss) ]][["events"]],p(ee,":chr",cc))
              message("----------");
              message("mrgCT.",ee," = ", sampCounts[[ss]][[p("mrgCT.",ee)]] )
              print( multidata.info[[ p(ss) ]][["events"]] )
          zzix <- which( ddax$rownum == dda$rownum[ idx[1] ])
          ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("START[",ee,"]"))
          if(DEBUG.MODE){ message("notes[i=",1,"]/[idx=",idx[1],"]/[ddax=",(zzix),"]=END") }
          for(i in seq_along(idx)[-1]){
            if(DEBUG.MODE){ message("i=",i) }
            if(DEBUG.MODE){ print(dda[idx[(i-1):(i)],c("sample","Chromosome","Start_bp","End_bp","event.type","rownum","notes")] )}
            if( dda$End_bp[idx[i-1]] + BUFFER.WINDOW < dda$Start_bp[idx[i]] || dda$rownum[idx[i-1]] + 1 < dda$rownum[idx[i]] ){
              if(DEBUG.MODE){ message("BREAK-step 1") }
              zzix <- which( ddax$rownum == dda$rownum[ idx[i-1] ])
              ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("END[",ee,"]"))
              if(DEBUG.MODE){ message("notes[i=",i-1,"]/[idx=",idx[i-1],"]/[ddax=",(zzix),"]=END") }

              idx.prev   <- 1:(i-1);
              zzix.prev  <- which( ddax$rownum %in% dda$rownum[ idx[idx.prev] ] )
              notes.prev <- ddax$notes[ zzix.prev ];
              ixx.start  <- max( which(grepl(p("START[",ee,"]"),notes.prev,fixed=T) ))
              idx.curr   <- ixx.start:(i-1);
              zzix.curr  <- which( ddax$rownum %in% dda$rownum[ idx[idx.curr] ])
              zzix.last  <- zzix.curr[[length(zzix.curr)]]
              is.clonal  <- any( f.na( ddax$event.clonality[zzix.curr] == "clonal" ))
              is.event   <- any( f.na( ddax$event.type[zzix.curr] == ee ))
          clonality.string <- ifelse(sum(kkxc) > 0,"clonal","nonclonal")

              if(is.clonal & is.event){
                sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- sampCounts[[ss]][[p("mrgCT.clonal.",ee)]]+1;
                multidata.info[[ p(ss) ]][["events.clonal"]] <- c(multidata.info[[ p(ss) ]][["events.clonal"]],p(ee,":chr",cc))
              }
              if(! is.event ){
                message("----------- NON EVENT: ",ss," / ",cc ," / ",ee);
                print(ddax[zzix.curr,])
                sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]-1;
                multidata.info[[ p(ss) ]][["events"]] <- multidata.info[[ p(ss) ]][["events"]][- length(multidata.info[[ p(ss) ]][["events"]])]
              } else {
                  multidata.info[[ p(ss) ]][["events.info"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- list(
                      event = ee,chrom=cc,
                      start = ddax$Start_bp[[ zzix.curr[[1]] ]],
                      end   = ddax$End_bp[[zzix.last]],
                      clonal = is.clonal,
                      idx    = p(zzix.curr,collapse=","),
                      N = length(zzix.curr)
                  )
                #multidata.info[[ p(ss) ]][["events.info"]][[length(multidata.info[[ p(ss) ]][["events.info"]])]] <- 
                #            p(multidata.info[[ p(ss) ]][["events.info"]][[length(multidata.info[[ p(ss) ]][["events.info"]])]],ddax$End_bp[[zzix.last]],":",clonality.string)
                  multidata.info[[ p(ss) ]][["events.idx"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- zzix.curr;
                  multidata.info[[ p(ss) ]][["events.data"]][[ length(multidata.info[[ p(ss) ]][["events.data"]])+1  ]] <- ddax[zzix.curr,,drop=F]
              }
              if(any( dda$event.type[idx[i:length(idx)]] == ee )){
                if(DEBUG.MODE){ message("BREAK-step 2") }
                sampCounts[[ss]][[p("mrgCT.",ee)]] <- sampCounts[[ss]][[p("mrgCT.",ee)]]+1;
                multidata.info[[ p(ss) ]][["events"]] <- c(multidata.info[[ p(ss) ]][["events"]],p(ee,":chr",cc))
                  message("----------");
                  message("mrgCT.",ee," = ", sampCounts[[ss]][[p("mrgCT.",ee)]] )
                  print( multidata.info[[ p(ss) ]][["events"]] )
  
                zzix <- which( ddax$rownum == dda$rownum[ idx[i] ])
                ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("START[",ee,"]"))
                if(DEBUG.MODE){ message("notes[i=",i,"]/[idx=",idx[i],"]/[ddax=",(zzix),"]=START") }
              }
            }  
          }
          zzix <- which( ddax$rownum == dda$rownum[ idx[length(idx)] ])
          ddax$notes[ zzix ] <- pd(ddax$notes[ zzix ],p("END[",ee,"]"))

              idx.prev   <- 1:length(idx);
              zzix.prev  <- which( ddax$rownum %in% dda$rownum[ idx[idx.prev] ] )
              notes.prev <- ddax$notes[ zzix.prev ];
              ixx.start  <- max( which(grepl(p("START[",ee,"]"),notes.prev,fixed=T) ))
              idx.curr   <- ixx.start:length(idx);
              zzix.curr  <- which( ddax$rownum %in% dda$rownum[ idx[idx.curr] ])
              zzix.last  <- zzix.curr[[length(zzix.curr)]]
              is.clonal  <- any( f.na( ddax$event.clonality[zzix.curr] == "clonal" ))
              clonality.string <- ifelse(is.clonal,"clonal","nonclonal")
              is.event   <- any( f.na( ddax$event.type[zzix.curr] == ee ))

              if(is.clonal){
                sampCounts[[ss]][[p("mrgCT.clonal.",ee)]] <- sampCounts[[ss]][[p("mrgCT.clonal.",ee)]]+1;
                multidata.info[[ p(ss) ]][["events.clonal"]] <- c(multidata.info[[ p(ss) ]][["events.clonal"]],p(ee,":chr",cc))
              }
              if(is.event){
                  multidata.info[[ p(ss) ]][["events.info"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- list(
                      event = ee,chrom=cc,
                      start = ddax$Start_bp[[ zzix.curr[[1]] ]],
                      end   = ddax$End_bp[[zzix.last]],
                      clonal = is.clonal,
                      idx    = p(zzix.curr,collapse=","),
                      N = length(zzix.curr)
                  )
                #multidata.info[[ p(ss) ]][["events.info"]][[length(multidata.info[[ p(ss) ]][["events.info"]])]] <- 
                #            p(multidata.info[[ p(ss) ]][["events.info"]][[length(multidata.info[[ p(ss) ]][["events.info"]])]],ddax$End_bp[[zzix.last]],":",clonality.string)
                  multidata.info[[ p(ss) ]][["events.idx"]][[ length(multidata.info[[ p(ss) ]][["events.idx"]])+1  ]] <- zzix.curr;
                  multidata.info[[ p(ss) ]][["events.data"]][[ length(multidata.info[[ p(ss) ]][["events.data"]])+1  ]] <- ddax[zzix.curr,,drop=F]
              }
          if(DEBUG.MODE){ message("notes[i=",length(idx),"]/[idx=",idx[length(idx)],"]/[ddax=",(zzix),"]=END") }
          multidata.info[[ p(ss) ]][["data"]] <- ddax;
        }
      }
    }



  }

  sampCounts[[ss]][[p("mrgCT.ALL")]] <- 0;
  for(ee in eventlist ){
    sampCounts[[ss]][[p("mrgCT.ALL")]] <- sampCounts[[ss]][[p("mrgCT.ALL")]] + sampCounts[[ss]][[p("mrgCT.",ee)]]
  }
}

#sampCounts.backup <- sampCounts
#multdata.info.backup <- multidata.info#
#
#table(
#  OLD = sapply( sampCounts.backup,'[[',"mrgCT.ALL"),
#  NEW = sapply( sampCounts,'[[',"mrgCT.ALL")
#)

################################################################################


###############################################
ss = "REBC-ADM2-TP-NB"

#DEBUGGING
# REBC-ADM3-TP-NB
# REBC-ADJH-TP-NT
# 
for( ss in samplist ){
   multidata.info[[ p(ss) ]]
   multidata.info[[ p(ss) ]][["events"]]
   sampCounts[[ss]]

   if( sampCounts[[ss]][["mrgCT.ALL"]] != length( multidata.info[[ p(ss) ]][["events"]] ) ){
     message("ERROR!!!");
   }
   message(ss," : mrgCT = ",sampCounts[[ss]][["mrgCT.ALL"]]," / ",length( multidata.info[[ p(ss) ]][["events"]] ) )
}
###################################


################################################################################




sampCountData <- do.call(rbind.data.frame,sampCounts);
setDT(sampCountData, keep.rownames = TRUE)[]
names(sampCountData)[[1]] <- "sample";
sampCountData.LABELED_REBC_ID <- data.frame( LABELED_REBC_ID = sapply(strsplit(sampCountData$sample,"-"),function(ss){
  p(ss[[1]],"-",ss[[2]]);
}))
sampCountData <- cbind.data.frame(sampCountData.LABELED_REBC_ID,sampCountData);





zzlist <- grepv("mrgCT.",names(sampCounts[[1]]))
zzlist <- zzlist[zzlist != "mrgCT.ALL"]
for( ss in names(multidata.info) ){
   ddax <- multidata.info[[ p(ss) ]][["data"]]
   ddad <- data.frame(SCNA.type = character(), CT = numeric());
   ddsc <- sampCounts[[ss]];
   #message("ss = ",ss);
   if(! is.null(ddax)){
      for(zz in zzlist){
        if( ddsc[[zz]] > 0 ){
          ddad <- rbind(ddad,data.frame(SCNA.type = zz, CT = ddsc[[zz]]));
        }
      }
   }
   multidata.info[[ p(ss) ]][["result"]] <- ddad
}



for( ss in samplist ){

   if(is.null(multidata.info[[ p(ss) ]])){
     ddad <- data.frame(SCNA.type = character(), CT = numeric());
     ddsc <- sampCounts[[ss]];
     multidata.info[[ss]] <- list();
     multidata.info[[ss]][["data"]] <- 
     for(zz in zzlist){
        if( ddsc[[zz]] > 0 ){
          ddad <- rbind(ddad,data.frame(SCNA.type = zz, CT = ddsc[[zz]]));
        }
     }
     kkn <-  dda$sample == ss & (dda$IS_SCNA == 1 | dda$event.type == "NA");
     ddaxn <- dda[kkn,]
     multidata.info[[ p(ss) ]][["data"]] <- ddaxn;
     multidata.info[[ p(ss) ]][["result"]] <- ddad
   }

}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


multidata.info.events.len <- sapply(sapply(multidata.info,'[[',"events"),length);
multidata.info.eventsINFO.len <- sapply(sapply(multidata.info,'[[',"events.info"),length);
multidata.info.eventsIDX.len <- sapply(sapply(multidata.info,'[[',"events.idx"),length);
multidata.info.eventsDATA.len <- sapply(sapply(multidata.info,'[[',"events.data"),length);


table(multidata.info.events.len,multidata.info.eventsINFO.len)
table(multidata.info.events.len,multidata.info.eventsINFO.len)
table(multidata.info.events.len,multidata.info.eventsIDX.len)
table(multidata.info.events.len,multidata.info.eventsDATA.len)

which(multidata.info.events.len !=  multidata.info.eventsINFO.len)
which(multidata.info.events.len !=  multidata.info.eventsIDX.len)
which(multidata.info.events.len !=  multidata.info.eventsDATA.len)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


keepcols.xx <- c(names(dda)[1:4],"Nstat","event.type","event.clonality","event.info","rownum","notes");
options(width=1000)

#----------- NON EVENT: REBC-ADM3-TP-NB / 11 / gain
#----------- NON EVENT: REBC-ADJH-TP-NT / 6 / gain
#----------- NON EVENT: REBC-ADJH-TP-NT / 6 / gain

samps.with.err <- c("REBC-ADM3-TP-NB","REBC-ADJH-TP-NT")

if(SINK) {
sink(p("res/temp.TEST.log"))
for(ss in samps.with.err ){
  cat("===========================================================================================\n")
  cat("===========================================================================================\n")
  print(ss);
  cat("\n");
  print( dda[dda$sample == ss,keepcols.xx] , width=1000000);
  cat("RESULT:\n");
  xkx <- t(sampCountData[sampCountData$sample == ss,names(sampCountData)[startsWith(names(sampCountData),"mrg")] ] )
  colnames(xkx) <- ""
  print(xkx)

  cat("FINAL RESULT:\n");
  xkx <- t(sampCountData.final[sampCountData$sample == ss,names(sampCountData)[startsWith(names(sampCountData),"mrg")] ] )
  colnames(xkx) <- ""
  print(xkx)

}
sink()
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

keepcols.xx <- c(names(dda)[1:4],"Nstat","event.type","event.clonality","event.info","rownum","notes");
options(width=1000)



names(sampCountData.B)


zzlist <- grepv("mrgCT.",names(sampCounts[[1]]))
for(zz in zzlist){
  message("\n",zz);
  print( table(A = sampCountData[[zz]],B=sampCountData.B[[zz]]) )
}
zzlist <- grepv("mrgCT.",names(sampCounts[[1]]))
for(zz in zzlist){
  message("\n",zz);
  print( table(A = sampCountData[[zz]] == sampCountData.B[[zz]]) )
}

flag.samp.issues <- list();
flag.samp.list <- c();
flag.samp.rebc <- c();

for(zz in zzlist){
  message("\n",zz);
  xxz <- sampCountData$sample[ (A = sampCountData[[zz]] != sampCountData.B[[zz]]) ]
  flag.samp.list <- unique(c(flag.samp.list,xxz))
  xxr <- sampCountData$LABELED_REBC_ID[ (A = sampCountData[[zz]] != sampCountData.B[[zz]]) ]
  flag.samp.rebc <- unique(c(flag.samp.rebc ,xxr))
  for(ss in xxz){
    flag.samp.issues[[ss]] <- c(flag.samp.issues[[ss]],zz)
  }
}


if(SINK) {

sink(p("temp.flagFile.untypedSpanEvents.",RUNVER,".log"))
for(ss in flag.samp.list ){
  cat("===========================================================================================\n")
  cat("===========================================================================================\n")
  print(ss);
  print(p("issues:",p(flag.samp.issues[[ss]],collapse=",")))
  print( multidata.info.B[[ss]][["data"]][,keepcols.xx] , width=1000000);
  cat("Result: BREAK AT UNTYPED SPANS:\n");
  print(multidata.info[[ss]][["result"]]  , width=1000000)
  cat("Result: MERGE ACROSS UNTYPED SPANS:\n");
  print(multidata.info.B[[ss]][["result"]]  , width=1000000)
  #cat("manual:\n");
  #print(t(ddc3[ f.na(ddc3$sample == ss), c("sample",grepv("_manual",names(ddc3))) ]))
}
sink()
}

DK.notes.spanMerge <- list(
   list(id="AF8T", merge=F,
        notes="IGNORE"),
   list(id="ADKY", merge=F,
        notes="Ignore, not in 354"),
   list(id="AF70", merge=F,
        notes="Ignore, not in 354"),
   list(id="ADM3", merge=F,
        notes="Complex Case 1"),
   list(id="AF65", merge=F,
        notes="Complex Case 2"),
   list(id="ACB3",merge=F,
        notes="NO MERGE"),
   list(id="ACBU",merge=F,
        notes="NO MERGE"),
   list(id="ACDX",merge=F,
        notes="NO MERGE"),
   list(id="ADJ9",merge=F,
        notes="NO MERGE"),
   list(id="ADJJ",merge=F,
        notes="NO MERGE"),
   list(id="ADLM",merge=F,
        notes="NO MERGE"),
   list(id="ADM0",merge=F,
        notes="NO MERGE"),
   list(id="ADMH",merge=F,
        notes="NO MERGE"),
   list(id="AF6I",merge=F,
        notes="NO MERGE"),
   list(id="ACCG",merge=T,
        notes="MERGE ME"),
   list(id="ADL2",merge=T,
        notes="MERGE ME"),
   list(id="ADJH",merge=T,
        notes="MERGE ME"),
   list(id="ADM3",merge=T,
        notes="MERGE ME")
)

DK.notes.spanMerge <- do.call(rbind.data.frame,DK.notes.spanMerge)
DK.notes.spanMerge$LABELED_REBC_ID <- p("REBC-",DK.notes.spanMerge$id);

DK.notes.spanMerge$foundInFlagList <- sapply( DK.notes.spanMerge$LABELED_REBC_ID , function(ss){
  ss %in% flag.samp.rebc
})

data.frame(
  LABELED_REBC_ID = flag.samp.rebc,
  foundInNoteList = sapply(flag.samp.rebc,function(ss){
    ss %in% DK.notes.spanMerge$LABELED_REBC_ID
  })
)

#7706 REBC-ADM3-TP-NB         11  29375409  29531275


swapToMerge <- DK.notes.spanMerge$LABELED_REBC_ID[DK.notes.spanMerge$merge]

sampEventData.final  <- lapply(multidata.info,'[[',"events");
sampEventData.clonal <- lapply(multidata.info,'[[',"events.clonal");
sampEventInfo        <- lapply(multidata.info,'[[',"events.info");
sampEventInfo.idx         <- lapply(multidata.info,'[[',"events.idx");
sampEventInfo.data        <- lapply(multidata.info,'[[',"events.data");

sampCountData.final <- sampCountData;

for(ssr in swapToMerge){
  ss <- sampCountData$sample[ sampCountData$LABELED_REBC_ID == ssr ];
  ix <- which(sampCountData$sample == ss)
  for(zz in names(sampCountData)){
    sampCountData.final[[zz]][[ix]] <- sampCountData.B[[zz]][[ix]]
  }
  sampEventData.final[[ss]]  <- multidata.info.B[[ss]][["events"]]
  sampEventData.clonal[[ss]] <- multidata.info.B[[ss]][["events.clonal"]]
  sampEventInfo[[ss]] <- multidata.info.B[[ss]][["events.info"]]
  sampEventInfo.idx[[ss]] <- multidata.info.B[[ss]][["events.idx"]]
  sampEventInfo.data[[ss]] <- multidata.info.B[[ss]][["events.data"]]
}

rownames(sampCountData.final) <- gsub("REBC-","",sampCountData.final$LABELED_REBC_ID)


###############################################################################################
###############################################################################################
###############################################################################################

ss = "REBC-ADMA-TP-NB"

sampEventInfo[[ss]]
sampEventInfo.data[[ss]]

for( ss in names(sampEventInfo)){
   ssi <- sampEventInfo[[ss]]
   for( i in seq_along(ssi)){
     ssi[[i]][["NXCCF"]] <- max( sampEventInfo.data[[ss]][[i]][["NXCCF"]], na.rm=T)
   }
   sampEventInfo[[ss]] <- ssi
}

###############################################################################################
###############################################################################################
###############################################################################################

sampCountData.final$fullChrom.base <- sapply(1:nrow(sampCountData.final),function(i){
   ss <- sampCountData.final$sample[[i]];
   ddfc <- multidata.info[[ss]][["fullchrom"]]
   if( length(ddfc) > 0){
     ddfc.isflag  <- unlist(sapply(ddfc,'[[',"flag"));
     ddfc.flagged <- ddfc[  ddfc.isflag]
     ddfc.noflag  <- ddfc[! ddfc.isflag]
     if(length(ddfc.noflag) == 0){
       ddfc.string <- "";
     } else {
       ddfc.string <- p(sort(sapply(ddfc.noflag,function(xx){
          p(xx[["chrom"]],":",xx[["ee"]])
       })),collapse=",")
     }
     ddfc.flagstring <- p((sapply(ddfc.flagged,function(xx){
        p(xx[["chrom"]],":",xx[["ee"]])
     })),collapse=",")
     ddfc.string;
   } else {
     ""
   }
})
sampCountData.final$fullChrom.flagged <- sapply(1:nrow(sampCountData.final),function(i){
   ss <- sampCountData.final$sample[[i]];
   ddfc <- multidata.info[[ss]][["fullchrom"]]
   if( length(ddfc) > 0){
     ddfc.isflag  <- unlist(sapply(ddfc,'[[',"flag"));
     ddfc.flagged <- ddfc[  ddfc.isflag]
     ddfc.noflag  <- ddfc[! ddfc.isflag]
     ddfc.string <- p((sapply(ddfc.noflag,function(xx){
        p(xx[["chrom"]],":",xx[["ee"]])
     })),collapse=",")
     ddfc.flagstring <- p((sapply(ddfc.flagged,function(xx){
        p(xx[["chrom"]],":",xx[["ee"]])
     })),collapse=",")
     ddfc.flagstring ;
   } else {
     ""
   }
})


manual.add.fullchrom <- c(
  c("ACF1:chr4"),
  c("ADK6:chr22"),
  c("ADL2:chr22"),
  c("ADLI:chr22")
)

multidata.info.fullchrom <- sapply(multidata.info,'[[',"fullchrom")
multidata.info.fullchrom <- sapply(names(multidata.info.fullchrom),function(ss){
  sampIdShort <- rownames(sampCountData.final)[sampCountData.final$sample == ss]
  zz = lapply(  multidata.info.fullchrom[[ss]], function(cdata){
    pairString <- p(sampIdShort,":chr",cdata$chrom)
    cdata$manualKeep <- pairString %in% manual.add.fullchrom;
    cdata$keep <- (! cdata$flag) | cdata$manualKeep
    cdata
  })
  zz[ order( sapply(zz,'[[',"chrom" ) ) ]
})
multidata.info.fullchrom.filt <- sapply(names(multidata.info.fullchrom),function(ss){
  sampIdShort <- rownames(sampCountData.final)[sampCountData.final$sample == ss]
  xx = multidata.info.fullchrom[[ ss ]]
  kp <- sapply(xx,'[[',"keep");
  if(length(kp) > 0){
    xx[kp]
  } else {
    list();
  }
})
multidata.info.fullchrom.dropped <- sapply(names(multidata.info.fullchrom),function(ss){
  sampIdShort <- rownames(sampCountData.final)[sampCountData.final$sample == ss]
  xx = multidata.info.fullchrom[[ ss ]]
  kp <- ! sapply(xx,'[[',"keep");
  if(length(kp) > 0){
    xx[kp]
  } else {
    list();
  }
})

##########
multidata.info.fullchrom.OLD <- sapply(multidata.info,'[[',"fullchrom.FULL")
multidata.info.fullchrom.OLD <- sapply(names(multidata.info.fullchrom.OLD),function(ss){
  sampIdShort <- rownames(sampCountData.final)[sampCountData.final$sample == ss]
  zz = lapply(  multidata.info.fullchrom.OLD[[ss]], function(cdata){
    pairString <- p(sampIdShort,":chr",cdata$chrom)
    cdata$manualKeep <- pairString %in% manual.add.fullchrom;
    cdata$keep <- (! cdata$flag) | cdata$manualKeep
    cdata
  })
  zz[ order( sapply(zz,'[[',"chrom" ) ) ]
})
multidata.info.fullchrom.OLD.filt <- sapply(names(multidata.info.fullchrom.OLD),function(ss){
  sampIdShort <- rownames(sampCountData.final)[sampCountData.final$sample == ss]
  xx = multidata.info.fullchrom.OLD[[ ss ]]
  kp <- sapply(xx,'[[',"keep");
  if(length(kp) > 0){
    xx[kp]
  } else {
    list();
  }
})

sampIdShort <- sapply(strsplit(manual.add.fullchrom,":"),'[[',1)


sampCountData.final$fullChrom.info <- sapply(1:nrow(sampCountData.final),function(i){
   ss <- sampCountData.final$sample[[i]];
   nn <- rownames(sampCountData.final)[[i]];
   ddfc <- multidata.info[[ss]][["fullchrom"]]
   if( length(ddfc) > 0){
     ddfc.isflag  <- unlist(lapply(ddfc,function(zz){
        ssee <- p(nn,":chr",zz$chrom);
        zz$flag && (! ssee %in% manual.add.fullchrom)
     }))
     ddfc.flagged <- ddfc[  ddfc.isflag]
     ddfc.noflag  <- ddfc[! ddfc.isflag]
     if(length(ddfc.noflag) == 0){
       ddfc.string <- "";
     } else {
       ddfc.string <- p(sort(sapply(ddfc.noflag,function(xx){
          p(xx[["chrom"]],":",xx[["ee"]])
       })),collapse=",")
     }
     ddfc.flagstring <- p(sapply(ddfc.flagged,function(xx){
        p(xx[["chrom"]],":",xx[["ee"]])
     }),collapse=",")
     ddfc.string;
   } else {
     ""
   }
})

as.character(sampCountData.final$fullChrom.info)

for(ee in eventlist ){
  sampCountData.final[[p("fullChrom.",ee)]] <- sapply(strsplit(sampCountData.final$fullChrom.info,","),function(xz){
    sum(grepl(p(ee),xz,fixed=T))
  })
  sampCountData.final[[p("fullChrom.clonal.",ee)]] <- sapply(strsplit(sampCountData.final$fullChrom.info,","),function(xz){
    sum(grepl(p("clonal.",ee),xz,fixed=T))
  })
}


sampCountData.final[sampCountData.final$fullChrom.info != sampCountData.final$fullChrom.base,]

sampCountData.final[355,]

multidata.info[["REBC-AF6W-TP-NT"]]

sampCountData.final[c("ADJH","ADM3","AF65"),]


multidata.info.Bx <- multidata.info.B
names(multidata.info.Bx) <- gsub("-[A-Z][A-Z]$","",gsub("-[A-Z][A-Z]-[A-Z][A-Z]$","",gsub("REBC-","",names(multidata.info.B))))

multidata.info.Bx[["AF65"]][["data"]][,keepcols.xx]

#write.table(multidata.info.Bx[["AF65"]][["data"]][,keepcols.xx],file="test.txt",row.names=F,col.names=T,sep='\t',quote=F)

####################################
########################################################################
########################################################################
########################################################################
########################################################################
####################################
multidata.info[[ss]][["fullchrom"]]

sampCountData.prefilter <- sampCountData.final

#sampEventInfo      #  <- lapply(multidata.info,'[[',"events.info");
#sampEventInfo.idx  #       <- lapply(multidata.info,'[[',"events.idx");
#sampEventInfo.data #       <- lapply(multidata.info,'[[',"events.data");
sampEventInfo.raw <- sampEventInfo 

sampEventInfo <- sapply(names(sampEventInfo),function(ss){
   fc      <- multidata.info.fullchrom.filt[[ss]];
   fc.old  <- multidata.info.fullchrom.OLD.filt[[ss]];
   fc.drop <- multidata.info.fullchrom.dropped[[ss]];

   se  <- sampEventInfo[[ss]]
   fcc <- sapply(fc,'[[',"chrom")
   fce <- sapply(fc,'[[',"event")
   fcc.old <- sapply(fc.old,'[[',"chrom")
   fce.old <- sapply(fc.old,'[[',"event")
   fcc.drop <- sapply(fc.drop ,'[[',"chrom")
   fce.drop <- sapply(fc.drop ,'[[',"event")
   lapply(se,function(see){
     see$fullChrom <- any( see$chrom == fcc & see$event == fce )
     see$fullChrom.FULL <- any( see$chrom == fcc.old & see$event == fce.old )
     see$fullChrom.DROP <- any( see$chrom == fcc.drop & see$event == fce.drop )
     see
   })
})
ss = "REBC-AF8W-TP-NB"

#multidata.info.fullchrom.filt[["REBC-AF92-TP-NB"]]


sampEventInfo <- sapply(names(sampEventInfo),function(ss){
   se <- sampEventInfo[[ss]];
   lapply(se,function(sx){
     centspan <- centromere.data[centromere.data$chr == sx$chrom,c("chromStart","chromEnd"),drop=T]
     sx$intersects.centromere <- sx$start <= centspan[[2]] && centspan[[1]] <= sx$end
     sx$inside.centromere     <- centspan[[1]] <= sx$start && sx$end <= centspan[[2]]
     sx$intersects.parm <- sx$start <= centspan[[2]]
     sx$inside.parm            <- sx$end   <= centspan[[2]]
     sx$mostly.inside.parm     <- sx$start <= centspan[[1]] && sx$end <= centspan[[2]] + 15000000
     sx$chrom.is.acrocentric  <- sx$chrom %in% acrocentric.chroms
     sx$filter.drop <- (! sx$fullChrom ) && (
                          sx$inside.centromere || 
                          ( sx$mostly.inside.parm && sx$chrom.is.acrocentric )
                       )        
     sx       
   })
})
sampEventInfo.dataset <- sampEventInfo.data
sampEventInfo.data <- lapply(names(sampEventInfo),function(ss){
   se <- sampEventInfo[[ss]];
   lapply(se,function(sx){
     c(list(sample = ss), sx)
   })
})
sampEventInfo.data <- unlist(sampEventInfo.data,recursive=F)
sampEventInfo.data <- do.call(rbind.data.frame,sampEventInfo.data)

if(SINK){
   write.table(sampEventInfo.data,file=p("eventTable.SCNA.SWH.",RUNVER,"m001.txt"),sep='\t',quote=F,row.names=F,col.names=T);
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


sampEventInfo.full <- sampEventInfo;
sampEventData.prefilter <- sampEventData.final
sampEventData.clonal.prefilter <- sampEventData.clonal;

sampEventInfo <- sapply(names(sampEventInfo.full),function(ss){
   se <- sampEventInfo[[ss]];
   se[ ! sapply(se,'[[',"filter.drop") ]
})

sampEventInfo.filtered <- sapply(names(sampEventInfo.full),function(ss){
   se <- sampEventInfo.full[[ss]];
   kp <- sapply(se,'[[',"filter.drop")
   if(length(kp) > 0){
     se[ kp ]
   } else {
     list()
   }
})

sampEventData.final <- sapply(names(sampEventInfo),function(ss){
   se <- sampEventInfo[[ss]];
   if(length(se) > 0){
     sapply(se,function(sx){
       p(sx$event,":chr",sx$chrom)
     })
   } else {
     character()
   }
})
sampEventData.clonal <- sapply(names(sampEventInfo),function(ss){
   se <- sampEventInfo[[ss]];
   kp <- sapply(se,'[[',"clonal")
   if( length(kp) > 0 && sum(kp) > 0){
     sort(sapply(se[kp],function(sx){
       p(sx$event,":chr",sx$chrom)
     }))
   } else {
     character();
   }
})
sampEventData.fullchrom <- sapply(names(sampEventInfo),function(ss){
   se <- sampEventInfo[[ss]];
   kp <- sapply(se,'[[',"fullChrom")
   if( length(kp) > 0 && sum(kp) > 0){
     sort(sapply(se[kp],function(sx){
       clonal.string <- ifelse(sx$clonal,"clonal.","");
       p(sx$chrom,":",clonal.string,sx$event)
     }))
   } else {
     character();
   }
})
sampEventData.fullchrom.FULL <- sapply(names(sampEventInfo),function(ss){
   se <- sampEventInfo[[ss]];
   kp <- sapply(se,'[[',"fullChrom.FULL")
   if( length(kp) > 0 && sum(kp) > 0){
     sort(sapply(se[kp],function(sx){
       clonal.string <- ifelse(sx$clonal,"clonal.","");
       p(sx$chrom,":",clonal.string,sx$event)
     }))
   } else {
     character();
   }
})
sampEventData.fullchrom.DROP <- sapply(names(sampEventInfo),function(ss){
   se <- sampEventInfo[[ss]];
   kp <- sapply(se,'[[',"fullChrom.DROP")
   if( length(kp) > 0 && sum(kp) > 0){
     sort(sapply(se[kp],function(sx){
       clonal.string <- ifelse(sx$clonal,"clonal.","");
       p(sx$chrom,":",clonal.string,sx$event)
     }))
   } else {
     character();
   }
})

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

scd.chromInfo <- (sampEventData.fullchrom);
chrom.list <- p("",c(1:22));

sampTallies.final <- data.frame(sample = sampCountData.final$sample );

ii=358
cc ="18"
ee = "CopyNeutralLOH"

for(cc in chromlist){
  for(ee in eventlist ){
    sampTallies.final[[p("fullChrom.",cc,".",ee)]] <- sapply(1:nrow(sampCountData.final),function(ii){
      p(cc,":",ee) %in% scd.chromInfo[[ii]] || p(cc,":clonal.",ee) %in% scd.chromInfo[[ii]]
    })
  }
}
for(cc in chromlist){
  for(ee in eventlist ){
    sampTallies.final[[p("events.",cc,".",ee)]] <- sapply( sampCountData.final[["sample"]],function(ss){
       sum( (p(ee,":chr",cc) == sampEventData.final[[ ss ]]  ) )
    })
  }
}
for(cc in chromlist){
  for(ee in eventlist ){
    sampTallies.final[[p("nonChrom.",cc,".",ee)]] <- sampTallies.final[[p("events.",cc,".",ee)]] - sampTallies.final[[p("fullChrom.",cc,".",ee)]]
  }
}
for(cc in chromlist){
  sampTallies.final[[p("fullChrom.",cc,".","ALL")]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(eventlist,function(ee){
        sampTallies.final[[p("fullChrom.",cc,".",ee)]][[ii]]
      }))
    })
  sampTallies.final[[p("events.",cc,".","ALL")]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(eventlist,function(ee){
        sampTallies.final[[p("events.",cc,".",ee)]][[ii]]
      }))
    })
  sampTallies.final[[p("nonChrom.",cc,".","ALL")]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(eventlist,function(ee){
        sampTallies.final[[p("nonChrom.",cc,".",ee)]][[ii]]
      }))
    })
}
for(ee in c(eventlist,"ALL") ){
  sampTallies.final[[p("fullChrom.","ALL",".",ee)]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(chromlist,function(cc){
        sampTallies.final[[p("fullChrom.",cc,".",ee)]][[ii]]
      }))
    })
  sampTallies.final[[p("events.","ALL",".",ee)]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(chromlist,function(cc){
        sampTallies.final[[p("events.",cc,".",ee)]][[ii]]
      }))
    })
  sampTallies.final[[p("nonChrom.","ALL",".",ee)]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(chromlist,function(cc){
        sampTallies.final[[p("nonChrom.",cc,".",ee)]][[ii]]
      }))
    })
}

#####################
#CLONAL
for(cc in chromlist){
  for(ee in eventlist ){
    sampTallies.final[[p("fullChrom.clonal.",cc,".",ee)]] <- sapply(1:nrow(sampCountData.final),function(ii){
      p(cc,":clonal.",ee) %in% scd.chromInfo[[ii]]
    })
  }
}
for(cc in chromlist){
  for(ee in eventlist ){
    sampTallies.final[[p("events.clonal.",cc,".",ee)]] <- sapply( sampCountData.final[["sample"]],function(ss){
       sum( (p(ee,":chr",cc) == sampEventData.clonal[[ ss ]]  ) )
    })
  }
}
for(cc in chromlist){
  for(ee in eventlist ){
    sampTallies.final[[p("nonChrom.clonal.",cc,".",ee)]] <- sampTallies.final[[p("events.clonal.",cc,".",ee)]] - sampTallies.final[[p("fullChrom.clonal.",cc,".",ee)]]
  }
}
for(cc in chromlist){
  sampTallies.final[[p("fullChrom.clonal.",cc,".","ALL")]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(eventlist,function(ee){
        sampTallies.final[[p("fullChrom.clonal.",cc,".",ee)]][[ii]]
      }))
    })
  sampTallies.final[[p("events.clonal.",cc,".","ALL")]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(eventlist,function(ee){
        sampTallies.final[[p("events.clonal.",cc,".",ee)]][[ii]]
      }))
    })
  sampTallies.final[[p("nonChrom.clonal.",cc,".","ALL")]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(eventlist,function(ee){
        sampTallies.final[[p("nonChrom.clonal.",cc,".",ee)]][[ii]]
      }))
    })
}
for(ee in c(eventlist,"ALL") ){
  sampTallies.final[[p("fullChrom.clonal.","ALL",".",ee)]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(chromlist,function(cc){
        sampTallies.final[[p("fullChrom.clonal.",cc,".",ee)]][[ii]]
      }))
    })
  sampTallies.final[[p("events.clonal.","ALL",".",ee)]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(chromlist,function(cc){
        sampTallies.final[[p("events.clonal.",cc,".",ee)]][[ii]]
      }))
    })
  sampTallies.final[[p("nonChrom.clonal.","ALL",".",ee)]] <- sapply(1:nrow(sampCountData.final),function(ii){
      sum(sapply(chromlist,function(cc){
        sampTallies.final[[p("nonChrom.clonal.",cc,".",ee)]][[ii]]
      }))
    })
}
#########################
#NONCLONAL

for(cc in c(chromlist,"ALL")){
  for(ee in c(eventlist,"ALL") ){
    for(ty in c("fullChrom","events","nonChrom")){
      sampTallies.final[[p(ty,".nonclonal.",cc,".",ee)]] <- sampTallies.final[[p(ty,".",cc,".",ee)]] - sampTallies.final[[p(ty,".clonal.",cc,".",ee)]]
    }
  }
}
#########################

sampCountData.new <- sampCountData[,c("LABELED_REBC_ID","sample",grepv("^rawCT",names(sampCountData.final)),grepv("^chrCT",names(sampCountData.final)) )]

sampCountData.new$events.list <- sapply(sampEventData.final,function(ev){
  p(ev,collapse=",");
})
sampCountData.new$events.clonal.list <- sapply(sampEventData.clonal,function(ev){
  p(ev,collapse=",");
})
sampCountData.new$fullchrom.events <- sapply(sampEventData.fullchrom,function(ev){
  p(gsub(".clonal","",ev),collapse=",");
})
sampCountData.new$fullchrom.events.clonal <- sapply(sampEventData.fullchrom,function(ev){
  p(grepv("clonal",ev),collapse=",");
})
sampCountData.new$fullchrom.flaggedAndManualDrop <- sapply(sampEventData.fullchrom.DROP,function(ev){
  p(ev,collapse=",");
})


options(width=250)

head(sampCountData.new)



grepv("mrgCT",names(sampCountData.final))

for(ee in eventlist){
  sampCountData.new[[p("mrgCT.",ee)]] <- sapply( sampEventData.final, function(se){
    sum(grepl(p("^",ee),se))
  })
  sampCountData.new[[p("mrgCT.clonal.",ee)]] <- sapply( sampEventData.clonal, function(se){
    sum(grepl(p("^",ee),se))
  })
}
sampCountData.new[[p("mrgCT.","ALL")]] <- sapply( sampEventData.final, function(se){
    length(se);
})
sampCountData.new[[p("mrgCT.clonal.","ALL")]] <- sapply( sampEventData.clonal, function(se){
    length(se);
})

#sampEventData.fullchrom

for(ee in eventlist){
  sampCountData.new[[p("fullChromCT.",ee)]] <- sapply( sampEventData.fullchrom, function(se){
    sum(grepl(p(ee),se))
  })
  sampCountData.new[[p("fullChromCT.clonal.",ee)]] <- sapply( sampEventData.fullchrom, function(se){
    sum(grepl(p("clonal.",ee),se))
  })
}
sampCountData.new[[p("fullChromCT.","ALL")]] <- sapply( sampEventData.fullchrom, function(se){
    length(se);
})
sampCountData.new[[p("fullChromCT.clonal.","ALL")]] <- sapply( sampEventData.fullchrom, function(se){
    sum(grepl(p("clonal"),se))
})


if(WRITETABLES){
  write.table(sampEventInfo.data,file=p("out/eventTable.SCNA.SWH.",RUNVER,"m001.txt"),sep='\t',quote=F,row.names=F,col.names=T);
  write.table(sampCountData.new ,file=p("out/simpleCountTable.SCNA.SWH.",RUNVER,"m10.txt"),row.names=F,col.names=T,sep='\t',quote=F)
  write.table(sampCountData.full ,file=p("out/superTable.SCNA.SWH.",RUNVER,"m10.txt"),row.names=F,col.names=T,sep='\t',quote=F)
  write.table(dda ,file=p("out/debug.dda.",RUNVER,"m10.txt"),row.names=F,col.names=T,sep='\t',quote=F)

}

###########################################################################
###########################################################################
###########################################################################
###########################################################################








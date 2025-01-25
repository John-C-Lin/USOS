# Determine subset of receptors from aircraft data
# Jan. 25th, 2025 by John C. Lin (John.Lin@utah.edu)

#################
fltdatdir <- "USOS_Aircraft_Data/"
plotTF <- FALSE

#criteria to create new receptor
DD.threshold <- 0.1  #[deg]
DZ.threshold <- 250  #[m]
DT.threshold <- 180  #[s]

finalresultname <- "USOS_selected_receptors.rds"
#################

receptors.all <- NULL
fltobjnames <- list.files(fltdatdir,".rds")
for(ff in 1:length(fltobjnames)){

  fltobjname <- fltobjnames[ff]
  fltlabel <- substring(fltobjname,1,nchar(fltobjname)-4)
  #fltobjname <- "20240726_RA_L1.rds"  #output from "readin_USOS_fltdat_ICARTT.r"
  print(paste("----------- Reading in:",fltobjname,"---------------"))
  fdat <- readRDS(paste0(fltdatdir,"/",fltobjname))


#remove NAs
sel<-!is.na(fdat[,"Latitude"])&!is.na(fdat[,"Longitude"])&!is.na(fdat[,"ALTAGL"])
fdat<-fdat[sel,]

YYYY <- substring(fltobjname,1,4)
start_of_year <- as.POSIXct(paste0(YYYY, "-01-01"), tz = "UTC")
# Add the day_of_year as days (subtract 1 because day 1 is Jan 1)
Time_UTC <- start_of_year + (floor(fdat$Day_of_year) - 1) * 86400 + fdat$StartTime_UTsec # 86400 seconds in a day
fdat <- data.frame(Time_UTC,fdat)

#temporally average lat/lon/AGL at MINUTE timescales (original data are in SECOND timescales)
#Time.min <- format(Time_UTC,format="%Y-%m-%d %H:%M")
#f <- function(x,ind){return(tapply(x,ind,mean,na.rm=T))}
#dat <- apply(fdat[,c("Latitude","Longitude","ALTAGL")],2,f,Time.min)
#Time.UTC <- strptime(rownames(dat),format="%Y-%m-%d %H:%M", tz="GMT")
#dat <- data.frame(Time.UTC,dat)

dat <- fdat[,c("Time_UTC","StartTime_UTsec","Day_of_year","Latitude","Longitude","ALTAGL")]
YYYYMMDDs <- unique(format(dat[,"Time_UTC"],"%Y%m%d"))


#for(i in 1:length(YYYYMMDDs)){
#for(i in 1:3){
#  sel<-YYYYMMDD%in%YYYYMMDDs[i]
#  print(paste("----- Processing:",YYYYMMDDs[i],"-----"))
#  if(sum(sel)<15){print("flight too short; skipping");runreceptorTF<-rep(FALSE,sum(sel));runreceptorTF.all<-c(runreceptorTF.all,runreceptorTF);next}
#  dat<-DAT[sel,]
  #Criteria for new receptor:
  #  a) whenever lat/lon change by DD.threshold [deg]
  #  b) AGL change by DZ.threshold [m]
  #  c) time change by DT.threshold [s]
  dx <- c(0,diff(dat[,"Longitude"]));dy <- c(0,diff(dat[,"Latitude"]))
  dx <- dx*cos(dat[,"Latitude"]*pi/180)   #apply cos(lat) factor for longitude
  dd <- sqrt(dx^2+dy^2)   #total distance [deg]
  dz <- c(0,diff(dat[,"ALTAGL"])) #[m]
  dt <- c(0,diff(as.numeric(dat$Time_UTC))) #[m]

  dd.c <- abs(cumsum(dd))
  xbreaks <- seq(dd.c[1],dd.c[length(dd.c)],DD.threshold)
  xcut <- cut(dd.c,breaks=xbreaks)
  nxcut <- as.numeric(xcut)
  ddTF <- !duplicated(nxcut)

  dz.c <- cumsum(abs(dz))
  xbreaks <- seq(dz.c[1],dz.c[length(dz.c)],DZ.threshold)
  xcut <- cut(dz.c,breaks=xbreaks)
  nxcut <- as.numeric(xcut)
  dzTF <- !duplicated(nxcut)

  dt.c <- cumsum(abs(dt))
  xbreaks <- seq(dt.c[1],dt.c[length(dt.c)],DT.threshold)
  xcut <- cut(dt.c,breaks=xbreaks)
  nxcut <- as.numeric(xcut)
  dtTF <- !duplicated(nxcut)

  # receptor "flag"--whether to run receptor or not
  runreceptorTF <- ddTF|dzTF|dtTF
  print(paste("number of receptors to run:",sum(runreceptorTF)))
if(plotTF){
  #X11(width=10,height=6);par(mfrow=c(1,2),cex.lab=1.3,cex.axis=1.3)
  plot(dat[,c("Longitude","Latitude")],pch=16,cex=1.2,main=fltobjname)
  points(dat[ddTF,c("Longitude","Latitude")],pch=16,cex=1.0,col="orange")
  points(dat[dzTF,c("Longitude","Latitude")],pch=17,cex=1.0,col="lightgreen")
  points(dat[dtTF,c("Longitude","Latitude")],pch=16,cex=1.0,col="blue")

  plot(dat[,c("Time_UTC","ALTAGL")],pch=16,cex=1.5,type="o",main=fltlabel)
  points(dat[dzTF,c("Time_UTC","ALTAGL")],pch=17,cex=1.3,col="lightgreen")
  points(dat[dtTF,c("Time_UTC","ALTAGL")],pch=17,cex=1.3,col="blue")
} #if(plotTF){

  print(paste("% of receptors selected:",signif(100*sum(runreceptorTF)/length(runreceptorTF),4),"%"))
  receptors <- dat[runreceptorTF,][-1,]
  receptors <- data.frame(fltlabel,receptors)

  # output 1 merged receptor list from all flight days
  receptors.all <- rbind(receptors.all,receptors)
  gc()
} # for(ff in 1:length(fltobjnames)){

saveRDS(receptors.all,file=finalresultname)
print(paste(finalresultname,"generated"))

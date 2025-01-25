# Reads in USOS aircraft data, which is in ICARTT format
# Jan. 25th, 2025 by John C. Lin (John.Lin@utah.edu)

#######################
datdir <- "USOS_Aircraft_Data/"
plotTF <- FALSE
#######################

xfiles <- list.files(path=datdir,pattern="USOS-ARL-Suite")
for(ff in 1:length(xfiles)){
xfile <- xfiles[ff] # e.g., "USOS-ARL-Suite_TwinOtter_20240719_RA_L1.ict"
xregexpr <- regexpr("TwinOtter",xfile)[1]
YYYYMMDD <- substring(xfile,xregexpr+10,xregexpr+17)
fltlabel <- substring(xfile,xregexpr+10,nchar(xfile)-4)  # e.g., "20240719_RA_L1"

nheaderlines <- scan(paste0(datdir,"/",xfile),nlines=1,sep=",")[1]
colnms <- scan(paste0(datdir,"/",xfile),nlines=1,sep=",",skip=nheaderlines-1,what="",strip.white=TRUE)
tmp <- read.csv(paste0(datdir,"/",xfile),skip=nheaderlines,col.names=colnms)
tmp[tmp==-9999] <- NA
dat <- tmp

# save data as RDS format
resultnm <- paste0(fltlabel,".rds")
saveRDS(dat,file=paste0(datdir,"/",resultnm))
print(paste(resultnm,"saved in",datdir))

if(plotTF){
hist(dat$ALTAGL,main=fltlabel)
plot(dat$Alt_Pressure ,dat$ALTAGL,pch=16,main=fltlabel)
plot(dat$StartTime_UTsec,dat$Alt_Pressure,pch=16,cex=0.5,main=fltlabel)
plot(dat$Longitude,dat$Latitude,pch=16,cex=0.5,main=fltlabel)
} # if(plotTF){

} # for(i in 1:length(xfiles)){
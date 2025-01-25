# Generates footprint figures from HRRR-STILT simulations for sites/cities in the CO2-USA merged GHG dataset, to examine COVID response
# Based on 'plot_foot_CH4_inversion_Uintah.r'
# January 3rd, 2022 by John C. Lin (John.Lin@utah.edu)

#########################
CITY <- "SLC"
SITE <- "UintahElementary"
outputdir <- paste0(SITE,"/by-id/")
pngfileTF <- FALSE; overwritepngfileTF <- TRUE
GoogleMapTF <- TRUE

#lat/lon range of plot
if(CITY=="BOSTON"){XLIMS <- c(-72.3, -70.8); YLIMS <- c(42.0, 43.0)}
if(CITY=="INDY"){ XLIMS <- c(-86.8, -85.6); YLIMS <- c(39.4, 40.6)}
if(CITY=="LA")  { XLIMS <- c(-119, -117); YLIMS <- c(33.5, 34.8)}
if(CITY=="NEC") { XLIMS <- c(-77.8, -76.2); YLIMS <- c(38.4, 39.6)}
#if(CITY=="SLC") { XLIMS <- c(-112.36, -111.7); YLIMS <- c(40.17, 41.4)}
if(CITY=="SLC") { XLIMS <- c(-113.1, -111.7); YLIMS <- c(40, 41.75)}
if(CITY=="TORONTO") {XLIMS <- c(-80.6,-78.6); YLIMS <- c(43.15, 44.5)}


#range of footprint plot (log10 limits)
ZLIMS<-c(-7,-1)
#########################

require(fields);require(ncdf4)
require(maps);require(RgoogleMaps)

#loop over the YYYYMMDD folders
YYYYMMDDs <- list.files(outputdir,pattern="202407")
YYYYMMDDs <- rev(YYYYMMDDs)
YYYYMMDDs <- sample(YYYYMMDDs,5)  #randomly sample a few to plot
for(i in 1:length(YYYYMMDDs)){
#for(i in c(1,50)){
#for(i in 13:length(YYYYMMDDs)){
  ncfilenm<-paste0(outputdir,"/",YYYYMMDDs[i],"/",YYYYMMDDs[i],"_foot.nc")
  if(!file.exists(ncfilenm)){print(paste("no footprint output for:",YYYYMMDDs[i]));next}

  #Generate footprint
  footfile<-nc_open(ncfilenm)
  lat<-ncvar_get(footfile,"lat");lon<-ncvar_get(footfile,"lon")
  foot.all<-ncvar_get(footfile,"foot")
  #sel.y<-lat>YLIMS[1]&lat<YLIMS[2];flat<-lat[sel.y]
  #sel.x<-lon>XLIMS[1]&lon<XLIMS[2];flon<-lon[sel.x]
  flat<-lat;flon<-lon
  footsum<-apply(foot.all,c(1,2),sum)   #sum over backward time
  #footsum<-footsum[sel.x,sel.y]
  foot.log<-log10(footsum);foot.log[footsum==0]<-NA

  #Read in trajectories 
  trajfile<-paste0(outputdir,"/",YYYYMMDDs[i],"/",YYYYMMDDs[i],"_traj.rds")
  if(!file.exists(trajfile)){print(paste("cannot find trajectory file",trajfile));next}
  dat.all<-readRDS(file=trajfile)
  receptor<-dat.all$receptor
  lon0<-receptor$long;lat0<-receptor$lati;agl0<-receptor$zagl
  receptor$run_time<-strftime(receptor$run_time,"%Y-%m-%d %H:%M",tz="UTC")
  dat<-dat.all$particle

  #Check whether PNG file exists or not
  if(pngfileTF){
    #pngfile <- paste0(outputdir,"/",YYYYMMDDs[i],"/",YYYYMMDDs[i],".png")
    pngfile <- paste0(YYYYMMDDs[i],".png")
    if(!overwritepngfileTF){if(file.exists(pngfile)){print(paste("PNG file already exists!; skip:",pngfile));next}}
  } #if(pngfileTF){

  lats<-dat$lati;lons<-dat$long
  #sel<-lats>YLIMS[1]&lats<YLIMS[2]
  #sel<-sel&(lons>XLIMS[1]&lons<XLIMS[2])
  #lats<-lats[sel];lons<-lons[sel]
if(!GoogleMapTF){
  dev.new();par(cex.axis=1.3,cex.lab=1.3,cex.main=1.1)
  image.plot(flon,flat,foot.log)
  #add trajectories
  points(lons,lats,pch=16,col="darkgray",cex=0.2)
  #add receptor location
  points(lon0,lat0,pch=17,col="black",cex=1.5)
  map("state",add=TRUE)  #;map("world",add=TRUE)
  map.cities(minpop=100000,cex=1.2)
  title(main=paste("Footprint (STILT model) at receptor (black triangle):\n",receptor$run_time,"UTC   ",
                     round(lon0,digits=3),"LON ",round(lat0,digits=3),"LAT ",round(agl0,digits=3),"m-AGL"),cex.main=1.1)
}else{
  #map defined by "zoom level"
  #  center<-c(mean(lats,na.rm=T),mean(lons,na.rm=T))
  #  center[1]<-center[1]-0.3;center[2]<-center[2]+0.2
  #  MyMap <- GetMap(center,zoom=9,maptype="terrain",frame=TRUE,GRAYSCALE=FALSE)
  #map defined by bounding box
  #bb <- qbbox(lat=lats, lon=lons)   #use trajectories to define map
  #bb <- qbbox(lat=flat, lon=flon)    #use footprint to define map
  bb <- qbbox(lat=YLIMS, lon=XLIMS)    #use footprint to define map
  MyMap <- GetMap.bbox(bb$lonR, bb$latR, maptype="terrain", frame=TRUE,GRAYSCALE=FALSE)
  PlotOnStaticMap(MyMap,mar=c(2,1,3,1))

  #Add trajectories
  #site.xy<- LatLon2XY.centered(MyMap, lat=as.numeric(dat$lati), lon=as.numeric(dat$long))
  site.xy<- LatLon2XY.centered(MyMap, lat=as.numeric(lats), lon=as.numeric(lons))
  #points(site.xy$newX,site.xy$newY,pch=16,col="black",cex=0.2)
  points(site.xy$newX,site.xy$newY,pch=16,col="darkgray",cex=0.2)

  #Add receptor location
  site.xy<- LatLon2XY.centered(MyMap, lat=as.numeric(lat0), lon=as.numeric(lon0))
  points(site.xy$newX,site.xy$newY,pch=17,col="black",cex=1.5)

  title(main=paste(SITE,"Site Footprint (STILT model) at receptor (black triangle):\n",receptor$run_time,"UTC   ",
                     round(lon0,digits=3),"LON ",round(lat0,digits=3),"LAT ",round(agl0,digits=3),"m-AGL"),cex.main=1.1)

  #Add footprint
  alpha<-0.45  #transparency
  lats.m<-matrix(rep(flat,length(flon)),byrow=T,ncol=length(flat))
  lons.m<-matrix(rep(flon,length(flat)),ncol=length(flat))
  xys<- LatLon2XY.centered(MyMap, lat=as.vector(lats.m), lon=as.vector(lons.m)) #converting lat/lon coordinates to map coordinates
  image.coords.x<-matrix(xys$newX,ncol=length(flat))
  image.coords.y<-matrix(xys$newY,ncol=length(flat))
  COLS<-c("white","darkgray","lightblue","RoyalBlue","darkgreen","yellow","orange","red") #colorscale
  colsc<-designer.colors(n=64,alpha=alpha,col=COLS)
  colsc.legend<-designer.colors(n=64,alpha=1.0,col=COLS)

  #generate image
  image(x=image.coords.x[,1],y=image.coords.y[1,],z=foot.log,
          col=colsc,ylim=c(MyMap$BBOX$ll[1],MyMap$BBOX$ur[1]),xlim=c(MyMap$BBOX$ll[2],MyMap$BBOX$ur[2]),add=T,zlim=ZLIMS)

  #Add legend 
  image.plot(x=image.coords.x[,1],y=image.coords.y[1,],z=foot.log,
        col=colsc.legend,ylim=c(MyMap$BBOX$ll[1],MyMap$BBOX$ur[1]),xlim=c(MyMap$BBOX$ll[2],MyMap$BBOX$ur[2]),add=T,
       legend.lab="log10(footprint)",legend.line=2,zlim=ZLIMS,legend.only=T,horizontal=T,legend.width=0.8)
} #if(!GoogleMapTF){

  if(pngfileTF){
    dev.copy(png,filename=pngfile);dev.off()
    print(paste(pngfile,"produced"))
  } #if(pngfileTF){
  nc_close(footfile)
  gc()
} #for(i in 1:length(YYYYMMDDs)){



#!/usr/bin/env Rscript
# STILT R Executable
# For documentation, see https://uataq.github.io/stilt/
# Ben Fasoli

print(paste("Starting time:",date()))
SITEs <- c("WBB","MtnMet","Syracuse","Hawthorne","Copperview","Herriman","RosePark","InlandPort","LakePark","TechCenter")
SITEs <- c(SITEs,"BryantMiddle","TheShop","EastHigh","Westminster","MtWire","UintahElementary") # these sites represent retroreflectors installed by NIST for dual-comb long-path instrument
for(ss in 1:length(SITEs)){
  SITE <- SITEs[ss]
  print(paste("--------------- SITE:",SITE,"---------------"))

  # User inputs ------------------------------------------------------------------
  project <- 'STILT_USOS'
  stilt_wd <- file.path('/uufs/chpc.utah.edu/common/home/u0791084/PROJECTS/USOS', project)
  #output_wd <- file.path(stilt_wd, 'out') 
  output_wd <- file.path('/uufs/chpc.utah.edu/common/home/u0791084/lin-group24/jcl/STILT_USOS_output',SITE) 
  lib.loc <- .libPaths()[1]

# Parallel simulation settings
n_cores <- 3
n_nodes <- 1
slurm   <- n_nodes > 1
slurm_options <- list(
  time      = '5:00:00',
  account   = 'lin-np',
  partition = 'lin-np'
)

# Simulation timing, yyyy-mm-dd HH:MM:SS (UTC)
YEARs <- (2024:2024)
#HHs <-  as.character(c(13:23,"00","01")) # subset of hours to simulate [UTC]
#HHs <-  as.character(c("00","06","12","18")) # subset of hours to simulate [UTC]
HHs <-  formatC(0:23,width=2,flag="0") # subset of hours to simulate [UTC]
run_times <- NULL
for(yy in 1:length(YEARs)){
  YEAR <- YEARs[yy]
  t_start <- paste0(YEAR,'-07-01 00:00:00')
  t_end   <- paste0(YEAR,'-08-31 23:00:00')
  run_times.sub <- seq(from = as.POSIXct(t_start, tz = 'UTC'),
                       to   = as.POSIXct(t_end, tz = 'UTC'), by   = 'hour')
  run_times <- c(run_times,run_times.sub)
} # for(yy in 1:length(YEARs)){
run_times <- as.POSIXct(as.numeric(run_times), tz='UTC', origin="1970-01-01")
run_times <- run_times[format(run_times,format="%H")%in%HHs]

# Receptor location(s)
# lati <- 40.5; long <- -112.0; zagl <- 5
if(SITE=="WBB"){
  lati <- 40.7634; long <- -111.848; zagl <- 36  
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="WBB"){
if(SITE=="MtnMet"){
  lati <- 40.7667; long <- -111.8284; zagl <- 4
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="MtnMet"){
if(SITE=="Syracuse"){
  lati <- 41.089; long <- -112.119; zagl <- 4
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="MtnMet"){
if(SITE=="Hawthorne"){
  lati <- 40.734477; long <- -111.872172; zagl <- 4    # Hawthorne Elementary School UDAQ monitoring site; inlet height guestimated by FRU height in Bares et al. [2019]
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="Hawthorne"){
if(SITE=="Copperview"){
  lati <- 40.59794; long <- -111.894; zagl <- 4    
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="Copperview"){
if(SITE=="Herriman"){
  lati <- 40.49639; long <- -112.036; zagl <- 4    
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="Herriman"){
if(SITE=="RosePark"){
  lati <- 40.79553; long <- -111.931; zagl <- 4    
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="RosePark"){
if(SITE=="InlandPort"){
  lati <- 40.80791; long <- -112.088; zagl <- 4    
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="InlandPort"){
if(SITE=="LakePark"){
  lati <- 40.7099; long <- -112.009; zagl <- 4    
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="LakePark"){
if(SITE=="TechCenter"){
  lati <- 40.77715; long <- -111.946; zagl <- 11   
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="TechCenter"){

# the following sites are for the retroreflectors installed by NIST for dual-comb long-path instrument
if(SITE=="BryantMiddle"){
  lati <- 40.76815; long <- -111.86923; zagl <- 15   
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="BryantMiddle"){
if(SITE=="TheShop"){
  lati <- 40.7605; long <- -111.88098; zagl <- 15   
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="TheShop"){
if(SITE=="EastHigh"){
  lati <- 40.75042; long <- -111.85516; zagl <- 10   
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="EastHigh"){
if(SITE=="Westminster"){
  lati <- 40.73258; long <- -111.85488; zagl <- 15   
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="Westminster"){
if(SITE=="MtWire"){
  lati <- 40.7706; long <- -111.79834; zagl <- 5   
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="MtWire"){
if(SITE=="UintahElementary"){
  lati <- 40.74201; long <- -111.846; zagl <- 10   
  receptors <- expand.grid(run_time = run_times, lati = lati, long = long, zagl = zagl, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
} # if(SITE=="UintahElementary"){

#only run those receptors when previous runs did NOT produce footprints
xfiles <- list.files(path=paste0(output_wd,"/footprints/"),full.names=FALSE)
YYYYMMDDHHmm <- format(receptors$run_time,format="%Y%m%d%H%M")
labels <- paste0(YYYYMMDDHHmm,"_",receptors$long,"_",receptors$lati,"_",receptors$zagl,"_foot.nc")
sel <- labels%in%xfiles
receptors <- receptors[!sel,]

# Footprint grid settings, must set at least xmn, xmx, ymn, ymx below
hnf_plume <- T
projection <- '+proj=longlat'
smooth_factor <- 1
time_integrate <- F
#xmn <- NA; xmx <- NA; ymn <- NA; ymx <- NA
#if(SITE=="Hawthorne") {xmn <- -112.36; xmx <- -111.7; ymn <- 40.17; ymx <- 41.4}
#if(TRUE) {xmn <- -112.36; xmx <- -111.7; ymn <- 40.17; ymx <- 41.4}    # SLV 
if(TRUE) {xmn <- -113.1; xmx <- -111.7; ymn <- 40; ymx <- 41.75}     # SLV + GSL + Utah Valley
xres <- 0.01
yres <- xres

# Meteorological data input
met_path           <- '/uufs/chpc.utah.edu/common/home/lin-group21/hrrr/hrrr/'
# met_file_format    <- '%Y%m%d.%Hz.hrrra'
met_file_format <- '%Y%m%d_%H'
met_subgrid_buffer <- 0.2
met_subgrid_enable <- TRUE 
met_subgrid_levels <- NA
n_met_min          <- 1

# Model control
n_hours       <- -24
numpar        <- 200
rm_dat        <- T
run_foot      <- TRUE
run_trajec    <- TRUE 
simulation_id <- NA
timeout       <- 3600
varsiwant     <- c('time', 'indx', 'long', 'lati', 'zagl', 'foot', 'mlht', 'dens',
                   'samt', 'sigw', 'tlgr')

# Transport and dispersion settings
capemin     <- -1
cmass       <- 0
conage      <- 48
cpack       <- 1
delt        <- 0
dxf         <- 1
dyf         <- 1
dzf         <- 0.01
efile       <- ''
emisshrs    <- 0.01
frhmax      <- 3
frhs        <- 1
frme        <- 0.1
frmr        <- 0
frts        <- 0.1
frvs        <- 0.01
hscale      <- 10800
ichem       <- 8
idsp        <- 2
initd       <- 0
k10m        <- 1
kagl        <- 1
kbls        <- 1
kblt        <- 5
kdef        <- 0
khinp       <- 0
khmax       <- 9999
kmix0       <- 150
kmixd       <- 3
kmsl        <- 0
kpuff       <- 0
krand       <- 4
krnd        <- 6
kspl        <- 1
kwet        <- 1
kzmix       <- 0
maxdim      <- 1
maxpar      <- numpar
mgmin       <- 10
mhrs        <- 9999
nbptyp      <- 1
ncycl       <- 0
ndump       <- 0
ninit       <- 1
nstr        <- 0
nturb       <- 0
nver        <- 0
outdt       <- 0
p10f        <- 1
pinbc       <- ''
pinpf       <- ''
poutf       <- ''
qcycle      <- 0
rhb         <- 80
rht         <- 60
splitf      <- 1
tkerd       <- 0.18
tkern       <- 0.18
tlfrac      <- 0.1
tout        <- 0
tratio      <- 0.75
tvmix       <- 1
veght       <- 0.5
vscale      <- 200
vscaleu     <- 200
vscales     <- -1
wbbh        <- 0
wbwf        <- 0
wbwr        <- 0
wvert       <- FALSE
w_option    <- 0
zicontroltf <- 0
ziscale     <- rep(list(rep(1, 24)), nrow(receptors))
z_top       <- 25000

# Transport error settings
horcoruverr <- NA
siguverr    <- NA
tluverr     <- NA
zcoruverr   <- NA

horcorzierr <- NA
sigzierr    <- NA
tlzierr     <- NA


# Interface to mutate the output object with user defined functions
before_trajec <- function() {output}
before_footprint <- function() {output}


# Source dependencies ----------------------------------------------------------
setwd(stilt_wd)
source('r/dependencies.r')


# Structure out directory ------------------------------------------------------
# Outputs are organized in three formats. by-id contains simulation files by
# unique simulation identifier. particles and footprints contain symbolic links
# to the particle trajectory and footprint files in by-id
#system(paste0('rm -r ', output_wd, '/footprints'), ignore.stderr = T)
#if (run_trajec) {
#  system(paste0('rm -r ', output_wd, '/by-id'), ignore.stderr = T)
#  system(paste0('rm -r ', output_wd, '/met'), ignore.stderr = T)
#  system(paste0('rm -r ', output_wd, '/particles'), ignore.stderr = T)
#}
for (d in c('by-id', 'particles', 'footprints')) {
  d <- file.path(output_wd, d)
  if (!file.exists(d))
    dir.create(d, recursive = T)
}


# Run trajectory simulations ---------------------------------------------------
stilt_apply(FUN = simulation_step,
            simulation_id = simulation_id,
            slurm = slurm, 
            slurm_options = slurm_options,
            n_cores = n_cores,
            n_nodes = n_nodes,
            before_footprint = list(before_footprint),
            before_trajec = list(before_trajec),
            lib.loc = lib.loc,
            capemin = capemin,
            cmass = cmass,
            conage = conage,
            cpack = cpack,
            delt = delt,
            dxf = dxf,
            dyf = dyf,
            dzf = dzf,
            efile = efile,
            emisshrs = emisshrs,
            frhmax = frhmax,
            frhs = frhs,
            frme = frme,
            frmr = frmr,
            frts = frts,
            frvs = frvs,
            hnf_plume = hnf_plume,
            horcoruverr = horcoruverr,
            horcorzierr = horcorzierr,
            hscale = hscale,
            ichem = ichem,
            idsp = idsp,
            initd = initd,
            k10m = k10m,
            kagl = kagl,
            kbls = kbls,
            kblt = kblt,
            kdef = kdef,
            khinp = khinp,
            khmax = khmax,
            kmix0 = kmix0,
            kmixd = kmixd,
            kmsl = kmsl,
            kpuff = kpuff,
            krand = krand,
            krnd = krnd,
            kspl = kspl,
            kwet = kwet,
            kzmix = kzmix,
            maxdim = maxdim,
            maxpar = maxpar,
            met_file_format = met_file_format,
            met_path = met_path,
            met_subgrid_buffer = met_subgrid_buffer,
            met_subgrid_enable = met_subgrid_enable,
            met_subgrid_levels = met_subgrid_levels,
            mgmin = mgmin,
            n_hours = n_hours,
            n_met_min = n_met_min,
            ncycl = ncycl,
            ndump = ndump,
            ninit = ninit,
            nstr = nstr,
            nturb = nturb,
            numpar = numpar,
            nver = nver,
            outdt = outdt,
            output_wd = output_wd,
            p10f = p10f,
            pinbc = pinbc,
            pinpf = pinpf,
            poutf = poutf,
            projection = projection,
            qcycle = qcycle,
            r_run_time = receptors$run_time,
            r_lati = receptors$lati,
            r_long = receptors$long,
            r_zagl = receptors$zagl,
            rhb = rhb,
            rht = rht,
            rm_dat = rm_dat,
            run_foot = run_foot,
            run_trajec = run_trajec,
            siguverr = siguverr,
            sigzierr = sigzierr,
            smooth_factor = smooth_factor,
            splitf = splitf,
            stilt_wd = stilt_wd,
            time_integrate = time_integrate,
            timeout = timeout,
            tkerd = tkerd,
            tkern = tkern,
            tlfrac = tlfrac,
            tluverr = tluverr,
            tlzierr = tlzierr,
            tout = tout,
            tratio = tratio,
            tvmix = tvmix,
            varsiwant = list(varsiwant),
            veght = veght,
            vscale = vscale,
            vscaleu = vscaleu,
            vscales = vscales,
            w_option = w_option,
            wbbh = wbbh,
            wbwf = wbwf,
            wbwr = wbwr,
            wvert = wvert,
            xmn = xmn,
            xmx = xmx,
            xres = xres,
            ymn = ymn,
            ymx = ymx,
            yres = yres,
            zicontroltf = zicontroltf,
            ziscale = ziscale,
            z_top = z_top,
            zcoruverr = zcoruverr)

print(paste("Ending time:",date()))

} # for(ss in 1:length(SITEs)){
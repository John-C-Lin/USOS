# Organizes STILT output by date
# Feb. 3rd, 2025 by John C. Lin (John.Lin@utah.edu)

#----------------#
outputdir <- "out/"
removedirTF <- FALSE
#----------------#

outdirs <- list.dirs(paste0(outputdir,"/by-id/"),full.names=FALSE) # e.g., "202407191453_-105.10216977_39.90150991_68"
# check to see whether there is valid trajectory output;  if not, then delete directory
trajfiles <- list.files(paste0(outputdir,"/particles/"))   # e.g., "202407191453_-105.10216977_39.90150991_68_traj.rds"
trajfiles <- substring(trajfiles,1,nchar(trajfiles)-9)   # remove '_traj.rds' from the end
SEL <- outdirs%in%trajfiles
print(paste("Number of receptors without trajectories:",sum(!SEL)))
if(removedirTF){
  files2remove <- paste0(outputdir,"/by-id/",outdirs[!SEL]) 
  print("removing directories without trajectories:")
  #file.remove(files2remove)
} # remove receptor directories that do not contain trajectories

# loop through each date, creating separate tar file
YYYYMMDDs <- substring(outdirs,1,8)
YYYYMMDDs <- unique(YYYYMMDDs[nchar(YYYYMMDDs)==8])
for(i in 1:length(YYYYMMDDs)){
  YYYYMMDD <- YYYYMMDDs[i]
  sel <- substring(outdirs,1,8)==YYYYMMDD
  tarfile <- paste0(YYYYMMDD,".tar.gz")
  command <- paste0("tar cfz ",tarfile," ",outputdir,"/by-id/",YYYYMMDD,"*")
  system(command)
  print(command)
  gc()
} # for(i in 1:length(YYYYMMDDs)){
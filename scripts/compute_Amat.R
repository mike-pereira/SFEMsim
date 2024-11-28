################################################################################

## PROJECTION MATRICES

################################################################################

#-------------------------------------------------------------------------------

#### IMPORT PACKAGES
library(SFEMsim)
library(Matrix)
library(parallel)
library(rgl)
library(Rvcg)


## Folder containing the data (nested sphere meshes)
datFolder="/home/mpereira/Documents/Work/data/Meshes/spheres/"


#--------------------------------------------------

##' Function to import a sphere mesh computed at a given discretization level
##' Arguments:
##' - ind : Discretization level: the sphere meshes are nested meshes and the argument ind indicates with nesting interation to consider
##' Returns : 
##' - A list contining the mesh information
importDatSphere<-function(ind=7){
  ### Import triangulation
  dat=vcgImport(paste0(datFolder,paste0("sphere_",ind,".ply")),updateNormals = F,readcolor = F,clean = F,silent = F)
  ### Nodes
  L=1
  nodeMat=L*t(dat$vb)[,-4]/mean(rowSums(t(dat$vb)[,-4]**2))**0.5 ## Node coordinates in Cartesian coordinates
  nodeMatSph=cbind(atan2((nodeMat[,1]**2+nodeMat[,2]**2)**0.5,nodeMat[,3]), atan2(nodeMat[,2],nodeMat[,1])) ## Node coordinates in Spherical coordinates
  ### Triangles
  triMat=t(dat$it)
  ## Number of nodes
  N=nrow(nodeMat)
  return(list(N=N,nodeMat=nodeMat,triMat=triMat,nodeMatSph=nodeMatSph,mesh=dat))
}


##' Function to project of a set of points in a set of mesh triangles and compute their barycentric coordinates
##' Arguments:
##' - pts : Matrix containing the coordinates of the points to project
##' - nodeMat : Matrix containing the coordinates of the mesh nodes
##' - triMat : Matrix containing the node indices of each mesh triangle
##' - triBar : Coordinates of the barycenter of each triangle
##' - n_cores : Number of cores to use for parallelism
##' Returns : 
##' - A matrix containing for each point the index of the mesh triangle on which it was projected, and its 3 barycentric coordinates
computeBarCoordFast<-function(pts,nodeMat,triMat,triBar,n_cores=2){
  nn=floor(nrow(pts)/n_cores)
  N=nrow(pts)
  nbtot=N %/% nn
  rmd = (nbtot * nn+1):N
  f <- function(i) {
    ## Function for the package SFEMsim used to compute the projection of (a set of) points
    ## on mesh triangle. Each point is projected to the triangle which barycenter is the closest 
    ## to the point, and then a ray tracing algorithm is used to compute the barycenteric coordinates of the points
    pr=sphBarCoordBT(pts[(i-1)*nn+1:nn,],nodeMat = nodeMat, 
                     triMat = triMat-1,triBar=triBar)
    return(pr)
  }
  res=mclapply(1:nbtot,f,mc.cores = n_cores)
  if((nbtot * nn+1)<N){
    res[[nbtot+1]]=sphBarCoordBT(pts[rmd,],nodeMat = nodeMat, triMat = triMat-1,triBar=triBar)
  }
  res=do.call(rbind,res)
  projPts=cbind(res[,1]+1,1-res[,3]-res[,4],res[,3],res[,4])
  return(projPts)
}



##' Function to round a value to the nearest multiple of an inverse power of 2
##' Arguments:
##' - val : Vector of values to round
##' - lvl_max : The values are rounded to the nearest multiple of 2^{-lvl_max} 
##' Returns : 
##' - A vector of rounded values
roundCoord<-function(val,lvl_max=8){
  
  val_tmp=val
  val_f=rep(-1,length(val_tmp))
  id_all=1:length(val_tmp)
  
  v_ref=seq(from=0,to=1,by=1/2**lvl_max)
  eps=1/2**(lvl_max)
  for (i in 1:length(v_ref)) {
    id=which((val_tmp >= (v_ref[i] - eps/2)) * (val_tmp < (v_ref[i] + eps/2)) ==1)
    if(length(id)>0){
      val_f[id_all[id]]=v_ref[i]
      val_tmp=val_tmp[-id]
      id_all=id_all[-id]
    }
  }
  
  return(val_f)
}

##' Function to compute a projection patrix between two nested sphere meshes
##' Arguments:
##' - nodeMat_fine,nodeMat_coarse, : Vector of values to round
##' - nodeMat_fine,nodeMat_coarse : Matrix containing the coordinates of the mesh nodes for the fine and coarse mesh
##' - triMat_coarse : Matrix containing the node indices of each  triangle of the coarse mesh
##' - lvl_max: The values of entries of the projection are rounded to the nearest multiple of 2^{-lvl_max} 
##' - n_cores : Number of cores to use for parallelism
##' Returns : 
##' - The projection matrix between the meshes
compute_Amat <- function(nodeMat_fine,nodeMat_coarse,triMat_coarse,n_cores=2,lvl_max=10) {
  

  triBar_coarse=(nodeMat_coarse[triMat_coarse[,1],]+nodeMat_coarse[triMat_coarse[,2],]+nodeMat_coarse[triMat_coarse[,3],])/3
  projPts = computeBarCoordFast(nodeMat_fine,nodeMat_coarse,triMat_coarse,triBar_coarse,n_cores)
 
  fb <- function(i) {
    idTri=triMat_coarse[projPts[i,1],]
    return(cbind(i,idTri,projPts[i,c(2,3,4)]))
  }
  Amat=mclapply(1:nrow(projPts),fb,mc.cores = n_cores)
  Amat=do.call(rbind,Amat)
  
  val=roundCoord(Amat[,3],lvl_max = lvl_max)
  return(sparseMatrix(i=Amat[,2],j=Amat[,1],x=val,dims = c(nrow(nodeMat_coarse),nrow(nodeMat_fine))))
  
}



#--------------------------------------------------
#--------------------------------------------------
#--------------------------------------------------


## Compute the projection matrices between two consecutive nested meshes for levels 2 to 9
A_list=list()
idFine_list=3:10
for(i in 1:length(idFine_list)){
  cat(i,"/",length(idFine_list),"\n")
  idFine=idFine_list[i]
  sphCoarse=importDatSphere(idFine-1)
  sphFine=importDatSphere(idFine)
  A_list[[i]]=compute_Amat(sphFine$nodeMat,sphCoarse$nodeMat, 
                           triMat = sphCoarse$triMat,n_cores = 15)
}

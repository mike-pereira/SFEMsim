################################################################################

## SCRIPT FOR SIMULATIONS ON THE SPHERE

################################################################################

#-------------------------------------------------------------------------------

#### IMPORT PACKAGES
library(SFEMsim)
library(Matrix)
library(rgl)
library(sf)
library(terra)
library(Rvcg)

## Folder containing the data (meshes, coastlines)
datFolder="/Users/erikjans/Downloads/data_fem/"
datFolder="/home/mpereira/Documents/Work/data/SFEM/"

#-------------------------------------------------------------------------------

## Function to plot borders
### Arguments
### geom : List of polygons defining the borders
### ... : Additional arguements to be passed to the lines3d function (eg. `linewidth=3`,...)
plotBorders<-function(geom,...){
  for(k in 1:length(geom)){
    u=geom[[k]][[1]]*pi/180 ## Extract vertex coordinates of polygon
    theta=u[,1] ## long
    delta=u[,2] ## lat
    lines3d(cos(delta)*cos(theta),cos(delta)*sin(theta),sin(delta),col="black",...) ## Plot 3D lines
  }
}

#-------------------------------------------------------------------------------


#### IMPORT COAST LINES
clV=vect(paste0(datFolder,"coast/WorldCoastlines.shp"))
coastLines=lapply(geom(clV,list=T), function(l){
  uu=as.matrix(l[[1]][[1]],byrow=T)
  list(uu)})


#-------------------------------------------------------------------------------

#### IMPORT SPHERE MESH

### Import triangulation
dat=vcgImport(paste0(datFolder,"mesh/sphere/sphere_7.ply"),updateNormals = F,readcolor = F,clean = F,silent = F)

### Nodes
L=1
nodeMat=L*t(dat$vb)[,-4]/mean(rowSums(t(dat$vb)[,-4]**2))**0.5 ## Node coordinates in Cartesian coordinates
nodeMatSph=cbind(atan2((nodeMat[,1]**2+nodeMat[,2]**2)**0.5,nodeMat[,3]), atan2(nodeMat[,2],nodeMat[,1])) ## Node coordinates in Spherical coordinates

### Triangles
triMat=t(dat$it)

## Number of nodes
N=nrow(nodeMat)

## Create sphere mesh for plots (with a radius smaller than 1, otherwise we can't see the borders)
mshplt=mesh3d(vertices = 0.999*t(nodeMat),triangles = t(triMat))


#-------------------------------------------------------------------------------

#### CREATE SELECTION OF POINTS IN LAND

## Compute coordinates of mesh nodes in spherical coordinates, using the same convention as the coast lines
nodeMatSphBorders=cbind(atan2(nodeMat[,2],nodeMat[,1]),atan2(nodeMat[,3],(nodeMat[,1]**2+nodeMat[,2]**2)**0.5))*180/pi

## Find points in land
mplg=st_multipolygon(coastLines) ## Polygons defining land
mpts=st_as_sf(as.data.frame(nodeMatSphBorders),coords = c("V1","V2")) ## Mesh nodes
inter_mpts=st_intersects(mpts,mplg) ## Find intersection
selPtInLand=(lengths(inter_mpts)==1) ## Binary variable of points in land

## Plot to check
inter_mpts=nodeMatSphBorders[selPtInLand,] ## Spherical coordinates of points in land
inter_pts=cbind(cos(inter_mpts[,2]*pi/180)*cos(inter_mpts[,1]*pi/180),
                cos(inter_mpts[,2]*pi/180)*sin(inter_mpts[,1]*pi/180),
                sin(inter_mpts[,2]*pi/180)) ## Cartesian coordinates of points in land
points3d(inter_pts,cex=0.1)

#-------------------------------------------------------------------------------

#### METRIC

## Compute coefficients of metric tensor in each triangle (Metric = PullBack of Euclidean metric)
triMetricMat=triMetricPull2d(nodeMat,triMat)

#-------------------------------------------------------------------------------

#### DIFFERENTIAL OPERATOR

## Expression of the function defining the diffusion directions <- Change here to get other diffusive behaviors
# f_expr=expression(sin(x[1])*sin(x[2]))
f_expr=expression(sin(theta)**2*(2*cos(theta)*cos(phi)+sin(phi)))


## Symbolic differentiation of the function (to compute function and its gradient)
f_symbol <- deriv(f_expr, namevec = c("theta", "phi"), function.arg = c("theta", "phi"))
## Compute the function values from symbolic expression: x =(theta, phi)
eval_f<-function(x){
  return(f_symbol(x[1],x[2]))
}
## Compute the gradient of the function (wrt spherical coordinates) from symbolic expression: x =(theta, phi)
grad_f<-function(x){
  return(attr(f_symbol(x[1],x[2]),"gradient"))
}

## Function that coefficients of the tensor associated with the diffusion matrix
GfuncSph=function(x){

  ## Compute gradient of f wrt to spherical coordinates
  v=grad_f(x)

  ## Compute the diffusion directions
  u=c(v[1],v[2])*sin(x[1])/((sin(x[1])**2*v[1]**2+1*v[2]**2))**0.5
  uorth=c(-v[2],(sin(x[1])**2)*v[1])*1.0/((sin(x[1])**2)*v[1]**2+v[2]**2)**0.5

  ## Weights associated to each direction
  rho=c(1,100)

  return((1/rho[1])*u%*%t(u)+(1/rho[2])*uorth%*%t(uorth))
}

## Compute coefficients of the diffusion tensor on each triangle
triDiffMat=triMetricSph(nodeMat, triMat, GfuncSph = GfuncSph)

#-------------------------------------------------------------------------------

####  POTENTIAL FUNCTION

## Values at the nodes
nodePotCoefs=(1+10000*(selPtInLand))/4

## Values at the triangles
triPotCoeffs=(nodePotCoefs[triMat[,1]]+nodePotCoefs[triMat[,2]]+nodePotCoefs[triMat[,3]])/3


#-------------------------------------------------------------------------------


####  FEM MATRICES

## Build FEM matrices
FEMatList=matFEM2d(nodeMat,triMat,triMetricMat,triPotCoeffs,triDiffMat)

## Print smallest and highest eigenvalues: if too big (>1e10), your diffusion has singularities
cat(c("Eigenvalue interval for Chebyshev approximation :", FEMatList$Eig))


#-------------------------------------------------------------------------------

#### SIMULATION

### Gamma function
nu=1
gamma<-function(x){
  x^(-(nu+1)/2)
}

## Simulation
set.seed(1234)
WN=matrix(rnorm(N),ncol=1)
z_sim=simPSD(gamma,FEMatList,WN=WN,tolCoefs = 10**(-3))


#-------------------------------------------------------------------------------

#### PLOT

## Values to plot
zplt=z_sim[,1]

## Plot simulation result
rgl.clear()
plotOnMesh(zplt,mshplt)
plotBorders(coastLines,lwd=2)
view3d(theta=0,phi=-80,zoom=0.65)

#-------------------------------------------------------------------------------

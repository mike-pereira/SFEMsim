################################################################################

## SCRIPT FOR SIMULATIONS ON THE BRAIN SURFACE

################################################################################

#-------------------------------------------------------------------------------

#### IMPORT PACKAGES
library(SFEMsim)
library(Matrix)
library(fields)
library(rgl)
library(Rvcg)

## Folder containing the data (meshes, coastlines)
datFolder="/Users/erikjans/Downloads/data_fem/"
datFolder="/home/mpereira/Documents/Work/data/SFEM/"

#-------------------------------------------------------------------------------

#### IMPORT SPHERE MESH

### Import triangulation
dat=vcgImport(paste0(datFolder,"mesh/brain.ply"),updateNormals = F,readcolor = F,clean = F,silent = F)

### Nodes
L=1
nodeMat=L*t(dat$vb)[,-4]/mean(rowSums(t(dat$vb)[,-4]**2))**0.5 ## Node coordinates in Cartesian coordinates

### Triangles
triMat=t(dat$it)

## Number of nodes
N=nrow(nodeMat)

## Create sphere mesh for plots (with a radius smaller than 1, otherwise we can't see the borders)
mshplt=mesh3d(vertices = 0.999*t(nodeMat),triangles = t(triMat))

#-------------------------------------------------------------------------------

#### METRIC

## Compute coefficients of metric tensor in each triangle (Metric = PullBack of Euclidean metric)
triMetricMat=triMetricPull2d(nodeMat,triMat)

#-------------------------------------------------------------------------------

#### DIFFERENTIAL OPERATOR

## Take A = Identity matrix
triDiffMat=triMetricMat

#-------------------------------------------------------------------------------

####  POTENTIAL FUNCTION

## Potential function
potf<-function(x){
  if (x[2]< -0.1){
    1.0
  }else{
    10000
  }
}

## Values at the triangles
triPotCoeffs=triPot(nodeMat,triMat,potf,coordSph = FALSE,useParallel = T)

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
view3d(theta=0,phi=-80,zoom=0.65)

#-------------------------------------------------------------------------------

################################################################################

## SCRIPT FOR 1D SIMULATIONS

################################################################################

#-------------------------------------------------------------------------------


#### IMPORT PACKAGES

library(SFEMsim)
library(Matrix)


#-------------------------------------------------------------------------------

#### CREATE CIRCLE MESH

# Number of nodes
N=10000

# Node coordinates
thetaN=seq(from=0,to=2*pi,length.out=N+1)[1:N]
rhoN=1/(1+0.25*cos(5*thetaN))**2
rhoN=rhoN/max(rhoN)
nodeMat=cbind(rhoN*cos(thetaN),rhoN*sin(thetaN))

# "Triangles" (or rather lines since we are 1D)
triMat=cbind(1:N,c(2:N,1))

# Centroids
triCenter=(nodeMat[triMat[,1],]+nodeMat[triMat[,2],])/2

## Plot surface
plot1d(nodeMat[,1],nodeMat,triMat,colorPlot=TRUE,ttl = "Plot of 1D the surface")

#-------------------------------------------------------------------------------

#### METRIC

## Compute coefficients of metric tensor in each triangle (Metric = PullBack of Euclidean metric)
triMetricMat=triMetricPull1d(nodeMat,triMat)


#-------------------------------------------------------------------------------

#### DIFFUSION MATRIX A

## Diffusion 1D "matrix" in circular coordinates
afunc_circ<-function(theta){
  (1+99*exp(-(theta-pi)**2/0.3**2))
  # (1+9*(1-exp(-(theta-pi)**2)))*300
}
plot(thetaN,afunc_circ(thetaN),type="l",
     main = "Diffusion coefficient",xlab = "Node circular coordinate", ylab="Value")

## Diffusion 1D "matrix" in Cartesian coordinates
afunc_cart<-function(x){
  theta=atan2(x[2],x[1])
  theta=theta+(2*pi)*(theta<0)
  return(afunc_circ(theta))
}

## Compute coefficients of the diffusion tensor on each segment
triDiffMat=as.matrix(apply(triCenter, 1, afunc_cart))


#-------------------------------------------------------------------------------

####  POTENTIAL FUNCTION V

## Potential function in circular coordinates
Vfunc_circ<-function(theta){
  return(10*(theta**0))
}
plot(thetaN,Vfunc_circ(thetaN),type="l",
     main = "Potential function V",xlab = "Node circular coordinate", ylab="Value")

## Potential function in Cartesian coordinates
Vfunc_cart<-function(x){
  theta=atan2(x[2],x[1])
  theta=theta+(2*pi)*(theta<0)
  return(Vfunc_circ(theta))
}

triPotCoeffs=as.matrix(apply(triCenter, 1, Vfunc_cart))

#-------------------------------------------------------------------------------


####  FEM MATRICES

## Build FEM matrices
FEMatList=matFEM1d(nodeMat,triMat,triMetricMat,triPotCoeffs,triDiffMat)

## Print smallest and highest eigenvalues: if too big (>1e10), your diffusion may have singularities
cat(c("Eigenvalue interval for Chebyshev approximation :", FEMatList$Eig))


#-------------------------------------------------------------------------------

#### SIMULATION

### Gamma function
nu=0.2
gamma<-function(x){
  sin(x)*x^(-(nu+1)/2)
}

## Simulation
set.seed(1234)
WN=matrix(rnorm(N),ncol=1)
z_sim=simPSD(gamma,FEMatList,WN=WN,tolCoefs = 10**(-4))


#-------------------------------------------------------------------------------

#### PLOT

## Values to plot
zplt=z_sim[,1]

## Plot simulation result as deformed surface
plot1d_deform(zplt,nodeMat,triMat,r=0.15,ttl = "Plot of the simulation as a deformed surface")



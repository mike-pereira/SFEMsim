################################################################################

## NUMERICAL EXPERIMENT: STRONG ERROR ON THE CIRCLE

################################################################################

#-------------------------------------------------------------------------------

#### IMPORT PACKAGES
library(SFEMsim)
library(Matrix)
library(Rvcg)
library(parallel)

#------------------------------------------------------

# On cluster
## export folder
export=T
expFolder="../export/"


#--------------------------------------------------

## Parameters

# Number of nodes of the fine mesh
N_fine=2**20

## Discretization levels
lvl_list=3:6

## Number of simulations
nsim=1

## Seed
seed=1234

## List of cases
alpha_l=c(0.5,1.05,1.5)

## Tolerance on coefficients
tolCoeffs=10**(-12)

## Use Cholesky
useChol=T

## Number of cores
mc.cores=min(nsim,15)

#-------------------------------------------------------------------------------

## Diffusion 1D "matrix" in circular coordinates
afunc_circ=NULL

## Cases
case=2
if(case==1){
  ## Potential function in circular coordinates
  Vfunc_circ=NULL
  ### Gamma function
  kappa=20
  gamma_alpha<-function(x,alpha,theta=0.9){
    (kappa**2+x)^(-alpha)
  }
}else{
  ## Potential function in circular coordinates
  Vmin=10**4
  Vfunc_circ<-function(theta){
    return(Vmin*(2*(theta<3*pi/2)*(theta>pi/2)+1))
  }
  ### Gamma function
  gamma_alpha<-function(x,alpha,theta=0.9){
    (Vmin**(alpha))*(x+cos(pi*theta)*sqrt(x))^(-alpha)
  }
}


#--------------------------------------------------

##' Function to create a coarse mesh from a fine mesh
##' Arguments:
##' - mesh_fine : A list containing the number of nodes, their coordinates and a table containing the adjacency information for the fine mesh
##' - lvl : Level or coarsifucation : we keep one point avery 2**lvl points
##' Returns : 
##' - A list containing the number of nodes, their coordinates, a table containing the adjacency information, and the projection matrix A for the coarse mesh
extract_coarse<-function(mesh_fine,lvl=1){
  l=2**lvl
  N_fine=mesh_fine$N
  id_coarse=seq(from=0,to=N_fine-1,by=l)
  N_coarse=length(id_coarse)
  nodeMat=mesh_fine$nodeMat[id_coarse+1,]
  theta=mesh_fine$theta[id_coarse+1]
  triMat=cbind(1:N_coarse,c(2:N_coarse,1))
  tp=NULL
  for (i in 0:(l-1)) {
    if(i==0){
      tp=rbind(tp,cbind((id_coarse/l)%%(N_coarse) + 1, (id_coarse)%%(N_fine) + 1, 1))
    }else{
      tp=rbind(tp,cbind((id_coarse/l)%%(N_coarse) + 1, (id_coarse + i)%%(N_fine) + 1, 1-i/l))
      tp=rbind(tp,cbind((id_coarse/l)%%(N_coarse) + 1, (id_coarse - i)%%(N_fine) + 1, 1-i/l))
    }
  }
  A=sparseMatrix(i=tp[,1],j=tp[,2],x=tp[,3])
  return(list(N=N_coarse,nodeMat=nodeMat,triMat=triMat,A=A,theta=theta,h=theta[2]-theta[1]))
}

##'Function to compute the FEM matrices for a mesh
##' Arguments:
##' - mesh : A list containing the number of nodes, their coordinates and a table containing the adjacency information for the fine mesh
##' Returns : 
##' - The same list with an added slots for the FEM matrices
addFEM<-function(mesh,afunc_circ,Vfunc_circ,computeChol = FALSE){
  #### METRIC
  ## Compute coefficients of metric tensor in each triangle (Metric = PullBack of Euclidean metric)
  triMetricMat=triMetricPull1d(mesh$nodeMat,mesh$triMat)
  
  # Centroids
  triCenter=(mesh$nodeMat[mesh$triMat[,1],]+mesh$nodeMat[mesh$triMat[,2],])/2
  
  #-------------------------------------------------------------------------------
  
  #### DIFFUSION MATRIX A
  if(is.null(afunc_circ)){
    triDiffMat=NULL
  }else{
    ## Diffusion 1D "matrix" in Cartesian coordinates
    afunc_cart<-function(x){
      theta=atan2(x[2],x[1])
      theta=theta+(2*pi)*(theta<0)
      return(afunc_circ(theta))
    }
    ## Compute coefficients of the diffusion tensor on each segment
    triDiffMat=as.matrix(apply(triCenter, 1, afunc_cart))
  }
  
  #-------------------------------------------------------------------------------
  
  ####  POTENTIAL FUNCTION V
  if(is.null(Vfunc_circ)){
    triPotCoeffs=NULL
  }else{
    ## Potential function in Cartesian coordinates
    Vfunc_cart<-function(x){
      theta=atan2(x[2],x[1])
      theta=theta+(2*pi)*(theta<0)
      return(Vfunc_circ(theta))
    }
    triPotCoeffs=as.matrix(apply(triCenter, 1, Vfunc_cart))
  }
  
  
  #-------------------------------------------------------------------------------
  
  
  ####  FEM MATRICES
  ## Build FEM matrices
  FEMatList=matFEM1d(mesh$nodeMat,mesh$triMat,triMetricMat,triPotCoeffs,triDiffMat,massLumping = FALSE)
  ## Print smallest and highest eigenvalues: if too big (>1e10), your diffusion may have singularities
  cat(c("Eigenvalue interval for Chebyshev approximation :", FEMatList$Eig))
  
  Cdh=Diagonal(mesh$N,1/diag(FEMatList$Scale))
  if(computeChol){
    Chol=expand(Cholesky(FEMatList$Mass))
  }else{
    Chol=NULL
  }
  
  ### Add result to list
  new_mesh=mesh
  new_mesh[["FEMatList"]]=FEMatList
  new_mesh[["Cdh"]]=Cdh
  new_mesh[["Chol"]]=Chol
  return(new_mesh)
}


#-------------------------------------------------------------------------------

#### CREATE FINE SPHERE MESH

# Node coordinates
thetaN=seq(from=0,to=2*pi,length.out=N_fine+1)[1:N_fine]
rhoN=1
nodeMat_fine=cbind(rhoN*cos(thetaN),rhoN*sin(thetaN))

# "Triangles" (or rather lines since we are 1D)
triMat_fine=cbind(1:N_fine,c(2:N_fine,1))

## Create a list containing all the info
sphFine=list(N=N_fine,nodeMat=nodeMat_fine,
             triMat=triMat_fine,theta=thetaN,h=thetaN[2]-thetaN[1])

## Compute FEM for fine mesh
sphFine=addFEM(sphFine,afunc_circ,Vfunc_circ,computeChol = useChol)

#-------------------------------------------------------------------------------

#### SIMULATION

## Create white noise at the fine level
set.seed(seed)
WN_fine=matrix(rnorm(nsim*sphFine$N),ncol=nsim)
z_sim_fine=array(dim=c(sphFine$N,nsim,length(alpha_l)))
err_list=array(dim=c(length(alpha_l),length(lvl_list)))
n_list=rep(0,length(lvl_list))
h_list=rep(0,length(lvl_list))
colnames(err_list)=paste0("lvl=",lvl_list)
rownames(err_list)=paste0("alpha=",alpha_l)
for(i in 1:length(alpha_l)){
  
  print(paste("Sampling fine level",": Case ",i,"/",length(alpha_l)))
  alpha=alpha_l[i]
  gamma<-function(x){
    gamma_alpha(x,alpha)
  }
  
  ## Compute fine samples
  z_sim_fine=mclapply(1:nsim,function(i){
    if(useChol){
      sim=simPSDChol(gamma,sphFine$FEMatList,WN=as.matrix(WN_fine[,i]),tolCoefs = tolCoeffs,verbose = F)
    }else{
      sim=simPSD(gamma,sphFine$FEMatList,WN=as.matrix(WN_fine[,i]),tolCoefs = tolCoeffs,verbose = F)
    }
    return(sim)
  },mc.cores=mc.cores)
  z_sim_fine=do.call(cbind,z_sim_fine)
  
  ## Compute simulations on coarse levels
  for(k in 1:length(lvl_list)){
    
    print(paste("Projecting level",k,"/",length(lvl_list)))
    lvl=lvl_list[k]
    
    ## Create coarse mesh from fine mesh
    sphCoarse=extract_coarse(sphFine,lvl=lvl)
    ## Compute FEM for coarse mesh
    sphCoarse=addFEM(sphCoarse,afunc_circ,Vfunc_circ,computeChol = useChol)
    n_list[k]=sphCoarse$N
    h_list[k]=sphCoarse$h
    
    #-------------------------------------------------------------------------------
    
    z_sim_coarse=mclapply(1:nsim,function(i){
      if(useChol){
        ## Project white noise at the coarse level (without the mass lumping approximation)
        #WN_coarse=as.matrix(solve(sphCoarse$Chol$L,sphCoarse$Chol$P %*% sphCoarse$A %*% t(sphFine$Chol$P) %*% sphFine$Chol$L %*% WN_fine))
        WN_coarse=projectNoise(sphFine$FEMatList$Mass,sphCoarse$FEMatList$Mass,sphCoarse$A,as.matrix(WN_fine[,i]))
        ## Compute samples
        z_tmp=simPSDChol(gamma,sphCoarse$FEMatList,WN=WN_coarse,tolCoefs = tolCoeffs,verbose = F)
      }else{
        ## Project white noise at the coarse level (using the mass lumping approximation)
        WN_coarse=as.matrix(sphCoarse$FEMatList$Scale %*% sphCoarse$A %*% sphFine$Cdh %*% as.matrix(WN_fine[,i]))
        ## Compute samples
        z_tmp=simPSD(gamma,sphCoarse$FEMatList,WN=WN_coarse,tolCoefs = tolCoeffs,verbose = F)
      }
      return(z_tmp)
    },mc.cores=mc.cores)
    z_sim_coarse=do.call(cbind,z_sim_coarse)
    
    ## reproject on fine mesh
    z_sim_coarse=(t(sphCoarse$A)%*%(z_sim_coarse))
    
    ## Compute error
    err=colSums((z_sim_coarse - z_sim_fine)* (sphFine$FEMatList$Mass %*% (z_sim_coarse - z_sim_fine)))
    err=mean(err)**0.5
    err_list[i,k]=err
    plot(z_sim_fine[,1],type="l",main=paste0("lvl=",lvl));points(z_sim_coarse[,1],type="l",col="red")
  }
}


#-------------------------------------------------------------------------------

## Save results
if(export){
  fileExp=paste0(expFolder,"results_circle_nstat.RData")
  # save.image(file=fileExp)
  save(list=c("err_list","n_list","h_list","alpha_l","lvl_list"),file=fileExp)
  print(paste0("Results saved at location :", fileExp))
}
#-------------------------------------------------------------------------------

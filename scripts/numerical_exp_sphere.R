################################################################################

## NUMERICAL EXPERIMENT: STRONG ERROR ON THE SPHERE

################################################################################

#-------------------------------------------------------------------------------

#### IMPORT PACKAGES
library(SFEMsim)
library(Matrix)
library(Rvcg)
library(parallel)

#------------------------------------------------------

# On cluster
## Load projection matrices (precomputed using the compute_Amat.R script)
load(paste0("../../../data/projMat/A_list.RData"))
## Folder containing the data (sphere meshes)
datFolder="../../../data/Meshes/spheres/"
## export folder
export=T
expFolder="../export/"


#--------------------------------------------------

## Parameters

## Level of discetization
idFine=10
lvl_list=5:(idFine-1)
lvl_list=idFine-1

## Number of simulations
nsim=1

## Seed
seed=1234

## List of cases
alpha_l=1/2+c(1/4,1/2+1/4,1+1/2+1/4)

## Tolerance on coefficients
tolCoeffs=10**(-12)

## Use Cholesky
useChol=F

## Number of cores
mc.cores=min(nsim,15)

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

##' Function to the mesh size h
##' Arguments:
##' - nodeMat : Matrix containing the coordinates of the mesh nodes
##' - triMat : Matrix containing the adjacency information of the mesh
##' Returns : 
##' - The mesh size h
compute_h<-function(nodeMat,triMat){
  x0=nodeMat[triMat[,1],]
  x1=nodeMat[triMat[,2],]
  x2=nodeMat[triMat[,3],]
  xx=cbind(x0,x1,x2)
  hh=t(apply(xx,1,function(v){
    c(norm(v[1:3]-v[7:9]),norm(v[4:6]-v[7:9]),norm(v[1:3]-v[4:6]))
  }))
  return(max(hh))
}


##' Function to compute the FEM matrices for a mesh
##' Arguments:
##' - mesh : A list containing the number of nodes, their coordinates and a table containing the adjacency information for the fine mesh
##' - computeChol : Whether to compute the Cholesky decomposition of the mass matrix
##' Returns : 
##' - The same list with an added slots for the FEM matrices
addFEM<-function(mesh,computeChol=TRUE){
  
  #### METRIC
  ## Compute coefficients of metric tensor in each triangle (Metric = PullBack of Euclidean metric)
  triMetricMat=triMetricPull2d(mesh$nodeMat,mesh$triMat)
  
  #-------------------------------------------------------------------------------
  
  #### DIFFERENTIAL OPERATOR 
  
  ## Expression of the function defining the diffusion directions <- Change here to get other diffusive behaviors
  # f_expr=expression(sin(x[1])*sin(x[2]))
  f_expr=expression(sin(theta)**2*(2*cos(theta)*cos(phi)+0*sin(phi)))
  
  
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
    
    # ## Gradient of f wrt to the metric (= G^{-1} d f / d x ), normalized
    # gv=c(v[1],(1/sin(x[1])**2)*v[2])*sin(x[1])/((sin(x[1])**2*v[1]**2+1*v[2]**2))**0.5
    # ## Orthogonal of the gradient, normalized
    # gv_orth=c(-v[2],v[1])*1.0/((sin(x[1])**2)*v[1]**2+v[2]**2)**0.5
    
    ## Compute the diffusion directions
    u=c(v[1],v[2])*sin(x[1])/((sin(x[1])**2*v[1]**2+1*v[2]**2))**0.5
    uorth=c(-v[2],(sin(x[1])**2)*v[1])*1.0/((sin(x[1])**2)*v[1]**2+v[2]**2)**0.5
    
    ## Weights associated to each direction
    rho=c(0.1+0.6/(1+exp(-cos(x[1])/0.25)),1)
    
    return(rho[1]*rho[2]*((1/rho[1])*u%*%t(u)+(1/rho[2])*uorth%*%t(uorth)))
  }
  
  ## Compute coefficients of the diffusion tensor on each triangle
  triDiffMat=triMetricSph(mesh$nodeMat, mesh$triMat, GfuncSph = GfuncSph)
  
  #-------------------------------------------------------------------------------
  
  ####  POTENTIAL FUNCTION
  potf<-function(x){
    500*(1+5*cos(pi*x[1])**2)
  }
  ## Potential function on each triangle
  triPotCoeffs=triPot(mesh$nodeMat,mesh$triMat,potf,coordSph = FALSE)
  
  #-------------------------------------------------------------------------------
  
  
  #-------------------------------------------------------------------------------
  
  ####  FEM MATRICES
  ## Build FEM matrices
  FEMatList=matFEM2d(mesh$nodeMat,mesh$triMat,triMetricMat,triPotCoeffs,triDiffMat,massLumping = FALSE)
  ## Print smallest and highest eigenvalues: if too big (>1e10), your diffusion may have singularities
  cat(c("Eigenvalue interval for Chebyshev approximation :", FEMatList$Eig))
  
  Cdh=Diagonal(mesh$N,1/diag(FEMatList$Scale))
  
  if(computeChol){
    Chol=expand(Cholesky(FEMatList$Mass))
  }else{
    Chol=NULL
  }
  
  h=compute_h(mesh$nodeMat,mesh$triMat)
  
  ### Add result to list
  new_mesh=mesh
  new_mesh[["FEMatList"]]=FEMatList
  new_mesh[["Cdh"]]=Cdh
  new_mesh[["Chol"]]=Chol
  new_mesh[["h"]]=h
  return(new_mesh)
}

#-------------------------------------------------------------------------------

#### CREATE FINE SPHERE MESH

# Import mesh
sphFine=importDatSphere(idFine)
sphFine=addFEM(sphFine,computeChol = useChol)

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
    # CREATE  COARSE MESH
    sphCoarse=importDatSphere(lvl)
    
    ## Compute projection matrix from precomputed quantities
    id_start=which(idFine_list==(lvl+1))
    id_end=which(idFine_list==(idFine))
    Amat=A_list[[id_start]]
    if((id_start+1)<= id_end){
      for (l in (id_start+1):id_end) {
        Amat=Amat%*%A_list[[l]]
      }
    }
    sphCoarse[["A"]]=Amat
    sphCoarse=addFEM(sphCoarse)
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
  }
}


#-------------------------------------------------------------------------------

## Save results
if(export){
  fileExp=paste0(expFolder,"results_nstat.RData")
  # save.image(file=fileExp)
  save(list=c("err_list","n_list","h_list","alpha_l","lvl_list"),file=fileExp)
  print(paste0("Results saved at location :", fileExp))
}
#-------------------------------------------------------------------------------

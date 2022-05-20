#' @title Distinct Lesion Parcellation using Center Detection and Mixture Modeling
#' @description This function finds distinct lesions based on a lesion probability map, and partitions voxels into individual lesion assignments.
#' @param probmap a 3D array or image of class \code{nifti}, containing the probability that each voxel is a lesion voxel
#' @param binmap a 3D array or \code{nifti} mask in which voxels are classified as either lesion voxels or not lesion voxels.
#' Note that mask should be in the same space as the probmap volume
#' @param smooth a number specifying the width (in voxels) of the desired smooth.
#' Note: For some probability maps with high confidence (many voxels at or near 1.0), increasing the smooth can be helpful.
#' @param minCenterSize an integer value representing the minimum number of connected voxels that can be considered a lesion center
#' @param parallel is a logical value that indicates whether the user's computer
#' is Linux or Unix (i.e. macOS), and should run the code in parallel
#' @param cores if parallel = TRUE, cores is an integer value that indicates how many cores
#' the function should be run on
#'
#' @importFrom ANTsRCore labelClusters
#' @importFrom extrantsr ants2oro oro2ants
#' @importFrom neurobase readnii writenii
#' @importFrom parallel mclapply
#' @importFrom pdist pdist
#' @importFrom plyr mapvalues
#' @importFrom mclust Mclust
#'
#' @return A list containing lesioncenters (a nifti file with labeled lesion centers),
#' lesionclusters.nn (a nifti file with labeled lesion clusters based on center detection and nearest neighbor assignment),
#' lesionclusters.gmm (a nifti file with labeled lesion clusters based on gaussian mixture modeling; if gmm=T),
#' and lesionclusters.cc (a nifti file with labeled lesion clusters based on naive connected components)
#' @examples \dontrun{
#' library(neurobase)
#' lesion.probs <- readnii('path/to/probabilitymap')
#' clusters <- lesionclusters(probmap = lesion.probs, binmap = lesion.probs>0.30,
#'                            parallel = TRUE, cores = 4) }
#' @export
lesionclusters=function(probmap,binmap,smooth=1.2,gmm=F,minCenterSize=10,
                        minLesionSize=1,parallel=F,cores=2,c3d_path=NULL){
  if(is.null(c3d_path)){
    if(file.exists("/Applications/ITK-SNAP.app/Contents/bin/c3d")){
      c3d_path="/Applications/ITK-SNAP.app/Contents/bin/c3d"
    }else{
      stop("Cannot find path to Convert3D\n
           If it is already installed via ITK-SNAP, please specify the correct path in the function call.\n
           If not, please download the software at http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.SNAP3.")
    }
  }


  ##################################
  ## Get connected components map ##
  ##################################

  scale=ceiling((1/mean(probmap@pixdim[2:4]))^3)
  clusmap=ants2oro(labelClusters(oro2ants(binmap),minClusterSize=minLesionSize*scale))

  #########################
  ## Do center detection ##
  #########################

  tempprob=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".nii.gz")
  writenii(probmap,tempprob)
  tempeigs=tempfile(pattern = "file", tmpdir = tempdir())
  system(paste0(c3d_path," ",tempprob," -smooth ", smooth, "vox -grad -oo ",paste0(tempeigs,"%02d.nii.gz")))

  tempx=tempfile(pattern = "file", tmpdir = tempdir())
  system(paste0(c3d_path," ",paste0(tempeigs,"00.nii.gz")," -grad -oo ",
                paste0(tempx,"%02d.nii.gz")))
  tempy=tempfile(pattern = "file", tmpdir = tempdir())
  system(paste0(c3d_path," ",paste0(tempeigs,"01.nii.gz")," -grad -oo ",
                paste0(tempy,"%02d.nii.gz")))
  tempz=tempfile(pattern = "file", tmpdir = tempdir())
  system(paste0(c3d_path," ",paste0(tempeigs,"02.nii.gz")," -grad -oo ",
                paste0(tempz,"%02d.nii.gz")))

  gxx=readnii(paste0(tempx,"00.nii.gz")); gxy=readnii(paste0(tempx,"01.nii.gz")); gxz=readnii(paste0(tempx,"02.nii.gz"))
  gyx=readnii(paste0(tempy,"00.nii.gz")); gyy=readnii(paste0(tempy,"01.nii.gz")); gyz=readnii(paste0(tempy,"02.nii.gz"))
  gzx=readnii(paste0(tempz,"00.nii.gz")); gzy=readnii(paste0(tempz,"01.nii.gz")); gzz=readnii(paste0(tempz,"02.nii.gz"))

  mask=binmap
  bigmat=cbind(as.vector(gxx[mask==1]),as.vector(gxy[mask==1]),as.vector(gxz[mask==1]),
               as.vector(gyx[mask==1]),as.vector(gyy[mask==1]),as.vector(gyz[mask==1]),
               as.vector(gzx[mask==1]),as.vector(gzy[mask==1]),as.vector(gzz[mask==1]))
  rm(gxx,gxy,gxz,gyx,gyy,gyz,gzx,gzy,gzz)

  biglist=split(bigmat,row(bigmat))
  biglist=lapply(biglist,matrix,nrow=3,byrow=T)
  rm(bigmat)

  getevals=function(matrix){
    thiseig=eigen(matrix)$values
    sort=order(abs(thiseig))
    return(thiseig[sort])
  }

  if(parallel==TRUE){
    result=matrix(unlist(mclapply(biglist,getevals,mc.cores=cores)),
                  ncol=3,byrow=T)
  }else if(parallel==FALSE){
    result=matrix(unlist(lapply(biglist,getevals)),ncol=3,byrow=T)
  }
  phes1=mask; phes1[mask==1]<-result[,1]
  phes2=mask; phes2[mask==1]<-result[,2]
  phes3=mask; phes3[mask==1]<-result[,3]

  centers=clusmap
  centers[centers!=0]<-1
  centers[phes1>0 | phes2>0 | phes3>0]<-0
  centers=ants2oro(labelClusters(oro2ants(centers),minClusterSize=minCenterSize*scale))
  centers[centers>0]=centers[centers>0]+max(as.vector(clusmap))


  ####################################
  ## Do nearest neighbor clustering ##
  ####################################

  uniqueclusts=unique(as.vector(clusmap))
  uniqueclusts=uniqueclusts[uniqueclusts!=0]
  nnOneLesion=function(x,clusmap,centers){
    assigndf=which(clusmap==x,arr.ind=T)
    centvals=unique(centers[assigndf])
    if(length(centvals)<=2){
      centvals[centvals==0]=x
      lesnum=max(centvals)
      assigndf=cbind(assigndf,lesnum)
    }else{
      assigned=which(clusmap==x & centers>0,arr.ind=T)
      unassigned=which(clusmap==x & centers==0,arr.ind=T)
      distmat=as.matrix(pdist(unassigned,assigned))
      closest=apply(distmat,1,which.min)
      unassigned=cbind(unassigned,centers[assigned[closest,]])
      assigned=cbind(assigned,centers[assigned])
      assigndf=rbind(assigned,unassigned)
    }
    return(assigndf)
  }
  if(parallel==TRUE){
    clusassignments=do.call(rbind,mclapply(uniqueclusts,nnOneLesion,
                                           clusmap,centers,mc.cores=cores))
  }else{
    clusassignments=do.call(rbind,lapply(uniqueclusts,nnOneLesion,clusmap,centers))
  }
  uniquelesions=clusassignments[,4] %>% as.vector() %>%
    table() %>% names() %>% as.numeric()
  clusassignments[,4]=mapvalues(clusassignments[,4],from=uniquelesions,
                                to=1:length(uniquelesions))
  nnmap=clusmap
  nnmap[clusassignments[,1:3]]=clusassignments[,4]

  centers[centers!=0]=mapvalues(centers[centers!=0],from=uniquelesions,
                                to=1:length(uniquelesions),
                                warn_missing=F)


  if(gmm==T){
    ####################################
    ## Do gaussian mixture clustering ##
    ####################################

    gmmOneLesion=function(x,clusmap,probmap,centers){
      coords=which(clusmap==x,arr.ind=T)
      assigndf=coords
      data=cbind(coords,probmap[coords]-min(probmap[coords]))
      data=as.data.frame(data)
      colnames(data)=c("x","y","z","p")
      centvals=unique(centers[coords])
      if(length(centvals)<=2){
        centvals[centvals==0]=x
        lesnum=max(centvals)
        assigndf=cbind(assigndf,lesnum)
      }else{
        nvox=nrow(coords); k=length(centvals)-1
        reprow=round(data$p*100,0)
        newdat=rep(1:nrow(data),reprow)
        #newdat=sample(1:nrow(data),size=nvox*10,replace=T,
        #              prob=(data$p)/sum(data$p))
        newdat=data[newdat,1:3]
        newdat$x=newdat$x+runif(nrow(newdat),-.5,.5)
        newdat$y=newdat$y+runif(nrow(newdat),-.5,.5)
        newdat$z=newdat$z+runif(nrow(newdat),-.5,.5)

        gmmauto <- NULL
        attempt <- 1
        while(is.null(gmmauto) && attempt <= 5) {
          attempt <- attempt + 1
          mink=max(1,round(k-.5*k,0))
          maxk=round(k+.5*k,0)
          gmmauto=try(
            Mclust(newdat,G=mink:maxk,nodelNames="VVV",verbose=F)
          )
        }

        mvn.pdf.i <- function(xi, mu, Sigma){
          1/sqrt( (2*pi)^length(xi) * det(Sigma) ) *
            exp(-(1/2) * t(xi - mu) %*% mvnorm.cov.inv(Sigma)
                %*% (xi - mu)  )
        }
        mvn.pdf <- function(X, mu, Sigma){
          apply(X, 1, function(xi) mvn.pdf.i(as.numeric(xi), mu, Sigma))
        }
        mvnorm.cov.inv <- function(Sigma) {
          # Eigendecomposition of covariance matrix
          E <- eigen(Sigma)
          Lambda.inv <- diag(E$values^-1)   # diagonal matrix
          Q <- E$vectors
          return(Q %*% Lambda.inv %*% t(Q))
        }

        mu=t(gmmauto$parameters$mean)
        cov=gmmauto$parameters$variance$sigma
        w=gmmauto$parameters$pro; k=length(w)
        mvn.c <- sapply(1:k, function(c) mvn.pdf(data[,1:3], mu[c,], cov[,, c]))
        gmm.softcluster <- t(w*t(mvn.c)) / rowSums(t(w*t(mvn.c)))
        gmm.cluster = apply(gmm.softcluster, 1, which.max)

        assigndf=cbind(assigndf,gmm.cluster)

      }
      return(assigndf)
    }
    if(parallel==TRUE){
      clusassignments=do.call(rbind,mclapply(uniqueclusts,gmmOneLesion,clusmap,
                                             probmap,centers,mc.cores=cores))
    }else{
      clusassignments=do.call(rbind,lapply(uniqueclusts,gmmOneLesion,clusmap,
                                           probmap,centers))
    }
    uniquelesions=clusassignments[,4] %>% as.vector() %>%
      table() %>% names() %>% as.numeric()
    clusassignments[,4]=mapvalues(clusassignments[,4],from=uniquelesions,
                                  to=1:length(uniquelesions))
    gmmmap=clusmap
    gmmmap[clusassignments[,1:3]]=clusassignments[,4]

    return(list(lesioncenters=centers, lesionclusters.nn=nnmap,
                lesionclusters.gmm=gmmmap, lesionclusters.cc=clusmap))
  }else{
    return(list(lesioncenters=centers, lesionclusters.nn=nnmap,
                lesionclusters.cc=clusmap))
  }

}

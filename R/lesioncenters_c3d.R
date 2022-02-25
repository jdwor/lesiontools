#' @title Distinct Lesion Centers using Convert3D gradients
#' @description This function finds the centers of distinct lesions based on a lesion probability map. The method is described in Dworkin et al., (2018).
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
#'
#' @return A list containing lesioncenters (a nifti file with labeled lesion centers) and lesioncount (an integer value representing the number of distinct lesions)
#' @examples \dontrun{
#' library(neurobase)
#' lesion.probs <- readnii('path/to/probabilitymap')
#' centers <- lesioncenters(probmap = lesion.probs, binmap = lesion.probs>0.30,
#'                          parallel = TRUE, cores = 4) }
#' @export
#' @references J.D. Dworkin, K.A. Linn, I. Oguz, G.M. Fleishman, R. Bakshi, G. Nair, P.A. Calabresi, R.G. Henry, J. Oh, N. Papinutto, D. Pelletier, W. Rooney, W. Stern, N.L. Sicotte, D.S. Reich, R.T. Shinohara. An automated statistical technique for counting distinct multiple sclerosis lesions. American Journal of Neuroradiology, 2018; 39, 626-633.
lesioncenters_c3d=function(probmap,binmap,smooth=1.2,minCenterSize=10,parallel=F,cores=2,c3d_path=NULL){
  if(is.null(c3d_path)){
    if(file.exists("/Applications/ITK-SNAP.app/Contents/bin/c3d")){
      c3d_path="/Applications/ITK-SNAP.app/Contents/bin/c3d"
    }else{
      stop("Cannot find path to Convert3D\n
           If it is already installed via ITK-SNAP, please specify the correct path in the function call.\n
           If not, please download the software at http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.SNAP3.")
    }
  }

  scale=ceiling((1/mean(probmap@pixdim[2:4]))^3)

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

  gxx=readnii(paste0(tempx,"00.nii.gz"))
  gxy=readnii(paste0(tempx,"01.nii.gz"))
  gxz=readnii(paste0(tempx,"02.nii.gz"))
  gyx=readnii(paste0(tempy,"00.nii.gz"))
  gyy=readnii(paste0(tempy,"01.nii.gz"))
  gyz=readnii(paste0(tempy,"02.nii.gz"))
  gzx=readnii(paste0(tempz,"00.nii.gz"))
  gzy=readnii(paste0(tempz,"01.nii.gz"))
  gzz=readnii(paste0(tempz,"02.nii.gz"))

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
  phes1=mask
  phes1[mask==1]<-result[,1]
  phes2=mask
  phes2[mask==1]<-result[,2]
  phes3=mask
  phes3[mask==1]<-result[,3]

  clusmap=ants2oro(labelClusters(oro2ants(binmap),minClusterSize=10*scale))

  les=clusmap
  les[les!=0]<-1
  les[phes1>0 | phes2>0 | phes3>0]<-0
  les=ants2oro(labelClusters(oro2ants(les),minClusterSize=minCenterSize*scale))

  return(list(lesioncenters=les,lesioncount=max(les)))
}

#' @title Distinct Lesion Centers
#' @description This function finds the centers of distinct lesions based on a lesion probability map. The method is described in Dworkin et al., (2018).
#' @param probmap a 3D array or image of class \code{nifti}, containing the probability that each voxel is a lesion voxel
#' @param binmap a 3D array or \code{nifti} mask in which voxels are classified as either lesion voxels or not lesion voxels.
#' Note that mask should be in the same space as the probmap volume
#' @param c3d a logical value reflecting whether or not the Convert3D imaging toolbox is installed.
#' @param minCenterSize an integer value representing the minimum number of connected voxels that can be considered a lesion center
#' @param radius an integer specifying radius of the neighborhood (in voxels) for which the hessian should be calculated.
#' @param parallel is a logical value that indicates whether the user's computer
#' is Linux or Unix (i.e. macOS), and should run the code in parallel
#' @param cores if parallel = TRUE, cores is an integer value that indicates how many cores
#' the function should be run on
#'
#' @importFrom ANTsRCore labelClusters
#' @importFrom extrantsr ants2oro oro2ants
#' @return A list containing lesioncenters (a nifti file with labeled lesion centers) and lesioncount (an integer value representing the number of distinct lesions)
#' @examples \dontrun{
#' library(neurobase)
#' lesion.probs <- readnii('path/to/probabilitymap')
#' centers <- lesioncenters(probmap = lesion.probs, binmap = lesion.probs>0.30,
#'                          parallel = TRUE, cores = 4) }
#' @export
#' @references J.D. Dworkin, K.A. Linn, I. Oguz, G.M. Fleishman, R. Bakshi, G. Nair, P.A. Calabresi, R.G. Henry, J. Oh, N. Papinutto, D. Pelletier, W. Rooney, W. Stern, N.L. Sicotte, D.S. Reich, R.T. Shinohara. An automated statistical technique for counting distinct multiple sclerosis lesions. American Journal of Neuroradiology, 2018; 39, 626-633.
lesioncenters=function(probmap,binmap,c3d=T,minCenterSize=10,radius=1,parallel=F,cores=2){
  scale=ceiling((1/mean(probmap@pixdim[2:4]))^3)
  if(c3d==T){
    if(file.exists("/Applications/ITK-SNAP.app/Contents/bin/c3d")){
      c3d_path="/Applications/ITK-SNAP.app/Contents/bin/c3d"
    }else{
      stop("Cannot find path to Convert3D\n
           If it is already installed via ITK-SNAP, please specify the correct path in the function call.\n
           If not, please download the software at http://www.itksnap.org/pmwiki/pmwiki.php?n=Downloads.SNAP3.")
    }
    tempprob=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".nii.gz")
    writenii(probmap,tempprob)
    tempeigs=tempfile(pattern = "file", tmpdir = tempdir())
    system(paste0(c3d_path," ",tempprob," -hesseig ",scale," -oo ",paste0(tempeigs,"%02d.nii.gz")))

    phes1=readnii(paste0(tempeigs,"00.nii.gz"))
    phes2=readnii(paste0(tempeigs,"01.nii.gz"))
    phes3=readnii(paste0(tempeigs,"02.nii.gz"))
  }else if(c3d==F){
    phes=hessian(probmap,mask=binmap,radius,parallel,cores)
    phes1=phes$eigval1
    phes2=phes$eigval2
    phes3=phes$eigval3
  }else{
    stop("'c3d' must be TRUE or FALSE")
  }
  clusmap=ants2oro(labelClusters(oro2ants(binmap),minClusterSize=20*scale))

  les=clusmap
  les[les!=0]<-1
  les[phes1>0 | phes2>0 | phes3>0]<-0
  les=ants2oro(labelClusters(oro2ants(les),minClusterSize=minCenterSize*scale))

  return(list(lesioncenters=les,lesioncount=max(les)))
}


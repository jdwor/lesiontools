#' @title Distinct Lesion Centers
#' @description This function finds the centers of distinct lesions based on a lesion probability map. The method is described in Dworkin et al., (2018).
#' @param probmap a 3D array or image of class \code{nifti}, containing the probability that each voxel is a lesion voxel
#' @param binmap a 3D array or \code{nifti} mask in which voxels are classified as either lesion voxels or not lesion voxels.
#' Note that mask should be in the same space as the probmap volume
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
lesioncenters=function(probmap,binmap,minCenterSize=10,radius=1,parallel=F,cores=2){
  scale=ceiling((1/mean(probmap@pixdim[2:4]))^3)
  phes=hessian(probmap,mask=binmap,radius,parallel,cores)
  clusmap=ants2oro(labelClusters(oro2ants(binmap),minClusterSize=20*scale))

  les=clusmap
  les[les!=0]<-1
  les[phes$eigval1>0 | phes$eigval2>0 | phes$eigval3>0]<-0
  les=ants2oro(labelClusters(oro2ants(les),minClusterSize=minCenterSize*scale))

  return(list(lesioncenters=les,lesioncount=max(les)))
}

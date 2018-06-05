#' @title Distance to Mask Boundary
#' @description This function finds the distance of each voxel to the nearest boundary in a given mask.
#' @param mask a 3D array or image of class \code{nifti}, containing a binary mask where 1 represents structure.
#'
#' @return A new image in which voxels have been assigned their distance to the nearest boundary.
#' @examples \dontrun{
#' library(neurobase)
#' lesion.mask <- readnii('path/to/mask')
#' dtb <- dtboundary(mask = lesion.mask) }
#' @export
dtboundary<-function(mask) {
  get.d.boundary.exact.balloon<-function(v,mask,d.max=30){
    if (mask[v[1],v[2],v[3]]==0) print('ERROR! - voxel outside of mask...')
    inf<-1000
    balloon.empty<-TRUE
    r<-1
    #expand balloon
    while(balloon.empty){
      balloon<-1-mask[(v[1]-r):(v[1]+r),(v[2]-r):(v[2]+r),(v[3]-r):(v[3]+r)]
      #If balloon had reached edge
      if (sum(balloon>0)){
        which.outside<-which(balloon>0,arr.ind=TRUE)
        d.out<-min(sqrt((which.outside[,1]-(r+1))^2+(which.outside[,2]-(r+1))^2+(which.outside[,3]-(r+1))^2))
        balloon.empty<-FALSE
      }
      else {
        if (r<=d.max) {
          r<-r+1
        }
        else {
          balloon.empty<-FALSE
          d.out<-inf
        }
      }
    }
    return(d.out)
  }
  which.mask.arrind<-which(mask>0,arr.ind=TRUE)
  #For each voxel in the mask
  min.d<-rep(0,dim(which.mask.arrind)[1])
  for (i in 1:(dim(which.mask.arrind)[1])) {
    #Get minimum distance to boundary
    min.d[i]<-get.d.boundary.exact.balloon(which.mask.arrind[i,],mask)
  }
  mask[mask>0]<-min.d
  return(mask)
}

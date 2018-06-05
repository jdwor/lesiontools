#' @title Register Label Image Based on Full Image
#' @description This function registers a full image to a fixed image, then applies the registration to a label or binary image in the same space as the full image.
#' @param fullimage an image of class \code{\link{nifti}}, in the same space as the label image, which will be registered to the fixed image.
#' @param labelimage an image of class \code{\link{nifti}} with limited structural information, in the same space as the full image, to which the full image registration will be applied.
#' @param fixedimage an image of class \code{\link{nifti}}, to which the other images will be registered.
#' @param typeofTransform the type of registration desired; this value is passed onto the registration function from extrantsr.
#' @param interpolator the type of interpolation desired; this value is passed onto the antsApplyTransforms function from ANTsRCore
#'
#' @importFrom ANTsRCore antsApplyTransforms
#' @importFRom extrantsr ants2oro oro2ants registration
#' @return A list containing image_reg (the registered version of fullimage) and label_reg (the registered version of labelimage).
#' library(neurobase)
#' library(fslr)
#' flair <- readnii('path/to/flair')
#' t1 <- readnii('path/to/t1')
#' tissue.class= <- fast(t1, opts='--nobias')
#' registered <- labelreg(t1,tissue.class,flair)
#' t1_reg = registered$image_reg
#' tissue.class_reg = registered$label_reg }
#' @export
labelreg=function(fullimage,labelimage,fixedimage,typeofTransform="Rigid",interpolator="lanczosWindowedSinc"){
  imtofix=registration(filename=fullimage,template.file=fixedimage,
                       typeofTransform=typeofTransform,remove.warp=FALSE,
                       outprefix="fun")
  labtofix<-antsApplyTransforms(fixed=oro2ants(fixedimage),moving=oro2ants(labelimage),
                                transformlist=imtofix$fwdtransforms,
                                interpolator=interpolator)
  return(list(image_reg=ants2oro(imtofix$outfile),label_reg=ants2oro(labtofix)))
}

#' @title Distinct Lesion Centers
#' @description This function finds the centers of distinct lesions based on a lesion probability map. The method is described in Dworkin (2018).
#' @param epi a T2*-EPI volume of class \code{\link{nifti}}.
#' @param t1 a T1-weighted volume of class \code{\link{nifti}}.
#' @param flair a T2-FLAIR volume of class \code{\link{nifti}}.
#' @param probmap an image of class \code{\link{nifti}}, containing the probability that each voxel
#' is a lesion voxel.
#' If a probability map is not included, the MIMoSA model will be applied (Valcarcel, 2017).
#' @param binmap a \code{\link{nifti}} mask in which voxels are classified as either lesion voxels
#' or not lesion voxels.
#' Note that mask should be in the same space as the probmap volume.
#' @param parallel is a logical value that indicates whether the user's computer
#' is Linux or Unix (i.e. macOS), and should run the code in parallel.
#' @param cores if parallel = TRUE, cores is an integer value that indicates how many cores
#' the function should be run on.
#' @param skullstrip a logical value reflecting whether or not the images have already been skull-stripped.
#' @param biascorrect a logical value reflecting whether or not the images have already been bias-corrected.
#' @param c3d a logical value reflecting whether or not the Convert3D imaging toolbox is installed.
#'
#' @importFrom ANTsRCore labelClusters
#' @importFrom neurobase niftiarr
#' @importFrom mimosa mimosa_data mimosa_model_No_PD_T2
#' @importFrom extrantsr ants2oro oro2ants bias_correct registration
#' @importFrom stats predict
#' @importFrom fslr fslsmooth fast fslerode
#' @return A list containing candidate.lesions (a nifti file with labeled lesions evaluated for CVS),
#' cvs.probmap (a nifti file in which candidate lesions are labeled with their CVS probability), and
#' cvs.biomarker (a numeric value representing the average CVS probability of a subject's lesions).
#' @examples \dontrun{
#' library(neurobase)
#' epi <- readnii('path/to/epi')
#' flair <- readnii('path/to/flair')
#' t1 <- readnii('path/to/t1')
#' cvs <- centralveins(epi = epi, t1 = t1, flair = flair
#'                          parallel = TRUE, cores = 4, c3d = T) }
#' @export
#' @references Alessandra Valcarcel (2017). mimosa: 'MIMoSA': A Method for Inter-Modal Segmentation Analysis. R package version 0.5.3. https://github.com/avalcarcel9/mimosa
centralveins=function(epi,t1,flair,probmap=NULL,binmap=NULL,parallel=F,
                      cores=2,skullstrip=F,biascorrect=F,c3d=F){
  if(biascorrect==F){
    epi=bias_correct(epi,correction="N4",reorient=F)
    t1=bias_correct(t1,correction="N4",reorient=F)
    flair=bias_correct(flair,correction="N4",reorient=F)
  }
  flair=registration(filename=flair,template.file=t1,typeofTransform="Rigid",
                     remove.warp=FALSE,outprefix="fun")
  flair=flair$outfile

  if(skullstrip==F){
    t1_ss=fslbet_robust(t1,correct=F)
    epi_ss=fslbet_robust(epi,correct=F)
    flair_ss=flair
    flair_ss[t1_ss==0]<-0
  }

  frangi=frangi(image=epi_ss,mask=epi_ss!=0,parallel=parallel,cores=cores,c3d)
  frangi[frangi<0]<-0

  regs=labelreg(epi,t1,frangi)
  epi_t1=regs$imagereg
  frangi_t1=regs$labelreg
  mask=t1!=0

  if(is.null(probmap) & is.null(binmap)){
    mimosa_data = mimosa_data(
      brain_mask = mask,
      FLAIR = flair_ss,
      T1 = t1_ss,
      normalize = 'Z',
      cores = cores,
      verbose = TRUE)

    mimosa_df = mimosa_data$mimosa_dataframe
    mimosa_cm = mimosa_data$top_voxels
    rm(mimosa_data)

    predictions = predict(mimosa_model_No_PD_T2,newdata = mimosa_df,type = 'response')
    probmap = niftiarr(brain_mask, 0)
    probmap[mimosa_cm == 1] = predictions
    probmap = fslsmooth(probability_map,sigma = 1.25,
                        mask = mask,retimg = TRUE,smooth_mask = TRUE)
  }

  les=lesioncenters(probmap,probmap>0.3,parallel=parallel,cores=cores)

  csf=fast(t1,opts='--nobias')
  csf[csf!=1]<-0
  csf=ants2oro(labelClusters(oro2ants(csf),minClusterSize=300))
  csf[csf>0]<-1
  csf=(csf!=1)
  csf=fslerode(csf, kopts = paste("-kernel boxv",3), verbose = TRUE)
  csf=(csf==0)

  lables=ants2oro(labelClusters(oro2ants(les),minClusterSize=27))
  for(j in lables){
    if(sum(csf[lables==j])>0){
      lables[lables==j]<-0
    }
  }
  les=lables>0

  dtb=dtboundary(les)

  lables=ants2oro(labelClusters(oro2ants(les),minClusterSize=27))
  probles=lables
  avprob=NULL
  maxles=max(as.vector(lables))
  for(j in 1:maxles){
    frangsub=frangi[lables==j]
    centsub=dtbmap[lables==j]
    coords=which(lables==j,arr.ind=T)
    prod=frangsub*centsub
    score=sum(prod)
    nullscores=NULL
    for(k in 1:1000){
      samp=sample(1:length(centsub))
      centsamp=centsub[samp]
      coordsamp=coords[samp,]
      sampprod=frangsub*centsamp
      sampscore=sum(sampprod)
      nullscores=c(nullscores,sampscore)
    }
    lesprob=sum(nullscores<score)/length(nullscores)
    avprob=c(avprob,lesprob)
    probles[lables==j]<-lesprob

    print(paste0("Done with ",i," lesion ",j," of ",maxles))
  }

  return(list(candidate.lesions=lables,cvs.probmap=probles,cvs.biomarker=mean(avprob)))
}

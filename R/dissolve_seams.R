#very similar code for these two --> should probably make less redundant

#' @export
dissolve_spilhaus_seam<-function(x,
                                 intersect_tol=0.01,
                                 snap_tol=0.05,
                                 buffer_tol=0.005,
                                 seam_res=1000,
                                 revalidate_geoms=TRUE){

  #for sorting back into proper order
  vals<-terra::values(x)
  if(!length(vals)){
    vals<-data.frame("TEMPORARY_ORDER"=seq_along(x))
  }else{
    vals[["TEMPORARY_ORDER"]]<-seq_along(x)
  }
  terra::values(x)<-vals

  #map one side of spilhaus border onto lonlat projection
  lim<-11825474
  lim.seq<-seq(-lim,lim,length.out=seam_res)
  seam.xy<-rbind(cbind(-lim,lim.seq),
                 cbind(lim.seq,lim))
  seam.lonlat<-rSpilhaus:::.spilhaus2lonlat(seam.xy[,1],seam.xy[,2])
  #fatten up to make sure intersection test gets everything
  seam<-terra::buffer(terra::vect(seam.lonlat,"lines"),intersect_tol)
  #find geometries overlapping with seam
  inds<-which(terra::is.related(x,seam,"intersects"))

  #for each problem geometry...
  #first diaggregate into individual pieces
  #then snap their borders together (within some tolerance)
  #then buffer the shapes slightly
  #(not ideal, but found this little bit of fudging necessary to cover all cases)
  #then finally reaggregate into a single, dissolving any overlapping borders
  #9/28-->calling aggregate without "by" arguments results in dropping field attributes
  #now fixed by manually re-adding fields to aggregated shapes
  for(i in inds){
    tmp<-terra::disagg(x[i])
    if(snap_tol>0){
      tmp<-terra::snap(tmp,tol=snap_tol)
    }
    tmp<-terra::aggregate(terra::buffer(tmp,width=buffer_tol))
    terra::values(tmp)<-terra::values(x)[i,,drop=FALSE]
    x<-rbind(x,tmp)
  }
  #eliminate the old geometry elements (R didn't seem to like directly replacing the indices for whatever reason)
  x<-x[-inds]

  x<-terra::sort(x,v="TEMPORARY_ORDER")
  vals<-terra::values(x)
  vals[["TEMPORARY_ORDER"]]<-NULL
  terra::values(x)<-vals

  if(revalidate_geoms){
    x<-terra::buffer(x,width=0)
  }

  #found that cropping resulted in more trouble than it was worth...
  # #crop to appropriate extent
  # rectangle<-cbind(c(-180,
  #                    -180,
  #                    180,
  #                    180),
  #                  c(-90,
  #                    90,
  #                    90,
  #                    -90))
  # terra::crop(x,terra::vect(rectangle,"polygons"))
  x
}

#buffer_tol of 500 may be overkill for many cases, but it was necessary for
#rnaturalearth land data at 50 m I needed to use 500 for seam around Bering strait,
#and 1200 for seam in Antarctica!
#is there a way to specifically target "cracks"? I think they're not snapping because they're not full borders...
#' @export
dissolve_IDL_seam<-function(x,
                            intersect_tol=1000,
                            snap_tol=5000,
                            buffer_tol=500,
                            seam_res=1000,
                            revalidate_geoms=TRUE){

  #for sorting back into proper order
  vals<-terra::values(x)
  if(!length(vals)){
    vals<-data.frame("TEMPORARY_ORDER"=seq_along(x))
  }else{
    vals[["TEMPORARY_ORDER"]]<-seq_along(x)
  }
  terra::values(x)<-vals

  #map -180/180 longitude line onto spilhaus projection
  seam.lonlat<-cbind(180,seq(-90,90,length.out=seam_res))
  seam.xy<-.lonlat2spilhaus(seam.lonlat[,1],seam.lonlat[,2])
  #split out separate components (thankfully one half always has positive y coords, the other always negative)
  seam.xy<-lapply(split(seam.xy,seam.xy[,2]>0),matrix,ncol=2)
  #fatten up to make sure intersection test gets everything
  seam<-terra::buffer(terra::vect(seam.xy,"lines"),intersect_tol)
  #find geometries overlapping with seam
  inds<-which(terra::is.related(x,seam,"intersects"))

  #for each problem geometry...
  #first diaggregate into individual pieces
  #then snap their borders together (within some tolerance)
  #then buffer the shapes slightly
  #(not ideal, but found this little bit of fudging necessary to cover all cases)
  #then finally reaggregate into a single, dissolving any overlapping borders
  #9/28-->calling aggregate without "by" arguments results in dropping field attributes
  #now fixed by manually re-adding fields to aggregated shapes
  for(i in inds){
    tmp<-terra::disagg(x[i])
    if(snap_tol>0){
      tmp<-terra::snap(tmp,tol=snap_tol)
    }
    tmp<-terra::aggregate(terra::buffer(tmp,width=buffer_tol))
    terra::values(tmp)<-terra::values(x)[i,,drop=FALSE]
    x<-rbind(x,tmp)
  }
  #eliminate the old geometry elements (R didn't seem to like directly replacing the indices for whatever reason)
  x<-x[-inds]

  x<-terra::sort(x,v="TEMPORARY_ORDER")
  vals<-terra::values(x)
  vals[["TEMPORARY_ORDER"]]<-NULL
  terra::values(x)<-vals

  if(revalidate_geoms){
    x<-terra::buffer(x,width=0)
  }

  #found that cropping resulted in more trouble than it was worth...
  # #crop to appropriate extent
  # lim<-11825474
  # rectangle<-cbind(c(-lim,
  #                    -lim,
  #                    lim,
  #                    lim),
  #                  c(-lim,
  #                    lim,
  #                    lim,
  #                    -lim))
  # terra::crop(x,terra::vect(rectangle,"polygons"))
  x
}

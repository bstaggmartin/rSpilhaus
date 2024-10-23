#' @export
patch_corners<-function(x,
                        NA_tol=1.1e6+15,
                        SA_tol=1.1e6+15,
                        Asia_tol=1.1e6+15,
                        revalidate_geoms=TRUE){
  lim<-11825474

  #for sorting back into proper order
  vals<-terra::values(x)
  if(!length(vals)){
    vals<-data.frame("TEMPORARY_ORDER"=seq_along(x))
  }else{
    vals[["TEMPORARY_ORDER"]]<-seq_along(x)
  }
  terra::values(x)<-vals

  if(NA_tol>0){
    NA.patch1<-terra::vect(
      cbind(c(-lim,seq(-lim,-lim+NA_tol,length.out=100)),
            c(lim,seq(lim-NA_tol,lim,length.out=100))),
      "polygons")
    inds<-which(terra::is.related(x,NA.patch1,"intersects"))
    if(length(inds)>0){
      for(i in inds){
        tmp<-terra::aggregate(terra::union(x[i],NA.patch1))
        terra::values(tmp)<-terra::values(x)[i,,drop=FALSE]
        x<-rbind(x,tmp)
      }
      x<-x[-inds]
    }
    NA.patch2<-terra::vect(
      cbind(c(lim,seq(lim,lim-NA_tol,length.out=100)),
            c(-lim,seq(-lim+NA_tol,-lim,length.out=100))),
      "polygons")
    inds<-which(terra::is.related(x,NA.patch2,"intersects"))
    if(length(inds)>0){
      for(i in inds){
        tmp<-terra::aggregate(terra::union(x[i],NA.patch2))
        terra::values(tmp)<-terra::values(x)[i,,drop=FALSE]
        x<-rbind(x,tmp)
      }
      x<-x[-inds]
    }
  }
  if(SA_tol>0){
    SA.patch<-
      terra::vect(
        cbind(c(-lim,seq(-lim,-lim+SA_tol,length.out=100)),
              c(-lim,seq(-lim+SA_tol,-lim,length.out=100))),
        "polygons")
    inds<-which(terra::is.related(x,SA.patch,"intersects"))
    if(length(inds)>0){
      for(i in inds){
        tmp<-terra::aggregate(terra::union(x[i],SA.patch))
        terra::values(tmp)<-terra::values(x)[i,,drop=FALSE]
        x<-rbind(x,tmp)
      }
      x<-x[-inds]
    }
  }
  if(Asia_tol>0){
    Asia.patch<-
      terra::vect(
        cbind(c(lim,seq(lim,lim-Asia_tol,length.out=100)),
              c(lim,seq(lim-Asia_tol,lim,length.out=100))),
        "polygons")
    inds<-which(terra::is.related(x,Asia.patch,"intersects"))
    if(length(inds)>0){
      for(i in inds){
        tmp<-terra::aggregate(terra::union(x[i],Asia.patch))
        terra::values(tmp)<-terra::values(x)[i,,drop=FALSE]
        x<-rbind(x,tmp)
      }
      x<-x[-inds]
    }
  }

  x<-terra::sort(x,v="TEMPORARY_ORDER")
  vals<-terra::values(x)
  vals[["TEMPORARY_ORDER"]]<-NULL
  terra::values(x)<-vals

  if(revalidate_geoms){
    terra::buffer(x,width=0)
  }else{
    x
  }
}

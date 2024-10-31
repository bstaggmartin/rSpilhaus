#very similar code for these two --> should probably make less redundant

#' @export
dissolve_spilhaus_seam<-function(x,
                                 seam_res=1000,
                                 seam_width=0.01,
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
  #extend to prevent border artifacts

  seam.lonlat<-lapply(split(seam.lonlat,seam.lonlat[,1]>0),matrix,ncol=2)
  min1.ind<-which.min(seam.lonlat[[1]][,1])
  min1<-seam.lonlat[[1]][min1.ind,]
  max2.ind<-which.max(seam.lonlat[[2]][,1])
  max2<-seam.lonlat[[2]][max2.ind,]
  if(min1.ind==1){
    seam.lonlat[[1]]<-rbind(c(max2[1]-360,max2[2]),
                            seam.lonlat[[1]])
  }else{
    seam.lonlat[[1]]<-rbind(seam.lonlat[[1]],
                            c(max2[1]-360,max2[2]))
  }
  if(max2.ind==1){
    seam.lonlat[[2]]<-rbind(c(min1[1]+360,min1[2]),
                            seam.lonlat[[2]])
  }else{
    seam.lonlat[[2]]<-rbind(seam.lonlat[[2]],
                            c(min1[1]+360,min1[2]))
  }

  #fatten up to make sure intersection test gets everything
  seam<-terra::buffer(terra::aggregate(terra::vect(seam.lonlat,"lines")),
                      seam_width,
                      capstyle="square")
  #find geometries overlapping with seam
  inds<-which(terra::is.related(x,seam,"intersects"))

  #new approach with convex hulls --> should be much cleaner...
  for(i in inds){
    uni<-terra::union(x[i],seam)
    tmp<-terra::disagg(uni[3])
    uni<-terra::aggregate(uni)

    # tmp.n<-rep(length(tmp),2)
    #need some step for remerging disaggregated bits by distance...
    #get distances and find the small ones...
    # dists<-dist(terra::geom(terra::centroids(tmp))[,c("x","y")])
    # tmp.grp<-which(.row(tmp.n)>.col(tmp.n),arr.ind=TRUE)[dists<2*seam_width,]
    #now have to set up to temporary grp variable...
    #tricky, but def doable...
    #I think this should always work?

    #better approach --> centroids of bounding boxes...
    if(length(tmp)>1){
      tmpy<-terra::geom(tmp)
      cents<-cbind((tapply(tmpy[,"x"],tmpy[,"geom"],min)+
                      tapply(tmpy[,"x"],tmpy[,"geom"],max))/2,
                   (tapply(tmpy[,"y"],tmpy[,"geom"],min)+
                      tapply(tmpy[,"y"],tmpy[,"geom"],max))/2)

      terra::values(tmp)<-
        data.frame(
          "TEMPORARY_GRP"=cutree(hclust(dist(cents)),h=2*seam_width)
        )
    }else{
      terra::values(tmp)<-data.frame("TEMPORARY_GRP"=1)
    }
    tmp<-terra::aggregate(terra::convHull(tmp,by="TEMPORARY_GRP"))

    #alpha hulls seem finicky...
    # tmp<-terra::aggregate(tmp,by="TEMPORARY_GRP")
    # terra::plot(terra::convHull(tmp[18]))
    # tmpy<-terra::geom(tmp[18])[,c("x","y")]
    # tmpy<-tmpy[!duplicated(tmpy),]
    # debug(alphahull::ashape)
    # alphahull::ashape(tmpy[,1],tmpy[,2],0.1)

    #I think if I take the difference between the hulls and the union of x and the seams...
    #in other words, just crop convex hulls with entire union of geom and seams!
    tmp<-terra::crop(tmp,uni)

    tmp<-terra::aggregate(terra::union(x[i],tmp))
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

#had an idea for new approach!
#find union of geoms with seams, split out intersection, than find convex hull of intersection and remerge

#' @export
dissolve_IDL_seam<-function(x,
                            seam_res=1000,
                            seam_width=1000,
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

  #map -180/180 longitude line onto spilhaus projection
  seam.lonlat<-cbind(180,seq(-90,90,length.out=seam_res))
  seam.xy<-rSpilhaus:::.lonlat2spilhaus(seam.lonlat[,1],seam.lonlat[,2])
  #split out separate components (thankfully one half always has positive y coords, the other always negative)
  seam.xy<-lapply(split(seam.xy,seam.xy[,2]>0),matrix,ncol=2)
  #extend to prevent border artifacts

  max1.ind<-which.max(seam.xy[[1]][,1])
  max1<-seam.xy[[1]][max1.ind,]
  max2.ind<-which.max(seam.xy[[2]][,2])
  max2<-seam.xy[[2]][max2.ind,]
  if(max1.ind==1){
    seam.xy[[1]]<-
      rbind(c(cos(pi/2)*(max2[1]-lim)-
                sin(pi/2)*(max2[2]-lim)+
                lim,
              sin(pi/2)*(max2[1]-lim)+
                cos(pi/2)*(max2[2]-lim)+
                lim),
            seam.xy[[1]])
  }else{
    seam.xy[[1]]<-
      rbind(seam.xy[[1]],
            c(cos(pi/2)*(max2[1]-lim)-
                sin(pi/2)*(max2[2]-lim)+
                lim,
              sin(pi/2)*(max2[1]-lim)+
                cos(pi/2)*(max2[2]-lim)+
                lim))
  }
  if(max2.ind==1){
    seam.xy[[2]]<-
      rbind(c(cos(-pi/2)*(max1[1]-lim)-
                sin(-pi/2)*(max1[2]-lim)+
                lim,
              sin(-pi/2)*(max1[1]-lim)+
                cos(-pi/2)*(max1[2]-lim)+
                lim),
            seam.xy[[2]])
  }else{
    seam.xy[[2]]<-
      rbind(seam.xy[[2]],
            c(cos(-pi/2)*(max1[1]-lim)-
                sin(-pi/2)*(max1[2]-lim)+
                lim,
              sin(-pi/2)*(max1[1]-lim)+
                cos(-pi/2)*(max1[2]-lim)+
                lim))
  }

  #fatten up to make sure intersection test gets everything
  seam<-terra::buffer(terra::aggregate(terra::vect(seam.xy,"lines")),
                      seam_width,
                      capstyle="square")
  #find geometries overlapping with seam
  inds<-which(terra::is.related(x,seam,"intersects"))

  #new approach with convex hulls --> should be much cleaner...
  for(i in inds){
    uni<-terra::union(x[i],seam)
    tmp<-terra::disagg(uni[3])
    uni<-terra::aggregate(uni)

    # tmp.n<-rep(length(tmp),2)
    #need some step for remerging disaggregated bits by distance...
    #get distances and find the small ones...
    # dists<-dist(terra::geom(terra::centroids(tmp))[,c("x","y")])
    # tmp.grp<-which(.row(tmp.n)>.col(tmp.n),arr.ind=TRUE)[dists<2*seam_width,]
    #now have to set up to temporary grp variable...
    #tricky, but def doable...
    #I think this should always work?

    #better approach --> centroids of bounding boxes...
    if(length(tmp)>1){
      tmpy<-terra::geom(tmp)
      cents<-cbind((tapply(tmpy[,"x"],tmpy[,"geom"],min)+
                      tapply(tmpy[,"x"],tmpy[,"geom"],max))/2,
                   (tapply(tmpy[,"y"],tmpy[,"geom"],min)+
                      tapply(tmpy[,"y"],tmpy[,"geom"],max))/2)

      terra::values(tmp)<-
        data.frame(
          "TEMPORARY_GRP"=cutree(hclust(dist(cents)),h=2*seam_width)
        )
    }else{
      terra::values(tmp)<-data.frame("TEMPORARY_GRP"=1)
    }
    tmp<-terra::aggregate(terra::convHull(tmp,by="TEMPORARY_GRP"))

    #I think if I take the difference between the hulls and the union of x and the seams...
    #in other words, just crop convex hulls with entire union of geom and seams!
    tmp<-terra::crop(tmp,uni)

    tmp<-terra::aggregate(terra::union(x[i],tmp))
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

# #old version
# dissolve_IDL_seam<-function(x,
#                             intersect_tol=1000,
#                             snap_tol=5000,
#                             buffer_tol=500,
#                             seam_res=1000,
#                             revalidate_geoms=TRUE){
#
#   #for sorting back into proper order
#   vals<-terra::values(x)
#   if(!length(vals)){
#     vals<-data.frame("TEMPORARY_ORDER"=seq_along(x))
#   }else{
#     vals[["TEMPORARY_ORDER"]]<-seq_along(x)
#   }
#   terra::values(x)<-vals
#
#   #map -180/180 longitude line onto spilhaus projection
#   seam.lonlat<-cbind(180,seq(-90,90,length.out=seam_res))
#   seam.xy<-.lonlat2spilhaus(seam.lonlat[,1],seam.lonlat[,2])
#   #split out separate components (thankfully one half always has positive y coords, the other always negative)
#   seam.xy<-lapply(split(seam.xy,seam.xy[,2]>0),matrix,ncol=2)
#   #fatten up to make sure intersection test gets everything
#   seam<-terra::buffer(terra::vect(seam.xy,"lines"),intersect_tol)
#   #find geometries overlapping with seam
#   inds<-which(terra::is.related(x,seam,"intersects"))
#
#   #for each problem geometry...
#   #first diaggregate into individual pieces
#   #then snap their borders together (within some tolerance)
#   #then buffer the shapes slightly
#   #(not ideal, but found this little bit of fudging necessary to cover all cases)
#   #then finally reaggregate into a single, dissolving any overlapping borders
#   #9/28-->calling aggregate without "by" arguments results in dropping field attributes
#   #now fixed by manually re-adding fields to aggregated shapes
#   for(i in inds){
#     tmp<-terra::disagg(x[i])
#     if(snap_tol>0){
#       tmp<-terra::snap(tmp,tol=snap_tol)
#     }
#     tmp<-terra::aggregate(terra::buffer(tmp,width=buffer_tol))
#     terra::values(tmp)<-terra::values(x)[i,,drop=FALSE]
#     x<-rbind(x,tmp)
#   }
#   #eliminate the old geometry elements (R didn't seem to like directly replacing the indices for whatever reason)
#   x<-x[-inds]
#
#   x<-terra::sort(x,v="TEMPORARY_ORDER")
#   vals<-terra::values(x)
#   vals[["TEMPORARY_ORDER"]]<-NULL
#   terra::values(x)<-vals
#
#   if(revalidate_geoms){
#     x<-terra::buffer(x,width=0)
#   }
#
#   #found that cropping resulted in more trouble than it was worth...
#   # #crop to appropriate extent
#   # lim<-11825474
#   # rectangle<-cbind(c(-lim,
#   #                    -lim,
#   #                    lim,
#   #                    lim),
#   #                  c(-lim,
#   #                    lim,
#   #                    lim,
#   #                    -lim))
#   # terra::crop(x,terra::vect(rectangle,"polygons"))
#   x
# }

# #old version
# dissolve_spilhaus_seam<-function(x,
#                                  intersect_tol=0.01,
#                                  snap_tol=0.05,
#                                  buffer_tol=0.005,
#                                  seam_res=1000,
#                                  revalidate_geoms=TRUE){
#
#   #for sorting back into proper order
#   vals<-terra::values(x)
#   if(!length(vals)){
#     vals<-data.frame("TEMPORARY_ORDER"=seq_along(x))
#   }else{
#     vals[["TEMPORARY_ORDER"]]<-seq_along(x)
#   }
#   terra::values(x)<-vals
#
#   #map one side of spilhaus border onto lonlat projection
#   lim<-11825474
#   lim.seq<-seq(-lim,lim,length.out=seam_res)
#   seam.xy<-rbind(cbind(-lim,lim.seq),
#                  cbind(lim.seq,lim))
#   seam.lonlat<-rSpilhaus:::.spilhaus2lonlat(seam.xy[,1],seam.xy[,2])
#   #fatten up to make sure intersection test gets everything
#   seam<-terra::buffer(terra::vect(seam.lonlat,"lines"),intersect_tol)
#   #find geometries overlapping with seam
#   inds<-which(terra::is.related(x,seam,"intersects"))
#
#   #for each problem geometry...
#   #first diaggregate into individual pieces
#   #then snap their borders together (within some tolerance)
#   #then buffer the shapes slightly
#   #(not ideal, but found this little bit of fudging necessary to cover all cases)
#   #then finally reaggregate into a single, dissolving any overlapping borders
#   #9/28-->calling aggregate without "by" arguments results in dropping field attributes
#   #now fixed by manually re-adding fields to aggregated shapes
#   for(i in inds){
#     tmp<-terra::disagg(x[i])
#     if(snap_tol>0){
#       tmp<-terra::snap(tmp,tol=snap_tol)
#     }
#     tmp<-terra::aggregate(terra::buffer(tmp,width=buffer_tol))
#     terra::values(tmp)<-terra::values(x)[i,,drop=FALSE]
#     x<-rbind(x,tmp)
#   }
#   #eliminate the old geometry elements (R didn't seem to like directly replacing the indices for whatever reason)
#   x<-x[-inds]
#
#   x<-terra::sort(x,v="TEMPORARY_ORDER")
#   vals<-terra::values(x)
#   vals[["TEMPORARY_ORDER"]]<-NULL
#   terra::values(x)<-vals
#
#   if(revalidate_geoms){
#     x<-terra::buffer(x,width=0)
#   }
#
#   #found that cropping resulted in more trouble than it was worth...
#   # #crop to appropriate extent
#   # rectangle<-cbind(c(-180,
#   #                    -180,
#   #                    180,
#   #                    180),
#   #                  c(-90,
#   #                    90,
#   #                    90,
#   #                    -90))
#   # terra::crop(x,terra::vect(rectangle,"polygons"))
#   x
# }

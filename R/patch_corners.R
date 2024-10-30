#2 key future improvements:
# - use tiling info for more accurate border extrapolation
#    - i.e., see what line SHOULD meet up with across border
# - use angle wrt origin so coordinates are distributed more evenly in corners
#adapting to prettified borders:
# - need to account for many, many more lines...tricky but not impossible
# - after thinking on it, probably best to use angle wrt origin to extrapolate
#   borders...otherwise they'll probably self-intersect and make a mess

.magic.corner.patcher<-function(x,
                                top,
                                left,
                                tol){
  vv<-hh<-11825474
  if(!top) hh<- -hh
  if(left) vv<- -vv

  #get coords
  coords<-x[,c("x","y")]
  nn<-nrow(coords)

  #find coords lying on focal line
  #tl --> y-(hh-tol)=x-vv --> y-hh+tol-x+vv
  #tr --> y-(hh-tol)=-(x-vv) --> y-hh+tol+x-vv
  #bl --> y-(hh+tol)=-(x-vv) --> y-hh-tol+x-vv
  #br --> y-(hh+tol)=x-vv --> y-hh-tol-x+vv
  inds<-coords[,2]-hh+
    (if(top) tol else -tol)+
    (if((top&left)|(!top&!left)) -coords[,1]+vv else coords[,1]-vv)!=
    0

  #parse out into individual blocks lying along focal line...
  runs<-rle(inds)
  #grab endpts of each block...
  ll<-c(0,cumsum(runs$lengths))[which(runs$values)]+1
  rr<-cumsum(runs$lengths)[runs$values]
  #correct cases with TRUE blocks bookending vector...I think this should work?
  if(ll[1]==1&rr[length(rr)]==nn){
    ll<-c(ll[length(ll)],ll[-c(1,length(ll))])
    rr<-rr[-length(rr)]
  }

  #get angles of endpts with respect to adjacent pts off focal line...
  ll.angs<-coords[ll,,drop=FALSE]-coords[(ll-2)%%nn+1,,drop=FALSE]
  ll.angs<-atan(ll.angs[,2]/ll.angs[,1])
  rr.angs<-coords[rr,,drop=FALSE]-coords[rr%%nn+1,,drop=FALSE]
  rr.angs<-atan(rr.angs[,2]/rr.angs[,1])

  #correcting angles when endpts touch border...
  #may have to play around with tolerances some more...
  ll.angs[abs(coords[ll,1])>(abs(vv)-100)]<-0
  rr.angs[abs(coords[rr,1])>(abs(vv)-100)]<-0
  ll.angs[abs(coords[ll,2])>(abs(vv)-100)]<-pi/2
  rr.angs[abs(coords[rr,2])>(abs(vv)-100)]<-pi/2

  #get indices between each endpt pair
  foo<-function(from,to,n){
    if(from>to){
      to<-n+to
    }
    (seq(from+1,length.out=to-from-1)-1)%%nn+1
  }
  tmp.seq<-seq_along(ll)
  mms<-lapply(tmp.seq,
              function(ii) foo(ll[ii],rr[ii],nn))

  #after much experimentation, best approach seems to be projecting
  #endpts to borders, then interpolating between them, then projecting
  #interpolations to border
  #otherwise there's too much ambiguity in whether slope angles should
  #rotate clockwise or counterclockwise
  ll.prj<-cbind(vv,
                tan(ll.angs)*(vv-coords[ll,1])+coords[ll,2])
  ll.prj.hh<-cbind((hh-coords[ll,2])/tan(ll.angs)+coords[ll,1],
                   hh)
  vv.dists<-sqrt(rowSums((ll.prj-coords[ll,,drop=FALSE])^2))
  hh.dists<-sqrt(rowSums((ll.prj.hh-coords[ll,,drop=FALSE])^2))
  tmp.inds<-hh.dists<vv.dists
  ll.prj[tmp.inds,]<-ll.prj.hh[tmp.inds,]

  rr.prj<-cbind(vv,
                tan(rr.angs)*(vv-coords[rr,1])+coords[rr,2])
  rr.prj.hh<-cbind((hh-coords[rr,2])/tan(rr.angs)+coords[rr,1],
                   hh)
  vv.dists<-sqrt(rowSums((rr.prj-coords[rr,,drop=FALSE])^2))
  hh.dists<-sqrt(rowSums((rr.prj.hh-coords[rr,,drop=FALSE])^2))
  tmp.inds<-hh.dists<vv.dists
  rr.prj[tmp.inds,]<-rr.prj.hh[tmp.inds,]

  #interpolation weights...
  ll.wgts<-lapply(tmp.seq,
                  function(ii)
                    sqrt(colSums((t(coords[mms[[ii]],,drop=FALSE])-coords[rr[ii],])^2)))
  rr.wgts<-lapply(tmp.seq,
                  function(ii)
                    sqrt(colSums((t(coords[mms[[ii]],,drop=FALSE])-coords[ll[ii],])^2)))
  #now just round them to edges!
  #now need 3 cases
  #if ll and rr lie on hh, round to hh
  #if rr and ll lie on vv, round to vv
  foo<-function(ii){
    tmp.prj<-(outer(ll.wgts[[ii]],ll.prj[ii,])+
                outer(rr.wgts[[ii]],rr.prj[ii,]))/
      rep(ll.wgts[[ii]]+rr.wgts[[ii]],2)

    # #I think a projection using the angle with the origin might be better, but this works for now
    # if(coords[ll[ii],2]==hh&coords[rr[ii],2]==hh){
    #   tmp.prj[,2]<-hh
    # }else if(coords[ll[ii],1]==vv&coords[rr[ii],1]==vv){
    #   tmp.prj[,1]<-vv
    # }else{
    #   dist.hh<-abs(tmp.prj[,2]-hh)
    #   dist.vv<-abs(tmp.prj[,1]-vv)
    #   tmp.prj[dist.hh<dist.vv,2]<-hh
    #   tmp.prj[dist.vv<dist.hh,1]<-vv
    # }
    #wait...tmp.prj will already be correct if both ll and rr lie on the same border
    #so you ONLY need to deal with case that they don't!
    #now using angle wrt origin for more even spacing of coords
    if(!(coords[ll[ii],2]==hh&coords[rr[ii],2]==hh)&
       !(coords[ll[ii],1]==vv&coords[rr[ii],1]==vv)){
      tmp.slopes<-tmp.prj[,2]/tmp.prj[,1]
      tmp.prj<-cbind(vv,
                     tmp.slopes*vv)
      tmp.inds<-abs(tmp.prj[,2])>abs(hh)
      tmp.prj[tmp.inds,]<-cbind(hh/tmp.slopes[tmp.inds],
                                hh)
    }
    tmp.prj
  }
  mm.prj<-lapply(tmp.seq,foo)

  #insert new coords back in
  coords[c(ll,unlist(mms,use.names=FALSE),rr),]<-rbind(ll.prj,
                                                       do.call(rbind,mm.prj),
                                                       rr.prj)
  x[,c("x","y")]<-coords

  #add in corner coordinate if it doesn't exist
  on.hh<-x[,"y"]==hh
  on.vv<-x[,"x"]==vv
  if(!any(on.hh&on.vv)){
    corner.check<-(on.hh[c(2:nn,1)]&on.vv)|(on.vv[c(2:nn,1)]&on.hh)
    if(any(corner.check)){
      tmp.inds<-which(corner.check) #should always be length 1
      tmp.tmp<-x[tmp.inds,,drop=FALSE]
      tmp.tmp[,c("x","y")]<-c(vv,hh)
      x<-rbind(x[seq_len(tmp.inds),,drop=FALSE],
               tmp.tmp,
               x[-seq_len(tmp.inds),,drop=FALSE])
    }
  }
  x
}

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
        tmp<-terra::aggregate(
          terra::union(
            x[i],
            terra::vect(
              .magic.corner.patcher(
                terra::geom(terra::union(x[i],NA.patch1)[3]),
                top=TRUE,left=TRUE,tol=NA_tol),"polygons"
            )
          )
        )
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
        tmp<-terra::aggregate(
          terra::union(
            x[i],
            terra::vect(
              .magic.corner.patcher(
                terra::geom(terra::union(x[i],NA.patch2)[3]),
                top=FALSE,left=FALSE,tol=NA_tol),"polygons"
            )
          )
        )
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
        tmp<-terra::aggregate(
          terra::union(
            x[i],
            terra::vect(
              .magic.corner.patcher(
                terra::geom(terra::union(x[i],SA.patch)[3]),
                top=FALSE,left=TRUE,tol=SA_tol),"polygons"
            )
          )
        )
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
        tmp<-terra::aggregate(
          terra::union(
            x[i],
            terra::vect(
              .magic.corner.patcher(
                terra::geom(terra::union(x[i],Asia.patch)[3]),
                top=TRUE,left=FALSE,tol=Asia_tol),"polygons"
            )
          )
        )
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

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
  #I think this should work, but it unfortunately might depend on whether...
  #...polygon is clockwise/counterclockwise rather than corner...
  # hh.ang<-if(top) pi/2 else -pi/2
  # vv.ang<-if(left) pi else 0

  ll.angs[abs(coords[ll,1])>(abs(vv)-100)]<-0
  rr.angs[abs(coords[rr,1])>(abs(vv)-100)]<-0
  #probs best to keep everything pi/2 and figure out best traversal later...
  ll.angs[abs(coords[ll,2])>(abs(vv)-100)]<-pi/2
  rr.angs[abs(coords[rr,2])>(abs(vv)-100)]<-pi/2
  # tmp.inds<-abs(coords[ll,2])>(abs(vv)-100)
  # ll.angs[tmp.inds&rr.angs<0]<- -pi/2
  # ll.angs[tmp.inds&rr.angs>0]<-pi/2
  # tmp.inds<-abs(coords[rr,2])>(abs(vv)-100)
  # rr.angs[tmp.inds&ll.angs<0]<- -pi/2
  # rr.angs[tmp.inds&ll.angs>0]<-pi/2

  #get angles for remaining points on focal line by interpolating between endpts
  #following function breaks if ll and rr are the same index...
  #probably don't want to do anything in that case? Will fix later
  foo<-function(from,to,n){
    if(from>to){
      to<-n+to
    }
    (seq(from+1,length.out=to-from-1)-1)%%nn+1
  }
  tmp.seq<-seq_along(ll)
  mms<-lapply(tmp.seq,
              function(ii) foo(ll[ii],rr[ii],nn))

  #10/29: alright, the above is fine
  #the problem is making sure you get the correct "angle gradient"
  #note: there's an "identity issue" of sorts--you don't need the entire unit circle of angles, just it's right side...
  #yeah, most appropriate to use atan rather than atan2 here...

  #I THINK weighted averages might be fine now???
  #Not quite...but I think you're getting close
  #Still have some cases where it seems like the angles shouldn't be going in the directions they are
  #Yeah, "minimum intervening angle" idea doesn't always work unfortunately
  #Having trouble finding systematic way to figure out "polarity" of angle traversal...
  #Can't tell if I'm just not finding the right criteria or if this is impossible to do simply...
  #alternate idea: project the corners (ll and rr)
  #then see if you can interpolate between, rounding to nearest edge...
  #probably the best way to do it...

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
  #otherwise, used edge proximity
  #could combine with earlier foo to get mm in the first place using this approach...
  foo<-function(ii){
    tmp.prj<-(outer(ll.wgts[[ii]],ll.prj[ii,])+
                outer(rr.wgts[[ii]],rr.prj[ii,]))/
      rep(ll.wgts[[ii]]+rr.wgts[[ii]],2)

    #I think a projection using the angle with the origin might be better, but this works for now
    if(coords[ll[ii],2]==hh&coords[rr[ii],2]==hh){
      tmp.prj[,2]<-hh
    }else if(coords[ll[ii],1]==vv&coords[rr[ii],1]==vv){
      tmp.prj[,1]<-vv
    }else{
      dist.hh<-abs(tmp.prj[,2]-hh)
      dist.vv<-abs(tmp.prj[,1]-vv)
      tmp.prj[dist.hh<dist.vv,2]<-hh
      tmp.prj[dist.vv<dist.hh,1]<-vv
    }
    tmp.prj
  }
  mm.prj<-lapply(tmp.seq,foo)

  # all.inds<-c(ll,unlist(mms,use.names=FALSE),rr)
  # all.pts<-coords[all.inds,,drop=FALSE]
  # all.angs<-tan(c(ll.angs,unlist(mm.angs,use.names=FALSE),rr.angs))
  coords[c(ll,unlist(mms,use.names=FALSE),rr),]<-rbind(ll.prj,
                                                       do.call(rbind,mm.prj),
                                                       rr.prj)


  #almost there! the last thing you need to do is somehow make it so that this always comes up with reasonable interpolation
  #depends on corner and orientation of ll vs rr...
  #idea: use angle with respect to origin???
  # foo<-function(ii){
  #   tmp.angs<-atan(coords[mms[[ii]],2]/coords[mms[[ii]],1])
  #   tmp.n<-length(tmp.angs)
  #
  #   (tmp.angs-min(tmp.angs[c(1,tmp.n)]))/diff(range(tmp.angs[c(1,tmp.n)]))*diff(range(ll.angs[ii],rr.angs[ii]))+min(ll.angs[ii],rr.angs[ii])
  #
  #
  #
  #
  #
  #
  #   #scale and shift to fit with actual ll and rr angs...
  #   #should be more robust
  #   #almost there, but need to get "flipping" baked in
  #
  #   # if(tmp.angs[tmp.n]>tmp.angs[1]&rr.angs[ii]<ll.angs[ii]){
  #   #   #things should go clockwise from ll to rr but are going counterclockwise
  #   #   rr.angs[ii]<-rr.angs[ii]+pi
  #   # }else  if(ll.angs[ii]<rr.angs[ii]){
  #   #   #things should go counter-clockwise from ll to rr but are going clockwise
  #   #   ll.angs[ii]<-ll.angs[ii]+pi
  #   # }
  #   #
  #   # (tmp.angs-ll.o.angs[ii])/
  #   #   (rr.o.angs[ii]-ll.o.angs[ii])*
  #   #   (rr.angs[ii]-ll.angs[ii])+
  #   #   ll.angs[ii]
  #   #
  #   # tan((tmp.angs-ll.o.angs[ii])/
  #   #   (rr.o.angs[ii]-ll.o.angs[ii])*
  #   #   (rr.angs[ii]-ll.angs[ii])+
  #   #   ll.angs[ii])
  # }
  # mm.angs<-lapply(tmp.seq,foo)
  # #put it all together...
  # all.inds<-unlist(mms,use.names=FALSE)
  # all.pts<-coords[all.inds,,drop=FALSE]
  # all.angs<-tan(unlist(mm.angs,use.names=FALSE))

  #not sure how well this will work...
  # tmp.prj<-cbind(vv,
  #                all.angs*(vv-all.pts[,1])+all.pts[,2])
  # tmp.prj.hh<-cbind((hh-all.pts[,2])/all.angs+all.pts[,1],
  #                   hh)
  # tmp.inds<-abs(tmp.prj[,2])>abs(hh)
  # tmp.prj[tmp.inds,]<-cbind((hh-all.pts[tmp.inds,2])/all.angs[tmp.inds]+all.pts[tmp.inds,1],
  #                           hh)

  # vv.dists<-sqrt(rowSums((tmp.prj-all.pts)^2))
  # hh.dists<-sqrt(rowSums((tmp.prj.hh-all.pts)^2))
  # tmp.inds<-hh.dists<vv.dists
  # tmp.prj[tmp.inds,]<-tmp.prj.hh[tmp.inds,]

  # coords[c(ll,unlist(mms,use.names=FALSE),rr),]<-rbind(ll.prj,
  #                                                      do.call(rbind,mm.prj),
  #                                                      rr.prj)
  #potential issue--whether point is directly projected precisely into corner depends on resolution
  #here's a potential fix?
  x[,c("x","y")]<-coords
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

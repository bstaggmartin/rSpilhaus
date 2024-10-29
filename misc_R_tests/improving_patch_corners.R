.magic.corner.patcher<-function(x,
                                top,
                                left,
                                tol){
  #let's try getting the angles between adjacent segments
  #hmmm...maybe find points that lie on patch line?
  tmpy<-testy[,c("x","y")]
  nn<-nrow(tmpy)
  # test.inds<-tmpy[,1]%in%terra::geom(NA.patch1)[,"x"]&
  #   tmpy[,2]%in%terra::geom(NA.patch1)[,"y"]
  #more robust test
  test.inds<-tmpy[,1]-tmpy[,2]+2*lim-NA_tol==0
  #yay! that seems to work here
  runs<-rle(test.inds)
  #grab endpts of each "TRUE" block...
  ll<-c(0,cumsum(runs$lengths))[which(runs$values)]+1
  rr<-cumsum(runs$lengths)[runs$values]
  #correct cases with TRUE blocks bookending...I think this should work?
  if(ll[1]==1&rr[length(rr)]==nn){
    ll<-c(ll[length(ll)],ll[-c(1,length(ll))])
    rr<-rr[-length(rr)]
  }
  ll.angs<-tmpy[ll,,drop=FALSE]-tmpy[(ll-2)%%nn+1,,drop=FALSE]
  ll.angs<-atan2(ll.angs[,2],ll.angs[,1])
  rr.angs<-tmpy[rr,,drop=FALSE]-tmpy[rr%%nn+1,,drop=FALSE]
  rr.angs<-atan2(rr.angs[,2],rr.angs[,1])
  #better angs for when endpts touch border...
  ll.angs[tmpy[ll,1]==-lim]<-0
  ll.angs[tmpy[ll,2]==lim]<- -pi/2
  rr.angs[tmpy[rr,1]==-lim]<-0
  rr.angs[tmpy[rr,2]==lim]<- -pi/2

  #breaks if ll and rr are the same index...
  #probably don't want to do anything in that case?
  foo<-function(from,to,n){
    if(from>to){
      to<-n+to
    }
    (seq(from+1,length.out=to-from-1)-1)%%nn+1
  }
  tmp.seq<-seq_along(ll)
  mms<-lapply(tmp.seq,
              function(ii) foo(ll[ii],rr[ii],nn))
  ll.wgts<-lapply(tmp.seq,
                  function(ii)
                    sqrt(colSums((t(tmpy[mms[[ii]],,drop=FALSE])-tmpy[rr[ii],])^2)))
  rr.wgts<-lapply(tmp.seq,
                  function(ii)
                    sqrt(colSums((t(tmpy[mms[[ii]],,drop=FALSE])-tmpy[ll[ii],])^2)))
  mm.angs<-lapply(tmp.seq,
                  function(ii)
                    (ll.wgts[[ii]]*ll.angs[ii]+rr.wgts[[ii]]*rr.angs[ii])/(ll.wgts[[ii]]+rr.wgts[[ii]]))
  #put it all together...
  all.inds<-c(ll,unlist(mms,use.names=FALSE),rr)
  all.pts<-tmpy[all.inds,,drop=FALSE]
  all.angs<-c(ll.angs,unlist(mm.angs,use.names=FALSE),rr.angs)
  tmp.prj<-cbind(-lim,
                 all.angs*(-lim-all.pts[,1])+all.pts[,2])
  tmp.inds<-tmp.prj[,2]>lim
  tmp.prj[tmp.inds,]<-cbind((lim-all.pts[tmp.inds,2])/all.angs[tmp.inds]+all.pts[tmp.inds,1],lim)
  tmpy[all.inds,]<-tmp.prj
  #potential issue--whether point is directly projected precisely into corner depends on resolution
  #here's a potential fix?
  testy[,c("x","y")]<-tmpy
  tt<-testy[,"y"]==lim
  ll<-testy[,"x"]==-lim
  if(!any(tt&ll)){
    corner.check<-(tt[c(2:nn,1)]&ll)|(ll[c(2:nn,1)]&tt)
    if(any(corner.check)){
      tmp.inds<-which(corner.check) #should always be length 1
      tmp.tmp<-testy[tmp.inds,,drop=FALSE]
      tmp.tmp[,c("x","y")]<-c(-lim,lim)
      testy<-rbind(testy[seq_len(tmp.inds),,drop=FALSE],
                   tmp.tmp,
                   testy[-seq_len(tmp.inds),,drop=FALSE])
    }
  }
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
        # tmp<-terra::aggregate(terra::union(x[i],NA.patch1))
        # terra::values(tmp)<-terra::values(x)[i,,drop=FALSE]
        # x<-rbind(x,tmp)

        #this all seems to be working!!! Though I'm sure there will be edge cases...
        #should break out into it's own function and figure out how to adapt to different corners...
        #other thing to consider: tiling--errors should be minimal...
        #...and this seems like a good starting point...
        #...BUT ideally should ensure things match up across border...

        testy<-terra::geom(terra::union(x[i],NA.patch1)[3]) #intersection of patch and geometry...
        #let's try getting the angles between adjacent segments
        #hmmm...maybe find points that lie on patch line?
        tmpy<-testy[,c("x","y")]
        nn<-nrow(tmpy)
        # test.inds<-tmpy[,1]%in%terra::geom(NA.patch1)[,"x"]&
        #   tmpy[,2]%in%terra::geom(NA.patch1)[,"y"]
        #more robust test
        test.inds<-tmpy[,1]-tmpy[,2]+2*lim-NA_tol==0
        #yay! that seems to work here
        runs<-rle(test.inds)
        #grab endpts of each "TRUE" block...
        ll<-c(0,cumsum(runs$lengths))[which(runs$values)]+1
        rr<-cumsum(runs$lengths)[runs$values]
        #correct cases with TRUE blocks bookending...I think this should work?
        if(ll[1]==1&rr[length(rr)]==nn){
          ll<-c(ll[length(ll)],ll[-c(1,length(ll))])
          rr<-rr[-length(rr)]
        }
        ll.angs<-tmpy[ll,,drop=FALSE]-tmpy[(ll-2)%%nn+1,,drop=FALSE]
        ll.angs<-atan2(ll.angs[,2],ll.angs[,1])
        rr.angs<-tmpy[rr,,drop=FALSE]-tmpy[rr%%nn+1,,drop=FALSE]
        rr.angs<-atan2(rr.angs[,2],rr.angs[,1])
        #better angs for when endpts touch border...
        ll.angs[tmpy[ll,1]==-lim]<-0
        ll.angs[tmpy[ll,2]==lim]<- -pi/2
        rr.angs[tmpy[rr,1]==-lim]<-0
        rr.angs[tmpy[rr,2]==lim]<- -pi/2

        #breaks if ll and rr are the same index...
        #probably don't want to do anything in that case?
        foo<-function(from,to,n){
          if(from>to){
            to<-n+to
          }
          (seq(from+1,length.out=to-from-1)-1)%%nn+1
        }
        tmp.seq<-seq_along(ll)
        mms<-lapply(tmp.seq,
                    function(ii) foo(ll[ii],rr[ii],nn))
        ll.wgts<-lapply(tmp.seq,
                        function(ii)
                          sqrt(colSums((t(tmpy[mms[[ii]],,drop=FALSE])-tmpy[rr[ii],])^2)))
        rr.wgts<-lapply(tmp.seq,
                        function(ii)
                          sqrt(colSums((t(tmpy[mms[[ii]],,drop=FALSE])-tmpy[ll[ii],])^2)))
        mm.angs<-lapply(tmp.seq,
                        function(ii)
                          (ll.wgts[[ii]]*ll.angs[ii]+rr.wgts[[ii]]*rr.angs[ii])/(ll.wgts[[ii]]+rr.wgts[[ii]]))
        #put it all together...
        all.inds<-c(ll,unlist(mms,use.names=FALSE),rr)
        all.pts<-tmpy[all.inds,,drop=FALSE]
        all.angs<-c(ll.angs,unlist(mm.angs,use.names=FALSE),rr.angs)
        tmp.prj<-cbind(-lim,
                       all.angs*(-lim-all.pts[,1])+all.pts[,2])
        tmp.inds<-tmp.prj[,2]>lim
        tmp.prj[tmp.inds,]<-cbind((lim-all.pts[tmp.inds,2])/all.angs[tmp.inds]+all.pts[tmp.inds,1],lim)
        tmpy[all.inds,]<-tmp.prj
        #potential issue--whether point is directly projected precisely into corner depends on resolution
        #here's a potential fix?
        testy[,c("x","y")]<-tmpy
        tt<-testy[,"y"]==lim
        ll<-testy[,"x"]==-lim
        if(!any(tt&ll)){
          corner.check<-(tt[c(2:nn,1)]&ll)|(ll[c(2:nn,1)]&tt)
          if(any(corner.check)){
            tmp.inds<-which(corner.check) #should always be length 1
            tmp.tmp<-testy[tmp.inds,,drop=FALSE]
            tmp.tmp[,c("x","y")]<-c(-lim,lim)
            testy<-rbind(testy[seq_len(tmp.inds),,drop=FALSE],
                         tmp.tmp,
                         testy[-seq_len(tmp.inds),,drop=FALSE])
          }
        }
        tmp<-terra::aggregate(terra::union(x[i],terra::vect(testy,"polygons")))
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

#' @export
expand_borders<-function(x,
                         amount=0.06,
                         buffer_tol=15,
                         prettify=FALSE,
                         frame=FALSE,
                         frame_tol=5,
                         revalidate_geoms=TRUE){

  if(amount>1|amount<0){
    stop("amount must be between 0 and 1")
  }
  if(amount<0.05&prettify&frame){
    stop("amount must be greater than 0.05 if prettify and frame are TRUE")
  }

  #for sorting back into proper order
  vals<-terra::values(x)
  if(!length(vals)){
    vals<-data.frame("TEMPORARY_ORDER"=seq_along(x))
  }else{
    vals[["TEMPORARY_ORDER"]]<-seq_along(x)
  }
  terra::values(x)<-vals

  geom.type<-terra::geomtype(x)
  lim<-11825474
  tmp<-terra::geom(x)
  rotate<-function(x,ang,org=c(0,0)){
    x[,c("x","y")]<-matrix(c(cos(ang)*(x[,"x"]-org[1])-
                               sin(ang)*(x[,"y"]-org[2])+
                               org[1],
                             sin(ang)*(x[,"x"]-org[1])+
                               cos(ang)*(x[,"y"]-org[2])+
                               org[2]),
                           ncol=2)
    x
  }
  tmp.prj<-terra::vect(do.call(rbind,list(rotate(tmp,pi,c(-lim,lim)),
                                          rotate(tmp,pi/2,c(-lim,-lim)),
                                          rotate(tmp,pi,c(-lim,-lim)),
                                          rotate(tmp,-pi/2,c(lim,lim)),
                                          rotate(tmp,0,c(0,0)),
                                          rotate(tmp,-pi/2,c(-lim,-lim)),
                                          rotate(tmp,pi,c(lim,lim)),
                                          rotate(tmp,pi/2,c(lim,lim)),
                                          rotate(tmp,pi,c(lim,-lim)))),
                       geom.type)

  nn<-length(x)
  out.prj<-NULL
  for(i in 1:nn){
    if(geom.type=="polygons"){
      tmp<-terra::aggregate(terra::buffer(tmp.prj[(0:8)*nn+i],width=buffer_tol))
    }else{
      tmp<-terra::aggregate(tmp.prj[(0:8)*nn+i])
    }
    terra::values(tmp)<-terra::values(x)[i,,drop=FALSE]
    if(is.null(out.prj)){
      out.prj<-tmp
    }else{
      out.prj<-rbind(out.prj,tmp)
    }
  }

  out.prj<-terra::sort(out.prj,v="TEMPORARY_ORDER")
  vals<-terra::values(out.prj)
  vals[["TEMPORARY_ORDER"]]<-NULL
  terra::values(out.prj)<-vals

  #get cropping polygons...
  if(prettify){
    #constraints lifted from ricardo's code...
    # cols<-rainbow(6)
    # terra::plot(out)
    # #lines
    # abline(1.089e7,-0.176,col=cols[1])
    # abline(-0.984e7/0.565,-1/0.565,col=cols[2])
    # abline(-1.378e7,0.46,col=cols[3])
    # abline(-1.274e7/0.172,1/0.172,col=cols[4])
    # abline(1e7,-0.5,col=cols[5])
    # abline(2.3e7,1,col=cols[6])
    # #horizontal/vertical lines
    # abline(v=-1.114e7,h=0.29e7,col=cols[1]) #interesects with line 2
    # abline(v=-1.17e7,h=0.39e7,col=cols[2]) #intersects with line 2
    # abline(v=0.295e7,h=-1.21e7,col=cols[3]) #intersect with line 3
    # abline(v=0.312e7,h=-1.2e7,col=cols[4]) #intersect with line 3
    # abline(v=0.4e7,h=-1.16e7,col=cols[5]) #intersects with line 3
    # abline(v=0.45e7,-1.11e7,col=cols[6]) #intersects with line 3

    #converting the above to a polygon...
    lim2<-1.1*lim
    cropper<-
      rbind(
        #left lim and line 6
        c(-lim2,2.3e7-lim2),
        #line 6 and line 1
        c((2.3e7-1.089e7)/(-0.176-1),1.089e7-0.176*(2.3e7-1.089e7)/(-0.176-1)),
        #line 1 and line 5
        c((1e7-1.089e7)/(-0.176+0.5),1.089e7-0.176*(1e7-1.089e7)/(-0.176+0.5)),
        #line 5 and right lim
        c(lim2,1e7-0.5*lim2),
        #right lim and line 4
        c(lim2,-1.274e7/0.172+1/0.172*lim2),
        #line 4 and line 3
        c((-1.274e7/0.172+1.378e7)/(0.46-1/0.172),-1.378e7+0.46*(-1.274e7/0.172+1.378e7)/(0.46-1/0.172)),
        #line 3 and horizontal line 6
        c((-1.11e7+1.378e7)/0.46,-1.11e7),
        #horizontal and vertical line 6
        c(0.45e7,-1.11e7),
        # #vertical line 6 and line 3
        # c(0.45e7,-1.378e7+0.46*0.45e7),
        # #line 3 and horizontal line 5
        # c((-1.16e7+1.378e7)/0.46,-1.16e7),
        #vertical line 6 and horizontal line 5
        c(0.45e7,-1.16e7),
        #horizontal and vertical line 5
        c(0.4e7,-1.16e7),
        #vertical line 5 and line 3
        c(0.4e7,-1.378e7+0.46*0.4e7),
        #line 3 and horizontal line 4
        c((-1.2e7+1.378e7)/0.46,-1.2e7),
        #horizontal and vertical line 4
        c(0.312e7,-1.2e7),
        # #vertical line 4 and line 3
        # c(0.312e7,-1.378e7+0.46*0.312e7),
        # #line 3 and horizontal line 3
        # c((-1.21e7+1.378e7)/0.46,-1.21e7),
        #vertical line 4 and horizontal line 3
        c(0.312e7,-1.21e7),
        #horizontal and vertical line 3
        c(0.295e7,-1.21e7),
        #vertical line 3 and line 3
        c(0.295e7,-1.378e7+0.46*0.295e7),
        #line 3 and bottom lim
        c((-lim2+1.378e7)/0.46,-lim2),
        #bottom lim and line 2
        c(lim2*0.565-0.984e7,-lim2),
        #line 2 and vertical line 1
        c(-1.114e7,-0.984e7/0.565+1.114e7/0.565),
        #vertical and horizontal line 1
        c(-1.114e7,0.29e7),
        #horizontal line 1 and line 2
        c(-0.565*0.29e7-0.984e7,0.29e7),
        #line 2 and vertical line 2
        c(-1.17e7,-0.984e7/0.565+1.17e7/0.565),
        #vertical and horizontal line 2
        c(-1.17e7,0.39e7),
        #horizontal line 2 and line 2
        c(-0.565*0.39e7-0.984e7,0.39e7),
        #line 2 and left lim
        c(-lim2,-0.984e7/0.565+lim2/0.565)
      )
  }else{
    wd<-amount*2*lim
    cropper<-cbind(c(-lim-wd,
                     -lim-wd,
                     lim+wd,
                     lim+wd),
                   c(-lim-wd,
                     lim+wd,
                     lim+wd,
                     -lim-wd))
  }

  #finally crop the dang thing
  out.prj<-terra::crop(out.prj,terra::vect(cropper,"polygons"))

  #adding frame if prettify is TRUE...
  if(frame&prettify&geom.type=="polygons"){
    wd<-amount*2*lim
    rectangle<-cbind(c(-lim-wd,
                       -lim-wd,
                       lim+wd,
                       lim+wd),
                     c(-lim-wd,
                       lim+wd,
                       lim+wd,
                       -lim-wd))
    rectangle<-terra::buffer(terra::erase(terra::vect(rectangle,"polygons"),
                                          terra::vect(cropper,"polygons")),
                             width=frame_tol)

    vals<-terra::values(out.prj)
    if(!length(vals)){
      vals<-data.frame("TEMPORARY_ORDER"=seq_along(out.prj))
    }else{
      vals[["TEMPORARY_ORDER"]]<-seq_along(out.prj)
    }
    terra::values(out.prj)<-vals

    inds<-which(terra::is.related(out.prj,rectangle,"intersects"))
    if(length(inds)>0){
      for(i in inds){
        tmp<-terra::aggregate(terra::union(out.prj[i],rectangle))
        terra::values(tmp)<-terra::values(out.prj)[i,,drop=FALSE]
        out.prj<-rbind(out.prj,tmp)
      }
      out.prj<-out.prj[-inds]
    }

    out.prj<-terra::sort(out.prj,v="TEMPORARY_ORDER")
    vals<-terra::values(out.prj)
    vals[["TEMPORARY_ORDER"]]<-NULL
    terra::values(out.prj)<-vals
  }

  if(revalidate_geoms){
    terra::buffer(out.prj,width=0)
  }else{
    out.prj
  }
}

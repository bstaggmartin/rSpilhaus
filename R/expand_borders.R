#TODO:
# - use new crack repair technique with convex hulls so buffering is no longer necessary
#    - 10/31 (oooh spooky) --> check!
# - adjust prettify border to work with 50 m res land from rnaturalearth (probs a good default)
# - adapt corner patching technique to extend diff pieces of geometry along border
#    - use angle wrt to origin to avoid self-intersection
#maybe use angle wrt to...
# abline(v=-3.75e6,h=2e6)
#ended up making it an argument "center"
#11/6 --> got it all working with no artifacts at any rnaturalearth resolution for now...
#looks ugly with multiple pieces of geometry and frame=TRUE, but it will suffice for now...
#merge_neg_space controls whether or not negative space in frame is merged with rest of frame
#(only does so with 1 piece of geometry--otherwise which one do we merge it with?)
#frame_tol controls how much frame is extended into geometries to detect/resolve intersections
#border_tol controls the distance to the frame border used to detect points of polygons lying on borders
#(and thus which points need to be extended to border)

#it all seems to be working reasonbaly well now, but should really clean up the function!

#' @export
expand_borders<-function(x,
                         amount=0.06,
                         seam_width=1000,
                         prettify=FALSE,
                         frame=FALSE,
                         dissolve=TRUE,
                         frame_tol=15,
                         merge_neg_space=TRUE,
                         neg_space_tol=5,
                         border_tol=0.1,
                         cent=c(-3e6,3e6),
                         revalidate_geoms=TRUE,
                         sample_method="simple",
                         patch_width=3,
                         patch_width_inc=0,
                         patch_method="modal",
                         max_patch_iter=1e4){

  if(amount>1|amount<0){
    stop("amount must be between 0 and 1")
  }
  if(amount<0.05&prettify&frame){
    stop("amount must be greater than 0.05 if prettify and frame are TRUE")
  }

  lim<-11825474
  ext.amount<-2*amount*lim

  if(inherits(x,"SpatVector")){
    #for sorting back into proper order
    vals<-terra::values(x)
    if(!length(vals)){
      vals<-data.frame("TEMPORARY_ORDER"=seq_along(x))
    }else{
      vals[["TEMPORARY_ORDER"]]<-seq_along(x)
    }
    terra::values(x)<-vals

    geom.type<-terra::geomtype(x)
    tmp.crop<-function(x,left=NA,top=NA,lab){
      if(is.na(left)){
        x0<- -lim
        x1<-lim
      }else if(left){
        x0<- -lim
        x1<- -lim+ext.amount
      }else{
        x0<-lim-ext.amount
        x1<-lim
      }
      if(is.na(top)){
        y0<- -lim
        y1<-lim
      }else if(top){
        y0<-lim-ext.amount
        y1<-lim
      }else{
        y0<- -lim
        y1<- -lim+ext.amount
      }
      out<-terra::crop(x,terra::vect(cbind(c(x0,x1),
                                           c(y0,y1)),
                                     "points"),
                       ext=TRUE)
      vals<-terra::values(out)$TEMPORARY_ORDER
      out<-terra::geom(out)
      cbind(out,
            "TEMPORARY_ORDER"=vals[out[,"geom"]],
            "TEMPORARY_LABEL"=if(nrow(out)) lab else NULL)
    }
    rotate<-function(x,ang,org=c(0,0)){
      if(nrow(x)>0&ang!=0){
        x[,c("x","y")]<-matrix(c(cos(ang)*(x[,"x"]-org[1])-
                                   sin(ang)*(x[,"y"]-org[2])+
                                   org[1],
                                 sin(ang)*(x[,"x"]-org[1])+
                                   cos(ang)*(x[,"y"]-org[2])+
                                   org[2]),
                               ncol=2)
      }
      x
    }
    #may be a more efficient way to do this, but fine for now...
    tmp.prj<-do.call(rbind,list(rotate(tmp.crop(x,left=TRUE,top=TRUE,lab=1),
                                       pi,c(-lim,lim)),
                                rotate(tmp.crop(x,left=NA,top=FALSE,lab=2),
                                       pi/2,c(-lim,-lim)),
                                rotate(tmp.crop(x,left=TRUE,top=FALSE,lab=3),
                                       pi,c(-lim,-lim)),
                                rotate(tmp.crop(x,left=FALSE,top=NA,lab=4),
                                       -pi/2,c(lim,lim)),
                                rotate(tmp.crop(x,lab=5),
                                       0,c(0,0)),
                                rotate(tmp.crop(x,left=TRUE,top=NA,lab=6),
                                       -pi/2,c(-lim,-lim)),
                                rotate(tmp.crop(x,left=FALSE,top=TRUE,lab=7),
                                       pi,c(lim,lim)),
                                rotate(tmp.crop(x,left=NA,top=TRUE,lab=8),
                                       pi/2,c(lim,lim)),
                                rotate(tmp.crop(x,left=FALSE,top=FALSE,lab=9),
                                       pi,c(lim,-lim))))
    #needed to fix this up to keep more craeful track of geometry...
    tmp.pos<-which(!duplicated(tmp.prj[,c("TEMPORARY_ORDER","TEMPORARY_LABEL")]))
    tmp.geoms<-inverse.rle(list(values=seq_along(tmp.pos),
                                lengths=c(tmp.pos[-1],nrow(tmp.prj)+1)-tmp.pos))
    geom.ids<-tmp.prj[tmp.pos,"TEMPORARY_ORDER"]
    tmp.prj<-tmp.prj[,1:(ncol(tmp.prj)-2),drop=FALSE]
    tmp.prj[,"geom"]<-tmp.geoms
    tmp.prj<-terra::vect(tmp.prj,
                         type=geom.type)
    seams<-c(-lim-ext.amount,lim+ext.amount)
    seams<-terra::vect(list(cbind(seams,-lim),
                            cbind(seams,lim),
                            cbind(-lim,seams),
                            cbind(lim,seams)),
                       "lines")
    seams<-terra::buffer(terra::aggregate(seams),
                         seam_width,
                         capstyle="square")

    nn<-length(x)
    out.prj<-vector("list",nn)
    for(i in 1:length(x)){
      tmp<-tmp.prj[geom.ids==i]
      if(length(tmp)>1){
        tmp<-tmp1<-terra::aggregate(tmp,dissolve=dissolve)
        if(dissolve&geom.type=="polygons"){
          if(terra::is.related(tmp1,seams,"intersects")){
            uni<-terra::union(tmp1,seams)
            tmp<-terra::disagg(uni[3])
            uni<-terra::aggregate(uni)
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
            tmp<-terra::aggregate(terra::union(tmp1,tmp))
          }
        }
      }
      terra::values(tmp)<-terra::values(x)[i,,drop=FALSE]
      out.prj[[i]]<-tmp
      # cat("\n",i)
    }
    out.prj<-do.call(rbind,out.prj)

    out.prj<-terra::sort(out.prj,v="TEMPORARY_ORDER")
    vals<-terra::values(out.prj)
    vals[["TEMPORARY_ORDER"]]<-NULL
    terra::values(out.prj)<-vals
  }else if(inherits(x,"SpatRaster")){

    geom.type<-NA

    dims<-dim(x)
    dims[c(1,2)]<-dims[c(2,1)]
    rez<-terra::res(x)

    tmp<-array(terra::values(x),
               dims)

    # nx<-floor(ext.amount/rez[1])
    # xseq<-seq_len(nx)
    # ny<-floor(ext.amount/rez[2])
    # yseq<-seq_len(ny)
    # new.ext<-c(-lim-(nx-0.5)*rez[1],
    #            lim+(nx-0.5)*rez[1],
    #            -lim-(ny-0.5)*rez[2],
    #            lim+(ny-0.5)*rez[2])

    #let's do what you did for the first time above..one at a time, but
    #"respecting" the funky order of the raster
    #so starting with bottom left corner...
    #ISSUE: differing resolution across x and y axes means that...
    #you can't just strictly mosaic things
    #suggests that you have to use an extract based method...

    # bl<-tmp[rev(yseq),dims[2]-xseq+1,,drop=FALSE]
    # bb<-aperm(tmp[,dims[2]-rev(xseq)+1,,drop=FALSE],
    #           c(2,1,3))

    xx<-terra::xFromCol(x,seq_len(dims[1]))
    yy<-terra::yFromRow(x,rev(seq_len(dims[2])))
    xxx<-seq(rez[1]/2,ext.amount,rez[1]) #new horizontal/x/row coords
    yyy<-seq(rez[2]/2,ext.amount,rez[2]) #new vertical/y/col coords
    nx<-length(xxx)
    ny<-length(yyy)

    # image(tmp[,500:1,1])
    # tmp.mat<-matrix(NA,600,ny)
    # vexp.tmp<-cbind(tmp.mat,tmp[,,1],tmp.mat)
    # tmp.mat<-matrix(NA,nx,500)
    # hexp.tmp<-rbind(tmp.mat,tmp[,,1],tmp.mat)
    # tmp.mat<-matrix(NA,nx,ny)

    #topleft
    tmp.coords<-expand.grid(-lim+xxx,lim-yyy)
    tl<-array(unlist(terra::extract(x,as.matrix(tmp.coords),method=sample_method),
                     use.names=FALSE),
              c(nx,ny,dims[3]))
    tl<-tl[rev(seq_len(nx)),rev(seq_len(ny)),,drop=FALSE]

    #left
    tmp.coords<-expand.grid(yy,-lim+xxx)
    ll<-array(unlist(terra::extract(x,as.matrix(tmp.coords),method=sample_method),
                     use.names=FALSE),
              c(dims[2],nx,dims[3]))
    ll<-aperm(ll,c(2,1,3))[rev(seq_len(nx)),rev(seq_len(dims[2])),,drop=FALSE]
    # image(rbind(ll[,,1],tmp[,,1])) #works

    # image(rbind(cbind(tl[,,1],ll[,,1],tmp.mat),
    #             vexp.tmp)) #works

    #bottomleft
    tmp.coords<-expand.grid(-lim+xxx,-lim+yyy)
    bl<-array(unlist(terra::extract(x,as.matrix(tmp.coords),method=sample_method),
                     use.names=FALSE),
              c(nx,ny,dims[3]))
    bl<-bl[rev(seq_len(nx)),,,drop=FALSE]


    #top
    tmp.coords<-expand.grid(lim-yyy,xx)
    tt<-array(unlist(terra::extract(x,as.matrix(tmp.coords),method=sample_method),
                     use.names=FALSE),
              c(ny,dims[1],dims[3]))
    tt<-aperm(tt,c(2,1,3))[,rev(seq_len(ny)),,drop=FALSE]
    # image(cbind(tt[,,1],tmp[,,1])) #works

    #bottom
    tmp.coords<-expand.grid(-lim+yyy,xx)
    bb<-array(unlist(terra::extract(x,as.matrix(tmp.coords),method=sample_method),
                     use.names=FALSE),
              c(ny,dims[1],dims[3]))
    bb<-aperm(bb,c(2,1,3))
    # image(cbind(tt[,,1],tmp[,,1],bb[,,1])) #works
    # image(cbind(rbind(tl[,,1],tt[,,1],tmp.mat),
    #             hexp.tmp,
    #             rbind(bl[,,1],bb[,,1],tmp.mat)),
    #       col=hcl.colors(100),
    #       breaks=seq(min(hexp.tmp,na.rm=TRUE),
    #                  max(hexp.tmp,na.rm=TRUE),
    #                  length.out=101)) #works

    #topright
    tmp.coords<-expand.grid(lim-xxx,lim-yyy)
    tr<-array(unlist(terra::extract(x,as.matrix(tmp.coords),method=sample_method),
                     use.names=FALSE),
              c(nx,ny,dims[3]))
    tr<-tr[,rev(seq_len(ny)),,drop=FALSE]

    #right
    tmp.coords<-expand.grid(yy,lim-xxx)
    rr<-array(unlist(terra::extract(x,as.matrix(tmp.coords),method=sample_method),
                     use.names=FALSE),
              c(dims[2],nx,dims[3]))
    rr<-aperm(rr,c(2,1,3))[,rev(seq_len(dims[2])),,drop=FALSE]
    # image(rbind(ll[,,1],tmp[,,1],rr[,,1])) #works

    #bottomright
    tmp.coords<-expand.grid(lim-xxx,-lim+yyy)
    br<-array(unlist(terra::extract(x,as.matrix(tmp.coords),method=sample_method),
                     use.names=FALSE),
              c(nx,ny,dims[3]))
    # image(rbind(cbind(tl[,,1],ll[,,1],tmp.mat),
    #             vexp.tmp,
    #             cbind(tr[,,1],rr[,,1],br[,,1]))) #works

    # image(rbind(cbind(tl[,,1],ll[,,1],bl[,,1]),
    #             cbind(tt[,,1],tmp[,,1],bb[,,1]),
    #             cbind(tr[,,1],rr[,,1],br[,,1])),
    #             col=hcl.colors(100),
    #             breaks=seq(min(hexp.tmp,na.rm=TRUE),
    #                        max(hexp.tmp,na.rm=TRUE),
    #                        length.out=101)) #seems to all be working...

    #now to put it all together...
    #first rbind top, middle, and bottom...
    ttt<-array(c(aperm(tl,c(2,3,1)),
                 aperm(tt,c(2,3,1)),
                 aperm(tr,c(2,3,1))),
               c(ny,dims[3],dims[1]+2*nx))
    mm<-array(c(aperm(ll,c(2,3,1)),
                aperm(tmp,c(2,3,1)),
                aperm(rr,c(2,3,1))),
              c(dims[2],dims[3],dims[1]+2*nx))
    bbb<-array(c(aperm(bl,c(2,3,1)),
                 aperm(bb,c(2,3,1)),
                 aperm(br,c(2,3,1))),
               c(ny,dims[3],dims[1]+2*nx))
    #now cbind left, middle, and right...
    out.prj<-x
    terra::nrow(out.prj)<-dims[2]+2*ny
    terra::ncol(out.prj)<-dims[1]+2*nx
    terra::ext(out.prj)<-c(-lim-xxx[nx]-rez[1]/2,
                           lim+xxx[nx]+rez[1]/2,
                           -lim-yyy[ny]-rez[2]/2,
                           lim+yyy[ny]+rez[2]/2)
    terra::values(out.prj)<-matrix(aperm(array(c(aperm(ttt,c(3,2,1)),
                                                 aperm(mm,c(3,2,1)),
                                                 aperm(bbb,c(3,2,1))),
                                               c(dims[1]+2*nx,dims[3],dims[2]+2*ny)),
                                         c(1,3,2)),
                                   ncol=dims[3])
    names(out.prj)<-names(x)
    terra::crs(out.prj)<-NULL
    terra::RGB(out.prj)<-terra::RGB(x)
    terra::units(out.prj)<-terra::units(x)

  }else{
    stop("Didn't recognize format of x--it should be a SpatVector or SpatRaster")
  }

  #get cropping polygons...
  if(prettify){
    # #constraints lifted from ricardo's code...
    # cols<-rainbow(6)
    # terra::plot(out.prj)
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
    # abline(v=0.447e7,h=-1.11e7,col=cols[6]) #intersects with line 3

    # #these should take care of any Alaskan weirdness...
    # #intersect with line 1...
    # abline(v=-8e6)
    # abline(h=11.9e6)
    # #horiztonal intersects with line 3, vertical with line 4
    # abline(h=-9.5e6)
    # abline(v=11.75e6)
    # #for central America...just modify line positioning a bit
    # # abline(v=0.295e7,h=-1.21e7,col=cols[3]) #intersect with line 3 --> change v to 0.293
    # abline(v=0.293e7)
    # # abline(v=0.312e7,h=-1.2e7,col=cols[4]) #intersect with line 3 --> change v to 0.3005, h to -1.1995
    # abline(v=0.3005e7)
    # # abline(v=0.4e7,h=-1.16e7,col=cols[5]) #intersects with line 3 --> change v to 0.39
    # abline(v=0.39e7)
    # # abline(v=0.45e7,-1.11e7,col=cols[6]) #intersects with line 3

    #need to alter first vertical line intersecting with line 3
    #0.45 --> 0.447
    #need to add corner around south central america for extra safety
    # abline(h=-12.53e6,v=2.825e6)

    #okay, last round I think...
    #there were some sneaky artifacts hiding about...
    # tmp.map<-terra::crop(out.prj,
    #                      terra::vect(cbind(c(-1.5e7,-1.5e7,1.5e7,1.5e7),
    #                                        c(-1.5e7,1.5e7,1.5e7,-1.5e7)),
    #                                  "polygons"))
    # terra::plot(tmp.map,
    # xlim=c(-1.5e7,1.5e7),ylim=c(-1.5e7,1.5e7)) #overall
    # xlim=c(1e6,6e6),ylim=c(-1.3e7,-1e7),asp=NA) #bottom central america
    # xlim=c(2.5e6,3.5e6),ylim=c(-1.26e7,-1.15e7),asp=NA) #bottom central america--zoomed
    # xlim=c(-5.1e6,-4.5e6),ylim=c(1.15e7,1.21e7),asp=NA) #top bering
    # xlim=c(1.15e7,1.21e7),ylim=c(-5.1e6,-4.5e6),asp=NA) #bottom bering
    # xlim=c(-1.2e7,-1.1e7),ylim=c(2.5e6,3.5e6),asp=NA) #top central america

    # lines(cropper)

    #fix top bering strait artifact--intersect with line 1
    # abline(h=1.176e7,v=-4.909e6)
    # abline(h=1.1725e7,v=-4.55e6)
    #lines to make bottom bering strait safe--intersect with line 4
    # abline(h=-4.55e6,v=1.1935e7)
    #fixing top central america is gonna be a bit rough...
    #2 new hor/ver lines to intersect with line 2 hor/ver line 1 and hor/ver line 2
    # abline(v=-1.159e7,h=3.01e6)
    # abline(h=2.93e6,v=-1.121e7)
    #and modify hor/ver line 1; this directly intersects with prev. hor/ver line
    # abline(h=2.88e6,v=-1.113e7)





    #converting the above to a polygon...
    lim2<-1.1*lim
    cropper<-
      rbind(
        #left lim and line 6
        c(-lim2,2.3e7-lim2),
        #line 6 and line 1
        c((2.3e7-1.089e7)/(-0.176-1),1.089e7-0.176*(2.3e7-1.089e7)/(-0.176-1)),
        #line 1 and NEW vertical line 1
        c(-8e6,1.089e7-0.176*-8e6),
        #NEW vertical and horizontal line 1
        c(-8e6,11.9e6),
        #NEW horizontal line 1 and line 1
        c((1.089e7-11.9e6)/0.176,11.9e6),
        #line 1 and NEW horizontal line 4
        c((1.089e7-1.176e7)/0.176,1.176e7),
        #NEW horizontal and vertical line 4
        c(-4.909e6,1.176e7),
        #NEW vertical line 4 and line 1
        c(-4.909e6,1.089e7-0.176*-4.909e6),
        #line 1 and NEW horizontal line 5
        c((1.089e7-1.1725e7)/0.176,1.1725e7),
        #NEW horizontal and vertical line 5
        c(-4.55e6,1.1725e7),
        #NEW vertical line 5 and line 1
        c(-4.55e6,1.089e7-0.176*-4.55e6),
        #line 1 and line 5
        c((1e7-1.089e7)/(-0.176+0.5),1.089e7-0.176*(1e7-1.089e7)/(-0.176+0.5)),
        #line 5 and right lim
        c(lim2,1e7-0.5*lim2),
        #right lim and line 4
        c(lim2,-1.274e7/0.172+lim2/0.172),
        #line 4 and NEW horizontal line 6
        c(-4.55e6*0.172+1.274e7,-4.55e6),
        #NEW horizontal and vertical line 6
        c(1.1935e7,-4.55e6),
        #NEW vertical line 6 and line 4
        c(1.1935e7,-1.274e7/0.172+1.1935e7/0.172),
        #line 4 and NEW vertical line 2
        c(11.75e6,-1.274e7/0.172+11.75e6/0.172),
        #NEW vertical and horizontal line 2
        c(11.75e6,-9.5e6),
        #NEW horizontal line 2 and line 3
        c((1.274e7-9.5e6*0.172),-9.5e6),
        # #line 4 and line 3
        # c((-1.274e7/0.172+1.378e7)/(0.46-1/0.172),-1.378e7+0.46*(-1.274e7/0.172+1.378e7)/(0.46-1/0.172)),
        #line 3 and horizontal line 6
        c((-1.11e7+1.378e7)/0.46,-1.11e7),
        #horizontal and vertical line 6
        c(0.447e7,-1.11e7),
        # #vertical line 6 and line 3
        # c(0.447e7,-1.378e7+0.46*0.447e7),
        # #line 3 and horizontal line 5
        # c((-1.16e7+1.378e7)/0.46,-1.16e7),
        #vertical line 6 and horizontal line 5
        c(0.447e7,-1.16e7),
        #horizontal and vertical line 5
        c(0.39e7,-1.16e7),
        #vertical line 5 and line 3
        c(0.39e7,-1.378e7+0.46*0.39e7),
        #line 3 and horizontal line 4
        c((-1.2e7+1.378e7)/0.46,-1.1995e7),
        #horizontal and vertical line 4
        c(0.3005e7,-1.1995e7),
        # #vertical line 4 and line 3
        # c(0.312e7,-1.378e7+0.46*0.312e7),
        # #line 3 and horizontal line 3
        # c((-1.21e7+1.378e7)/0.46,-1.21e7),
        #vertical line 4 and horizontal line 3
        c(0.3005e7,-1.21e7),
        #horizontal and vertical line 3
        c(0.293e7,-1.21e7),
        #vertical line 3 and line 3
        c(0.293e7,-1.378e7+0.46*0.295e7),
        #line 3 and NEW vertical line 3
        c(2.825e6,-1.378e7+0.46*2.825e6),
        #NEW vertical and horizontal line 3
        c(2.825e6,-12.53e6),
        #NEW horizontal line 3 and line3
        c((-12.53e6+1.378e7)/0.46,-12.53e6),
        #line 3 and bottom lim
        c((-lim2+1.378e7)/0.46,-lim2),
        #bottom lim and line 2
        c(lim2*0.565-0.984e7,-lim2),
        #line 2 and vertical line 1
        c(-1.113e7,-0.984e7/0.565+1.113e7/0.565),
        #vertical and horizontal line 1
        c(-1.113e7,0.288e7),
        # #horizontal line 1 and line 2
        # c(-0.565*0.288e7-0.984e7,0.288e7),
        #horizontal line 1 and NEW vertical line 8
        c(-1.121e7,0.288e7),
        #NEW vertical and horizontal line 8
        c(-1.121e7,2.93e6),
        #NEW horizontal line 8 and line 2
        c(-0.565*2.93e6-0.984e7,2.93e6),
        #line 2 and NEW horizontal line 7
        c(-0.565*3.01e6-0.984e7,3.01e6),
        #NEW horizontal and vertical line 7
        c(-1.159e7,3.01e6),
        #NEW vertical line 7 and line 2
        c(-1.159e7,-0.984e7/0.565+1.159e7/0.565),
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
  if(is.na(geom.type)){
    out.prj<-terra::crop(out.prj,terra::vect(cropper,"polygons"),mask=TRUE)
  }else{
    out.prj<-terra::crop(out.prj,terra::vect(cropper,"polygons"))
  }


  #adding frame if prettify is TRUE...
  if(frame&prettify){
    wd<-amount*2*lim
    rectangle<-cbind(c(-lim-wd,
                       -lim-wd,
                       lim+wd,
                       lim+wd),
                     c(-lim-wd,
                       lim+wd,
                       lim+wd,
                       -lim-wd))
    if(is.na(geom.type)){
      out.prj<-terra::crop(out.prj,terra::vect(rectangle,"polygons"),extend=TRUE)
      foc.nas<-
        seq_len(terra::ncell(out.prj))[
          -terra::cells(out.prj,terra::vect(cropper,"polygons"))[,2]
          ]

      #idea! --> alternate na.policy to progressively blur?
      tmp<-terra::focal(out.prj,patch_width,fun=patch_method,na.policy="only")
      counter<-1
      while(any(is.na(terra::values(tmp)[foc.nas,]))&counter<max_patch_iter){
        patch_width<-patch_width+patch_width_inc
        tmp<-terra::focal(tmp,patch_width,fun=patch_method,na.policy="only")
        counter<-counter+1
      }
      terra::values(out.prj)[foc.nas,]<-terra::values(tmp)[foc.nas,]
      terra::RGB(out.prj)<-terra::RGB(x)
      terra::units(out.prj)<-terra::units(x)
    }else{
      rect1<-terra::erase(terra::vect(rectangle,"polygons"),
                          terra::vect(cropper,"polygons"))
      rectangle<-terra::buffer(rect1,
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
        bl<-lt<-tr<-rb<-TRUE
        for(i in inds){
          tmp<-terra::union(out.prj[i],rectangle)
          tmp<-terra::geom(tmp[length(tmp)])
          coords<-terra::vect(tmp[,c("x","y"),drop=FALSE],"points")

          #border detection not perfect...buffering helps, but doesn't seem to get all cases...
          on.border<-terra::is.related(terra::buffer(coords,border_tol),
                                       terra::vect(rbind(cropper,cropper[1,]),"lines"),
                                       "intersects")

          # terra::plot(coords)
          # terra::plot(terra::vect(cropper,"lines"),add=TRUE)
          # points(coords[on.border,],col="red")

          if(any(on.border)){
            foc.coords<-tmp[on.border,c("x","y"),drop=FALSE]
            angs<-atan2(foc.coords[,2],foc.coords[,1])

            slices<-findInterval(angs,c(-pi,-3*pi/4,-pi/4,pi/4,3*pi/4,pi))
            #1,5-->left lim
            #2-->bottom lim
            #3-->right lim
            #4-->top lim

            angs<-(foc.coords[,2]-cent[2])/(foc.coords[,1]-cent[1])

            for(j in unique(slices)){
              tmp.inds<-slices==j
              # y=angs*x-angs*-3.75e6+2e6
              # y=angs*(x+3.75e6)+2e6
              foc.coords[tmp.inds,]<-switch(j,
                                            cbind(-lim-wd,
                                                  angs[tmp.inds]*(-lim-wd-cent[1])+cent[2]),
                                            cbind((-lim-wd-cent[2])/angs[tmp.inds]+cent[1],
                                                  -lim-wd),
                                            cbind(lim+wd,
                                                  angs[tmp.inds]*(lim+wd-cent[1])+cent[2]),
                                            cbind((lim+wd-cent[2])/angs[tmp.inds]+cent[1],
                                                  lim+wd),
                                            cbind(-lim-wd,
                                                  angs[tmp.inds]*(-lim-wd-cent[1])+cent[2]))
            }
            tmp[on.border,c("x","y")]<-foc.coords

            #insert corners...
            #technically don't need to check for corners after they're found first time...
            #just a way to optimize later on
            #also would probably be best to do a loop
            if(bl|lt|tr|rb){
              b.checks<-tmp[,c("y","x","y","x"),drop=FALSE]==do.call(cbind,
                                                                     lapply(c(-lim-wd,
                                                                              -lim-wd,
                                                                              lim+wd,
                                                                              lim+wd),
                                                                            rep,each=nrow(tmp)))
              colnames(b.checks)<-c("bb","ll","tt","rr")
              tmp<-cbind(tmp,
                         b.checks)
              for(j in unique(tmp[,"part"])){
                if(bl|lt|tr|rb){
                  p.inds.tf<-tmp[,"part"]==j
                  p.inds<-which(p.inds.tf)
                  nn<-length(p.inds)
                  if(bl){
                    if(any(tmp[p.inds,"bb"]&tmp[p.inds,"ll"])&bl){
                      bl<-FALSE
                    }else{
                      corner.check<-(tmp[p.inds[c(2:nn,1)],"bb"]&tmp[p.inds,"ll"])|
                        (tmp[p.inds[c(2:nn,1)],"ll"]&tmp[p.inds,"bb"])
                      if(any(corner.check)){
                        if(sum(corner.check)>1) stop()
                        tmp.inds<-p.inds[corner.check] #should always be length 1
                        tmp.tmp<-tmp[tmp.inds,,drop=FALSE]
                        tmp.tmp[,c("x","y","bb","ll","tt","rr")]<-c(-lim-wd,-lim-wd,1,1,0,0)
                        tmp<-rbind(tmp[seq_len(tmp.inds),,drop=FALSE],
                                   tmp.tmp,
                                   tmp[-seq_len(tmp.inds),,drop=FALSE])
                        p.inds.tf<-append(p.inds.tf,TRUE,tmp.inds)
                        p.inds<-which(p.inds.tf)
                        nn<-nn+1
                        bl<-FALSE
                      }
                    }
                  }
                  if(lt){
                    if(any(tmp[p.inds,"ll"]&tmp[p.inds,"tt"])){
                      lt<-FALSE
                    }else{
                      corner.check<-(tmp[p.inds[c(2:nn,1)],"ll"]&tmp[p.inds,"tt"])|
                        (tmp[p.inds[c(2:nn,1)],"tt"]&tmp[p.inds,"ll"])
                      if(any(corner.check)){
                        if(sum(corner.check)>1) stop()
                        tmp.inds<-p.inds[corner.check] #should always be length 1
                        tmp.tmp<-tmp[tmp.inds,,drop=FALSE]
                        tmp.tmp[,c("x","y","bb","ll","tt","rr")]<-c(-lim-wd,lim+wd,0,1,1,0)
                        tmp<-rbind(tmp[seq_len(tmp.inds),,drop=FALSE],
                                   tmp.tmp,
                                   tmp[-seq_len(tmp.inds),,drop=FALSE])
                        p.inds.tf<-append(p.inds.tf,TRUE,tmp.inds)
                        p.inds<-which(p.inds.tf)
                        nn<-nn+1
                        lt<-FALSE
                      }
                    }
                  }
                  if(tr){
                    if(any(tmp[p.inds,"tt"]&tmp[p.inds,"rr"])){
                      tr<-FALSE
                    }else{
                      corner.check<-(tmp[p.inds[c(2:nn,1)],"tt"]&tmp[p.inds,"rr"])|
                        (tmp[p.inds[c(2:nn,1)],"rr"]&tmp[p.inds,"tt"])
                      if(any(corner.check)){
                        if(sum(corner.check)>1) stop()
                        tmp.inds<-p.inds[corner.check] #should always be length 1
                        tmp.tmp<-tmp[tmp.inds,,drop=FALSE]
                        tmp.tmp[,c("x","y","bb","ll","tt","rr")]<-c(lim+wd,lim+wd,0,0,1,1)
                        tmp<-rbind(tmp[seq_len(tmp.inds),,drop=FALSE],
                                   tmp.tmp,
                                   tmp[-seq_len(tmp.inds),,drop=FALSE])
                        p.inds.tf<-append(p.inds.tf,TRUE,tmp.inds)
                        p.inds<-which(p.inds.tf)
                        nn<-nn+1
                        tr<-FALSE
                      }
                    }
                  }
                  if(rb){
                    if(any(tmp[p.inds,"rr"]&tmp[p.inds,"bb"])){
                      rb<-FALSE
                    }else{
                      corner.check<-(tmp[p.inds[c(2:nn,1)],"rr"]&tmp[p.inds,"bb"])|
                        (tmp[p.inds[c(2:nn,1)],"bb"]&tmp[p.inds,"rr"])
                      if(any(corner.check)){
                        if(sum(corner.check)>1) stop()
                        tmp.inds<-p.inds[corner.check] #should always be length 1
                        tmp.tmp<-tmp[tmp.inds,,drop=FALSE]
                        tmp.tmp[,c("x","y","bb","ll","tt","rr")]<-c(lim+wd,-lim-wd,1,0,0,1)
                        tmp<-rbind(tmp[seq_len(tmp.inds),,drop=FALSE],
                                   tmp.tmp,
                                   tmp[-seq_len(tmp.inds),,drop=FALSE])
                        p.inds.tf<-append(p.inds.tf,TRUE,tmp.inds)
                        p.inds<-which(p.inds.tf)
                        nn<-nn+1
                        rb<-FALSE
                      }
                    }
                  }
                }
              }
              tmp<-tmp[,1:5,drop=FALSE]
            }
            tmp<-terra::aggregate(terra::union(out.prj[i],terra::buffer(terra::vect(tmp,"polygons"),0)))
          }else{
            tmp<-out.prj[i]
          }
          terra::values(tmp)<-terra::values(out.prj)[i,,drop=FALSE]
          out.prj<-rbind(out.prj,tmp)

        }
        out.prj<-out.prj[-inds]

        #deal with negative space
        neg.space<-terra::union(rect1,terra::aggregate(out.prj))[1]
        if(length(out.prj)==1&merge_neg_space){
          tmp<-terra::aggregate(terra::union(out.prj,terra::buffer(neg.space,neg_space_tol)))
          terra::values(tmp)<-terra::values(out.prj)
          out.prj<-tmp
        }else{
          vals<-terra::values(out.prj)[length(out.prj)+1,]
          vals[["TEMPORARY_ORDER"]]<-length(out.prj)+1
          terra::values(neg.space)<-vals
          out.prj<-rbind(out.prj,neg.space)
        }
      }

      out.prj<-terra::sort(out.prj,v="TEMPORARY_ORDER")
      vals<-terra::values(out.prj)
      vals[["TEMPORARY_ORDER"]]<-NULL
      terra::values(out.prj)<-vals
    }
  }

  if(!is.na(geom.type)&revalidate_geoms){
    .revalidate.geoms(out.prj)
  }else{
    out.prj
  }
}

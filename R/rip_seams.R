#' @export
rip_spilhaus_seam<-function(x,
                            NA_tol=1.1e6,
                            SA_tol=1.1e6,
                            Asia_tol=1.1e6,
                            seam_res=1000){

  #argument warnings
  if(seam_res<600){
    warning("seam resolutions below 600 usually result in reprojection errors",
            " in converting longtidue-latitude to Spilhaus coordinates")
  }
  if(NA_tol<1e6|SA_tol<1e6|Asia_tol<1e6){
    warning("corner tolerances below 1 million (1e6) usually result in reprojection errors")
  }

  #max coordinate in spilhaus projection
  lim<-11825474

  #convert spilhaus border to lonlat
  lim.seq<-seq(-lim,lim,length.out=seam_res)
  seam.xy<-rbind(cbind(-lim,lim.seq),
                 cbind(lim.seq,lim),
                 cbind(lim,-lim.seq),
                 cbind(-lim.seq,-lim))
  seam.lonlat<-.spilhaus2lonlat(seam.xy[,1],seam.xy[,2])

  #eliminate failed conversions
  seam.lonlat<-seam.lonlat[complete.cases(seam.lonlat),]

  #make corner patches
  NA.patch<-.spilhaus2lonlat(c(seq(-lim,-lim+NA_tol,length.out=100),
                               seq(lim,lim-NA_tol,length.out=100)),
                             c(seq(lim-NA_tol,lim,length.out=100),
                               seq(-lim+NA_tol,-lim,length.out=100)))
  SA.patch<-.spilhaus2lonlat(seq(-lim,-lim+SA_tol,length.out=100),
                             seq(-lim+SA_tol,-lim,length.out=100))
  Asia.patch<-.spilhaus2lonlat(seq(lim-Asia_tol,lim,length.out=100),
                               seq(lim,lim-Asia_tol,length.out=100))

  #rip along IDL
  seam.lonlat<-lapply(split(seam.lonlat,seam.lonlat[,1]>0),matrix,ncol=2)

  ##the following section simply extends each seam slightly past IDL##
  #find westernmost coord pair in west hemisphere seam
  nn<-nrow(seam.lonlat[[1]])
  min1<-which.min((seam.lonlat[[1]][,1]+seam.lonlat[[1]][c(nn,1:(nn-1)),1])/2)
  min1<-c(min1-2,min1-1)%%nn+1
  #find easternmost coord pair in east hemisphere seam
  nn<-nrow(seam.lonlat[[2]])
  max2<-which.max((seam.lonlat[[2]][,1]+seam.lonlat[[2]][c(nn,1:(nn-1)),1])/2)
  max2<-c(max2-2,max2-1)%%nn+1
  #graft easternmost coord pair of east hemisphere seam to...
  #...western edge of west hemisphere seam
  tmp.xx<-seam.lonlat[[1]][,1]
  tmp.yy<-seam.lonlat[[1]][,2]
  if(min1[1]<min1[2]){
    tmp.xx<-append(tmp.xx,seam.lonlat[[2]][max2[2],1]-360,min1[1])
    tmp.yy<-append(tmp.yy,seam.lonlat[[2]][max2[2],2],min1[1])
    tmp.xx<-append(tmp.xx,seam.lonlat[[2]][max2[1],1]-360,min1[2])
    tmp.yy<-append(tmp.yy,seam.lonlat[[2]][max2[1],2],min1[2])
  }else{
    tmp.xx<-append(tmp.xx,seam.lonlat[[2]][max2[1],1]-360,min1[2]-1)
    tmp.yy<-append(tmp.yy,seam.lonlat[[2]][max2[1],2],min1[2]-1)
    tmp.xx<-append(tmp.xx,seam.lonlat[[2]][max2[2],1]-360,min1[1]+1)
    tmp.yy<-append(tmp.yy,seam.lonlat[[2]][max2[2],2],min1[1]+1)
  }
  tmp.seam<-cbind(tmp.xx,tmp.yy)
  #graft westernmost coord pair of west hemisphere seam to...
  #...eastern edge of east hemisphere seam
  tmp.xx<-seam.lonlat[[2]][,1]
  tmp.yy<-seam.lonlat[[2]][,2]
  if(max2[1]<max2[2]){
    tmp.xx<-append(tmp.xx,seam.lonlat[[1]][min1[2],1]+360,max2[1])
    tmp.yy<-append(tmp.yy,seam.lonlat[[1]][min1[2],2],max2[1])
    tmp.xx<-append(tmp.xx,seam.lonlat[[1]][min1[1],1]+360,max2[2])
    tmp.yy<-append(tmp.yy,seam.lonlat[[1]][min1[1],2],max2[2])
  }else{
    tmp.xx<-append(tmp.xx,seam.lonlat[[1]][min1[1],1]+360,max2[2]-1)
    tmp.yy<-append(tmp.yy,seam.lonlat[[1]][min1[1],2],max2[2]-1)
    tmp.xx<-append(tmp.xx,seam.lonlat[[1]][min1[2],1]+360,max2[1]+1)
    tmp.yy<-append(tmp.yy,seam.lonlat[[1]][min1[2],2],max2[1]+1)
  }
  seam.lonlat[[1]]<-tmp.seam
  seam.lonlat[[2]]<-cbind(tmp.xx,tmp.yy)

  #convert seams to polygons, merge with corner patches, use to erase x
  #buffer with 0 to prevent invalid geometry
  seam<-terra::buffer(terra::vect(seam.lonlat,
                                  "polygons"),
                      width=0)
  seam<-rbind(terra::aggregate(terra::union(seam[1],
                                            terra::vect(list(NA.patch,SA.patch),
                                                        "polygons"))),
              terra::aggregate(terra::union(seam[2],
                                            terra::vect(Asia.patch,
                                                        "polygons"))))
  terra::erase(x,seam)

}

#' @export
rip_IDL_seam<-function(x,
                       seam_width=0.005,
                       seam_res=2000){
  if(seam_width<0.001){
    warning("seam widths below 0.001 usually result in reprojection errors",
            " in converting from Spilhaus to longitude-latitude coordinates")
  }
  if(seam_res<1200){
    warning("seam resolutions below 1200 usually result in reprojection errors",
            " in converting from Spilhaus to longitude-latitude coordinates")
  }

  #max coordinate in spilhaus projection
  lim<-11825474

  seam.lonlat<-cbind(180-seam_width/2,seq(-90,90,length.out=seam_res))
  seam.lonlat<-rbind(seam.lonlat,-seam.lonlat)

  seam.xy<-rSpilhaus:::.lonlat2spilhaus(seam.lonlat[,1],seam.lonlat[,2])

  #doesn't seem like failed conversion are a problem here...

  #convenient line to rip seam
  seam.xy<-lapply(split(seam.xy,seam.xy[,2]>0),matrix,ncol=2)

  ##the following section simply extends each seam slightly past spilhaus borders##
  #find rightmost coord pair in lower seam
  nn<-nrow(seam.xy[[1]])
  max1<-which.max((seam.xy[[1]][,1]+seam.xy[[1]][c(nn,1:(nn-1)),1])/2)
  max1<-c(max1-2,max1-1)%%nn+1
  #find uppermost coord pair in upper seam
  nn<-nrow(seam.xy[[2]])
  max2<-which.max((seam.xy[[2]][,2]+seam.xy[[2]][c(nn,1:(nn-1)),2])/2)
  max2<-c(max2-2,max2-1)%%nn+1
  #graft rightmost coord pair of lower seam to...
  #...upper edge of upper seam
  tmp.xx<-seam.xy[[1]][,1]
  tmp.yy<-seam.xy[[1]][,2]
  if(max1[1]<max1[2]){
    tmp.xx<-append(tmp.xx,
                   cos(pi/2)*(seam.xy[[2]][max2[2],1]-lim)-
                     sin(pi/2)*(seam.xy[[2]][max2[2],2]-lim)+
                     lim,
                   max1[1])
    tmp.yy<-append(tmp.yy,
                   sin(pi/2)*(seam.xy[[2]][max2[2],1]-lim)+
                     cos(pi/2)*(seam.xy[[2]][max2[2],2]-lim)+
                     lim,
                   max1[1])
    tmp.xx<-append(tmp.xx,
                   cos(pi/2)*(seam.xy[[2]][max2[1],1]-lim)-
                     sin(pi/2)*(seam.xy[[2]][max2[1],2]-lim)+
                     lim,
                   max1[2])
    tmp.yy<-append(tmp.yy,
                   sin(pi/2)*(seam.xy[[2]][max2[1],1]-lim)+
                     cos(pi/2)*(seam.xy[[2]][max2[1],2]-lim)+
                     lim,
                   max1[2])
  }else{
    tmp.xx<-append(tmp.xx,
                   cos(pi/2)*(seam.xy[[2]][max2[1],1]-lim)-
                     sin(pi/2)*(seam.xy[[2]][max2[1],2]-lim)+
                     lim,
                   max1[2]-1)
    tmp.yy<-append(tmp.yy,
                   sin(pi/2)*(seam.xy[[2]][max2[1],1]-lim)+
                     cos(pi/2)*(seam.xy[[2]][max2[1],2]-lim)+
                     lim,
                   max1[2]-1)
    tmp.xx<-append(tmp.xx,
                   cos(pi/2)*(seam.xy[[2]][max2[2],1]-lim)-
                     sin(pi/2)*(seam.xy[[2]][max2[2],2]-lim)+
                     lim,
                   max1[1]+1)
    tmp.yy<-append(tmp.yy,
                   sin(pi/2)*(seam.xy[[2]][max2[2],1]-lim)+
                     cos(pi/2)*(seam.xy[[2]][max2[2],2]-lim)+
                     lim,
                   max1[1]+1)
  }
  tmp.seam<-cbind(tmp.xx,tmp.yy)
  #graft westernmost coord pair of west hemisphere seam to...
  #...eastern edge of east hemisphere seam
  tmp.xx<-seam.xy[[2]][,1]
  tmp.yy<-seam.xy[[2]][,2]

  if(max2[1]<max2[2]){
    tmp.xx<-append(tmp.xx,
                   cos(-pi/2)*(seam.xy[[1]][max1[2],1]-lim)-
                     sin(-pi/2)*(seam.xy[[1]][max1[2],2]-lim)+
                     lim,
                   max2[1])
    tmp.yy<-append(tmp.yy,
                   sin(-pi/2)*(seam.xy[[1]][max1[2],1]-lim)+
                     cos(-pi/2)*(seam.xy[[1]][max1[2],2]-lim)+
                     lim,
                   max2[1])
    tmp.xx<-append(tmp.xx,
                   cos(-pi/2)*(seam.xy[[1]][max1[1],1]-lim)-
                     sin(-pi/2)*(seam.xy[[1]][max1[1],2]-lim)+
                     lim,
                   max2[2])
    tmp.yy<-append(tmp.yy,
                   sin(-pi/2)*(seam.xy[[1]][max1[1],1]-lim)+
                     cos(-pi/2)*(seam.xy[[1]][max1[1],2]-lim)+
                     lim,
                   max2[2])
  }else{
    tmp.xx<-append(tmp.xx,
                   cos(-pi/2)*(seam.xy[[1]][max1[1],1]-lim)-
                     sin(-pi/2)*(seam.xy[[1]][max1[1],2]-lim)+
                     lim,
                   max2[2]-1)
    tmp.yy<-append(tmp.yy,
                   sin(-pi/2)*(seam.xy[[1]][max1[1],1]-lim)+
                     cos(-pi/2)*(seam.xy[[1]][max1[1],2]-lim)+
                     lim,
                   max2[2]-1)
    tmp.xx<-append(tmp.xx,
                   cos(-pi/2)*(seam.xy[[1]][max1[2],1]-lim)-
                     sin(-pi/2)*(seam.xy[[1]][max1[2],2]-lim)+
                     lim,
                   max2[1]+1)
    tmp.yy<-append(tmp.yy,
                   sin(-pi/2)*(seam.xy[[1]][max1[2],1]-lim)+
                     cos(-pi/2)*(seam.xy[[1]][max1[2],2]-lim)+
                     lim,
                   max2[1]+1)
  }
  seam.xy[[1]]<-tmp.seam
  seam.xy[[2]]<-cbind(tmp.xx,tmp.yy)

  #convert seams to polygons, use to erase x
  #doesn't look like invalid geometry is an issue here
  seam<-terra::vect(seam.xy,
                    "polygons")
  terra::erase(x,seam)

}

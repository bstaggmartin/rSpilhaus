.insert.cols<-function(x,insert,holder,def.nm){
  check<-try(.check.cols(x,holder),silent=TRUE)
  if(inherits(check,"try-error")){
    colnames(insert)<-if(is.character(x)) x else def.nm
    holder<-cbind(holder,insert)
  }else{
    holder[,check]<-insert
  }
  holder
}

.check.cols<-function(x,dat,nice.nm=NULL){
  in.x<-x
  if(is.character(x)){
    x<-match(x,colnames(dat))
  }
  if(is.na(x)|x<1|x>ncol(dat)){
    stop("could not find column for ",nice.nm," (",in.x,")")
  }
  x
}

#possibly hacky fix because running into trouble reprojecting anything RIGHT at border...
#figured out the issue--buffering was causing things to expand past border
#I don't think this is necessary anymore, but I'll keep it for extra safety...
#actually, still seems necessary because of NA border...
.clamp.borders<-function(x,lims=c(-11825474,11825474)){
  x[x<lims[1]]<-lims[1]
  x[x>lims[2]]<-lims[2]
  x
}

#' @export
spilhaus2lonlat<-function(dat,...){
  UseMethod("spilhaus2lonlat")
}

# dat<-rnaturalearth::ne_load(scale=50,type="NE2_50M_SR",category="raster",destdir="misc_R_tests/rasters")
# x_res<-500;y_res<-500;method<-"simple"
#' @export
#' @method spilhaus2lonlat SpatRaster
spilhaus2lonlat.SpatRaster<-function(dat,
                                     lon_res=720,
                                     lat_res=360,
                                     sample_method="simple",
                                     border_tol=5e4, #seems to be a robust default...
                                     patch_method="modal",
                                     patch_width=3){

  xx<-seq(-180+180/lon_res,180-180/lon_res,length.out=lon_res)
  yy<-seq(-90+90/lat_res,90-90/lat_res,length.out=lat_res)

  tmp<-expand.grid("lon"=xx,"lat"=rev(yy))
  tmp<-lonlat2spilhaus(tmp,lon="lon",lat="lat",x="x",y="y")

  #hmmm...cells only works if method="simple"
  #need another test for bilinear...
  vals<-terra::extract(dat,as.matrix(tmp[,c("x","y")]),
                       method=sample_method)
                       #,cells=TRUE)
  # if(!is.null(vals$cell)){
  #   foc.nas<-is.na(vals$cell)
  #   vals<-vals[-1]
  # }else{
  #   lim<-11825474
  #   ltmp<-tmp[,c("x","y")]<=-lim+border_tol
  #   htmp<-tmp[,c("x","y")]>=lim-border_tol
  #   foc.nas<-ltmp[,1]|ltmp[,2]|htmp[,1]|htmp[,2]
  # }

  out<-dat
  terra::nrow(out)<-lat_res
  terra::ncol(out)<-lon_res
  terra::ext(out)<-c(-180,180,
                     -90,90)
  terra::values(out)<-matrix(unlist(vals,use.names=FALSE),
                             ncol=ncol(vals),
                             dimnames=list(NULL,colnames(vals)))
  names(out)<-colnames(vals)
  terra::crs(out)<-NULL

  #fill in seam...
  #the seam turned out to be an error resulting from incorrect extent specification!
  #filling now unnecessary
  # tmp<-terra::focal(out,patch_width,fun=patch_method,na.policy="only")
  # while(any(is.na(terra::values(tmp)[foc.nas,]))){
  #   tmp<-terra::focal(tmp,patch_width,fun=patch_method,na.policy="only")
  # }
  # terra::values(out)[foc.nas,]<-terra::values(tmp)[foc.nas,]

  #hmmm...seams seem to present a problem
  #perhaps extend borders slightly?

  #looks like these should be done after the last step...
  terra::RGB(out)<-terra::RGB(dat)
  terra::units(out)<-terra::units(dat)

  #might have accidentally discarded some important attributes
  #but this seems like a good start

  out
}

#' @export
#' @method spilhaus2lonlat SpatVector
spilhaus2lonlat.SpatVector<-function(dat,
                                     xlim=c(-11825474,11825474),
                                     ylim=xlim,
                                     rip=TRUE,
                                     rip_seam_res=2000,
                                     rip_seam_width=0.005,
                                     dissolve=TRUE,
                                     dissolve_seam_res=1000,
                                     dissolve_seam_width=2*rip_seam_width,
                                     revalidate_geoms=TRUE){
  if(rip){
    dat<-rip_IDL_seam(dat,
                      seam_width=rip_seam_width,
                      seam_res=rip_seam_res)
  }

  tmp<-spilhaus2lonlat(terra::geom(dat),
                       lon="x",lat="y",
                       x="x",y="y",
                       xlim=xlim,ylim=ylim)
  dat<-terra::vect(tmp,
                   type=terra::geomtype(dat),
                   atts=terra::values(dat))

  if(is.null(revalidate_geoms)) revalidate_geoms<-TRUE
  revalidate_geoms[is.na(revalidate_geoms)]<-TRUE
  revalidate_geoms<-rep(revalidate_geoms,length.out=2)
  if(revalidate_geoms[1]){
    dat<-.revalidate.geoms(dat)
  }

  if(dissolve){
    dat<-dissolve_spilhaus_seam(dat,
                                seam_res=dissolve_seam_res,
                                seam_width=dissolve_seam_width,
                                revalidate_geoms=revalidate_geoms[2])
  }

  dat
}

#' @export
#' @method spilhaus2lonlat data.frame
spilhaus2lonlat.data.frame<-function(dat,
                                     x="x",y="y",
                                     lon="lon",lat="lat",
                                     xlim=c(-11825474,11825474),
                                     ylim=xlim){
  x<-.check.cols(x,dat,"Spilhaus x coordinates")
  y<-.check.cols(y,dat,"Spilhaus y coordinates")

  tmp<-.spilhaus2lonlat(.clamp.borders(dat[[x]],xlim),
                        .clamp.borders(dat[[y]],ylim))

  dat[[lon]]<-tmp[,1]
  dat[[lat]]<-tmp[,2]
  dat
}

#' @export
#' @method spilhaus2lonlat matrix
spilhaus2lonlat.matrix<-function(dat,
                                 x="x",y="y",
                                 lon="lon",lat="lat",
                                 xlim=c(-11825474,11825474),
                                 ylim=xlim){
  if(!is.numeric(dat)){
    stop("dat must be a numeric matrix")
  }
  x<-.check.cols(x,dat,"Spilhaus x coordinates")
  y<-.check.cols(y,dat,"Spilhaus y coordinates")

  tmp<-.spilhaus2lonlat(.clamp.borders(dat[,x],xlim),
                        .clamp.borders(dat[,y],ylim))

  dat<-.insert.cols(lon,tmp[,1,drop=FALSE],dat,"lon")
  .insert.cols(lat,tmp[,2,drop=FALSE],dat,"lat")
}

.spilhaus2lonlat<-function(spilhaus_x, spilhaus_y) {
  # constants
  e = sqrt(0.00669438)
  lat_center_deg = -49.56371678
  lon_center_deg = 66.94970198
  azimuth_deg = 40.17823482

  # parameters derived from constants
  lat_center_rad = lat_center_deg * pi / 180
  lon_center_rad = lon_center_deg * pi / 180
  azimuth_rad = azimuth_deg * pi / 180
  conformal_lat_center = -pi / 2 + 2 * atan(
    tan(pi/4 + lat_center_rad/2) *
      ((1 - e * sin(lat_center_rad)) / (1 + e * sin(lat_center_rad))) ^ (e / 2)
  )
  alpha = -asin(cos(conformal_lat_center) * cos(azimuth_rad))
  lambda_0 = lon_center_rad + atan2(tan(azimuth_rad), -sin(conformal_lat_center))
  beta = pi + atan2(-sin(azimuth_rad), -tan(conformal_lat_center))

  adams_x = (spilhaus_y - spilhaus_x) * sqrt(2) / 2
  adams_y = - (spilhaus_y + spilhaus_x) * sqrt(2) / 2

  adams_ws2 = "+proj=adams_ws2 +no_defs +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
  projection_fun = function(x, y) {
    tryCatch(
      sf::sf_project(from=adams_ws2, to=sf::st_crs(4326), pts=c(x, y)),
      error=function(e) c(NA, NA)
    )
  }
  projected = sf::sf_project(from=adams_ws2, to=sf::st_crs(4326), pts=cbind(adams_x, adams_y), keep = TRUE, warn = FALSE)
  lon_s = projected[,1]
  lat_s = projected[,2]

  #transformed coords in radians
  lon_s_rad = lon_s * pi / 180
  lat_s_rad = lat_s * pi / 180

  # conformal latitude
  lat_c = asin(sin(alpha) * sin(lat_s_rad) + cos(alpha) * cos(lat_s_rad) * cos(lon_s_rad - beta))

  # longitude, in radians
  lon = lambda_0 + atan2(
    cos(lat_s_rad) * sin(lon_s_rad - beta),
    sin(alpha) * cos(lat_s_rad) * cos(lon_s_rad - beta) - cos(alpha) * sin(lat_s_rad)
  )

  # latitude (iterative formula from https://mathworld.wolfram.com/ConformalLatitude.html)
  lat = lat_c
  for (i in 0:9) {
    lat = -0.5 * pi + 2 * atan(
      tan(pi / 4 + lat_c / 2) *
        ((1 + e * sin(lat)) / (1 - e * sin(lat))) ^ (e / 2)
    )
  }

  # coordinates in degrees
  longitude = ((lon * 180 / pi + 180) %% 360) - 180
  latitude = lat * 180 / pi

  return(cbind(longitude, latitude))
}

#' @export
lonlat2spilhaus<-function(dat,...){
  UseMethod("lonlat2spilhaus")
}

#' @export
#' @method lonlat2spilhaus SpatRaster
lonlat2spilhaus.SpatRaster<-function(dat,
                                     x_res=500,
                                     y_res=500,
                                     sample_method="simple",
                                     patch_method="modal",
                                     patch_width=3,
                                     max_patch_iter=500){

  lim<-11825474
  xx<-seq(-lim+lim/x_res,lim-lim/x_res,length.out=x_res)
  yy<-seq(-lim+lim/y_res,lim-lim/y_res,length.out=y_res)

  tmp<-expand.grid("x"=xx,"y"=rev(yy))
  tmp<-spilhaus2lonlat(tmp,x="x",y="y",lon="lon",lat="lat")

  #use cell index to determine which NAs are reprojection vs. "true" NAs
  #hmmm...cells only works if method="simple"
  #need another test for bilinear...
  vals<-terra::extract(dat,as.matrix(tmp[,c("lon","lat")]),
                       method=sample_method,cells=TRUE)
  if(!is.null(vals$cell)){
    foc.nas<-is.na(vals$cell)
    vals<-vals[-1]
  }else{
    foc.nas<-is.na(tmp[,"lon"])|is.na(tmp[,"lat"])
  }


  out<-dat
  terra::nrow(out)<-y_res
  terra::ncol(out)<-x_res
  terra::ext(out)<-c(-lim,lim,
                     -lim,lim)
  terra::values(out)<-matrix(unlist(vals,use.names=FALSE),
                             ncol=ncol(vals),
                             dimnames=list(NULL,colnames(vals)))
  names(out)<-colnames(vals)
  terra::crs(out)<-NULL

  #fill in corners...
  tmp<-terra::focal(out,patch_width,fun=patch_method,na.policy="only")
  counter<-1
  while(any(is.na(terra::values(tmp)[foc.nas,]))&counter<max_patch_iter){
    tmp<-terra::focal(tmp,patch_width,fun=patch_method,na.policy="only")
    counter<-counter+1
  }
  terra::values(out)[foc.nas,]<-terra::values(tmp)[foc.nas,]

  #looks like these should be done after the last step...
  terra::RGB(out)<-terra::RGB(dat)
  terra::units(out)<-terra::units(dat)

  #might have accidentally discarded some important attributes
  #but this seems like a good start

  out
}

#' @export
#' @method lonlat2spilhaus SpatVector
lonlat2spilhaus.SpatVector<-function(dat,
                                     rip=TRUE,
                                     rip_seam_res=1000,
                                     NA_tol=1.1e6,
                                     SA_tol=1.1e6,
                                     Asia_tol=1.1e6,
                                     patch=TRUE,
                                     NA_patch_tol=NA_tol+15,
                                     SA_patch_tol=SA_tol+15,
                                     Asia_patch_tol=Asia_tol+15,
                                     dissolve=TRUE,
                                     dissolve_seam_res=1000,
                                     dissolve_seam_width=1000,
                                     revalidate_geoms=TRUE){
  if(rip){
    dat<-rip_spilhaus_seam(x=dat,
                           NA_tol=NA_tol,
                           SA_tol=SA_tol,
                           Asia_tol=Asia_tol,
                           seam_res=rip_seam_res)
  }

  tmp<-lonlat2spilhaus(terra::geom(dat),
                       lon="x",lat="y",
                       x="x",y="y")
  dat<-terra::vect(tmp,
                   type=terra::geomtype(dat),
                   atts=terra::values(dat))

  if(is.null(revalidate_geoms)) revalidate_geoms<-TRUE
  revalidate_geoms[is.na(revalidate_geoms)]<-TRUE
  revalidate_geoms<-rep(revalidate_geoms,length.out=3)
  if(revalidate_geoms[1]){
    dat<-.revalidate.geoms(dat)
  }

  if(patch){
    dat<-patch_corners(dat,
                       NA_tol=NA_patch_tol,
                       SA_tol=SA_patch_tol,
                       Asia_tol=Asia_patch_tol,
                       revalidate_geoms=revalidate_geoms[2])
  }

  if(dissolve){
    dat<-dissolve_IDL_seam(dat,
                           seam_res=dissolve_seam_res,
                           seam_width=dissolve_seam_width,
                           revalidate_geoms=revalidate_geoms[3])
  }

  dat
}

#' @export
#' @method lonlat2spilhaus data.frame
lonlat2spilhaus.data.frame<-function(dat,
                                     lon="lon",lat="lat",
                                     x="x",y="y"){
  lon<-.check.cols(lon,dat,"latitude coordinates")
  lat<-.check.cols(lat,dat,"longitude coordinates")

  tmp<-.lonlat2spilhaus(dat[[lon]],dat[[lat]])

  dat[[x]]<-tmp[,1]
  dat[[y]]<-tmp[,2]
  dat
}

#' @export
#' @method lonlat2spilhaus matrix
lonlat2spilhaus.matrix<-function(dat,
                                 lon="lon",lat="lat",
                                 x="x",y="y"){
  if(!is.numeric(dat)){
    stop("dat must be a numeric matrix")
  }
  lon<-.check.cols(lon,dat,"latitude coordinates")
  lat<-.check.cols(lat,dat,"longitude coordinates")

  tmp<-.lonlat2spilhaus(dat[,lon],dat[,lat])

  dat<-.insert.cols(x,tmp[,1,drop=FALSE],dat,"x")
  .insert.cols(y,tmp[,2,drop=FALSE],dat,"y")
}

.lonlat2spilhaus<-function(longitude, latitude){

  # constants (https://github.com/OSGeo/PROJ/issues/1851)
  e = sqrt(0.00669438)
  lat_center_deg = -49.56371678
  lon_center_deg = 66.94970198
  azimuth_deg = 40.17823482

  # parameters derived from constants
  lat_center_rad = lat_center_deg * pi / 180
  lon_center_rad = lon_center_deg * pi / 180
  azimuth_rad = azimuth_deg * pi / 180
  conformal_lat_center = -pi / 2 + 2 * atan(
    tan(pi/4 + lat_center_rad/2) *
      ((1 - e * sin(lat_center_rad)) / (1 + e * sin(lat_center_rad))) ^ (e / 2)
  )
  alpha = -asin(cos(conformal_lat_center) * cos(azimuth_rad))
  lambda_0 = lon_center_rad + atan2(tan(azimuth_rad), -sin(conformal_lat_center))
  beta = pi + atan2(-sin(azimuth_rad), -tan(conformal_lat_center))

  # coordinates in radians
  lon = longitude * pi / 180
  lat = latitude * pi / 180

  # conformal latitude, in radians
  lat_c = -pi / 2 + 2 * atan(
    tan(pi/4 + lat/2) * ((1 - e * sin(lat)) / (1 + e * sin(lat))) ^ (e / 2)
  )

  # transformed lat and lon, in degrees
  lat_s = 180 / pi * asin(sin(alpha) * sin(lat_c) - cos(alpha) * cos(lat_c) * cos(lon - lambda_0))
  lon_s = 180 / pi * (
    beta + atan2(
      cos(lat_c) * sin(lon - lambda_0),
      (sin(alpha) * cos(lat_c) * cos(lon - lambda_0) + cos(alpha) * sin(lat_c))
    )
  )

  # projects transformed coordinates onto plane (Adams World in a Square II)
  adams_ws2 = "+proj=adams_ws2 +no_defs +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
  projected = sf::sf_project(from=sf::st_crs(4326), to=adams_ws2, pts=cbind(lon_s, lat_s))
  adams_x = projected[,1]
  adams_y = projected[,2]
  spilhaus_x = -(adams_x + adams_y) / sqrt(2)
  spilhaus_y = (adams_x - adams_y) / sqrt(2)

  return(cbind(spilhaus_x, spilhaus_y)) #, adams_x, adams_y, lon_s, lat_s))
}

#' @export
make_graticules<-function(lon=30,lat=30,
                          lon_start=0,lat_start=0,
                          lon_res=100,lat_res=100){
  lons<-c(rev(seq(lon_start,-180,by=-lon)),
          seq(lon_start,180,by=lon)[-1])
  lats<-c(rev(seq(lat_start,-90,by=-lat)),
          seq(lat_start,90,by=lat)[-1])
  fine.lons<-seq(-180,180,length.out=lon_res)
  fine.lats<-seq(-90,90,length.out=lat_res)
  lon.lines<-lapply(lons,function(ii) cbind(ii,fine.lats))
  lat.lines<-lapply(lats,function(ii) cbind(fine.lons,ii))
  out<-terra::vect(c(lon.lines,lat.lines),"lines")
  terra::values(out)<-data.frame("GRATICULE_TYPE"=rep(c("lon","lat"),
                                                      c(length(lon.lines),
                                                        length(lat.lines))),
                                 "GRATICULE_DEGREE"=c(lons,lats))
  out
}

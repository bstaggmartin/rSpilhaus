library(rSpilhaus)

gen_graticules<-function(lon=30,lat=30,
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
terra::plot(gen_graticules())

#figure out if these functions are "good enough"...
#so graticules and other line-based features don't need special treatment...
#nope, need special treatment for line corner patching and dissolution...
#not too big a deal, but will be annoying

#probably best approach is to simply make special patch_corners/dissolve functions for
#line data...
#dissolution is likely not necessary, though corner patching definitely is
grats<-lonlat2spilhaus(gen_graticules(lon=5,lat=5))


terra::plot(spilhaus2lonlat(grats)) #works without error, but quite messy
#collapse of graticules at poles and artifacts at IDL--pretty unavoidable I reckon
#can increase rip_seam_width to ensure IDL is simply erased
terra::plot(grats)
terra::plot(expand_borders(grats,amount=0.2))
abline(h=c(-11825474,11825474),
       v=c(-11825474,11825474))

coasts<-terra::vect(rnaturalearth::ne_download(type="coastline",category="physical"))
lakes<-terra::as.lines(terra::vect(rnaturalearth::ne_download(type="lakes",category="physical")))
coasts<-rbind(coasts,lakes)
coast2spil<-lonlat2spilhaus(coasts)
terra::plot(coast2spil)
terra::plot(grats,add=TRUE)

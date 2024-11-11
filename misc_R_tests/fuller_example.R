rm(list=ls())
library(rSpilhaus)
library(rnaturalearth)
library(terra)

####LOADING IN STUFF TO PLAY AROUND WITH####

lands<-vect(ne_download(type="land",category="physical"))
countries<-vect(ne_download(type="countries"))
lakes<-vect(ne_download(type="lakes",category="physical"))
coasts<-vect(ne_download(type="coastline",category="physical"))
depth.vec<-setNames(c(10:1*1000,200,0),
                    LETTERS[1:12])
depths<-setNames(vector("list",length(depth.vec)),rev(depth.vec))
for(i in names(depth.vec)){
  depths[[as.character(depth.vec[i])]]<-
    vect(ne_download(type=paste("bathymetry",
                                i,
                                depth.vec[i],
                                sep="_"),
                     scale=10,
                     category="physical"))
}
grats<-make_graticules(lon=15,lat=15)

####CONVERSION####

#let's make an "over-the-top" map of countries (including lakes), depths, and graticules...

#let's say we have some point data too
#don't sweat the next 4 steps
#just a quick way to generate uniformly distributed points on a sphere...
dat<-matrix(rnorm(3000),ncol=3)
dat<-dat/sqrt(rowSums(dat^2))
dat<-cbind("lon"=atan2(dat[,2],dat[,1]),
           "lat"=atan(dat[,3]/sqrt(rowSums(dat[,-3]^2))))
dat<-dat*180/pi
#now let's erase anything on land because we're talking MARINE bio here
dat<-dat[!is.related(vect(dat,"points"),lands,"within"),]
dat<-lonlat2spilhaus(dat,
                     lon="lon",lat="lat",
                     x="x",y="y")

#erase lakes from countries and reproject...
countries<-erase(countries,lakes)
spil.countries<-lonlat2spilhaus(countries)

#add lake coasts to coastlines and reproject...
coastlines<-rbind(coasts,as.lines(lakes))
spil.coastlines<-lonlat2spilhaus(coastlines)

#simplifying the depth data so it's easier to work with...
spil.depths<-lapply(depths,simplifyGeom,tolerance=0.1)
#reproject and combine all the depth data
spil.depths<-do.call(rbind,unname(lapply(spil.depths,lonlat2spilhaus)))
#make a color scale for depth values...
depth.cols<-hcl.colors(length(depth.vec)+1,"Blues")[-(length(depth.vec)+1)]
depth.cols<-depth.cols[match(values(spil.depths)$depth,depth.vec)]

#reproject graticules
spil.grats<-lonlat2spilhaus(grats)

#put it all together!
plot(spil.depths,
     col=depth.cols,
     border=NA)
plot(spil.countries,
     col="darkseagreen3",
     border="darkseagreen2",
     add=TRUE)
plot(spil.coastlines,
     col="darkseagreen3",
     add=TRUE)
points(dat[,c("x","y")],
       pch=16,cex=0.5)
plot(spil.grats,
     col="gray30",
     add=TRUE)

####TILING####

#note that tiling the depths data takes a hot sec...
exp.spil.depths<-expand_borders(spil.depths)
exp.spil.countries<-expand_borders(spil.countries)
exp.spil.coastlines<-expand_borders(spil.coastlines)
exp.spil.grats<-expand_borders(spil.grats)
exp.dat<-expand_borders(vect(dat[,c("x","y")],"points"))

#even the line and point data are properly tiled!
plot(exp.spil.depths,
     col=depth.cols,
     border=NA)
plot(exp.spil.countries,
     col="darkseagreen3",
     border="darkseagreen2",
     add=TRUE)
plot(exp.spil.coastlines,
     col="darkseagreen3",
     add=TRUE)
plot(exp.dat,
     add=TRUE,
     pch=16,cex=0.5)
plot(exp.spil.grats,
     col="gray30",
     add=TRUE)

####PRETTIFICATION####

#you really only want frame=TRUE for elements representing land
#frame=TRUE adds a nice solid box around the map edges
#this tends to look ugly if you join it with elements in the ocean
#note that frame has no effect on line/point data though!

#again, this takes a bit...
exp.spil.depths<-expand_borders(spil.depths,
                                prettify=TRUE,frame=FALSE)

#best to use simple land (rather than individual countries) here
#you could use countries, but it gets ugly due to many individual polygons!
#note that here frame is TRUE
exp.spil.lands<-expand_borders(lonlat2spilhaus(lands),
                               prettify=TRUE,frame=TRUE)

#no need to specify frame as FALSE for the following because line/point data...
exp.spil.coastlines<-expand_borders(spil.coastlines,
                                    prettify=TRUE)
exp.dat<-expand_borders(vect(dat[,c("x","y")],"points"),
                        prettify=TRUE)

#for graticules, I recommend either prettifying and erasing with land...
exp.spil.grats<-expand_borders(spil.grats,
                               prettify=TRUE)
exp.spil.grats<-erase(exp.spil.grats,exp.spil.lands)
#OR just tiling it again to cover everything...
# exp.spil.grats<-expand_borders(spil.grats)

#time to put it all together!
plot(exp.spil.depths,
     col=depth.cols,
     border=NA)
plot(exp.spil.lands,
     col="darkseagreen3",
     border=NA,
     add=TRUE)
#no benefit to plotting coastlines separately now...but whatevs
plot(exp.spil.coastlines,
     col="darkseagreen4",
     add=TRUE)
plot(exp.dat,
     add=TRUE,
     pch=16,cex=0.5)
plot(exp.spil.grats,
     col="gray30",
     add=TRUE)

####DEMONSTRATING CORNER ARTIFACTS####

#messing up polygons
bad.spil.countries<-lonlat2spilhaus(countries,NA_tol=5e5)
plot(bad.spil.countries)

#messing up lines
bad.spil.grats<-lonlat2spilhaus(make_graticules(lon=5,lat=5),NA_tol=5e5)
plot(bad.spil.grats)

#thankfully doesn't really affect point data
#well, other than erasing anything in the corners...
#but this the spilhaus projection...
#so you really shouldn't have any point data in the corners!

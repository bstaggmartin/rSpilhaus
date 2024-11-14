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
     col="darkseagreen4",
     add=TRUE)
points(dat[,c("x","y")],
       pch=16,cex=0.1)
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
     col="darkseagreen4",
     add=TRUE)
plot(exp.dat,
     add=TRUE,
     pch=16,cex=0.2)
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
                                prettify=TRUE)

#best to use simple land (rather than individual countries) here
#you could use countries, but it gets ugly due to many individual polygons!
#note that here frame is TRUE
exp.spil.lands<-expand_borders(lonlat2spilhaus(lands),
                               prettify=TRUE,frame=TRUE)

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
     pch=16,cex=0.2)
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

####MESSING AROUND WITH RASTERS####

#be sure to change destdir for your own setup here!!!
ne_download(scale=50,type="HYP_50M_SR",category="raster",destdir="misc_R_tests/rasters")
#annoyingly, had to manually package the downloaded files up for this to work...
land.rast<-
  ne_load(scale=50,type="HYP_50M_SR",category="raster",destdir="misc_R_tests/rasters")
ne_download(scale=50,type="OB_50M",category="raster",destdir="misc_R_tests/rasters")
water.rast<-
  ne_load(scale=50,type="OB_50M",category="raster",destdir="misc_R_tests/rasters")

#initial tests with land data...
spil.land.rast<-lonlat2spilhaus(land.rast)
plot(spil.land.rast)

#this seems to give a more "natural" result...
spil.land.rast<-lonlat2spilhaus(land.rast,
                                sample_method="bilinear",
                                patch_method="mean")
plot(spil.land.rast)

#how does it work with NA cells?
new.land.rast<-crop(land.rast,lands,mask=TRUE,touches=FALSE)
plot(new.land.rast,colNA="lightblue1")
spil.land.rast<-lonlat2spilhaus(new.land.rast,
                                sample_method="bilinear",
                                patch_method="mean")
plot(spil.land.rast,colNA="lightblue1") #nice

#testing expand_borders
plot(expand_borders(spil.land.rast),colNA="lightblue1")
plot(expand_borders(spil.land.rast,prettify=TRUE,frame=TRUE),colNA="lightblue1")
#trying to find prettier ways to fill in frame...
#this doesn't really work because most cells are NA in the beginning...
#oooh, maybe a weighting function will make it better?
# foo<-function(x){
#   if(sum(is.na(x))>length(x)/2){
#     NA
#   }else{
#     mean(x,na.rm=TRUE)
#   }
# }
#playing around with kernel types definitely helps...
#but still some ways to go
#perhaps a final bluring step? Seems to work well!
#A little fiddly (and modal now seems broken with more precise kernel specification)
#But overall not too bad
#Probably make it so that patch_kernel=NA defaults to no kernel spec

#Oooh, got a new one--progressive blurring?
plot(expand_borders(spil.land.rast,prettify=TRUE,frame=TRUE,
                    sample_method="bilinear",
                    patch_width=3,
                    patch_method="mean",
                    patch_width_inc=0,
                    patch_kernel=NA,
                    final_blur=0),
     colNA="lightblue1")
#funky-looking...
#will need to play around with this more in the future...
#I think the best strategy would be to keep blurring previous cells in each patching iteration...
#but terra doesn't yet support limiting such blurring operations to particular cells (unless those cells are NA)
#I can think of some potential workarounds with masking, but it doesn't seem worth it at the moment

#oh well, moving on for now
#erasing land from water?
new.water.rast<-crop(water.rast,
                     erase(terra::vect(cbind(c(-180,-180,180,180),
                                             c(-90,90,90,-90)),
                                       "polygons"),
                           lands),
                     mask=TRUE)
plot(new.water.rast,colNA="darkseagreen3")

#found it was actually best to use full water raster and plot land over it...
#this avoids awkward cell gaps in the final plot
spil.water.rast<-lonlat2spilhaus(water.rast,
                                 sample_method="bilinear",
                                 patch_method="mean")
plot(spil.water.rast)

#do a quick mockup with prettification...
plot(expand_borders(spil.water.rast,prettify=TRUE,
                    sample_method="bilinear"))
plot(expand_borders(spil.land.rast,prettify=TRUE,frame=TRUE,
                    sample_method="bilinear",
                    patch_width=10,patch_width_inc=2,
                    patch_method="mean"),
     add=TRUE,
     bgalpha=0)
plot(exp.spil.grats,
     add=TRUE,
     col="black")
#this honestly looks pretty bad, but you get the idea...
#(projecting country borders beyond prettification boundaries is just tricky...)
plot(expand_borders(spil.countries,prettify=TRUE,frame=TRUE),
     add=TRUE,
     border="gray90")

#how about one with tiling?
plot(expand_borders(spil.water.rast,
                    sample_method="bilinear"))
plot(expand_borders(spil.land.rast,
                    sample_method="bilinear"),
     add=TRUE,
     bgalpha=0)
plot(expand_borders(spil.countries),
     add=TRUE,
     border="gray90")
plot(expand_borders(spil.grats),
     add=TRUE,
     col="gray50")

#now for funsies --> one REALLY zoomed out
plot(expand_borders(spil.water.rast,
                    sample_method="bilinear",
                    amount=0.3))
plot(expand_borders(spil.land.rast,
                    sample_method="bilinear",
                    amount=0.3),
     add=TRUE,
     bgalpha=0)
plot(expand_borders(spil.countries,
                    amount=0.3),
     add=TRUE,
     border="gray90")
plot(expand_borders(lonlat2spilhaus(make_graticules(lon=10,lat=10)),
                    amount=0.3),
     add=TRUE,
     col="gray50")
#so delightfully odd-looking...

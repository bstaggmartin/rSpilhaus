library(rSpilhaus)

#download/aggregate land data (helpful for testing few large geoms)
og<-terra::vect(rnaturalearth::ne_download(scale=50,
                                           type="land",
                                           category="physical"))
og<-terra::aggregate(og)
#OR get country map (helpful for testing many small geoms)
# og<-terra::vect(rnaturalearth::countries110)
terra::plot(og)

#convert to spilhaus projection
og2spil<-lonlat2spilhaus(og)
terra::plot(og2spil)

#expand borders via tiling
tiled<-expand_borders(og2spil)
terra::plot(tiled)

#expand borders via "prettifying"
#I lifted the prettification border directly from Ricardo's code, but:
# - minor artifacts around Alaska towards top left
# - minor artifiacts in Central America towards bottom right
#Will have to further refine cropping polygon in the future
#10/30 update: working better now! Still get some tiny artifacts around...
#...top Alaska and both Central Americas at higher resolutions
#(probably just have to tweak existing coordinates a little bit)
prettified<-expand_borders(og2spil,prettify=TRUE,frame=TRUE)
terra::plot(prettified)

#can now convert back to lonlat coordinates quite easily
og2spil2og<-spilhaus2lonlat(og2spil)
terra::plot(og2spil2og)

#testing out the same with marine ecoregions
ecor<-terra::vect("C:/Users/bruce/Documents/spilhaus/Marine_Ecoregions_Of_the_World__MEOW_.shp")
ecor<-terra::project(ecor,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
terra::plot(ecor)

ecor2spil<-lonlat2spilhaus(ecor)
terra::plot(ecor2spil)

ecor.tiled<-expand_borders(ecor2spil)
terra::plot(ecor.tiled)

ecor.prettified<-expand_borders(ecor2spil,prettify=TRUE,frame=FALSE)
terra::plot(ecor.prettified)

ecor2spil2ecor<-spilhaus2lonlat(ecor2spil)
terra::plot(ecor2spil2ecor)

#doing a quick mock up of a proper map...ish
terra::plot(ecor.prettified,
            col=hcl.colors(20,"dark2"),
            lwd=1)
terra::plot(expand_borders(og2spil,
                           prettify=TRUE,
                           frame=TRUE,
                           amount=0.1),
            col="black",
            border=NA,
            add=TRUE)

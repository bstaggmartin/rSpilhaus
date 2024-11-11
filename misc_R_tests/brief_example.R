library(rSpilhaus)

#download/aggregate land data (helpful for testing few large geoms)
# og<-terra::vect(rnaturalearth::ne_download(scale=110,
                                           # type="land",
                                           # category="physical"))
# og<-terra::aggregate(og)
#OR get country map (helpful for testing many small geoms)
og<-terra::vect(rnaturalearth::countries110)
terra::plot(og)

lakes<-terra::vect(rnaturalearth::ne_download(scale=110,
                                              type="lakes",
                                              category="physical"))
og<-terra::erase(og,lakes)
# lakes2spil<-lonlat2spilhaus(lakes)
# og2spil<-terra::erase(og2spil,lakes2spil)
# terra::plot(og2spil)

#convert to spilhaus projection
og2spil<-lonlat2spilhaus(og)
terra::plot(og2spil)

#expand borders via tiling
tiled<-expand_borders(og2spil,amount=0.5)
terra::plot(tiled)

#expand borders via "prettifying"
#I lifted the prettification border directly from Ricardo's code, but:
# - minor artifacts around Alaska towards top left
# - minor artifiacts in Central America towards bottom right
#Will have to further refine cropping polygon in the future
#10/30 update: working better now! Still get some tiny artifacts around...
#...top Alaska and both Central Americas at higher resolutions
#(probably just have to tweak existing coordinates a little bit)
# terra::plot(expand_borders(og2spil,prettify=TRUE,frame=TRUE,cent=c(-3e6,3e6)))
# debug(expand_borders)
prettified<-expand_borders(og2spil,prettify=TRUE,frame=TRUE,amount=0.051)
#still getting minor artifact with land res of 50...
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

#so prettify totally breaks here for some reason...seems related quadrant detection...
#not a big deal here, so leave for later
ecor.prettified<-expand_borders(ecor2spil,prettify=TRUE,frame=FALSE)
terra::plot(ecor.prettified)

ecor2spil2ecor<-spilhaus2lonlat(ecor2spil)
terra::plot(ecor2spil2ecor)

#doing a quick mock up of a proper map...ish
terra::plot(ecor.prettified,
            col=hcl.colors(20,"dynamic"),
            lwd=1)
terra::plot(expand_borders(og2spil,
                           prettify=TRUE,
                           frame=TRUE,
                           amount=0.1),
            col="black",
            border=NA,
            add=TRUE)

.revalidate.geoms<-function(x){
  if(terra::geomtype(x)=="polygons"){
    terra::buffer(x,0)
  }else{
    x
  }
}

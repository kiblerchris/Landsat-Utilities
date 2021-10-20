preprocessLandsat <- function(inputPath, 
                              maskShapefile = NA, 
                              allBands = TRUE, 
                              maskClouds = TRUE, 
                              writeFile = TRUE, 
                              toMemory = FALSE,
                              fileFormat = "GTiff"){
  
  require(raster)
  
  setwd(inputPath)
  
  ### Iterate over every image subdirectory
  
  for (dir in list.dirs(path = ".", full.names = FALSE)[nchar(list.dirs(path = ".", full.names = FALSE)) > 0]){
    
    setwd(paste(inputPath, "/", dir, sep = ""))
    
    ### Create Raster Brick
    
    #List image files that are actually Landsat bands
    band_list <- list.files(path = ".", glob2rx("*band*.tif$"))
    
    #Drop aerosol band for Landsat 8
    if(substr(band_list[1], 1,4) == "LC08"){
      band_list <- band_list[-1]
      print("Recognized Landsat 8")
    } else if(substr(band_list[1], 1,4) %in% c("LT05", "LE07")){
      print("Recognized Landsat 5-7")
      } else {stop("Sensor Not Recognized")}
    
    print(paste("Processing ", substr(band_list[1], 18, 25), sep = ""))
    
    print(paste("Bands included in file:", band_list, sep = " "))
    
    #Merge bands into a brick
    image_brick <- stack(band_list)
    
    #Set NA value for output
    NAvalue(image_brick) <- -9999
    
    ### Crop and Mask Images Using Shapefile
    
    if(!is.na(maskShapefile)){
      
      polygon <- shapefile(maskShapefile)
      
      image_brick <- crop(image_brick, polygon)
      image_brick <- mask(image_brick, polygon)
      
      print("Crop and Mask using Shapefile")
      
    }
    
    ### Overwrite Invalid Pixel Values
    
    vals <- values(image_brick)
    vals[vals < 0] <- NA
    vals[vals > 10000] <- NA
    
    #Set all bands to NA if any band is NA for a pixel
    if(allBands == TRUE){

      vals[rowSums(is.na(vals)) > 0, ] <- NA
      
    }
    
    image_brick <- setValues(image_brick, vals)
    
    ### Mask Clouds
    
    if(maskClouds == TRUE){
      
      #Load pixel QA band
      pixel_qa_band_name <- list.files(path = ".", glob2rx("*pixel_qa*.tif$"))
      pixel_qa_band <- raster(pixel_qa_band_name)
      
      if(!is.na(maskShapefile)){
        
        #Crop to same extent as image brick
        pixel_qa_band <- crop(pixel_qa_band, polygon)
        pixel_qa_band <- mask(pixel_qa_band, polygon)
        
      }
      
      if(substr(band_list[1], 1,4) %in% c("LT05", "LE07")){
        
        #Mask medium and high confidence clouds
        image_brick[pixel_qa_band %in% c(130, 132, 136, 144, 160, 176, 224)] <- NA
        
        print("QA Mask for Landsat 5-7")
        
      } else if(substr(band_list[1], 1,4) == "LC08"){
        
        #Medium and high confidence clouds and high confidence cirrus
        image_brick[pixel_qa_band %in% c(386, 388, 392, 400, 416, 432, 898, 900, 904, 928, 944, 480, 992, 
                                         834, 836, 840, 848, 864, 880, 898, 900, 904, 912, 928, 944, 992)] <- NA
        
        print("QA Mask for Landsat 8")
        
      } else{warning("No QA Mask Found")}
      
    }
    
    ### Export Output
    
    setwd(inputPath)
    
    if(writeFile == TRUE){
      
      writeRaster(image_brick, 
                  substr(band_list[1], 1, nchar(band_list[1]) - 10), 
                  format = fileFormat, 
                  NAflag = -9999)
      
    }
    
    if(toMemory == TRUE){
      
      assign(substr(band_list[1], 1, nchar(band_list[1]) - 10), 
             image_brick,
             envir = .GlobalEnv)
      
    }
    
    print("Image completed")
    
  }
}

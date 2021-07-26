prepCollection2 <- function(inputPath,
                            dropCoastal = TRUE,
                            maskShapefile = NA,
                            maskValidRange = TRUE,
                            allBands = TRUE, 
                            maskClouds = TRUE, 
                            writeFile = TRUE, 
                            toMemory = FALSE,
                            fileFormat = "GTiff"){
  
  require(raster)
  
  setwd(inputPath)
  
  ### Iterate over every image subdirectory
  
  for (dir in list.dirs(path = ".", full.names = FALSE, recursive = FALSE)){
    
    setwd(paste(inputPath, "/", dir, sep = ""))
    
    ### Create Raster Brick
  
    #Find bands in directory
    band_list <- list.files(path = ".", glob2rx("*SR_B*.TIF$"))
    
    print(paste0("Processing: ", paste(unlist(strsplit(band_list[1], "_"))[c(1, 4)], collapse = "_")))
    
    if(dropCoastal == TRUE){
      
      if(substr(band_list[1], 1,4) == "LC08"){
        band_list <- band_list[-1]
        print("Recognized Landsat 8 and dropped coastal band")
      } else if(substr(band_list[1], 1,4) %in% c("LT05", "LE07")){
        print("Recognized Landsat 5-7")
      } else {stop("Sensor not recognized")}
      
    }
  
    #Print band list
    print("Bands included in file:")
    print(band_list)
    
    #Merge bands into a brick
    image_brick <- brick(stack(band_list))
    
    #Set NA value for output
    NAvalue(image_brick) <- -9999
    
    print("Created brick")
    
    #Apply correction coefficients for Collection 2
    image_brick <- calc(image_brick, fun = function(x){x * 0.0000275 + (-0.2)})
    
    print("Calculated reflectance values")
    
    ### Overwrite Invalid Pixel Values
    
    if(maskValidRange == TRUE){
      
      vals <- values(image_brick)
      vals[vals < 0] <- NA
      vals[vals > 1] <- NA
      
      #Set all bands to NA if any band is NA for a pixel
      if(allBands == TRUE){
        
        vals[rowSums(is.na(vals)) > 0, ] <- NA
        
      }
      
      image_brick <- setValues(image_brick, vals)
      
      print("Masked erroneous values")
      
      if(!is.na(maskShapefile)){
        
        polygon <- shapefile(maskShapefile)
        
        image_brick <- crop(image_brick, polygon)
        image_brick <- mask(image_brick, polygon)
        
        print("Cropped and masked using shapefile")
        
      }
      
    }
    
    if(maskClouds == TRUE){
      
      #Load pixel QA band
      pixel_qa_band_name <- list.files(path = ".", glob2rx("*QA_PIXEL*.TIF$"))
      pixel_qa_band <- raster(pixel_qa_band_name)
      
      print("Loaded QA band")
      
      if(!is.na(maskShapefile)){
        
        #Crop to same extent as image brick
        pixel_qa_band <- crop(pixel_qa_band, polygon)
        pixel_qa_band <- mask(pixel_qa_band, polygon)
        
      }
      
      #Mask fill values
      image_brick[bitwAnd(bitwShiftL(1, 0), values(pixel_qa_band)) != 0] <- NA
      
      #Mask high confidence cirrus (for Landsat 8 only)
      image_brick[bitwAnd(bitwShiftL(1, 2), values(pixel_qa_band)) != 0] <- NA
      
      #Mask high confidence clouds
      image_brick[bitwAnd(bitwShiftL(1, 3), values(pixel_qa_band)) != 0] <- NA
      
      #Mask high confidence cloud shadows
      image_brick[bitwAnd(bitwShiftL(1, 4), values(pixel_qa_band)) != 0] <- NA
      
      print("Masked fill values, high confidence clouds, and cloud shadows")
      
    }
  
    rm(pixel_qa_band)
    
    ### Export Output
    
    setwd(inputPath)
    
    if(writeFile == TRUE){
      
      writeRaster(image_brick, 
                  filename = paste(unlist(strsplit(band_list[1], "_"))[c(1, 4)], collapse = "_"), 
                  format = fileFormat, 
                  NAflag = -9999)
      
    }
    
    if(toMemory == TRUE){
      
      assign(paste(unlist(strsplit(band_list[1], "_"))[c(1, 4)], collapse = "_"), 
             image_brick,
             envir = .GlobalEnv)
      
    }
    
    print("Image completed")
  
  }
}
  
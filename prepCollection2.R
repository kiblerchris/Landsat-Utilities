######################################################
######################################################
### Landsat Collection 2, Level 2 Preprocessing
### Code by Christopher Kibler
### Last Updated September 30, 2021
### Please report any bugs to kibler@ucsb.edu
######################################################
######################################################

#How to use:
# 1. Download Landsat 4-8 Collection 2, Level 2 products from USGS Earth Explorer.
# 2. Move the zip file(s) into a new directory. 
# 3. Unzip the files so each image has its own subdirectory. DO NOT CHANGE THE DEFAULT FILE NAMES.
# 4. Each subdirectory should contain all of the individual band files for the image.
# 5. Set the R working directory to the main directory that contains the image subdirectories.
# 6. Run the code below. The code is very verbose and should describe all of the processes.

#Arguments:
# inputPath: The file path for the main directory that contains the image subdirectories
# dropCoastal: Drop the Landsat 8 coastal/aerosol band (band 1)
# maskShapefile: File path for an optional shapefile (.shp) that will be used to mask the images
# maskValidRange: Set all surface reflectance values outside of 0-1 range to NA
# allBands: If a value in any band is outside of the 0-1 range, set all bands for that pixel to NA
# maskClouds: Mask high confidence clouds (all sensors) and cirrus (Landsat 8) using pixel QA band
# maskSaturated: Mask saturated pixels using the radiometric QA band
# writeFile: Write output image to current working directory
# toMemory: Assign output image to a new object in memory
# fileFormat: Output format for image file. See ?writeRaster for options.

prepCollection2 <- function(inputPath,
                            dropCoastal = TRUE,
                            maskShapefile = NA,
                            maskValidRange = TRUE,
                            allBands = TRUE, 
                            maskClouds = TRUE, 
                            maskSaturated = TRUE,
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
    
    image_id <- band_list[1]
    image_attributes <- strsplit(image_id, ".", fixed = TRUE)[[1]][1]
    image_attributes <- as.list(strsplit(image_attributes, "_")[[1]])
    names(image_attributes) <- c("sensor", 
                                 "correction_level", 
                                 "pathrow", 
                                 "acquisition_date", 
                                 "production_date",
                                 "collection", 
                                 "tier",
                                 "product", 
                                 "first_band")
    
    print(image_attributes)
    
    print(paste("Processing:", image_attributes$sensor, image_attributes$acquisition_date, sep = " "))
    
    if(dropCoastal == TRUE){
      
      if(image_attributes$sensor == "LC08"){
        band_list <- band_list[-1]
        print("Recognized Landsat 8 and dropped coastal band")
      } else if(image_attributes$sensor %in% c("LT05", "LE07")){
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
    
    ### Mask to Shapefile Boundaries
    
    if(!is.na(maskShapefile)){
      
      polygon <- shapefile(maskShapefile)
      
      image_brick <- crop(image_brick, polygon)
      image_brick <- mask(image_brick, polygon)
      
      print("Cropped and masked using shapefile")
      
    }
    
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
      
      print("Masked valid range")
      
    }
    
    if(maskClouds == TRUE){
      
      #Load pixel QA band
      pixel_qa_band_name <- list.files(path = ".", glob2rx("*QA_PIXEL*.TIF$"))
      pixel_qa_band <- raster(pixel_qa_band_name)
      
      print("Loaded pixel QA band")
      
      if(!is.na(maskShapefile)){
        
        #Crop to same extent as image brick
        pixel_qa_band <- crop(pixel_qa_band, polygon)
        pixel_qa_band <- mask(pixel_qa_band, polygon)
        
      }
      
      #Mask fill values
      image_brick[bitwAnd(bitwShiftL(1, 0), values(pixel_qa_band)) != 0] <- NA
      
      if(image_attributes$sensor == "LC08"){
        
        #Mask high confidence cirrus (for Landsat 8 only)
        image_brick[bitwAnd(bitwShiftL(1, 2), values(pixel_qa_band)) != 0] <- NA
        
        print("Masked cirrus")
      }
      
      #Mask high confidence clouds
      image_brick[bitwAnd(bitwShiftL(1, 3), values(pixel_qa_band)) != 0] <- NA
      
      #Mask high confidence cloud shadows
      image_brick[bitwAnd(bitwShiftL(1, 4), values(pixel_qa_band)) != 0] <- NA
      
      print("Masked fill values, high confidence clouds, and cloud shadows")
      
      rm(pixel_qa_band_name, pixel_qa_band)
      
    }
  
    if(maskSaturated == TRUE){
      
      radsat_qa_band_name <- list.files(path = ".", glob2rx("*QA_RADSAT*.TIF$"))
      radsat_qa_band <- raster(radsat_qa_band_name)
      
      print("Loaded saturation QA band")
      
      if(!is.na(maskShapefile)){
        
        #Crop to same extent as image brick
        radsat_qa_band <- crop(radsat_qa_band, polygon)
        radsat_qa_band <- mask(radsat_qa_band, polygon)
        
      }
      
      image_brick[radsat_qa_band != 0] <- NA
      
      print("Masked saturated values")
      
      rm(radsat_qa_band_name, radsat_qa_band)
      
    }
    
    ### Export Output
    
    setwd(inputPath)
    
    if(writeFile == TRUE){
      
      writeRaster(image_brick, 
                  filename = paste(image_attributes$sensor, image_attributes$acquisition_date, sep = "_"), 
                  format = fileFormat, 
                  NAflag = -9999)
      
    }
    
    if(toMemory == TRUE){
      
      assign(paste(image_attributes$sensor, image_attributes$acquisition_date, sep = "_"), 
             image_brick,
             envir = .GlobalEnv)
      
    }
    
    print("Image completed")
  
  }
}

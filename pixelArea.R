pixelArea <- function(number_of_pixels, output_unit = "km2"){
  
  require(measurements)
  
  m2 <- 30 * 30 * number_of_pixels
  
  output <- conv_unit(m2, from = "m2", to = output_unit)
  
  print(paste(output, output_unit, sep = " "))
  
  return(output)
  
}

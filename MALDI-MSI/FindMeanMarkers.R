# FindMeanMarkers function 
FindMeanMarkers <- function(object, features = NULL, ident.1 = NULL, ident.2 = NULL) {

  data <- GetAssayData(object = object, layer = "data") 
  if (is.null(features)) {
    features <- rownames(object) 
  }

  mean.fxn <- function(x) {
    return(rowMeans(x, na.rm = TRUE))
  }
  
  cells.1 <- WhichCells(object = object, idents = ident.1) 
  cells.2 <- WhichCells(object = object, idents = ident.2) 
  
  results <- data.frame(
    Feature = features,
    Mean_Ident1 = NA,
    Mean_Ident2 = NA,
    FC = NA
  )

  for (i in seq_along(features)) {
    feature <- features[i]
    
    data.1 <- mean.fxn(data[feature, cells.1, drop = FALSE] )
    data.2 <- mean.fxn(data[feature, cells.2, drop = FALSE] )

    FC <- data.1 / data.2

    results[i, "Mean_Ident1"] <- data.1
    results[i, "Mean_Ident2"] <- data.2
    results[i, "FC"] <- FC
  }
  return(results)
}



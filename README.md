The connect_trees function is used to spatially link repeated forest inventory data derived from point clouds.   

The inventory datasets must contain the fields X, Y, and BHD (breast height diameter) in centimeters (.csv, UTF8).   

The output of the function is an adjusted version of the new inventory, in which trees have been spatially matched (connected) to those from the previous inventory. 

# Usage:

result <- connect_trees(
  inventur_alt_path = "path/to/Inventur_alt.csv",
  inventur_neu_path = "path/to/Inventur_neu.csv"
)

transformed <- result$Inventur_neu_transformed
head(transformed)

## Plotting

plot(transformed$X, transformed$Y, pch = 19, col = "black",cex = 0.7)
points(result$Inventur_alt_filtered$X, result$Inventur_alt_filtered$Y, pch = 20, col = "red", cex = 0.5)
 

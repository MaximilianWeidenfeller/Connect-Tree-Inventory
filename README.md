Die Funktion connect_trees dient dazu, serielle Inventurdaten aus Punktwolken von Wäldern räumlich zu verbinden. 
Die Inventurdaten müssen  X, Y, und BHD in cm enthalten. Das Ergebnis der Funktion ist eine angepasste Version der neuen Inventur.

Anwendung:

result <- connect_trees(
  inventur_alt_path = "path/to/Inventur_alt.csv",
  inventur_neu_path = "path/to/Inventur_neu.csv"
)

transformed <- result$Inventur_neu_transformed
head(transformed)

# Optional: plotten
plot(transformed$X, transformed$Y, pch = 19, col = "black",cex = 0.7)
points(result$Inventur_alt_filtered$X, result$Inventur_alt_filtered$Y, pch = 20, col = "red", cex = 0.5)
 

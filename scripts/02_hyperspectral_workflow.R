############################################################
# HYPERSPECTRAL WORKFLOW
# Goal: Dimensionality reduction + Trait prediction
############################################################

# ---- STEP 1: Load packages ----
library(terra)
library(sf)
library(ranger)
library(pls)

# ---- STEP 2: Load hyperspectral cube ----
hs <- rast("hyperspectral_cube.tif")

# ---- STEP 3: Load field sample points with trait values ----
samples <- st_read("sample_points.shp")
samples_v <- vect(samples)

# ---- STEP 4: Extract spectral signatures ----
spectra <- extract(hs, samples_v)

# ---- STEP 5: Prepare matrix for modelling ----
spectra_matrix <- as.matrix(spectra[,-1])
trait <- samples$trait_value  # Replace with real trait

# ---- STEP 6: PCA (optional dimensionality reduction) ----
pca <- prcomp(spectra_matrix, scale.=TRUE)
plot(pca$x[,1], pca$x[,2], main="PCA of Spectra")

# ---- STEP 7: Random Forest model ----
df <- data.frame(trait=trait, spectra_matrix)
rf_model <- ranger(trait ~ ., data=df, importance="impurity")

print(rf_model)

# ---- STEP 8: Predict trait map ----
pred_map <- predict(hs, rf_model)
plot(pred_map, main="Predicted Trait Map")

############################################################

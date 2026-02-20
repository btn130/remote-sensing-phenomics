############################################################
# MULTISPECTRAL PHENOMICS WORKFLOW (Sentinel-2 example)
# Goal: Extract vegetation indices + plot-level traits
############################################################

# ---- STEP 1: Load required packages ----
library(terra)
library(sf)

# ---- STEP 2: Load multispectral bands ----
# Replace with your actual band paths
nir  <- rast("B08.tif")      # NIR band
red  <- rast("B04.tif")      # Red band
rededge <- rast("B05.tif")   # Red edge band

# ---- STEP 3: Compute vegetation indices ----
ndvi <- (nir - red) / (nir + red)
ndre <- (nir - rededge) / (nir + rededge)

# ---- STEP 4: Visualize indices ----
plot(ndvi, main="NDVI")
plot(ndre, main="NDRE")

# ---- STEP 5: Load plot shapefile ----
plots <- st_read("plots.shp")
plots_v <- vect(plots)

# ---- STEP 6: Extract plot-level traits ----
traits <- extract(c(ndvi, ndre), plots_v, fun=mean, na.rm=TRUE)

# ---- STEP 7: Save extracted traits ----
write.csv(traits, "multispectral_traits.csv", row.names=FALSE)

############################################################

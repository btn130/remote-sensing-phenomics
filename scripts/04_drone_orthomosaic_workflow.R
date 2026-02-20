############################################################
# DRONE ORTHOMOSAIC WORKFLOW
# Goal: Plot-level trait extraction
############################################################

# ---- STEP 1: Load packages ----
library(terra)
library(sf)

# ---- STEP 2: Load drone orthomosaic ----
ortho <- rast("drone_multispectral.tif")

# ---- STEP 3: Compute NDVI (example: NIR=4, Red=3) ----
ndvi <- (ortho[[4]] - ortho[[3]]) / (ortho[[4]] + ortho[[3]])
plot(ndvi, main="Drone NDVI")

# ---- STEP 4: Create simple plot grid ----
ext <- ext(ndvi)
plots_grid <- as.polygons(ext, n=100)

# ---- STEP 5: Extract plot-level NDVI ----
traits <- extract(ndvi, plots_grid, fun=mean, na.rm=TRUE)

# ---- STEP 6: Save traits ----
write.csv(traits, "drone_traits.csv", row.names=FALSE)

############################################################

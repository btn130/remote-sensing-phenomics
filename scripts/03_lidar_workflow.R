############################################################
# LIDAR WORKFLOW
# Goal: Canopy height + structural traits
############################################################

# ---- STEP 1: Load packages ----
library(lidR)
library(terra)

# ---- STEP 2: Read LiDAR point cloud ----
las <- readLAS("lidar.laz")

# ---- STEP 3: Generate Digital Terrain Model (DTM) ----
dtm <- grid_terrain(las, res=1, algorithm=knnidw())

# ---- STEP 4: Normalize heights ----
las_norm <- normalize_height(las, dtm)

# ---- STEP 5: Create Canopy Height Model (CHM) ----
chm <- grid_canopy(las_norm, res=1, p2r())
plot(chm, main="Canopy Height Model")

# ---- STEP 6: Extract canopy structural metrics ----
metrics <- grid_metrics(
  las_norm,
  ~list(
    mean_h = mean(Z),
    p95 = quantile(Z, 0.95),
    cover = sum(Z > 2)/length(Z)
  ),
  res=5
)

plot(metrics)

############################################################

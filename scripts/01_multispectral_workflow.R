############################################################
# 01_multispectral_sentinel2.R
# MULTISPECTRAL PHENOMICS WORKFLOW (Sentinel-2 L2A)
# - Download via STAC (Earth-Search)
# - Cloud mask using SCL
# - Create monthly median NDVI + NDRE composite
# - Extract plot-level traits (mean/SD)
############################################################

# ---- STEP 0: Packages ----
pkgs <- c("rstac","sf","terra","dplyr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, require, character.only = TRUE))

# ---- STEP 1: User settings ----
stac_url <- "https://earth-search.aws.element84.com/v1"
collection <- "sentinel-2-l2a"

# AOI polygon (recommended). Use your Bhutan boundary or field boundary.
aoi_path <- "data/shapefiles/DzongLatLong.shp"          # <- put a small AOI here
#sf::st_layers(aoi_path) # run this to see what layer is present
aoi_layer <- "DzongLatLong"                   # <- layer name inside gpkg

# Optional plot polygons (for plot-level phenomics traits)
plots_path <- "data/plots.gpkg"      # <- optional
plots_layer <- "plots"

# Target month
YEAR <- 2021
MONTH <- 4

# Keep only one dzongkhag (CHANGE NAME HERE)
dz <- "Paro"   # e.g., "Bumthang", "Monggar", etc.

# Output folders
dir.create("data/raw/s2", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/rasters", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

out_ndvi <- sprintf("outputs/rasters/NDVI_%s_%04d_%02d_median.tif", dz, YEAR, MONTH)
out_ndre <- sprintf("outputs/rasters/NDRE_%s_%04d_%02d_median.tif", dz, YEAR, MONTH)
out_csv  <- sprintf("outputs/tables/plot_traits_%s_%04d_%02d.csv", dz, YEAR, MONTH)

# ---- STEP 2: Read AOI (must be polygon) and prep bbox for STAC ----
# Read AOI and subset ONE dzongkhag ----
aoi_sf <- sf::st_read(aoi_path, quiet = TRUE) |> 
  sf::st_make_valid() |> 
  sf::st_transform(4326)

# Keep only one district for faster processing (CHANGE NAME HERE)
# Your column appears to be DzgName
aoi_sf <- aoi_sf[aoi_sf$DzgName == dz, ]

if (nrow(aoi_sf) == 0) stop("Dzongkhag not found in DzgName: ", dz)

# Dissolve to single polygon (helps crop/mask)
aoi_sf <- sf::st_union(aoi_sf) |> sf::st_as_sf()

bb <- sf::st_bbox(aoi_sf)
bbox <- c(bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"])

# ---- STEP 3: Month datetime range (split into halves to reduce 502 risk) ----
date_start <- as.Date(sprintf("%04d-%02d-01", YEAR, MONTH))
date_end   <- as.Date(seq(date_start, by="month", length.out=2)[2] - 1)

ranges <- c(
  sprintf("%sT00:00:00Z/%sT23:59:59Z", format(date_start, "%Y-%m-%d"),
          format(min(date_start + 14, date_end), "%Y-%m-%d")),
  sprintf("%sT00:00:00Z/%sT23:59:59Z", format(date_start + 15, "%Y-%m-%d"),
          format(date_end, "%Y-%m-%d"))
)

# ---- STEP 4: Robust STAC request with retries ----
post_request_retry <- function(req, tries = 8, base_sleep = 4) {
  for (i in seq_len(tries)) {
    out <- try(rstac::post_request(req), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    message("STAC request failed (try ", i, "/", tries, "). Retrying...")
    Sys.sleep(base_sleep * (2^(i - 1)) + runif(1, 0, 2))
  }
  stop("STAC request failed after retries.")
}

get_items_bbox <- function(datetime_range, limit = 200) {
  req <- rstac::stac(stac_url) |>
    rstac::stac_search(
      collections = collection,
      bbox        = as.numeric(bbox),
      datetime    = datetime_range,
      limit       = limit
    )
  Sys.sleep(runif(1, 1, 4))
  post_request_retry(req)
}

# ---- STEP 5: Download item metadata (not the imagery) ----
items1 <- get_items_bbox(ranges[1], limit = 80) # this limit can increase to 200 if whole of Bhutan
items2 <- get_items_bbox(ranges[2], limit = 80)
features <- c(items1$features, items2$features)
stopifnot(length(features) > 0)

# ---- STEP 6: (Optional) filter by cloud cover ----
max_cloud <- 60
features <- Filter(function(f) {
  cc <- f$properties$`eo:cloud_cover`
  if (is.null(cc)) return(TRUE)
  return(cc <= max_cloud)
}, features)

stopifnot(length(features) > 0)

# ---- STEP 7: Cloud-mask + NDVI/NDRE per scene (stream COGs directly; no huge downloads) ----
get_href <- function(f, key) {
  h <- f$assets[[key]]$href
  if (!is.null(h)) return(h)
  stop("Missing asset '", key, "'. Available: ", paste(names(f$assets), collapse=", "))
}

# SCL classes to mask (cloud/shadow/cirrus/snow/no-data)
bad_scl <- c(0,1,3,8,9,10,11)

aoi_v <- terra::vect(aoi_sf)

ndvi_list <- list()
ndre_list <- list()

for (i in seq_along(features)) {
  f <- features[[i]]
  
  # Earth-Search S2 assets are named: red, nir, rededge1, scl (as you saw)
  red_href     <- get_href(f, "red")
  nir_href     <- get_href(f, "nir")
  rededge_href <- get_href(f, "rededge1")
  scl_href     <- get_href(f, "scl")
  
  red <- terra::rast(red_href)
  nir <- terra::rast(nir_href)
  re1 <- terra::rast(rededge_href)
  scl <- terra::rast(scl_href)
  
  # Project AOI to raster CRS for correct crop/mask
  aoi_r <- terra::project(aoi_v, terra::crs(nir))
  aoi_ext <- terra::ext(aoi_r)
  
  # Skip if this tile doesn’t overlap AOI (prevents crop error)
  if (!terra::relate(terra::ext(nir), aoi_ext, "intersects")) next
  
  # Crop by extent (fast), then mask by polygon (exact)
  red <- terra::mask(terra::crop(red, aoi_ext), aoi_r)
  nir <- terra::mask(terra::crop(nir, aoi_ext), aoi_r)
  re1 <- terra::mask(terra::crop(re1, aoi_ext), aoi_r)
  scl <- terra::mask(terra::crop(scl, aoi_ext), aoi_r)
  
  # Align grids
  if (!isTRUE(all.equal(terra::res(red), terra::res(nir)))) red <- terra::resample(red, nir, "bilinear")
  if (!isTRUE(all.equal(terra::res(re1), terra::res(nir)))) re1 <- terra::resample(re1, nir, "bilinear")
  if (!isTRUE(all.equal(terra::res(scl), terra::res(nir)))) scl <- terra::resample(scl, nir, "near")
  
  # NDVI / NDRE
  ndvi <- (nir - red) / (nir + red)
  ndre <- (nir - re1) / (nir + re1)
  
  # Cloud mask using SCL
  ndvi <- terra::mask(ndvi, scl, maskvalues = bad_scl, updatevalue = NA)
  ndre <- terra::mask(ndre, scl, maskvalues = bad_scl, updatevalue = NA)
  
  ndvi_list[[length(ndvi_list)+1]] <- ndvi
  ndre_list[[length(ndre_list)+1]] <- ndre
  
  message("Processed scene ", i, "/", length(features))
}

if (length(ndvi_list) == 0) {
  stop("No valid scenes overlapped AOI for this month.")
}

# ---- STEP 8: Monthly composite (median is robust) ----
ndvi_stack <- terra::rast(ndvi_list)
ndre_stack <- terra::rast(ndre_list)

ndvi_month <- terra::app(ndvi_stack, median, na.rm = TRUE)
ndre_month <- terra::app(ndre_stack, median, na.rm = TRUE)

# ---- STEP 9: Save rasters ----
terra::writeRaster(ndvi_month, out_ndvi, overwrite = TRUE)
terra::writeRaster(ndre_month, out_ndre, overwrite = TRUE)

# ---- STEP 10: Plot-level traits (optional) ----
if (file.exists(plots_path)) {
  plots_sf <- sf::st_read(plots_path, layer = plots_layer, quiet = TRUE) |> sf::st_make_valid()
  plots_v  <- terra::vect(plots_sf) |> terra::project(terra::crs(ndvi_month))
  
  trait_mean <- terra::extract(c(ndvi_month, ndre_month), plots_v, fun = mean, na.rm = TRUE)
  trait_sd   <- terra::extract(c(ndvi_month, ndre_month), plots_v, fun = sd,   na.rm = TRUE)
  
  out <- cbind(
    trait_mean,
    ndvi_sd = trait_sd[, names(ndvi_month)],
    ndre_sd = trait_sd[, names(ndre_month)]
  )
  write.csv(out, out_csv, row.names = FALSE)
}

# ---- STEP 11: Coverage metric (how much is cloud-free) ----
cov_pct <- terra::global(!is.na(ndvi_month), "mean", na.rm = TRUE)[1,1] * 100
message(sprintf("Monthly NDVI clear-sky coverage: %.1f%%", cov_pct))
writeLines(capture.output(sessionInfo()), 
           sprintf("outputs/tables/sessionInfo_%04d_%02d.txt", YEAR, MONTH))


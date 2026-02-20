\# 🌿 Remote Sensing Phenomics  

\*\*Author:\*\* Sangay Dorji  

PhD (Climate \& Agricultural Modelling), University of Queensland  



---



\## 📌 Project Overview



This repository demonstrates applied remote sensing workflows for plant phenotyping (phenomics), including:



\- Multispectral vegetation indices (Sentinel-2)

\- Hyperspectral trait modelling

\- LiDAR canopy structure analysis

\- Drone-based plot-level trait extraction



The objective is to showcase reproducible workflows for extracting plant structural and biochemical traits from remote sensing imagery.



---



\## 🛰️ 1. Multispectral Workflow (Sentinel-2)



\### Objectives

\- Compute NDVI and NDRE

\- Generate phenology time series

\- Extract plot-level vegetation traits



\### Methods

1\. Surface reflectance preprocessing

2\. Vegetation index computation

3\. Plot extraction using shapefiles

4\. Statistical summarisation



\### Outputs

\- `NDVI\_map.tif`

\- `plot\_traits.csv`

\- Seasonal phenology curve



---



\## 🌈 2. Hyperspectral Trait Prediction



\### Objectives

\- Extract spectral signatures (100+ bands)

\- Perform dimensionality reduction (PCA)

\- Predict plant traits using machine learning (RF / PLSR)



\### Methods

\- Spectral extraction from hyperspectral cube

\- PCA analysis

\- Random Forest regression

\- Full-scene prediction mapping



\### Outputs

\- Trait prediction map

\- Model performance metrics

\- Variable importance plot



---



\## 🌳 3. LiDAR Structural Metrics



\### Objectives

\- Generate canopy height model (CHM)

\- Extract structural traits (mean height, P95, canopy cover)



\### Methods

\- Point cloud normalization

\- Terrain modelling

\- Grid-based canopy metrics



\### Outputs

\- CHM raster

\- Structural metrics map

\- Height distribution analysis



---



\## 🚁 4. Drone Phenomics Workflow



\### Objectives

\- Extract plot-level vegetation indices

\- Compute canopy cover

\- Generate trait tables for statistical analysis



\### Methods

\- Orthomosaic loading

\- NDVI computation

\- Grid-based plot segmentation

\- Trait extraction



---



\## 📊 Technologies Used



\- R

\- terra

\- sf

\- lidR

\- ranger

\- pls



---



\## 🔁 Reproducibility



To run the workflows:



```bash

Rscript scripts/01\_multispectral\_workflow.R



\## Data Access



Sentinel-2 data can be downloaded from:

https://earth-search.aws.element84.com



Hyperspectral data (NEON):

https://data.neonscience.org



LiDAR data:

https://lpdaac.usgs.gov/products/gedi02\_a\_v002/



This repository does not include raw imagery due to file size limitations.

All data sources are publicly available and links are provided above.






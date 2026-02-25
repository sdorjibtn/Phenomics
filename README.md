\# Phenomics Pipeline – UAV LiDAR + Sentinel-2 Integration



Author: Sangay Dorji



This repository demonstrates a phenomics workflow integrating:



\- Sentinel-2 multispectral data (NDVI, NDRE)

\- UAV LiDAR canopy metrics (zq95, CHM)

\- Random Forest modelling (ranger)

\- Spatial prediction and validation



\## Study Area

Gatton, Queensland, Australia



\## Workflow

1\. Sentinel-2 monthly composites

2\. LiDAR canopy metric extraction

3\. CRS alignment (WGS84 UTM56S → GDA94 MGA56)

4\. Pixel sampling

5\. Random Forest regression

6\. Variable importance analysis

7\. Full-area trait prediction



\## Model Performance (Example)

\- RMSE ≈ 1.6 m

\- R² ≈ 0.60–0.64 (5-fold CV)



\## Skills Demonstrated

\- Remote sensing preprocessing

\- CRS handling (EPSG 28356 / 32756)

\- Raster alignment

\- Machine learning (Random Forest)

\- Cross-validation

\- Reproducible workflow


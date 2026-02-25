#!/usr/bin/env Rscript
# ============================================================
# Hyperspectral Workflow (EnMAP L2A) | Laptop + Bunya
# ------------------------------------------------------------
# Author: Sangay Dorji
#
# What this script does:
#  1) Loads an EnMAP L2A hyperspectral cube (surface reflectance)
#  2) Clips to an AOI (bbox in WGS84)
#  3) Applies optional quality/cloud mask (if a mask raster exists)
#  4) Computes key hyperspectral vegetation indices (NDVI, NDRE, PRI)
#  5) Runs PCA (first 3 components) to summarize spectra
#  6) Exports GeoTIFFs + quicklook PNGs + AOI mean spectra CSV
#
# Notes:
#  • EnMAP is hyperspectral (hundreds of bands). Be mindful of memory.
#  • Use env vars for paths so this runs unchanged on laptop and Bunya.
#  • This is analysis-ready once you have EnMAP files locally.
#
# Recommended env vars:
#  - ENMAP_ROOT    : project root folder
#  - ENMAP_CUBE    : path to the EnMAP L2A reflectance cube GeoTIFF
#  - ENMAP_MASK    : optional path to a cloud/quality mask (0/1 or classes)
#  - ENMAP_OUTDIR  : output folder (optional)
#  - ENMAP_TMPDIR  : temp folder (recommended on HPC scratch)
#  - TERRA_THREADS : threads
#  - REMOVE_TMP_DIR: TRUE/FALSE
#
# Example (Windows PowerShell):
#  $env:ENMAP_ROOT="C:\Users\dorji\Documents\Phenomics\hyperspectral\enmap"
#  $env:ENMAP_CUBE="C:\...\ENMAP01_L2A_SPECTRAL_IMAGE.TIF"
#  $env:ENMAP_MASK="C:\...\CLOUD_MASK.TIF"   # optional
#  $env:TERRA_THREADS="4"
#  Rscript 01_enmap_hyperspectral_workflow.R
#
# Example (Bunya bash):
#  export ENMAP_ROOT=/scratch/user/s4597131/phenomics/hyperspectral/enmap
#  export ENMAP_CUBE=/scratch/user/s4597131/phenomics/hyperspectral/enmap/data/ENMAP_L2A_SPECTRAL_IMAGE.tif
#  export ENMAP_MASK=/scratch/user/s4597131/phenomics/hyperspectral/enmap/data/CLOUD_MASK.tif
#  export ENMAP_TMPDIR=/scratch/user/s4597131/phenomics/tmp_enmap
#  export TERRA_THREADS=8
#  Rscript 01_enmap_hyperspectral_workflow.R
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

# ============================================================
# 0) HELPERS
# ============================================================

msg <- function(...) cat(sprintf("[%s] %s\n",
                                 format(Sys.time(), "%F %T"),
                                 paste(..., collapse = " ")))

pick_band_by_wavelength <- function(wl, target_nm) {
  # Returns index of wavelength closest to target_nm
  which.min(abs(wl - target_nm))
}

clip_to_aoi <- function(r, aoi_vect_4326) {
  aoi_proj <- terra::project(aoi_vect_4326, terra::crs(r))
  aoi_proj <- terra::buffer(aoi_proj, width = 0)
  terra::mask(terra::crop(r, aoi_proj, snap = "out"), aoi_proj)
}

# ============================================================
# 1) PATHS + SETTINGS (Laptop + Bunya)
# ============================================================

default_root <- if (.Platform$OS.type == "windows") {
  "C:/Users/dorji/Documents/Phenomics/hyperspectral/enmap"
} else {
  file.path("hyperspectral", "enmap")
}

project_root <- Sys.getenv("ENMAP_ROOT", unset = default_root)

# Required: reflectance cube path (EnMAP L2A spectral image)
default_cube <- file.path(project_root, "data", "ENMAP_L2A_SPECTRAL_IMAGE.tif")
cube_path <- Sys.getenv("ENMAP_CUBE", unset = default_cube)

# Optional: cloud/quality mask (depends on your download contents)
mask_path <- Sys.getenv("ENMAP_MASK", unset = "")  # empty = no mask applied

default_out <- file.path(project_root, "outputs", "enmap_products")
out_dir <- Sys.getenv("ENMAP_OUTDIR", unset = default_out)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

default_tmp <- file.path(out_dir, "tmp")
tmp_dir <- Sys.getenv("ENMAP_TMPDIR", unset = default_tmp)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = tmp_dir)

threads <- as.integer(Sys.getenv("TERRA_THREADS", unset = "2"))
threads <- max(1L, threads)
terraOptions(threads = threads)

remove_tmp_dir <- as.logical(Sys.getenv("REMOVE_TMP_DIR", unset = "FALSE"))

# AOI bbox (WGS84) — edit for your site (example: Gatton)
aoi_sfc <- st_as_sfc(st_bbox(c(
  xmin = 152.25,
  xmax = 152.45,
  ymin = -27.65,
  ymax = -27.45
), crs = 4326))
aoi_sf <- st_sf(geometry = aoi_sfc)
aoi_vect_4326 <- terra::vect(aoi_sf)

# PCA settings (reduce memory)
pca_nbands_max <- 120   # sample this many bands for PCA (hyperspectral has lots)
pca_ncells     <- 15000 # sample pixels for PCA model fit

msg("Project root:", normalizePath(project_root, winslash = "/", mustWork = FALSE))
msg("Cube:", cube_path)
msg("Mask:", ifelse(mask_path == "", "(none)", mask_path))
msg("Out dir:", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
msg("Tmp dir:", normalizePath(tmp_dir, winslash = "/", mustWork = FALSE))
msg("terra threads:", threads)

if (!file.exists(cube_path)) stop("EnMAP cube not found: ", cube_path)

# ============================================================
# 2) LOAD + CLIP
# ============================================================

msg("Reading EnMAP cube (this may take time) ...")
cube <- rast(cube_path)

# Clip early to reduce memory
msg("Clipping cube to AOI ...")
cube_aoi <- clip_to_aoi(cube, aoi_vect_4326)

# Optional mask
if (mask_path != "" && file.exists(mask_path)) {
  msg("Applying mask:", mask_path)
  m <- rast(mask_path)
  m <- clip_to_aoi(m, aoi_vect_4326)
  m <- resample(m, cube_aoi[[1]], method = "near")
  
  # IMPORTANT:
  # Mask conventions differ by provider. Adjust this logic to your mask!
  # Option A: mask is 1 = clear, 0 = cloud -> keep where m==1
  # Option B: mask classes -> define "clear classes"
  #
  # Below assumes 1 = clear, 0 = not clear:
  clear <- m == 1
  cube_aoi <- mask(cube_aoi, clear, maskvalues = FALSE)
} else {
  msg("No mask applied (mask missing or ENMAP_MASK not set).")
}

# ============================================================
# 3) BAND/WAVELENGTH HANDLING
# ============================================================

# If wavelengths are embedded, terra sometimes stores them as "wavelength" in metadata.
# We try to retrieve them; if not available, we use band indices you can edit.
wl <- tryCatch({
  w <- wavelengths(cube_aoi)
  if (is.null(w) || length(w) == 0) NULL else as.numeric(w)
}, error = function(e) NULL)

if (!is.null(wl)) {
  msg("Wavelengths found. Range:", paste(range(wl, na.rm = TRUE), collapse = " - "), "nm")
  i_red   <- pick_band_by_wavelength(wl, 665)
  i_nir   <- pick_band_by_wavelength(wl, 842)
  i_re705 <- pick_band_by_wavelength(wl, 705)
  i_r531  <- pick_band_by_wavelength(wl, 531)
  i_r570  <- pick_band_by_wavelength(wl, 570)
} else {
  msg("No wavelengths metadata found. Using fallback band indices (EDIT if needed).")
  # Fallback indices (you MUST adjust to your EnMAP product band order if needed)
  i_red   <- 50
  i_nir   <- 90
  i_re705 <- 75
  i_r531  <- 35
  i_r570  <- 40
}

red   <- cube_aoi[[i_red]]
nir   <- cube_aoi[[i_nir]]
re705 <- cube_aoi[[i_re705]]
r531  <- cube_aoi[[i_r531]]
r570  <- cube_aoi[[i_r570]]

# ============================================================
# 4) HYPERSPECTRAL INDICES
# ============================================================

msg("Computing indices ...")

# NDVI
ndvi <- (nir - red) / (nir + red)

# NDRE (using ~705 nm red-edge)
ndre <- (nir - re705) / (nir + re705)

# PRI (Photochemical Reflectance Index)
# PRI = (R531 - R570) / (R531 + R570)
pri <- (r531 - r570) / (r531 + r570)

names(ndvi) <- "NDVI"
names(ndre) <- "NDRE"
names(pri)  <- "PRI"

# Save indices
writeRaster(ndvi, file.path(out_dir, "ENMAP_NDVI.tif"), overwrite = TRUE)
writeRaster(ndre, file.path(out_dir, "ENMAP_NDRE.tif"), overwrite = TRUE)
writeRaster(pri,  file.path(out_dir, "ENMAP_PRI.tif"),  overwrite = TRUE)
msg("Saved indices GeoTIFFs.")

# ============================================================
# 5) PCA (DIMENSION REDUCTION FOR MAPPING)
# ============================================================

msg("Preparing PCA (sampling pixels + bands) ...")

nb <- nlyr(cube_aoi)

# Choose a subset of bands (evenly spaced) to keep PCA lighter
k <- min(pca_nbands_max, nb)
band_idx <- unique(round(seq(1, nb, length.out = k)))
cube_sub <- cube_aoi[[band_idx]]

# Sample pixels for PCA model fit
samp <- spatSample(
  cube_sub,
  size = pca_ncells,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE,
  warn = FALSE
)

if (nrow(samp) < 1000) {
  msg("Warning: very few valid pixels for PCA (", nrow(samp), "). PCA may be unstable.")
}

# Fit PCA on sampled spectra
pca <- prcomp(samp, center = TRUE, scale. = TRUE)

# Project full raster (PC1-3) using PCA loadings
# We standardize layers using PCA center/scale before multiplying by rotation
msg("Projecting PC1-3 to rasters ...")
center <- pca$center
scalev <- pca$scale
rot <- pca$rotation[, 1:3, drop = FALSE]

# Standardize raster layers
cube_std <- (cube_sub - center) / scalev

# Compute PCs: PC = sum_i (layer_i * loading_i)
pc1 <- app(cube_std, fun = function(...) { as.matrix(cbind(...)) %*% rot[,1] })
pc2 <- app(cube_std, fun = function(...) { as.matrix(cbind(...)) %*% rot[,2] })
pc3 <- app(cube_std, fun = function(...) { as.matrix(cbind(...)) %*% rot[,3] })

pc <- c(pc1, pc2, pc3)
names(pc) <- c("PC1","PC2","PC3")

writeRaster(pc, file.path(out_dir, "ENMAP_PCA_PC1_PC3.tif"), overwrite = TRUE)
msg("Saved PCA raster (PC1-PC3).")

# ============================================================
# 6) AOI MEAN SPECTRUM (PORTFOLIO OUTPUT)
# ============================================================

msg("Extracting AOI mean spectrum ...")
mean_spec <- global(cube_aoi, fun = mean, na.rm = TRUE)
mean_spec <- as.numeric(mean_spec[1, ])

spec_df <- data.frame(
  band = seq_along(mean_spec),
  wavelength_nm = if (!is.null(wl)) wl else NA_real_,
  reflectance_mean = mean_spec
)

write.csv(spec_df, file.path(out_dir, "AOI_mean_spectrum.csv"), row.names = FALSE)
msg("Saved AOI_mean_spectrum.csv")

# ============================================================
# 7) QUICKLOOK PNGs (GOOD FOR README)
# ============================================================

png(file.path(out_dir, "quicklook_NDVI.png"), width = 1400, height = 1000, res = 150)
plot(ndvi, main = "EnMAP NDVI (AOI)")
dev.off()

png(file.path(out_dir, "quicklook_NDRE.png"), width = 1400, height = 1000, res = 150)
plot(ndre, main = "EnMAP NDRE (AOI)")
dev.off()

png(file.path(out_dir, "quicklook_PRI.png"), width = 1400, height = 1000, res = 150)
plot(pri, main = "EnMAP PRI (AOI)")
dev.off()

png(file.path(out_dir, "quicklook_PC1.png"), width = 1400, height = 1000, res = 150)
plot(pc[["PC1"]], main = "EnMAP PCA - PC1 (AOI)")
dev.off()

msg("Saved quicklook PNGs in:", out_dir)

# ============================================================
# 8) CLEAN UP
# ============================================================

terra::tmpFiles(remove = TRUE)

if (isTRUE(remove_tmp_dir)) {
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
  msg("Removed tmp directory:", tmp_dir)
} else {
  msg("Kept tmp directory (set REMOVE_TMP_DIR=TRUE to delete):", tmp_dir)
}

msg("DONE.")
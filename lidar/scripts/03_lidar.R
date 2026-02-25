

#!/usr/bin/env Rscript
# ============================================================
# LiDAR Crop Phenomics Workflow (lidR + terra) | Laptop + Bunya
# ------------------------------------------------------------
# Author: Sangay Dorji
#
# What this script does:
#  1) (Optional) Unzips LiDAR delivery ZIPs into a single folder
#  2) Builds a LAScatalog from LAZ/LAS tiles
#  3) Uses an existing DTM mosaic to normalize heights
#  4) Generates CHM + structural metrics (veg-friendly thresholds)
#  5) Mosaics chunk outputs into single rasters
#  6) Saves GeoTIFFs + quicklook PNGs (for GitHub README)
#
# GitHub notes:
#  - Do NOT commit raw LAZ/ZIPs. Commit code + small samples only.
#  - Use env vars for paths so the same script runs on laptop + HPC.
#
# Recommended env vars:
#  - LIDAR_ROOT   : project root folder
#  - LIDAR_ZIPDIR : folder containing delivery ZIPs (optional)
#  - LIDAR_DTM    : path to DTM GeoTIFF
#  - LIDAR_OUTDIR : output folder (optional)
#  - LIDAR_TMPDIR : temp folder (recommended on HPC scratch)
#  - TERRA_THREADS: number of threads for terra
#  - REMOVE_TMP_DIR: TRUE/FALSE (optional clean up)
#
# Example (Windows PowerShell):
#  $env:LIDAR_ROOT="C:\Users\dorji\Documents\Phenomics\lidar"
#  $env:LIDAR_ZIPDIR="C:\Users\dorji\Downloads\DATA_1018584\QLD Government\Point Clouds\AHD"
#  $env:LIDAR_DTM="C:\Users\dorji\Documents\Phenomics\lidar\data\lidar_mosaic\DTM_5m.tif"
#  $env:TERRA_THREADS="4"
#  Rscript 01_lidar_phenomics.R
#
# Example (Bunya bash):
#  export LIDAR_ROOT=/scratch/user/s4597131/phenomics/lidar
#  export LIDAR_ZIPDIR=/scratch/user/s4597131/phenomics/lidar_zips
#  export LIDAR_DTM=/scratch/user/s4597131/phenomics/lidar/data/lidar_mosaic/DTM_5m.tif
#  export LIDAR_TMPDIR=/scratch/user/s4597131/phenomics/tmp_lidar
#  export TERRA_THREADS=8
#  Rscript 01_lidar_phenomics.R
# ============================================================

suppressPackageStartupMessages({
  library(lidR)
  library(terra)
})

options(lidR.raster.backend = "terra")

# ============================================================
# 0) PATHS + SETTINGS (Laptop + HPC)
# ============================================================

msg <- function(...) cat(sprintf("[%s] %s\n",
                                 format(Sys.time(), "%F %T"),
                                 paste(..., collapse = " ")))

# --- Project root ---
default_root <- if (.Platform$OS.type == "windows") {
  "C:/Users/dorji/Documents/Phenomics/lidar"
} else {
  file.path("lidar") # safe default if not set; override via env var
}
project_root <- Sys.getenv("LIDAR_ROOT", unset = default_root)

# --- Optional ZIP folder (only needed if do_unzip = TRUE) ---
# If you don't set it, the script won't try to unzip.
zip_dir <- Sys.getenv("LIDAR_ZIPDIR", unset = "")

# --- Data folder containing LAZ/LAS tiles (inside project) ---
laz_dir <- file.path(project_root, "data")

# --- DTM mosaic path ---
# You can override with env var LIDAR_DTM
default_dtm <- file.path(project_root, "data", "lidar_mosaic", "DTM_5m.tif")
dtm_path <- Sys.getenv("LIDAR_DTM", unset = default_dtm)

# --- Output folder ---
default_out <- file.path(project_root, "outputs", "lidar_products_gatton")
out_dir <- Sys.getenv("LIDAR_OUTDIR", unset = default_out)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Temp folder (recommended on HPC scratch) ---
default_tmp <- file.path(out_dir, "tmp")
tmp_dir <- Sys.getenv("LIDAR_TMPDIR", unset = default_tmp)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = tmp_dir)

# --- Threads ---
threads <- as.integer(Sys.getenv("TERRA_THREADS", unset = "2"))
threads <- max(1L, threads)
terraOptions(threads = threads)

# --- Clean-up toggle ---
remove_tmp_dir <- as.logical(Sys.getenv("REMOVE_TMP_DIR", unset = "FALSE"))

# --- Processing settings ---
do_unzip <- FALSE  # set TRUE only once if you need to extract delivery ZIPs
res_m <- 2         # 2 m for veg phenomics; use 5 m if noisy
chunk_size_m   <- 1000
chunk_buffer_m <- 30
overwrite <- TRUE

msg("Project root:", normalizePath(project_root, winslash = "/", mustWork = FALSE))
msg("LAZ/LAS dir:", normalizePath(laz_dir, winslash = "/", mustWork = FALSE))
msg("DTM path:", dtm_path)
msg("Output dir:", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
msg("Temp dir:", normalizePath(tmp_dir, winslash = "/", mustWork = FALSE))
msg("terra threads:", threads)

# ============================================================
# 1) OPTIONAL: UNZIP DELIVERY FILES (run once)
# ============================================================

if (isTRUE(do_unzip)) {
  if (zip_dir == "") stop("do_unzip=TRUE but LIDAR_ZIPDIR not set. Provide zip_dir.")
  msg("Unzipping LiDAR ZIPs from:", zip_dir)
  
  dir.create(laz_dir, recursive = TRUE, showWarnings = FALSE)
  
  zips <- list.files(zip_dir, pattern = "\\.zip$", full.names = TRUE, recursive = TRUE)
  if (length(zips) == 0) stop("No .zip files found in: ", zip_dir)
  
  for (zf in zips) {
    msg("  ->", basename(zf))
    unzip(zf, exdir = laz_dir, overwrite = TRUE)
  }
  msg("Unzip complete. Extracted into:", normalizePath(laz_dir, winslash = "/", mustWork = FALSE))
}

# ============================================================
# 2) BASIC CHECKS
# ============================================================

if (!dir.exists(laz_dir)) stop("LAZ/LAS folder not found: ", laz_dir)
if (!file.exists(dtm_path)) stop("DTM not found: ", dtm_path)

las_files <- list.files(
  laz_dir,
  pattern = "\\.(laz|las)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)
if (length(las_files) == 0) stop("No LAZ/LAS files found in: ", laz_dir)

msg("Found", length(las_files), "LiDAR tiles.")
dtm <- rast(dtm_path)

# ============================================================
# 3) BUILD LASCATALOG
# ============================================================

ctg <- readLAScatalog(laz_dir)
opt_chunk_size(ctg)   <- chunk_size_m
opt_chunk_buffer(ctg) <- chunk_buffer_m
opt_progress(ctg)     <- TRUE
# Optional: store intermediate chunk outputs (helpful for debugging)
# opt_output_files(ctg) <- file.path(out_dir, "chunks", "{ID}")

msg("Catalog ready. Chunk size:", chunk_size_m, "m | buffer:", chunk_buffer_m, "m")

# ============================================================
# 4) CHUNK PROCESSOR (called by catalog_apply)
# ============================================================

process_chunk <- function(cluster, dtm, res_m) {
  
  las <- tryCatch(readLAS(cluster), error = function(e) NULL)
  if (is.null(las) || is.empty(las)) return(NULL)
  
  # Crop DTM to chunk extent for speed
  dtm_c <- tryCatch(crop(dtm, ext(las)), error = function(e) NULL)
  if (is.null(dtm_c)) return(NULL)
  
  # Normalize heights (Z becomes above-ground height)
  las_n <- tryCatch(normalize_height(las, dtm_c), error = function(e) NULL)
  if (is.null(las_n) || is.empty(las_n)) return(NULL)
  
  # Remove negative heights
  las_n <- filter_poi(las_n, Z >= 0)
  if (is.empty(las_n)) return(NULL)
  
  # Canopy surface (after normalization, DSM ~ CHM)
  dsm <- rasterize_canopy(
    las_n,
    res = res_m,
    algorithm = p2r(subcircle = min(0.2, res_m / 2))
  )
  chm_r <- if (inherits(dsm, "SpatRaster")) dsm else rast(dsm)
  
  # Veg-friendly structural metrics (vegetables often < 1.5 m)
  metrics_fun <- ~list(
    zmean   = mean(Z),
    zsd     = sd(Z),
    zq75    = quantile(Z, 0.75),
    zq90    = quantile(Z, 0.90),
    zq95    = quantile(Z, 0.95),
    cover02 = sum(Z > 0.2) / length(Z),
    cover05 = sum(Z > 0.5) / length(Z),
    cover10 = sum(Z > 1.0) / length(Z)
  )
  
  met <- grid_metrics(las_n, metrics_fun, res = res_m)
  met_r <- if (inherits(met, "SpatRaster")) met else rast(met)
  names(met_r) <- c("zmean","zsd","zq75","zq90","zq95","cover02","cover05","cover10")
  
  list(chm = chm_r, met = met_r)
}

# ============================================================
# 5) RUN PROCESSING ACROSS CATALOG
# ============================================================

msg("Running catalog_apply() ...")
out_list <- catalog_apply(ctg, process_chunk, dtm = dtm, res_m = res_m)

out_list <- Filter(Negate(is.null), out_list)
if (length(out_list) == 0) stop("No chunk outputs. Check CRS, DTM overlap, and input coverage.")

msg("Chunks processed:", length(out_list))

# ============================================================
# 6) MOSAIC CHUNK OUTPUTS
# ============================================================

msg("Mosaicking CHM ...")
chm_list <- lapply(out_list, \(x) x$chm)
chm_list <- Filter(\(x) inherits(x, "SpatRaster"), chm_list)
if (length(chm_list) == 0) stop("No CHM rasters found.")

chm_m <- do.call(terra::mosaic, c(chm_list, list(fun = "max")))

msg("Mosaicking metrics ...")
met_list <- lapply(out_list, \(x) x$met)
met_list <- Filter(\(x) inherits(x, "SpatRaster"), met_list)
if (length(met_list) == 0) stop("No metric rasters found.")

mosaic_mean <- function(layer) {
  lst <- lapply(met_list, \(r) r[[layer]])
  do.call(terra::mosaic, c(lst, list(fun = "mean")))
}
mosaic_max <- function(layer) {
  lst <- lapply(met_list, \(r) r[[layer]])
  do.call(terra::mosaic, c(lst, list(fun = "max")))
}

# Mean mosaics
zmean_m   <- mosaic_mean("zmean")
zsd_m     <- mosaic_mean("zsd")
cover02_m <- mosaic_mean("cover02")
cover05_m <- mosaic_mean("cover05")
cover10_m <- mosaic_mean("cover10")

# Max mosaics
zq75_m <- mosaic_max("zq75")
zq90_m <- mosaic_max("zq90")
zq95_m <- mosaic_max("zq95")

met_m <- c(zmean_m, zsd_m, zq75_m, zq90_m, zq95_m, cover02_m, cover05_m, cover10_m)
names(met_m) <- c("zmean","zsd","zq75","zq90","zq95","cover02","cover05","cover10")

# ============================================================
# 7) SAVE OUTPUTS + QUICKLOOKS
# ============================================================

chm_out <- file.path(out_dir, sprintf("CHM_res%sm.tif", res_m))
met_out <- file.path(out_dir, sprintf("metrics_res%sm.tif", res_m))

writeRaster(chm_m, chm_out, overwrite = overwrite)
writeRaster(met_m, met_out, overwrite = overwrite)

png(file.path(out_dir, "CHM.png"), width = 1400, height = 1000, res = 150)
plot(chm_m, main = sprintf("CHM (above-ground height), res=%sm", res_m))
dev.off()

png(file.path(out_dir, "metrics_zq95.png"), width = 1400, height = 1000, res = 150)
plot(met_m[["zq95"]], main = "zq95 (95th percentile height)")
dev.off()

png(file.path(out_dir, "metrics_cover05.png"), width = 1400, height = 1000, res = 150)
plot(met_m[["cover05"]], main = "cover05 (fraction of points > 0.5 m)")
dev.off()

msg("CHM saved:", chm_out)
msg("Metrics saved:", met_out)
msg("Figures saved in:", out_dir)
msg("DONE.")

# ============================================================
# 8) CLEAN UP TEMP FILES
# ============================================================

terra::tmpFiles(remove = TRUE)

if (isTRUE(remove_tmp_dir)) {
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
  msg("Removed tmp directory:", tmp_dir)
} else {
  msg("Kept tmp directory (set REMOVE_TMP_DIR=TRUE to delete):", tmp_dir)
}


# Prepare DTM from lidar
#!/usr/bin/env Rscript
# ============================================================
# Build DTM from LAScatalog (Professional Workflow)
# ------------------------------------------------------------
# - Processes all LiDAR tiles using chunking
# - Writes temporary chunk rasters
# - Mosaics into a single DTM
# - Automatically removes tmp folder after success
# ============================================================

suppressPackageStartupMessages({
  library(lidR)
  library(terra)
})

options(lidR.raster.backend = "terra")

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

# Change only this root path for laptop vs Bunya
project_root <- "~/Phenomics/lidar"

laz_dir <- file.path(project_root, "data")                 # 104 .laz files
out_dir <- file.path(project_root, "data", "lidar_mosaic")
tmp_dir <- file.path(out_dir, "tmp")

res_m <- 5
chunk_size_m   <- 1000
chunk_buffer_m <- 30

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir,  recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(sprintf("[%s] %s\n",
                                 format(Sys.time(), "%F %T"),
                                 paste(..., collapse = " ")))

# ------------------------------------------------------------
# 1) Create LAScatalog
# ------------------------------------------------------------
ctg <- readLAScatalog(laz_dir)
opt_chunk_size(ctg)   <- chunk_size_m
opt_chunk_buffer(ctg) <- chunk_buffer_m
opt_progress(ctg)     <- TRUE

# IMPORTANT: write chunk rasters to temp folder
opt_output_files(ctg) <- file.path(tmp_dir, "dtm_{XLEFT}_{YBOTTOM}")

msg("Processing LAScatalog...")
msg("Resolution:", res_m, "m")

# ------------------------------------------------------------
# 2) Rasterize terrain (writes chunk GeoTIFFs to tmp/)
# ------------------------------------------------------------
dtm_out <- rasterize_terrain(ctg, res = res_m, algorithm = tin())

# ------------------------------------------------------------
# 3) Mosaic chunk rasters
# ------------------------------------------------------------
dtm_files <- list.files(tmp_dir, pattern = "\\.tif$", 
                        full.names = TRUE, recursive = TRUE)

if (length(dtm_files) == 0) {
  stop("No DTM GeoTIFFs written to tmp_dir:\n", tmp_dir)
}

msg("Mosaicking", length(dtm_files), "chunks...")

dtm_list <- lapply(dtm_files, terra::rast)

# Mean in overlap zones (safe for terrain)
dtm <- do.call(terra::mosaic, c(dtm_list, list(fun = "mean")))

dtm_path <- file.path(out_dir, sprintf("DTM_%sm.tif", res_m))
writeRaster(dtm, dtm_path, overwrite = TRUE)

msg("DTM saved:", dtm_path)

# ------------------------------------------------------------
# 4) Clean up temporary files
# ------------------------------------------------------------
msg("Cleaning temporary files...")

try({
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
}, silent = TRUE)

msg("Temporary folder removed.")
msg("DONE.")

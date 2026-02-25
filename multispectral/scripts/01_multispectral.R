# ============================================================
# Sentinel-2 L2A Monthly Multispectral Workflow (R | Laptop + Bunya)
# ------------------------------------------------------------
# Author: Sangay Dorji
# Sentinel raw has not been saved due to file size
# Workflow:
#  1) STAC search Sentinel-2 L2A over AOI (WGS84)
#  2) Read COG assets directly (no manual downloads)
#  3) Cloud-mask using SCL
#  4) Compute NDVI + NDRE per scene
#  5) Create MONTHLY median composites (20 m grid)
#  6) Save GeoTIFFs + AOI mean time series CSV + plot + summary
#
# Key improvements:
#  • Works on Windows laptop and Linux/HPC without editing code:
#      - Set OUTDIR via env var S2_OUTDIR (recommended on Bunya)
#      - Otherwise uses a sensible default (Windows path if on Windows)
#  • Safer temp handling:
#      - Uses out_dir/tmp by default
#      - Uses scratch temp on HPC if S2_TMPDIR is set
#      - Cleans terra tmp files; optionally removes tmp folder
#  • Robust per-scene processing:
#      - Skips scenes that fail to read/download
#      - Avoids crashing entire run due to one bad scene
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(rstac)
})

# ============================================================
# 1) SETTINGS
# ============================================================

year      <- 2023
max_items <- 120

# ---- Output directory ----
# Preferred: set env var S2_OUTDIR (works everywhere)
#   Windows (PowerShell):
#     $env:S2_OUTDIR="C:\Users\dorji\Documents\Phenomics\multispectral\outputs"
#   Bunya (bash):
#     export S2_OUTDIR=/scratch/user/s4597131/phenomics/outputs/multispectral
#
# If not set, choose a platform-specific default.
default_outdir <- if (.Platform$OS.type == "windows") {
  "C:/Users/dorji/Documents/Phenomics/multispectral/outputs"
} else {
  file.path("outputs", "sentinel2_demo")
}
out_dir <- Sys.getenv("S2_OUTDIR", unset = default_outdir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Temp directory ----
# Preferred on HPC: set S2_TMPDIR to a fast scratch location.
default_tmpdir <- file.path(out_dir, "tmp")
tmp_dir <- Sys.getenv("S2_TMPDIR", unset = default_tmpdir)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = tmp_dir)

# ---- Threads ----
# Set TERRA_THREADS on HPC to match cpus-per-task, otherwise default to 2.
threads <- as.integer(Sys.getenv("TERRA_THREADS", unset = "2"))
threads <- max(1L, threads)
terraOptions(threads = threads)

# ---- Remove tmp folder at end? ----
# Keep FALSE during development (useful for debugging).
remove_tmp_dir <- as.logical(Sys.getenv("REMOVE_TMP_DIR", unset = "FALSE"))

# ---------------------------
# AOI (WGS84 bbox) - edit bbox
# ---------------------------
# Example: Gatton farm area
aoi_sfc <- st_as_sfc(st_bbox(c(
  xmin = 152.25,
  xmax = 152.45,
  ymin = -27.65,
  ymax = -27.45
), crs = 4326))

aoi_sf        <- st_sf(geometry = aoi_sfc) # for STAC search
aoi_vect_4326 <- terra::vect(aoi_sf)       # for terra crop/mask

# ---------------------------
# SCL cloud mask classes (good pixels)
# ---------------------------
good_scl <- c(4, 5, 6, 7)  # vegetation, bare soil, water, unclassified

# ============================================================
# 2) HELPERS
# ============================================================

msg <- function(...) {
  cat(sprintf("[%s] %s\n",
              format(Sys.time(), "%F %T"),
              paste(..., collapse = " ")))
}

get_asset_href <- function(feature, asset_name) {
  href <- feature$assets[[asset_name]]$href
  if (is.null(href)) stop("Missing asset: ", asset_name)
  href
}

# CRS-safe clip: reproject AOI to raster CRS before crop/mask
clip_to_aoi <- function(r, aoi_wgs_vect) {
  aoi_proj <- terra::project(aoi_wgs_vect, terra::crs(r))
  aoi_proj <- terra::buffer(aoi_proj, width = 0) # avoid edge precision issues
  terra::mask(terra::crop(r, aoi_proj, snap = "out"), aoi_proj)
}

safe_rast <- function(href) {
  # Read a raster safely; return NULL if it fails
  tryCatch(rast(href), error = function(e) NULL)
}

# ============================================================
# 3) STAC SEARCH
# ============================================================

datetime_range <- sprintf(
  "%d-01-01T00:00:00Z/%d-12-31T23:59:59Z",
  year, year
)

msg("Output dir:", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
msg("Temp dir:", normalizePath(tmp_dir, winslash = "/", mustWork = FALSE))
msg("terra threads:", threads)
msg("Searching STAC Sentinel-2 L2A:", datetime_range)

items <- stac("https://earth-search.aws.element84.com/v1") |>
  stac_search(
    collections = "sentinel-2-l2a",
    intersects  = aoi_sf,
    datetime    = datetime_range,
    limit       = max_items
  ) |>
  post_request()

n_items <- length(items$features)
msg("Scenes matched:", n_items)
if (n_items == 0) stop("No scenes found. Try expanding AOI or increasing max_items.")

# Scene dates -> months
item_dt   <- vapply(items$features, \(f) f$properties$datetime, character(1))
item_date <- as.Date(substr(item_dt, 1, 10))
item_ym   <- format(item_date, "%Y-%m")
months    <- sort(unique(item_ym))
msg("Months detected:", paste(months, collapse = ", "))

# Time-series output table
monthly_ts <- data.frame(
  year_month = character(),
  ndvi_mean  = numeric(),
  ndre_mean  = numeric()
)

# ============================================================
# 4) PROCESS MONTH BY MONTH
# ============================================================

for (ym in months) {
  
  idx <- which(item_ym == ym)
  msg("Processing month:", ym, "| scenes:", length(idx))
  
  ndvi_list <- list()
  ndre_list <- list()
  
  for (j in idx) {
    
    f <- items$features[[j]]
    
    # Read assets (COG URLs)
    red_href <- get_asset_href(f, "red")
    nir_href <- get_asset_href(f, "nir")
    re1_href <- get_asset_href(f, "rededge1")
    scl_href <- get_asset_href(f, "scl")
    
    red <- safe_rast(red_href)
    nir <- safe_rast(nir_href)
    re1 <- safe_rast(re1_href)
    scl <- safe_rast(scl_href)
    
    if (any(vapply(list(red, nir, re1, scl), is.null, logical(1)))) {
      msg("  Skipping scene (failed to read one or more assets).")
      next
    }
    
    # Clip to AOI (CRS-safe)
    red <- clip_to_aoi(red, aoi_vect_4326)
    nir <- clip_to_aoi(nir, aoi_vect_4326)
    re1 <- clip_to_aoi(re1, aoi_vect_4326)
    scl <- clip_to_aoi(scl, aoi_vect_4326)
    
    # Use 20 m grid = re1; resample red/nir to 20 m, SCL nearest
    red20 <- resample(red, re1, method = "bilinear")
    nir20 <- resample(nir, re1, method = "bilinear")
    scl20 <- resample(scl, re1, method = "near")
    
    # Cloud mask using SCL good classes
    good <- scl20 %in% good_scl
    red20 <- mask(red20, good, maskvalues = FALSE)
    nir20 <- mask(nir20, good, maskvalues = FALSE)
    re1   <- mask(re1,   good, maskvalues = FALSE)
    
    # Indices
    ndvi <- (nir20 - red20) / (nir20 + red20)
    ndre <- (nir20 - re1)   / (nir20 + re1)
    
    ndvi_list[[length(ndvi_list) + 1]] <- ndvi
    ndre_list[[length(ndre_list) + 1]] <- ndre
  }
  
  if (length(ndvi_list) == 0 || length(ndre_list) == 0) {
    msg("  No valid scenes for month:", ym, "- skipping composite.")
    next
  }
  
  # Monthly median composite
  ndvi_med <- app(rast(ndvi_list), median, na.rm = TRUE)
  ndre_med <- app(rast(ndre_list), median, na.rm = TRUE)
  
  names(ndvi_med) <- paste0("NDVI_", ym)
  names(ndre_med) <- paste0("NDRE_", ym)
  
  # Save rasters
  out_ndvi <- file.path(out_dir, paste0("NDVI_", ym, "_20m.tif"))
  out_ndre <- file.path(out_dir, paste0("NDRE_", ym, "_20m.tif"))
  writeRaster(ndvi_med, out_ndvi, overwrite = TRUE)
  writeRaster(ndre_med, out_ndre, overwrite = TRUE)
  
  # AOI mean time series
  ndvi_mean <- as.numeric(global(ndvi_med, mean, na.rm = TRUE)[1, 1])
  ndre_mean <- as.numeric(global(ndre_med, mean, na.rm = TRUE)[1, 1])
  
  monthly_ts <- rbind(monthly_ts, data.frame(
    year_month = ym,
    ndvi_mean  = ndvi_mean,
    ndre_mean  = ndre_mean
  ))
  
  msg("Saved:", basename(out_ndvi), "and", basename(out_ndre))
}

# ============================================================
# 5) SAVE TIME SERIES + PLOT
# ============================================================

ts_csv <- file.path(out_dir, paste0("AOI_monthly_timeseries_", year, ".csv"))
write.csv(monthly_ts, ts_csv, row.names = FALSE)
msg("Saved time series:", basename(ts_csv))

png(file.path(out_dir, paste0("AOI_monthly_timeseries_", year, ".png")),
    width = 1200, height = 700, res = 150)

par(mar = c(5, 5, 3, 5))
x <- seq_len(nrow(monthly_ts))

plot(x, monthly_ts$ndvi_mean, type = "b", pch = 16,
     xlab = "Month", ylab = "NDVI (AOI mean)",
     main = paste0("AOI Monthly NDVI & NDRE (", year, ")"))
axis(1, at = x, labels = monthly_ts$year_month, las = 2, cex.axis = 0.8)

par(new = TRUE)
plot(x, monthly_ts$ndre_mean, type = "b", pch = 17,
     axes = FALSE, xlab = "", ylab = "")
axis(4)
mtext("NDRE (AOI mean)", side = 4, line = 3)

legend("topleft", legend = c("NDVI", "NDRE"), pch = c(16, 17), bty = "n")
dev.off()

# Simple phenology metric: peak NDVI month
if (nrow(monthly_ts) > 0) {
  peak_i     <- which.max(monthly_ts$ndvi_mean)
  peak_month <- monthly_ts$year_month[peak_i]
  writeLines(paste("Peak NDVI month:", peak_month),
             con = file.path(out_dir, "phenology_summary.txt"))
  msg("Peak NDVI month:", peak_month)
} else {
  msg("No monthly stats computed (monthly_ts is empty).")
}

msg("DONE. Outputs in:", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

# ============================================================
# 6) CLEAN UP TEMP FILES
# ============================================================

terra::tmpFiles(remove = TRUE)

if (isTRUE(remove_tmp_dir)) {
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
  msg("Removed tmp directory:", tmp_dir)
} else {
  msg("Kept tmp directory (set REMOVE_TMP_DIR=TRUE to delete):", tmp_dir)
}

#!/usr/bin/env Rscript
# ============================================================
# Phenomics RF Demo (Gatton, QLD): Predict LiDAR canopy trait from Sentinel-2
# ------------------------------------------------------------
# Author: Sangay Dorji
#
# What this script does:
#   1) Loads Sentinel-2 monthly NDVI + NDRE composites (GeoTIFFs)
#   2) Loads LiDAR target (zq95 from metrics_res2m.tif OR CHM_res2m.tif)
#   3) Aligns CRS + extent (GDA94 MGA56 vs WGS84 UTM56S handled)
#   4) Samples pixels to build training table (no field data required)
#   5) Trains Random Forest regression with K-fold CV
#   6) Predicts full-area trait map + residual map
#   7) Exports CSV metrics, variable importance, GeoTIFFs, and PNG quicklooks
#
# Works on:
#   - Laptop (Windows)
#   - Bunya / Linux HPC
#
# Folder assumptions (recommended):
#   PHENO_ROOT/
#     multispectral/outputs/   (your NDVI/NDRE monthly GeoTIFFs)
#     lidar/outputs/lidar_products_gatton/  (CHM_res2m.tif, metrics_res2m.tif)
#     modelling/outputs/       (created by script)
#
# Optional environment variables:
#   PHENO_ROOT    : project root path
#   S2_DIR        : override Sentinel output directory
#   LIDAR_DIR     : override LiDAR product directory
#   PHENO_OUTDIR  : override output directory
#   PHENO_TMPDIR  : override temp directory
#   TERRA_THREADS : terra threads (e.g., 4 laptop, 8 HPC)
#   TARGET        : "zq95" (default) or "chm"
#   N_SAMPLES     : number of training samples (default 15000)
#   K_FOLDS       : CV folds (default 5)
#   NUM_TREES     : RF trees (default 700)
#   MTRY_FRAC     : mtry fraction (default 0.35)
#   MIN_NODE      : min node size (default 5)
#
# Example (Windows PowerShell):
#   $env:PHENO_ROOT="C:\Users\dorji\Documents\Phenomics"
#   $env:TERRA_THREADS="4"
#   $env:TARGET="zq95"
#   Rscript 03_rf_predict_lidar_trait_from_s2.R
#
# Example (Bunya bash):
#   export PHENO_ROOT=/scratch/user/s4597131/phenomics
#   export TERRA_THREADS=8
#   export TARGET=zq95
#   Rscript 03_rf_predict_lidar_trait_from_s2.R
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(ranger)
})

msg <- function(...) cat(sprintf("[%s] %s\n",
                                 format(Sys.time(), "%F %T"),
                                 paste(..., collapse = " ")))

# ============================================================
# 0) PATHS + SETTINGS
# ============================================================

default_root <- if (.Platform$OS.type == "windows") {
  "C:/Users/dorji/Documents/Phenomics"
} else {
  "phenomics"
}
project_root <- Sys.getenv("PHENO_ROOT", unset = default_root)

s2_dir <- Sys.getenv("S2_DIR", unset = file.path(project_root, "multispectral", "outputs"))
lidar_dir <- Sys.getenv("LIDAR_DIR", unset = file.path(project_root, "lidar", "outputs", "lidar_products_gatton"))

out_dir <- Sys.getenv("PHENO_OUTDIR", unset = file.path(project_root, "modelling", "outputs"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tmp_dir <- Sys.getenv("PHENO_TMPDIR", unset = file.path(out_dir, "tmp"))
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

threads <- max(1L, as.integer(Sys.getenv("TERRA_THREADS", unset = "2")))
terraOptions(tempdir = tmp_dir, threads = threads)

use_target <- tolower(Sys.getenv("TARGET", unset = "zq95"))  # "zq95" or "chm"
n_samples  <- as.integer(Sys.getenv("N_SAMPLES", unset = "15000"))
k_folds    <- as.integer(Sys.getenv("K_FOLDS", unset = "5"))
num_trees  <- as.integer(Sys.getenv("NUM_TREES", unset = "700"))
mtry_frac  <- as.numeric(Sys.getenv("MTRY_FRAC", unset = "0.35"))
min_node   <- as.integer(Sys.getenv("MIN_NODE", unset = "5"))

msg("Project root:", normalizePath(project_root, winslash="/", mustWork=FALSE))
msg("S2 dir:", normalizePath(s2_dir, winslash="/", mustWork=FALSE))
msg("LiDAR dir:", normalizePath(lidar_dir, winslash="/", mustWork=FALSE))
msg("Out dir:", normalizePath(out_dir, winslash="/", mustWork=FALSE))
msg("Tmp dir:", normalizePath(tmp_dir, winslash="/", mustWork=FALSE))
msg("terra threads:", threads)
msg("Target:", use_target)
msg("Sampling:", n_samples, "| CV folds:", k_folds)

# ============================================================
# 1) LOAD SENTINEL-2 PREDICTORS (NDVI + NDRE monthly)
# ============================================================

ndvi_files <- list.files(s2_dir, pattern = "NDVI_monthly_median_.*_20m\\.tif$", full.names = TRUE, ignore.case = TRUE)
ndre_files <- list.files(s2_dir, pattern = "NDRE_monthly_median_.*_20m\\.tif$", full.names = TRUE, ignore.case = TRUE)

if (length(ndvi_files) == 0) stop("No NDVI GeoTIFFs found in: ", s2_dir)
if (length(ndre_files) == 0) stop("No NDRE GeoTIFFs found in: ", s2_dir)

ndvi_files <- sort(ndvi_files)
ndre_files <- sort(ndre_files)

msg("Found NDVI:", length(ndvi_files))
msg("Found NDRE:", length(ndre_files))

ndvi <- rast(ndvi_files)
ndre <- rast(ndre_files)

# Name layers WITHOUT hyphens (important!)
ndvi_month <- sub(".*_(\\d{4}-\\d{2})_20m.*", "\\1", basename(ndvi_files))
ndre_month <- sub(".*_(\\d{4}-\\d{2})_20m.*", "\\1", basename(ndre_files))
ndvi_month <- gsub("-", "_", ndvi_month)
ndre_month <- gsub("-", "_", ndre_month)

names(ndvi) <- paste0("NDVI_", ndvi_month)
names(ndre) <- paste0("NDRE_", ndre_month)

# Align all predictors to a single reference grid to avoid geometry mismatches
ref_s2 <- ndvi[[1]]
ndvi_a <- resample(ndvi, ref_s2, method = "bilinear")
ndre_a <- resample(ndre, ref_s2, method = "bilinear")

e_s2 <- intersect(ext(ndvi_a), ext(ndre_a))
if (is.null(e_s2)) stop("NDVI and NDRE do not overlap (unexpected). Check inputs.")
ndvi_a <- crop(ndvi_a, e_s2)
ndre_a <- crop(ndre_a, e_s2)

pred <- c(ndvi_a, ndre_a)
names(pred) <- make.names(names(pred), unique = TRUE)

if (!hasValues(pred[[1]])) stop("Predictor stack has no values after alignment. Check NDVI/NDRE grids.")
msg("Predictor layers:", nlyr(pred))

# ============================================================
# 2) LOAD LiDAR TARGET (zq95 from metrics OR CHM)
# ============================================================

metrics_path <- file.path(lidar_dir, "metrics_res2m.tif")
chm_path     <- file.path(lidar_dir, "CHM_res2m.tif")

if (use_target == "zq95") {
  if (!file.exists(metrics_path)) stop("Missing metrics_res2m.tif in: ", lidar_dir)
  met <- rast(metrics_path)
  if (!("zq95" %in% names(met))) stop("zq95 not found in metrics_res2m.tif. Layers: ", paste(names(met), collapse=", "))
  y <- met[["zq95"]]
  names(y) <- "zq95"
} else if (use_target == "chm") {
  if (!file.exists(chm_path)) stop("Missing CHM_res2m.tif in: ", lidar_dir)
  y <- rast(chm_path)
  names(y) <- "chm"
} else {
  stop("TARGET must be 'zq95' or 'chm'")
}

if (!hasValues(y)) stop("LiDAR target raster has no values: ", use_target)
msg("Loaded LiDAR target:", names(y))

# ============================================================
# 3) ALIGN CRS + OVERLAP (GDA94 MGA56 vs WGS84 UTM56S)
#    - Project Sentinel predictors to LiDAR CRS
#    - Resample LiDAR target to Sentinel grid (in LiDAR CRS)
# ============================================================

pred_p <- project(pred, crs(y), method = "bilinear")
names(pred_p) <- make.names(names(pred_p), unique = TRUE)

ref <- pred_p[[1]]

y20 <- resample(y, ref, method = "bilinear")

if (!hasValues(y20)) {
  stop(
    "Resampled LiDAR target has no values.\n",
    "CRS(pred): ", crs(pred[[1]]), "\n",
    "CRS(y): ", crs(y), "\n"
  )
}

e <- intersect(ext(ref), ext(y20))
if (is.null(e)) stop("No overlap between predictors and target after projection.\n",
                     "ext(ref): ", as.character(ext(ref)), "\n",
                     "ext(y20): ", as.character(ext(y20)))

pred_p <- crop(pred_p, e)
y20    <- crop(y20,    e)

mask_ref <- !is.na(y20)

pred_m <- mask(pred_p, mask_ref)
y_m    <- mask(y20,    mask_ref)

# Ensure names match model training table later
names(pred_m) <- make.names(names(pred_m), unique = TRUE)

if (!hasValues(pred_m) || !hasValues(y_m)) stop("Aligned rasters have no values (check NAs/coverage).")
msg("Alignment OK. Final predictors:", nlyr(pred_m))

# ============================================================
# 4) BUILD TRAINING TABLE (AUTO-SAMPLING)
# ============================================================

msg("Sampling", n_samples, "pixels ...")

stack_all <- c(y_m, pred_m)

samp <- spatSample(stack_all, size = n_samples, method = "random",
                   na.rm = TRUE, as.df = TRUE, warn = FALSE)

if (nrow(samp) < 2000) stop("Too few samples after NA removal: ", nrow(samp))

colnames(samp)[1] <- "target"
samp <- samp[complete.cases(samp), , drop = FALSE]

# Sanitize names for ranger formula AND keep them aligned with raster layer names
names(samp) <- make.names(names(samp), unique = TRUE)
names(samp)[1] <- "target"

# Make sure the predictor names used in samp match pred_m
pred_names <- names(pred_m)
samp_names <- setdiff(names(samp), "target")

# If they differ, force samp predictor order/names to match rasters
# (this prevents "independent variables not found" at prediction)
if (!all(samp_names %in% pred_names)) {
  missing_in_raster <- setdiff(samp_names, pred_names)
  stop("Some training predictors are missing in raster stack pred_m:\n", paste(missing_in_raster, collapse = ", "))
}

# Reorder samp columns: target first, then predictors in raster order
samp <- samp[, c("target", pred_names), drop = FALSE]

msg("Training rows used:", nrow(samp))
write.csv(samp, file.path(out_dir, "training_table_sample.csv"), row.names = FALSE)

# ============================================================
# 5) RANDOM FOREST REGRESSION + K-FOLD CV
# ============================================================

set.seed(42)
fold_id <- sample(rep(1:k_folds, length.out = nrow(samp)))

p <- ncol(samp) - 1
mtry <- max(1, floor(mtry_frac * p))

msg("RF settings:", "trees=", num_trees, "mtry=", mtry, "min.node.size=", min_node)

cv <- vector("list", k_folds)

for (k in 1:k_folds) {
  tr <- samp[fold_id != k, , drop = FALSE]
  te <- samp[fold_id == k, , drop = FALSE]
  
  rf <- ranger(
    target ~ .,
    data = tr,
    num.trees = num_trees,
    mtry = mtry,
    min.node.size = min_node,
    importance = "permutation",
    seed = 42
  )
  
  pr <- predict(rf, te)$predictions
  rmse <- sqrt(mean((pr - te$target)^2))
  r2   <- 1 - sum((pr - te$target)^2) / sum((te$target - mean(te$target))^2)
  
  cv[[k]] <- data.frame(fold = k, rmse = rmse, r2 = r2, n_test = nrow(te))
  msg("  Fold", k, "RMSE=", round(rmse, 3), "R2=", round(r2, 3))
}

cv_df <- do.call(rbind, cv)
write.csv(cv_df, file.path(out_dir, "cv_metrics.csv"), row.names = FALSE)

# Train final model
rf_final <- ranger(
  target ~ .,
  data = samp,
  num.trees = num_trees,
  mtry = mtry,
  min.node.size = min_node,
  importance = "permutation",
  seed = 42
)

saveRDS(rf_final, file.path(out_dir, "rf_model.rds"))

imp <- sort(rf_final$variable.importance, decreasing = TRUE)
imp_df <- data.frame(variable = names(imp), importance = as.numeric(imp))
write.csv(imp_df, file.path(out_dir, "variable_importance.csv"), row.names = FALSE)

# ============================================================
# 6) PREDICT FULL MAP (names MUST match)
# ============================================================

# Ensure raster names match model variable names
model_vars <- rf_final$forest$independent.variable.names
names(pred_m) <- make.names(names(pred_m), unique = TRUE)

if (!all(model_vars %in% names(pred_m))) {
  stop("Model expects variables not found in pred_m:\n",
       paste(setdiff(model_vars, names(pred_m)), collapse = ", "))
}

# Predict using only the required layers (and in correct order)
pred_for_map <- pred_m[[model_vars]]

msg("Predicting full map ...")
pred_map <- terra::predict(pred_for_map, rf_final, type = "response")

names(pred_map) <- paste0("pred_", names(y_m))

out_map <- file.path(out_dir, paste0("rf_pred_", names(y_m), "_from_s2.tif"))
writeRaster(pred_map, out_map, overwrite = TRUE)

# Residual map (Observed - Predicted)
resid_map <- y_m - pred_map
names(resid_map) <- "residual_obs_minus_pred"
writeRaster(resid_map, file.path(out_dir, "rf_residual_map.tif"), overwrite = TRUE)

# ============================================================
# 7) QUICKLOOK FIGURES
# ============================================================

png(file.path(out_dir, "quicklook_observed_target.png"), width = 1400, height = 1000, res = 150)
plot(y_m, main = paste("Observed LiDAR target (on Sentinel grid):", names(y_m)))
dev.off()

png(file.path(out_dir, "quicklook_predicted_target.png"), width = 1400, height = 1000, res = 150)
plot(pred_map, main = paste("RF Predicted:", names(pred_map)))
dev.off()

png(file.path(out_dir, "quicklook_residual.png"), width = 1400, height = 1000, res = 150)
plot(resid_map, main = "Residual (Observed - Predicted)")
dev.off()

png(file.path(out_dir, "variable_importance.png"), width = 1200, height = 900, res = 150)
par(mar=c(6,10,3,2))
top <- head(imp_df, 20)
barplot(rev(top$importance), names.arg = rev(top$variable), horiz = TRUE, las = 1,
        main = "RF Variable Importance (Top 20)")
dev.off()

msg("Saved prediction map:", out_map)
msg("Saved outputs in:", out_dir)

# Cleanup
terra::tmpFiles(remove = TRUE)
msg("DONE.")

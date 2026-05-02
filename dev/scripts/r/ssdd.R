# Structure Separation Distance Density (SSDD)
# =============================================
# R translation of SSDD.ipynb
#
# Computes building-level SSDD metrics for Wildland-Urban Interface (WUI) analysis.
#
# Two components:
#   SD (Structure Density)  -- kernel density + basal area fraction, blended via alpha_D
#   SS (Separation)         -- inverse-distance + orientation-weighted inverse-distance, blended via alpha_S
#   SSDD                    -- linear blend of SD and SS via beta
#
# Dependencies:
#   sf       >= 1.0-0  (st_minimum_rotated_rectangle requires sf >= 1.0)
#   ggplot2
#   pbapply            (progress-bar lapply; install with install.packages("pbapply"))
#
# Usage:
#   Edit the PARAMETERS section below, then source() or Rscript this file.


# =========================
# PARAMETERS
# =========================

INPUT_PATH   <- "/path/to/Buildings.shp"   # SHP / GPKG / GeoJSON
INPUT_LAYER  <- NULL                         # e.g. "buildings" for GPKG; NULL for SHP

OUT_DIR      <- "/path/to/output"
RUN_NAME     <- "my_run"

TARGET_EPSG  <- 3310   # NAD83 / California Albers (equal-area, meters)

# ---------- SD knobs ----------
r_D            <- 100.0    # SD buffer radius (m)
alpha_D        <- 0.5      # 1 = kernel density only, 0 = basal area only
kernel_type    <- "quartic"
weight_by_area <- FALSE    # TRUE => weight neighbor contributions by footprint area

# ---------- SS knobs ----------
r_S         <- 50.0    # SS neighbor search radius (m)
epsilon     <- 0.5     # distance floor to avoid division-by-zero (m)
sigma_theta <- 15.0    # orientation tolerance (deg); smaller => orientation matters more
alpha_S     <- 0.5     # 1 = distance only, 0 = orientation-weighted distance only

# ---------- Master blend ----------
beta <- 0.5   # 1 = SD only, 0 = SS only

# ---------- Normalization ----------
NORM_METHOD <- "robust"   # "minmax" or "robust"
P_LOW  <- 2
P_HIGH <- 98

# ---------- QA ----------
MAKE_QA_PLOTS <- TRUE


# =========================
# SETUP
# =========================

library(sf)
library(ggplot2)
library(pbapply)


# =========================
# UTILITY FUNCTIONS
# =========================

norm_vec <- function(x, method = "minmax", p_low = 2, p_high = 98) {
  # Normalize a numeric vector to [0, 1].
  x <- as.numeric(x)
  if (method == "minmax") {
    lo <- min(x, na.rm = TRUE)
    hi <- max(x, na.rm = TRUE)
    if (hi == lo) return(rep(0, length(x)))
    (x - lo) / (hi - lo)
  } else if (method == "robust") {
    lo <- as.numeric(quantile(x, p_low  / 100, na.rm = TRUE))
    hi <- as.numeric(quantile(x, p_high / 100, na.rm = TRUE))
    if (hi == lo) return(rep(0, length(x)))
    x_clip <- pmax(pmin(x, hi), lo)
    (x_clip - lo) / (hi - lo)
  } else {
    stop(paste("Unknown normalization method:", method))
  }
}

quartic_kernel <- function(u) {
  # K(u) = (1 - u^2)^2 for u in [0, 1], else 0.
  ifelse(u < 0 | u > 1, 0.0, (1.0 - u^2)^2)
}

kernel_value <- function(u, kernel = "quartic") {
  if (kernel == "quartic") return(quartic_kernel(u))
  stop(paste("Unknown kernel type:", kernel))
}

dominant_orientation_degrees <- function(poly_geom) {
  # Angle of the longest edge of the minimum rotated rectangle, in [0, 180).
  # Requires sf >= 1.0-0.
  mrr    <- st_minimum_rotated_rectangle(st_sfc(poly_geom, crs = NA_crs_))
  coords <- st_coordinates(mrr)[1:5, c("X", "Y")]   # 4 corners + closing point

  lengths <- numeric(4)
  dxs     <- numeric(4)
  dys     <- numeric(4)
  for (k in seq_len(4)) {
    dx        <- coords[k + 1, "X"] - coords[k, "X"]
    dy        <- coords[k + 1, "Y"] - coords[k, "Y"]
    lengths[k] <- sqrt(dx^2 + dy^2)
    dxs[k]    <- dx
    dys[k]    <- dy
  }

  best  <- which.max(lengths)
  angle <- (atan2(dys[best], dxs[best]) * 180 / pi) %% 180
  angle
}

angle_difference_deg <- function(a, b) {
  # Acute angular difference between two orientations in [0, 90].
  diff <- abs(a - b) %% 180
  diff <- pmin(diff, 180 - diff)
  pmin(diff, 90)
}

orientation_factor <- function(theta_deg, sigma_theta) {
  # g(theta) = exp(-(theta / sigma)^2).
  if (sigma_theta <= 0) return(rep(1, length(theta_deg)))
  exp(-((theta_deg / sigma_theta)^2))
}

qa_hist <- function(x, title, xlab = "") {
  df <- data.frame(x = x[!is.na(x)])
  print(
    ggplot(df, aes(x = x)) +
      geom_histogram(bins = 40) +
      labs(title = title, x = xlab, y = "Count")
  )
}


# =========================
# 1. LOAD, CLEAN, REPROJECT
# =========================

t0 <- proc.time()
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading buildings...\n")
if (is.null(INPUT_LAYER)) {
  bld <- st_read(INPUT_PATH, quiet = TRUE)
} else {
  bld <- st_read(INPUT_PATH, layer = INPUT_LAYER, quiet = TRUE)
}

# Keep only valid polygon geometries
bld <- bld[!st_is_empty(bld), ]
bld <- bld[st_geometry_type(bld) %in% c("POLYGON", "MULTIPOLYGON"), ]
bld <- st_buffer(bld, dist = 0)   # fix common geometry issues

cat(sprintf("  Input CRS  : %s\n", st_crs(bld)$input))
bld <- st_transform(bld, crs = TARGET_EPSG)
cat(sprintf("  Analysis CRS: EPSG:%d\n", TARGET_EPSG))

if (!"bld_id" %in% names(bld)) {
  bld$bld_id <- seq_len(nrow(bld)) - 1L
}
bld$bld_area <- as.numeric(st_area(bld))
bld$rep_pt   <- st_point_on_surface(st_geometry(bld))
cat(sprintf("  Loaded %s buildings\n", format(nrow(bld), big.mark = ",")))

# Geometry arrays (sfc) for efficient indexing
polys    <- st_geometry(bld)           # polygon sfc
pts_geom <- bld$rep_pt                 # point sfc

# sf data frames used for st_distance / st_intersects
polys_sf <- st_as_sf(data.frame(idx = seq_len(nrow(bld))), geometry = polys)
pts_sf   <- st_as_sf(data.frame(idx = seq_len(nrow(bld))), geometry = pts_geom)


# =========================
# 2. DOMINANT ORIENTATION
# =========================

cat("Computing building orientations...\n")
bld$phi_deg <- pbsapply(seq_len(nrow(bld)), function(i) {
  dominant_orientation_degrees(polys[[i]])
})


# =========================
# 3. PRE-COMPUTE SPATIAL CANDIDATE LISTS
# Using st_intersects() with buffered geometries to find candidates efficiently.
# sf uses an R-tree index internally.
# =========================

cat("Building spatial candidate lists...\n")

# SD candidate points (rep_pt within r_D of each rep_pt)
pts_buf_sd   <- st_buffer(pts_sf,   dist = r_D)
cand_pts_sd  <- st_intersects(pts_buf_sd, pts_sf,   sparse = TRUE)

# SD candidate polygons (polygon intersecting r_D buffer of each polygon)
polys_buf_sd  <- st_buffer(polys_sf, dist = r_D)
cand_polys_sd <- st_intersects(polys_buf_sd, polys_sf, sparse = TRUE)

# SS candidate polygons (polygon intersecting r_S buffer of each polygon)
polys_buf_ss  <- st_buffer(polys_sf, dist = r_S)
cand_polys_ss <- st_intersects(polys_buf_ss, polys_sf, sparse = TRUE)


# =========================
# 4. SD: KERNEL DENSITY (KD) + BASAL AREA FRACTION (BA)
# =========================

compute_KD <- function(i) {
  # KD_i = (1 / (pi * r_D^2)) * sum_j  w_j * K(dist(ci, cj) / r_D)
  idxs <- setdiff(cand_pts_sd[[i]], i)
  if (length(idxs) == 0L) return(0)

  dists <- as.numeric(st_distance(pts_sf[i, ], pts_sf[idxs, ]))
  keep  <- dists <= r_D
  dists <- dists[keep]
  j_idx <- idxs[keep]
  if (length(dists) == 0L) return(0)

  u <- dists / r_D
  k <- kernel_value(u, kernel = kernel_type)
  w <- if (weight_by_area) bld$bld_area[j_idx] else rep(1.0, length(j_idx))
  sum(w * k) / (pi * r_D^2)
}

compute_BA <- function(i) {
  # BA_i = sum_j area(P_j ∩ buffer(P_i, r_D)) / area(buffer(P_i, r_D))
  win      <- st_buffer(polys_sf[i, ], dist = r_D)
  win_area <- as.numeric(st_area(win))
  if (win_area <= 0) return(0)

  idxs <- cand_polys_sd[[i]]
  if (length(idxs) == 0L) return(0)

  # Intersect candidate polygons with the window; sum resulting areas
  candidates_sf <- polys_sf[idxs, ]
  suppressWarnings({
    clipped <- st_intersection(candidates_sf, win)
  })
  if (nrow(clipped) == 0L) return(0)

  inter_area <- sum(as.numeric(st_area(clipped)), na.rm = TRUE)
  inter_area / win_area
}

cat("Computing SD components...\n")
cat("  KD (kernel density)...\n")
bld$KD_raw <- unlist(pblapply(seq_len(nrow(bld)), compute_KD))

cat("  BA (basal area fraction)...\n")
bld$BA_raw <- unlist(pblapply(seq_len(nrow(bld)), compute_BA))

bld$KD <- norm_vec(bld$KD_raw, method = NORM_METHOD, p_low = P_LOW, p_high = P_HIGH)
bld$BA <- norm_vec(bld$BA_raw, method = NORM_METHOD, p_low = P_LOW, p_high = P_HIGH)
bld$SD <- alpha_D * bld$KD + (1 - alpha_D) * bld$BA
cat(sprintf("  SD computed — mean=%.3f, sd=%.3f\n", mean(bld$SD), sd(bld$SD)))


# =========================
# 5. SS: DISTANCE PROXY (DP) + ORIENTATION PROXY (OP)
# =========================

compute_SS_terms <- function(i) {
  # Returns named vector c(dp, op, m).
  # DP_raw_i = mean(1 / (d_ij + eps))
  # OP_raw_i = mean(g(theta_ij) / (d_ij + eps))
  # d_ij = wall-to-wall distance between polygons i and j.
  idxs <- setdiff(cand_polys_ss[[i]], i)
  if (length(idxs) == 0L) return(c(dp = 0, op = 0, m = 0L))

  dij  <- as.numeric(st_distance(polys_sf[i, ], polys_sf[idxs, ]))
  keep <- dij <= r_S
  dij  <- dij[keep]
  j_idx <- idxs[keep]
  if (length(dij) == 0L) return(c(dp = 0, op = 0, m = 0L))

  inv    <- 1 / (dij + epsilon)
  phi_j  <- bld$phi_deg[j_idx]
  theta  <- angle_difference_deg(bld$phi_deg[i], phi_j)
  orient <- orientation_factor(theta, sigma_theta)

  c(dp = mean(inv), op = mean(orient * inv), m = length(dij))
}

cat("Computing SS components...\n")
cat("  SS terms (distance + orientation)...\n")
ss_list <- pblapply(seq_len(nrow(bld)), compute_SS_terms)
ss_mat  <- do.call(rbind, ss_list)

bld$DP_raw       <- ss_mat[, "dp"]
bld$OP_raw       <- ss_mat[, "op"]
bld$SS_neighbors <- as.integer(ss_mat[, "m"])

bld$DP <- norm_vec(bld$DP_raw, method = NORM_METHOD, p_low = P_LOW, p_high = P_HIGH)
bld$OP <- norm_vec(bld$OP_raw, method = NORM_METHOD, p_low = P_LOW, p_high = P_HIGH)
bld$SS <- alpha_S * bld$DP + (1 - alpha_S) * bld$OP
cat(sprintf("  SS computed — mean=%.3f, sd=%.3f\n", mean(bld$SS), sd(bld$SS)))


# =========================
# 6. SSDD + PROVENANCE FIELDS
# =========================

bld$SSDD      <- beta * bld$SD + (1 - beta) * bld$SS
bld$SSDD_geom <- (pmax(bld$SD, 1e-9)^beta) * (pmax(bld$SS, 1e-9)^(1 - beta))
cat(sprintf("  SSDD computed — mean=%.3f, sd=%.3f\n", mean(bld$SSDD), sd(bld$SSDD)))

# Store run parameters on every feature for self-describing outputs
bld$SD_r_m   <- as.numeric(r_D)
bld$SS_r_m   <- as.numeric(r_S)
bld$alpha_D  <- as.numeric(alpha_D)
bld$alpha_S  <- as.numeric(alpha_S)
bld$beta     <- as.numeric(beta)
bld$sigma_th <- as.numeric(sigma_theta)
bld$eps_m    <- as.numeric(epsilon)
bld$kernel   <- kernel_type
bld$w_area   <- weight_by_area
bld$norm     <- NORM_METHOD
bld$p_low    <- as.integer(P_LOW)
bld$p_high   <- as.integer(P_HIGH)


# =========================
# 7. QA PLOTS
# =========================

if (MAKE_QA_PLOTS) {
  qa_hist(bld$KD_raw, "KD_raw distribution",              "KD_raw")
  qa_hist(bld$BA_raw, "BA_raw (basal fraction) distribution", "BA_raw")
  qa_hist(bld$SD,     "SD distribution",                  "SD")
  qa_hist(bld$DP_raw, "DP_raw distribution",              "DP_raw")
  qa_hist(bld$OP_raw, "OP_raw distribution",              "OP_raw")
  qa_hist(bld$SS,     "SS distribution",                  "SS")
  qa_hist(bld$SSDD,   "SSDD distribution",                "SSDD")
}


# =========================
# 8. SAVE OUTPUTS
# =========================

run_tag <- sprintf(
  "SDr%d_SSr%d_aD%.2f_aS%.2f_b%.2f_sig%.0f_eps%.1f_norm%s",
  as.integer(r_D), as.integer(r_S),
  alpha_D, alpha_S, beta,
  sigma_theta, epsilon,
  NORM_METHOD
)

out_gpkg <- file.path(OUT_DIR, sprintf("%s_SSDD_%s.gpkg",            RUN_NAME, run_tag))
out_csv  <- file.path(OUT_DIR, sprintf("%s_SSDD_%s.csv",             RUN_NAME, run_tag))
out_txt  <- file.path(OUT_DIR, sprintf("%s_SSDD_%s_RunSummary.txt",  RUN_NAME, run_tag))

# Drop helper geometry column before saving
out_sf        <- bld
out_sf$rep_pt <- NULL

# GeoPackage (preserves geometry; preferred over Shapefile)
st_write(out_sf, out_gpkg, layer = "buildings_ssdd", driver = "GPKG",
         delete_dsn = TRUE, quiet = TRUE)

# CSV (attribute table, no geometry)
csv_cols <- c(
  "bld_id",
  "KD_raw", "BA_raw", "KD", "BA", "SD",
  "DP_raw", "OP_raw", "DP", "OP", "SS", "SS_neighbors",
  "SSDD", "SSDD_geom",
  "SD_r_m", "SS_r_m",
  "alpha_D", "alpha_S", "beta",
  "sigma_th", "eps_m",
  "kernel", "w_area",
  "norm", "p_low", "p_high"
)
csv_cols_existing <- intersect(csv_cols, names(out_sf))
out_tbl <- st_drop_geometry(out_sf)[, csv_cols_existing]
write.csv(out_tbl, out_csv, row.names = FALSE)

# Run summary text
elapsed_sec <- as.numeric((proc.time() - t0)["elapsed"])
writeLines(
  c(
    "SSDD run summary",
    "----------------",
    sprintf("RUN_NAME: %s", RUN_NAME),
    sprintf("INPUT: %s (layer=%s)", INPUT_PATH,
            ifelse(is.null(INPUT_LAYER), "NULL", INPUT_LAYER)),
    sprintf("N buildings: %s", format(nrow(bld), big.mark = ",")),
    sprintf("CRS (analysis): EPSG:%d", TARGET_EPSG),
    "",
    "Parameters:",
    sprintf("  SD radius r_D (m)    : %s", r_D),
    sprintf("  SS radius r_S (m)    : %s", r_S),
    sprintf("  alpha_D              : %s", alpha_D),
    sprintf("  alpha_S              : %s", alpha_S),
    sprintf("  beta                 : %s", beta),
    sprintf("  sigma_theta (deg)    : %s", sigma_theta),
    sprintf("  epsilon (m)          : %s", epsilon),
    sprintf("  kernel               : %s", kernel_type),
    sprintf("  weight_by_area       : %s", weight_by_area),
    sprintf("  normalization        : %s (P_LOW=%d, P_HIGH=%d)",
            NORM_METHOD, P_LOW, P_HIGH),
    "",
    "Outputs:",
    sprintf("  GPKG: %s (layer=buildings_ssdd)", out_gpkg),
    sprintf("  CSV : %s", out_csv),
    "",
    sprintf("Elapsed seconds: %.2f", elapsed_sec)
  ),
  out_txt
)

cat(sprintf("\nWrote: %s\n", out_gpkg))
cat(sprintf("Wrote: %s\n", out_csv))
cat(sprintf("Wrote: %s\n", out_txt))
cat(sprintf("\nDone in %.1fs\n", elapsed_sec))

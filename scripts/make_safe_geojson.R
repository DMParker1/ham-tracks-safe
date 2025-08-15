(library(sf)

# ---- user knobs ----
home_lat <- 33.6; home_lon <- -117.8   # EDIT: your QTH (approx is fine)
geofence_m <- 5000                         # drop points within this radius
min_age_days <- 14                          # only include tracks older than this
add_jitter <- TRUE                         # randomize remaining points slightly
jitter_min_m <- 800; jitter_max_m <- 800

in_dir  <- "data/raw"
out_geo <- "data/processed/tracks_safe.geojson"

# ---- helpers ----
hav <- function(lat1, lon1, lat2, lon2){
  R <- 6371000
  toRad <- pi/180
  dlat <- (lat2-lat1)*toRad; dlon <- (lon2-lon1)*toRad
  a <- sin(dlat/2)^2 + cos(lat1*toRad)*cos(lat2*toRad)*sin(dlon/2)^2
  2*R*asin(sqrt(a))
}

dest <- function(lat, lon, brg_deg, dist_m){
  R <- 6371000; brg <- brg_deg*pi/180
  lat1 <- lat*pi/180; lon1 <- lon*pi/180; dR <- dist_m/R
  lat2 <- asin(sin(lat1)*cos(dR) + cos(lat1)*sin(dR)*cos(brg))
  lon2 <- lon1 + atan2(sin(brg)*sin(dR)*cos(lat1), cos(dR)-sin(lat1)*sin(lat2))
  c(lat2*180/pi, ((lon2*180/pi + 540) %% 360) - 180)
}

read_any_gpx_points <- function(path){
  for (layer in c("track_points","tracks","waypoints")){
    try({
      g <- suppressWarnings(st_read(path, layer = layer, quiet = TRUE))
      if (nrow(g) > 0) return(g)
    }, silent = TRUE)
  }
  st_as_sf(data.frame(), coords = c("lon","lat"), crs = 4326)
}

# ---- load & combine ----
files <- list.files(in_dir, pattern = "\\.gpx$", full.names = TRUE)
if (length(files) == 0) stop("No GPX files found in data/raw")

pts_list <- lapply(files, read_any_gpx_points)
g <- do.call(rbind, pts_list)
if (nrow(g) == 0) stop("No points found in GPX layers")

# Try to normalize time column if present
time_col <- intersect(names(g), c("time","Name","name","desc"))
if ("time" %in% names(g)) {
  g$when <- as.POSIXct(g$time, tz = "UTC")
} else {
  g$when <- NA
}

# ---- geofence removal ----
coords <- st_coordinates(st_transform(g, 4326))
d <- hav(coords[,2], coords[,1], home_lat, home_lon)
keep <- d > geofence_m
g <- g[keep, , drop = FALSE]

# ---- delay filter ----
if (!all(is.na(g$when))) {
  cutoff <- Sys.time() - min_age_days*24*3600
  g <- g[ which(g$when < cutoff), , drop = FALSE]
}

# ---- jitter (optional) ----
if (add_jitter && nrow(g) > 0) {
  coords <- st_coordinates(st_transform(g, 4326))
  dist <- runif(nrow(g), jitter_min_m, jitter_max_m)
  brg  <- runif(nrow(g), 0, 360)
  jittered <- t(mapply(dest, coords[,2], coords[,1], brg, dist))
  g <- st_set_geometry(g, st_sfc(lapply(seq_len(nrow(g)), function(i)
    st_point(c(jittered[i,2], jittered[i,1]))), crs = 4326))
}

# ---- write GeoJSON ----
dir.create(dirname(out_geo), showWarnings = FALSE, recursive = TRUE)
st_write(g, out_geo, append = FALSE, quiet = TRUE)
cat("Wrote:", out_geo, "with", nrow(g), "points\n")

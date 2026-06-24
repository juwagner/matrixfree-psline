# ------------------------------------------------------------------------------
# Prepare LUCAS dataset
# ------------------------------------------------------------------------------

load("data/LUCAS.SOIL_corr.Rdata")

# ------------------------------------------------------------------------------
# Map spectrals to spectral bands
bands <- list(
  B2 = c(493 - 32, 493 + 32),      # blue
  B3 = c(560 - 17.5, 560 + 17.5),  # green
  B4 = c(665 - 15,   665 + 15),    # red
  B5 = c(705 - 7.5,  705 + 7.5),   # red edge 1
  B6 = c(740 - 7.5,  740 + 7.5),   # red edge 2
  B7 = c(783 - 10,   783 + 10),    # red edge 3
  B8 = c(842 - 57.5, 842 + 57.5),  # NIR
  B8a = c(865 - 10,   865 + 10),   # narrow NIR
  B9 = c(945 - 10,   945 + 10),    # water vapor
  B11 = c(1610 - 45,  1610 + 45),  # SWIR1
  B12 = c(2190 - 90,  2190 + 90)   # SWIR2
)

spc <- LUCAS.SOIL$spc
wl <- as.numeric(colnames(spc))

aggregate_band <- function(spc, wl, lower, upper) {
  idx <- wl >= lower & wl <= upper
  rowMeans(spc[, idx])
}

df <- as.data.frame(
  lapply(bands, function(b) {
    aggregate_band(spc, wl, b[1], b[2])
  })
)

# ------------------------------------------------------------------------------
# Assemble and clean data (non-unique ID)
LUCAS_agg <- data.frame(
  ID = LUCAS.SOIL$sample.ID,
  df,
  date = LUCAS.SOIL$date,
  lat = LUCAS.SOIL$GPS_LAT,
  long = LUCAS.SOIL$GPS_LONG,
  OC = LUCAS.SOIL$OC
)

LUCAS_agg <- LUCAS_agg[!duplicated(LUCAS_agg$ID), ]

# ------------------------------------------------------------------------------
# Store data
save(LUCAS_agg, file = "data/LUCAS_agg.Rdata")

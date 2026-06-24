# ------------------------------------------------------------------------------
# Load LUCAS dataset
# ------------------------------------------------------------------------------

library(sp)
library(BBmisc)

# ------------------------------------------------------------------------------
# raw data

file_path <- "data/LUCAS_agg.Rdata"
if (file.exists(file_path)) {
  load("data/LUCAS_agg.Rdata")
} else {
  source("src/utils/prepare_lucas_data.R")
}

X <- as.matrix(LUCAS_agg[c("B2", "B3", "B4")])  # selected spectral bands
U <- cbind(LUCAS_agg$lat, LUCAS_agg$long)       # geo coordinates
y <- LUCAS_agg$OC                               # organic carbonate
n <- length(y)

# ------------------------------------------------------------------------------
# data preprocessing
indicator <- which(U[,1]>85 | U[,2]>85)
X <- normalize(X[-indicator,], margin=2)
U <- normalize(U[-indicator,], margin=2)
y <- normalize(y[-indicator], method="range")

# ------------------------------------------------------------------------------
# cleanup
rm(LUCAS_agg, indicator)
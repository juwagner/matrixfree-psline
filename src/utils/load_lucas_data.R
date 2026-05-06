# ------------------------------------------------------------------------------
# Load LUCAS dataset
# ------------------------------------------------------------------------------

library(sp)
library(BBmisc)

# ------------------------------------------------------------------------------
# raw data
load("data/lucas_data.Rdata")

X <- as.matrix(LUCAS.SOIL$spc[2:4])                   # selected spectral bands
U <- cbind(LUCAS.SOIL$GPS_LAT, LUCAS.SOIL$GPS_LONG)   # geo coordinates
y <- LUCAS.SOIL$OC                                    # organic carbonate
n <- length(y)

# ------------------------------------------------------------------------------
# data preprocessing
indicator <- which(U[,1]>85 | U[,2]>85)
X <- normalize(X[-indicator,], margin=2)
U <- normalize(U[-indicator,], margin=2)
y <- normalize(y[-indicator], method="range")

# ------------------------------------------------------------------------------
# cleanup
rm(LUCAS.SOIL, indicator)
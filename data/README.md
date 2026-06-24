

The data used in the repo are from the Land Use and Coverage Area frame Survey (LUCAS) 2009/2012.
The data are freely available and can be downloaded from the official European Soil Data Centre (ESDAC) after prior registration.
(See https://esdac.jrc.ec.europa.eu/projects/lucas).

Load the downloaded CSV file into R and convert it into a DataFrame.
Store the DataFrame into "data/LUCAS_raw.Rdata".

The workflow performs the following steps:
- Loads the raw LUCAS spectral data
- Extracts the spectral measurements (spc)
- Transforms the high-resolution spectra (~4200 wavelengths) into aggregated spectral bands (Sentinel‑2 bands)
- Creates a reduced dataset containing only the relevant variables
- Stores the processed dataset into "data/LUCAS_agg.Rdata"

The aggregate dataset is used for further computation.
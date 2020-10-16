#require(remotes)
#library(tidyverse) # install.packages("tidyverse")
#library(janitor) # install.packages("janitor")
#library(scales)
#library(mgcv)  # remotes::install_github("ps")
#library(EpiNow2) # remotes::install_github("epiforecasts/EpiNow2")
# will ask to install RTools from Rtools35.exe first
#library(frs)     # remotes::install_github("ellisp/frs-r-package/pkg", force=T)
#library(patchwork)
#library(glue) # install.packages("glue") remove.packages("glue")
#library(surveillance) # for backprojNP()
#library(ggplot2)
#library(future)
# Load the package required to read JSON files from web.
library(httr) 
library(jsonlite) 
#library(rlist) 
#library(dplyr) 
#library(data.table)


#### input arg
defaultArgs <- list (
  outFile = NULL
)

args <- R.utils::commandArgs(trailingOnly = TRUE,
                             asValues = TRUE ,
                             defaults = defaultArgs)


# data read.
# jsonURL <- "https://opendata.arcgis.com/datasets/b913e9591eae4912b33dc5b4e88646c5_10.geojson"
jsonURL <- "https://opendata.arcgis.com/datasets/89d7a90aafa24519847c89b249be96ca_13.geojson"

jsonResponse<-GET(jsonURL)
http_type(jsonResponse)
jsonResponseText <- content(jsonResponse, as = "text") #JSON response structured into raw data
dataJSON = fromJSON(jsonResponseText)
data_WI = dataJSON$features
data_WI = do.call("rbind", data_WI)
# Convert JSON file to a data frame.
data_WI_check <- as.data.frame(data_WI)
###############################################

saveRDS(data_WI_check, args$outFile)
q()

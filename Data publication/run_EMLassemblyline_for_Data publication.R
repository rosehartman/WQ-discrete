# This script executes an EMLassemblyline workflow.

# Initialize workspace --------------------------------------------------------

# Update EMLassemblyline and load

#remotes::install_github("EDIorg/EMLassemblyline")
library(EMLassemblyline)
library(dplyr)
library(readr)

# Define paths for your metadata templates, data, and EML

root<-"Data publication"

path_templates <- file.path(root, "metadata_templates")
path_data <- file.path(root, "data_objects")
path_eml <- file.path(root, "eml")


# Add data ----------------------------------------------------------------

data<-discretewq::wq()%>%
  mutate(Date=as.character(Date, format="%Y-%m-%d"),
         Datetime=as.character(Datetime),
         Notes=stringr::str_replace_all(Notes, '"', "'"),# Replace full quotes with single quotes to avoid data read errors
         Notes=stringr::str_replace_all(Notes, stringr::fixed('\n'), " "), # Replace line breaks with space to avoid data read errors
  Notes=stringr::str_replace_all(Notes, stringr::fixed('\r'), " "))%>% # Replace line breaks with space to avoid data read errors
  select(Source, Station=StationID, Latitude, Longitude, Field_coords, Date, Datetime, Depth, 
         Sample_depth_surface, Sample_depth_bottom, Tide, Temperature, Temperature_bottom, 
         Conductivity, Salinity, Secchi, Microcystis, Chlorophyll, Notes)

write_csv(data, file.path(path_data, "Delta_Integrated_WQ.csv"), )
# Create metadata templates ---------------------------------------------------

# Below is a list of boiler plate function calls for creating metadata templates.
# They are meant to be a reminder and save you a little time. Remove the 
# functions and arguments you don't need AND ... don't forget to read the docs! 
# E.g. ?template_core_metadata

# Create core templates (required for all data packages)

EMLassemblyline::template_core_metadata(
  path = path_templates,
  license = "CCBY",
  file.type = ".docx")

# Create provenance template

EMLassemblyline::template_provenance(
  path=path_templates,
)

# Create table attributes template (required when data tables are present)

EMLassemblyline::template_table_attributes(
  path = path_templates,
  data.path = path_data,
  data.table = c("Delta_Integrated_WQ.csv", "Delta_Integrated_WQ_metadata.csv"))

# Create categorical variables template (required when attributes templates
# contains variables with a "categorical" class)


#### NEED TO DO THIS AND BELOW ###
EMLassemblyline::template_categorical_variables(
  path = path_templates, 
  data.path = path_data)

# Create geographic coverage (required when more than one geographic location
# is to be reported in the metadata).

EMLassemblyline::template_geographic_coverage(
  path = path_templates, 
  data.path = path_data, 
  data.table = "Delta_Integrated_WQ.csv", 
  lat.col = "Latitude",
  lon.col = "Longitude",
  site.col = "Station")

# Make EML from metadata templates --------------------------------------------

# Once all your metadata templates are complete call this function to create 
# the EML.

EMLassemblyline::make_eml(
  path = path_templates,
  data.path = path_data,
  eml.path = path_eml, 
  dataset.title = "Six decades (1959-2020) of water quality in the upper San Francisco Estuary: an integrated database of 11 discrete monitoring surveys in the Sacramento San Joaquin Delta, Suisun Bay, and Suisun Marsh", 
  temporal.coverage = c("1959-06-13", "2020-07-10"), 
  maintenance.description = "ongoing", 
  data.table = c("Delta_Integrated_WQ.csv", "Delta_Integrated_WQ_metadata.csv"), 
  data.table.description = c("Integrated water quality database", "Information on each survey included in the integrated database"),
  data.table.quote.character=c('"','"'),
  user.id = "sbashevkin",
  user.domain = "EDI", 
  package.id = "edi.731.1")

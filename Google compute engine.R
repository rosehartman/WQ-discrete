require(googleComputeEngineR)
require(googleCloudStorageR)

#gcs_create_bucket(name="discretewq",
#                  projectId="models-277122",
#                  location="US-WEST1",
#                  storageClass = "STANDARD")

# Upload necessary files with something like:
# gcs_upload("Discrete Temp Data.Rds")
# NOTE: Seems to work better with objects that can be "loaded" (i.e. those saved with `save` not `saveRDS`)

vm <- gce_vm(template = "rstudio", zone="us-west1-a",
             name = "rstudio-server",
             username = "", password = "",  predefined_type = "e2-highmem-16",
             dynamic_image = "gcr.io/gcer-public/persistent-rstudio",
             disk_size_gb=100)
gce_set_metadata(list(GCS_SESSION_BUCKET = "discretewq"), vm)

# To save all files from VM, run in cloud: googleCloudStorageR::gcs_save_all(bucket="discretewq")
# Load files from local machine with something like:
# gcs_load("model.Rds")

# Also can save outputs with something like:
# gcs_upload(model, object_function = function(input, output) save(input, file=output), name="model.Rds")


# To add a project
# 1) Create Rstudio project in the VM
# 2) Create file within the project folder named "_gcssave.yaml" with the following 2 lines:

# bucket: discretewq
# loaddir:
  
# 3) Create .Rprofile file within project folder with the following 15 lines:

#   .First <- function(){
#     cat("\n# Welcome Sam! Today is ", date(), "\n")
#     
#     ## will look for download if GCS_SESSION_BUCKET env arg set
#     googleCloudStorageR::gcs_first()
#   }
# 
# 
# .Last <- function(){
#   # will only upload if a _gcssave.yaml in directory with bucketname
#   googleCloudStorageR::gcs_last()
#   message("\nGoodbye Sam at ", date(), "\n")
# }
# 
# message("n*** Successfully loaded .Rprofile ***n")

# 4) exit project and stop VM
# 5) At next startup of VM or new VM, just create a project with same name and files will be added. 


# Download results --------------------------------------------------------

gcs_get_object("modellc4_predictions.Rds", saveToDisk="Temperature smoothing model/modellc4_predictions.Rds", overwrite = T)
#gcs_get_object("modellc4.Rds", saveToDisk="Temperature smoothing model/modellc4.Rds")

#gcs_get_object("CC_gam1_predictions.Rds", saveToDisk="Temperature smoothing model/CC_gam1_predictions.Rds")
#gcs_get_object("CC_gam1.Rds", saveToDisk="Temperature smoothing model/CC_gam1.Rds")

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
             username = "", password = "",  predefined_type = "e2-highmem-8",
             dynamic_image = "gcr.io/gcer-public/persistent-rstudio",
             disk_size_gb=100)
gce_set_metadata(list(GCS_SESSION_BUCKET = "discretewq"), vm)

# Create project named "WQ-discrete-cloud"
# Saved objects should be automatically backed up to the bucket

# Load files from local machine with something like:
# gcs_load("model.Rds")

# Also can save outputs with something like:
# gcs_upload(model, object_function = function(input, output) save(input, file=output), name="model.Rds")
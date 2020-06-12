require(googleComputeEngineR)

vm <- gce_vm(template = "rstudio", zone="us-west1-a",
             name = "rstudio-server",
             username = "", password = "",  predefined_type = "e2-highmem-16",
             disk_size_gb=100)

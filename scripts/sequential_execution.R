# Be sure that working dir is set to the project director setwd("path/to/your/project")

# If any modifications are needed to the pipeline after running the shiny app, must remember to switch working dir or you will get issues with file saving and data imports.

setwd("~/RNAseq/endocrine-ER-stress")

# Sequential execution
source("scripts/00_Data_Preprocessing.R")
source("scripts/01_Master_RNAseq_Pipeline.R")
source("scripts/02_Comparative_Analysis_Pipeline.R")
source("scripts/03_Granular_Analysis_DotPlots.R")
source("scripts/04_Custom_Unique_genesets_DotPlots.R")
source("scripts/05_GSEA_Analysis_Unique_Genes.R")

setwd("RNAseq_Shiny_App/")  # switch working dir to the shiny app dir

source("scripts/06_prepare_data_for_shiny.R")

source("app.R") # execute app.R to use the app

# source("deploy.R") # this deploys the app.R to shinyapps.io

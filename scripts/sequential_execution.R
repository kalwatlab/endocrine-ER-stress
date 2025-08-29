# Be sure that working dir is set to the project directory e.g. setwd("path/to/your/project")

setwd("~/RNAseq/endocrine-ER-stress")

# Sequential execution
source("scripts/00_Data_Preprocessing.R")
source("scripts/01_Master_RNAseq_Pipeline.R")
source("scripts/02_Comparative_Analysis_Pipeline.R")
source("scripts/03_Granular_Analysis_DotPlots.R")
source("scripts/04_Custom_Unique_genesets_DotPlots.R")
source("scripts/05_GSEA_Analysis_Unique_Genes.R")
source("scripts/06_prepare_data_for_shiny.r")

setwd("RNAseq_Shiny_App/")
source("app.R") # source app.R to load everything - be sure to switch working dir to the shiny app directory before running the app
shiny::runApp("app.R") # run the app

# source("deploy.R") # this deploys the app.R to shinyapps.io

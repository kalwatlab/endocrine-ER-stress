# Deployment script for RNAseq Shiny Browser

library(rsconnect)

# Check files exist
if (!file.exists("app.R")) {
  stop("app.R not found!")
}

if (!file.exists("merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt")) {
  stop("Data file not found!")
}

# Get file size
data_size <- file.size("merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt") / 1024^2
cat("Data file size:", round(data_size, 2), "MB\n")

if (data_size > 100) {
  warning("Data file may be too large for free tier!")
}

# Deploy
cat("Deploying app...\n")
deployApp(
  appName = "endocrine-ER-stress",
  appTitle = "Endocrine ER Stress RNAseq Browser",
  forceUpdate = TRUE
)

cat("Deployment complete!\n")
cat("App URL: https://diabetes-detectives.shinyapps.io/endocrine-ER-stress/\n")

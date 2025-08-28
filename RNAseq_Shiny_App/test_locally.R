# Test the Shiny app locally

# For testing with sample data (faster):
if (file.exists("merged edgeR and TPMs/merged_edgeR_TPMs_all_lines_SAMPLE.txt")) {
  cat("Testing with sample data...\n")
  file.rename("merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt",
              "merged edgeR and TPMs/merged_edgeR_TPMs_all_lines_FULL.txt")
  file.rename("merged edgeR and TPMs/merged_edgeR_TPMs_all_lines_SAMPLE.txt",
              "merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt")
}

# Run the app
shiny::runApp()

# To switch back to full data:
# file.rename('merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt',
#             'merged edgeR and TPMs/merged_edgeR_TPMs_all_lines_SAMPLE.txt')
# file.rename('merged edgeR and TPMs/merged_edgeR_TPMs_all_lines_FULL.txt',
#             'merged edgeR and TPMs/merged_edgeR_TPMs_all_lines.txt')

# RNAseq Shiny Browser

## Quick Start

1. **Test locally with sample data:**
   ```r
   source('test_locally.R')
   ```

2. **Test with full data:**
   ```r
   shiny::runApp()
   ```

3. **Deploy to shinyapps.io:**
   ```r
   source('deploy.R')
   ```

## Files

- `app.R` - Main Shiny application
- `merged edgeR and TPMs/` - Data directory
  - `merged_edgeR_TPMs_all_lines.txt` - Full dataset
  - `merged_edgeR_TPMs_all_lines_SAMPLE.txt` - Sample dataset for testing
- `deploy.R` - Deployment script
- `test_locally.R` - Local testing script

## Data Size

- Full dataset: 9.75 MB
- Sample dataset: 0.23 MB

## Notes

- Free shinyapps.io tier limited to 1GB RAM and 25 active hours/month
- Use sample dataset for development and testing
- Deploy with full dataset for production

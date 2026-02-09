# SQPAD Analysis

This folder contains scripts to reproduce results from Lele and Solymos.

The R Markdown file was used for producing the Online Supplement.

The `00-functions-sqpad.R` file contains functions for fitting and 
summarizing SQPAD models to data. Note that SQPAD is available as part
of the detect R package (<https://github.com/psolymos/detect>).

The bSims R package was used for conducting _in-silico_ simulations
(<https://github.com/psolymos/bSims>).

Files:

```text
├── README.md                          # readme
├── 00-cloud-vm.md                     # run scripts on server
├── 00-functions-sqpad.R               # core SQPAD functions
├── 01-functions-helpers.R             # helper functions
├── 02-bsims-simulations.R             # script for simulations
├── 03-bsims-estimation.R              # script for estimation
├── 04-field-data.R                    # field data analysis
├── 05-summaries.R                     # produce summaries
├── 06-visualization.R                 # produce plots
├── estimates_bsims_mc.parquet         # saved simulation results
├── estimates_field_data_mc.parquet    # saved field data results
├── sqpad-appendix.pdf                 # PDF appendix with math
├── sqpad-supplement.pdf               # PDF appendix with results
└── sqpad-supplement.Rmd               # R markdown file for PDF
```

## Combining PDF files

Combine the appendix and supplement into a single PDF file:

```R
qpdf::pdf_combine(
    input = c(
        "analysis/sqpad-appendix.pdf", 
        "analysis/sqpad-supplement.pdf"),
    output = "paper/sqpad-online-supplement.pdf")
```

# SQPAD Analysis

This folder contains scripts to reproduce results from Lele and Solymos.

The R Markdown file was used for producing the Online Supplement.

The `00-functions-sqpad.R` file contains functions for fitting and 
summarizing SQPAD models to data. Note that SQPAD is available as part
of the detect R package (<https://github.com/psolymos/detect>).

The bSims R package was used for conducting in-silico simulations
(<https://github.com/psolymos/bSims>).

## Combining PDF files

Combine the appendix and supplement into a single PDF file:

```R
qpdf::pdf_combine(
    input = c(
        "analysis/sqpad-appendix.pdf", 
        "analysis/sqpad-supplement.pdf"),
    output = "paper/sqpad-online-supplement.pdf")
```

## Analysis for IST publication

![](9_misc/logos/InSilicoTreatment_Logo_Colour.png){width=25%}

This repository contains the analysis on public and in-house datasets necessary to reproduce the results in the IST (In Silico Treatment) publication.
IST is a computational method to assess the translatability between animal models and human disease data at the pathway and gene level. 
The manuscript and the package vignette will guide the user on IST's capabilities and user manual.

### Repo structure

The organisation is:

* `0_rawdata`: internal/external gene expression data (quantified, differential analysis)
* `1_preprocessing`: data fetching, wrangling and preparation for IST runs
* `2_IPF`: IPF animal models and treatments
* `3_NASH`: NASH animal models and treatments
* `9_misc`: miscellaneous supporting material (figures, logos, etc)
* `9_temp`: folder for testing code, not included in publication

Each report comes with a companion folder with the exported results.
If the report name is `3_report.Rmd`, the output folder will be `3_report_output/`, so it is always caught by `.gitignore`.

Important file locations are defined in the `config.yml` file and fecthed using the `config` R package.

### Requirements

The `IST` R package, also published with the manuscript.

Package dependencies are described in the `renv.lock` file using the `renv` R package via `renv::snapshot()`.

### Reproduce analysis

The `Makefile` contains automated instructions to rerun the whole repo. 
You can reproduce regenerate the whole analysis by running:

```
make all
```

To run specific parts from the tree structure, there are more specific commands:

```
make preprocessing
make ipf
make nash
```


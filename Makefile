SHELL=/bin/bash
RSCRIPT=Rscript

# you may need some lines nashdescriptive
# to add other environmental variables, 
# use a controlled R package environment, etc

# tasks: 3 steps
all: preprocessing ipf nash

# rmd files
%.html: %.Rmd; $(RSCRIPT) -e "require(rmarkdown); render('$<');"

# Step 1 - data prep
preprocessing: genemaps ipfhuman ipfsig nashhuman nashsig genesets

genemaps: 1_preprocessing/0_genemapping.Rmd 1_preprocessing/0_genemapping.html
ipfhuman: 1_preprocessing/1_IPF_human.Rmd 1_preprocessing/1_IPF_human.html
ipfsig: 1_preprocessing/2_IPF_signatures.Rmd 1_preprocessing/2_IPF_signatures.html
nashhuman: 1_preprocessing/3_NASH_human.Rmd 1_preprocessing/3_NASH_human.html
nashsig: 1_preprocessing/4_NASH_signatures.Rmd 1_preprocessing/4_NASH_signatures.html
genesets: 1_preprocessing/5_genesets.Rmd 1_preprocessing/5_genesets.html

# Step 2 - IPF analysis
ipf: ipfdescriptive ipfgsea ipfist

ipfdescriptive: 2_IPF/1_descriptive.Rmd 2_IPF/1_descriptive.html
ipfgsea: 2_IPF/2_GSEA.Rmd 2_IPF/2_GSEA.html
ipfist: 2_IPF/3_IST_disease.Rmd 2_IPF/3_IST_disease.html

# Step 3 - NASH analysis
nash: nashdescriptive nashgsea nashist

nashdescriptive: 3_NASH/1_descriptive.Rmd 3_NASH/1_descriptive.html
nashgsea: 3_NASH/2_GSEA.Rmd 3_NASH/2_GSEA.html
nashist: 3_NASH/3_IST_disease.Rmd 3_NASH/3_IST_disease.html

# Exploring pig single cell immune tissue Shiny App
This Shiny app allows interactive exploration of scRNA-seq data across four porcine immune tissues (Bone marrow, spleen, thymus, lymph node).

# This interactive app features: 
Home Page: Includes ready-to-use and preprocessed .cloupe files, and cell-type-specific differential gene expression analysis.
Gene Expression Page: Visualize expression of selected or typed-in genes across tissues using cell annotations, UMAP, and violin plots.
Reference Mapping Page: Compare immune cell types from the reference (Bone Marrow) to predicted cell types in other tissues.

** Large data files (.rds, .cloupe) are excluded from this repository
The deployable feature of the app expects data files(.rds) n the data/ directory- if you want to run this app locally, please download and place the .rds files under data/ directory. Link to exsisting .rds file: (AgBioCommons?)

# Run ShinyPIGGI locally:
To run this app locally, follow the steps below:

## 1. R version checks ##
 Check you R version in R Studio using command: ```r R.version.string ```
 This app requires R version >= 4.3.0

 ## 2. Git Clone the repository ##
 git clone https://github.com/kapoormuskan/Pig_Immune_Tissue_ShinyApp.git
 cd Pig_Immune_Tissue_ShinyApp

 ## 3. Install required packages for app by running the following commands: ""
 ```r
setwd("path/to/Pig_Immune_Tissue_ShinyApp")  # replace this with your path
packages <- readLines("r-requirements.txt")
install.packages(packages, repos = "https://cloud.r-project.org")
```

## 4. Launch the app inside R studio itself ##
```r
shiny::runApp('.', launch.browser = TRUE)
```

## ShinyPIGGI citation ##


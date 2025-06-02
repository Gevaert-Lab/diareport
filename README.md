# 🎯  diareport R package  version of DIA-Report

Welcome to **DIA-Report**! Your ultimate buddy for DE reports on DIA proteomics data. 🧑‍🔬🔬

With **DIA-Report**, you can generate advanced differential expression analysis reports based on **DIA-NN results** with just one magical command. ✨

The statistical analysis is powered by **MSQrob2** and **Qfeatures** packages, allowing you to filter your data with various parameters and perform differential expression on proteomics data. The high number of parameters allows you to customize the analysis to your needs (normalization methods, missing value filtering, aggregation methods) and avoid relying solely on black-box statistical methods. 🎛️🧪


## 📋 Requirements

-   [R](https://www.r-project.org/) (v.4.5.0). It should work with older versions but has not been fully tested.
-   [Bioconductor](https://www.bioconductor.org/install/): `install.packages("BiocManager")`
-   Install phantomjs in your R installation: `webshot::install_phantomjs()`
-   [Quarto](https://quarto.org/docs/download/) for creating HTML reports

## 🛠️ Installation

1.  Install R and Bioconductor:

    ``` r
    install.packages("BiocManager")
    ```

2.  Install phantomjs:

    ``` r
    webshot::install_phantomjs()
    ```

2.  Install devtools:

    ``` r
    install.packages("devtools")
      ```  

3.  Install Quarto:

    Follow the instructions on the [Quarto website](https://quarto.org/docs/download/).
    
    
## How to install the R package

1. Clone the Github repository 

2. Use `devtools::install()` to install the diareport package in your system.


## 📂 Quarto Templates Availabale 

The repository contains two Quarto templates:

-   `Template_DIA-NN_dev.qmd`: Template for protein-level analysis report
-   `Template_DIA-NN_peptide_dev.qmd`: Template for peptide-level analysis report

## How to run an analsysis 

To run an analysis you  can use :
`render_dia_report(params, template_file, report_target_folder, output_filename)`

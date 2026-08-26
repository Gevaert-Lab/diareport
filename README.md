
[![Build and Push Docker Image](https://github.com/Gevaert-Lab/diareport/actions/workflows/release.yml/badge.svg?branch=main)](https://github.com/Gevaert-Lab/diareport/actions/workflows/release.yml)
# 🎯 diareport R Package

Welcome to **DiaReport** — your ultimate companion for differential expression (DE) reports on DIA proteomics data. 🧑‍🔬🔬

With **DiaReport**, you can generate advanced differential expression analysis reports based on **DIA-NN results** with just one magical command. ✨

The statistical analysis is powered by the **MSqrob2** and **QFeatures** packages, enabling you to filter your data with various parameters and perform differential expression on proteomics data. The broad selection of parameters allows you to customize the analysis to your needs (normalization methods, missing value filtering, aggregation methods) and avoid reliance on black-box statistical methods. 🎛️🧪

------------------------------------------------------------------------


## 🛠️ Requirements Installation

1.  **Install R  version >= 4.4.0**   

2.  **Install Bioconductor dependencies**   
    ``` r 
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
 
        BiocManager::install(c("QFeatures", "SummarizedExperiment", "MsCoreUtils", "msqrob2"))
    ```
3.  **Install R Package:** `install.packages(c("devtools", "webshot2", "chromote"))`

4.  **Install Quarto (version > 1.8.25):**

    Follow the instructions on the [Quarto website](https://quarto.org/docs/download/).

5.  **Browser Requirement (for Report Exports)**
     `webshot2` is used to capture static snapshots of interactive plots. It requires a *Chromium-based browser* (Google Chrome, Microsoft Edge, or Chromium) to be installed on your system.
     
     **Windows Users**: If you encounter a "path to Chrome" error, you must manually point R to your browser executable. You can add this to your .Renviron file or run it at the start of your session: `Sys.setenv(CHROMOTE_CHROME = "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe")`
------------------------------------------------------------------------

## 🛠 How to Install the diareport Package

**For development use:**

1.  Clone the GitHub repository.

2.  Use `devtools::install()` to install the diareport package on your system.

**For end user:**

``` r
devtools::install_github('Gevaert-Lab/diareport')
```

------------------------------------------------------------------------

## 📂 Quarto Templates Available

The repository contains six Quarto templates:

-   `dea_base_protein`: Template for protein-level analysis report
-   `dea_base_peptide`: Template for precursor-level analysis reports (No PTMs)
-   `dea_partincluded_protein`: Template for protein-level analysis report including Absent-from-DE analysis
-   `dea_partincluded_peptide`: Template for precursor-level analysis reports including Absent-from-DE analysis (No PTMs)
-   `dea_ev_protein`: Template for precursor-level analysis reports for Extra Vesicles (EV) experiments
-   `dea_ev_peptide`: Template for precursor-level analysis reports for Extra Vesicles (EV) experiments

Specify the template name via the `template_file` parameter in the `render_dia_report` function

------------------------------------------------------------------------

## 🚀 How to Run an Analysis: Base Protein Report ( `dea_base_protein`)

**Important:**  
The output folder must be writable. The function does not overwrite existing
reports unless explicitly allowed; if a file with the same name already exists,
rendering will fail.

1.  Download [DIA_data.zip](https://raw.githubusercontent.com/Gevaert-Lab/DIA-Report/main/example_report/DIA_data.zip), which includes the DIA-NN report and EDF file from [Staes, An, et al.](https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00048?ref=PDF).
2.  Unzip it into a folder called `../path/DIA_data/`.
3.  Use the following code to run the analysis, adjusting the paths accordingly 

**Remark:**  Always use the *full path* to indicate a folder or file.

``` r
library(diareport)

report_target_folder  <- '../path/UPS_spike'
template_file = "dea_base_protein"
output_filename = "Output_report.html"

params <- list(
  title =  "Quantitative Benchmarking DIA",
  subtitle = 'DE Analysis',
  author= 'Your Name',
  description= 'This experiment is designed for DIA benchmarking of different quantitative workflows',
  input_file= '../path/DIA_data/report.tsv',
  design_file = '../path/DIA_data/annotation_DIA.csv',
  folder_prj = report_target_folder,
  contrast= 'Group',
  aggr_method= 'medianPolish',
  normalization ='center.median',
  formula = '~ -1 + Group',
  Proteotypic = TRUE,
  pep_per_prot=2,
  nNonZero= 30,
  FC_thr = 1,
  comparisons= c('GroupA - GroupB', 'GroupA - GroupC','GroupB - GroupC'),
  quantitative_features= 'Precursor.Quantity',
  filtering_contaminant= FALSE,
  filtPerGroup= 'all',
  wildstr_run = 'DIA_Yeast_UPS2_',
  mbr= TRUE,
  DIANN_ver2= FALSE,
  comparison_label  = c("A - B","A - C","B - C")
)

# To run the analysis:
render_dia_report(params, template_file, report_target_folder, output_filename)
```
*Remark*: The fully rendered  HTML report related UPS2 dataset is available online: https://zenodo.org/records/18701402

------------------------------------------------------------------------


## 🚀 How to Run an Analysis: EV Protein Report ( `dea_ev_protein`)


1. Extra Vesicle Proteomics experiments is included in the package  [Pauwels et al.](https://doi.org/10.1002/jev2.70103).

The EV dataset shipped with the package is provided for demonstration and vignette purposes only.

**Remark:**  Always use the *full path* to indicate a folder or file.
**Note:** Ensure `report_target_folder` is writable and does not already contain
an `EV_report.html` file.

``` r
library(diareport)

report_target_folder  <- 'your_path/EV_proteomics'
template_file = "dea_ev_protein"
output_filename = "EV_report.html"

params <- list(
  title =  "EV Proteomics: Sample Preparation Comparison Methods",
  subtitle = 'DE Analysis',
  author= 'Your Name',
  description= 'Same EV material prepared with different protocols: ultracentrifugation (UC) and ultrafiltration with 96-well plate (UF96)',
  input_file =  system.file("extdata/input_data/report.parquet", package = "diareport"),
  design_file = system.file("extdata/input_data/experiment_design.csv", package = "diareport"),
  folder_prj = report_target_folder,
  contrast= 'Group',
  aggr_method= 'medianPolish',
  normalization ='center.median',
  formula = '~ -1 + Group',
  Proteotypic = TRUE,
  pep_per_prot= 2,
  nNonZero= 30,
  FC_thr = 1,
  comparisons= c('GroupUF96 - GroupUC'),
  quantitative_features= 'Precursor.Quantity',
  filtPerGroup= 'at_least_one',
  wildstr_run = 'CMB-',
  mbr= TRUE,
  DIANN_ver2= TRUE,
  comparison_label  = c('UF96 - UC'),
  filtering_contaminant = TRUE,
  contaminant_str = '_BOVIN',
  filt_NaNDE = FALSE
)

# To run the analysis:
render_dia_report(params, template_file, report_target_folder, output_filename)
```

*Remark*: The fully rendered  HTML report is available online: https://zenodo.org/records/18754965

------------------------------------------------------------------------


## ⚙️ Input Parameters

The parameters must be specified in a YAML file. You can find an example in `parameters.yaml`. The parameter descriptions are as follows:

-   **title**: Title of the report
-   **subtitle**: Subtitle
-   **author**: Author name
-   **description**: Description of the experiment
-   **input_file**: Path to the DIA-NN report file (*tsv \| parquet*)
-   **design_file**: Path to the experiment design file
-   **folder_prj**: Path to the root output folder for the analysis
-   **formula**: Formula used in the linear model
-   **contrast**: Name of the column in the experiment design file used in the model (default: Group)
-   **aggr_method**: Summarization method used [*medianPolish, robustSummary(), colMeans(), colMedians(), base::colSums()*]
-   **normalization**: Normalization method [*sum, max, center.mean, center.median, div.mean, div.median, diff.median, quantiles, quantiles.robust, vsn, medianofRatios, LogMedian*]
-   **Proteotypic**: Include only proteotypic peptides (Boolean: TRUE / FALSE)
-   **pep_per_prot**: Number of peptides per protein
-   **nNonZero**: Minimum percentage of samples with non-missing values (used with `filtPerGroup`)
-   **comparisons**: List of comparisons to use
-   **FC_thr**: log2FC threshold (default 1)
-   **adjpval_thr**: Statistical threshold for significant hits (adj.P-value, default 0.05)
-   **ensembl_annotation**: Path to the file containing additional annotations for the proteins
-   **ensembl_col**: List of annotation fields to add; the last column must be a gene symbol
-   **filtering_contaminant**: Enable/disable filtering of contaminants
-   **contaminant_str**: String marking contaminants in the FASTA file (e.g., *Cont*)
-   **cofounder_list**: List of confounder names for analysis
-   **PCA_comparison**: List of confounder names for PCA plots (e.g., *Group-ConfA*)
-   **quantitative_features**: Quantitative feature to use [`Precursor.Quantity` (default option), `Precursor.Normalized`],
-   **filtPerGroup**: Filtering of NaN values based on `nNonZero`, applied to at least one group (`at_least_one`), all groups (`all`), or across all samples (empty string `''`)
-   **mbr**: If MBR is used in DIA-NN, Lib.Q-value and Lib.PG.Q-values are used to select precursors instead of Global.Q-value & Global.PG.Q-value (Boolean: TRUE / FALSE)
-   **wildstr_run**: Wildcard string for run file identification (default: CMB-)
-   **DIANN_ver2**: Set to TRUE if DIA-NN results were generated with version \> 2; otherwise, FALSE
-   **comparison_label**: List of comparison labels without variable names (e.g., *GroupA - GroupB → A - B*)
-   **keep_design_order**: If TRUE, keeps the sample order as in the design file (default: FALSE)

Parameters can also be provided as an R list (see `example_report/Run_DIAReportR.R`) if rendering the report from R.

**Note:** With DIANN v2+, the report is saved in **parquet** format, and the path in `input_file` **must** also include the **'report.protein_description.tsv'** file. If this file is not found, an exception is thrown.

------------------------------------------------------------------------

## 📝 Experiment Design File (EDF)

The experiment design file is a CSV that must include the following columns:

-   **Sample**: Sample name (used in all plots, should be meaningful and not too long)
-   **Run**: Raw file name *without file extension (mzML/.d/.raw)*
-   **Group**: Groups in the experiment (e.g., Cancer/control, mutation/WT)
-   **Replicate**: Label for the sample replicates

Example:

| Sample | Run | Group | Replicate |
|--------------|-------------------------------|-------------|-------------|
| B000250_ratio01_DIA | B000250_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio01_DIA | A | 1 |
| B000254_ratio02_DIA | B000254_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio02_DIA | B | 1 |
| B000258_ratio04_DIA | B000258_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio04_DIA | C | 1 |
| B000262_ratio08_DIA | B000262_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio08_DIA | D | 1 |
| B000266_ratio10_DIA | B000266_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio10_DIA | E | 1 |

**Note:** The names in the "Run" column must match the run names in the DIA-NN report (see the "run" column in the DIA-NN report).

**Note 2:** Matching between "Run" and "Sample" is based on the "Run" column, excluding the file extension. This applies to DIA-NN results from both older and newer versions (≥2).

The EDF file may also include confounder columns, which can be used in confounder analysis, PCA plots, and as fixed effects in the linear model.

------------------------------------------------------------------------


## 🐳 Docker version 

`diareport` is available as a Docker image based on R/RStudio, with all required dependencies pre-installed.

Before proceeding, ensure you have Docker installed and running on your machine. If you haven't installed it yet, you can download it from the official Docker website.

**Quick Start**

- Create a docker-compose.yml file in your working directory with the following content:``` yaml
``` yaml
version: "3.9"
services:
  rstudio:
    image: aarge/gevaertlab_diareport:latest   
    container_name: rstudio-diareport
    environment:
      - PASSWORD=changeme
      - CHROMOTE_CHROME=/usr/bin/google-chrome
    ports:
      - "8787:8787"
    volumes:
      - "./:/home/rstudio/work"
    shm_size: '1gb'
```
**Note on File Access**: 
The line - "./:/home/rstudio/work" mounts your current local directory to the container. 
This allows the RStudio instance inside the container to access your local files. Any scripts or data you save into the /home/rstudio/work directory while in RStudio will be immediately accessible on your local machine, and your work will persist even after you stop or remove the container.

- Run the following command in your terminal:
```bash
docker compose up -d
```
- Once the container is running, open your web browser and navigate to:
    - 👉 [http://localhost:8787](http://localhost:8787)  
    - **Username:** `rstudio`  
    - **Password:** `changeme` (or whatever was set in `docker-compose.yml`)  


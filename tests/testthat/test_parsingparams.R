library(testthat)
library(SummarizedExperiment)

test_that("parsing_base", {

  params <- list()
  params$design_file<-  testthat::test_path("ExperimentalDesignCMB-1536.csv")
  params$input_file <- testthat::test_path("dfmsrob_wide.Rds")
  params$mbr <- TRUE
  params$description<- ""
  params$title<- ""
  params$subtitle<- ""
  params$author<- ""
  params$aggr_method<- "medianPolish"
  params$normalization<- "quantiles"
  params$filtering_contaminant <- FALSE
  params$contaminant_str <- ''
  params$ensembl_annotation <- ""
  params$ensembl_col<- ""
  params$confounder_list<- ""
  params$pep_per_prot<- 3
  params$nNonZero <- 30
  params$filtPerGroup <- 'all'
  params$DIANN_ver2 <- TRUE
  params$Proteotypic <- TRUE
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$comparison_label <- c('A vs. B')
  params$keep_design_order <- FALSE
  params$contrast <- 'Group'
  params$formula <- '~ -1 + Group'
  params$FC_thr <- 2
  params$adjpval_thr <- 0.05


  res <- validate_params(params )

  # Perform the tests
  expect_true(res)

} )



test_that("parsing_I", {

  params <- list()
  params$design_file<-  testthat::test_path("ExperimentalDesignCMB-1536.csv")
  params$input_file <- testthat::test_path("dfmsrob_wide.Rds")
  params$mbr <- TRUE
  params$description<- ""
  params$title<- ""
  params$subtitle<- ""
  params$author<- ""
  params$aggr_method<- "goofy"
  params$normalization<- "quantiles"
  params$filtering_contaminant <- FALSE
  params$contaminant_str <- ''
  params$folder_prj<- ""
  params$ensembl_annotation <- ""
  params$ensembl_col<- ""
  params$confounder_list<- ""
  params$pep_per_prot<- 3
  params$nNonZero <- 30
  params$filtPerGroup <- 'all'
  params$DIANN_ver2 <- TRUE
  params$Proteotypic <- TRUE
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$comparison_label <- c('A vs. B')
  params$keep_design_order <- FALSE
  params$contrast <- 'Group'
  params$formula <- '~ -1 + Group'
  params$FC_thr <- 2
  params$adjpval_thr <- 0.05


  #res <- validate_params(params )

  # Perform the tests
  expect_error(validate_params(params ), 'Aggregation method must be one of: medianPolish, RobustSummary, colMeans, colMedians.')

} )



test_that("parsing_II", {

  params <- list()
  params$design_file<-  testthat::test_path("ExperimentalDesignCMB-1536.csv")
  params$input_file <- testthat::test_path("XXX.csv")
  params$mbr <- TRUE
  params$description<- ""
  params$title<- ""
  params$subtitle<- ""
  params$author<- ""
  params$aggr_method<- "medianPolish"
  params$normalization<- "quantiles"
  params$filtering_contaminant <- FALSE
  params$contaminant_str <- ''
  params$folder_prj<- ""
  params$ensembl_annotation <- ""
  params$ensembl_col<- ""
  params$confounder_list<- ""
  params$pep_per_prot<- 3
  params$nNonZero <- 30
  params$filtPerGroup <- 'all'
  params$DIANN_ver2 <- TRUE
  params$Proteotypic <- TRUE
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$comparison_label <- c('A vs. B')
  params$keep_design_order <- FALSE
  params$contrast <- 'Group'
  params$formula <- '~ -1 + Group'
  params$FC_thr <- 2
  params$adjpval_thr <- 0.05


  #res <- validate_params(params )

  # Perform the tests
  expect_error(validate_params(params ), "Input file does not exist or is not specified.")

} )



test_that("parsing_III", {

  params <- list()
  params$design_file<-  testthat::test_path("ExperimentalDesignCMB-1536.csv")
  params$input_file <- testthat::test_path("dfmsrob_wide.Rds")
  params$mbr <- TRUE
  params$description<- ""
  params$title<- ""
  params$subtitle<- ""
  params$author<- ""
  params$aggr_method<- "medianPolish"
  params$normalization<- "quantiles"
  params$filtering_contaminant<- FALSE
  params$contaminant_str<- ""
  params$folder_prj<- ""
  params$ensembl_annotation <- testthat::test_path("dfmsrob_wide.Rds")
  params$ensembl_col<- ""
  params$confounder_list<- ""
  params$pep_per_prot<- 3
  params$nNonZero <- 30
  params$filtPerGroup <- 'all'
  params$DIANN_ver2 <- TRUE
  params$Proteotypic <- TRUE
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$comparison_label <- c('A vs. B')
  params$keep_design_order <- FALSE
  params$contrast <- 'Group'
  params$formula <- '~ -1 + Group'
  params$FC_thr <- 2
  params$adjpval_thr <- 0.05

  #res <- validate_params(params )

  # Perform the tests
  expect_error(validate_params(params ), "When provide an ensembl_annotation file, you must provide a non-empty 'ensembl_col'")

} )



test_that("parsing_IV", {

  params <- list()
  params$design_file<-  testthat::test_path("ExperimentalDesignCMB-1536.csv")
  params$input_file <- testthat::test_path("dfmsrob_wide.Rds")
  params$mbr <- TRUE
  params$description<- ""
  params$title<- ""
  params$subtitle<- ""
  params$author<- ""
  params$aggr_method<- "medianPolish"
  params$normalization<- "quantiles"
  params$filtering_contaminant<- FALSE
  params$contaminant_str<- ""
  params$folder_prj<- ""
  params$ensembl_annotation <- ""
  params$ensembl_col<- c('A','B','hgnc_symbol')
  params$confounder_list<- ""
  params$pep_per_prot<- 3
  params$nNonZero <- 30
  params$filtPerGroup <- 'all'
  params$DIANN_ver2 <- TRUE
  params$Proteotypic <- TRUE
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$comparison_label <- c('A vs. B')
  params$keep_design_order <- FALSE
  params$contrast <- 'Group'
  params$formula <- '~ -1 + Group'
  params$FC_thr <- 2
  params$adjpval_thr <- 0.05

  #res <- validate_params(params )

  # Perform the tests
  expect_error(validate_params(params ),
               "Please provide an ensembl_annotation file, you  provided a non-empty 'ensembl_col' with a non valid 'ensembl_annotation' file")

} )




test_that("parsing_V", {

  params <- list()
  params$design_file<-  testthat::test_path("ExperimentalDesignCMB-1536.csv")
  params$input_file <- testthat::test_path("dfmsrob_wide.Rds")
  params$mbr <- TRUE
  params$description<- ""
  params$title<- ""
  params$subtitle<- ""
  params$author<- ""
  params$aggr_method<- "medianPolish"
  params$normalization<- "quantiles"
  params$filtering_contaminant<- TRUE
  params$contaminant_str<- ""
  params$folder_prj<- ""
  params$ensembl_annotation <- ""
  params$ensembl_col<- ""
  params$confounder_list<- ""
  params$pep_per_prot<- 3
  params$nNonZero <- 30
  params$filtPerGroup <- 'all'
  params$DIANN_ver2 <- TRUE
  params$Proteotypic <- TRUE
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$comparison_label <- c('A vs. B')
  params$keep_design_order <- FALSE
  params$contrast <- 'Group'
  params$formula <- '~ -1 + Group'
  params$FC_thr <- 2
  params$adjpval_thr <- 0.05
  params$contaminant_str <- ''

  #res <- validate_params(params )

  # Perform the tests
  expect_error(validate_params(params ),
               "When 'filtering_contaminant' is TRUE, you must provide a non-empty 'contaminant_str'")
} )


test_that("parsing_VI", {

  params <- list()
  params$design_file<-  testthat::test_path("ExperimentalDesignCMB-1536.csv")
  params$input_file <- testthat::test_path("dfmsrob_wide.Rds")
  params$mbr <- TRUE
  params$description<- ""
  params$title<- ""
  params$subtitle<- ""
  params$author<- ""
  params$aggr_method<- "medianPolish"
  params$normalization<- "quantiles"
  params$filtering_contaminant<- FALSE
  params$contaminant_str <- "cont__"
  params$folder_prj<- ""
  params$ensembl_annotation <- ""
  params$ensembl_col<- ""
  params$confounder_list<- ""
  params$pep_per_prot<- 3
  params$nNonZero <- 30
  params$filtPerGroup <- 'all'
  params$DIANN_ver2 <- TRUE
  params$Proteotypic <- TRUE
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$comparison_label <- c('A vs. B')
  params$keep_design_order <- FALSE
  params$contrast <- 'Group'
  params$formula <- '~ -1 + Group'
  params$FC_thr <- 2
  params$adjpval_thr <- 0.05

  #res <- validate_params(params )

  # Perform the tests
  expect_error(validate_params(params ),
               "If 'filtering_contaminant' is FALSE , you are not allowed to  provide a 'contaminant_str'")
} )



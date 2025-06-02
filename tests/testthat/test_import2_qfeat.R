library(testthat)
library(dplyr)
library(SummarizedExperiment)


test_that("qfeat_base", {

   params <- list()
   params$design_file<-  testthat::test_path("ExperimentalDesignCMB-1536.csv")
   params$mbr <- 'TRUE'
   params$quantitative_features <- 'Precursor.Quantity'
   params$wildstr_run  <- 'CMB-'
   params$comparisons <- c('GroupA - GroupB')
   params$keep_design_order <- FALSE
   lst_wide_columns <- c('Run', 'Precursor.Id', 'Modified.Sequence', 'Stripped.Sequence', 'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'Proteotypic')
   min_col_need_design <- c("Sample","Run", "Group", "Replicate")
   colname__ <- c("Precursor.Id" , "Modified.Sequence","Stripped.Sequence","Protein.Group",
                  "Protein.Ids","Protein.Names","Genes","Proteotypic","First.Protein.Description")

   dfMsqrob <- readRDS(  testthat::test_path( "dfmsrob_wide.Rds"))

   L <- readLines(params$design_file, n = 1)
   if (grepl(";", L)) design <- read.csv2(params$design_file) else design <- read.csv( params$design_file)
   out <- import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = colname__  )
  # Perform the tests
  #print(import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = lst_wide_columns  )$error)
  expect_equal(out$status,0) # Ensure the data is loaded
  expect_equal(dim(out$q_feat[['precursor']])[1],25801)
  expect_equal(dim(out$q_feat[['precursor']])[2],6)

} )


test_that("qf_run_name_wrong", {
  # Define the function to test

  params <- list()

  params$design_file<-   testthat::test_path("ExperimentalDesignCMB-1536_wrongA.csv")
  params$mbr <- 'TRUE'
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$keep_design_order <- FALSE
  lst_wide_columns <- c('Run', 'Precursor.Id', 'Modified.Sequence', 'Stripped.Sequence', 'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'Proteotypic')
  min_col_need_design <- c("Sample","Run", "Group", "Replicate")
  colname__ <- c("Precursor.Id" , "Modified.Sequence","Stripped.Sequence","Protein.Group",
                 "Protein.Ids","Protein.Names","Genes","Proteotypic","First.Protein.Description")
  #
  dfMsqrob <- readRDS(  testthat::test_path( "dfmsrob_wide.Rds"))
  #
  L <- readLines(params$design_file, n = 1)
  if (grepl(";", L)) design <- read.csv2(params$design_file) else design <- read.csv( params$design_file)
  out <- import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = colname__  )
  # Perform the tests
  #print(import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = lst_wide_columns  )$error)
  expect_equal(out$status,1) # Ensure the data is loaded

})

test_that("qf_run_less_sample_design", {
  # Define the function to test

  params <- list()
  params$design_file<-  testthat::test_path("ExperimentalDesignCMB-1536_wrongB.csv")

  params$mbr <- 'TRUE'
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$keep_design_order <- FALSE
  lst_wide_columns <- c('Run', 'Precursor.Id', 'Modified.Sequence', 'Stripped.Sequence', 'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'Proteotypic')
  min_col_need_design <- c("Sample","Run", "Group", "Replicate")
  colname__ <- c("Precursor.Id" , "Modified.Sequence","Stripped.Sequence","Protein.Group",
                 "Protein.Ids","Protein.Names","Genes","Proteotypic","First.Protein.Description")
  #
  dfMsqrob <- readRDS(  testthat::test_path( "dfmsrob_wide.Rds"))
  #
  L <- readLines(params$design_file, n = 1)
  if (grepl(";", L)) design <- read.csv2(params$design_file) else design <- read.csv( params$design_file)
  out <- import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = colname__  )
  # Perform the tests
  #print(import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = lst_wide_columns  )$error)
  expect_equal(out$status,0) # Ensure the data is loaded
  expect_equal(dim(out$q_feat[['precursor']])[2],5)

})




test_that("qf_run_less_sample_diann", {
  # Define the function to test
  params <- list()
  params$design_file<-    testthat::test_path( "ExperimentalDesignCMB-1536.csv")
  params$mbr <- 'TRUE'
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$keep_design_order <- FALSE
  lst_wide_columns <- c('Run', 'Precursor.Id', 'Modified.Sequence', 'Stripped.Sequence', 'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'Proteotypic')
  min_col_need_design <- c("Sample","Run", "Group", "Replicate")
  colname__ <- c("Precursor.Id" , "Modified.Sequence","Stripped.Sequence","Protein.Group",
                 "Protein.Ids","Protein.Names","Genes","Proteotypic","First.Protein.Description")
  #
  dfMsqrob <- readRDS(  testthat::test_path( "dfmsrob_wide.Rds"))
  dfMsqrob <- dfMsqrob %>% select(- `CMB-1536DIA_UPRC_9606X01280_2p_Ih25cm2_PM2_CMB-1536_KRGEV_10` )

  #
  L <- readLines(params$design_file, n = 1)
  if (grepl(";", L)) design <- read.csv2(params$design_file) else design <- read.csv( params$design_file)
  out <- import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = colname__  )
  # Perform the tests
  expect_equal(out$status,1) # Ensure the data is loaded

})

test_that("qf_wrong_confounder", {
  params <- list()

  params$design_file<-   testthat::test_path( "ExperimentalDesignCMB-1536.csv")
  params$mbr <- 'TRUE'
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$confounder_list <- c('goofy','Age')
  dfMsqrob <- readRDS(  testthat::test_path( "dfmsrob_wide.Rds"))
  min_col_need_design <- c("Sample","Run", "Group", "Replicate")
  colname__ <- c("Precursor.Id" , "Modified.Sequence","Stripped.Sequence","Protein.Group",
                 "Protein.Ids","Protein.Names","Genes","Proteotypic","First.Protein.Description")
  L <- readLines(params$design_file, n = 1)
  if (grepl(";", L)) design <- read.csv2(params$design_file) else design <- read.csv(params$design_file)
  out <- import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = colname__  )
  expect_equal(out$status,1)
  #expect_equal(out$error, 'confounder values are not present in design file')
} )

# Run all tests
#test_dir(".", reporter = "summary")

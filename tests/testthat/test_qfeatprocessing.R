library(testthat)
library(dplyr)
library(SummarizedExperiment)




test_that("addrow_information", {

  params <- list()
  params$design_file<-  testthat::test_path('ExperimentalDesignCMB-1536.csv')
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
  dfMsqrob <- readRDS( testthat::test_path( "dfmsrob_wide.Rds"))
  #
  L <- readLines(params$design_file, n = 1)
  if (grepl(";", L)) design <- read.csv2(params$design_file) else design <- read.csv( params$design_file)
  out <- import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = colname__  )
  g <- add_rowdata_detection(out$q_feat,  design , assay = 'precursor' )

  # Perform the tests
  expect_true(all(  c('percA','percB','pep_per_prot','ZeroA','ZeroB') %in% colnames( as.data.frame(   rowData(g[['precursor']]) )  ) ))

})



test_that("filteringNA_qfeat", {

  params <- list()
  params$design_file<- testthat::test_path('ExperimentalDesignCMB-1536.csv')
  params$mbr <- 'TRUE'
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$filtPerGroup <- ''
  params$Proteotypic <- TRUE
  params$filtering_contaminant <- FALSE
  params$pep_per_prot <- 2
  params$nNonZero <- 30
  params$keep_design_order <- FALSE
  lst_wide_columns <- c('Run', 'Precursor.Id', 'Modified.Sequence', 'Stripped.Sequence', 'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'Proteotypic')
  min_col_need_design <- c("Sample","Run", "Group", "Replicate")
  colname__ <- c("Precursor.Id" , "Modified.Sequence","Stripped.Sequence","Protein.Group",
                 "Protein.Ids","Protein.Names","Genes","Proteotypic","First.Protein.Description")
  #
  dfMsqrob <-  readRDS( testthat::test_path( "dfmsrob_wide.Rds"))
  #
  L <- readLines(params$design_file, n = 1)
  if (grepl(";", L)) design <- read.csv2(params$design_file) else design <- read.csv( params$design_file)
  out <- import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = colname__  )

  g <- add_rowdata_detection(out$q_feat,  design , assay = 'precursor' )
  res_filt  <-  filteringNA_qfeat(pe_ = g, params, design)


  expect_equal( dim(res_filt$q_feat[['precursor']])[1]  ,  19965    )

})

test_that("procesProtein_qfeat", {

  params <- list()
  params$design_file<- testthat::test_path('ExperimentalDesignCMB-1536.csv')
  params$mbr <- 'TRUE'
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$filtPerGroup <- ''
  params$Proteotypic <- TRUE
  params$filtering_contaminant <- FALSE
  params$pep_per_prot <- 2
  params$nNonZero <- 30
  params$keep_design_order <- FALSE
  lst_wide_columns <- c('Run', 'Precursor.Id', 'Modified.Sequence', 'Stripped.Sequence', 'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'Proteotypic')
  min_col_need_design <- c("Sample","Run", "Group", "Replicate")
  colname__ <- c("Precursor.Id" , "Modified.Sequence","Stripped.Sequence","Protein.Group",
                 "Protein.Ids","Protein.Names","Genes","Proteotypic","First.Protein.Description")
  #
  dfMsqrob <-  readRDS( testthat::test_path( "dfmsrob_wide.Rds"))
  #
  L <- readLines(params$design_file, n = 1)
  if (grepl(";", L)) design <- read.csv2(params$design_file) else design <- read.csv( params$design_file)
  out <- import2_qfeature (dfMsqrob, design, params, min_col_need_design, diann_colname = colname__  )

  g <- add_rowdata_detection(out$q_feat,  design , assay = 'precursor' )
  res_filt  <-  filteringNA_qfeat(pe_ = g, params, design)
  suppressWarnings(
  processing_res <-  processing_qfeat_protein(   res_filt$q_feat , params,  MsCoreUtils::medianPolish )
  )
  expect_in( names(processing_res$q_feat)  ,  c('precursor', 'precursorLog', 'precursorNorm', 'proteinRS')    )
  expect_equal( dim(processing_res$q_feat[['proteinRS']])[1] , 2454 )

})

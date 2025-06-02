library(testthat)
library(msqrob2)
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)


test_that("msqrob2_base", {

  params <- list()
  params$design_file<-  testthat::test_path("ExperimentalDesignCMB-1536.csv")
  params$mbr <- 'TRUE'
  params$quantitative_features <- 'Precursor.Quantity'
  params$wildstr_run  <- 'CMB-'
  params$comparisons <- c('GroupA - GroupB')
  params$comparison_label <- c('A vs. B')

  params$keep_design_order <- FALSE
  params$contrast <- 'Group'
  params$formula <- '~ -1 + Group'
  params$FC_thr <- 2
  params$adjpval_thr <- 0.05
  params$ensembl_annotation <- ''

  lst_wide_columns <- c('Run', 'Precursor.Id', 'Modified.Sequence', 'Stripped.Sequence', 'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'Proteotypic')
  min_col_need_design <- c("Sample","Run", "Group", "Replicate")
  colname__ <- c("Precursor.Id" , "Modified.Sequence","Stripped.Sequence","Protein.Group",
                 "Protein.Ids","Protein.Names","Genes","Proteotypic","First.Protein.Description")

  pe <- readRDS(  testthat::test_path( "qfeat_proc&filt.Rds"))

  #L <- readLines(params$design_file, n = 1)
  #if (grepl(";", L)) design <- read.csv2(params$design_file) else design <- read.csv( params$design_file)
  suppressWarnings(

  out <- msqrob_model (pe, params,layer = 'proteinRS')
  )

  # Perform the tests
  expect_equal(length(out$de_comp),1) # Ensure the data is loaded
  expect_equal(names(out$de_comp) , 'A vs. B')
  expect_true(is_ggplot(out$de_comp$`A vs. B`$volcano))
  expect_equal(out$status , 0)
  expect_equal(  dim(out$de_comp$`A vs. B`$toptable)[1],1486)

} )

#' @author Andrea Argentini
#' @title  step_read_diann
#'
#' @description
#' Step function that reads the DIANN report and populates the workflow input list
#' with the design and dfMsqrob objects.
#'
#' @param input A list containing at least params_report.
#' @return The updated input list with design and dfMsqrob.
step_read_diann <- function(input) {
  res <- read_DIANN_report(params = input$params_report)
  if (res$status == 1) stop(res$error)
  input$design <- res$df_design
  input$dfMsqrob <- res$w_diann
  input
}

#' @author Andrea Argentini
#' @title  step_import_qfeat
#'
#' @description
#' Step function that wraps the import2_qfeature function. It processes input data and
#' adds the q_feat object to the input list.
#'
#' @param input A list containing dfMsqrob, design, and params_report.
#' @return The updated input list with the pe (q_feat) object.
step_import_qfeat <- function(input) {
  min_info_design <- c("Sample", "Run", "Group", "Replicate")
  diann_colname <- c("Precursor.Id", "Modified.Sequence", "Stripped.Sequence",
                     "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                     "Proteotypic", "First.Protein.Description")
  res <- import2_qfeature(input$dfMsqrob, input$design, params = input$params_report,
                          min_info_design, diann_colname = diann_colname)
  if (res$status == 1) stop(res$error)
  input$pe <- res$q_feat
  input
}

#' @author Andrea Argentini
#' @title  step_add_rowdata_precursor
#'
#' @description
#' Step function that adds row data detection at the precursor level using add_rowdata_detection.
#'
#' @param input A list containing pe and design.
#' @return The updated input list with pe modified at precursor level.
step_add_rowdata_precursor <- function(input) {
  input$pe <- add_rowdata_detection(input$pe, input$design, assay = 'precursor')
  input
}

#' @author Andrea Argentini
#' @title  step_filter_na
#'
#' @description
#' Step function that applies filteringNA_qfeat to filter missing values in the qfeature object.
#'
#' @param input A list containing pe, params_report, and design.
#' @return The updated input list with filtered pe.
step_filter_na <- function(input) {
  res <- filteringNA_qfeat(input$pe, params = input$params_report, input$design)
  if (res$status == 1) stop(res$error)
  input$pe <- res$q_feat
  input
}

#' @author Andrea Argentini
#' @title  step_processing_protein
#'
#' @description
#' Step function that aggregates qfeature object to protein level using the specified aggregation method.
#'
#' @param input A list containing pe, params_report, and aggr_method_f.
#' @return The updated input list with pe aggregated at protein level.
step_processing_protein <- function(input) {
  aggr_method_f <- input$aggr_method_f
  res <- processing_qfeat_protein(input$pe, params = input$params_report, aggr_method_f)
  if (res$status == 1) stop(res$error)
  input$pe <- res$q_feat
  input
}


#' @author Andrea Argentini
#' @title  step_processing_protein
#'
#' @description
#' Step function that aggregates precursor to peptide sequence level
#'
#' @param input A list containing pe, params_report, and aggr_method_f.
#' @return The updated input list with pe aggregated at protein level.
step_processing_peptide<- function(input) {
  aggr_method_f <- input$aggr_method_f
  res <- processing_qfeat_peptide(input$pe, params = input$params_report, aggr_method_f)
  if (res$status == 1) stop(res$error)
  input$pe <- res$q_feat
  input
}

#' @author Andrea Argentini
#' @title  step_add_rowdata_protein
#'
#' @description
#' Step function that adds row data detection at the proteinRS level using add_rowdata_detection.
#'
#' @param input A list containing pe and design.
#' @return The updated input list with pe modified at proteinRS level.
step_add_rowdata_protein <- function(input) {
  input$pe <- add_rowdata_detection(input$pe, input$design, assay = 'proteinRS')
  input
}


#' @author Andrea Argentini
#' @title  step_add_rowdata_protein
#'
#' @description
#' Step function that adds row data detection at the proteinRS level using add_rowdata_detection.
#'
#' @param input A list containing pe and design.
#' @return The updated input list with pe modified at proteinRS level.
step_add_rowdata_peptide <- function(input) {
  input$pe <- add_rowdata_detection(input$pe, input$design, assay = 'peptideNorm')
  input
}

#' @author Andrea Argentini
#' @title  step_add_ensembl
#'
#' @description
#' Step function that adds Ensembl annotation to the qfeature object if annotation is specified in params_report.
#'
#' @param input A list containing pe and params_report.
#' @return The updated input list with Ensembl annotation (if applied).
step_add_ensembl <- function(input) {
  if (!is.null(input$params_report$ensembl_annotation) && input$params_report$ensembl_annotation != '') {
    res <- add_ensembl_annotation(input$pe, params = input$params_report)
    if (res$status == 0) {
      input$pe <- res$q_feat
    } else {
      stop(res$error)
    }
  }
  input
}

#' @author Andrea Argentini
#' @title  step_msqrob_de
#'
#' @description
#' Step function that runs msqrob_model and adds differential expression results to the input list.
#'
#' @param input A list containing pe and params_report.
#' @return The updated input list with de_comparison and updated pe.
step_msqrob_de <- function(input) {
  print(input$layer)
  res <- msqrob_model(input$pe, input$params_report,input$layer )
  if (res$status == 0) {
    input$pe <- res$q_feat
    input$de_comparison <- res$de_comp
  } else {
    stop(res$error)
  }
  input
}

#' @author Andrea Argentini
#' @title  step_partial_result
#'
#' @description
#' Step function for pariati data
#'
#' @param input A list containing pe and params_report.
#' @return The updated input list with de_comparison and updated pe.
step_partial_result <- function(input) {
  print(input$layer)
  res <- partially_present(input$pe, input$params_report,input$layer )
  if (res$status == 0) {
    input$part_item <- res$part_item
    input$part_value <- res$part_value
  } else {
    stop(res$error)
  }
  input
}


#' @author Andrea Argentini
#' @title  run_workflow_pipeline
#'
#' @description
#' Executes a sequence of workflow step functions in order, passing the result of each step
#' as input to the next. Each step function should accept a single list argument (the current
#' workflow state) and return a list with updated state. This function is suitable for modular,
#' extensible workflow pipelines, such as those used in DIA report processing.
#'
#' @param initial_input A named list containing the initial state or data required by the pipeline.
#' @param steps A list of step functions, each accepting a single list argument and returning a list.
#' @param ... Additional arguments to pass to each step function.
#'
#' @return The final state list returned by the last step function in the pipeline.
#'
#' @details
#' If any step function returns a list with a `status` element set to 1, the pipeline
#' stops with an error using the `error` element of that list.
#'

run_workflow_pipeline <- function(initial_input, steps, ...) {
  result <- initial_input
  for (step in steps) {
    result <- step(result, ...)
    # Optionally, check for errors or status here if your functions return status lists
    if (is.list(result) && !is.null(result$status) && result$status == 1) {
      stop(result$error)
    }
  }
  return(result)
}

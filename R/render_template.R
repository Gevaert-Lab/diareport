#' Validate file name
#'
#' @param filename The file name to be validated.
#' @return TRUE if the file name is valid, otherwise stops with an error message.
#' @importFrom assertthat assert_that is.string
validate_filename <- function(filename) {
  # Define invalid characters for file names
  invalid_chars <- "[<>:\"/\\|?*]"

  # Check if filename is a string
  assertthat::assert_that(assertthat::is.string(filename), msg = "filename must be a string.")

  # Check if filename length is less than 40 characters
  if (nchar(filename) > 40) {
    stop("The file name must be less than 40 characters.")
  }

  # Check if filename contains invalid characters
  if (grepl(invalid_chars, filename)) {
    stop("The file name contains invalid characters. Invalid characters are: <>:\"/\\|?*")
  }

  TRUE
}

#' Validate report folder path
#'
#' @param report_folder The folder path to be validated.
#' @return TRUE if the folder path is valid and writable, otherwise stops with an error message.
#' @importFrom assertthat assert_that is.writeable
validate_folder <- function(report_folder) {
  # Define invalid characters for Windows file system
  invalid_chars <- "[<>:\"/\\|?*]"

  # # Check if the folder path contains invalid characters
  # if (grepl(invalid_chars, report_folder)) {
  #   stop("The folder path contains invalid characters. Invalid characters are: <>:\"/\\|?*")
  # }

  if (!dir.exists(file.path(report_folder))) {
    dir.create(file.path( report_folder),recursive = TRUE)
  }
  dir.create(file.path( report_folder, "Result"),recursive = TRUE)

  # Check if the folder path is writable
  assertthat::assert_that(assertthat::is.writeable(report_folder), msg = "The folder path is not writable.")

  TRUE
}



#' Validate template parameter
#'
#' @param template The template parameter to be validated.
#' @return TRUE if the template is valid, otherwise stops with an error message.
#' @importFrom assertthat assert_that is.string
validate_template <- function(template) {
  # Define the list of valid templates
  valid_templates <- c( "Template_DIA-NN_peptide_dev.qmd",
                        "Template_DIA-NN_dev.qmd",
                        "Template_DIA-NN_dev_A.qmd",
                        "Template_DIA-NN_dev_EV.qmd",
                        "Template_DIA-NN_peptide_dev_A.qmd")

  # Check if template is a string
  assertthat::assert_that(assertthat::is.string(template), msg = "template must be a string.")

  # Check if template belongs to the list of valid templates
  if (!template %in% valid_templates) {
    stop("Invalid template. The template must be one of the following: ", paste(valid_templates, collapse = ", "))
  }

  TRUE
}


#' Validate parameters for the DIA-NN report
#'
#' @param params A list of parameters to be validated.
#' @return TRUE if all parameters are valid, otherwise stops with an error message.
#' @importFrom assertthat assert_that is.string
validate_params <- function(params) {

  # is_empty <- function(x) {
  #   return(length(x) == 0 || (length(x) == 1 && x == ''))
  # }

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  is_empty <- function(x) {
    is.character(x) && length(x) == 1 && x == ''
  }

  check_formula <- function(x) {
    tryCatch({
      as.formula(x)
      TRUE
    }, error = function(e) {
      FALSE
    })
  }

  check_path <- function(x) {
    if (length(x) == 1 && x == '' ){
      TRUE
    }else{
      file.exists(x)
    }
  }
  requirements <- list(
    input_file = list(
      type = "string",
      check = function(x) file.exists(x),
      msg = "Input file does not exist or is not specified."
    ),
    design_file = list(
      type = "string",
      check = function(x) file.exists(x),
      msg = "Design file does not exist or is not specified."
    ),
    description= list(
      type = "string"
    ),
    title= list(
      type = "string"
    ),
    subtitle = list(
      type = "string"
    ),
    author = list(
      type = "string"
    ),
    contrast = list(
      type = "string",
      check = function(x) !is_empty(x),
      msg = "Contrast must not be empty."
    ),
    comparisons = list(
      type = "character",
      check = function(x) length(x) >= 1 && all(x != ''),
      msg = "Comparisons must contain at least one value."
    ),
    comparison_label = list(
      type = "character",
      check = function(x) length(x) >= 1 && all(x != ''),
      msg = "Comparison label must contain at least one value."
    ),
    formula = list(
      type = "string",
      check = check_formula,
      msg = "Formula is invalid, please check its syntax."
    ),
    aggr_method = list(
      type = "string",
      check = function(x) x %in% c("medianPolish", "RobustSummary", "colMeans", "colMedians"),
      msg = "Aggregation method must be one of: medianPolish, RobustSummary, colMeans, colMedians."
    ),
    normalization = list(
      type = "string",
      check = function(x) x %in%  c("sum", "max", "center.mean", "center.median", "div.mean", "div.median", "diff.meda", "quantiles","quantiles.robust","vsn"),
      msg = "Normalization must be a recognized method among: [sum max center.mean center.median div.mean div.median diff.meda quantiles quantiles.robust vsn]"
    ),
    FC_thr = list(
      type = "numeric",
      check = function(x) x > 0,
      msg = "FC_thr must be a positive number."
    ),
    pep_per_prot= list(
      type = "numeric",
      check = function(x) x > 0,
      msg = "pep_per_prot must be a positive number."
    ),
    nNonZero= list(
      type = "numeric",
      check = function(x) x > 0,
      msg = "nNonZero must be a positive number."
    ),
    adjpval_thr = list(
      type = "numeric",
      check = function(x) x > 0 && x <= 1,
      msg = "adjpval_thr must be between 0 and 1."
    ),
    filtPerGroup = list(
      type = "string",
      check = function(x) x %in% c('','all','at_least_one'),
      msg = "filtPerGroup must be '', 'all', or 'at_least_one'."
    ),
    wildstr_run = list(
      type = 'string',
      check = function(x) !is_empty(x),
      msg = "wildstr_run must not be empty."
    ),
    DIANN_ver2= list(
      type= 'logical',
      check = function(x) is.logical(x),
      msg = "DIANN_ver2 must be logical (TRUE or FALSE)."
    ),
    keep_design_order= list(
      type= 'logical',
      check = function(x) is.logical(x),
      msg = "keep_design_order must be logical (TRUE or FALSE)."
    ),
    mbr = list(
      type= 'logical',
      check = function(x) is.logical(x),
      msg = "mbr must be logical (TRUE or FALSE)."
    ),
    Proteotypic = list(
      type= 'logical',
      check = function(x) is.logical(x),
      msg = "Proteotypic must be logical (TRUE or FALSE)."
    ),
    ensembl_annotation  = list(
      type= 'string',
      check= function(x) check_path(x),
      msg = "Ensembl annotation file does not exist (if specified)."
    ),
    ensembl_col = list(
      type = 'character',
      check = function(x) any(x %in% c("hgnc_symbol") ) ||  params$ensembl_col == "" ,
      msg = "hgnc_symbol columns must be indicated"
    ),
    filtering_contaminant= list(type= 'logical',
      check = function(x) is.logical(x),
      msg = "filtering_contaminant must be logical (TRUE or FALSE)." ),
    contaminant_str = list(type= 'string'),
    confounder_list = list(type = 'character'),
    filt_NaNDE= list(
      type= 'logical',
      check = function(x) is.logical(x),
      msg = "Include NaN from DE analysis (TRUE or FALSE)."
    )
    # Add more as needed
  )


  for (p in names(requirements)) {
    val <- params[[p]]
    req <- requirements[[p]]
    # Type check

    if (req$type == "string") {
      assertthat::assert_that(assertthat::is.string(val), msg = paste0("'", p, "' must be a string."))
    } else if (req$type == "numeric") {
      assertthat::assert_that(is.numeric(val), msg = paste0("'", p, "' must be numeric."))
    } else if (req$type == "logical") {
      assertthat::assert_that(is.logical(val), msg = paste0("'", p, "' must be logical (TRUE/FALSE)."))
    } else if (req$type == "character") {
      assertthat::assert_that(is.character(val), msg = paste0("'", p, "' must be a character vector."))
    } # Add more types as needed

    # Value check (if provided)
    if (!is.null(req$check)) {
      assertthat::assert_that(req$check(val), msg = req$msg %||% paste0("Invalid value for '", p, "'.") )
    }
  }
  # to be added
  if (!is.null(params$filtering_contaminant) && isTRUE(params$filtering_contaminant)) {
    if (is.null(params$contaminant_str) || params$contaminant_str == "") {
      stop("When 'filtering_contaminant' is TRUE, you must provide a non-empty 'contaminant_str'")
    }
  }

  if ( ! is.null(params$contaminant_str) && ! is_empty(params$contaminant_str) ) {
    if (is.null(params$filtering_contaminant) || params$filtering_contaminant == FALSE) {
      stop("If 'filtering_contaminant' is FALSE , you are not allowed to  provide a 'contaminant_str'")
    }
  }


  if (!is.null(params$ensembl_annotation) && (! is_empty(params$ensembl_annotation)) ) {
    if (is.null(params$ensembl_col) || all(params$ensembl_col == "")) {
      stop("When provide an ensembl_annotation file, you must provide a non-empty 'ensembl_col'")
    }
  }


  if (!is.null(params$ensembl_col) && (! is_empty(params$ensembl_col)) ) {
    if (is.null(params$ensembl_annotation) || params$ensembl_annotation == "") {
      stop("Please provide an ensembl_annotation file, you  provided a non-empty 'ensembl_col' with a non valid 'ensembl_annotation' file")
    }
  }

  TRUE
}

#' @title merge_default_parameters
#'
#' @param params_int parameters
#' @return merged parameters
#' @importFrom yaml read_yaml

merge_default_parameters <- function  ( params_int  ){

  yaml_path <- system.file("config", "default_parameter.yaml", package = "diareport")
  default_p <- read_yaml(yaml_path )

  miss <- base::setdiff(names(default_p$params), names(params_int))
   for (a in miss) {
     params_int[[a]] <- default_p$params[[a]] }

  return (params_int)
}


#'@author andrea Argentini
#'@title Render a DIA-NN report using a Quarto template
#'
#' @param params_report parameters
#' @param template template name
#' @param report_folder description
#' @param report_filename output report file
#' @details
#' The `params` list must contain the following elements:
#' \describe{
#'   \item{\code{report_target_folder}}{The target folder where the report and parameters will be saved.}
#'   \item{\code{input_file}}{The path to the input file.}
#'   \item{\code{design_file}}{The path to the design file.}
#'   \item{\code{title}}{The title of the report.}
#'   \item{\code{subtitle}}{The subtitle of the report.}
#'   \item{\code{author}}{The author of the report.}
#'   \item{\code{description}}{A description of the report.}
#'   \item{\code{contrast}}{The contrast for the analysis.}
#'   \item{\code{aggr_method}}{The aggregation method.}
#'   \item{\code{normalization}}{The normalization method.}
#'   \item{\code{formula}}{The formula for the analysis.}
#'   \item{\code{Proteotypic}}{Whether to use proteotypic peptides.}
#'   \item{\code{pep_per_prot}}{The minimum number of peptides per protein.}
#'   \item{\code{nNonZero}}{The minimum number of non-zero values.}
#'   \item{\code{confounder_list}}{A vector of confounders.}
#'   \item{\code{PCA_comparison}}{A vector of comparisons for PCA plots.}
#'   \item{\code{comparisons}}{A vector of comparisons for differential analysis.}
#'   \item{\code{FC_thr}}{The fold change threshold.}
#'   \item{\code{filtering_contaminant}}{Whether to filter contaminants.}
#'   \item{\code{quantitative_features}}{The quantitative feature to use from DIA-NN output.}
#'   \item{\code{filtPerGroup}}{Whether to filter per group.}
#'   \item{\code{mbr}}{Whether MBR was used in DIA-NN.}
#'   \item{\code{template_file}}{The name of the Quarto template file.}
#'   \item{\code{output_filename}}{The name of the output HTML file.}
#'   \item{\code{DIANN_ver2}}{Whether to read DIANN version > 2 output.}
#' }
#'
#' @return The full path to the rendered report.
#'
#' @importFrom quarto quarto_render
#' @importFrom fs file_move
#' @importFrom logger log_info log_threshold log_appender log_formatter INFO appender_console appender_file
#' @importFrom yaml as.yaml
#' @importFrom utils modifyList
#' @importFrom assertthat assert_that is.string
#' @importFrom matrixStats colMedians
#' @importFrom MsCoreUtils medianPolish  robustSummary

#' @export
render_dia_report <- function(params_report, template, report_folder, report_filename ) {

  # Validate parameters


  validate_template( template)
  validate_folder(report_folder)
  validate_filename( filename = report_filename)

  params_report <- merge_default_parameters(params_report)

  validate_params(params_report)

  # other set up
  log_threshold(logger::INFO)
  log_appender(logger::appender_console)


  log_formatter(logger::formatter_glue)
  # to be done better
  if (params_report$aggr_method  == 'medianPolish'){
    aggr_method_f <- MsCoreUtils::medianPolish
  }
  if (params_report$aggr_method  == 'RobustSummary'){
    aggr_method_f <-  MsCoreUtils::robustSummary
  }
  if (params_report$aggr_method  == 'colMeans'){
    aggr_method_f <- colMeans
  }
  if (params_report$aggr_method  == 'colMeadians'){
    aggr_method_f <- matrixStats::colMedians
  }
  # base protein
  if (template == 'Template_DIA-NN_dev.qmd' ) {
    logfile <- file.path(report_folder, "logfile_protein.log")
    file.create(logfile)  # This will truncate/overwrite the file
    log_appender(logger::appender_file(logfile ), index = 2)
    workflow_steps <- list(
      step_read_diann,
      step_import_qfeat,
      step_add_rowdata_precursor,
      step_filter_na,
      step_processing_protein,
      step_add_rowdata_protein,
      step_add_ensembl,
      step_msqrob_de
    )
    initial_input <- list(params_report = params_report, aggr_method_f = aggr_method_f,layer='proteinRS')

  }

  if (template == 'Template_DIA-NN_dev_EV.qmd' ) {
    logfile <- file.path(report_folder, "logfile_protein.log")
    file.create(logfile)  # This will truncate/overwrite the file
    log_appender(logger::appender_file(logfile ), index = 2)
    workflow_steps <- list(
      step_read_diann,
      step_import_qfeat,
      step_add_rowdata_precursor,
      step_filter_na_ev,
      step_processing_protein_ev,
      step_add_rowdata_protein,
      step_add_ensembl,
      step_msqrob_de,
      step_partial_result
    )
    initial_input <- list(params_report = params_report, aggr_method_f = aggr_method_f,layer='proteinRS')

  }

  # vA peptide
  if (template == 'Template_DIA-NN_peptide_dev_A.qmd'   ){
    logfile <- file.path(report_folder, "logfile_peptide.log")
    file.create(logfile)  # This will truncate/overwrite the file
    log_appender(logger::appender_file(logfile ), index = 2)
    workflow_steps <- list(
      step_read_diann,
      step_import_qfeat,
      step_add_rowdata_precursor,
      step_filter_na,
      step_processing_peptide,
      step_add_rowdata_peptide,
      step_msqrob_de,
      step_partial_result
    )
    params_report$part_item <- ''
    params_report$part_value <- ''
    initial_input <- list(params_report = params_report, aggr_method_f =  base::colSums,layer='peptideNorm')

  }
  # vA protein
   if (template == 'Template_DIA-NN_dev_A.qmd'   ){
     logfile <- file.path(report_folder, "logfile_protein.log")
     file.create(logfile)  # This will truncate/overwrite the file
     log_appender(logger::appender_file(logfile ), index = 2)
     workflow_steps <- list(
       step_read_diann,
       step_import_qfeat,
       step_add_rowdata_precursor,
       step_filter_na,
       step_processing_protein,
       step_add_rowdata_protein,
       step_add_ensembl,
       step_msqrob_de,
       step_partial_result
     )
     params_report$part_item <- ''
     params_report$part_value <- ''
     initial_input <- list(params_report = params_report, aggr_method_f = aggr_method_f,layer='proteinRS')


   }
  # base peptide
  if (template == "Template_DIA-NN_peptide_dev.qmd"){

    logfile <- file.path(report_folder, "logfile_peptide.log")
    file.create(logfile)  # This will truncate/overwrite the file
    log_appender(logger::appender_file(logfile ), index = 2)

    workflow_steps <- list(
      step_read_diann,
      step_import_qfeat,
      step_add_rowdata_precursor,
      step_filter_na,
      step_processing_peptide,
      step_add_rowdata_peptide,
      step_msqrob_de
    )

    initial_input <- list(params_report = params_report, aggr_method_f =  base::colSums,layer='peptideNorm')

    }


  # Prepare aggr_method_f in your input if you want to select it dynamically before running pipeline
 # initial_input <- list(params_report = params_report, aggr_method_f = aggr_method_f)

  result <- run_workflow_pipeline(initial_input, workflow_steps)

  # Find the template file within the package
  template_source_folder <- system.file("quarto_template", package = "diareport")
  if (template_source_folder == "") {
    stop("Template folder not found in the package.")
  }

  # Create a unique temporary working directory
  temp_work_dir <- file.path(tempdir(), paste0("quarto_temp_", Sys.getpid()))
  dir.create(temp_work_dir, recursive = TRUE, showWarnings = FALSE)
  log_info('Temp folder created : {temp_work_dir}')

  # Copy the entire template folder content to the temporary directory
  # This copies all files and subfolders (e.g., resource folders with JS/CSS files)
  success <- file.copy(from = template_source_folder,
                       to = temp_work_dir,
                       recursive = TRUE)
  if (!success) {
    stop("Failed to copy the template folder to the temporary directory.")
  }
  log_info('Copy Template file ...done')
  # Construct the path to the copied template file in the temp directory.
  # Assumes that the template file is directly inside the copied folder.
  temp_template_path <- file.path(temp_work_dir, basename(template_source_folder), template)
  if (!file.exists(temp_template_path)) {
    stop("Template file not found in the temporary directory: ", temp_template_path)
  }
  saveRDS(result$pe, file.path(temp_work_dir,basename(template_source_folder), 'qf_in.RDS'  ))
  saveRDS(result$de_comparison, file.path(temp_work_dir,basename(template_source_folder), 'DEcomp_in.RDS'  ))

  params_report$qf_obj <-   file.path(temp_work_dir,basename(template_source_folder),'qf_in.RDS'  )
  params_report$de_obj <-   file.path(temp_work_dir,basename(template_source_folder),'DEcomp_in.RDS'  )


  if (template == "Template_DIA-NN_dev_A.qmd"  | template == 'Template_DIA-NN_peptide_dev_A.qmd'| template == 'Template_DIA-NN_dev_EV.qmd' ){
    saveRDS(result$part_item, file.path(temp_work_dir,basename(template_source_folder), 'part_item.RDS'  ))
    saveRDS(result$part_value, file.path(temp_work_dir,basename(template_source_folder), 'part_value.RDS'  ))


    params_report$part_item <-   file.path(temp_work_dir,basename(template_source_folder),'part_item.RDS'  )
    params_report$part_value <-   file.path(temp_work_dir,basename(template_source_folder),'part_value.RDS'  )

  }



  path <- file.path(temp_work_dir, basename(template_source_folder))

  tryCatch({
    withr::with_dir(path, {
      quarto_render(
        input = temp_template_path,
        output_format = "html",
        output_file = report_filename,
        execute_params = params_report,

        quarto_args = c("--output-dir", path)
      )
    })
  }, error = function(e) {
    print("Error in Quarto rendering:")
    print(e$message)
    print("Cleaning Temp folder")
    unlink(temp_work_dir, recursive = TRUE)
    stop(e)
  })
  #})

  resource_folder_name <- paste0(tools::file_path_sans_ext(report_filename), "_files")
  rendered_report_path <- file.path(path, report_filename)


  if (!dir.exists(report_folder)) {
    dir.create(report_folder, recursive = TRUE)
  }


  log_info('Copying rendered html report ...')
  # Copy the rendered HTML report to the target folder
  file.copy(from = rendered_report_path, to = file.path(report_folder, report_filename), overwrite = TRUE)
  # If a resource folder was generated, copy it as well
  temp_resource_path <- file.path(temp_work_dir, resource_folder_name)
  if (dir.exists(temp_resource_path)) {
    file.copy(from = temp_resource_path,
              to = file.path(report_folder, resource_folder_name),
              recursive = TRUE, overwrite = TRUE)
  }

  log_info('Cleaning temp folder ...')
  # Optionally, remove the temporary working directory to clean up
  unlink(temp_work_dir, recursive = TRUE)



  log_info("Report generated successfully at: {report_folder}")
  #return(report_folder)
}

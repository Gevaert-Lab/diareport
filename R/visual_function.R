

#' @author andrea argentini
#' @title render_child
#' @description
#'  This function allows to render other template.Rmd in the main quarto document
#' @param data  DE result for each comparison
#' @param path PAth where to store the result
#' @param pe Q feature object, in case this data is required
#' @param label layer name of layer in the qfeature object
#' @param sample_rel list of sample names related to the comparison
#' @param template name of the template .Rms file (contrast,heatmap etc)
#' @return none
#' @export
render_child <- function(data, path, pe, label , sample_rel,  template) {
  if (missing(sample_rel)  ){
    # _templateBArPlot.Rmd _templateContrast _templatePval
    res = knitr::knit_child(
      text = xfun::read_utf8( template),
      envir = rlang::env(data = data, pe = pe,  label = label,  path = path),
      quiet = TRUE
    )
    cat(res, sep = '\n')
    cat("\n")
  } else{
    # _templateHeatmap.Rmd
    res = knitr::knit_child(
      text = xfun::read_utf8( template),
      envir = rlang::env(data = data, pe = pe, label = label , sample_rel = sample_rel, params =params, path = path),
      quiet = TRUE
    )
    cat(res, sep = '\n')
    cat("\n")
  }
}


#' Generate PCA Plots
#'
#' This function generates PCA plots from a list of confounders or variables
#' in the `colData` of a QFeatures object. It supports combinations of
#' one or two variables for color and shape, with specific type combinations:
#' - Character | Factor and Character | Factor
#' - Numerical and Character | Factor
#' - Character | Factor (single variable)
#'
#' The function returns a list of ggplot objects and saves the plots as PDFs.
#'
#' @param pe A QFeatures object containing the data for PCA.
#' @param params A list containing parameters for the run, including `PCA_comparison`
#'   (a character vector of variables to use for PCA) and `folder_prj` (the output folder).
#' @param layer The layer of the QFeatures object to use for PCA.
#' @param test_skip Logical. If `TRUE`, skips saving plots as PDF files.
#' @return A list containing:
#'   - `plots`: A list of ggplot objects for each PCA plot.
#'   - `msgs`: A list of messages indicating the type of PCA performed.
#'   - `status`: A list of statuses (0 for success, 1 for invalid combinations).
#' @importFrom ggplot2 ggplot geom_point xlab ylab ggtitle labs
#' @importFrom dplyr %>%
#' @importFrom SummarizedExperiment assay
#' @importFrom QFeatures filterNA
#' @importFrom stats prcomp
#' @importFrom logger log_info
#' @importFrom grDevices dev.off pdf
#' @importFrom scales percent
#' @export

generate_pca_plots <- function( pe, params, layer,  test_skip= FALSE ) {
  # Define a fixed color palette

  var_topca <- params$PCA_comparison
  output_plot<- list()
  output_msg<- list()
  output_state <- list()
  for (i in seq_along(var_topca)) { # Iterate using index
    v <- var_topca[i] # Access element by index
    log_info(v)

    comparisonPCA <- unlist(strsplit(var_topca[var_topca == v], '-', fixed = TRUE))

    prcompPe <- pe[[layer]] %>%
      filterNA() %>%
      assay() %>%
      t() %>%
      prcomp()

    if (length(comparisonPCA) == 1) {
      # Handle single variable case
      log_info('PCA single variable')
      single_comp <- comparisonPCA[1]
      log_info(single_comp)

      pca_ <- ggplot(data = data.frame(prcompPe$x, SampleName= colData(pe)[['SampleName']],
                                       single_comp = colData(pe)[[single_comp]])  ) +
        ggtitle(paste0("PCA by ", single_comp)) +

        geom_point(aes(x = PC1, y = PC2, colour = factor(single_comp),
                       text = paste("Sample:", SampleName)), size = 3 ) +
        xlab(paste("PC1", percent(summary(prcompPe)$importance[,"PC1"][[2]], accuracy = 0.1))) +
        ylab(paste("PC2", percent(summary(prcompPe)$importance[,"PC2"][[2]], accuracy = 0.1))) +
        labs(colour = single_comp)

      #plot(pca_)
      output_plot[[i]] <- pca_ # Assign plot to list element
      output_msg[[i]] <- 'PCA single variable'
      output_state[[i]] <- 0

      log_info(file.path(params$folder_prj, "Result"))
      if ( test_skip == FALSE){
        pdf(file = file.path(params$folder_prj, "Result", paste0("PCA by ", v, ".pdf")), paper = "a4")
        plot(pca_)

        invisible(dev.off())
      }


    } else if (length(comparisonPCA) == 2) {
      # Handle two variables case
      first_comp <- comparisonPCA[1]
      second_comp <- comparisonPCA[2]
      log_info('2 variables -> Same type ')


      if (class(colData(pe)[[first_comp]]) == class(colData(pe)[[second_comp]])) {
        # Both variables are of the same type
        if (is.character(colData(pe)[[first_comp]]) | is.factor(colData(pe)[[first_comp]])  ) {
          # Both are character variables
          pca_ <- ggplot(data = data.frame(prcompPe$x,SampleName= colData(pe)[['SampleName']],
                                           first_comp = colData(pe)[[first_comp]],
                                           second_comp = colData(pe)[[second_comp]]
          )) +
            ggtitle(paste0("PCA by ", v)) +
            geom_point(aes(x = PC1, y = PC2, colour = factor(first_comp), shape = factor(second_comp),
                           text = paste("Sample:",SampleName)), size = 3) +
            xlab(paste("PC1", percent(summary(prcompPe)$importance[,"PC1"][[2]], accuracy = 0.1))) +
            ylab(paste("PC2", percent(summary(prcompPe)$importance[,"PC2"][[2]], accuracy = 0.1))) +
            labs(colour = first_comp, shape = second_comp)

          output_plot[[i]] <- pca_ # Assign plot to list element
          output_msg[[i]] <- 'PCA two variables are both categorical'
          output_state[[i]] <- 0

          log_info(file.path(params$folder_prj, "Result"))
          if ( test_skip == FALSE){
            pdf(file = file.path(params$folder_prj, "Result", paste0("PCA by ", v, ".pdf")), paper = "a4")
            plot(pca_)

            invisible(dev.off())
          }
        } else {
          # Both are numeric variables
          log_info('I m breaking')
          output_plot[[i]] <- NULL
          output_msg[[i]] <- 'PCA two variable are both numeric NOT ALLOWED '
          output_state[[i]] <- 1
          break
        }
      } else {
        # Variables are of different types
        if (is.numeric(colData(pe)[[first_comp]])) {
          num_comp <- first_comp
          chr_comp <- second_comp
        } else {
          num_comp <- second_comp
          chr_comp <- first_comp
        }
        log_info(paste0("numerical ", num_comp))
        log_info(paste0("categorical ", chr_comp))

        pca_ <- ggplot(data = data.frame(prcompPe$x,SampleName= colData(pe)[['SampleName']],
                                         num_comp = colData(pe)[[num_comp]],
                                         chr_comp = colData(pe)[[chr_comp]]
        )) +
          ggtitle(paste0("PCA by ", v)) +
          geom_point(aes(x = PC1, y = PC2, colour = num_comp,
                         shape = factor(chr_comp),
                         text = paste("Sample:", SampleName )), size = 3 ) +
          xlab(paste("PC1", percent(summary(prcompPe)$importance[,"PC1"][[2]], accuracy = 0.1))) +
          ylab(paste("PC2", percent(summary(prcompPe)$importance[,"PC2"][[2]], accuracy = 0.1))) +
          labs(colour = num_comp, shape = chr_comp)

        #plot(pca_)
        output_plot[[i]] <- pca_ # Assign plot to list element
        output_msg[[i]] <- 'PCA two variables are numeric and categorical '
        output_state[[i]] <- 0

        log_info(file.path(params$folder_prj, "Result"))
        if ( test_skip == FALSE){
          pdf(file = file.path(params$folder_prj, "Result", paste0("PCA by ", v, ".pdf")), paper = "a4")
          plot(pca_)

          invisible(dev.off())
        }
      }
    }
  }
  return ( list(plots = output_plot, msgs=  output_msg,  status = output_state) )
}

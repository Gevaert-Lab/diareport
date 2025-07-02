
#' @author Andrea Argentini
#' @title check_columns_presence
#' @description This function checks the present of a set of features in the columns name
#' of the input data frame.
#' @param df input data frame where features needs to be checked
#' @param min_features data frame containing experiment design data
#' @return status  int 0 / 1 error found
#' @return type_raw format file of the raw file detected ,
#' @return error error message
#' @importFrom logger log_info

check_columns_presence  <- function  ( df, min_features){
  status <- 0
  #type_raw <- NA
  error <- ''


  if  ( ! all( min_features %in% colnames(df)) == TRUE){
    log_info(length(colnames(df) ))
    #cat ( 'Design file not recognized. It should contains at least the following columns:',min_col_need_design,sep='\n\n' )
    #error <-  paste( c('Design file not recognized. It should contains at least the following columns:', paste(min_col_need_design,sep=' ')) ,sep=',' )
    error <- 'Placeholder'
    status <- 1
    return(list(status=status ,error=error))
  }


  return(list(status=status, error= ''))

}

#' @author Andrea Argentini
#' @title  add_ensembl_annotation
#' @description
#' This function adds Ensembl annotation to the assay proteinRS in the q-feature object.
#' Ensembl annotation file is csv file, that must include hgnc_symbol columns. Annotation columns indicated in
#' ensembl_col are included in the rowData(pe_[['proteinRS']]), joining the gene symbol present in DIA-NN report.
#' hgnc_symbol information should be the same of the Genes column present in the DIA-NN, otherwise an error is throw.
#'
#' @param pe_ q_feature object where to add rowdata
#' @param params EDF data frame
#' @return status int 0 non error  / 1 error found
#' @return q_feat q-feature created / modified
#' @return error error message
#' @importFrom dplyr  %>% distinct select left_join join_by
#' @importFrom utils capture.output read.csv
#' @importFrom logger log_info
#' @importFrom SummarizedExperiment rowData
add_ensembl_annotation  <-function (pe_, params ) {
  ## read it
  log_info('Reading Annotation file ...')
  ensembl_db_table <-read.csv(params$ensembl_annotation)
  ## make it simple
  res_check <- check_columns_presence ( ensembl_db_table, min_features = c('hgnc_symbol')  )

  if (res_check$status == 0) {
    log_info('Joining Annotation information into QFeature ...')

    ensembl_db_table <- ensembl_db_table %>%
        distinct(ensembl_gene_id , .keep_all=TRUE)  %>%
        select( params$ensembl_col, hgnc_symbol  ) %>%
        distinct( hgnc_symbol, .keep_all=TRUE)
    temp <- as.data.frame(rowData(pe_[['proteinRS']]))  %>%
            left_join( ensembl_db_table , join_by(Genes ==  hgnc_symbol))

    if (  length(intersect(temp$Genes, ensembl_db_table$hgnc_symbol)) == 0 ){
      return ( list( status= 1 , error= 'Join based hgnc_symbol and Genes did not work. Check your Ids or gene symbol'))

    }


    if (dim(temp)[1] != dim(as.data.frame(rowData(pe_[['proteinRS']])) )[1] ){

      return ( list( status= 1 , error= 'Join based hgnc_symbol and Genes did not work. Check you Ids or gene symbol'))

    }else{
      # adding information specified in  params$ensembl_col
      ## hgnc_symbol should not included
      for (i in  1:(length( params$ensembl_col)) ){
        ann_col <- params$ensembl_col[i]

        rowData(pe_[['proteinRS']])[[ann_col]] <- temp[[ann_col]]
      }}
  }else {
    msg <-  capture.output(cat ( 'Annotation file not recognized. It should contains at least the following column: hgnc_symbol' ))
    return ( list( status= 1 , error= msg))
  }
  return ( list( status= 0 , error= '', q_feat = pe_))
}


#' @author Andrea Argentini
#' @title  filteringNA_qfeat_protein
#' @description
#' This function adds number and percentage of non missing value for each group in assay specified.
#' Moreover it adds computed the number of peptide per protein and the total number of non missing values.
#'
#' @param pe_ q_feature object
#' @param params EDF data frame
#' @param design EDF data frame
#' @return status: int 0 non error  / 1 error found
#' @return q_feat: q-feature created/ modified
#' @return error: error message
#' @importFrom SummarizedExperiment colData
#' @importFrom logger log_info
#' @importFrom  stats as.formula
#' @importFrom QFeatures filterFeatures VariableFilter
#'
filteringNA_qfeat <- function(pe_ , params, design){
  group_val <- design %>% distinct(Group) %>%
    arrange(Group) %>%  pull()
  error <- ''
  status <- 0
  size <- dim(colData(pe_))[1]

  tryCatch( expr = {

    if (params$filtPerGroup != '') {
      log_info('filtering per group')

      ## or
      if (params$filtPerGroup == 'at_least_one') {
        log_info('filtering per group criteria: in at least one group ')
        perc_group_val <- paste0("perc", group_val)
        formula_condition <- paste(paste0(perc_group_val, " >= ", params$nNonZero), collapse = " | ")
        dynamic_formula <- as.formula(paste0("~ ", formula_condition))
      }
      # and
      if (params$filtPerGroup == 'all') {
        log_info('filtering per group criteria: for all the groups ')
        perc_group_val <- paste0("perc", group_val)
        formula_condition <- paste(paste0(perc_group_val, " >= ", params$nNonZero), collapse = " & ")
        dynamic_formula <- as.formula(paste0("~ ", formula_condition))
      }

      pe_ <-  filterFeatures(pe_,dynamic_formula)

    }else{
      log_info('filtering across all samples')
      formula <- as.formula( paste0("~ nNonZero >= ", round(size  * ( params$nNonZero / 100))) )
      pe_ <- filterFeatures(pe_, formula)
      #pe_ <- filterFeatures(pe_, ~ nNonZero >= round(size  * ( params$nNonZero / 100))   )

    }

  },error = function(err){
    print(paste("Q-feature Filtering NaN :  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  } )

  tryCatch( expr = {
    if (params$Proteotypic){
      log_info('Proteotypic filtering')
      formula <- as.formula( paste0("~ Proteotypic == 1") )
      pe_ <- filterFeatures(pe_, formula )
      #pe_ <- filterFeatures(pe_, ~ Proteotypic == 1)

    }

  },error = function(err){
    print(paste("Q-feature Protetypic :  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  } )


  tryCatch( expr = {
    if (params$filtering_contaminant){
      log_info('Filtering contaminants')

      pe_ <- filterFeatures ( pe_ ,VariableFilter("Protein.Ids", params$contaminant_str, "contains", not=TRUE))
    }

  },error = function(err){
    print(paste("Q-feature Contaminant :  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  })

  tryCatch( expr = {

    log_info('Filtering only n peptides per protein')
    formula <- as.formula( paste0("~   pep_per_prot >= ", params$pep_per_prot) )
    pe_ <- filterFeatures(pe_, formula )

  },error = function(err){
    print(paste("Q-feature Num peptide Protein :  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  })


  return( list(error= '', status= 0,q_feat = pe_ ))
}

#' @author Andrea Argentini
#' @title  msqrob_model
#' @description
#' This function perform DE Analysis with the input formula and the comparison of interest
#'
#' @param pe_ q_feature object where to perfomr DE analsysis
#' @param params parameters
#' @param layer name of the layer of the Qfeature obj
#' @return error: modified q-feature object
#' @importFrom SummarizedExperiment rowData assay colData
#' @importFrom logger log_info
#' @importFrom msqrob2 msqrob getCoef makeContrast hypothesisTest
#' @importFrom dplyr left_join select group_by summarise distinct n

msqrob_model <- function(pe_, params, layer ){


  tryCatch( expr = {


    pe_ <- msqrob(object = pe_, i = layer,
                  formula = as.formula( params$formula)  ,ridge = FALSE, overwrite = TRUE)
    contrast_list <- paste0(params$comparisons, "=0")


    coef <- rowData(pe_[[layer]])$msqrobModels[[1]] %>% getCoef %>% names
    log_info('Model fitted ...')
    if (is.null(coef)) {
      coef <- rowData(pe_[[layer]])$msqrobModels[[2]] %>% getCoef %>% names
      #log_info(getCoef(rowData(pe[['proteinRS']])$msqrobModels[[2]]))
      if (is.null(coef)){
        return( list(error= 'Model is not able to converge', status= 1,q_feat =NULL,de_comp= NULL ))
      }
      getCoef(rowData(pe_[[layer]])$msqrobModels[[2]])


    }else{
      #log_info(getCoef(rowData(pe[[layer]])$msqrobModels[[1]]))
      ## print coefficent
      getCoef(rowData(pe_[[layer]])$msqrobModels[[1]])
    }
    log_info('Making contrast & testing')
    L <- makeContrast(contrast_list, parameterNames = coef)
    pe_ <- hypothesisTest(object = pe_, i = layer, contrast = L , overwrite=TRUE)

    res_DE <-  lapply(params$comparisons, dep_volcano, data= pe_,  p=params , layer= layer )

    # params$comparison
    names(res_DE) <-  params$comparison_label

    return (list(error= '', status= 0,q_feat = pe_,de_comp= res_DE ))

  },error = function(err){
    print(paste("Msqrob modeling :  ",err))
    return( list(error= err, status= 1,q_feat =NULL,de_comp= NULL ))
  } )


}

#' @author Andrea Argentini
#' @title  partial_present
#' @description
#' This function perform DE Analysis with the input formula and the comparison of interest
#' This work only with one VARIABLE Design --> ~ Group
#' @param pe_ q_feature object where to perfomr DE analsysis
#' @param params parameters
#' @param layer name of the layer of the Qfeature obj
#' @return error: modified q-feature object
#' @importFrom SummarizedExperiment rowData assay
#' @importFrom logger log_info
#' @importFrom dplyr left_join select  filter bind_rows

partially_present <- function( pe_, params , layer ){
  #  params$comparisons

  tryCatch( expr = {

  if (layer == 'proteinRS'){
    sele_ <- c('Protein.Group','Protein.Ids','Genes')
  }else{
    sele_ <- c( 'Stripped.Sequence', 'Protein.Group','Protein.Ids','Genes')
  }

  params_C <-gsub("Group", "", params$comparisons)
  part_list <- list()
  matrix_list <- list()

  for ( cmp in params_C){
    log_info(paste0('Computing Partial Analysis for ',cmp))
    #print(cmp)
    #cmp <-gsub("Group", "", cmp)
    contrast_values_parsed <- diareport::parse_comparison(cmp,'Group',pe_)
    #print(contrast_values_parsed)
    filt_val <- diareport::select_samples_comparison(contrast_values_parsed, pe_ , 'Group')
    #print(filt_val)
    term <- str_split(cmp, " - ")
    #print(term[[1]])
    log_info(paste0('SelectingPartial Analysis Hits ',cmp))
    percterm <- paste0('perc', unlist(term))
    #print(percterm[[1]])
    #print( percterm[[2]])
    log_info(paste(percterm,collapse = ' '))
    pp_ <-  rowData(pe_[[layer]] )  %>%
       as.data.frame()   %>%
       filter(.data[[percterm[[1]]]] > 1 &  .data[[percterm[[2]]]] == 0 ) %>%
       select ( sele_,.data[[percterm[[1]]]] , .data[[percterm[[2]]]] )

    pp__ <-   rowData(pe_[[layer]] )  %>%
      as.data.frame()   %>%
      filter(.data[[percterm[[1]]]] == 0 &  .data[[percterm[[2]]]] > 1 ) %>%
      select ( sele_,.data[[percterm[[1]]]] , .data[[percterm[[2]]]] )

    final <- dplyr::bind_rows(pp_, pp__)
    log_info(dim(final)[1])
    log_info(dim(final)[2])

    if (layer == 'proteinRS'){
      rownames(final)  <- final$Protein.Ids
    }else{
      rownames(final)  <- final$Stripped.Sequence
    }
    log_info( paste(filt_val,collapse = ' '))
    #print(final %>% head())
    log_info( paste0(rownames(final)[1:10],collapse = ' ' ))

    part_list[[length(part_list) + 1]]  <- final
    ## only the selected Id.
    log_info(paste0('Selecting  Partial Analysis Values ',cmp))

    matrix_list[[length(matrix_list) + 1]] <-  assay(pe_[[layer]])[rownames(final),filt_val]

  }
  names(part_list) <-  params_C
  names(matrix_list) <-  params_C
  return ( list( error= '', status= 0 , part_item =  part_list , part_value = matrix_list ) )
  } ,error = function(err){
    print(paste("Partially Present Analysis :  ",err))
    return( list(error= err, status= 1,part_item =NULL,part_value= NULL ))
  } )
}


#' @author Andrea Argentini
#' @title DEP_volcano
#' @description
#' This function computes volcano plot and return the toptable for each comparison
#' Remark : Model result are supposed to be in proteinRS layer.
#' @param label contrast name
#' @param data Qfeatures object
#' @param params document parameters
#' @param layer QFeature object
#' @return toptable result as dataframe
#' @return volcano volcano ggplot object
#' @return volcano2file volcano ggplot annotated
#' @importFrom SummarizedExperiment rowData
#' @importFrom logger log_info
#' @importFrom ggplot2 ggplot geom_point theme_minimal ylab geom_vline geom_hline aes scale_color_manual ggtitle
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join join_by case_when mutate
#' @importFrom utils head
#' @importFrom ggrepel geom_text_repel
dep_volcano <- function ( label, data  ,params,layer){

  cmp = label

  all_res <-  rowData(data[[layer]])[[label]]
  perc_field <- rowData(data[[layer]]) %>% colnames() %>%  stringr::str_subset('perc')

  if (layer == 'proteinRS') {
    log_info(paste0(cmp,' :Starting fetching Top Table (protein) ...'))

    # protein case
    all_res <- all_res %>% rownames_to_column(var = 'Uniprot_id' )
    if (   all( ! params$ensembl_annotation == '' ))  {
      temp <- as.data.frame(rowData(data[[layer]])) %>% rownames_to_column('Uniprot_id') %>% dplyr::select(Uniprot_id,Genes, Protein.Names, perc_field, head(params$ensembl_col,-1) )

    }else{
      #perc_field <- rowData(data[['proteinRS']]) %>% colnames() %>%  stringr::str_subset('perc')
      temp <- as.data.frame(rowData(data[[layer]])) %>% rownames_to_column('Uniprot_id') %>% dplyr::select(Uniprot_id,Genes, Protein.Names, perc_field )
    }

    all_res <-  all_res %>% left_join( temp, by=join_by(Uniprot_id))

  }else{
    log_info(paste0(cmp,' :Starting fetching Top Table (peptide) ...'))

    ## peptide case
    all_res <- all_res %>% rownames_to_column(var = 'precursor_id' )
    temp <- as.data.frame(rowData(data[['peptideNorm']])) %>% rownames_to_column('precursor_id') %>% select(precursor_id,Genes, Protein.Ids, Protein.Names,perc_field )

    all_res <-  all_res %>% left_join( temp, by=join_by(precursor_id))


  }
  # previous condition -> more complex head(params$ensembl_col,-1) %in%  names(rowData(pe[['proteinRS']])))


  #all_res <-  all_res %>% left_join( temp, by=join_by(Uniprot_id))
  if (params$filt_NaNDE == TRUE) {
  log_info(paste0(cmp,' :#by MSqrob: ', dim(all_res)[1]))
  all_res <- all_res[ ! is.na(all_res$adjPval),]
  log_info(paste0(cmp,' :#by MSqrob after p-adj Null filt.: ', dim(all_res)[1]))
  }else{
    log_info(paste0(cmp,' :#by MSqrob: ', dim(all_res)[1]))
  }

  all_res$differential_expressed <- "NO"
  all_res$differential_expressed[all_res$logFC >= params$FC_thr & all_res$adjPval < params$adjpval_thr] <- "UP"
  all_res$differential_expressed[all_res$logFC <= - params$FC_thr & all_res$adjPval <  params$adjpval_thr] <- "DOWN"

  log_info(paste0(cmp,' Making Volcano plot ...'))

  if ( ! params$ensembl_annotation == '') {
    ## adding ensemble annotation
    p1 <- ggplot(data = all_res , aes(x = logFC, y = -log10(pval) ,col=differential_expressed ,
                                      text = sprintf("Protein_name: %s <br> Gene_symbol: %s  <br> Chromosome name: %s",   all_res$Protein.Names, all_res$Genes,all_res$chromosome_name)   )  )  +
      geom_point() +
      theme_minimal() +
      #geom_text_repel() +
      geom_vline(xintercept = c(- params$FC_thr, params$FC_thr),col="grey") +
      geom_hline(yintercept = -log10(params$adjpval_thr),col="grey") +
      scale_color_manual(values=c("DOWN"="blue","NO"="black", "UP"="red"))+
      ggtitle(paste0("Volcano ",cmp) )

    DEall <- all_res[!is.na(all_res$adjPval) ,append(c('Uniprot_id',  "Protein.Names" , "Genes", "adjPval","pval","logFC","differential_expressed",perc_field),head(params$ensembl_col,-1) ) ]

  }else{
   ##

    p1 <- ggplot(data = all_res , aes(x = logFC, y = -log10(pval) ,col=differential_expressed ,
                                      text = sprintf("Protein_name: %s <br> Gene_symbol: %s", all_res$Protein.Names, all_res$Genes)   )  )  +
      geom_point() +
      theme_minimal() +
      #geom_text_repel() +
      geom_vline(xintercept = c(- params$FC_thr, params$FC_thr),col="grey") +
      geom_hline(yintercept = -log10(params$adjpval_thr),col="grey") +
      scale_color_manual(values=c("DOWN"="blue","NO"="black", "UP"="red"))+
      ggtitle(paste0("Volcano ",cmp) )

    #perc_field <- rowData(data[['proteinRS']]) %>% colnames() %>%  stringr::str_subset('perc')


    if (layer == 'proteinRS') {
      DEall <- all_res[!is.na(all_res$adjPval) , c('Uniprot_id',  "Protein.Names" , "Genes", "adjPval","pval","logFC", "differential_expressed",perc_field)]

    } else {
      DEall <- all_res[!is.na(all_res$adjPval) , c('precursor_id',  "Protein.Names" , "Genes", "adjPval","pval","logFC", "differential_expressed",perc_field)]

    }
  }
  ## volcano annotate with gene name

  if (layer == 'proteinRS') {
    log_info(paste0(cmp,' preparing xport final result (protein)...'))
    all_res_file <- all_res  %>% mutate( label_DE = case_when( differential_expressed == 'UP' ~ Genes,
                                                               differential_expressed == 'DOWN' ~ Genes ,
                                                               TRUE  ~ NA  ))
  } else {
    log_info(paste0(cmp,' preparing xport final result (peptide)...'))

    all_res_file <- all_res %>%  mutate( Gene_v =  case_when( str_detect(Genes, ";") ~ str_split(Genes, ";", simplify = TRUE)[, 1],
                                                             TRUE ~ Genes)) %>%
                          mutate( label_DE = case_when( differential_expressed == 'UP' ~ paste(Gene_v,precursor_id,sep='_'),
                                    differential_expressed == 'DOWN' ~  paste(Gene_v,precursor_id,sep='_') ,
                                    TRUE  ~ NA  ))
  }

  log_info(paste0(cmp,' preparing annotated volcano plot ...'))

  p_toFile <- ggplot(data = all_res_file , aes(x = logFC, y = -log10(pval) ,col=differential_expressed ,
                                               label= label_DE  )  )  +
    geom_point() +
    geom_text_repel() +
    geom_vline(xintercept = c(- params$FC_thr, params$FC_thr),col="grey") +
    geom_hline(yintercept = -log10(params$adjpval_thr),col="grey") +
    scale_color_manual(values=c("DOWN"="blue","NO"="black", "UP"="red"))+
    ggtitle(paste0("Volcano ",cmp) )


  return ( list( toptable =DEall , volcano = p1, volcano2file =p_toFile ) )
}



#' @author Andrea Argentini
#' @title  add_rowdata_detection
#' @description
#' This function adds number and percentage of non missing value for each group in assay specified.
#' Moreover it adds computed the number of peptide per protien and the total number of non missing values.
#'
#' @param pe_ q_feature object where to add rowdata
#' @param design EDF data frame
#' @param assay assay name where to add the statistics
#' @return pe_: modified q-feature object
#' @importFrom SummarizedExperiment rowData assay colData
#' @importFrom logger log_info
#' @importFrom SummarizedExperiment "rowData<-"
#' @importFrom dplyr left_join select group_by summarise distinct n

add_rowdata_detection <- function ( pe_ , design , assay ){
  if (assay == 'precursor' | assay == 'peptideNorm'){
    log_info(paste0('Computing extra info (nNonZero) at ', assay, ' ...'))
    tmp <-  rowData(pe_[[assay]])
    tmp$nNonZero <- rowSums(!is.na(assay(pe_[[assay]])))
    rowData(pe_[[assay]]) <- tmp
    #rowData(pe_[[assay]])$nNonZero <- rowSums(!is.na(assay(pe_[[assay]])))
    log_info(paste0('Computing extra info (pep_per_prot) at ', assay, ' ...'))

    rowData(pe_[[assay]])$pep_per_prot <-
      left_join(rowData(pe_[[assay]]) %>% as.data.frame %>% dplyr::select(Protein.Ids),
                rowData(pe_[[assay]]) %>% as.data.frame %>% dplyr::group_by(Protein.Ids) %>%
                  summarise(pep_per_prot = length(unique(Stripped.Sequence))))$pep_per_prot

    ## add statistics per group

  }
  ## statistics on percetange of detection
  group_val <- design %>% distinct(Group) %>%
    arrange(Group) %>%  pull()
  group_size <- design  %>% group_by(Group) %>% summarise(n_=n()) %>%
    arrange(Group) %>%  pull(n_)
  names(group_size) <- group_val

  log_info(paste0('Computing  % not missing at ', assay, ' ...'))
  for ( v in group_val){
    log_info(v)
    val_count = paste0('Zero',v)
    rowData(pe_[[assay]])[[val_count]] <- assay(pe_[[assay]])[,rownames(colData(pe_)[colData(pe_)$Group == v,]),drop = FALSE]  %>%
      is.na  %>% rowSums()
    ##
    val_perc = paste0('perc',v)
    rowData(pe_[[assay]])[[val_perc]] <- ((group_size[[v]] - rowData(pe_[[assay]])[[val_count]] ) / group_size[[v]]) *  100
  }

  return ( pe_ )
}





#' @author Andrea Argentini
#' @title  processing_qfeat_protein
#'
#' @description
#' This function applies the following steps from precursor assays:
#' 1. Log2 transformation of the percursor intensities (assay name: precursorLog)
#' 2. Normalization of log2transformed intensities based on params$normalization method (assay name: precursorNorm).
#' 3. Protein summarization based on  aggr_method_f method (assay name: proteinRS)
#'
#' @param pe_ q_feature object where to add rowdata
#' @param params EDF data frame
#' @param aggr_mth_fun function for protein summarization
#' @return status int 0 non error  / 1 error found
#' @return q_feat q-feature created / modified
#' @return error error message
#' @importFrom QFeatures logTransform normalize aggregateFeatures
#' @importFrom logger log_info

processing_qfeat_protein <- function(pe_ , params, aggr_mth_fun ){
  error <- ''
  status <- 0
  log_info(paste0('Assays in q-feat object: ', paste(names(pe_), collapse = ", ")) )

  log_info('Intensity log tranformation')
  if ( ! params$normalization == 'vsn'){
  tryCatch( expr = {
    pe_ <- logTransform(pe_, base = 2, i = "precursor",
                       name = "precursorLog")

  },error = function(err){
    print(paste("Q-feature Log-trasformation:  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  } )
  }

  # pe <- logTransform(pe, base = 2, i = "precursor",
  #                    name = "precursorLog")


  log_info('Normalization')

  tryCatch( expr = {

    if  (  params$normalization == 'vsn'){
      pe_ <- normalize(pe_,  method = params$normalization, i = "precursor",
                       name = "precursorNorm")

    }else{
      pe_ <- normalize(pe_,  method = params$normalization, i = "precursorLog",
                       name = "precursorNorm")
    }
    #pe_ <- normalize(pe_,  method = params$normalization, i = "precursorLog",
    #                           name = "precursorNorm")

  },error = function(err){
    print(paste("Q-feature Normalization:  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  } )

  #MsCoreUtils::medianPolish()
  log_info('Summarization at protein level')

  tryCatch( expr = {
    pe_ <- aggregateFeatures(pe_, i = "precursorNorm",
                            fcol = "Protein.Ids",
                            name = "proteinRS",
                            fun = aggr_mth_fun,
                            na.rm = TRUE)

  },error = function(err){
    print(paste("Q-feature Summarization:  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  } )



  log_info(paste0('Assays in q-feat object: ', paste(names(pe_), collapse = ", ")) )

  return( list(error= '', status= 0, q_feat = pe_ ))

}




#' @author Andrea Argentini
#' @title  processing_qfeat_peptide
#'
#' @description
#' This function applies the following steps from precursor assay:
#' 1. Summarize precursor to peptide (stripped sequence) intensities using sum function (assay name PeptideRawSum)
#' 1. Log2 transformation of the peptide intensities (assay name: peptideLog)
#' 2. Normalization of log2transformed peptide intensities based on params$normalization method (assay name: peptideNorm).
#'
#' @param pe_ q_feature object where to add rowdata
#' @param params EDF data frame
#' @param aggr_mth_fun function for protein summarization
#' @return status: int 0 non error  / 1 error found
#' @return q_feat: q-feature created / modified
#' @return error: error message
#' @importFrom QFeatures logTransform normalize aggregateFeatures infIsNA
#' @importFrom logger log_info

processing_qfeat_peptide <- function(pe_ , params, aggr_mth_fun ){
  error <- ''
  status <- 0
  log_info(paste0('Assays in q-feat object: ', paste(names(pe_), collapse = ", ")) )

  log_info('Summarization Precursor -> Peptide (Sum)')


  tryCatch( expr = {
    pe_ <- aggregateFeatures(pe_, i = "precursor",
                             fcol = "Stripped.Sequence",
                             name = "PeptideRawSum",
                             fun = aggr_mth_fun ,
                             # slower but better than medianPolish
                             na.rm = TRUE)

  },error = function(err){
    print(paste("Q-feature Summarization Precursor -> Peptide:  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  } )

  if ( ! params$normalization == 'vsn'){

  log_info('Intensity log tranformation')

  tryCatch( expr = {
    pe_ <- logTransform(pe_, base = 2, i = "PeptideRawSum",
                        name = "peptideLog")
    pe_ <- infIsNA(pe_, i='peptideLog')

  },error = function(err){
    print(paste("Q-feature Log-trasformation:  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  } )
  }

  log_info('Normalization')


  tryCatch( expr = {

    if  (  params$normalization == 'vsn'){
      pe_ <- QFeatures::normalize(pe_,  method = params$normalization, i = "PeptideRawSum",
                                  name = "peptideNorm")
    }else{
      pe_ <- QFeatures::normalize(pe_,  method = params$normalization, i = "peptideLog",
                                  name = "peptideNorm")
    }

  },error = function(err){
    print(paste("Q-feature Normalization:  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  } )

  log_info(paste0('Assays in q-feat object: ', paste(names(pe_), collapse = ", ")) )

  return( list(error= '', status= 0, q_feat = pe_ ))

}







#' @author Andrea
#' @title dfTowideDIANN_20
#' @description
#' This function after some quality check filters DIA-NN result using Q-value and
#' and PG Qvalues.
#' It makes a wide version of the DIA-NN result using the precursorquant information
#' @param data frame containing the DIA-NN report data
#' @param precursorquan columns to use to pivot into a wide format
#' @param mbr  mbr parameter
#' @param wide_colums List of the columns that should be attached to the result
#' @return status : A data frame of DIA-NN result in a wide format
#' @importFrom dplyr  %>% distinct select left_join join_by filter .data
#' @importFrom  tidyr pivot_wider

dfTowideDIANN_20 <- function(data, precursorquan, mbr, wide_colums) {

  if (mbr == FALSE) {

    data %>%
      filter(
        Global.PG.Q.Value <= 0.01 &
          Global.Q.Value <= 0.01 &
          Precursor.Id != "" &
          .data[[precursorquan]] > 0
      ) %>%
      dplyr::select(
        wide_colums,
        .data[[precursorquan]]
      ) %>%
      tidyr::pivot_wider(
        names_from = Run,
        values_from = .data[[precursorquan]]
      )

  }else{   data %>%
      filter(
        Lib.PG.Q.Value <= 0.01 &
          Lib.Q.Value <= 0.01 &
          Precursor.Id != "" &
          .data[[precursorquan]] > 0
      ) %>%
      dplyr::select(
        wide_colums,
        .data[[precursorquan]]
      ) %>%
      tidyr::pivot_wider(
        names_from = Run,
        values_from = .data[[precursorquan]]
      ) }

}


#' @author Andrea Argentini
#' @title  read_DIANN_report
#'
#' @description
#' bla bla
#'
#' @param params q_feature
#' @return status int 0 non error  / 1 error found
#' @return error error message is any
#' @return w_diann diann data in a wide format
#' @return df_design design data as data.frame
#' @importFrom arrow  read_parquet
#' @importFrom utils read.csv2 read.csv
#' @importFrom dplyr %>% mutate select  left_join join_by  rename all_of everything
#' @importFrom QFeatures logTransform normalize aggregateFeatures
#' @importFrom logger log_info
#' @importFrom stringr str_split

read_DIANN_report <- function( params ){
## Describe
  # adding simple function used only inside this part
  detect_file_in_folder <- function(folder_path, pattern) {
    file <- list.files(folder_path, pattern = pattern); if(length(file) > 0) file[1] else "" }

  tryCatch( expr = {

    if (params$DIANN_ver2){
      log_info('Reading DIA-NN PARQUET format...')
      data <- read_parquet(params$input_file)
      lst_wide_columns <- c('Run', 'Precursor.Id', 'Modified.Sequence', 'Stripped.Sequence', 'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'Proteotypic')

    }else{
      log_info('Reading DIA-NN TSV format...')
      data <- read.csv(params$input_file,sep='\t')
      lst_wide_columns <- c('Run', 'Precursor.Id', 'Modified.Sequence', 'Stripped.Sequence', 'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'Proteotypic','First.Protein.Description')

    }
  },error = function(err){
    print(paste("Reading DIANN report :  ",err))
    return( list(error= err, status= 1, w_diann =NULL, df_design = NULL  ))
  } )


   checkDIANN <- check_columns_presence(data , min_features = append(lst_wide_columns, params$quantitative_features))
   if (checkDIANN$status == 1 ) {
     checkDIANN$error <- paste( 'DIA-NN report not recognized. It should contains at least the following columns:',params$quantitative_features,sep='\n')
     return( list(error= checkDIANN$error, status= 1, w_diann =NULL, df_design = NULL  ))

   }

   dfMsqrob <- dfTowideDIANN_20( data, precursorquan = params$quantitative_features,
                                mbr = params$mbr,
                                wide_colums = lst_wide_columns)


   if (params$DIANN_ver2){
     log_info('Retriving Protein description (ver > 2.0)...')

     check_protein_desc <-  detect_file_in_folder(file.path(dirname(params$input_file)),
                                                  'report.protein_description.tsv')
     #log_info(check_protein_desc)
     #log_info( file.path(dirname(params$input_file),check_protein_desc))
     if (check_protein_desc != ""){
       prt_desc <- read.csv( file.path(dirname(params$input_file),check_protein_desc),sep="\t"  )

     }else{
       msg <- paste0('report.protein_description.tsv',  'not found in folder ->',
                     dirname(params$input_file),'\n',' Check your input path !' )
       return( list(error= msg, status= 1, w_diann =NULL, df_design = NULL ))


     }
     # join the  protein description information use Protein.Ids not Protein.Names
     dfMsqrob<- dfMsqrob %>% mutate(app = sapply(str_split(Protein.Ids, ";"), function(x) ifelse(length(x) > 0, x[1], NA))) %>%
       left_join( prt_desc %>% select(Protein.Id,Description) ,
                  join_by(app == Protein.Id)) %>%
       select(- app) %>% rename( First.Protein.Description =  Description   )

     diann_colname <- c("Precursor.Id" , "Modified.Sequence","Stripped.Sequence","Protein.Group",
                        "Protein.Ids","Protein.Names","Genes","Proteotypic","First.Protein.Description")
     #  order the columns in right ways: 1-10 precursor feature columns + 11:  quantitative feature per sample

     dfMsqrob <- dfMsqrob %>% select(all_of(diann_colname), everything())

   }



   tryCatch( expr = {
     L <- readLines(params$design_file, n = 1)
     if (grepl(";", L)) design <- read.csv2(params$design_file) else design <- read.csv(params$design_file)

     ## exit without any error
     return( list(error= '', status= 0, w_diann =dfMsqrob, df_design = design ))
   },error = function(err){
     print(paste("Reading Design file :  ",err))
     return( list(error= err, status= 1, w_diann =dfMsqrob, df_design = NULL ))

   } )




}


#' @author Andrea Argentini
#' @title  import2_qfeature
#' @description
#' This function create an q-feature obj, containing the DIA-NN precursor intensities in the assay precursor.
#' ColData is also set with the experiment design information. Several sanity checks are performed like:
#' 1)  columns present in design file
#' 2)  check number and  sample name found in the DIA-NN report and in the design file
#' 3) check confounder variable present in the design file
#' Finally name matching between run in DIA-report andin the design file is done,
#' using params$wildstr_run to select the run name among the the other information columns
#'
#' @param diaNN_data data frame containing the DIA-NN in a wide format
#' @param design design information as data frame
#' @param params list of parameters
#' @param min_col_need_design list of the mandatory fields in EDF file
#' @param diann_colname columns names
#' @return status int 0 non error  / 1 error found
#' @return q_feat q-feature created/ modified
#' @return error error message
#' @importFrom arrow  read_parquet
#' @importFrom utils read.csv
#' @importFrom dplyr %>% mutate select  left_join join_by  rename arrange
#' @importFrom QFeatures readQFeatures
#' @importFrom logger log_info
#' @importFrom stringr str_detect str_which
#' @importFrom tibble tibble
#' @importFrom SummarizedExperiment 'colData<-'
import2_qfeature <- function (diaNN_data, design, params, min_col_need_design, diann_colname ){

  is_empty <- function(contrast) {
    return(length(contrast) == 0 || (length(contrast) == 1 && contrast == ''))
  }

  log_info('Check EDF file information ...')
  result_check <- check_columns_presence (design, min_features =  min_col_need_design )

  if (result_check$status == 1) {
    msg <- paste('Design file not recognized. It should contains at least the following columns:',min_col_need_design,sep='\n' )
    return( list(error= msg, status= result_check$status, q_feat =NULL ))
  }

  log_info('Check Sample names and their number in DIA-NN and EDF ...')

  checkLength<- check_length_design_data(diaNN_data, design)

  if (checkLength$status==1){

    return( list(error= checkLength$error, status= checkLength$status ,q_feat =NULL ))
  }
  ##  show a more message for this special case
  if (checkLength$status==2){
    diaNN_data <- checkLength$data_
    log_info(checkLength$message)
    #checkLength$message
  }

  ## first check confounders in design file
  log_info('Check confounder in EDF ...')
  if ( ! is_empty(params$confounder_list)  ) {
    check_confounder_list <- check_columns_presence( df = design,  min_features =  params$confounder_list)
    #check_confounder_list<-checkConfounder(confounder= params$confounder_list, colsDesign=colnames(design))
    if (check_confounder_list$status == 1){
      msg <- paste('Confounder not found in the design file -> ',params$confounder_list,sep='\n' )

      return( list(error= msg, status= check_confounder_list$status,q_feat =NULL ))

    }
  }

  log_info('Check comparisons with EDF ...')
  var2check <- 'Group'
  if (! is_empty(params$confounder_list) ) {
    var2check <- append(var2check,params$confounder_list)
  }

  checkVar_res <-checkVariables(inputParams =params$comparisons ,
                                dfDesign = design, variables= var2check)
  if (checkVar_res$status==1){
    #stop(checkVar_res$error)
    return( list(error= checkVar_res$error, status= checkVar_res$status,q_feat = NULL ))
  }

  #  params$wildstr_run
  log_info('Matching DIANN sample name with  EDF ...')

  samplenames <- tibble(
    base_name_sample = names(diaNN_data)[str_which(names(diaNN_data) ,   params$wildstr_run )  ]
  )

  samplenames <- samplenames %>% left_join(  design %>% dplyr::select(Run, Sample) , join_by(base_name_sample ==  Run) )

  names(diaNN_data)[str_which(names(diaNN_data),  params$wildstr_run ) ] <- samplenames$Sample

  if ( !is.null(params$keep_design_order) && (params$keep_design_order == TRUE) ) {
    ## order sample to follow the original design file order
    diaNN_data <- diaNN_data[, c(diann_colname, design$Sample) ]
  }

  tryCatch( expr = {
    log_info('Creating Qfeature object ....')
    pe <- readQFeatures( diaNN_data,
                         fnames = "Precursor.Id",
                         quantCols =  str_detect(names(diaNN_data), paste( diann_colname , collapse = "|"), negate=TRUE) ,
                         name = "precursor")

    if ( !is.null(params$keep_design_order) && (params$keep_design_order == FALSE)) {
      design <- design %>%  arrange(factor(Sample, levels = samplenames$Sample))

    }
    log_info('Setting Qfeature ColData ....')

    # add Coldata
    SummarizedExperiment::colData(pe)$Group <- factor(design$Group)
    SummarizedExperiment::colData(pe)$SampleName <- design$Sample
    #group column from design is now Group in coData(pe)
    SummarizedExperiment::colData(pe)$Replicate <- factor(design$Replicate)
    # setting confounder
    #custom_col <- setdiff(colnames(design),min_col_need_design )
    if (! is_empty(params$confounder_list) ){
        log_info('Importing ConFounder in Q-feature ColData ....')
        for (col_add in params$confounder_list){

          if (is.character(design[[col_add]])) {
            log_info(paste('Type Chr - ', col_add))
            colData(pe)[[col_add]]<- as.factor(design[[col_add]])
            }else{
             log_info(paste('Type Num - ', col_add))
             colData(pe)[[col_add]] <-  as.numeric(design[[col_add]])
            }
          }
    }
    return( list(error= '', status= 0,q_feat =pe ))

  }, error = function(err){
    print(paste("Q-feature err:  ",err))
    return( list(error= err, status= 1,q_feat =NULL ))
  } )
}



#' @author Andrea Argentini
#' @title check_columns_presence
#' @description This function checks some main sanity controls in the design file
#' e.g. Name of columns, min available columns.
#' @param df data frame containing experiment design data
#' @param min_features  vector with the column names to be checked
#' @return status : int 0 / 1 error found
#' @return type_raw: format file of the raw file detected ,
#' @return error: error message

check_columns_presence  <- function  ( df, min_features){
  status <- 0
  #type_raw <- NA
  error <- ''


  if  ( ! all( min_features %in% colnames(df)) == TRUE){
    log_info(length(colnames(df) ))
    #cat ( 'Design file not recognized. It should contains at least the following columns:',min_col_need_design,sep='\n\n' )
    #error <-  paste( c('Design file not recognized. It should contains at least the following columns:', paste(min_col_need_design,sep=' ')) ,sep=',' )
    error <- 'Placeholder'
    status <- 1
    return(list(status=status ,error=error))
  }


  return(list(status=status, error= ''))

}

#' @author Andrea Argentini
#' @title check_length_design_data
#' @description
#' This function checks the names and the number of samples in DIA-NN report and experiment design data,
#' if DIA-NN report has more samples than design file, only sample present in design file are kept.
#'Remark : Model result are supposed to be in proteinRS layer.
#' @param data_ data frame containing the DIA-NN report data
#' @param design data frame containing experiment design data
#' @return status  int 0 / 1: error found, 2: samples in DIA-NN report are more than samples in experiment design data
#' @return error error message
#' @return message: message returned if data frame containing the DIA-NN report data is filtered
#' @importFrom dplyr  %>%  select pull
check_length_design_data  <- function  (data_ , design){
  status <- 0
  error <- ''
  message <- ''

  data_sample <- colnames(data_)[10:length(colnames(data_))]
  ## filename does not exist
  d_sample <- design %>% dplyr::select(Run) %>% pull()

  if (length(data_sample) < length(d_sample)){

    error <- paste0('Number of samples in the design file and in DIA-NN does not match \n  Samples detected in DIANN :\n',  paste(unlist(data_sample), collapse = "\n")  ,  '\n Samples detected in EDF',paste(unlist(d_sample), collapse = "\n")  , '\n')
    status <- 1
    return(list(status=status,error=error,message=message))
  }

  if (length(data_sample[!d_sample %in% data_sample]) >= 1){
    error <- 'Samples in the design file and in DIA-NN do not match'
    status <- 1
    return(list(status=status,error=error,message=message))
  }
  ### pay attention here
  if (length(data_sample) > length(d_sample)){
    status <- 2
    dfSample<- data_[,(colnames(data_)%in% d_sample)]
    ## "First.Protein.Description" on hold for the moment
    df  <- cbind(data_[, c("Precursor.Id" , "Modified.Sequence","Stripped.Sequence","Protein.Group",
                        "Protein.Ids","Protein.Names","Genes","Proteotypic","First.Protein.Description")], dfSample)
    message <- 'Number of samples in DIA-NN is bigger than number of samples in design file.'
    return(list(status=status, error = error, message=message , data_ = df))
  }else{
    return(list(status=status,error=error,message=message))
  }

  #no error exit
  #str_match(data_sample,'\\..*')
  #type_raw <- str_match(data_sample,'\\..*')[,1][1]


}


#' @author Andrea Argentini
#' @title checkVariables
#' @description
#' This function implement sanity check if the values are correct with respect
#' to the design file
#' @param inputParams comparisons groups in input parameter list
#' @param dfDesign data frame containing experiment design data
#' @param variables data frame containing experiment design data
#' @return status  int 0 / 1 error message
#' @return error
#' @importFrom utils capture.output

checkVariables <- function(inputParams, dfDesign, variables) {
  # Initialize lists to collect status and errors
  statusList <- list()
  errorList <- list()

  # 1. Split by '-'
  terms <- unlist(strsplit(inputParams, ' - ', fixed = TRUE))

  # Function to process each term (A and B from the prompt)
  processTerm <- function(term) {
    # Remove parentheses
    term <- gsub("[()]", "", term)
    ## changed to deal with just one variable and not interacted terms
    if (str_detect( term,"[:+]")){
      # 2. Split by ':' and '+'
      parts <- unlist(strsplit(term, "[:+]"))
      values <- trimws(parts)
    }else{
      values <- trimws(term)
    }
    values <- values[values != ""]  # Remove empty strings

    return(values)
  }

  # Process each term
  extractedValues <- unlist(lapply(terms, processTerm))

  extractedValues <- unique(extractedValues)

  for (variable in variables) {

    extracted <- gsub(variable, "", extractedValues[grepl(variable, extractedValues,ignore.case = TRUE)], ignore.case = TRUE)

    # Check if the extracted values are present in the corresponding column in dfDesign
    if (!all(extracted %in% unique(dfDesign[[variable]])) == TRUE) {
      error <- capture.output(cat(paste(variable, ' values in comparison do not match', variable, 'values in design file')))
      status <- 1
    } else {
      error <- ""
      status <- 0
    }

    # Append the status and error to the lists
    statusList[[variable]] <- status
    errorList[[variable]] <- error
  }

  # Aggregate the results
  overallStatus <- ifelse(any(unlist(statusList) == 1), 1, 0)
  overallError <- paste(unlist(errorList), collapse = "\n")

  return(list(status = overallStatus, error = overallError))
}



#' @author Andrea Argentini
#' @title check_dependencies
#' @description
#' This function gets a vector with names of packages and it
#' Check if packages are available, if not it installs them.
#' Otherwise, it loads the  requiredpackages.
#' @param required_packages list with packages to be installed
#' @return none
#' @importFrom utils install.packages
#' @importFrom BiocManager install
#' @export
check_dependencies = function(required_packages = required_packages){
  suppressPackageStartupMessages(
    for(i in required_packages){
      # require returns TRUE invisibly if it was able to load package
      if(! require(i, character.only = TRUE, quietly = TRUE)){
        #  If package was not able to be loaded then re-install
        tryCatch(install.packages(i , dependencies = TRUE), error = function(e) { NULL })
        tryCatch(BiocManager::install(i), error = function(e) { NULL })
        require(i, character.only = TRUE, quietly = TRUE)
      }
    }
  )

}



#' @author Andrea Argentini
#' @title parse_comparison
#' @description
#' This function split the input  comparison label and extract variables and
#' their validity with respect to the ColData() of the current Q-feature object
#' @param input_comparison input comparison label with possible
#' @param variable_names list of the variable name in the formula
#' @param pe  Q-feature object related to the analysis.
#' @return list with left and right parsed value of the comparison
#' @export
parse_comparison <- function(input_comparison, variable_names, pe) {
  # Check if variable_names is empty
  if (length(variable_names) == 0 || (length(variable_names) == 1 && variable_names == '')) {
    stop("The variable_names vector is empty.")
  }

  # Split the string by '-'
  split_strings <- strsplit(input_comparison, " - ")[[1]]

  # Check if the split resulted in two parts
  if (length(split_strings) != 2) {
    stop("Input string does not contain exactly one '-' separator.")
  }

  # Extract substrings A and B
  A <- trimws(split_strings[1]) # Remove leading and trailing whitespaces
  B <- trimws(split_strings[2]) # Remove leading and trailing whitespaces

  # Check if either A or B is empty
  if (A == "" || B == "") {
    stop("One of the left or right substrings is empty.")
  }

  # Function to check if values belong to colData(pe)
  check_values <- function(sub_string, variable_names, pe) {
    # Split the substring by '+'
    terms <- strsplit(sub_string, " \\+ ")[[1]]
    values <- list()

    for (variable in variable_names) {
      values[[variable]] <- c() # Initialize an empty vector for the variable

      for (term in terms) {

        if (length(term) > 0) {
          # Check if the value belongs to colData(pe)
          if (term %in% colData(pe)[[variable]]) {
            values[[variable]] <- c(values[[variable]], term)
          }
        }
      }

      # If no valid value found, set it to NA
      if (length(values[[variable]]) == 0) {
        values[[variable]] <- NA
      }
    }
    values
  }

  # Check values from A and B
  values_A <- check_values(A, variable_names, pe)
  values_B <- check_values(B, variable_names, pe)

  # Return the checked values as a list
  list(Aleft = values_A, Bright = values_B)
}



#' @author Andrea Argentini
#' @title select_samples_comparison
#' @description
#' This function takes the text parsed from the comparison label and select the right.
#' sample in each contrast
#' @param test_parsing input comparison label, text is already parsed
#' @param pe Q-feature object related to the analysis.
#' @param variable_names list of the variable name in the formula
#' @return list with left and right parsed value of the comparison
#' @export

select_samples_comparison <- function(test_parsing, pe, variable_names) {
  sample_sel <- character(0)  # Initialize as character vector

  # Helper function to get samples for a given term (Aleft or Bright)
  get_samples <- function(term, pe, variable_names) {
    # Check if variable_names is valid
    if (is.null(variable_names) || length(variable_names) == 0) {
      return(character(0)) # Return empty character vector if variable_names is missing or empty
    }

    # Create a logical vector to store the combined conditions
    condition <- rep(TRUE, nrow(colData(pe)))

    # Iterate over each variable name and add the condition
    for (variable in variable_names) {
      if (is.null(term[[variable]]) || is.na(term[[variable]])) {
        return(character(0)) # Return empty character vector if variable is missing or NA in the term
      }
      condition <- condition & (colData(pe)[[variable]] == term[[variable]])
    }

    # Select samples that match all conditions
    samples <- rownames(colData(pe)[condition, , drop = FALSE])
    return(samples)
  }

  # Get samples for A left and Bright
  samples_A <- get_samples(test_parsing$Aleft, pe, variable_names)
  samples_B <- get_samples(test_parsing$Bright, pe, variable_names)

  # Combine the samples and return unique values
  sample_sel <- unique(c(samples_A, samples_B))

  return(sample_sel)
}


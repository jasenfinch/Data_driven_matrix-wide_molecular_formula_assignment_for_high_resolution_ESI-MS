
techniques <- c(
  'FIE-HRMS',
  'RP-UHPLC-HRMS'
) 

brachy_targets <- purrr::map(
  techniques,
  ~{
    target_name <- .x %>% 
      stringr::str_replace_all('-','_') %>% 
      paste0('brachy_',.) %>% 
      rlang::sym()
    
    if (stringr::str_detect(.x,'LC')){
      pre_treatment_parameters <- metabolyseR::analysisParameters('pre-treatment')
      metabolyseR::changeParameter(pre_treatment_parameters,'RSDthresh') <- 20
    } else {
      pre_treatment_parameters <- NULL
    }
    
    list( 
      hrmtargets::tar_input_file_path(!!target_name,
                                      mzML_files = metaboData::filePaths(.x,
                                                                         'BdistachyonEcotypes',
                                                                         ask = FALSE),
                                      sample_info = metaboData::runinfo(.x,
                                                                        'BdistachyonEcotypes',
                                                                        ask = FALSE)),
      hrmtargets::tar_spectral_processing(!!target_name,
                                          type = .x %>%
                                            stringr::str_remove('UHP'),
                                          export_path = NULL),
      hrmtargets::tar_pre_treatment(!!target_name,
                                    parameters = pre_treatment_parameters,
                                    export_path = NULL),
      hrmtargets::tar_mf_assignment(
        !!target_name,
        parameters = assignments::assignmentParameters(.x %>%
                                                      stringr::str_remove('UHP')),
        export_path = NULL
      ),
      
      ## Match assignments to B. distachyon specific KEGG compounds
      tar_KEGG_matches(!!target_name,
                       KEGG_IRs_brachy)
    )
  }
)

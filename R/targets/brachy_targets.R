
techniques <- c(
  'FIE-HRMS',
  'RP-UHPLC-HRMS'
) 

assignment_parameters <- list(
  `FIE-HRMS` = assignments::assignmentParameters("FIE-HRMS"),
  `RP-UHPLC-HRMS` = assignments::assignmentParameters("RP-LC-HRMS")
)

assignments::adducts(assignment_parameters$`FIE-HRMS`)$n <- c("[M-H]1-","[M+Cl]1-","[M+Cl37]1-","[M+K-2H]1-","[M-2H]2-","[2M-H]1-")
assignments::adducts(assignment_parameters$`RP-UHPLC-HRMS`)$n <- c("[M-H]1-","[M+Cl]1-","[M+Cl37]1-","[M+K-2H]1-","[M-2H]2-","[2M-H]1-")

assignment_parameters <- assignment_parameters %>% 
  purrr::map(
    ~{
      assignments::isotopes(.x) <- c('13C','18O')
      .x
    }
  )

assignment_parameters <- assignment_parameters %>% 
  purrr::map(
    ~{
      assignments::ppm(.x) <- 4
      .x
    }
  )

brachy_targets <- purrr::map2(
  techniques,
  assignment_parameters,
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
        parameters = .y,
        separate = TRUE,
        export_path = NULL
      ),
      
      ## Match assignments to B. distachyon specific KEGG compounds
      tar_KEGG_matches(!!target_name,
                       KEGG_IRs_brachy)
    )
  }
)

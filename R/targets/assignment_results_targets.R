
assignment_results_targets <- list(
  tar_target(
    assignment_results,
    list(
      `spiked-urine_FIE-HRMS` = spiked_urine_results_molecular_formula_assignment,
      `brachy_FIE-HRMS` = brachy_FIE_HRMS_results_molecular_formula_assignment,
      `brachy_C18-UHPLC-HRMS` = brachy_RP_UHPLC_HRMS_results_molecular_formula_assignment
    ) %>%   map(~tibble(
      n_features = ncol(.x@data),
      n_correlations = nrow(.x@correlations),
      n_relationships = nrow(.x@relationships),
      n_assignments = .x %>% 
        assignments::assignments() %>% 
        nrow(),
      n_MFs = .x %>% 
        assignments::summariseAssignments() %>% 
        nrow()
    )) %>% 
      bind_rows(.id = 'matrix') %>% 
      tidyr::separate(matrix,
                      c('matrix','technique'),sep = '_')
  ),
  
  tar_target(
    assignment_matches,
    list(
      `spiked-urine_FIE-HRMS` = spiked_urine_KEGG_matches,
      `brachy_FIE-HRMS` = brachy_FIE_HRMS_KEGG_matches,
      `brachy_C18-UHPLC-HRMS` = brachy_RP_UHPLC_HRMS_KEGG_matches
    ) %>% 
      map(~tibble(
        KEGG_MF_matches = .x$MF %>% 
          unique() %>% 
          length(),
        KEGG_compound_matches = nrow(.x)
      )) %>% 
      bind_rows(.id = 'matrix') %>% 
      mutate(
        matched_standards = c(spiked_urine_correct_assignments$MF %>% 
                                unique() %>% 
                                length(),
                              rep(NA,2))
      ) %>% 
      tidyr::separate(matrix,
                      c('matrix','technique'),sep = '_')
  )
)
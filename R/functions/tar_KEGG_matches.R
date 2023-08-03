
tar_KEGG_matches <- function(name,db){
  envir <- tar_option_get("envir")
  tidy_eval <- tar_option_get("tidy_eval")
  
  name <- tar_deparse_language(rlang::enexpr(name))
  assignments_name <- paste0(name,'_results_molecular_formula_assignment') %>% 
    rlang::sym()
  IRs_name <- rlang::enexpr(db)
  
  results_name <- paste0(name,'_KEGG_matches')
  
  results_expr <- rlang::expr(!!assignments_name %>% 
                                assignments::assignments() %>% 
                                dplyr::select(MF,Adduct) %>% 
                                dplyr::distinct() %>% 
                                dplyr::rowwise() %>% 
                                dplyr::group_split() %>% 
                                purrr::map_dfr(
                                  ~KEGGmatches(.x$MF,.x$Adduct,!!IRs_name)
                                ) %>% 
                                dplyr::distinct() %>% 
                                dplyr::rowwise() %>% 
                                dplyr::mutate(MF = cheminf::smilesToMF(SMILES)) %>% 
                                dplyr::ungroup()
  )
  
  command_results <- tar_tidy_eval(
    results_expr,
    envir = envir,
    tidy_eval = tidy_eval
  )
  
  target_results <- tar_target_raw(
    results_name,
    command_results
  )
  
  return(list(target_results))
}

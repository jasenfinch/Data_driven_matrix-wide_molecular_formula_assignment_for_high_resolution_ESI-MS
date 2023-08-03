
KEGGmatches <- function(mf,adduct,IRs){
  rule_expr <- mzAnnotation::adduct_rules() %>% 
    dplyr::filter(Name == adduct) %>% 
    .$Rule %>% 
    rlang::parse_expr()
  
  IRs %>%
    cheminf::filterMF(mf) %>% 
    cheminf::filterIP(rule = !!rule_expr) %>% 
    cheminf::entries() %>% 
    mutate(Adduct = adduct)
}

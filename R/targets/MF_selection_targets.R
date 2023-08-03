
KEGG_alternative_MFs_ppm <- 6

MF_selection_targets <- list(
  
  ## Calculate golden rules for KEGG MFs
  KEGG_compounds_valid_MFs_golden_rules = tar_target(
    KEGG_compounds_valid_MFs_golden_rules,
    KEGG_compounds_valid_MFs$MF %>% 
      mzAnnotation::goldenRules()
  ),
  
  ## Calculate KEGG MF plausibilities
  KEGG_compound_plausibilities = tar_target(
    KEGG_compound_plausibilities,
    KEGG_compounds_valid_MFs_golden_rules %>% 
      mzAnnotation::goldenRulesScore()
  ),
  
  ## Generate alternative KEGG MF plausibilities
  KEGG_alternative_MFs = tar_target(
    KEGG_alternative_MFs,
    KEGG_compounds_valid_MFs %>% 
      dplyr::rowwise() %>% 
      dplyr::group_split() %>% 
      furrr::future_map(~ mzAnnotation::ipMF(.x$M,
                                             adduct = 'M',
                                             ppm = KEGG_alternative_MFs_ppm) %>% 
                          dplyr::mutate(Rank = rank(100 - `Plausibility (%)`,
                                                    ties.method = 'min'),
                                        `KEGG MF` = .x$MF),
                        .options = furrr::furrr_options(seed = TRUE))
  ),
  
  ## Retrieve the MF plausibility ranks of the KEGG compounds
  KEGG_compound_MF_ranks = tar_target(
    KEGG_compound_MF_ranks,
    KEGG_alternative_MFs %>% 
      furrr::future_map(~{
        .x %>% 
          dplyr::mutate(N = nrow(.)) %>% 
          dplyr::filter(MF == `KEGG MF`) %>% 
          dplyr::select(MF,`Measured M`,`Plausibility (%)`,Rank,N)
      }) %>% 
      dplyr::bind_rows()
  ),
  
  ## Plot the KEGG compound MF ranks against the exact mass
  KEGG_compound_MF_ranks_mass_plot = tar_target(
    KEGG_compound_MF_ranks_mass_plot,
    KEGG_compound_MF_ranks %>% 
      ggplot2::ggplot(ggplot2::aes(x = Rank,
                                   y = `Measured M`)) +
      ggplot2::geom_point(alpha = 0.3,size = 1) +
      jfmisc::theme_neat() +
      ggplot2::labs(x = '*P<sub>MF</sub>* rank',
                    y = 'Monoisotopic mass') +
      theme(axis.title.x = ggtext::element_markdown())
  ),
  
  ## Plot the KEGG compound MF ranks against the number of generated MFs
  KEGG_compound_MF_ranks_n_plot = tar_target(
    KEGG_compound_MF_ranks_n_plot,
    KEGG_compound_MF_ranks %>% 
      ggplot2::ggplot(ggplot2::aes(x = Rank,
                                   y = N)) +
      ggplot2::geom_point(alpha = 0.3,size = 1) +
      jfmisc::theme_neat() +
      ggplot2::labs(x = '*P<sub>MF</sub>* rank',
                    y = 'Number of generated molecular formulas') +
      theme(axis.title.x = ggtext::element_markdown(face = 'bold'))
  )
)


search_ppm <- 4

assignment_parameters <- {
  p <- assignments::assignmentParameters("FIE-HRMS") 
  assignments::adducts(p)$n <- c("[M-H]1-","[M+Cl]1-","[M+Cl37]1-","[M+K-2H]1-","[M-2H]2-","[2M-H]1-")
  assignments::adducts(p)$p <- c("[M+Na]1+","[M+H]1+","[M+K]1+","[M+K41]1+","[M+2H]2+","[2M+H]1+")
  assignments::isotopes(p) <- c("13C","18O","34S")
  assignments::ppm(p) <- 4
  p
}

spiked_urine_example_feature <- 'n80.97488' 


matrix_info <- tibble::tibble(
  sample = c(
    'standards',
    'spiked_urine'
  ),
  tag = c(
    'standards-mix',
    'spiked-urine'
  ),
  sample_info = c(
    c(
      'sample_information.csv',
      'spiked_urine_runinfo.csv'
    )
  ),
  QCidx = c(
    'D3_mix',
    'QC'
  )
)

assignment_outcomes <- list(
  standards = c(
    'No correlations',
    'No relevant relationships',
    'No relevant relationships',
    'Eliminated, unassigned',
    'Alternative adduct, alternative MF',
    'Matching adduct, MF outside top 3',
    'Matching adduct, matching MF'
  ),
  spiked_urine = c(
    'No correlations',
    'No relevant relationships',
    'No relevant relationships',
    'Eliminated, unassigned',
    'Alternative adduct, alternative MF',
    'Matching adduct, MF outside top 3',
    'Matching adduct, alternative MF',
    'Matching adduct, matching MF'
  ) 
)

standards_processing_targets <- matrix_info %>% 
  dplyr::rowwise() %>% 
  dplyr::group_split() %>% 
  purrr::map(
    ~{
      list(
        ## download mzML files
        hrmtargets::tar_input_piggyback(
          !!rlang::sym(.x$sample),
          tag = !!.x$tag,
          repo = 'jasenfinch/Data_driven_matrix-wide_molecular_formula_assignment_for_high_resolution_ESI-MS-devel',
          sample_info_file = .x$sample_info),
        
        ## Perform spectral binning
        hrmtargets::tar_spectral_processing(
          !!rlang::sym(.x$sample),
          type = 'FIE-HRMS',
          export_path = NULL),
        
        ## pre-treat spiked_urine features
        hrmtargets::tar_pre_treatment(
          !!rlang::sym(.x$sample),
          parameters = metabolyseR::analysisParameters(
            'pre-treatment') %>% 
            {
              metabolyseR::parameters(.,
                                      'pre-treatment') <- metabolyseR::preTreatmentParameters(
                                        list(
                                          occupancyFilter = 'maximum',
                                          impute = 'all',
                                          QC = 'RSDfilter',
                                          transform = 'TICnorm'
                                        )) 
              metabolyseR::changeParameter(.,'QCidx') <- .x$QCidx
              .
            },
          plots = 'unsupervised_RF',
          export_path = NULL),
        
        ## Assign molecular formulas
        hrmtargets::tar_mf_assignment(
          !!rlang::sym(.x$sample),
          parameters = assignment_parameters,
          separate = TRUE,
          export_path = NULL,
          resources = tar_resources(
            future = tar_resources_future(
              plan = future::tweak(
                jfmisc::suitableParallelStrategy(),
                workers = jfmisc::suitableParallelWorkers(
                  RAM_per_worker = 8
                )
              )
            )
          )
        )
      )   
    }
  )

standards_compound_targets <- list(
  ## chemical standards compound info file
  tarchetypes::tar_file(
    standards_compound_info_file,
    {
      'spiked_chemical_standards_info.csv' %>%
        piggyback::pb_download(dest = 'data',
                               tag = 'spiked-urine')

      'data/spiked_chemical_standards_info.csv'
    }
  ),

  ## load chemical standards compound info
  tar_target(
    standards_compound_info,
    readr::read_csv(standards_compound_info_file) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(SMILES = cheminf::convert(
        InChI,
        'INCHI',
        'SMILES'),
        MF = cheminf::smilesToMF(SMILES))
  ),

  ## calculate chemical standards compound descriptors
  tar_target(
    standards_compound_db,
    cheminf::metaboliteDB(
      standards_compound_info %>%
        dplyr::rename(ID = CID,
                      NAME = Name)
    )
  ),

  ## calculate possible adducts for the standards
  tar_target(
    standards_compound_adducts,
    standards_compound_db %>%
      cheminf::entries() %>%
      .$ID %>%
      furrr::future_map(cheminf::calcAdducts,
                        db = standards_compound_db) %>%
      purrr::set_names(standards_compound_db %>%
                         cheminf::entries() %>%
                         .$ID) %>%
      dplyr::bind_rows(.id = 'ID') %>%
      dplyr::left_join(mzAnnotation::adduct_rules() %>%
                         dplyr::select(Name,Nelec),
                       by = c('Adduct' = 'Name')) %>%
      dplyr::mutate(ID = as.numeric(ID),
                    Mode = replace(Nelec,
                                   Nelec > 0,
                                   'n') %>%
                      replace(Nelec < 0,
                              'p')) %>%
      dplyr::select(-Nelec,-MF) %>%
      dplyr::filter(Possible == TRUE,
                    Adduct %in%
                      mzAnnotation::adduct_rules()$Name[mzAnnotation::adduct_rules()$Default == 1]) %>%
      dplyr::left_join(standards_compound_db %>%
                         cheminf::entries(),
                       by = 'ID') %>%
      dplyr::select(ID,NAME,InChI,MF,Adduct,Mode,`m/z`) %>%
      dplyr::filter(
        Adduct %in% {assignment_parameters %>%
            assignments::adducts() %>% 
            purrr::flatten_chr()})
  ),

  ## calculate possible spiked compound isotopes
  tar_target(
    standards_compound_isotopes,
    assignment_parameters %>%
      assignments::isotopes() %>%
      purrr::map_dfr(~{
        isotope <- .x
        standards_compound_adducts %>%
          dplyr::group_split(MF) %>%
          furrr::future_map(~{
            if (mzAnnotation::isotopePossible(.x$MF[1],isotope = isotope)){
              .x %>%
                dplyr::rowwise() %>%
                dplyr::mutate(Isotope = isotope,
                              `m/z` = MF %>%
                                mzAnnotation::calcAccurateMass() %>%
                                mzAnnotation::calcMZ(adduct = Adduct,
                                                     isotope = isotope))
            }
          },.options = furrr::furrr_options(seed = TRUE))
      })
  ),

  ## combine all possible PIPs for the chemical standards
  tar_target(
    standards_compound_PIPs,
    dplyr::bind_rows(
      standards_compound_adducts,
      standards_compound_isotopes
    )
  )
)

standards_evaluation_targets <- matrix_info %>%
  dplyr::rowwise() %>%
  dplyr::group_split() %>%
  purrr::map(
    ~{
      list(
        tar_target_raw(
          paste0(.x$sample,'_features'),
          rlang::expr(
            !!rlang::sym(paste0(.x$sample,'_results_molecular_formula_assignment')) %>%
              assignments::featureData() %>%
              {
                tibble::tibble(
                  Feature = colnames(.)
                )
              } %>%
              dplyr::mutate(
                Mode = stringr::str_sub(Feature,1,1),
                mz = stringr::str_remove_all(Feature,'[:alpha:]') %>%
                  as.numeric()
              )
          )
        ),

        tar_target_raw(
          paste0(.x$sample,'_compound_matches'),
          rlang::expr(
            standards_compound_PIPs %>%
              dplyr::rowwise() %>%
              dplyr::group_split() %>%
              purrr::map_dfr(~{
                ppm_range <- mzAnnotation::ppmRange(.x$`m/z`,
                                                    assignments::ppm(assignment_parameters))

                matches <- !!rlang::sym(paste0(.x$sample,'_features')) %>%
                  dplyr::filter(Mode == .x$Mode,
                                mz > ppm_range$lower,
                                mz < ppm_range$upper) %>%
                  dplyr::mutate(`PPM error` = mzAnnotation::ppmError(mz,
                                                                     .x$`m/z`) %>%
                                  abs())

                if(nrow(matches) > 0){
                  dplyr::bind_cols(dplyr::select(.x,
                                                 ID:Adduct,Isotope,
                                                 `Theoretical m/z` = `m/z`),
                                   matches)
                } else {
                  NULL
                }
              }) %>%
              dplyr::bind_rows()
          )
        ),

        tar_target_raw(
          paste0(.x$sample,'_feature_matches'),
          rlang::expr(
            !!rlang::sym(paste0(.x$sample,'_compound_matches')) %>%
              dplyr::select(
                Feature,
                Mode,mz,
                MF,
                Isotope,
                Adduct,
                `PPM error`) %>%
              dplyr::distinct() %>%
              dplyr::arrange(MF,Adduct,Isotope) %>%
              dplyr::group_split(MF) %>%
              purrr::map_dfr(~{
                if (any(is.na(.x$Isotope))){
                  return(.x)
                } else {
                  return(NULL)
                }
              })
          )
        ),

        tar_target_raw(
          paste0(.x$sample,'_correlations'),
          rlang::expr(
            !!rlang::sym(paste0(.x$sample,'_results_molecular_formula_assignment')) %>%
              assignments::correlations() %>%
              dplyr::select(
                contains('Feature')
              ) %>%
              {
                dplyr::bind_rows(
                  dplyr::select(.,Feature = Feature1),
                  dplyr::select(.,Feature = Feature2),
                ) %>%
                  dplyr::distinct()
              } %>%
              dplyr::mutate(
                correlations = TRUE
              )
          )
        ),

        tar_target_raw(
          paste0(.x$sample,'_relationships'),
          rlang::expr(
            !!rlang::sym(paste0(.x$sample,'_results_molecular_formula_assignment')) %>%
              assignments::relationships() %>%
              dplyr::filter(
                dplyr::if_all(
                  dplyr::contains('Transformation'),
                  ~is.na(.x)
                ),
                coefficient > 0
              ) %>%
              dplyr::select(
                contains(
                  c(
                    'Feature',
                    'Adduct',
                    'Isotope')
                )
              ) %>%
              {
                dplyr::bind_rows(
                  dplyr::select(
                    .,
                    Feature = Feature1,
                    Adduct = Adduct1,
                    Isotope = Isotope1
                  ),
                  dplyr::select(
                    .,
                    Feature = Feature2,
                    Adduct = Adduct2,
                    Isotope = Isotope2
                  ),
                ) %>%
                  dplyr::distinct()
              } %>%
              dplyr::mutate(
                relationships = TRUE
              )
          )
        ),

        tar_target_raw(
          paste0(.x$sample,'_assignments'),
          rlang::expr(
            !!rlang::sym(paste0(.x$sample,'_results_molecular_formula_assignment')) %>%
              assignments::assignments() 
          )
        ),
        
        tar_target_raw(
          paste0(.x$sample,'_correct_assignments'),
          rlang::expr(
            !!rlang::sym(paste0(.x$sample,'_assignments')) %>%
              inner_join(!!rlang::sym(paste0(.x$sample,'_feature_matches')) %>% 
                           select(Feature,MF,Adduct,Isotope), 
                         by = c("Feature", 
                                "Isotope", 
                                "Adduct", 
                                "MF"))
          )
        ),
        
        tar_target_raw(
          paste0(.x$sample,'_assignments_AI'),
          rlang::expr(
            !!rlang::sym(paste0(.x$sample,'_assignments')) %>%
              dplyr::filter(
                !stringr::str_detect(
                  Iteration,
                  'T'
                )
              ) %>%
              dplyr::select(
                Feature,
                Adduct,
                Isotope,
                MF
              )
          )
        ),

        tar_target_raw(
          paste0(.x$sample,'_assignment_outcomes'),
          rlang::expr(
            {
              assignment_outcomes <- !!rlang::sym(paste0(.x$sample,'_feature_matches')) %>%
                dplyr::select(
                  MF,
                  Adduct,
                  Isotope,
                  Feature
                ) %>%
                dplyr::distinct() %>%
                dplyr::left_join(
                  !!rlang::sym(paste0(.x$sample,'_correlations')),
                  by = 'Feature'
                ) %>%
                dplyr::left_join(
                  !!rlang::sym(paste0(.x$sample,'_relationships')),
                  by = c(
                    'Feature',
                    'Adduct',
                    'Isotope'
                  )
                ) %>%
                dplyr::left_join(
                  !!rlang::sym(paste0(.x$sample,'_assignments_AI')) %>%
                    dplyr::select(
                      Feature
                    ) %>%
                    dplyr::mutate(
                      assigned = TRUE
                    ),
                  by = c(
                    'Feature'
                  )
                ) %>%
                dplyr::left_join(
                  !!rlang::sym(paste0(.x$sample,'_assignments_AI')) %>%
                    dplyr::select(
                      Feature,
                      Isotope,
                      Adduct
                    ) %>%
                    dplyr::mutate(
                      matching_adduct = TRUE
                    ),
                  by = c(
                    'Feature',
                    'Isotope',
                    'Adduct'
                  )
                ) %>%
                dplyr::left_join(
                  !!rlang::sym(paste0(.x$sample,'_assignments_AI')) %>%
                    dplyr::mutate(
                      matching_assignment = TRUE
                    ),
                  by = c(
                    'Feature',
                    'Adduct',
                    'Isotope',
                    'MF'
                  )
                )

              mf_top_ranked <- assignment_outcomes %>%
                dplyr::filter(
                  dplyr::if_all(
                    correlations:matching_adduct,
                    ~.x == TRUE)
                ) %>%
                dplyr::mutate(
                  `m/z` = Feature %>%
                    stringr::str_remove_all(
                      '[:alpha:]'
                    ) %>%
                    as.numeric()
                ) %>%
                dplyr::rowwise() %>%
                dplyr::group_split() %>%
                furrr::future_map_dfr(
                  ~{
                    mfs <- mzAnnotation::ipMF(
                      mz = .x$`m/z`,
                      adduct = .x$Adduct,
                      isotope = .x$Isotope,
                      ppm = assignments::ppm(assignment_parameters)
                    ) %>%
                      dplyr::slice(1:3)

                    .x %>%
                      dplyr::mutate(
                        mf_top_3 = MF %in% mfs$MF
                      )
                  },
                  seed = TRUE
                ) %>%
                dplyr::select(
                  MF:Feature,
                  mf_top_3
                )

              assignment_outcomes %>%
                dplyr::left_join(
                  mf_top_ranked,
                  by = c(
                    'MF',
                    'Adduct',
                    'Isotope',
                    'Feature'
                  )
                ) %>%
                dplyr::relocate(
                  mf_top_3,
                  .before = matching_assignment
                ) %>%
                dplyr::mutate(
                  dplyr::across(
                    correlations:matching_assignment,
                    ~replace(
                      .x,
                      is.na(.x),
                      FALSE
                    )
                  )
                )
            }
          )
        ),

        tar_target_raw(
          paste0(.x$sample,'_assignment_outcomes_non_iso'),
          rlang::expr(
            !!rlang::sym(paste0(.x$sample,'_assignment_outcomes')) %>%
              dplyr::filter(
                is.na(Isotope),
                !Adduct %in% c(
                  '[M+Cl37]1-',
                  '[M+K41]1+'
                )
              )
          )
        ),

        tar_target_raw(
          paste0(.x$sample,'_assignment_outcomes_summary'),
          rlang::expr(
            !!rlang::sym(paste0(.x$sample,'_assignment_outcomes_non_iso')) %>%
              dplyr::count(
                correlations,
                relationships,
                assigned,
                matching_adduct,
                mf_top_3,
                matching_assignment
              ) %>%
              dplyr::mutate(
                `Assignment outcome` = factor(
                  !!assignment_outcomes[[.x$sample]],
                  levels = c(
                    'No correlations',
                    'No relevant relationships',
                    'Eliminated, unassigned',
                    'Alternative adduct, alternative MF',
                    'Matching adduct, MF outside top 3',
                    'Matching adduct, alternative MF',
                    'Matching adduct, matching MF'
                  )
                )
              ) %>%
              dplyr::group_by(
                `Assignment outcome`
              ) %>%
              dplyr::summarise(
                n = sum(n)
              ) %>%
              dplyr::mutate(
                `%` = n / sum(n) * 100
              ) %>%
              dplyr::rename(
                `# IPs` = n
              )
          )
        )
      )
    }
  )

standards_figures_targets <- list(
  ## Plot component solutions for example feature
  tar_target(
    spiked_urine_feature_solutions_plot,
    {
      pl <- plotFeatureSolutions(
        spiked_urine_results_molecular_formula_assignment,
        tibble::tibble(
          component = c(43,423,1316),
          border = c('red',rep('black',2))
        ),
        spiked_urine_example_feature,
        'A&I1'
      ) &
        guides(
          edge_colour = ggraph::guide_edge_colorbar(
            title = 'correlation coefficient',
            title.position = 'top',
            title.hjust = 0.5,
            barheight = 0.75,
            barwidth = 8,
            frame.colour = 'black',
            direction = 'horizontal'
          )
        )

      legend <- ggpubr::get_legend(pl)

      {pl &
          theme(legend.position = 'none')} /
        legend +
        patchwork::plot_layout(heights = c(6,1))
    }
  ),
  
  ## Match assignments to human specific KEGG compounds
  tar_KEGG_matches(spiked_urine,
                   KEGG_IRs_human),
  
  tar_target(
    assignment_outcomes_summary,
    list(
      standards = standards_assignment_outcomes_summary,
      spiked_urine = spiked_urine_assignment_outcomes_summary
    ) %>% 
      dplyr::bind_rows(
        .id = 'sample'
      ) %>% 
      dplyr::select(-`%`) %>% 
      tidyr::spread(
        sample,`# IPs`,
        fill = 0
      )
  )
)

standards_targets <- list(
  standards_processing_targets,
  standards_compound_targets,
  standards_evaluation_targets,
  standards_figures_targets
)
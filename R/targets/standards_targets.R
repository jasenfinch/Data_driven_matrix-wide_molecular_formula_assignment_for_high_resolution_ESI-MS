
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

standards_targets <- list(
  standards_processing_targets
)

# standards_evaluation_targets <- list(
#   
# )
# 
# standards_targets <- list(
#   ## spiked_urine compound info file
#   tarchetypes::tar_file(
#     spiked_urine_compound_info_file,
#     {
#       'spiked_chemical_standards_info.csv' %>% 
#         piggyback::pb_download(dest = 'data',
#                                tag = 'spiked-urine')
#       
#       'data/spiked_chemical_standards_info.csv' 
#     }
#   ),
#   
#   ## load spiked_urine compound info
#   tar_target(
#     spiked_urine_compound_info,
#     readr::read_csv(spiked_urine_compound_info_file) %>% 
#       dplyr::rowwise() %>% 
#       dplyr::group_split() %>% 
#       purrr::map_dfr(~ .x %>% 
#                        dplyr::mutate(SMILES = cheminf::convert(
#                          InChI,
#                          'INCHI',
#                          'SMILES'),
#                          MF = cheminf::smilesToMF(SMILES)))
#   ),
#   
#   ## calculate spiked compound descriptors
#   tar_target(
#     spiked_urine_compound_db,
#     cheminf::metaboliteDB(
#       spiked_urine_compound_info %>% 
#         dplyr::rename(ID = CID,
#                       NAME = Name)
#     )
#   ),
#   
#   ## calculate possible spiked compound adducts
#   tar_target(
#     spiked_urine_compound_adducts,
#     spiked_urine_compound_db %>%
#       cheminf::entries() %>%
#       .$ID %>%
#       purrr::map(cheminf::calcAdducts,
#                  db = spiked_urine_compound_db) %>%
#       set_names(spiked_urine_compound_db %>%
#                   cheminf::entries() %>%
#                   .$ID) %>%
#       dplyr::bind_rows(.id = 'ID') %>%
#       dplyr::left_join(mzAnnotation::adduct_rules() %>% 
#                          dplyr::select(Name,Nelec),
#                        by = c('Adduct' = 'Name')) %>%
#       dplyr::mutate(ID = as.numeric(ID),
#                     Mode = replace(Nelec,
#                                    Nelec > 0,
#                                    'n') %>% 
#                       replace(Nelec < 0,
#                               'p')) %>%
#       dplyr::select(-Nelec,-MF) %>%
#       dplyr::filter(Possible == TRUE,
#                     Adduct %in% 
#                       mzAnnotation::adduct_rules()$Name[mzAnnotation::adduct_rules()$Default == 1]) %>%
#       dplyr::left_join(spiked_urine_compound_db %>%
#                          cheminf::entries(),
#                        by = 'ID') %>%
#       dplyr::select(ID,NAME,InChI,MF,Adduct,Mode,`m/z`) %>%
#       dplyr::filter(
#         Adduct %in% {assignments::adducts(spiked_urine_results_molecular_formula_assignment) %>% 
#             purrr::flatten_chr()})
#   ),
#   
#   ## calculate possible spiked compound isotopes
#   tar_target(
#     spiked_urine_compound_isotopes,
#     spiked_urine_results_molecular_formula_assignment %>%
#       assignments::isotopes() %>% 
#       map_dfr(~{
#         isotope <- .x
#         spiked_urine_compound_adducts %>% 
#           group_split(MF) %>% 
#           furrr::future_map(~{
#             if (mzAnnotation::isotopePossible(.x$MF[1],isotope = isotope)){
#               .x %>% 
#                 rowwise() %>% 
#                 mutate(Isotope = isotope,
#                        `m/z` = MF %>% 
#                          mzAnnotation::calcAccurateMass() %>% 
#                          mzAnnotation::calcMZ(adduct = Adduct,
#                                               isotope = isotope))
#             }
#           },.options = furrr::furrr_options(seed = TRUE))    
#       })
#   ),
#   
#   ## combine possible spiked compound PIP
#   tar_target(
#     spiked_urine_compound_PIPs,
#     bind_rows(
#       spiked_urine_compound_adducts,
#       spiked_urine_compound_isotopes
#     )
#   ),
#   
#   tar_target(
#     spiked_urine_assignments,
#     spiked_urine_results_molecular_formula_assignment %>% 
#       assignments::assignments()
#   ),
#   
#   ## Match assignments to human specific KEGG compounds
#   tar_KEGG_matches(spiked_urine,
#                    KEGG_IRs_human),
#   
#   ## Extract the spiked urine features
#   tar_target(
#     spiked_urine_features,
#     spiked_urine_results_pre_treatment %>% 
#       {
#         tibble::tibble(
#           Feature = metabolyseR::features(.,
#                                           type = 'pre-treated')
#         ) 
#       } %>% 
#       dplyr::mutate(
#         Mode = stringr::str_sub(Feature,1,1),
#         mz = stringr::str_remove_all(Feature,'[:alpha:]') %>%
#           as.numeric()
#       )
#   ),
#   
#   ## Directly match spiked compound putative adducts to features
#   tar_target(
#     spiked_urine_compound_matches,
#     spiked_urine_compound_PIPs %>%
#       dplyr::rowwise() %>% 
#       dplyr::group_split() %>% 
#       purrr::map_dfr(~{
#         ppm_range <- mzAnnotation::ppmRange(.x$`m/z`,
#                                             spiked_urine_search_ppm)
#         
#         matches <- spiked_urine_features %>%
#           dplyr::filter(Mode == .x$Mode,
#                         mz > ppm_range$lower,
#                         mz < ppm_range$upper) %>%
#           dplyr::mutate(`PPM error` = mzAnnotation::ppmError(mz,
#                                                              .x$`m/z`) %>%
#                           abs())
#         
#         if(nrow(matches) > 0){
#           dplyr::bind_cols(dplyr::select(.x,
#                                          ID:Adduct,Isotope,
#                                          `Theoretical m/z` = `m/z`),
#                            matches)
#         } else {
#           NULL
#         }
#       }) %>%
#       dplyr::bind_rows()
#   ),
#   
#   ## Reduce compound matches to MF feature matches
#   tar_target(
#     spiked_urine_feature_matches,
#     spiked_urine_compound_matches %>% 
#       select(Feature,Mode,mz,MF,Isotope,Adduct,`PPM error`) %>% 
#       distinct() %>% 
#       arrange(MF,Adduct,Isotope) %>% 
#       group_split(MF) %>% 
#       map_dfr(~{
#         if (any(is.na(.x$Isotope))){
#           return(.x)
#         } else {
#           return(NULL)
#         }
#       })
#   ),
#   
#   ## Identify correctly assigned features
#   spiked_urine_correct_assignments = tar_target(
#     spiked_urine_correct_assignments,
#     spiked_urine_assignments %>% 
#       inner_join(spiked_urine_feature_matches %>% 
#                    select(Feature,MF,Adduct,Isotope), 
#                  by = c("Feature", 
#                         "Isotope", 
#                         "Adduct", 
#                         "MF"))
#   ),
#   
#   ## Identify assigned PIPs
#   tar_target(spiked_urine_assigned,
#              spiked_urine_feature_matches %>% 
#                inner_join(spiked_urine_assignments %>% 
#                             select(Feature),
#                           by = 'Feature')),
#   
#   ## Identify unassigned PIPs
#   tar_target(
#     spiked_urine_unassigned,
#     spiked_urine_feature_matches %>% 
#       anti_join(spiked_urine_assignments %>% 
#                   select(Feature),
#                 by = 'Feature')),
#   
#   ## Identify unassigned features that had correlations
#   tar_target(
#     spiked_urine_unassigned_correlated,
#     spiked_urine_unassigned %>% 
#       inner_join(spiked_urine_results_molecular_formula_assignment %>% 
#                    assignments::correlations() %>% 
#                    {bind_rows(
#                      select(.,Feature = Feature1),
#                      select(.,Feature = Feature2)
#                    )} %>% 
#                    distinct(),
#                  by = 'Feature'
#       ) 
#   ),
#   ## Identify unassigned features that did not have correlations
#   tar_target(
#     spiked_urine_unassigned_not_correlated,
#     spiked_urine_unassigned %>% 
#       anti_join(spiked_urine_results_molecular_formula_assignment %>% 
#                   assignments::correlations() %>% 
#                   {bind_rows(
#                     select(.,Feature = Feature1),
#                     select(.,Feature = Feature2)
#                   )} %>% 
#                   distinct(),
#                 by = 'Feature'
#       ) 
#   ),
#   
#   ## Identify unassigned features were correlated but did not have relationships
#   tar_target(
#     spiked_urine_unassigned_without_relationships,
#     spiked_urine_unassigned_correlated %>% 
#       anti_join(spiked_urine_results_molecular_formula_assignment %>% 
#                   assignments::relationships() %>% 
#                   dplyr::filter(
#                     is.na(Transformation1),
#                     is.na(Transformation2)) %>% 
#                   {bind_rows(
#                     select(.,Feature = Feature1),
#                     select(.,Feature = Feature2)
#                   )} %>% 
#                   distinct(),
#                 by = 'Feature')
#   ),
#   
#   
#   ## Identify unassigned features that had relationships
#   tar_target(
#     spiked_urine_unassigned_with_relationships,
#     spiked_urine_results_molecular_formula_assignment %>% 
#       assignments::relationships() %>% 
#       dplyr::filter(
#         is.na(Transformation1),
#         is.na(Transformation2),
#         Feature1 %in% 
#           spiked_urine_unassigned_correlated$Feature |
#           Feature2 %in% 
#           spiked_urine_unassigned_correlated$Feature) %>% 
#       {bind_rows(
#         select(.,Feature = Feature1),
#         select(.,Feature = Feature2)
#       )} %>% 
#       distinct() %>% 
#       dplyr::filter(Feature %in% 
#                       spiked_urine_unassigned_correlated$Feature)
#   ),
#   
#   ## Identify unassigned features were correlated and had correct relationships 
#   tar_target(
#     spiked_urine_unassigned_with_correct_relationships,
#     spiked_urine_feature_matches %>% 
#       inner_join(spiked_urine_results_molecular_formula_assignment %>% 
#                    assignments::relationships() %>% 
#                    dplyr::filter(
#                      is.na(Transformation1),
#                      is.na(Transformation2)) %>% 
#                    inner_join(spiked_urine_unassigned_correlated %>% 
#                                 select(Feature,Isotope,Adduct) %>% 
#                                 distinct(),
#                               by = c('Feature1' = 'Feature',
#                                      'Isotope1' = 'Isotope',
#                                      'Adduct1' = 'Adduct')) %>% 
#                    inner_join(spiked_urine_unassigned_correlated %>% 
#                                 select(Feature,Isotope,Adduct) %>% 
#                                 distinct(),
#                               by = c('Feature2' = 'Feature',
#                                      'Isotope2' = 'Isotope',
#                                      'Adduct2' = 'Adduct')) %>% 
#                    {bind_rows(
#                      select(.,
#                             Feature = Feature1,
#                             Isotope = Isotope1,
#                             Adduct = Adduct1),
#                      select(.,
#                             Feature = Feature2,
#                             Isotope = Isotope2,
#                             Adduct = Adduct2)
#                    )} %>% 
#                    distinct(),
#                  by = c('Feature','Isotope','Adduct')
#       )
#   ),
#   
#   ## Identify incorrectly assigned features
#   spiked_urine_incorrect_assignments = tar_target(
#     spiked_urine_incorrect_assignments,
#     spiked_urine_assignments %>% 
#       anti_join(spiked_urine_correct_assignments %>% 
#                   select(Feature), 
#                 by = "Feature") %>% 
#       inner_join(spiked_urine_feature_matches %>% 
#                    select(Feature), 
#                  by = "Feature") %>% 
#       distinct()
#   ),
#   
#   ## Identify correct relationships for the incorrectly assigned features
#   tar_target(
#     spiked_urine_incorrect_assigned_correct_relationships,
#     bind_rows(
#       spiked_urine_results_molecular_formula_assignment %>% 
#         assignments::relationships() %>% 
#         dplyr::filter(is.na(Transformation1),
#                       is.na(Transformation2)) %>% 
#         inner_join(spiked_urine_feature_matches %>% 
#                      select(Feature,Isotope,Adduct),
#                    by = c('Feature1' = 'Feature',
#                           'Isotope1' = 'Isotope',
#                           'Adduct1' = 'Adduct')) %>% 
#         select(Feature = Feature1,
#                Isotope = Isotope1,
#                Adduct = Adduct1),
#       spiked_urine_results_molecular_formula_assignment %>% 
#         assignments::relationships() %>% 
#         dplyr::filter(is.na(Transformation1),
#                       is.na(Transformation2)) %>% 
#         inner_join(spiked_urine_feature_matches %>% 
#                      select(Feature,Isotope,Adduct),
#                    by = c('Feature2' = 'Feature',
#                           'Isotope2' = 'Isotope',
#                           'Adduct2' = 'Adduct'))%>% 
#         select(Feature = Feature2,
#                Isotope = Isotope2,
#                Adduct = Adduct2)
#     ) %>% 
#       bind_rows() %>% 
#       distinct() %>% 
#       inner_join(spiked_urine_incorrect_assignments %>% 
#                    select(Feature) %>% 
#                    distinct(),
#                  by = 'Feature')
#   ),
#   
#   ## Identify incorrectly assigned features that did not have correct relationships
#   tar_target(
#     spiked_urine_incorrect_assigned_incorrect_relationships,
#     spiked_urine_feature_matches %>% 
#       inner_join(
#         bind_rows(
#           spiked_urine_results_molecular_formula_assignment %>% 
#             assignments::relationships() %>% 
#             dplyr::filter(is.na(Transformation1),
#                           is.na(Transformation2)) %>% 
#             inner_join(spiked_urine_incorrect_assignments %>% 
#                          select(Feature) %>% 
#                          distinct(),
#                        by = c('Feature1' = 'Feature')) %>% 
#             select(Feature = Feature1,
#                    Isotope = Isotope1,
#                    Adduct = Adduct1),
#           spiked_urine_results_molecular_formula_assignment %>% 
#             assignments::relationships() %>% 
#             dplyr::filter(is.na(Transformation1),
#                           is.na(Transformation2)) %>% 
#             inner_join(spiked_urine_incorrect_assignments %>% 
#                          select(Feature) %>% 
#                          distinct(),
#                        by = c('Feature2' = 'Feature'))%>% 
#             select(Feature = Feature2,
#                    Isotope = Isotope2,
#                    Adduct = Adduct2)
#         ) %>% 
#           distinct() %>% 
#           anti_join(spiked_urine_incorrect_assigned_correct_relationships %>% 
#                       select(Feature) %>% 
#                       distinct(),
#                     by = 'Feature') %>% 
#           select(Feature) %>% 
#           distinct(),
#         by = 'Feature')
#   ),
#   
#   ## Identify correct components for the incorrectly assigned features
#   tar_target(
#     spiked_urine_incorrect_assigned_correct_components,
#     spiked_urine_feature_matches %>% 
#       inner_join(spiked_urine_incorrect_assignments %>% 
#                    select(Feature),
#                  by = 'Feature') %>% 
#       inner_join(spiked_urine_incorrect_assigned_correct_relationships, 
#                  by = c("Feature", "Isotope", "Adduct")) %>% 
#       rowwise() %>% 
#       group_split() %>% 
#       map_dfr(~{
#         assignments::featureComponents(
#           spiked_urine_results_molecular_formula_assignment,
#           .x$Feature,
#           type = 'all'
#         ) %>% 
#           inner_join(.x %>% 
#                        select(Feature,MF,Isotope,Adduct),
#                      by = c('Feature','MF','Isotope','Adduct'))
#       })
#   ),
#   
#   ## Identify incorrectly assigned features without correct components
#   tar_target(
#     spiked_urine_incorrect_assigned_incorrect_components,
#     spiked_urine_incorrect_assigned_correct_relationships %>%
#       anti_join(spiked_urine_incorrect_assigned_correct_components %>%
#                   select(Feature,Isotope,Adduct),
#                 by = c("Feature", "Isotope", "Adduct"))
#   ),
#   
#   ## Retieved putative compound matches for incorrectly assigned features
#   spiked_urine_incorrect_assignment_compounds = tar_target(
#     spiked_urine_incorrect_assignment_compounds,
#     spiked_urine_results_molecular_formula_assignment %>% 
#       assignments::assignments() %>% 
#       dplyr::filter(Feature %in% spiked_urine_feature_matches$Feature) %>% 
#       dplyr::filter(!(MF %in% {spiked_urine_compound_db %>% 
#           cheminf::descriptors() %>% 
#           .$MF})) %>% 
#       dplyr::select(Feature,`Assigned MF` = MF) %>% 
#       dplyr::left_join(spiked_urine_feature_matches, 
#                        by = "Feature")
#   ),
#   
#   ## Identify standard compound MFs that it was possible for the algorithm to assign
#   tar_target(
#     spike_urine_standards_MF_assignment_possible,
#     bind_rows(
#       spiked_urine_unassigned_with_correct_relationships,
#       spiked_urine_feature_matches %>% 
#         inner_join(spiked_urine_incorrect_assigned_correct_relationships, 
#                    by = c("Feature", "Isotope", "Adduct")),
#       spiked_urine_correct_assignments
#     ) %>% 
#       select(MF) %>% 
#       distinct()
#   ),
#   
#   ## Plot component solutions for example feature
#   tar_target(
#     spiked_urine_feature_solutions_plot,
#     {
#       pl <- plotFeatureSolutions(
#         spiked_urine_results_molecular_formula_assignment,
#         tibble::tibble(
#           component = c(72,868,2370),
#           border = c('red',rep('black',2))
#         ),
#         spiked_urine_example_feature,
#         'A&I1'
#       ) &
#         guides(
#           edge_colour = ggraph::guide_edge_colorbar(
#             title = 'correlation coefficient',
#             title.position = 'top',
#             title.hjust = 0.5,
#             barheight = 0.75,
#             barwidth = 8,
#             frame.colour = 'black',
#             direction = 'horizontal'
#           )
#         )
#       
#       legend <- ggpubr::get_legend(pl)
#       
#       {pl &
#           theme(legend.position = 'none')} /
#         legend +
#         patchwork::plot_layout(heights = c(6,1))
#     }
#   )
# )

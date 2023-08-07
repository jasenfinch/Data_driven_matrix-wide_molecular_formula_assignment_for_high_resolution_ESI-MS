
manuscript_targets <- list(
  ## Approach plot
  tar_target(
    approach,
    assignmentApproach()
  ),
  
  ## Additional file info
  tarchetypes::tar_render(
    additional_files,
    'manuscript/additional_files.Rmd',
    output_file = 'Additional_files.pdf',
    output_dir = 'exports/manuscript/supplementary_materials'
  ),
  
  ## Supplementary materials and methods
  tarchetypes::tar_render(
    supplementary_methods,
    'manuscript/supplementary_methods.Rmd',
    output_file = 'Additional file 2.pdf',
    output_dir = 'exports/manuscript/supplementary_materials'
  ),
  
  ## Supplementary tables
  tarchetypes::tar_render(
    supplementary_tables,
    'manuscript/supplementary_tables.Rmd',
    output_file = 'Additional file 3.pdf',
    output_dir = 'exports/manuscript/supplementary_materials'
  ),
  
  ## Supplementary table S1 of spiked compounds
  tarchetypes::tar_file(
    table_s1,
    spiked_urine_compound_info %>% 
      dplyr::mutate(MF = SMILES %>% 
                      cheminf::chemicalDescriptors() %>% 
                      .$MF) %>% 
      dplyr::select(Name,InChI,`Pubchem CID` = CID,MF) %>% 
      jfmisc::exportCSV(
        file = 'exports/manuscript/supplementary_materials/Additional file 1.csv')
  ),
  
  ## Supplementary table S5 of all assignments
  tarchetypes::tar_file(
    table_s6,
    bind_rows(
      list(
        `Human_urine_FIE-HRMS` = spiked_urine_results_molecular_formula_assignment %>% 
          assignments::assignments(),
        `Brachypodium distachyon_leaf_FIE-HRMS` = brachy_FIE_HRMS_results_molecular_formula_assignment %>% 
          assignments::assignments(),
        `Brachypodium distachyon_leaf_C18-UHPLC-HRMS` = brachy_RP_UHPLC_HRMS_results_molecular_formula_assignment %>% 
          assignments::assignments()
      ),
      .id = 'Type'
    )  %>%
      tidyr::separate(
        Type,
        c('Organism',
          'Matrix',
          'Technique'),
        sep = '_'
      ) %>% 
      left_join(
        bind_rows(
          list(
            `Human_urine_FIE-HRMS` = spiked_urine_KEGG_matches,
            `Brachypodium distachyon_leaf_FIE-HRMS` = brachy_FIE_HRMS_KEGG_matches,
            `Brachypodium distachyon_leaf_C18-UHPLC-HRMS` = brachy_RP_UHPLC_HRMS_KEGG_matches
          ),
          .id = 'Type'
        )  %>%
          tidyr::separate(
            Type,
            c('Organism',
              'Matrix',
              'Technique'),
            sep = '_'
          ) %>% 
          select(-InChI,-SMILES) %>% 
          group_split(Organism,Matrix,Technique,Adduct,MF) %>% 
          map_dfr(~.x %>% 
                    mutate(`KEGG ID` = paste0(ID,collapse = '; ')) %>% 
                    .[1,]
          ) %>% 
          select(-ID),
        by = c('Organism','Matrix','Technique','Adduct','MF')
      ) %>% 
      select(Organism:Technique,
             `Ionisation mode` = Mode,
             `Retention time` = RetentionTime,
             `m/z` = `Measured m/z`,
             `Molecular formula` = MF,
             Isotope,
             Adduct,
             `PPM error`,
             `MF Plausibility (%)`,
             `KEGG ID`) %>% 
      jfmisc::exportCSV(
        file = 'exports/manuscript/supplementary_materials/Additional file 4.csv'
      )
  ),
  
  ## Supplementary table S6 directly matched chemical standard ionisation products
  tarchetypes::tar_file(
    table_s7,
    spiked_urine_feature_matches %>% 
      select(MF,
             Isotope,
             Adduct,
             `m/z` = mz,
             `Ionisation mode` = Mode,
             `PPM error`) %>% 
      jfmisc::exportCSV(
        file = 'exports/manuscript/supplementary_materials/Additional file 5.csv'
      )
  ),
  
  ## Supplementary table S7 assigned chemical standards
  tarchetypes::tar_file(
    table_s8,
    spiked_urine_correct_assignments %>% 
      select(MF,
             Isotope,
             Adduct,
             `m/z` = `Measured m/z`,
             `Ionisation mode` = Mode,
             `PPM error`) %>% 
      jfmisc::exportCSV(
        file = 'exports/manuscript/supplementary_materials/Additional file 6.csv'
      )
  ),
  
  ## zip supplementary material
  tarchetypes::tar_file(
    supplementary_zip,
    {
      
      zipfile <- 'exports/manuscript/supplementary_materials.zip'
      
      if (file.exists(zipfile)) {
        unlink(zipfile) 
      }
      
      zip(zipfile,
          c(supplementary_methods[1],
            supplementary_tables[1],
            table_s1,
            table_s6,
            table_s7,
            table_s8),
          flags = '-r9Xj')
      
      zipfile
    }
  ),
  
  
  ## render manuscript
  tarchetypes::tar_render(
    manuscript,
    "manuscript/manuscript.Rmd",
    output_dir = "exports/manuscript",
    quiet = TRUE,output_format = "all"
  )
)


KEGG_compound_IRs <- tarchetypes::tar_map(
  values = tibble(
    organism = c('human','brachy'),
    KEGG_id = c('hsa','bdi')
  ),
  ## Identify KEGG compounds found organism
  tar_target(
    KEGG_IDs,
    {
      enzymes <- KEGGREST::keggLink(KEGG_id,'enzyme') %>%
        names()
      KEGGREST::keggLink('compound','enzyme') %>%
        {tibble::tibble(Enzyme = names(.),
                        Compound = .)} %>%
        dplyr::filter(Enzyme %in% enzymes) %>%
        dplyr::select(Compound) %>%
        dplyr::distinct() %>%
        tibble::deframe() %>%
        stringr::str_remove_all('cpd:')
    }
  ),
  ## Filter InChIs to only include specified organism
  tar_target(
    specific_KEGG_InChIs,{
      future::plan('multicore',workers = 16)
      
      specific_KEGG_InChIs <- KEGG_InChIs %>% 
        dplyr::filter(ID %in% KEGG_IDs) %>% 
        tidyr::drop_na() %>% 
        dplyr::rowwise() %>% 
        group_split() %>% 
        furrr::future_map_dfr(~{
          .x %>%
            dplyr::mutate(
              SMILES = cheminf::convert(InChI,
                                        'INCHI',
                                        'SMILES'))
        })
      
      future::plan('multicore',workers = 64)
      
      specific_KEGG_InChIs
    }
  ),
  
  ## Calculate the chemical descriptors for orgainism specific KEGG compounds
  tar_target(
    KEGG_IRs,
    specific_KEGG_InChIs %>%
      dplyr::mutate(NAME = "") %>% 
      cheminf::metaboliteDB(.)
  ),
  
  
  names = 'organism'
)

KEGG_compound_targets <- list(
  ## Retrieve KEGG compounds
  tar_target(
    KEGG_compounds,
    KEGGREST::keggList('compound') %>% 
      names() %>% 
      stringr::str_remove('cpd:') %>% 
      tibble::tibble(ENTRY = .) %>% 
      split(rep(1:ceiling(nrow(.)/10), each=10, length.out=nrow(.))) %>%
      purrr::map_dfr(~{
        .x$ENTRY %>% 
          KEGGREST::keggGet() %>% 
          purrr::map_dfr(~{
            purrr::map_dfc(.x,~ paste(.x,collapse = ';;;'))
          })
      }
      )
  ), 
  
  ## Collate KEGG molecular formulas
  tar_target(
    KEGG_compounds_MFs,
    KEGG_compounds %>% 
      dplyr::rowwise() %>% 
      dplyr::select(MF = FORMULA) %>% 
      dplyr::distinct()
  ),
  
  ## Remove invalid KEGG moleuclar formulas
  tar_target(
    KEGG_compounds_valid_MFs,
    KEGG_compounds_MFs %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(valid_MF = validMF(MF)) %>% 
      dplyr::filter(valid_MF == TRUE) %>%
      dplyr::select(-valid_MF) %>% 
      dplyr::mutate(M = mzAnnotation::calcAccurateMass(MF)) %>% 
      dplyr::filter(M >= 50,
                    M <= 1000) %>% 
      dplyr::rowwise() %>% 
      targets::tar_group(),
    iteration = 'group'
  ),
  
  tarchetypes::tar_file(
    KEGG_InChIs_file,
    {
      piggyback::pb_download('KEGG_InChI.csv',
                             dest = 'data/',
                             tag = 'KEGG-InChI')
      'data/KEGG_InChI.csv' 
    }),
  
  ## Parse KEGG InChIs
  tar_target(
    KEGG_InChIs,
    readr::read_csv(KEGG_InChIs_file)
  ),
  
  KEGG_compound_IRs
)


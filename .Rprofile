source("renv/activate.R")

renv::restore(prompt = FALSE)

source("R/utils.R")

pacman::p_load(magrittr,targets,purrr,dplyr,ggraph,
               kableExtra,jfmisc,tidyselect,
               install = FALSE)

"R/functions/" %>%
  list.files(full.names = TRUE) %>%
  purrr::walk(source)   

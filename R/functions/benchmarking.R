memUsage <- function(feature_data,max_cors,n_cores,interval){
  
  if (n_cores > 1) {
    future::plan('multicore',workers = n_cores)  
  } else {
    future::plan('sequential')
  }
  
  parameters <- assignments::assignmentParameters('FIE-HRMS')
  parameters@correlations_parameters$maxCor <- max_cors
  
  mem_usage <- profvis::profvis({
    assignment <- assignments::assignMFs(
      feature_data,
      parameters) 
  },
  interval = interval) %>%
    {
      .$x$message$prof 
    } %>%
    tibble::as_tibble() %>%
    {
      tibble::tibble(
        min = min(.$memalloc),
        max = max(.$memalloc)) %>%
        dplyr::mutate(
          `RAM usage (MB)` = max - min,
          workers = n_cores,
          `# correlations` = max_cors) %>%
        dplyr::select(-min,-max)
    }
  
  return(mem_usage)
}

procBenchmark <- function(feature_data,max_cors,n_cores,reps){
  if (n_cores > 1) {
    future::plan('multicore',workers = n_cores)  
  } else {
    future::plan('sequential')
  }
  
  parameters <- assignments::assignmentParameters('FIE-HRMS')
  parameters@correlations_parameters$maxCor <- max_cors
  
  proc_time <- rbenchmark::benchmark(
    {
      assignment <- assignments::assignMFs(
        feature_data,
        parameters) 
    },
    replications = reps,
    columns = c('replications',
                'elapsed')
  ) %>%
    tibble::as_tibble() %>% 
    dplyr::rename(`time (seconds)` = elapsed) %>% 
    dplyr::mutate(workers = n_cores,
           `# correlations` = max_cors)
  
  return(proc_time)
}

n_correlations <- seq(3,5,1) %>% 
  {10^. * 2.5}
n_cores <- seq(32,96,32)

benchmark_intervals <- tidyr::expand_grid(
  n_correlations = n_correlations,
  n_cores = n_cores
)

interval <- 0.5
n_reps <- 1

benchmarking_branches <- tarchetypes::tar_map(
  values = benchmark_intervals,
  targets::tar_target(
    processing_time,
    procBenchmark(
      spiked_urine_results_pre_treatment,
      n_correlations,
      n_cores,
      reps = n_reps
    )
  )
)

benchmarking_targets <- list(
  benchmarking_branches,
  
  tarchetypes::tar_combine(
    processing_time_results,
    benchmarking_branches,
    command = dplyr::bind_rows(!!!.x,.id = 'combination')
  ),
  
  targets::tar_target(
    processing_time_plot,
    processing_time_results %>% 
      mutate(`Time (minutes)` = `time (seconds)` / 60,
             workers = factor(workers),
             `# correlations` = as.integer(`# correlations`)
      ) %>% 
      ggplot(aes(
        x = `# correlations`,
        y = `Time (minutes)`,
        fill = workers
      )) +
      geom_line(linetype = 2) +
      geom_point(shape = 21,
                 size = 3) +
      scale_x_log10(labels = scales::comma,breaks = n_correlations) +
      ggthemes::scale_fill_ptol() +
      jfmisc::theme_neat() +
      labs(x = 'No. correlations',
           fill = 'CPU workers')
  ),
  
  targets::tar_target(
    cpu_information,
    minder::cpuInfo()
  ),
  
  targets::tar_target(
    memory_information,
    minder::memoryInfo() 
  )
)


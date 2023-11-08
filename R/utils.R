## Resolve conflicts
# conflict_prefer(quiet = TRUE)

# options(future.debug = TRUE)

## Set targets options
targets::tar_option_set(
  error = 'continue',
  memory = 'transient',
  garbage_collection = TRUE
)

## Parallel backend
metaboMisc::suitableParallelPlan(
  workers = jfmisc::suitableParallelWorkers(
    RAM_per_worker = 12
  )
)

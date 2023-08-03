# Data-driven, matrix-wide molecular formula assignment for ultra-high resolution ESI-MS

[![Docker](https://github.com/jasenfinch/Data_driven_matrix-wide_molecular_formula_assignment_for_high_resolution_ESI-MS/workflows/Docker/badge.svg)](https://github.com/jasenfinch/Data_driven_matrix-wide_molecular_formula_assignment_for_high_resolution_ESI-MS/actions)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

This is the code and analysis repository for the article:

Finch, J.P., Wilson, T., Lyons, L., Phillips, H., Beckmann, M. and Draper, J., 2023. Data-driven, matrix-wide molecular formula assignment for ultra-high resolution ESI-MS]

All code is written in R and the [targets](https://docs.ropensci.org/targets/) package has been used for workflow management.
The [renv](https://github.com/rstudio/renv) package has been used to ensure a reproducible R environment.

## Compile the manuscript

### Using docker

The manuscript can be compiled using a pre-built docker image, directly from GitHub:

``` sh
docker run -v $(pwd):/home/Data_driven_matrix-wide_molecular_formula_assignment_for_high_resolution_ESI-MS ghcr.io/jasenfinch/molecular-formula-assignment:latest
```

### Locally

To generate the manuscript, simply [clone the repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository), open the R console, set the working directory to the repository clone using `setwd()` and run the command `targets::tar_make()`.

## Generated by projecttemplates (https://jasenfinch.github.io/projecttemplates/) v0.5.0
FROM ghcr.io/jasenfinch/hrm-docker:2023-05-07

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y cmake

RUN Rscript -e "install.packages(c('renv'), repos = c(CRAN = 'https://cloud.r-project.org'))"

WORKDIR /home/Data_driven_matrix-wide_molecular_formula_assignment_for_high_resolution_ESI-MS
    
COPY renv.lock renv.lock

RUN Rscript -e "renv::consent(provided = TRUE); renv::restore(prompt = FALSE)"

ENTRYPOINT ["Rscript","run.R"]

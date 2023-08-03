plotFeatureSolutions <- function(assignment,components,feature,iteration){
  components %>%
    rowwise() %>% 
    group_split() %>% 
    map(~{
      component_stats <- assignments:::component(
        assignment,
        .x$component,
        iteration,
        type = 'all'
      ) %>% 
        assignments::nodes() %>% 
        dplyr::filter(Feature == feature) %>% 
        mutate(across(where(is.numeric),~signif(.x,digits = 3)))
      
      assignments::plotComponent(
        assignment,
        .x$component,
        iteration,
        type = 'all',
        highlight = feature,
        border = .x$border,
        label_size = 2.5,
        axis_offset = 0.2
      ) +
        labs(  caption = glue::glue('
          P<sub>c</sub> = {component_stats$`Component Plausibility`};
          Degree = {component_stats$Degree};
          AIS<sub>c</sub> = {component_stats$AIS};
          
          P<sub>MF</sub> = {component_stats$`MF Plausibility (%)`}%;
          Δ ppm = {component_stats$`PPM error`}')
        )
    }) %>% 
    patchwork::wrap_plots() +
    patchwork::plot_layout(guides = 'collect')
}

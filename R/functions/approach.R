
assignmentApproach <- function(){
  nodes <- tibble::tibble(
    title = c(
      '',
      'Correlation analysis',
      'Relationship calculation',
      'Adduct and isotopic assignment iterations',
      'Transformation assignment iterations',
      ''),
    body = c(
      '*m/z* feature intensity matrix',
      'Create a network graph of linked features.',
      'Calculation of isotopic, adduct and transformation relationships between correlated features.',
      'Component subgraphs are extracted based on molecular formulas generated from ionisation products of related features. Feature assignments are selected by sequential elimination from components based on component plausibility scores.',
      'Identify features that are linked to assigned nodes by correlations with a transformation relationship. Multiple iterations are performed to allow possible assignments to propagate through the correlation network as more features are assigned.',
      'Assigned *m/z* features'
    ),
    x = 1,
    y = c(5,4.3,3.5,2.5,1.25,0.25)
  ) %>% 
    dplyr::mutate(label = body %>% 
             stringr::str_wrap(65) %>% 
             stringr::str_replace_all('\\n','<br>') %>% 
             {glue::glue('**{title}**<br>{.}')} %>% 
               stringr::str_replace_all(stringr::coll('****<br>'),'')
             )
  
  edges <- tibble::tibble(
    x = c(1),
    xend = c(1),
    y = c(4.85,4.07,3.22,2.07,0.76),
    yend = c(4.53,3.79,2.93,1.75,0.39)
  )
  
  ggplot2::ggplot() +
    ggtext::geom_richtext(
      data = nodes,
      ggplot2::aes(x = x,
                   y = y,
                   label = label),
      size = 3) +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(
        x = x,
        xend = xend,
        y = y,
        yend = yend
      ),
      arrow = ggplot2::arrow(length = ggplot2::unit(2,'mm')),
      size = 0.5
    ) +
    # ggplot2::theme_bw()
    ggraph::theme_graph()
}

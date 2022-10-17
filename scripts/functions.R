
# Functions ----

source(here::here("report_input", "custom_settings.R"))

## WebGestalt functions for PPI-networks ----

### filter number of terms according to FDR status ----
WGR_filterFDR_ppi <- function(WGR_prepped_df, n_terms = 20){
  if(is.null(WGR_prepped_df)){
    return(NULL)
  }
  
  unique_terms <- WGR_prepped_df %>% distinct(description) %>% nrow()
  unique_terms_FDRthreshold <- WGR_prepped_df %>% 
    # distinct(description, p_threshold) %>% 
    filter(p_threshold == TRUE) %>% 
    nrow()
  to_filter <- case_when(unique_terms > n_terms & unique_terms_FDRthreshold >= n_terms ~ "FDR_threshold", 
                         unique_terms > n_terms & unique_terms_FDRthreshold < n_terms ~ "n_min_FDR",
                         TRUE ~ "no_filter"
  )
  
  # If there is more than n_terms terms and min n_terms of these are above the FDR threshold, only above FDR threshold is kept
  if(to_filter == "FDR_threshold"){
    WGR_prepped_df_filtered <- WGR_prepped_df %>% filter(p_threshold == TRUE)}
  
  # If there is more than n_terms terms and none of these are significant, n_terms of these are kept
  if(to_filter == "n_min_FDR"){
    WGR_prepped_df_filtered <- WGR_prepped_df %>% ungroup() %>% slice_min(FDR, n = n_terms)
  }
  
  # If we have less than n_terms terms, all of them are kept
  if(to_filter == "no_filter"){
    WGR_prepped_df_filtered <- WGR_prepped_df
  }
  
  return(WGR_prepped_df_filtered)
}


### make figures ----
WGR_colplot_ppi <- function(WGR_prepped_df, n_rows){
  if(n_rows == 0){
    enrichResult_message <- "Unfortunately, no enriched terms were found for the selected targets."
    return(enrichResult_message)
  } 
  
  WGR_prepped_df <- WGR_prepped_df %>% 
    rowwise() %>% 
    mutate(description = ifelse(str_length(description) > 80,
                                str_wrap(description, width = plyr::round_any(str_length(description)/2, 10, f = ceiling)), 
                                description)) %>% 
    ungroup()
  
  description_order <- WGR_prepped_df %>% group_by(description) %>% summarise(enrichmentRatio_max = max(enrichmentRatio)) %>% arrange(desc(enrichmentRatio_max)) %>% pull(description)
  
  plot_data <- WGR_prepped_df %>% 
    mutate(description = fct_relevel(description, description_order)) %>%
    group_by(description, myhit_name) %>% 
    filter(enrichmentRatio == max(enrichmentRatio)) %>% 
    ungroup() %>% 
    arrange(description, desc(enrichmentRatio)) %>% 
    mutate(order = desc(row_number())) 
  
  plot <- plot_data %>%   
    ggplot(aes(y = enrichmentRatio, x = order)) + 
    geom_bar(position = "identity", stat = "identity", aes(fill = factor(p_threshold, c(TRUE, FALSE)))) +
    coord_flip() +
    facet_grid(rows = vars(description), scales = "free_y", space = "free", switch = "y") + 
    scale_x_continuous(
      breaks = plot_data$order,
      labels = plot_data$myhit_name,    
      expand = c(0,0)) +
    labs(x = NULL, y = "Enrichment ratio") + 
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          legend.box = "vertical"
          ) + 
    scale_fill_manual(breaks = c(TRUE, FALSE),
                      labels = c("FDR < 0.05", "FDR \u2265 0.05"),
                      values = c(figure_colors$p_threshold_true, figure_colors$p_threshold_false),
                      name = NULL,
                      drop = FALSE
                      )
  
  return(plot)
  
  }
 

## WebGestalt functions for shared features ----

### run WGR ----

WebGestaltR_safe <- safely(WebGestaltR)

WGR_run_shared <- function(interest_data, ref_data, database){
  interest_file <- interest_data %>% distinct(uniprot)
  ref_file <- ref_data %>% distinct(uniprot)
  enrichResult <- NULL
  
  enrichResult <- WebGestaltR_safe(enrichMethod = "ORA", 
                              organism = "hsapiens",
                              enrichDatabase = database,
                              interestGene = interest_file$uniprot, interestGeneType = "uniprotswissprot",
                              referenceGene = ref_file$uniprot, referenceGeneType = "uniprotswissprot", 
                              fdrThr = 1, minNum = 10, fdrMethod = "BH",
                              isOutput = FALSE, 
                              projectName = NULL)  
  
  return(enrichResult)
  
}


### prep WGR result ----
WGR_prep_shared <- function(WGR_result){
  if(is.null(WGR_result)){
    return(NULL) 
    }
  
  enrichResult_prep <- WGR_result %>% 
      mutate(userId = str_split(userId, "\\;")) %>% 
      unnest(userId) %>% 
      left_join(my_hits_bmk %>% select(target, uniprot), ., by = c("uniprot" = "userId")) %>% 
      filter(!is.na(description)) %>% 
      group_by(description, link) %>% 
      summarise(targets = str_c(target, collapse = ", ") %>% str_wrap(width = 20), 
                enrichmentRatio = first(enrichmentRatio),
                FDR = first(FDR)) %>% 
      mutate(FDR_threshold = ifelse(FDR < 0.05, TRUE, FALSE)) 
  return(enrichResult_prep)
  
}


### filter number of terms according to FDR status ----
WGR_filterFDR_shared <- function(WGR_prepped_df, n_terms = 20){
  if(is.null(WGR_prepped_df)){
    return(NULL)
  }
  
  unique_terms <- WGR_prepped_df %>% distinct(description) %>% nrow()
  unique_terms_FDRthreshold <- WGR_prepped_df %>% distinct(description, FDR_threshold) %>% filter(FDR_threshold == TRUE) %>% nrow()
  to_filter <- case_when(unique_terms > n_terms & unique_terms_FDRthreshold >= n_terms ~ "FDR_threshold", 
                         unique_terms > n_terms & unique_terms_FDRthreshold < n_terms ~ "n_min_FDR",
                         TRUE ~ "no_filter"
                         )
  
  # If there is more than n_terms terms and min n_terms of these are above the FDR threshold, only above FDR threshold is kept
  if(to_filter == "FDR_threshold"){
    WGR_prepped_df_filtered <- WGR_prepped_df %>% filter(FDR_threshold == TRUE)}
  
  # If there is more than n_terms terms and none of these are significant, n_terms of these are kept
  if(to_filter == "n_min_FDR"){
    WGR_prepped_df_filtered <- WGR_prepped_df %>% ungroup() %>% slice_min(FDR, n = n_terms)
  }
  
  # If we have less than n_terms terms, all of them are kept
  if(to_filter == "no_filter"){
    WGR_prepped_df_filtered <- WGR_prepped_df
  }
  
  return(WGR_prepped_df_filtered)
}


### calculate figure height WGR ----
WGR_figureheight_shared <- function(WGR_prepped_df){ 
  if(is.null(WGR_prepped_df)){
    return(1)
  }
  
  unique_goterms <- WGR_prepped_df %>% distinct(description) %>% nrow()
  max_rows <- WGR_prepped_df %>% summarise(count = str_count(targets, "\n") + 1) %>% ungroup() %>% slice_max(count) %>% distinct(count) %>% pull(count)
  goterm_height <- unique_goterms * 0.3 * max_rows + 1.5
  
  return(goterm_height)
}


### make figures ---- 
WGR_colplot_shared <- function(WGR_prepped_df){
  if(is.null(WGR_prepped_df)){
    enrichResult_message <- "Unfortunately, no enriched terms were found for the selected targets."
    return(enrichResult_message)}
  
  WGR_prepped_df <- WGR_prepped_df %>% 
    rowwise() %>% 
    mutate(description = ifelse(str_length(description) > 60,
                                str_wrap(description, width = plyr::round_any(str_length(description)/2, 10, f = ceiling)), 
                                description)) %>% 
    ungroup()
  
  enrichResult_plot <-
    WGR_prepped_df %>%
    ungroup() %>% 
    mutate(description = fct_reorder(description, enrichmentRatio)) %>% 
    ggplot(aes(x = enrichmentRatio, y = description)) +
    labs(x = "Enrichment ratio", y = NULL) +
    geom_col(aes(alpha = factor(FDR_threshold, c(TRUE, FALSE)), fill = factor(FDR_threshold, c(TRUE, FALSE)))) +
    scale_alpha_manual(breaks = c(TRUE, FALSE),
                       labels = c("FDR < 0.05", "FDR \u2265 0.05"),
                       values = c(1,0.75),
                       drop = FALSE,
                       name = NULL
                       ) +
    scale_fill_manual(breaks = c(TRUE, FALSE),
                      labels = c("FDR < 0.05", "FDR \u2265 0.05"),
                      values = c(figure_colors$p_threshold_true, figure_colors$p_threshold_false),
                      name = NULL,
                      drop = FALSE,
                      ) +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.box = "vertical",
          text = element_text(size = 12),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.y = unit(1.5, "lines")
          ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE),
           alpha = guide_legend(nrow = 1, byrow = TRUE)) + 
    coord_cartesian(expand = FALSE)  
  
  
  enrichResult_text <- 
    WGR_prepped_df %>%
    ungroup() %>% 
    mutate(description = fct_reorder(description, enrichmentRatio)) %>% 
    ggplot(aes(y = description)) +
    geom_text(x = 0.5, 
              aes(label = targets, 
                  color = factor(FDR_threshold, c(TRUE, FALSE))), 
              size = 4,
              vjust = 0.45
              ) +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.box = "vertical",
          plot.margin = margin(t = 0, b = 0),
          panel.spacing.y = unit(1.5, "lines")
          ) + 
    scale_color_manual(breaks = c(TRUE, FALSE),
                       labels = c("FDR < 0.05", "FDR \u2265 0.05"),
                       values = c(figure_colors$p_threshold_true, figure_colors$p_threshold_false),
                       name = NULL,
                       drop = FALSE
                       ) +
    scale_alpha_manual(breaks = c(TRUE, FALSE),
                       values = c(1, 0.75),
                       ) +
    guides(color = "none",
           alpha = "none")
  
  enrichResult_figure <- enrichResult_plot + enrichResult_text + 
    plot_layout(widths = c(2.25,1)) &
    theme(legend.position = "top")
    
  return(enrichResult_figure)
  }




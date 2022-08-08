# settings ----- 

## Thresholds ---- 

# P-value threshold - Should we make an option to change this when making the report?
p_value_threshold <- 0.05

# Minimum number of hits - filtering on p-value will occur if there are more than this number of significant hits.
min_num_hits <- 10


# Set plot theme ----
theme_set(theme_classic() +
            theme(strip.background = element_blank(),
                  strip.text.y = element_text(angle = 0),
            )

)



## FIGURE COLORS ----

# doesn't work as lists...
about_text <- list(
  # Significance
  text_p_threshold_true = "blue",
  text_p_threshold_false = "gray",
  
  color_p_threshold_true = "#00688B",
  color_p_threshold_false = "gray",
  
  # diabetes colors
  text_diabetes_fraction_diabetes = "dark blue",
  text_diabetes_fraction_allother = "light blue",
  
  color_diabetes_fraction_diabetes = "#1F78B4",
  color_diabetes_fraction_allother = "#A6CEE3"
  
)

figure_colors <- list(
  
  # general 
  p_threshold_true = "#00688B", # #00688B = dark blue
  p_threshold_false = "gray",
  
  omics_legend_overwrite = "gray",
  
  # volcano plots w/ best results
  background_points = "lightgray",

  
  # diabetes fractions of co-mentions figure
  diabetes_fractions_diabetes = "#1F78B4", # #1F78B4 = darker blue
  diabetes_fractions_allother = "#A6CEE3", # #A6CEE3 = lighter blue
  
  # PPI networks
  prim_node_text = "#00688B", # #00688B = dark blue
  sec_node_text = "gray25",
  sec_node_confint_true = "#00688B",
  sec_node_confint_false = "gray",
  edges = "lightgray"
)

## TABLE COLORS ----
table_colors <- list(  
  
  # coloring text according to significance in tables
  text_p_treshold_true = "#00688B",  # #00688B = dark blue
  text_p_treshold_false = "gray",
  
  # color bars - lighter color
  color_bar_p_threshold_true = "#C5F4FF",  # #C5F4FF = lighter blue
  color_bar_p_threshold_false = "#DFDFE2" # #DFDFE2 = lighter grey 
)
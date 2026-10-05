
rm(list=ls())

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(scales)
library(tidyverse)

# Load data

df <- read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/01_raw_inputs/Biomphalaria_upd_complete_oct26.csv", sep = ",", header = TRUE)


table(df$country)

# ---------------------------
## EXPLORATORY ANALYSIS
# ---------------------------

# 1. Histogram per country (considering nature of the variables)

# Variable order 
ordered_vars <- c(
  "elev","TPI","SS","flow_acc",
  "temp_modis","pp_oneweek","NDVI_modis","NDVI_q90_t",
  "prop_water","prop_trees","prop_crops","prop_built",
  "pop_rate_change","sampling_experience","time_sin","time_cos"
)

# Which unit to use per variable 
site_level_vars <- c("elev","TPI","SS","flow_acc","pop_rate_change",
                     "prop_water","prop_trees","prop_crops","prop_built")
visit_level_vars <- setdiff(ordered_vars, site_level_vars)

#  Groups for color coding
var_groups <- tibble(
  var = ordered_vars,
  group = c(
    rep("Geomorphology and hydrology", 4),
    rep("Climate and vegetation",      4),
    rep("Land cover",               4),
    rep("Anthropogenic and detection",4)
  )
)

# Colorblind-friendly palette
pal <- c(
  "Geomorphology and hydrology" = "#1B333D",
  "Climate and vegetation"      = "#759FBC",
  "Land cover"                  = "#bfd3c1",
  "Anthropogenic and detection" = "#FFE7D8"
)

# Labels (edit freely)
var_labels <- tibble(
  var = ordered_vars,
  label = c(
    "Elevation", "Topographic index", "SSF", "Flow accumulation",
    "Temperature", "Precipitation (1-week)", "NDVI (MODIS)", "NDVI (Planet)",
    "Prop. water", "Prop. trees", "Prop. crops", "Prop. built",
    "Pop. rate change", "Sampling experience", "Sampling time (sin)", "Sampling time (cos)"
  )
)


# helpers
parse_num <- function(x) readr::parse_number(as.character(x))
first_non_na <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else x[1] }
median_or_na <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else median(x) }

make_plot_one_country <- function(df, ctry, bins = 30) {
  
  df_country <- df %>%
    filter(country == ctry) %>%
    mutate(across(all_of(ordered_vars), parse_num))
  
  df_site <- df_country %>%
    group_by(site_uid) %>%
    summarise(
      elev = first_non_na(elev),
      TPI  = first_non_na(TPI),
      SS   = first_non_na(SS),
      flow_acc = first_non_na(flow_acc),
      pop_rate_change = first_non_na(pop_rate_change),
      
      prop_water = median_or_na(prop_water),
      prop_trees = median_or_na(prop_trees),
      prop_crops = median_or_na(prop_crops),
      prop_built = median_or_na(prop_built),
      .groups = "drop"
    )
  
  d_site <- df_site %>%
    pivot_longer(cols = all_of(site_level_vars), names_to = "var", values_to = "value")
  
  d_visit <- df_country %>%
    pivot_longer(cols = all_of(visit_level_vars), names_to = "var", values_to = "value")
  
  d <- bind_rows(d_site, d_visit) %>%
    left_join(var_groups, by = "var") %>%
    left_join(var_labels, by = "var") %>%
    mutate(label = factor(label, levels = var_labels$label)) %>%
    filter(!is.na(value) & is.finite(value))
  
  ggplot(d, aes(x = value, fill = group)) +
    geom_histogram(aes(y = after_stat(ncount)),
                   bins = bins, color = "white", linewidth = 0.12, alpha = 0.92) +
    facet_wrap(~ label, ncol = 4, scales = "free_x", drop = FALSE) +
    scale_fill_manual(values = pal, guide = guide_legend(title = NULL)) +
    labs(title = paste0("Predictor distributions: ", ctry), x = NULL, y = "Relative frequency") +
    theme_classic(base_size = 10.5) +
    theme(
      legend.position = "bottom",
      axis.line = element_line(linewidth = 0.25, colour = "grey35"),
      axis.ticks = element_line(linewidth = 0.25, colour = "grey35"),
      strip.background = element_blank(),
      strip.text = element_text(face = "plain", size = 9),
      plot.title = element_text(face = "bold", size = 12)
    )
}

# ---- run for all countries (UG, DRC, CHAD) ----
countries <- c("UG","DRC","CHAD")

p_ug   <- make_plot_one_country(df, "UG")
p_drc  <- make_plot_one_country(df, "DRC")
p_chad <- make_plot_one_country(df, "CHAD")

# Save (PNG )
ggsave("APPENDIX_dist_UG_Oct26.png",   p_ug,   width = 180, height = 130, units = "mm", dpi = 300)
ggsave("APPENDIX_dist_DRC_Oct26.png",  p_drc,  width = 180, height = 130, units = "mm", dpi = 300)
ggsave("APPENDIX_dist_CHAD_Oct26.png", p_chad, width = 180, height = 130, units = "mm", dpi = 300)

# 2. Transformations 

# For flow accumulation and precipitation 

df$flow_acc = log(df$flow_acc+1)
df$pp_oneweek = log(df$pp_oneweek+1)

# Proportion - median 

df_site <- df %>%
  group_by(site_uid) %>%
  summarise(
    prop_water_med = median_or_na(prop_water),
    prop_trees_med = median_or_na(prop_trees),
    prop_crops_med = median_or_na(prop_crops),
    prop_built_med = median_or_na(prop_built),
    .groups = "drop"
  )

df = left_join(df, df_site, by = 'site_uid')

## after transformation

p_ug   <- make_plot_one_country(df, "UG")
p_drc  <- make_plot_one_country(df, "DRC")
p_chad <- make_plot_one_country(df, "CHAD")

# Save (PNG )
ggsave("APPENDIX_dist_UG_t_Oct26.png",   p_ug,   width = 180, height = 130, units = "mm", dpi = 300)
ggsave("APPENDIX_dist_DRC_t_Oct26.png",  p_drc,  width = 180, height = 130, units = "mm", dpi = 300)
ggsave("APPENDIX_dist_CHAD_t_Oct26.png", p_chad, width = 180, height = 130, units = "mm", dpi = 300)


# 3. Scaling 

rm(list = setdiff(ls(), c('df')))

# Define predictor columns
predictors <- c("NDVI_q90_t", "pp_oneweek", "temp_modis", "elev", "SS", "flow_acc", 
                "pop_rate_change", "prop_water_med", "prop_trees_med", "prop_crops_med", "prop_built_med",
                'sampling_experience', 'time_sin', 'time_cos', 'NDVI_modis', 'TPI')  #  names

# 3.1 Compute scaling parameters by country
scaling_params <- df %>%
  group_by(country) %>%
  summarise(across(all_of(predictors), 
                   list(mean = ~ mean(.x, na.rm = TRUE), 
                        sd   = ~ sd(.x,   na.rm = TRUE)))) %>%
  ungroup()

# Write to CSV
write.csv(scaling_params, "~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/02_clean_data/scaling_parameters_by_country_Oct26.csv", row.names = FALSE)

# 3.2 Apply scaling to data
df_scaled <- df %>%
  group_by(country) %>%
  mutate(across(all_of(predictors),
                ~ (.x - mean(.x, na.rm=TRUE)) / sd(.x, na.rm=TRUE))) %>%
  ungroup()

# 4. Correlation pairwise 

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# -------------------------------
# Helper for prettier variable labels
# -------------------------------
# -------------------------------
# Predictor labels
# -------------------------------

var_labels <- c(
  "elev"                = "Elevation",
  "TPI"                 = "Topographic index (TI)",
  "SS"                  = "Slope steepness factor (SSF)",
  "flow_acc"            = "Flow accumulation",
  "temp_modis"          = "Temperature (MODIS)",
  "pp_oneweek"          = "Precipitation (1-week)",
  "NDVI_modis"          = "NDVI (MODIS)",
  "NDVI_q90_t"          = "NDVI (Planet)",
  "prop_water_med"      = "Proportion water",
  "prop_trees_med"      = "Proportion trees",
  "prop_crops_med"      = "Proportion crops",
  "prop_built_med"      = "Proportion built",
  "pop_rate_change"     = "Population rate change",
  "sampling_experience" = "Sampling experience",
  "time_sin"            = "Sampling time (sin)",
  "time_cos"            = "Sampling time (cos)"
)

pretty_lab <- function(x) {
  
  out <- unname(var_labels[x])
  
  # Fallback for any variable not included above
  missing <- is.na(out)
  
  out[missing] <- x[missing] |>
    str_replace_all("_", " ") |>
    str_squish() |>
    str_to_sentence()
  
  out
}

# Keep a fixed predictor order across countries for comparability
pred_keep <- predictors[predictors %in% names(df)]

for (ctry in unique(df$country)) {
  
  df_ctry <- df %>%
    filter(country == ctry)
  
  mat <- df_ctry %>%
    dplyr::select(all_of(pred_keep)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
    dplyr::select(where(~ !all(is.na(.)))) %>%
    dplyr::select(where(~ sum(!is.na(.)) > 2))
  
  if (nrow(mat) < 3 || ncol(mat) < 2) next
  
  # Correct Pearson correlation
  corr_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
  
  # Fixed order for between-country comparison
  ord <- colnames(corr_mat)
  
  # Convert matrix to long format and keep lower triangle only
  corr_df <- as.data.frame(as.table(corr_mat), stringsAsFactors = FALSE) %>%
    rename(var1 = Var1, var2 = Var2, r = Freq) %>%
    mutate(i = match(var1, ord),
           j = match(var2, ord)) %>%
    filter(i > j) %>%
    mutate(
      var1 = factor(pretty_lab(var1), levels = rev(pretty_lab(ord))),
      var2 = factor(pretty_lab(var2), levels = pretty_lab(ord)),
      # Only print stronger correlations to reduce clutter
      lab = ifelse(abs(r) >= 0.5, sprintf("%.2f", r), "")
    )
  
  p <- ggplot(corr_df, aes(x = var2, y = var1, fill = r)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = lab), size = 2.8, family = "sans") +
    scale_fill_gradient2(
      low = "#3F618C",
      mid = "white",
      high = "#A65E58",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson's r"
    ) +
    coord_equal() +
    labs(
      x = NULL,
      y = NULL,
      title = ctry
    ) +
    theme_minimal(base_size = 9, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(b = 6)),
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    )
  # Vector export for manuscript
  ggsave(
    filename = paste0("corrplot_", ctry, ".pdf"),
    plot = p,
    width = 170 / 25.4,
    height = 140 / 25.4,
    device = cairo_pdf
  )
  
  # High-resolution raster export if needed
  ggsave(
    filename = paste0("corrplot_", ctry, ".png"),
    plot = p,
    width = 170 / 25.4,
    height = 140 / 25.4,
    dpi = 600
  )
}


## 5. Collinearity with VIF 
library(car)

# Extract only the fixed effect part of the model
# It's possible to extract the fixed effects part of the formula and refit without random effects for VIF calculation

write.csv(df_scaled, "~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/02_clean_data/df_final_scaled_Bio_upd_Oct26.csv", row.names = FALSE)

df_scaled = read.csv('~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/02_clean_data/df_final_scaled_Bio_upd_Oct26.csv')

## Calculate VIF 

calc_vif <- function(data, predictors) {
  
  dat <- data %>%
    dplyr::select(all_of(predictors)) %>%
    tidyr::drop_na()
  
  out <- sapply(predictors, function(v) {
    
    others <- setdiff(predictors, v)
    
    form <- reformulate(
      others,
      response = v
    )
    
    mod <- lm(
      form,
      data = dat
    )
    
    r2 <- summary(mod)$r.squared
    
    1 / (1 - r2)
  })
  
  tibble::tibble(
    Variable = names(out),
    VIF = as.numeric(out)
  ) %>%
    arrange(desc(VIF))
}

predictors <- c(
  "NDVI_q90_t",
  "pp_oneweek",
  "temp_modis",
  "elev",
  "TPI",
  "SS",
  "flow_acc",
  "pop_rate_change",
  "prop_water_med",
  "prop_trees_med",
  "prop_crops_med",
  "prop_built_med",
  "sampling_experience",
  "time_sin",
  "time_cos",
  "NDVI_modis"
)

vif_UG <- calc_vif(
  df_scaled %>% filter(country == "UG"),
  predictors
)

vif_DRC <- calc_vif(
  df_scaled %>% filter(country == "DRC"),
  predictors
)

vif_CHAD <- calc_vif(
  df_scaled %>% filter(country == "CHAD"),
  predictors
)

vif_UG
vif_DRC
vif_CHAD

library(dplyr)
library(tidyr)

# Add country labels
vif_UG <- vif_UG %>%
  mutate(Country = "Uganda")

vif_DRC <- vif_DRC %>%
  mutate(Country = "DRC")

vif_CHAD <- vif_CHAD %>%
  mutate(Country = "Chad")

# Combine
vif_table <- bind_rows(
  vif_DRC,
  vif_UG,
  vif_CHAD
) %>%
  mutate(
    Variable = dplyr::recode(
      Variable,
      "pp_oneweek"          = "Precipitation (1-week)",
      "temp_modis"          = "Temperature",
      "NDVI_q90_t"          = "NDVI (Planet)",
      "elev"                = "Elevation",
      "flow_acc"            = "Flow accumulation",
      "SS"                  = "Slope steepness factor",
      "TPI"                 = "Topographic wetness index (TI)",
      "pop_rate_change"     = "Population rate change",
      "prop_water_med"      = "Proportion water",
      "prop_trees_med"      = "Proportion trees",
      "prop_crops_med"      = "Proportion crops",
      "prop_built_med"      = "Proportion built",
      "sampling_experience" = "Sampling experience",
      "time_sin"            = "Sampling time (sin)",
      "time_cos"            = "Sampling time (cos)",
      "NDVI_modis"          = "NDVI (MODIS)"
    )
  ) %>%
  dplyr::select(Variable, Country, VIF) %>%
  pivot_wider(
    names_from = Country,
    values_from = VIF
  ) %>%
  dplyr::select(
    Variable,
    DRC,
    Uganda,
    Chad
  ) %>%
  mutate(
    across(
      c(DRC, Uganda, Chad),
      ~ round(.x, 2)
    )
  )

vif_table

# No VIF exceeding 10, so variables are retained 

variable_order <- c(
  "Precipitation (1-week)",
  "Temperature",
  "NDVI (Planet)",
  "Elevation",
  "Flow accumulation",
  "Slope steepness factor",
  "Topographic wetness index (TI)",
  "Population rate change",
  "Proportion water",
  "Proportion trees",
  "Proportion crops",
  "Proportion built",
  "Sampling experience",
  "Sampling time (sin)",
  "Sampling time (cos)",
  "NDVI (MODIS)"
)

vif_table <- vif_table %>%
  mutate(
    Variable = factor(
      Variable,
      levels = variable_order
    )
  ) %>%
  arrange(Variable)

vif_table

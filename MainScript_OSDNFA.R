
rm(list=ls())

#############

## Figure 4 
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpattern)

extract_fixed_terms <- function(f) {
  if (is.na(f) || f == "") return(character(0))
  
  rhs <- str_split_fixed(f, "~", 2)[2] |>
    str_remove_all("\\([^\\)]*\\|[^\\)]*\\)") |>   # remove random effects
    str_remove_all("offset\\([^\\)]*\\)") |>
    str_squish()
  
  terms <- str_split(rhs, "\\+", simplify = FALSE)[[1]] |> str_trim()
  terms[terms != "" & !terms %in% c("1", "0", "-1")]
}

extract_sig_terms <- function(s) {
  if (is.na(s) || s == "") return(character(0))
  out <- str_split(s, ",", simplify = FALSE)[[1]] |> str_trim()
  out[out != ""]
}

summarise_selected_models <- function(df, country) {
  df_sel <- df %>%
    filter(AIC - min(AIC, na.rm = TRUE) <= 2) %>%
    mutate(model_id = row_number())
  
  if (nrow(df_sel) == 0) stop(paste("No AIC<2 models for", country))
  
  var_long <- df_sel %>%
    transmute(model_id, Variable = lapply(as.character(Formula), extract_fixed_terms)) %>%
    unnest(Variable) %>%
    distinct()
  
  sig_long <- df_sel %>%
    transmute(model_id, Variable = lapply(as.character(Significant_Vars), extract_sig_terms)) %>%
    unnest(Variable) %>%
    distinct()
  
  sig_counts <- inner_join(var_long, sig_long, by = c("model_id", "Variable")) %>%
    dplyr::count(Variable, name = "SigCount")
  
  var_long %>%
    dplyr::count(Variable, name = "Count") %>%
    left_join(sig_counts, by = "Variable") %>%
    mutate(
      SigCount   = coalesce(SigCount, 0L),
      Country    = country,
      n_models   = nrow(df_sel),
      Percent    = 100 * Count / n_models,
      SigPercent = 100 * SigCount / n_models
    )
}

# Files
files <- c(
  Uganda = "~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/04_model_outputs/Models_output_UG_upd_Oct26.csv",
  DRC    = "~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/04_model_outputs/Models_output_DRC_upd_Oct26.csv",
  Chad   = "~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/04_model_outputs/Models_output_CHAD_upd_Oct26.csv"
)

UG_models   <- read.csv(files["Uganda"])
DRC_models  <- read.csv(files["DRC"])
CHAD_models <- read.csv(files["Chad"])

# Variable labels
var_lut <- c(
  "pp_oneweek"          = "Precipitation (1-week)",
  "temp_modis"          = "Temperature (MODIS)",
  "NDVI_modis"          = "NDVI (MODIS)",
  "NDVI_q90_t"          = "NDVI (Planet)",
  "flow_acc"            = "Flow accumulation",
  "SS"                  = "Slope steepness factor (SSF)",
  "elev"                = "Elevation",
  "pop_rate_change"     = "Population rate change",
  "prop_water_med"      = "Proportion of water",
  "prop_trees_med"      = "Proportion of trees",
  "prop_crops_med"      = "Proportion of crops",
  "prop_built_med"      = "Proportion of built",
  "sampling_experience" = "Sampling experience",
  "time_sin"            = "Sampling time (sin)",
  "time_cos"            = "Sampling time (cos)"
)

pal <- c("Uganda" = "#111A26", "DRC" = "#3F618C", "Chad" = "#6883A6")

# Build dataset
plot_df <- bind_rows(lapply(names(files), function(ctry) {
  summarise_selected_models(read.csv(files[[ctry]]), ctry)
})) %>%
  mutate(
    Variable_label = dplyr::recode(as.character(Variable), !!!var_lut, .default = as.character(Variable))
  )

# Order variables
var_order <- plot_df %>%
  group_by(Variable_label) %>%
  summarise(score = max(Percent, na.rm = TRUE), .groups = "drop") %>%
  arrange(score) %>%
  pull(Variable_label)

# Ensure all Variable x Country combinations exist
plot_df <- plot_df %>%
  complete(
    Variable_label,
    Country = c("Uganda", "DRC", "Chad"),
    fill = list(Count = 0, SigCount = 0, Percent = 0, SigPercent = 0)
  ) %>%
  mutate(
    Country = factor(Country, levels = c("Uganda", "DRC", "Chad")),
    Variable_label = factor(Variable_label, levels = var_order),
    sig_flag = SigPercent > 0
  )

pos <- position_dodge2(width = 0.82, preserve = "single")

p4a <- ggplot(plot_df, aes(x = Variable_label)) +
  geom_col(
    aes(y = Percent, fill = Country, group = Country),
    position = pos,
    width = 0.78,
    color = "grey15",
    linewidth = 0.35
  ) +
  geom_col_pattern(
    aes(
      y = SigPercent,
      group = Country,
      colour = Country,
      pattern = sig_flag
    ),
    position = pos,
    width = 0.78,
    fill = "transparent",
    pattern = "stripe",
    pattern_angle = 45,
    pattern_alpha = 0.4,
    pattern_density = 0.2,
    pattern_size = 0.6,
    pattern_spacing = 0.04,
    pattern_frequency = 2,
    pattern_colour = "grey90",
    linewidth = 0.55,
    show.legend = FALSE
  ) +
  scale_pattern_manual(values = c(`FALSE` = "none", `TRUE` = "stripe")) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Percentage of AIC<2 models (%)",
    fill = "Country"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.y = element_text(size = 8)
  )

print(p4a)

#### Select best model 

library(dplyr)

select_best_row <- function(df,
                            sin_var = "time_sin",
                            cos_var = "time_cos",
                            aic_window = 2,
                            formula_col = "Formula") {
  
  cand <- df %>%
    filter(AIC - min(AIC, na.rm = TRUE) < aic_window)
  
  has_sin <- grepl(paste0("\\b", sin_var, "\\b"), cand[[formula_col]])
  has_cos <- grepl(paste0("\\b", cos_var, "\\b"), cand[[formula_col]])
  
  cand_ok <- cand %>%
    mutate(has_sin = has_sin, has_cos = has_cos) %>%
    filter((has_sin & has_cos) | (!has_sin & !has_cos))
  
  # If enforcing the rule leaves nothing, fall back to the original candidate set
  if (nrow(cand_ok) == 0) cand_ok <- cand
  
  cand_ok %>%
    slice_max(order_by = R2_m, n = 1, with_ties = FALSE)
}

UG_best_row   <- select_best_row(UG_models,   sin_var="time_sin", cos_var="time_cos")
DRC_best_row  <- select_best_row(DRC_models,  sin_var="time_sin", cos_var="time_cos")
CHAD_best_row <- select_best_row(CHAD_models, sin_var="time_sin", cos_var="time_cos")

best_rows <- bind_rows(
  mutate(UG_best_row,   Country = "Uganda"),
  mutate(DRC_best_row,  Country = "DRC"),
  mutate(CHAD_best_row, Country = "Chad")
)

best_rows[, c("Country","Formula","AIC","R2_m","R2_c")]

library(lme4)

fit_from_string <- function(formula_string, data, family = binomial()) {
  f <- as.formula(formula_string)
  glmer(f, data = data, family = family,
        control = glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 2e5)))
}

df_all <- read.csv('~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/02_clean_data/df_final_scaled_Bio_upd_Oct26.csv')

UG_data   <- df_all %>% filter(country == "UG")
DRC_data  <- df_all %>% filter(country == "DRC")
CHAD_data <- df_all %>% filter(country == "CHAD")

UG_best_model   <- fit_from_string(UG_best_row$Formula[1],   UG_data)
DRC_best_model  <- fit_from_string(DRC_best_row$Formula[1],  DRC_data)
CHAD_best_model <- fit_from_string(CHAD_best_row$Formula[1], CHAD_data)

library(broom.mixed)
library(tidyr)

tidy_or <- function(mod, country_label){
  broom.mixed::tidy(mod, conf.int = TRUE) %>%
    filter(effect == "fixed", term != "(Intercept)") %>%
    mutate(
      Country = country_label,
      odds_ratio = exp(estimate),
      lower_ci   = exp(conf.low),
      upper_ci   = exp(conf.high)
    ) %>%
    dplyr::select(Country, term, odds_ratio, lower_ci, upper_ci)
}

or_UG   <- tidy_or(UG_best_model,   "Uganda")
or_DRC  <- tidy_or(DRC_best_model,  "DRC")
or_CHAD <- tidy_or(CHAD_best_model, "Chad")

or_all <- bind_rows(or_UG, or_DRC, or_CHAD)

all_terms <- sort(unique(or_all$term))

or_all_aligned <- or_all %>%
  group_by(Country) %>%
  complete(term = all_terms,
           fill = list(odds_ratio = NA, lower_ci = NA, upper_ci = NA)) %>%
  ungroup()

or_all_aligned <- or_all_aligned %>%
  mutate(term_label = dplyr::recode(term, !!!var_lut, .default = term))

# Order predictors by (say) Uganda OR magnitude, or just alphabetically:
term_order <- or_all_aligned %>%
  pull(term_label) %>%
  unique() %>%
  sort(decreasing=T)

or_all_aligned$term_label <- factor(or_all_aligned$term_label, levels = term_order)

library(ggplot2)

# Palette
pal_country <- c("Uganda"="#111A26","DRC"="#3F618C","Chad"="#6883A6")

# Optional: distinct shapes per country (helps if colors are similar)
#shape_country <- c("Chad" = 16, "DRC" = 16, "Uganda" = 16)  # all circles (clean)
shape_country <- c("Chad" = 16, "DRC" = 17, "Uganda" = 15)

pos <- position_dodge(width = 0.55)

p <- ggplot(or_all_aligned, aes(x = term_label, y = odds_ratio,
                                color = Country, shape = Country)) +
  geom_pointrange(
    aes(ymin = lower_ci, ymax = upper_ci),
    position = pos,
    linewidth = 0.9
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  scale_y_log10(
    breaks = c(0.1, 0.3, 1, 3, 10, 30, 100),
    labels = c("0.1","0.3","1","3","10","30","100")
  ) +
  scale_color_manual(values = pal_country) +
  scale_shape_manual(values = shape_country) +
  coord_flip() +
  labs(
    x = "Predictor variables",
    y = "Odds Ratio (95% CI)",
    color = "Best model",
    shape = "Best model"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10)
  )

print(p)

## Odd ratios original values

library(dplyr)
library(tidyr)

# 1) Load the UN-SCALED data (or whatever file has original units)
df_orig <- read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/01_raw_inputs/Biomphalaria_upd_complete_oct26.csv")

df_orig <- df_orig %>%
  mutate(site_uid = paste(country, site_id_renum, sep = "_"))

# Proportion - median 

median_or_na <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else median(x) }

df_median <- df_orig %>%
  group_by(site_uid) %>%
  summarise(
    prop_water_med = median_or_na(prop_water),
    prop_trees_med = median_or_na(prop_trees),
    prop_crops_med = median_or_na(prop_crops),
    prop_built_med = median_or_na(prop_built),
    .groups = "drop"
  )

df_orig = left_join(df_orig, df_median, by = 'site_uid')

# Matching labels
df_orig <- df_orig %>%
  mutate(Country = dplyr::recode(country, UG = "Uganda", DRC = "DRC", CHAD = "Chad"))

# 2) SD per term per country (original units)
terms_needed <- unique(or_all_aligned$term) %>%
  sort(decreasing=T)

# One table containing:
# - SD in original units
# - SD on the scale actually used in the model
# - country-specific median in original units

# ============================================================
# Reference values on ORIGINAL scale for transformed variables
# ============================================================

ref_long <- df_orig %>%
  group_by(Country) %>%
  summarise(
    pp_oneweek = median(pp_oneweek, na.rm = TRUE),
    flow_acc   = median(flow_acc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    -Country,
    names_to = "term",
    values_to = "reference_value"
  )

# --------------------------------------------------
# SD per term per country
# pp_oneweek and flow_acc were log1p transformed
# before scaling in the models
# --------------------------------------------------

df_for_sd <- df_orig %>%
  mutate(
    pp_oneweek = log1p(pp_oneweek),
    flow_acc   = log1p(flow_acc)
  )

sd_long <- df_for_sd %>%
  group_by(Country) %>%
  summarise(
    across(
      all_of(terms_needed),
      ~ sd(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    -Country,
    names_to = "term",
    values_to = "sd_orig"
  )

sd_long

# 3) Choose meaningful increments Δ (same across countries, per predictor)
delta_lut <- c(
  temp_modis       = 1,     #  °C
  pp_oneweek       = 1,    #  mm
  NDVI_q90_t       = 0.01,   # 0.1 NDVI
  NDVI_modis       = 0.01,
  prop_water_med   = 0.01,   # percentage points if 0-1
  prop_trees_med   = 0.01,
  prop_crops_med   = 0.01,
  prop_built_med   = 0.01,
  pop_rate_change  = 1,     # 1 unit (define what 1 means here)
  flow_acc         = 10000,  # pick a sensible unit
  SS               = 0.1,
  TPI              = 1,
  sampling_experience = 2,
  time_sin         = 1,
  time_cos         = 1
)

#Back-transform: current ORs are per 1 SD.
or_plot <- or_all_aligned %>%
  left_join(sd_long, by = c("Country", "term")) %>%
  mutate(
    delta = ifelse(term %in% names(delta_lut), delta_lut[term], 1),
    
    # convert existing ORs back to beta (log-odds per 1 SD)
    beta_sd     = log(odds_ratio),
    beta_low_sd = log(lower_ci),
    beta_high_sd= log(upper_ci),
    
    # OR for Δ original units
    odds_ratio = exp(beta_sd     * (delta / sd_orig)),
    lower_ci   = exp(beta_low_sd * (delta / sd_orig)),
    upper_ci   = exp(beta_high_sd* (delta / sd_orig))
  )


# Optional: label terms with increment info for clarity
or_plot <- or_plot %>% mutate(term_label = paste0(term_label, " (per ", delta, ")"))

or_plot <- or_all_aligned %>%
  left_join(sd_long, by = c("Country", "term")) %>%
  left_join(ref_long, by = c("Country", "term")) %>%
  mutate(
    
    delta = ifelse(
      term %in% names(delta_lut),
      delta_lut[term],
      1
    ),
    
    # Recover coefficient on standardized model scale
    beta_sd      = log(odds_ratio),
    beta_low_sd  = log(lower_ci),
    beta_high_sd = log(upper_ci),
    
    # Convert requested ORIGINAL-unit increment
    # into an increment on the standardized MODEL scale
    delta_z = case_when(
      
      term %in% c("pp_oneweek", "flow_acc") ~
        (
          log1p(reference_value + delta) -
            log1p(reference_value)
        ) / sd_orig,
      
      TRUE ~ delta / sd_orig
    ),
    
    # OR for requested original-unit increment
    odds_ratio = exp(beta_sd      * delta_z),
    lower_ci   = exp(beta_low_sd  * delta_z),
    upper_ci   = exp(beta_high_sd * delta_z)
  )

or_plot <- or_plot %>% mutate(term_label = paste0(term_label, " (per ", delta, ")"))


p4b <- ggplot(or_plot, aes(x = term_label, y = odds_ratio,
                           color = Country, shape = Country)) +
  geom_pointrange(
    aes(ymin = lower_ci, ymax = upper_ci),
    position = pos,
    linewidth = 0.9
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  scale_y_log10(
    breaks = c(0.1, 0.3, 1, 3, 10, 30, 100),
    labels = c("0.1","0.3","1","3","10","30","100")
  ) +
  scale_color_manual(values = pal_country) +
  scale_shape_manual(values = shape_country) +
  coord_flip() +
  labs(
    x = "Predictor variables",
    y = "Odds Ratio (95% CI)",
    color = "Best model",
    shape = "Best model"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8)
  )



print(p4b)

## figure 4 combined
## Warning appears because some models don't have those predictors!

library(grid)
library(patchwork)
library(svglite)

theme_fig4 <- theme_classic(base_size = 9) +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(size = 9, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    legend.title = element_text(size = 9, face = "bold", colour = "black"),
    legend.text = element_text(size = 8, colour = "black"),
    legend.key = element_blank(),
    legend.background = element_blank(),
    panel.grid.minor = element_blank(),
    plot.tag = element_text(size = 10, face = "plain", colour = "black")
  )

pos_bar <- position_dodge2(width = 0.80, preserve = "single")

p4a <- ggplot(plot_df, aes(x = Variable_label)) +
  geom_col(
    aes(y = Percent, fill = Country, group = Country),
    position = pos_bar,
    width = 0.72,
    color = "grey20",
    linewidth = 0.22
  ) +
  geom_col_pattern(
    aes(
      y = SigPercent,
      group = Country,
      colour = Country,
      pattern = sig_flag
    ),
    position = pos_bar,
    width = 0.72,
    fill = "transparent",
    pattern = "stripe",
    pattern_angle = 45,
    pattern_alpha = 0.30,
    pattern_density = 0.2,
    pattern_size = 0.22,
    pattern_spacing = 0.03,
    pattern_colour = "grey85",
    linewidth = 0.22,
    show.legend = FALSE
  ) +
  scale_pattern_manual(values = c(`FALSE` = "none", `TRUE` = "stripe")) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Percentage of AIC<2 models (%)"
  ) +
  theme_fig4 +
  theme(panel.grid.major.x = element_line(color = "grey85", linewidth = 0.35),
        legend.position = "none",
        plot.margin = margin(5.5, 4, 5.5, 5.5)
  )

pos_or <- position_dodge(width = 0.55)

p4b <- ggplot(or_plot, aes(x = term_label, y = odds_ratio,
                           colour = Country, shape = Country)) +
  geom_linerange(
    aes(ymin = lower_ci, ymax = upper_ci),
    position = pos_or,
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  geom_point(
    position = pos_or,
    size = 2.2,
    stroke = 0.2,
    na.rm = TRUE
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "#8B0000"
  ) +
  scale_y_log10(
    breaks = c(0.1, 0.3, 1, 3, 10, 30, 100),
    labels = c("0.1", "0.3", "1", "3", "10", "30", "100")
  ) +
  scale_colour_manual(values = pal_country) +
  scale_shape_manual(values = shape_country) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Odds ratio (95% CI)",
    colour = "Country",
    shape = "Country"
  ) +
  guides(
    colour = guide_legend(nrow = 1, override.aes = list(size = 2.5)),
    shape  = guide_legend(nrow = 1)
  ) +
  theme_fig4 +
  theme(
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.25),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(0.15, "cm"),
    plot.margin = margin(5.5, 5.5, 5.5, 4)
  )

fig4 <- p4a + p4b +
  plot_layout(widths = c(1.00, 1.08)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag.position = c(0.01, 0.99)
  )

fig4

ggsave(
  filename = "Figure4_complete.svg",
  plot = fig4,
  width = 18,
  height = 13,
  units = "cm",
  device = svglite::svglite,
  bg = "transparent"
)

## Test everything with the best model 

library(sf)
library(dplyr)
library(tidyr)
library(spdep)
library(DHARMa)
library(ggplot2)
library(patchwork)
library(scales)
library(grid)
library(svglite)

# --------------------------------------------------
# 1) Helpers
# --------------------------------------------------

country_lut <- c(UG = "Uganda", DRC = "DRC", CHAD = "Chad")
country_levels <- c("Chad", "DRC", "Uganda")

utm_epsg_from_lonlat <- function(lon, lat) {
  lon0 <- median(lon, na.rm = TRUE)
  lat0 <- median(lat, na.rm = TRUE)
  zone <- floor((lon0 + 180) / 6) + 1
  if (lat0 >= 0) 32600 + zone else 32700 + zone
}

prep_country <- function(df_scaled, cc) {
  df_c <- df_scaled %>%
    filter(country == cc) %>%
    mutate(eventDate = as.POSIXct(eventDate, tz = "UTC"))
  
  epsg <- utm_epsg_from_lonlat(df_c$lon, df_c$lat)
  
  pts <- st_as_sf(df_c, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  pts_utm <- st_transform(pts, epsg)
  xy <- st_coordinates(pts_utm)
  
  pts_utm %>%
    mutate(
      x_utm = xy[, 1],
      y_utm = xy[, 2]
    ) %>%
    st_drop_geometry()
}

to_z <- function(u) qnorm(pmin(pmax(u, 1e-6), 1 - 1e-6))

attach_dharma <- function(dat, mod, country_label) {
  r <- simulateResiduals(mod, plot = FALSE)
  
  dat %>%
    mutate(
      Country = country_label,
      date    = as.Date(eventDate),
      resi    = r$scaledResiduals - 0.5,
      resid_z = to_z(r$scaledResiduals)
    )
}

moran_by_week_distance <- function(df,
                                   dist_vec = 1:10,
                                   min_points = 5,
                                   date_col = "eventDate",
                                   x_m = "x_utm",
                                   y_m = "y_utm",
                                   res_col = "resi") {
  
  df <- df %>%
    mutate(
      week_relative = as.integer(difftime(.data[[date_col]],
                                          min(.data[[date_col]], na.rm = TRUE),
                                          units = "weeks")),
      x_km = .data[[x_m]] / 1000,
      y_km = .data[[y_m]] / 1000
    )
  
  weeks <- sort(unique(df$week_relative))
  
  out <- expand_grid(distance = dist_vec, week = weeks) %>%
    mutate(
      moran_I = NA_real_,
      p_value = NA_real_,
      npoints = NA_integer_,
      n_links = NA_integer_
    )
  
  for (dist in dist_vec) {
    for (wk in weeks) {
      subset_df <- df %>% filter(week_relative == wk)
      npts <- nrow(subset_df)
      
      idx <- which(out$distance == dist & out$week == wk)
      out$npoints[idx] <- npts
      
      if (npts < min_points) next
      
      coords <- cbind(subset_df$x_km, subset_df$y_km)
      nb <- dnearneigh(coords, 0, dist, longlat = FALSE)
      n_links <- sum(card(nb))
      out$n_links[idx] <- n_links
      
      if (n_links == 0) next
      
      lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
      mt <- tryCatch(
        moran.test(subset_df[[res_col]], lw, zero.policy = TRUE),
        error = function(e) NULL
      )
      
      if (is.null(mt)) next
      
      out$moran_I[idx] <- unname(mt$estimate[["Moran I statistic"]])
      out$p_value[idx] <- mt$p.value
    }
  }
  
  out
}

make_lag_pairs <- function(dat, max_lag = 300) {
  dat %>%
    arrange(site_uid, date) %>%
    group_by(site_uid) %>%
    mutate(
      resid_prev = lag(resid_z),
      date_prev  = lag(date),
      gap_days   = as.integer(date - date_prev),
      abs_diff   = abs(resid_z - resid_prev)
    ) %>%
    ungroup() %>%
    filter(!is.na(gap_days), gap_days >= 1, gap_days <= max_lag, !is.na(abs_diff))
}

theme_heat <- theme_classic(base_size = 9) +
  theme(
    text = element_text(colour = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    strip.text = element_text(face = "bold", size = 10, colour = "black"),
    strip.background = element_blank(),
    axis.title = element_text(size = 9, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 9, colour = "black"),
    legend.text  = element_text(size = 8, colour = "black"),
    legend.key = element_blank(),
    panel.spacing.x = unit(4, "mm"),
    panel.spacing.y = unit(2, "mm"),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

# --------------------------------------------------
# 2) Data preparation
# --------------------------------------------------

df_scaled <- read.csv("~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/02_clean_data/df_final_scaled_Bio_upd_Oct26.csv")

country_data <- list(
  UG   = prep_country(df_scaled, "UG"),
  DRC  = prep_country(df_scaled, "DRC"),
  CHAD = prep_country(df_scaled, "CHAD")
)

country_models <- list(
  UG   = UG_best_model,
  DRC  = DRC_best_model,
  CHAD = CHAD_best_model
)

set.seed(123)
country_resid <- lapply(names(country_data), function(cc) {
  attach_dharma(country_data[[cc]], country_models[[cc]], country_lut[[cc]])
})
names(country_resid) <- names(country_data)

# --------------------------------------------------
# 3) Spatial autocorrelation heatmap
# --------------------------------------------------

dist_vec <- 1:10

results <- bind_rows(lapply(names(country_resid), function(cc) {
  moran_by_week_distance(country_resid[[cc]], dist_vec = dist_vec) %>%
    mutate(Country = country_lut[[cc]])
})) %>%
  mutate(
    Country = factor(Country, levels = country_levels),
    week_f = factor(week, levels = sort(unique(week)))
  ) %>%
  group_by(Country) %>%
  mutate(
    q = p.adjust(p_value, method = "BH"),
    sig_fdr = q < 0.05
  ) %>%
  ungroup()

lim_I <- max(abs(results$moran_I), na.rm = TRUE)

week_breaks <- levels(results$week_f)
if (length(week_breaks) > 10) {
  week_breaks <- week_breaks[seq(1, length(week_breaks), by = 4)]
}

p_spatial <- ggplot(results, aes(x = distance, y = week_f, fill = moran_I)) +
  geom_tile(width = 0.98, height = 0.98) +
  geom_point(
    data = subset(results, sig_fdr),
    aes(x = distance, y = week_f),
    inherit.aes = FALSE,
    shape = 15, size = 0.45, colour = "grey15", alpha = 0.55
  ) +
  facet_wrap(~ Country, nrow = 1) +
  scale_x_continuous(
    breaks = c(2, 4, 6, 8, 10),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    breaks = week_breaks,
    expand = c(0, 0)
  ) +
  scale_fill_gradient2(
    low = "#6883A6",
    mid = "white",
    high = "#B35C4A",
    midpoint = 0,
    limits = c(-lim_I, lim_I),
    oob = squish,
    na.value = "grey94",
    name = "Moran's I"
  ) +
  guides(
    fill = guide_colorbar(
      barheight = unit(28, "pt"),
      frame.colour = "grey70",
      ticks.colour = "grey40"
    )
  ) +
  labs(
    x = "Distance threshold (km)",
    y = "Relative week"
  ) +
  theme_heat + 
  theme (legend.position = 'bottom')

# --------------------------------------------------
# 4) Temporal lag heatmap
# --------------------------------------------------

max_lag_1d <- 30

heat_1d <- bind_rows(lapply(names(country_resid), function(cc) {
  make_lag_pairs(country_resid[[cc]], max_lag = 300)
})) %>%
  filter(gap_days <= max_lag_1d) %>%
  group_by(Country, site_uid, lag_day = gap_days) %>%
  summarise(
    mean_absdiff   = mean(abs_diff, na.rm = TRUE),
    median_absdiff = median(abs_diff, na.rm = TRUE),
    n_pairs        = n(),
    .groups = "drop"
  ) %>%
  filter(n_pairs >= 1)

site_rank_1d <- heat_1d %>%
  filter(lag_day <= 14) %>%
  group_by(Country, site_uid) %>%
  summarise(score = median(median_absdiff, na.rm = TRUE), .groups = "drop") %>%
  arrange(Country, score) %>%
  transmute(site_key = paste(Country, site_uid, sep = "||")) %>%
  pull(site_key)

heat_1d <- heat_1d %>%
  mutate(
    Country = factor(Country, levels = country_levels),
    site_key = factor(paste(Country, site_uid, sep = "||"), levels = site_rank_1d)
  )

cap1 <- quantile(heat_1d$mean_absdiff, 0.95, na.rm = TRUE)

heat_1d <- heat_1d %>%
  mutate(
    mean_absdiff_cap = pmin(mean_absdiff, cap1),
    alpha_n = rescale(pmin(n_pairs, 10), to = c(0.45, 1))
  )

p_temporal <- ggplot(heat_1d, aes(x = site_key, y = lag_day, fill = mean_absdiff_cap)) +
  geom_tile(aes(alpha = alpha_n), width = 0.98, height = 0.98) +
  facet_grid(. ~ Country, scales = "free_x", space = "free_x") +
  scale_alpha_identity() +
  scale_y_continuous(
    breaks = c(1, 3, 5, 7, 10, 14, 21, 28, 30),
    limits = c(1, max_lag_1d),
    expand = c(0, 0)
  ) +
  scale_fill_viridis_c(
    option = "C",
    direction = -1,
    limits = c(min(heat_1d$mean_absdiff_cap, na.rm = TRUE), cap1),
    oob = squish,
    name = "Mean |Δ residual|"
  ) +
  guides(
    fill = guide_colorbar(
      barheight = unit(28, "pt"),
      frame.colour = "grey70",
      ticks.colour = "grey40"
    )
  ) +
  labs(
    x = "Sites",
    y = "Gap between visits (days)"
  ) +
  theme_heat +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_text(margin = margin(t = 6)),
    legend.position = 'bottom'
  )

# --------------------------------------------------
# 5) Combined figure
# --------------------------------------------------

fig_auto <- p_temporal + p_spatial +
  plot_layout(widths = c(1, 1.35)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(size = 10, face = "plain", colour = "black"),
    plot.tag.position = c(0.01, 0.99)
  )

fig_auto

ggsave(
  filename = "Figure5_complete.svg",
  plot = fig_auto,
  width = 18,
  height = 13,
  units = "cm",
  device = svglite::svglite,
  bg = "transparent"
)

### Variation for DRC with spatial term 

library(mgcv)
library(dplyr)
library(ggplot2)
library(DHARMa)
library(scales)
library(grid)

# -----------------------------
# 1) Fit DRC model with spatial smooth
# -----------------------------
UG   <- prep_country(df_scaled, "UG")
DRC  <- prep_country(df_scaled, "DRC")
CHAD <- prep_country(df_scaled, "CHAD")

DRC_spdat <- DRC %>%
  dplyr::mutate(
    bio_pres = as.integer(as.character(bio_pres)),
    site_uid = factor(site_uid),
    x_utm = as.numeric(x_utm),
    y_utm = as.numeric(y_utm)
  )

n_xy <- n_distinct(DRC_spdat$x_utm, DRC_spdat$y_utm)
k_xy <- min(50, n_xy - 1)

m_drc_gam <- bam(
  bio_pres ~ pp_oneweek + temp_modis + NDVI_q90_t +
     prop_water_med + prop_trees_med + prop_built_med +
    sampling_experience + time_sin + time_cos + NDVI_modis+
    s(x_utm, y_utm, bs = "tp", k = k_xy) +
    s(site_uid, bs = "re"),
  data = DRC_spdat,
  family = binomial(),
  method = "fREML"
)

DRC_best_model_spatial <- m_drc_gam

# -----------------------------
# 2) Attach DHARMa residuals
# -----------------------------
set.seed(123)

UG_M1   <- attach_dharma(UG,        UG_best_model,           "Uganda")
DRC_M2  <- attach_dharma(DRC_spdat, DRC_best_model_spatial,  "DRC")
CHAD_M1 <- attach_dharma(CHAD,      CHAD_best_model,         "Chad")

# -----------------------------
# 3) Moran's I by week x distance
# -----------------------------
dist_vec <- 1:10

UG_res   <- moran_by_week_distance(UG_M1,   dist_vec = dist_vec)
DRC_res  <- moran_by_week_distance(DRC_M2,  dist_vec = dist_vec)
CHAD_res <- moran_by_week_distance(CHAD_M1, dist_vec = dist_vec)

plot_df_sp <- bind_rows(
  mutate(UG_res,   Country = "Uganda"),
  mutate(DRC_res,  Country = "DRC"),
  mutate(CHAD_res, Country = "Chad")
) %>%
  mutate(
    Country = factor(Country, levels = c("Chad", "DRC", "Uganda")),
    week_f  = factor(week, levels = sort(unique(week)))
  ) %>%
  group_by(Country) %>%
  mutate(
    q = p.adjust(p_value, method = "BH"),
    sig_fdr = q < 0.05
  ) %>%
  ungroup()

# symmetric colour scale
lim_I2 <- max(abs(plot_df_sp$moran_I), na.rm = TRUE)

# fewer week labels, same style as before
week_breaks2 <- levels(plot_df_sp$week_f)
if (length(week_breaks2) > 10) {
  week_breaks2 <- week_breaks2[seq(1, length(week_breaks2), by = 4)]
}

# -----------------------------
# 4) Plot with same style as previous figure
# -----------------------------
p_spatial_DRC <- ggplot(plot_df_sp, aes(x = distance, y = week_f, fill = moran_I)) +
  geom_tile(width = 0.98, height = 0.98) +
  geom_point(
    data = subset(plot_df_sp, sig_fdr),
    aes(x = distance, y = week_f),
    inherit.aes = FALSE,
    shape = 15,
    size = 0.45,
    colour = "grey15",
    alpha = 0.55
  ) +
  facet_wrap(~ Country, nrow = 1) +
  scale_x_continuous(
    breaks = c(2, 4, 6, 8, 10),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    breaks = week_breaks2,
    expand = c(0, 0)
  ) +
  scale_fill_gradient2(
    low = "#6883A6",
    mid = "white",
    high = "#B35C4A",
    midpoint = 0,
    limits = c(-lim_I2, lim_I2),
    oob = squish,
    na.value = "grey94",
    name = "Moran's I"
  ) +
  guides(
    fill = guide_colorbar(
      barheight = unit(28, "pt"),
      frame.colour = "grey70",
      ticks.colour = "grey40"
    )
  ) +
  labs(
    x = "Distance threshold (km)",
    y = "Relative week"
  ) +
  theme_heat +
  theme(
    legend.position = "bottom",
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

p_spatial_DRC

ggsave(
  filename = "Figure_appendix_DRC_spatial_Moran.svg",
  plot = p_spatial_DRC,
  width = 18,
  height = 15,
  units = "cm",
  device = svglite::svglite,
  bg = "transparent"
)

# ============================================================
# Compare DRC coefficients:
# non-spatial GLMM vs spatial GAM
# ============================================================

library(broom.mixed)
library(tibble)

# Non-spatial GLMM
coef_nonspatial <- broom.mixed::tidy(
  DRC_best_model,
  effects = "fixed",
  conf.int = TRUE
) %>%
  filter(term != "(Intercept)") %>%
  transmute(
    term,
    beta_nonspatial = estimate,
    low_nonspatial = conf.low,
    high_nonspatial = conf.high,
    sig_nonspatial = conf.low > 0 | conf.high < 0
  )

# Spatial GAM: parametric/fixed effects
gam_par <- as.data.frame(
  summary(DRC_best_model_spatial)$p.table
) %>%
  rownames_to_column("term")

coef_spatial <- gam_par %>%
  filter(term != "(Intercept)") %>%
  transmute(
    term,
    beta_spatial = Estimate,
    se_spatial = `Std. Error`,
    low_spatial = beta_spatial - 1.96 * se_spatial,
    high_spatial = beta_spatial + 1.96 * se_spatial,
    sig_spatial = low_spatial > 0 | high_spatial < 0
  )

coef_compare <- full_join(
  coef_nonspatial,
  coef_spatial,
  by = "term"
) %>%
  mutate(
    absolute_change =
      beta_spatial - beta_nonspatial,
    
    percent_change =
      100 *
      abs(beta_spatial - beta_nonspatial) /
      abs(beta_nonspatial),
    
    same_direction =
      sign(beta_nonspatial) ==
      sign(beta_spatial),
    
    significance_same =
      sig_nonspatial ==
      sig_spatial
  )

coef_compare

# Proper manuscript labels
term_labels <- c(
  "pp_oneweek"          = "Precipitation (1-week)",
  "temp_modis"          = "Temperature (MODIS)",
  "NDVI_q90_t"          = "NDVI (Planet)",
  "prop_water_med"      = "Proportion of water",
  "prop_trees_med"      = "Proportion of trees",
  "prop_built_med"      = "Proportion of built area",
  "sampling_experience" = "Sampling experience",
  "time_sin"            = "Sampling time (sin)",
  "time_cos"            = "Sampling time (cos)",
  "NDVI_modis"          = "NDVI (MODIS)"
)

coef_long <- coef_compare %>%
  mutate(
    term_label = dplyr::recode(
      term,
      !!!term_labels
    )
  ) %>%
  dplyr::select(
    term_label,
    beta_nonspatial,
    beta_spatial
  ) %>%
  pivot_longer(
    cols = c(
      beta_nonspatial,
      beta_spatial
    ),
    names_to = "Model",
    values_to = "beta"
  ) %>%
  mutate(
    Model = dplyr::recode(
      Model,
      beta_nonspatial = "Non-spatial GLMM",
      beta_spatial = "Spatial sensitivity"
    )
  )

fig_DRC_spatial_rev = ggplot(
  coef_long,
  aes(
    x = beta,
    y = term_label,
    colour = Model
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_point(
    size = 2.5,
    position = position_dodge(width = 0.35)
  ) +
  theme_classic() +
  labs(
    x = "Standardised coefficient",
    y = NULL,
    colour = NULL
  )

ggsave(
  filename = "FigureDRCspatial_reviewer.svg",
  plot = fig_DRC_spatial_rev,
  width = 18,
  height = 13,
  units = "cm",
  device = svglite::svglite,
  bg = "transparent"
)

#### Sampling experience and time-corrected experience

library(dplyr)
library(lme4)
library(broom.mixed)
library(ggplot2)

df_all <- df_scaled

# If jul_dat is not already present in df_scaled, create it
if (!"jul_dat" %in% names(df_all)) {
  
  df_all$jul_dat <- julian.Date(
    as.Date(df_all$eventDate)
  )
  
  # Standardise within country, not across all countries
  df_all <- df_all %>%
    group_by(country) %>%
    mutate(
      jul_dat = (
        jul_dat - mean(jul_dat, na.rm = TRUE)
      ) / sd(jul_dat, na.rm = TRUE)
    ) %>%
    ungroup()
}

# --------------------------------------------------
# 2. Residualise sampling experience against julian date
# --------------------------------------------------

df_all <- df_all %>%
  group_by(country) %>%
  group_modify(
    ~ {
      
      m_exp <- lm(
        sampling_experience ~ jul_dat,
        data = .x,
        na.action = na.exclude
      )
      
      .x %>%
        mutate(
          exp_resid = residuals(m_exp)
        )
    }
  ) %>%
  ungroup()

df_all %>%
  group_by(country) %>%
  summarise(
    r_raw = cor(
      sampling_experience,
      jul_dat,
      use = "complete.obs"
    ),
    r_resid = cor(
      exp_resid,
      jul_dat,
      use = "complete.obs"
    ),
    .groups = "drop"
  )

ggplot(
  df_all,
  aes(
    x = jul_dat,
    y = exp_resid
  )
) +
  geom_point(
    alpha = 0.35,
    size = 1
  ) +
  facet_wrap(
    ~ country,
    scales = "free_x"
  ) +
  labs(
    x = "Julian date",
    y = "Time-adjusted sampling experience"
  ) +
  theme_bw()

# --------------------------------------------------
# 3. Replace sampling experience by exp_resid
#    and add jul_dat
# --------------------------------------------------

make_time_adjusted_formula <- function(formula_string) {
  
  f <- as.formula(formula_string)
  
  f2 <- update(
    f,
    . ~ . -
      sampling_experience +
      exp_resid +
      jul_dat
  )
  
  paste(
    deparse(f2),
    collapse = " "
  )
}

UG_f_with   <- UG_best_row$Formula[1]
DRC_f_with  <- DRC_best_row$Formula[1]
CHAD_f_with <- CHAD_best_row$Formula[1]

UG_f_all <- make_time_adjusted_formula(
  UG_f_with
)

DRC_f_all <- make_time_adjusted_formula(
  DRC_f_with
)

CHAD_f_all <- make_time_adjusted_formula(
  CHAD_f_with
)

UG_f_all
DRC_f_all
CHAD_f_all

# --------------------------------------------------
# 4. Country datasets
# --------------------------------------------------

UG_data   <- df_all %>% filter(country == "UG")
DRC_data  <- df_all %>% filter(country == "DRC")
CHAD_data <- df_all %>% filter(country == "CHAD")


# --------------------------------------------------
# 5. Fit models
# --------------------------------------------------

fit_from_string <- function(
    formula_string,
    data
) {
  
  glmer(
    as.formula(formula_string),
    data = data,
    family = binomial(),
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 2e5
      )
    )
  )
}

UG_mod_all <- fit_from_string(
  UG_f_all,
  UG_data
)

DRC_mod_all <- fit_from_string(
  DRC_f_all,
  DRC_data
)

CHAD_mod_all <- fit_from_string(
  CHAD_f_all,
  CHAD_data
)


# --------------------------------------------------
# 6. Extract coefficients
# --------------------------------------------------

tidy_beta <- function(
    mod,
    country_label
) {
  
  broom.mixed::tidy(
    mod,
    effects = "fixed",
    conf.int = TRUE
  ) %>%
    filter(
      term != "(Intercept)"
    ) %>%
    mutate(
      Country = country_label
    ) %>%
    dplyr::select(
      Country,
      term,
      estimate,
      conf.low,
      conf.high,
      p.value
    )
}

beta_all <- bind_rows(
  tidy_beta(
    UG_mod_all,
    "Uganda"
  ),
  tidy_beta(
    DRC_mod_all,
    "DRC"
  ),
  tidy_beta(
    CHAD_mod_all,
    "Chad"
  )
)

beta_time_exp <- beta_all %>%
  filter(
    term %in% c(
      "exp_resid",
      "jul_dat"
    )
  ) %>%
  mutate(
    OR = exp(estimate),
    OR_low = exp(conf.low),
    OR_high = exp(conf.high)
  )

beta_time_exp

# ============================================================
# ORIGINAL-SCALE INFORMATION FOR TIME-ADJUSTED MODELS
# ============================================================

library(dplyr)
library(tidyr)

# Raw data
df_orig_exp <- read.csv(
  "~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/01_raw_inputs/Biomphalaria_upd_complete_oct26.csv"
)

df_orig_exp <- df_orig_exp %>%
  mutate(
    site_uid = paste(
      country,
      site_id_renum,
      sep = "_"
    ),
    Country = dplyr::recode(
      country,
      UG = "Uganda",
      DRC = "DRC",
      CHAD = "Chad"
    ),
    
    # Raw Julian date in DAYS
    jul_dat_raw = as.numeric(
      julian.Date(
        as.Date(eventDate)
      )
    )
  )

median_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else median(x)
}

df_median_exp <- df_orig_exp %>%
  group_by(site_uid) %>%
  summarise(
    prop_water_med = median_or_na(prop_water),
    prop_trees_med = median_or_na(prop_trees),
    prop_crops_med = median_or_na(prop_crops),
    prop_built_med = median_or_na(prop_built),
    .groups = "drop"
  )

df_orig_exp <- df_orig_exp %>%
  left_join(
    df_median_exp,
    by = "site_uid"
  )

# ============================================================
# SDs corresponding to the scale entering the models
# ============================================================

terms_for_sd <- c(
  "temp_modis",
  "SS",
  "prop_water_med",
  "prop_trees_med",
  "prop_built_med",
  "pp_oneweek",
  "NDVI_q90_t",
  "NDVI_modis",
  "flow_acc",
  "elev"
)

df_for_sd_exp <- df_orig_exp %>%
  mutate(
    pp_oneweek = log1p(pp_oneweek),
    flow_acc   = log1p(flow_acc)
  )

sd_long_exp <- df_for_sd_exp %>%
  group_by(Country) %>%
  summarise(
    across(
      all_of(terms_for_sd),
      ~ sd(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    -Country,
    names_to = "term",
    values_to = "sd_model"
  )

ref_log_exp <- df_orig_exp %>%
  group_by(Country) %>%
  summarise(
    pp_oneweek = median(
      pp_oneweek,
      na.rm = TRUE
    ),
    flow_acc = median(
      flow_acc,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(
    -Country,
    names_to = "term",
    values_to = "reference_value"
  )

julian_sd <- df_orig_exp %>%
  group_by(Country) %>%
  summarise(
    julian_sd_days = sd(
      jul_dat_raw,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

delta_lut_exp <- c(
  temp_modis       = 1,       # +1 °C
  SS               = 0.1,
  prop_water_med   = 0.01,    # +1 percentage point
  prop_trees_med   = 0.01,
  prop_built_med   = 0.01,
  pp_oneweek       = 1,       # +1 mm
  NDVI_q90_t       = 0.05,
  NDVI_modis       = 0.05,
  flow_acc         = 10000,
  elev             = 50
)

# ============================================================
# CONVERT TIME-ADJUSTED MODEL COEFFICIENTS TO ORs
# ============================================================

beta_plot <- beta_all %>%
  left_join(
    sd_long_exp,
    by = c("Country", "term")
  ) %>%
  left_join(
    ref_log_exp,
    by = c("Country", "term")
  ) %>%
  left_join(
    julian_sd,
    by = "Country"
  ) %>%
  mutate(
    
    # --------------------------------------------------------
    # Change expressed on the model's standardised scale
    # --------------------------------------------------------
    
    delta_z = case_when(
      
      # Precipitation and flow accumulation:
      # log1p-transformed before scaling
      term %in% c(
        "pp_oneweek",
        "flow_acc"
      ) ~
        (
          log1p(
            reference_value +
              delta_lut_exp[term]
          ) -
            log1p(reference_value)
        ) / sd_model,
      
      # Julian date was standardised within country.
      # Convert 30 calendar days to SD units.
      term == "jul_dat" ~
        30 / julian_sd_days,
      
      # Residualised sampling experience is already
      # expressed on the standardised experience scale
      term == "exp_resid" ~
        0.1,
      
      # Cyclic terms: report a +0.2 change
      term %in% c(
        "time_sin",
        "time_cos"
      ) ~
        0.2,
      
      # Ordinary predictors
      term %in% names(delta_lut_exp) ~
        delta_lut_exp[term] /
        sd_model,
      
      TRUE ~ NA_real_
    ),
    
    OR = exp(
      estimate * delta_z
    ),
    
    OR_low = exp(
      conf.low * delta_z
    ),
    
    OR_high = exp(
      conf.high * delta_z
    )
  )

beta_plot <- beta_plot %>%
  mutate(
    term_label = case_when(
      
      term == "temp_modis" ~
        "Temperature (+1 °C)",
      
      term == "SS" ~
        "Slope steepness factor (+0.1)",
      
      term == "time_sin" ~
        "Sampling time (sin) (+0.2)",
      
      term == "time_cos" ~
        "Sampling time (cos) (+0.2)",
      
      term == "prop_water_med" ~
        "Proportion of water (+1 pp)",
      
      term == "prop_trees_med" ~
        "Proportion of trees (+1 pp)",
      
      term == "prop_built_med" ~
        "Proportion of built area (+1 pp)",
      
      term == "pp_oneweek" ~
        "Precipitation (1-week) (+1 mm)",
      
      term == "NDVI_q90_t" ~
        "NDVI (Planet) (+0.05)",
      
      term == "NDVI_modis" ~
        "NDVI (MODIS) (+0.05)",
      
      term == "jul_dat" ~
        "Julian date (+30 days)",
      
      term == "flow_acc" ~
        "Flow accumulation (+10,000)",
      
      term == "exp_resid" ~
        "Time-adjusted sampling experience (+0.1)",
      
      term == "elev" ~
        "Elevation (+50 m)",
      
      TRUE ~ term
    )
  )

term_order_exp <- c(
  "Temperature (+1 °C)",
  "Slope steepness factor (+0.1)",
  "Sampling time (sin) (+0.2)",
  "Sampling time (cos) (+0.2)",
  "Proportion of water (+1 pp)",
  "Proportion of trees (+1 pp)",
  "Proportion of built area (+1 pp)",
  "Precipitation (1-week) (+1 mm)",
  "NDVI (Planet) (+0.05)",
  "NDVI (MODIS) (+0.05)",
  "Julian date (+30 days)",
  "Flow accumulation (+10,000)",
  "Time-adjusted sampling experience (+0.1)",
  "Elevation (+50 m)"
)

beta_plot <- beta_plot %>%
  mutate(
    Country = factor(
      Country,
      levels = c(
        "Uganda",
        "DRC",
        "Chad"
      )
    ),
    term_label = factor(
      term_label,
      levels = rev(term_order_exp)
    )
  )

pal_country <- c(
  "Uganda" = "#111A26",
  "DRC"    = "#3F618C",
  "Chad"   = "#6883A6"
)

shape_country <- c(
  "Uganda" = 15,
  "DRC" = 17,
  "Chad" = 16
)

pos_exp <- position_dodge(
  width = 0.55
)

p_exp_adjusted <- ggplot(
  beta_plot,
  aes(
    x = term_label,
    y = OR,
    colour = Country,
    shape = Country
  )
) +
  geom_linerange(
    aes(
      ymin = OR_low,
      ymax = OR_high
    ),
    position = pos_exp,
    linewidth = 0.75,
    na.rm = TRUE
  ) +
  geom_point(
    position = pos_exp,
    size = 2.2,
    na.rm = TRUE
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "#8B0000",
    linewidth = 0.45
  ) +
  scale_y_log10() +
  scale_colour_manual(
    values = pal_country
  ) +
  scale_shape_manual(
    values = shape_country
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Odds ratio (95% CI)",
    colour = "Country",
    shape = "Country"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    panel.grid.major.x = element_line(
      colour = "grey85",
      linewidth = 0.35
    ),
    panel.grid.major.y = element_line(
      colour = "grey90",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.y = element_text(
      size = 8
    )
  )

p_exp_adjusted

ggsave(
  filename = "FigA4.svg",
  plot = p_exp_adjusted,
  width = 18,
  height = 13,
  units = "cm",
  device = svglite::svglite,
  bg = "transparent"
)




######################
## Extra figures


##### 

library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(scales)
library(patchwork)

# ============================================================
# 1. CURRENT Chad snail observations
# ============================================================

df_chad <- read.csv(
  "~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/01_raw_inputs/Biomphalaria_upd_complete_oct26.csv"
) %>%
  filter(country == "CHAD") %>%
  mutate(
    date = as.Date(eventDate),
    week = floor_date(date, "week", week_start = 1)
  )

chad_week_prev <- df_chad %>%
  group_by(week) %>%
  summarise(
    prev = mean(bio_pres, na.rm = TRUE),
    n_surveys = n(),
    .groups = "drop"
  )

# ============================================================
# 2. DAILY CHIRPS RAINFALL
# ============================================================

chirps_raw <- read.csv(
  "~/Documents/Documents/Chad/Chad_files/CHIRPS_nearest_CH.csv"
)

# Check names once
names(chirps_raw)

rain_cols <- grep(
  "^X[0-9]{8}_precipitation$",
  names(chirps_raw),
  value = TRUE
)

chirps_long <- chirps_raw %>%
  pivot_longer(
    cols = all_of(rain_cols),
    names_to = "date_col",
    values_to = "pp"
  ) %>%
  mutate(
    date = as.Date(
      str_extract(date_col, "[0-9]{8}"),
      format = "%Y%m%d"
    ),
    week = floor_date(
      date,
      "week",
      week_start = 1
    )
  )

chad_week_rain <- chirps_long %>%
  group_by(week) %>%
  summarise(
    precipitation = median(
      pp,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

theme_ts <- theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_prev <- ggplot(
  chad_week_prev,
  aes(
    x = week,
    y = prev
  )
) +
  geom_line(
    linewidth = 0.8,
    colour = "#111A26"
  ) +
  geom_point(
    size = 1.8,
    colour = "#111A26"
  ) +
  scale_y_continuous(
    labels = percent_format(
      accuracy = 1
    )
  ) +
  labs(
    y = "Snail prevalence"
  ) +
  theme_ts


p_rain <- ggplot(
  chad_week_rain,
  aes(
    x = week,
    y = precipitation
  )
) +
  geom_line(
    linewidth = 0.8,
    colour = "#3F618C"
  ) +
  labs(
    y = "Weekly precipitation (mm)"
  ) +
  theme_ts


p_effort <- ggplot(
  chad_week_prev,
  aes(
    x = week,
    y = n_surveys
  )
) +
  geom_col(
    fill = "#6883A6",
    width = 6
  ) +
  labs(
    x = "Week",
    y = "Number of reports"
  ) +
  theme_ts


p_full <- (
  p_prev /
    p_rain /
    p_effort
) +
  plot_layout(
    heights = c(
      1,
      1,
      0.9
    )
  )

p_full


ggsave(
  filename = "FigA5.svg",
  plot = p_full,
  width = 18,
  height = 13,
  units = "cm",
  device = svglite::svglite,
  bg = "transparent"
)




########
#### Reports per month 
library(dplyr)
library(lubridate)
library(ggplot2)

reports_month <- df_all %>%
  filter(country %in% c("UG", "DRC")) %>%
  dplyr::mutate(
    eventDate = as.Date(eventDate),
    month_date = floor_date(eventDate, "month"),
    country = dplyr::recode(country, 'UG' = "Uganda", 'DRC' = "DRC")
  ) %>%
  dplyr::count(country, month_date, name = "n_reports") %>%
  arrange(country, month_date) %>%
  group_by(country) %>%
  mutate(time_index = row_number()) %>%
  ungroup()

reportspermonth = ggplot(reports_month, aes(x = month_date, y = n_reports, colour = country)) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_point(size = 2.2, alpha = 0.95) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.9,
    linetype = "22"
  ) +
  scale_x_date(
    date_labels = "%b %Y",
    date_breaks = "3 months",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_colour_manual(
    values = c("DRC" = "#3F618C", "Uganda" = "#111A26")
  ) +
  labs(
    x = "Month",
    y = "Number of reports",
    colour = "Country"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 15),
    panel.border = element_blank()
  )


ggsave(
  filename = "FigA6.svg",
  plot = reportspermonth,
  width = 18,
  height = 13,
  units = "cm",
  device = svglite::svglite,
  bg = "transparent"
)


## Correlation between sampling_experience and Julian date

library(dplyr)
library(lubridate)
library(ggplot2)

pal_country <- c(
  "Uganda" = "#111A26",
  "DRC"    = "#3F618C",
  "Chad"   = "#6883A6"
)


df <- df_all %>%
  mutate(
    eventDate = as.Date(eventDate),
    jul_dat = julian.Date(as.Date(df_all$eventDate)),
    country = case_when(
      country %in% c("UG", "Uganda")   ~ "Uganda",
      country %in% c("DRC", "Congo")   ~ "DRC",
      country %in% c("CHAD", "Chad")   ~ "Chad",
      TRUE ~ as.character(country)
    )
  ) %>%
  filter(country %in% c("Uganda", "DRC", "Chad")) %>%
  filter(!is.na(jul_dat), !is.na(sampling_experience))

cor_stats <- df %>%
  group_by(country) %>%
  summarise(
    r = cor(jul_dat, sampling_experience, method = "pearson"),
    p_value = cor.test(jul_dat, sampling_experience, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      "r = ", round(r, 2),
      "\np = ", signif(p_value, 2)
    )
  )

cor_stats

label_pos <- df %>%
  group_by(country) %>%
  summarise(
    x = min(jul_dat, na.rm = TRUE) + 5,
    y = max(sampling_experience, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(cor_stats, by = "country")

cor_exp_jul = ggplot(df, aes(x = jul_dat, y = sampling_experience, colour = country)) +
  geom_point(alpha = 0.45, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5) +
  facet_wrap(~ country, scales = "free") +
  geom_text(
    data = label_pos,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    colour = "black"
  ) +
  scale_colour_manual(values = pal_country) +
  labs(
    x = "Julian day",
    y = "Sampling experience"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

cor_exp_jul <- cor_exp_jul &
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7.5),
    strip.text = element_text(size = 8.5, face = "bold"),
    plot.title = element_text(size = 8.5, face = "bold"),
    plot.margin = margin(3, 3, 3, 3)
  )

ggsave(
  filename = "FigA3.svg",
  plot = cor_exp_jul,
  width = 18,
  height = 13,
  units = "cm",
  device = svglite::svglite,
  bg = "transparent"
)

# ============================================================
# VALIDATION
# CURRENT Oct26 models and data
# ============================================================

library(dplyr)
library(purrr)
library(lme4)
library(ggplot2)
library(patchwork)
library(pROC)
library(tidyr)

# df_scaled
# UG_best_row
# DRC_best_row
# CHAD_best_row
# already exist


# ============================================================
# 1. METRICS
# ============================================================

brier <- function(y, p) {
  mean((y - p)^2, na.rm = TRUE)
}

auc_roc <- function(y, p) {
  
  ok <- !is.na(y) & !is.na(p)
  
  y <- y[ok]
  p <- p[ok]
  
  if (length(unique(y)) < 2) {
    return(NA_real_)
  }
  
  as.numeric(
    pROC::auc(
      pROC::roc(
        y,
        p,
        quiet = TRUE,
        direction = "<"
      )
    )
  )
}


# ============================================================
# 2. CURRENT COUNTRY DATA + FORMULAS
# ============================================================

UG_data <- df_scaled %>%
  filter(country == "UG")

DRC_data <- df_scaled %>%
  filter(country == "DRC")

CHAD_data <- df_scaled %>%
  filter(country == "CHAD")


UG_f <- UG_best_row$Formula[1]
DRC_f <- DRC_best_row$Formula[1]
CHAD_f <- CHAD_best_row$Formula[1]

# ============================================================
# 3. KNOWN-SITE CV
# Each site contributes observations to train and test
# ============================================================

make_within_site_folds <- function(
    df,
    site_col = "site_uid",
    k = 5,
    seed = 1
) {
  
  set.seed(seed)
  
  df %>%
    group_by(
      .data[[site_col]]
    ) %>%
    mutate(
      fold = sample(
        rep(
          1:k,
          length.out = n()
        )
      )
    ) %>%
    ungroup()
}


cv_known_site <- function(
    df,
    formula_str,
    k = 5,
    seed = 1,
    site_col = "site_uid",
    y_col = "bio_pres"
) {
  
  df_f <- make_within_site_folds(
    df,
    site_col = site_col,
    k = k,
    seed = seed
  )
  
  fold_results <- map(
    1:k,
    function(f) {
      
      train <- df_f %>%
        filter(fold != f)
      
      test <- df_f %>%
        filter(fold == f)
      
      mod <- glmer(
        as.formula(formula_str),
        data = train,
        family = binomial(),
        control = glmerControl(
          optimizer = "bobyqa",
          optCtrl = list(
            maxfun = 2e5
          )
        )
      )
      
      # Known sites:
      # use site random intercepts estimated
      # from training observations
      p <- predict(
        mod,
        newdata = test,
        type = "response"
      )
      
      # Prevalence-only benchmark
      p_null <- rep(
        mean(
          train[[y_col]],
          na.rm = TRUE
        ),
        nrow(test)
      )
      
      pred_df <- tibble(
        fold = f,
        site_uid = test[[site_col]],
        y = test[[y_col]],
        p = p,
        p_null = p_null
      )
      
      B_model <- brier(
        pred_df$y,
        pred_df$p
      )
      
      B_null <- brier(
        pred_df$y,
        pred_df$p_null
      )
      
      metrics <- tibble(
        fold = f,
        n_test = nrow(test),
        n_sites = n_distinct(
          test[[site_col]]
        ),
        prevalence = mean(
          test[[y_col]],
          na.rm = TRUE
        ),
        Brier = B_model,
        Brier_null = B_null,
        Brier_skill =
          1 - B_model / B_null,
        AUC = auc_roc(
          pred_df$y,
          pred_df$p
        )
      )
      
      list(
        metrics = metrics,
        preds = pred_df
      )
    }
  )
  
  metrics <- bind_rows(
    map(
      fold_results,
      "metrics"
    )
  )
  
  preds <- bind_rows(
    map(
      fold_results,
      "preds"
    )
  )
  
  B_model <- brier(
    preds$y,
    preds$p
  )
  
  B_null <- brier(
    preds$y,
    preds$p_null
  )
  
  pooled <- tibble(
    Brier = B_model,
    Brier_null = B_null,
    Brier_skill =
      1 - B_model / B_null,
    AUC = auc_roc(
      preds$y,
      preds$p
    )
  )
  
  list(
    metrics = metrics,
    preds = preds,
    pooled = pooled
  )
}

# ============================================================
# 4. UNSEEN-SITE CV
# Entire sites are held out
# ============================================================

make_site_folds <- function(
    df,
    site_col = "site_uid",
    y_col = "bio_pres",
    k = 5,
    seed = 1
) {
  
  set.seed(seed)
  
  site_table <- df %>%
    group_by(
      .data[[site_col]]
    ) %>%
    summarise(
      
      # Helps distribute sites containing
      # presences across folds
      any_presence = as.integer(
        any(
          .data[[y_col]] == 1,
          na.rm = TRUE
        )
      ),
      
      n_obs = n(),
      
      .groups = "drop"
    ) %>%
    group_by(any_presence) %>%
    mutate(
      fold = sample(
        rep(
          1:k,
          length.out = n()
        )
      )
    ) %>%
    ungroup()
  
  df %>%
    left_join(
      site_table %>%
        dplyr::select(
          all_of(site_col),
          fold
        ),
      by = site_col
    )
}


cv_unseen_site <- function(
    df,
    formula_str,
    k = 5,
    seed = 1,
    site_col = "site_uid",
    y_col = "bio_pres"
) {
  
  df_f <- make_site_folds(
    df,
    site_col = site_col,
    y_col = y_col,
    k = k,
    seed = seed
  )
  
  fold_results <- map(
    1:k,
    function(f) {
      
      train <- df_f %>%
        filter(fold != f)
      
      test <- df_f %>%
        filter(fold == f)
      
      mod <- glmer(
        as.formula(formula_str),
        data = train,
        family = binomial(),
        control = glmerControl(
          optimizer = "bobyqa",
          optCtrl = list(
            maxfun = 2e5
          )
        )
      )
      
      # Held-out sites have no estimated random intercept.
      # Predict from fixed effects only.
      p <- predict(
        mod,
        newdata = test,
        type = "response",
        re.form = NA,
        allow.new.levels = TRUE
      )
      
      p_null <- rep(
        mean(
          train[[y_col]],
          na.rm = TRUE
        ),
        nrow(test)
      )
      
      pred_df <- tibble(
        fold = f,
        site_uid = test[[site_col]],
        y = test[[y_col]],
        p = p,
        p_null = p_null
      )
      
      B_model <- brier(
        pred_df$y,
        pred_df$p
      )
      
      B_null <- brier(
        pred_df$y,
        pred_df$p_null
      )
      
      metrics <- tibble(
        fold = f,
        n_test = nrow(test),
        n_sites = n_distinct(
          test[[site_col]]
        ),
        n_presence = sum(
          test[[y_col]] == 1,
          na.rm = TRUE
        ),
        prevalence = mean(
          test[[y_col]],
          na.rm = TRUE
        ),
        Brier = B_model,
        Brier_null = B_null,
        Brier_skill =
          1 - B_model / B_null,
        AUC = auc_roc(
          pred_df$y,
          pred_df$p
        )
      )
      
      list(
        metrics = metrics,
        preds = pred_df
      )
    }
  )
  
  metrics <- bind_rows(
    map(
      fold_results,
      "metrics"
    )
  )
  
  preds <- bind_rows(
    map(
      fold_results,
      "preds"
    )
  )
  
  B_model <- brier(
    preds$y,
    preds$p
  )
  
  B_null <- brier(
    preds$y,
    preds$p_null
  )
  
  pooled <- tibble(
    Brier = B_model,
    Brier_null = B_null,
    Brier_skill =
      1 - B_model / B_null,
    AUC = auc_roc(
      preds$y,
      preds$p
    )
  )
  
  list(
    metrics = metrics,
    preds = preds,
    pooled = pooled
  )
}

# ============================================================
# 5. RUN BOTH VALIDATION SCHEMES
# ============================================================

# Known sites
known_ug <- cv_known_site(
  UG_data,
  UG_f,
  k = 5,
  seed = 1
)

known_drc <- cv_known_site(
  DRC_data,
  DRC_f,
  k = 5,
  seed = 1
)

known_chad <- cv_known_site(
  CHAD_data,
  CHAD_f,
  k = 5,
  seed = 1
)

# Unseen sites
unseen_ug <- cv_unseen_site(
  UG_data,
  UG_f,
  k = 5,
  seed = 1
)

unseen_drc <- cv_unseen_site(
  DRC_data,
  DRC_f,
  k = 5,
  seed = 1
)

unseen_chad <- cv_unseen_site(
  CHAD_data,
  CHAD_f,
  k = 5,
  seed = 1
)

# ============================================================
# 6. POOLED VALIDATION SUMMARY
# ============================================================

validation_summary <- bind_rows(
  
  known_ug$pooled %>%
    mutate(
      Country = "Uganda",
      Validation = "Known sites"
    ),
  
  unseen_ug$pooled %>%
    mutate(
      Country = "Uganda",
      Validation = "Unseen sites"
    ),
  
  known_drc$pooled %>%
    mutate(
      Country = "DRC",
      Validation = "Known sites"
    ),
  
  unseen_drc$pooled %>%
    mutate(
      Country = "DRC",
      Validation = "Unseen sites"
    ),
  
  known_chad$pooled %>%
    mutate(
      Country = "Chad",
      Validation = "Known sites"
    ),
  
  unseen_chad$pooled %>%
    mutate(
      Country = "Chad",
      Validation = "Unseen sites"
    )
  
) %>%
  dplyr::select(
    Country,
    Validation,
    Brier,
    Brier_null,
    Brier_skill,
    AUC
  )

validation_summary

write.csv(
  validation_summary,
    "/Users/noni/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/07_tables/validation_summary.csv",
  row.names = FALSE
)

validation_summary %>%
  mutate(
    across(
      c(
        Brier,
        Brier_null,
        Brier_skill,
        AUC
      ),
      ~ round(.x, 3)
    )
  )

# ============================================================
# 7. FIGURE DATA
# ============================================================

fold_metrics <- bind_rows(
  
  known_ug$metrics %>%
    mutate(
      Country = "Uganda",
      Validation = "Known sites"
    ),
  
  unseen_ug$metrics %>%
    mutate(
      Country = "Uganda",
      Validation = "Unseen sites"
    ),
  
  known_drc$metrics %>%
    mutate(
      Country = "DRC",
      Validation = "Known sites"
    ),
  
  unseen_drc$metrics %>%
    mutate(
      Country = "DRC",
      Validation = "Unseen sites"
    ),
  
  known_chad$metrics %>%
    mutate(
      Country = "Chad",
      Validation = "Known sites"
    ),
  
  unseen_chad$metrics %>%
    mutate(
      Country = "Chad",
      Validation = "Unseen sites"
    )
)

unseen_preds <- bind_rows(
  
  unseen_ug$preds %>%
    mutate(
      Country = "Uganda"
    ),
  
  unseen_drc$preds %>%
    mutate(
      Country = "DRC"
    ),
  
  unseen_chad$preds %>%
    mutate(
      Country = "Chad"
    )
)

cal_unseen <- unseen_preds %>%
  group_by(Country) %>%
  mutate(bin = ntile(p, 5)) %>%
  group_by(Country, bin) %>%
  summarise(
    p_mean = mean(p, na.rm = TRUE),
    y_mean = mean(y, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    se = sqrt((y_mean * (1 - y_mean)) / n),
    y_low = pmax(0, y_mean - 1.96 * se),
    y_high = pmin(1, y_mean + 1.96 * se)
  )

cal_limits <- cal_unseen %>%
  group_by(Country) %>%
  summarise(
    lim = max(c(p_mean, y_high), na.rm = TRUE) * 1.08,
    .groups = "drop"
  )

cal_unseen <- cal_unseen %>%
  left_join(cal_limits, by = "Country")

country_levels <- c("Uganda", "DRC", "Chad")

fold_metrics <- fold_metrics %>%
  mutate(
    Country = factor(Country, levels = country_levels),
    Validation = factor(Validation, levels = c("Known sites", "Unseen sites"))
  )

cal_unseen <- cal_unseen %>%
  mutate(
    Country = factor(Country, levels = country_levels)
  )


pal <- c(
  "Uganda" = "#111A26",
  "DRC"    = "#3F618C",
  "Chad"   = "#6883A6"
)

# Colours for the two prediction targets in panel (a)
validation_pal <- c(
  "Known sites"  = "#5B5F97",
  "Unseen sites" = "#C97A40"
)

theme_val <- theme_classic(base_size = 9) +
  theme(
    text = element_text(colour = "black"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.tag = element_text(size = 10),
    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank()
  )


# ============================================================
# Prepare pooled Brier-skill values
# ============================================================

pooled_skill <- validation_summary %>%
  mutate(
    Country = factor(
      Country,
      levels = c("Uganda", "DRC", "Chad")
    ),
    Validation = factor(
      Validation,
      levels = c("Known sites", "Unseen sites")
    )
  )

fold_metrics <- fold_metrics %>%
  mutate(
    Country = factor(
      Country,
      levels = c("Uganda", "DRC", "Chad")
    ),
    Validation = factor(
      Validation,
      levels = c("Known sites", "Unseen sites")
    )
  )

pos_val <- position_dodge(width = 0.40)


# ============================================================
# (a) Brier skill
# small = individual folds
# large = pooled CV result
# ============================================================

p_skill <- ggplot() +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.45
  ) +
  
  # Individual folds
  geom_point(
    data = fold_metrics,
    aes(
      x = Country,
      y = Brier_skill,
      colour = Validation,
    ),
    position = pos_val,
    size = 1.5,
    alpha = 0.30
  ) +
  
  # Pooled result
  geom_point(
    data = pooled_skill,
    aes(
      x = Country,
      y = Brier_skill,
      colour = Validation,
    ),
    position = pos_val,
    size = 3.4,
    stroke = 0.7
  ) +
  
  scale_colour_manual(
    values = validation_pal
  ) +
  
  scale_shape_manual(
    values = c(
      "Known sites" = 16,
      "Unseen sites" = 17
    )
  ) +
  
  labs(
    x = NULL,
    y = "Brier skill score",
    colour = "Prediction target",
  ) +
  
  theme_val +
  
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )


# ============================================================
# (b) Calibration — unseen sites only
# ============================================================
p_cal <- ggplot(
  cal_unseen,
  aes(
    x = p_mean,
    y = y_mean
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = y_low,
      ymax = y_high,
      colour = Country
    ),
    width = 0,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  geom_point(
    aes(
      size = n,
      colour = Country
    ),
    alpha = 0.95,
    show.legend = c(
      colour = FALSE,
      size = TRUE
    )
  ) +
  facet_wrap(
    ~ Country,
    nrow = 1,
    scales = "free"
  ) +
  scale_colour_manual(
    values = pal
  ) +
  scale_size_continuous(
    name = "Observations per bin",
    range = c(2.5, 6),
    breaks = c(250, 500, 750)
  ) +
  guides(
    size = guide_legend(
      override.aes = list(
        colour = "grey45",
        alpha = 1
      )
    )
  ) +
  labs(
    x = "Mean predicted probability",
    y = "Observed frequency"
  ) +
  theme_val +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "plain",
      size = 9
    )
  )

# ============================================================
# Combined figure
# ============================================================


fig_val <- p_skill / p_cal +
  plot_layout(
    heights = c(0.8, 1.05)
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag.position = c(0.01, 0.99)
  )

fig_val

ggsave(
  filename = "Figure6_validation.svg",
  plot = fig_val,
  width = 18,
  height = 15,
  units = "cm",
  device = svglite::svglite,
  bg = "transparent"
)


## Figure 3 , presence / predictions

## proportion pres / abs 

library(dplyr)
library(ggplot2)
library(patchwork)

# lme4::predict() includes random effects by default.
# For fixed-effects only, use: predict(mod, type = "response", re.form = NA)
add_phat <- function(dat, mod, country_label) {
  dat %>%
    mutate(
      Country = country_label,
      y = bio_pres,
      p_hat = predict(mod, type = "response")
    )
}

UG_best_model   <- fit_from_string(UG_best_row$Formula[1],   UG_data)
DRC_best_model  <- fit_from_string(DRC_best_row$Formula[1],  DRC_data)
CHAD_best_model <- fit_from_string(CHAD_best_row$Formula[1], CHAD_data)

pred_all <- bind_rows(
  add_phat(UG_data,   UG_best_model,   "Uganda"),
  add_phat(DRC_data,  DRC_best_model,  "DRC"),
  add_phat(CHAD_data, CHAD_best_model, "Chad")
)

site_stats <- pred_all %>%
  filter(!is.na(site_uid), !is.na(y), !is.na(p_hat)) %>%
  group_by(Country, site_uid) %>%
  summarise(
    n_obs = n(),
    P     = sum(y == 1, na.rm = TRUE),
    A     = sum(y == 0, na.rm = TRUE),
    bal   = (P - A) / (P + A),
    p_med = median(p_hat, na.rm = TRUE),
    p_q05 = quantile(p_hat, 0.05, na.rm = TRUE),
    p_q25 = quantile(p_hat, 0.25, na.rm = TRUE),
    p_q75 = quantile(p_hat, 0.75, na.rm = TRUE),
    p_q95 = quantile(p_hat, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_obs >= 10) %>%
  mutate(
    visits_bin = cut(
      n_obs,
      breaks = c(10, 20, 50, 100, Inf),
      labels = c("10–19", "20–49", "50–99", "≥100"),
      right = FALSE
    )
  ) %>%
  group_by(Country) %>%
  arrange(bal, .by_group = TRUE) %>%
  mutate(site_key = factor(site_uid, levels = site_uid)) %>%
  ungroup()

alpha_vals  <- c("10–19" = 0.30, "20–49" = 0.55, "50–99" = 0.80, "≥100" = 1.00)
bal_breaks  <- c(-1, 0, 1)
prob_breaks <- c(0, 0.5, 1)
prob_labels <- c("0", "0.5", "1.0")

make_pair <- function(cty, dat = site_stats) {
  d <- filter(dat, Country == cty)
  
  p1 <- ggplot(d, aes(x = bal, y = site_key)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_col(fill = "azure3", width = 0.65) +
    scale_x_continuous(limits = c(-1, 1), breaks = bal_breaks) +
    labs(title = cty, x = "(P − A)/(P + A)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  p2 <- ggplot(d, aes(y = site_key)) +
    geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.3) +
    geom_linerange(aes(xmin = p_q05, xmax = p_q95, alpha = visits_bin), linewidth = 0.35) +
    geom_linerange(aes(xmin = p_q25, xmax = p_q75, alpha = visits_bin), linewidth = 1.05) +
    geom_point(aes(x = p_med, alpha = visits_bin), size = 1.2) +
    scale_x_continuous(limits = c(0, 1), breaks = prob_breaks, labels = prob_labels) +
    scale_alpha_manual(values = alpha_vals, name = "Visits per site", drop = FALSE) +
    labs(x = "Pred. pres.", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  p1 + p2 + plot_layout(widths = c(1, 1))
}

p_row <- wrap_plots(
  lapply(c("Chad", "DRC", "Uganda"), make_pair),
  ncol = 3
) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.location = "plot"
  )

p_row

ggsave(
  filename = "Figure_3_complete.svg",
  plot = p_row,
  width = 18,
  height = 11,
  units = "cm",
  device = svglite::svglite,
  bg = "transparent"
)

# ============================================================
# Models with recordedBy
# ============================================================

UG_models_ID <- read.csv(
  "~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/04_model_outputs/Models_output_UG_upd_Oct26_ID.csv"
)

DRC_models_ID <- read.csv(
  "~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/04_model_outputs/Models_output_DRC_upd_Oct26_ID.csv"
)

CHAD_models_ID <- read.csv(
  "~/OneDrive - KU Leuven/UG_DRC_CHAD/PaperRevision_Oct26/04_model_outputs/Models_output_CHAD_upd_Oct26_ID.csv"
)


# ============================================================
# Helper: extract fixed-effect structure only
# ============================================================

get_fixed_structure <- function(formula_string) {
  
  f <- as.formula(formula_string)
  
  # Remove random-effect terms from formula
  fixed_f <- lme4::nobars(f)
  
  # Extract fixed-effect term labels
  terms_fixed <- attr(
    terms(fixed_f),
    "term.labels"
  )
  
  # Sort so term order in the formula does not matter
  paste(
    sort(terms_fixed),
    collapse = " + "
  )
}

UG_target_structure <- get_fixed_structure(
  UG_best_row$Formula[1]
)

DRC_target_structure <- get_fixed_structure(
  DRC_best_row$Formula[1]
)

CHAD_target_structure <- get_fixed_structure(
  CHAD_best_row$Formula[1]
)

UG_target_structure
DRC_target_structure
CHAD_target_structure

UG_same_ID <- UG_models_ID %>%
  mutate(
    fixed_structure = sapply(
      Formula,
      get_fixed_structure
    )
  ) %>%
  filter(
    fixed_structure == UG_target_structure
  )


DRC_same_ID <- DRC_models_ID %>%
  mutate(
    fixed_structure = sapply(
      Formula,
      get_fixed_structure
    )
  ) %>%
  filter(
    fixed_structure == DRC_target_structure
  )


CHAD_same_ID <- CHAD_models_ID %>%
  mutate(
    fixed_structure = sapply(
      Formula,
      get_fixed_structure
    )
  ) %>%
  filter(
    fixed_structure == CHAD_target_structure
  )

UG_same_ID %>%
  dplyr::select(Formula, AIC, R2_m, R2_c)

DRC_same_ID %>%
  dplyr::select(Formula, AIC, R2_m, R2_c)

CHAD_same_ID %>%
  dplyr::select(Formula, AIC, R2_m, R2_c)


## ## Chad subset - to see experience 

CHAD_data <- df_all %>% filter(country == "CHAD")

CHAD_reduced = CHAD_data %>% filter(eventDate < '2024-07-01')

CHAD_reduced_model   <- fit_from_string(CHAD_best_row$Formula[1], CHAD_reduced)

pvals <- summary(CHAD_reduced_model)$coefficients[, "Pr(>|z|)"]

or_CHAD <- tidy_or(CHAD_best_model, "Chad")

or_CHAD_1 <- tidy_or(CHAD_reduced_model, "Chad")


## With reduced dataset 

CHAD_reduced <- CHAD_data %>%
  filter(eventDate < '2024-07-01') %>%
  mutate(     jul_dat_z = as.numeric(scale(jul_dat))
    )

# residualise WITHIN the reduced dataset
exp_mod_reduced <- lm(
  sampling_experience ~ jul_dat_z,
  data = CHAD_reduced
)

CHAD_reduced$exp_resid <- residuals(exp_mod_reduced)

CHAD_reduced_time_model <- glmer(
  bio_pres ~ pp_oneweek +
    NDVI_q90_t +
    prop_trees_med +
    exp_resid +
    jul_dat +
    NDVI_modis +
    (1 | site_uid),
  data = CHAD_reduced,
  family = binomial(),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

summary(CHAD_reduced_time_model)

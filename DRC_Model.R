library(dplyr)
library(lme4)
library(MuMIn)
library(future)
library(future.apply)
library(parallelly)

# ============================================================
# DATA
# ============================================================

df_scaled <- read.csv(
  "~/OneDrive - KU Leuven/UG_DRC_CHAD/Data/df_final_scaled_Bio_upd 2.csv"
)

# This is DRC, so I would rename it to avoid confusion
drc_data_model <- df_scaled %>%
  filter(country == "DRC")


# ============================================================
# PREDICTORS
# ============================================================

fixed_effects <- c(
  "pp_oneweek",
  "temp_modis",
  "NDVI_q90_t",
  "elev",
  "flow_acc",
  "SS",
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


# ============================================================
# IMPORTANT: USE THE SAME OBSERVATIONS FOR EVERY MODEL
# ============================================================

all_model_vars <- c(
  "bio_pres",
  "site_uid",
  fixed_effects
)

n_before <- nrow(drc_data_model)

drc_data_model <- drc_data_model[
  complete.cases(drc_data_model[, all_model_vars]),
]

cat(
  "Rows before:", n_before,
  "\nRows used in all models:", nrow(drc_data_model),
  "\n"
)


# ============================================================
# MODEL CONTROL
# ============================================================

control_params <- glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 5e5)
)


# ============================================================
# GENERATE ALL COMBINATIONS
# ============================================================

generate_combinations <- function(effects) {
  
  unlist(
    lapply(
      seq_along(effects),
      function(k) combn(effects, k, simplify = FALSE)
    ),
    recursive = FALSE
  )
}

combinations <- generate_combinations(fixed_effects)

length(combinations)
# 32767


# ============================================================
# NULL MODEL
# Fit this ONCE for R2 calculation
# ============================================================

null_model <- glmer(
  bio_pres ~ 1 + (1 | site_uid),
  data = drc_data_model,
  family = binomial,
  control = control_params
)


# ============================================================
# PARALLEL SETUP
# ============================================================

available <- parallelly::availableCores()

workers <- max(
  1,
  min(6, available - 1)
)

cat(
  "Available cores:", available,
  "\nUsing workers:", workers,
  "\n"
)

plan(
  multisession,
  workers = workers
)


# ============================================================
# FUNCTION FOR ONE MODEL
# ============================================================

fit_one_model <- function(vars) {
  
  fixed_formula <- paste(
    vars,
    collapse = " + "
  )
  
  formula <- as.formula(
    paste(
      "bio_pres ~",
      fixed_formula,
      "+ (1 | site_uid)"
    )
  )
  
  # --------------------------------
  # Fit model
  # --------------------------------
  
  model <- tryCatch(
    glmer(
      formula,
      data = drc_data_model,
      family = binomial,
      control = control_params
    ),
    error = function(e) NULL
  )
  
  
  # --------------------------------
  # If model fails
  # --------------------------------
  
  if (is.null(model)) {
    
    return(
      data.frame(
        Formula = paste(deparse(formula), collapse = " "),
        AIC = NA_real_,
        R2_m = NA_real_,
        R2_c = NA_real_,
        Significant_Vars = NA_character_,
        Converged = FALSE,
        Singular = NA
      )
    )
  }
  
  
  # --------------------------------
  # AIC
  # --------------------------------
  
  model_AIC <- AIC(model)
  
  
  # --------------------------------
  # R2
  #
  # IMPORTANT:
  # reuse the already fitted null model
  # --------------------------------
  
  r2_values <- suppressWarnings(
    MuMIn::r.squaredGLMM(
      model,
      null = null_model
    )
  )
  
  R2_m <- r2_values[1, "R2m"]
  R2_c <- r2_values[1, "R2c"]
  
  
  # --------------------------------
  # P-values
  # --------------------------------
  
  coef_table <- coef(summary(model))
  
  pvals <- coef_table[, "Pr(>|z|)"]
  
  # Remove intercept because normally we don't
  # want "(Intercept)" reported as an explanatory variable
  pvals <- pvals[
    names(pvals) != "(Intercept)"
  ]
  
  significant_vars <- names(
    pvals[pvals < 0.05]
  )
  
  
  # --------------------------------
  # Convergence information
  # --------------------------------
  
  convergence_message <-
    model@optinfo$conv$lme4$messages
  
  converged <- is.null(
    convergence_message
  )
  
  singular <- lme4::isSingular(
    model,
    tol = 1e-4
  )
  
  
  # --------------------------------
  # Return one row
  # --------------------------------
  
  data.frame(
    Formula =
      paste(deparse(formula), collapse = " "),
    
    AIC =
      model_AIC,
    
    R2_m =
      R2_m,
    
    R2_c =
      R2_c,
    
    Significant_Vars =
      if (
        length(significant_vars) > 0
      ) {
        paste(
          significant_vars,
          collapse = ", "
        )
      } else {
        "None"
      },
    
    Converged =
      converged,
    
    Singular =
      singular,
    
    stringsAsFactors = FALSE
  )
}


# ============================================================
# RUN ALL 32,767 MODELS IN PARALLEL
# ============================================================

library(progressr)

handlers(
  handler_progress(type = "time")
)


with_progress({
  
  p <- progressor(
    steps = length(combinations)
  )
  
  results_list <- future_lapply(
    seq_along(combinations),
    
    function(i) {
      
      result <- fit_one_model(
        combinations[[i]]
      )
      
      p(
        message = paste0(
          "Model ",
          i,
          " / ",
          length(combinations)
        )
      )
      
      result
    },
    
    future.packages = c(
      "lme4",
      "MuMIn",
      "progressr"
    ),
    
    future.scheduling = 2
  )
})

# ============================================================
# COMBINE RESULTS
# ============================================================

results_DRC <- bind_rows(
  results_list
)


# ============================================================
# ADD DELTA AIC + SORT
# ============================================================

results_DRC <- results_DRC %>%
  mutate(
    delta_AIC = AIC - min(
      AIC,
      na.rm = TRUE
    )
  ) %>%
  arrange(AIC)


# ============================================================
# RETURN R TO NORMAL
# ============================================================

plan(sequential)


# ============================================================
# SAVE
# ============================================================

write.csv(
  results_DRC,
  "~/Library/CloudStorage/OneDrive-KULeuven/UG_DRC_CHAD/Models_output_DRC_upd_Oct26.csv",
  row.names = FALSE
)

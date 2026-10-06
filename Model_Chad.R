df_scaled = read.csv('df_final_scaled_Bio_upd.csv')


chad_data_model = df_scaled %>% filter(country=='CHAD')

control_params <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5e5))

# Full formula 

full_formula <- bio_pres ~ pp_oneweek + temp_modis + NDVI_q90_t + 
  elev + flow_acc + SS + pop_rate_change + prop_water_med +
  prop_trees_med + prop_crops_med + prop_built_med + sampling_experience + 
  time_sin + time_cos + NDVI_modis

# Extract all fixed effect terms from the full model (without the interaction term)
fixed_effects <- c("pp_oneweek", 'temp_modis','NDVI_q90_t',
                   'elev','flow_acc', 'SS','pop_rate_change', 'prop_water_med',
                   'prop_trees_med','prop_crops_med', 'prop_built_med','sampling_experience', 
                   'time_sin','time_cos','NDVI_modis')

# Function to generate all subsets of a set (without the empty set)
generate_combinations <- function(effects) {
  unlist(lapply(1:length(effects), function(k) combn(effects, k, simplify = FALSE)), recursive = FALSE)
}


# Generate all combinations of the fixed effects
combinations <- generate_combinations(fixed_effects)

# Initialize storage for AIC, R2 values, and formulas
results <- data.frame(Formula = character(length(combinations)),
                      AIC = numeric(length(combinations)),
                      R2_m = numeric(length(combinations)),
                      R2_c = numeric(length(combinations)),
                      Significant_Vars = character(length(combinations)),  # Store significant variables
                      stringsAsFactors = FALSE)

library(lme4)
library(MuMIn)
library(future.apply)
library(progressr)
library(cli)

# 1. SET THREAD LIMITS BEFORE SPAWNING WORKERS
# This ensures all 16 background sessions inherit 1-thread limit for linear algebra
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# 2. Use 'multisession' across 16 background R processes
plan(multisession, workers = 16)

# 3. Configure progressr to display live models/second ({cli::pb_rate})
handlers(
  handler_cli(
    format = "{cli::pb_spin} [{cli::pb_bar}] {cli::pb_percent} | {cli::pb_current}/{cli::pb_total} | {cli::pb_rate} | ETA: {cli::pb_eta}",
    clear = FALSE
  )
)

# 4. Worker function
fit_single_model <- function(comb, data_obj, ctrl) {
  fixed_formula <- paste(comb, collapse = " + ")
  form <- as.formula(paste("bio_pres ~", fixed_formula, "+ (1 | recordedBy)"))
  
  res <- tryCatch({
    model <- glmer(form, data = data_obj, family = binomial, control = ctrl)
    
    r2_values <- r.squaredGLMM(model)
    pvals <- summary(model)$coefficients[, "Pr(>|z|)"]
    significant_vars <- names(pvals)[pvals < 0.05]
    sig_str <- if (length(significant_vars) > 0) paste(significant_vars, collapse = ", ") else "None"
    
    data.frame(
      Formula = paste(deparse(form), collapse = " "),
      AIC = AIC(model),
      R2_m = r2_values[1],
      R2_c = r2_values[3],
      Significant_Vars = sig_str,
      Status = "Success",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      Formula = paste(deparse(form), collapse = " "),
      AIC = NA,
      R2_m = NA,
      R2_c = NA,
      Significant_Vars = NA,
      Status = paste("Error:", e$message),
      stringsAsFactors = FALSE
    )
  })
  
  return(res)
}

# 5. Run parallel loop with live speed tracking
with_progress({
  p <- progressor(along = combinations)
  
  results_list <- future_lapply(
    combinations, 
    function(comb) {
      res <- fit_single_model(comb, data_obj = chad_data_model, ctrl = control_params)
      p() # Updates count & speed calculation
      return(res)
    }, 
    future.seed = TRUE,
    future.packages = c("lme4", "MuMIn"),
    future.chunk.size = 200 # Minimizes communication lag between cores
  )
})

# Combine results and clean up workers
results <- do.call(rbind, results_list)
plan(sequential)

write.csv(results,'Models_output_CHAD_upd_Oct26_ID.csv')

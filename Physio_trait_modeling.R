library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(forcats)
library(lme4)
library(lmerTest)
library(arm)  # for rescale

all <- read.csv("acclim_data_all.csv") #dataset contains all individual names and experimental information

mod <- nls(Final_V02max~a*(Mass^b), data=all, start=list(a=1, b=1)) #correct for non-linear mass vs Vo2max relationship
summary(mod)
resid <- resid(mod)
all$FinalVO2max_residuals <- resid

####function for creating data summary df and plotting, rescaling raw data for comparartive purposes, and fitting GLS/PGLS models####
analyze_variable_test <- function(data, varname) {
  
  # Calculate summary statistics for the given variable
  var_sum <- data_summary(data, varname, groupnames = c("Treatment", "Strain"))
  
  # Create a plot with summary statistics
  var_plot <- var_sum %>%
    dplyr::mutate(Strain = forcats::fct_relevel(Strain, "PO", "LS", "LN", "ML", "AZ", "ME")) %>%
    ggplot(aes(x = as.factor(Treatment), y = !!sym(varname), group = Strain, color = Strain)) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = !!sym(varname) - std, ymax = !!sym(varname) + std),
                  width = .3, position = position_dodge(0.05), linewidth = .5) +
    scale_color_manual(values = color_frame$color) +
    scale_x_discrete(labels = c("1000" = "Low", "3000" = "Med", "4500" = "High")) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab(varname) +
    xlab("Hypoxia Treatment by Taxa") +
    facet_grid(. ~ Strain, scales = "free")
  # Fit GLS models for each strain and extract coefficient estimates
  gls_data <- data %>%
    dplyr::mutate(
      !!varname := arm::rescale(.data[[varname]]),# rescale response (mean=0, SD=0.5)
      TreatVal  = as.numeric(as.character(Treatment))# numeric predictor for slope
    )
  
  gls_model_df <- gls_data %>%
    dplyr::group_by(Strain) %>%
    dplyr::summarise(
      gls_model = list(
        nlme::gls(stats::reformulate("TreatVal", response = varname),  # slope on numeric TreatVal
                  weights     = nlme::varIdent(form = ~ 1 | Treatment),              # allow var by Treatment levels
                  na.action   = na.omit,
                  correlation = nlme::corCompSymm(form = ~ 1 | Family),
                  data        = dplyr::pick(dplyr::everything())
        )
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      estimate  = purrr::map_dbl(gls_model, ~ unname(coef(.x)["TreatVal"])),
      std_error = purrr::map_dbl(gls_model, ~ summary(.x)$tTable["TreatVal","Std.Error"]),
      t_value   = purrr::map_dbl(gls_model, ~ summary(.x)$tTable["TreatVal","t-value"]),
      p_value   = purrr::map_dbl(gls_model, ~ summary(.x)$tTable["TreatVal","p-value"])
    ) %>%
    dplyr::mutate(
      elev   = dplyr::case_when(
        Strain == "AZ" ~ 2700,
        Strain == "ML" ~ 1500,
        Strain == "ME" ~ 4350,
        Strain %in% c("LS", "LN") ~ 430,
        Strain == "PO" ~ 0
      ),
      Strain = forcats::fct_relevel(Strain, "PO", "LS", "LN", "ML", "AZ", "ME")
    )
  
  #Slope vs elevation (PGLS)
  slope_model   <- nlme::gls(estimate ~ elev, data = gls_model_df,
                             correlation = ape::corPagel(1, CAE_cytb_tree, form = ~ Strain))
  slope_summary <- summary(slope_model)
  slope_estimate <- stats::coef(slope_model)[2]
  p_val <- slope_summary$tTable[2, 4]
  
  fit_data <- data.frame(
    elev = gls_model_df$elev,
    fitted_values = stats::predict(slope_model)
  )
  
  var_slope_plot <- gls_model_df %>%
    ggplot(aes(x = elev, y = estimate, color = Strain)) +
    geom_point(size = 5) +
    geom_line(data = fit_data, aes(x = elev, y = fitted_values),
              color = "black", linewidth = 0.5, linetype = "dashed") +
    scale_color_manual(values = color_frame$color) +
    theme_bw() +
    theme(legend.position = "none") +
    ylim(-0.00016, 0.0003) +
    labs(y = "Flexibility", x = "Native Elevation (m a.s.l.)") +
    annotate("text", x = 2500, y = max(gls_model_df$estimate, na.rm = TRUE) - 0.1*max(gls_model_df$estimate),
             label = paste0("Slope = ", round(slope_estimate, 10), "\np-val = ", round(p_val, 10)),
             hjust = 0, size = 3, color = "black")
  
  list(
    var_plot      = var_plot,
    var_slope_plot = var_slope_plot,
    gls_results   = as.data.frame(gls_model_df) %>% dplyr::select(-gls_model),
    slope_results = slope_summary
  )
}

#example useage on hemoglobin concentration as a trait of interest
hb_result_rescale <- analyze_variable_test(data = all, varname = "hb")
hb_result_rescale$var_plot
hb_result_rescale$var_slope_plot
hb_result_rescale$gls_results
hb_result_rescale$slope_results
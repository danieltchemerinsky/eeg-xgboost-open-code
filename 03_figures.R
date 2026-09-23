# ==============================================================================
# PUBLICATION FIGURES
# ==============================================================================
# 
# Figure functions. Sourcing this file defines them; it draws nothing.
# Run 01_loocv_analysis.R first to create the objects they need.
#
# REQUIRED OBJECTS (from main analysis):
# - loo_predictions: matrix of class probabilities
# - actual_outcomes: vector of true class labels (0,1,2,3 for CPC 1,2,3,5)
# - train_data: predictor data frame
# - train_label: outcome labels
# - params: XGBoost parameters
# - nrounds: number of boosting rounds
# - full: full dataset with Patient and Outcome columns
# - table1, table2, auc_results: performance metrics
#
# ==============================================================================

# Load required packages
library(ggplot2)
library(pROC)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(xgboost)

# Check if patchwork is installed, if not provide instructions
if(!requireNamespace("patchwork", quietly = TRUE)) {
  message("Please install patchwork: install.packages('patchwork')")
}

# ==============================================================================
# NATURE-STYLE THEME
# ==============================================================================

theme_nature <- function(base_size = 8) {
  # Use sans-serif font (works on all systems)
  # On Windows: uses Arial, on Mac: Helvetica, on Linux: DejaVu Sans
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      # Panel
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Axes
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.1, "cm"),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size + 1, face = "bold"),
      
      # Legend
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.position = "bottom",
      
      # Plot
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0, color = "grey40"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Colorblind-friendly palette
nature_colors <- c(
  "blue" = "#2166AC",
  "red" = "#B2182B", 
  "green" = "#1B7837",
  "orange" = "#E66101",
  "purple" = "#762A83",
  "grey" = "#525252",
  "black" = "#000000"
)

# ==============================================================================
# FIGURE 2: FEATURE IMPORTANCE (FOREST PLOT WITH 95% CI)
# ==============================================================================

generate_figure2 <- function(all_importance = NULL, train_data = NULL, 
                              train_label = NULL, params = NULL, nrounds = NULL) {
  
  cat("Generating Figure 2: Feature Importance with 95% CI...\n")
  
  # If all_importance from LOO-CV is provided, use it
  if(!is.null(all_importance)) {
    
    # Calculate summary statistics from LOO folds
    summary_stats <- all_importance %>%
      dplyr::select(Feature, Fold, Gain) %>%
      tidyr::complete(Feature, Fold = 1:max(Fold), fill = list(Gain = 0)) %>%
      dplyr::group_by(Feature) %>%
      dplyr::summarise(
        Median = median(Gain),
        Lower_95 = quantile(Gain, 0.025),
        Upper_95 = quantile(Gain, 0.975),
        Mean = mean(Gain),
        .groups = "drop"
      ) %>%
      dplyr::arrange(desc(Mean)) %>%
      head(12)
    
  } else {
    # The bars in Figure 2 are the spread of each feature's Gain across the
    # leave-one-out folds, so they cannot be drawn without `all_importance`.
    # A single model fitted to all patients gives one Gain per feature and no
    # spread at all; drawing bars around it would imply information the figure
    # does not have.
    stop("`all_importance` not found. Run 01_loocv_analysis.R first: Figure 2 ",
         "needs the per-fold feature importance it produces.")
  }
  
  # Clean feature names
  summary_stats <- summary_stats %>%
    dplyr::mutate(
      Feature_Clean = Feature,
      Feature_Clean = gsub("\\.Post$", "", Feature_Clean),
      Feature_Clean = gsub("\\.Post\\.", ".", Feature_Clean),
      Feature_Clean = gsub("Sound", "Ear", Feature_Clean),
      Feature_Clean = gsub("Stern", "Sternum", Feature_Clean),
      Feature_Clean = gsub("LeftArm", "Left Arm", Feature_Clean),
      Feature_Clean = gsub("RightArm", "Right Arm", Feature_Clean),
      Feature_Clean = gsub("LeftLeg", "Left Leg", Feature_Clean),
      Feature_Clean = gsub("\\.", " · ", Feature_Clean)
    )
  
  # Reorder for plotting
  summary_stats$Feature_Clean <- factor(
    summary_stats$Feature_Clean, 
    levels = rev(summary_stats$Feature_Clean)
  )
  
  # Create forest plot
  fig2 <- ggplot(summary_stats, aes(x = Median, y = Feature_Clean)) +
    # Spread across folds. geom_errorbarh() was deprecated in ggplot2 4.0.0;
    # the horizontal form is now geom_errorbar() with orientation = "y".
    geom_errorbar(
      aes(xmin = Lower_95, xmax = Upper_95),
      orientation = "y",
      width = 0.25,
      color = "#525252",  # grey
      linewidth = 0.6
    ) +
    # Median points
    geom_point(
      size = 2.5,
      color = "#B2182B",  # red
      shape = 16
    ) +
    # Scales
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.1)),
      breaks = seq(0, 0.15, 0.03)
    ) +
    labs(
      x = "Gain (Relative Contribution)",
      y = NULL
    ) +
    theme_nature(base_size = 9) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 8)
    )
  
  return(fig2)
}


# ==============================================================================
# FIGURE 3: ONE-VS-ALL ROC CURVES
# ==============================================================================

generate_figure3 <- function(loo_predictions, actual_outcomes, n_patients = 114) {
  
  cat("Generating Figure 3: One-vs-All ROC Curves with 95% CI...\n")
  
  # Calculate ROC for each CPC class
  roc_cpc1 <- roc(response = (actual_outcomes == 0), 
                  predictor = loo_predictions[, 1], quiet = TRUE)
  ci_cpc1 <- ci.auc(roc_cpc1, method = "delong")
  
  roc_cpc2 <- roc(response = (actual_outcomes == 1), 
                  predictor = loo_predictions[, 2], quiet = TRUE)
  ci_cpc2 <- ci.auc(roc_cpc2, method = "delong")
  
  roc_cpc5 <- roc(response = (actual_outcomes == 3), 
                  predictor = loo_predictions[, 4], quiet = TRUE)
  ci_cpc5 <- ci.auc(roc_cpc5, method = "delong")
  
  # Print AUC values
  cat("\n  Updated AUC Values:\n")
  cat(sprintf("  CPC 1 vs All: %.3f (95%% CI: %.2f-%.2f)\n", 
              auc(roc_cpc1), ci_cpc1[1], ci_cpc1[3]))
  cat(sprintf("  CPC 2 vs All: %.3f (95%% CI: %.2f-%.2f)\n", 
              auc(roc_cpc2), ci_cpc2[1], ci_cpc2[3]))
  cat(sprintf("  CPC 5 vs All: %.3f (95%% CI: %.2f-%.2f)\n", 
              auc(roc_cpc5), ci_cpc5[1], ci_cpc5[3]))
  
  # Create legend labels with CIs
  label_cpc1 <- sprintf("CPC 1 vs All, AUC: %.2f (%.2f-%.2f)",
                        auc(roc_cpc1), ci_cpc1[1], ci_cpc1[3])
  label_cpc2 <- sprintf("CPC 2 vs All, AUC: %.2f (%.2f-%.2f)",
                        auc(roc_cpc2), ci_cpc2[1], ci_cpc2[3])
  label_cpc5 <- sprintf("CPC 5 vs All, AUC: %.2f (%.2f-%.2f)",
                        auc(roc_cpc5), ci_cpc5[1], ci_cpc5[3])
  
  # Combine ROC data - sort by specificity descending for proper line connection
  roc_data <- bind_rows(
    data.frame(
      Specificity = roc_cpc1$specificities,
      Sensitivity = roc_cpc1$sensitivities,
      Class = label_cpc1
    ) %>% arrange(desc(Specificity), Sensitivity),
    data.frame(
      Specificity = roc_cpc2$specificities,
      Sensitivity = roc_cpc2$sensitivities,
      Class = label_cpc2
    ) %>% arrange(desc(Specificity), Sensitivity),
    data.frame(
      Specificity = roc_cpc5$specificities,
      Sensitivity = roc_cpc5$sensitivities,
      Class = label_cpc5
    ) %>% arrange(desc(Specificity), Sensitivity)
  )
  
  # Set factor order for legend
  roc_data$Class <- factor(roc_data$Class, 
                           levels = c(label_cpc1, label_cpc2, label_cpc5))
  
  # Create plot using geom_path (connects points in data order)
  fig3 <- ggplot(roc_data, aes(x = Specificity, y = Sensitivity, 
                                color = Class, linetype = Class)) +
    # Reference diagonal
    geom_abline(intercept = 1, slope = 1, 
                linetype = "dotted", color = "grey60", linewidth = 0.4) +
    # ROC curves - geom_path connects in data order (which we sorted)
    geom_step(linewidth = 0.8, direction = "vh") +
    # Scales
    scale_x_reverse(
      name = paste0("Specificity (N = ", n_patients, ")"),
      limits = c(1, 0),
      breaks = seq(1, 0, -0.2),
      expand = c(0.02, 0.02)
    ) +
    scale_y_continuous(
      name = "Sensitivity",
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = c(0.02, 0.02)
    ) +
    scale_color_manual(
      values = c("#000000", "#B2182B", "#2166AC")  # black, red, blue
    ) +
    scale_linetype_manual(values = c("dashed", "dotdash", "solid")) +
    coord_equal() +
    theme_nature(base_size = 9) +
    theme(
      # Anchor the legend's bottom-right corner to the bottom-right of the
      # panel. Positioning by centre (the previous c(0.65, 0.22)) pushed the
      # long CI labels off the right edge, clipping the intervals mid-number.
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.key.width = unit(0.75, "cm"),
      # Was 13pt, which overrode the theme and left the legend larger than the
      # axis labels. Keep it just below the 9pt base size.
      legend.text = element_text(size = 7),
      legend.spacing.y = unit(0.02, "cm"),
      legend.margin = margin(1, 1, 1, 1)
    ) +
    guides(color = guide_legend(ncol = 1), 
           linetype = guide_legend(ncol = 1))
  
  return(fig3)
}


# ==============================================================================
# FIGURE 4: Clean Version - Good vs Poor Outcome ROC (Two Panels)
# ==============================================================================

generate_figure4_clean <- function(roc_114_react, roc_114_all,
                                   roc_103_react, roc_103_all,
                                   sep_sens, sep_spec) {

  library(ggplot2)
  library(pROC)
  library(patchwork)

  # Print SSEP values for debugging
  cat("SSEP values:\n")
  cat("  Sensitivity:", sep_sens, "\n")
  cat("  Specificity:", sep_spec, "\n")

  # Helper to extract and prepare ROC data with proper step function
  get_roc_df <- function(roc_obj) {
    # Get unique thresholds and corresponding sens/spec
    df <- data.frame(
      spec = roc_obj$specificities,
      sens = roc_obj$sensitivities
    )
    # Remove duplicates and sort by specificity (descending) then sensitivity
    df <- df[!duplicated(df), ]
    df <- df[order(-df$spec, df$sens), ]
    return(df)
  }

  # Extract ROC data
  df_114_r <- get_roc_df(roc_114_react)
  df_114_a <- get_roc_df(roc_114_all)
  df_103_r <- get_roc_df(roc_103_react)
  df_103_a <- get_roc_df(roc_103_all)

  # Get AUC values with CI
  auc_114_r <- sprintf("%.2f", as.numeric(auc(roc_114_react)))
  auc_114_a <- sprintf("%.2f", as.numeric(auc(roc_114_all)))
  auc_103_r <- sprintf("%.2f", as.numeric(auc(roc_103_react)))
  auc_103_a <- sprintf("%.2f", as.numeric(auc(roc_103_all)))

  # --- Panel A: N = 114 ---
  panel_a <- ggplot() +
    # Reference diagonal
    geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),
                 linetype = "dotted", color = "gray50", linewidth = 0.5) +
    # Reactivity curve (black dashed) - use geom_step for proper ROC appearance
    geom_step(data = df_114_r, aes(x = spec, y = sens),
              color = "black", linetype = "dashed", linewidth = 0.8, direction = "vh") +
    # All predictors curve (blue solid)
    geom_step(data = df_114_a, aes(x = spec, y = sens),
              color = "#2166AC", linetype = "solid", linewidth = 0.8, direction = "vh") +
    scale_x_reverse(limits = c(1, 0)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Specificity (N = 114)", y = "Sensitivity", tag = "A") +
    # Manual legend using annotate - positioned to avoid cropping
    annotate("segment", x = 0.55, xend = 0.40, y = 0.20, yend = 0.20,
             color = "black", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = 0.38, y = 0.20,
             label = paste0("Reactivity (", auc_114_r, ")"),
             hjust = 0, size = 2.8) +
    annotate("segment", x = 0.55, xend = 0.40, y = 0.12, yend = 0.12,
             color = "#2166AC", linetype = "solid", linewidth = 0.8) +
    annotate("text", x = 0.38, y = 0.12,
             label = paste0("All Predictors (", auc_114_a, ")"),
             hjust = 0, size = 2.8) +
    theme_minimal(base_size = 11, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 9, color = "black"),
      plot.tag = element_text(size = 14, face = "bold"),
      aspect.ratio = 1
    )

  # --- Panel B: N = 103 with SSEP ---
  panel_b <- ggplot() +
    # Reference diagonal
    geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),
                 linetype = "dotted", color = "gray50", linewidth = 0.5) +
    # Reactivity curve (black dashed)
    geom_step(data = df_103_r, aes(x = spec, y = sens),
              color = "black", linetype = "dashed", linewidth = 0.8, direction = "vh") +
    # All predictors curve (blue solid)
    geom_step(data = df_103_a, aes(x = spec, y = sens),
              color = "#2166AC", linetype = "solid", linewidth = 0.8, direction = "vh") +
    # SSEP point (green) - placed at correct coordinates
    geom_point(aes(x = sep_spec, y = sep_sens),
               color = "#228B22", size = 2) +
    scale_x_reverse(limits = c(1, 0)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Specificity (N = 103)", y = "Sensitivity", tag = "B") +
    # Manual legend using annotate - positioned to avoid cropping
    annotate("segment", x = 0.55, xend = 0.40, y = 0.24, yend = 0.24,
             color = "black", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = 0.38, y = 0.24,
             label = paste0("Reactivity (", auc_103_r, ")"),
             hjust = 0, size = 2.8) +
    annotate("segment", x = 0.55, xend = 0.40, y = 0.16, yend = 0.16,
             color = "#2166AC", linetype = "solid", linewidth = 0.8) +
    annotate("text", x = 0.38, y = 0.16,
             label = paste0("All Predictors (", auc_103_a, ")"),
             hjust = 0, size = 2.8) +
    annotate("point", x = 0.475, y = 0.08, color = "#228B22", size = 3) +
    annotate("text", x = 0.38, y = 0.08, label = "SSEP", hjust = 0, size = 2.8) +
    theme_minimal(base_size = 11, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 9, color = "black"),
      plot.tag = element_text(size = 14, face = "bold"),
      aspect.ratio = 1
    )

  # Combine
  combined <- panel_a + panel_b + plot_layout(ncol = 2)

  return(combined)
}

cat("Clean figure 4 generator loaded.\n\n")
cat("Usage:\n")
cat("fig4 <- generate_figure4_clean(roc_114_react, roc_114_all,\n")
cat("                               roc_103_react, roc_103_all,\n")
cat("                               sep_sens, sep_spec)\n")
cat("print(fig4)\n")
# ==============================================================================
# FIGURE S1: PARTIAL DEPENDENCE PLOTS
# ==============================================================================

generate_figure_s1 <- function(train_data, train_label, params, nrounds, 
                                top_features = NULL) {
  
  cat("Generating Figure S1: Partial Dependence Plots...\n")
  
  # Train final model
  dtrain_full <- xgb.DMatrix(data = as.matrix(train_data), label = train_label, missing = NA)
  set.seed(42)
  final_model <- xgb.train(
    params = params,
    data = dtrain_full,
    nrounds = nrounds,
    verbose = 0
  )
  
  # Get top features if not provided
  if(is.null(top_features)) {
    importance <- xgb.importance(model = final_model)
    top_features <- importance$Feature[1:4]
  }
  
  pdp_list <- list()
  
  for(i in 1:4) {
    raw_name <- top_features[i]
    
    # Clean name
    clean_name <- raw_name
    clean_name <- gsub("\\.Post$", "", clean_name)
    clean_name <- gsub("\\.Post\\.", ".", clean_name)
    clean_name <- gsub("Sound", "Ear", clean_name)
    clean_name <- gsub("\\.", ".", clean_name)
    
    # Check if feature exists
    if(!raw_name %in% names(train_data)) {
      cat("  Warning: Feature", raw_name, "not found. Skipping.\n")
      next
    }
    
    # Grid values
    feat_vals <- train_data[[raw_name]]
    xlims <- quantile(feat_vals, probs = c(0.01, 0.99), na.rm = TRUE)
    grid_vals <- seq(xlims[1], xlims[2], length.out = 50)
    
    # Calculate PDP
    pdp_data <- data.frame()
    
    for(j in seq_along(grid_vals)) {
      temp_data <- train_data
      temp_data[[raw_name]] <- grid_vals[j]
      
      preds_flat <- predict(final_model, as.matrix(temp_data))
      preds <- matrix(preds_flat, ncol = 4, byrow = TRUE)
      
      pdp_data <- bind_rows(pdp_data, data.frame(
        x = grid_vals[j],
        Good = mean(preds[, 1]),
        Intermediate = mean(preds[, 2] + preds[, 3]),
        Poor = mean(preds[, 4])
      ))
    }
    
    # Center data
    pdp_data <- pdp_data %>%
      mutate(
        Good = Good - mean(Good),
        Intermediate = Intermediate - mean(Intermediate),
        Poor = Poor - mean(Poor)
      ) %>%
      pivot_longer(cols = c(Good, Intermediate, Poor),
                   names_to = "Outcome", values_to = "Influence")
    
    pdp_data$Outcome <- factor(pdp_data$Outcome, 
                               levels = c("Good", "Intermediate", "Poor"))
    
    # Create plot
    pdp_list[[i]] <- ggplot(pdp_data, aes(x = x, y = Influence,
                                           color = Outcome, linetype = Outcome)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3) +
      geom_line(linewidth = 0.7) +
      scale_color_manual(
        values = c("Good" = "#1B7837",
                   "Intermediate" = "#E66101",
                   "Poor" = "#B2182B"),
        labels = c("Good (CPC 1)", "Interm. (CPC 2-3)", "Poor (CPC 5)")
      ) +
      scale_linetype_manual(
        values = c("dashed", "solid", "dotdash"),
        labels = c("Good (CPC 1)", "Interm. (CPC 2-3)", "Poor (CPC 5)")
      ) +
      labs(
        x = clean_name,
        y = "Relative Influence",
        title = paste0("Rank ", i, ": ", clean_name)
      ) +
      theme_nature(base_size = 8) +
      theme(
        legend.position = "none",  # No legend in individual panels
        plot.title = element_text(size = 9, face = "bold")
      )
  }
  
  # Combine with patchwork - no legend (will be explained in figure caption)
  fig_s1 <- (pdp_list[[1]] | pdp_list[[2]]) / (pdp_list[[3]] | pdp_list[[4]])

  return(fig_s1)
}


# ==============================================================================
# MAIN EXECUTION FUNCTION
# ==============================================================================

generate_all_figures <- function(output_dir = "eeg-xgboost-open-code/figures") {
  
  cat("\n")
  cat("================================================================\n")
  cat("GENERATING PUBLICATION-QUALITY FIGURES\n")
  cat("================================================================\n\n")
  
  # Check required objects
  required <- c("loo_predictions", "actual_outcomes", "train_data", 
                "train_label", "params", "nrounds", "full")
  
  missing <- required[!sapply(required, exists, envir = .GlobalEnv)]
  if(length(missing) > 0) {
    stop(paste("Missing required objects:", paste(missing, collapse = ", "),
               "\nRun 01_loocv_analysis.R first."))
  }
  
  # Create output directory
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  n_patients <- nrow(full)
  
  # --------------------------------------------------------------------------
  # FIGURE 2
  # --------------------------------------------------------------------------
  cat("\n--- Figure 2 ---\n")
  
  # Check if all_importance exists (from LOO-CV with importance tracking)
  if(exists("all_importance", envir = .GlobalEnv)) {
    fig2 <- generate_figure2(all_importance = all_importance)
  } else {
    fig2 <- generate_figure2(train_data = train_data, train_label = train_label,
                             params = params, nrounds = nrounds)
  }
  
  ggsave(file.path(output_dir, "Figure2_FeatureImportance.pdf"), fig2,
         width = 180, height = 120, units = "mm", dpi = 300)
  ggsave(file.path(output_dir, "Figure2_FeatureImportance.png"), fig2,
         width = 180, height = 120, units = "mm", dpi = 300)
  cat("  Saved: Figure2_FeatureImportance.pdf/png\n")
  
  # --------------------------------------------------------------------------
  # FIGURE 3
  # --------------------------------------------------------------------------
  cat("\n--- Figure 3 ---\n")
  
  fig3 <- generate_figure3(loo_predictions, actual_outcomes, n_patients)
  
  ggsave(file.path(output_dir, "Figure3_ROC_OneVsAll.pdf"), fig3,
         width = 100, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(output_dir, "Figure3_ROC_OneVsAll.png"), fig3,
         width = 100, height = 100, units = "mm", dpi = 300)
  cat("  Saved: Figure3_ROC_OneVsAll.pdf/png\n")
  
  # --------------------------------------------------------------------------
  # FIGURE S1
  # --------------------------------------------------------------------------
  cat("\n--- Figure S1 ---\n")
  
  # Which four features to show.
  #
  # Ranked by mean Gain across the leave-one-out folds - the same basis as
  # Figure 2, so the two figures agree on the order. Ranks 3 and 4 are separated
  # by only a few per cent, so a ranking taken from a single refit is not stable
  # enough to settle their order; the fold means are deterministic and need no
  # refit.
  if (exists("all_importance", envir = .GlobalEnv) &&
      nrow(get("all_importance", envir = .GlobalEnv)) > 0) {

    top_features <- get("all_importance", envir = .GlobalEnv) %>%
      dplyr::select(Feature, Fold, Gain) %>%
      tidyr::complete(Feature, Fold = 1:max(Fold), fill = list(Gain = 0)) %>%
      dplyr::group_by(Feature) %>%
      dplyr::summarise(Mean = mean(Gain), .groups = "drop") %>%
      dplyr::arrange(desc(Mean)) %>%
      head(4) %>%
      dplyr::pull(Feature)

    cat("  Features ranked by mean gain across the leave-one-out folds.\n")

  } else {
    warning("all_importance not found: falling back to a single-model ranking, ",
            "which is not stable across sessions. Run 01_loocv_analysis.R first.")
    dtrain_full <- xgb.DMatrix(data = as.matrix(train_data), label = train_label,
                               missing = NA)
    set.seed(42)
    final_model <- xgb.train(params = params, data = dtrain_full,
                             nrounds = nrounds, verbose = 0)
    top_features <- xgb.importance(model = final_model)$Feature[1:4]
  }

  cat("  Figure S1 features: ", paste(top_features, collapse = ", "), "\n", sep = "")

  fig_s1 <- generate_figure_s1(train_data, train_label, params, nrounds, top_features)
  
  ggsave(file.path(output_dir, "FigureS1_PartialDependence.pdf"), fig_s1,
         width = 180, height = 150, units = "mm", dpi = 300)
  ggsave(file.path(output_dir, "FigureS1_PartialDependence.png"), fig_s1,
         width = 180, height = 150, units = "mm", dpi = 300)
  cat("  Saved: FigureS1_PartialDependence.pdf/png\n")
  
  # --------------------------------------------------------------------------
  # SUMMARY
  # --------------------------------------------------------------------------
  cat("\n================================================================\n")
  cat("FIGURE GENERATION COMPLETE\n")
  cat("================================================================\n")
  cat("\nOutput files in:", normalizePath(output_dir), "\n")
  cat("\nFigure 4 is drawn separately: run 02_roc_reactivity_vs_all_predictors.R\n")
  cat("for the N=114 and N=103 ROC objects, then call generate_figure4_clean().\n")
  
  return(list(fig2 = fig2, fig3 = fig3, fig_s1 = fig_s1))
}


# ==============================================================================
# QUICK RUN (uncomment to execute)
# ==============================================================================
# figures <- generate_all_figures()



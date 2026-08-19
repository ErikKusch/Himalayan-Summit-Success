#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  ANALYSIS & VISUALIZATIONS: Statistical analysis and visualization of
#'      - changes in mean conditions at each peak
#'      - weather stability at each peak
#'  DEPENDENCIES:
#'  - "Himalayan Summit Success - DATA.R"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
source("Himalayan Summit Success - DATA.r")

Dir <- file.path(Dir.Exports, "Gradual Changes in the Himalayan Seasons")
if (!dir.exists(Dir)) dir.create(Dir, recursive = TRUE)

# PREPARATION =============================================================
peaks_ts_df <- Peaks_ts_DAILY_df

trend_summary <- function(df) {
    df <- df %>% filter(!is.na(Year), !is.na(Value))

    out <- tibble::tibble(
        n = nrow(df),
        slope_per_year = NA_real_,
        slope_low_per_year = NA_real_,
        slope_high_per_year = NA_real_,
        slope_per_decade = NA_real_,
        slope_low_per_decade = NA_real_,
        slope_high_per_decade = NA_real_,
        p_value = NA_real_,
        pearson_r = NA_real_
    )

    if (out$n < 3 || dplyr::n_distinct(df$Year) < 2 || dplyr::n_distinct(df$Value) < 2) {
        return(out)
    }

    fit <- lm(Value ~ Year, data = df)
    fit_coef <- summary(fit)$coefficients["Year", ]
    fit_ci <- confint(fit)["Year", ]
    cor_res <- cor.test(df$Year, df$Value, method = "pearson")

    out$slope_per_year <- unname(fit_coef["Estimate"])
    out$slope_low_per_year <- unname(fit_ci[1])
    out$slope_high_per_year <- unname(fit_ci[2])
    out$slope_per_decade <- out$slope_per_year * 10
    out$slope_low_per_decade <- out$slope_low_per_year * 10
    out$slope_high_per_decade <- out$slope_high_per_year * 10
    out$p_value <- unname(fit_coef["Pr(>|t|)"])
    out$pearson_r <- unname(cor_res$estimate)
    out
}

# PANEL A: MEAN CHANGE ====================================================
#' 1. calculate the linear trend for each variable, season and peak for both mean and AR1
mainA_peak_data <- peaks_ts_df %>%
    filter(Variable %in% climate_vars$name) %>%
    filter(PeakID %in% eightks_sf$PKNAME) %>%
    mutate(Year = as.numeric(format(Date, "%Y"))) %>%
    left_join(climate_vars, by = c("Variable" = "name")) %>%
    mutate(
        Variable = ifelse(is.na(display), Variable, display),
        Season = factor(Season,
            levels = c("Pre", "Post"),
            labels = c("Pre-Monsoon", "Post-Monsoon")
        )
    ) %>%
    group_by(Variable, Season, PeakID, Year) %>%
    summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

mainB_peak_data <- pred_df %>%
    filter(Variable %in% climate_vars$name) %>%
    filter(PeakID %in% eightks_sf$PKNAME) %>%
    left_join(climate_vars, by = c("Variable" = "name")) %>%
    mutate(
        Season = factor(Season,
            levels = c("Pre", "Post"),
            labels = c("Pre-Monsoon", "Post-Monsoon")
        ),
        Variable = ifelse(is.na(display), Variable, display)
    ) %>%
    group_by(Variable, Season, PeakID, Year) %>%
    summarise(Value = mean(AR1, na.rm = TRUE), .groups = "drop")

mainAB_trend_summary <- dplyr::bind_rows(
    mainA_peak_data %>% mutate(Metric = "Mean"),
    mainB_peak_data %>% mutate(Metric = "AR1")
) %>%
    group_by(Variable, Season, PeakID, Metric) %>%
    group_modify(~ trend_summary(.x)) %>%
    ungroup()

#' 2. summarise for all peaks to get an overall trend for each variable and season, for both mean and AR1, insert PeakID = "Overall" instead of original PeakID
mainAB_trend_summaryOverall <- mainAB_trend_summary %>%
    mutate(PeakID = "Overall") %>%
    group_by(Variable, Season, PeakID, Metric) %>%
    summarise(
        n = sum(n),
        slope_per_year = mean(slope_per_year, na.rm = TRUE),
        slope_low_per_year = mean(slope_low_per_year, na.rm = TRUE),
        slope_high_per_year = mean(slope_high_per_year, na.rm = TRUE),
        slope_per_decade = mean(slope_per_decade, na.rm = TRUE),
        slope_low_per_decade = mean(slope_low_per_decade, na.rm = TRUE),
        slope_high_per_decade = mean(slope_high_per_decade, na.rm = TRUE),
        p_value = mean(p_value, na.rm = TRUE),
        pearson_r = mean(pearson_r, na.rm = TRUE),
        .groups = "drop"
    )
trends_df <- rbind(mainAB_trend_summaryOverall, mainAB_trend_summary) %>%
    arrange(Variable, Season, Metric, PeakID)
## rename Mean to "Mean Conditions" and AR1 to "Weather Stability" in the Metric column
trends_df$Metric <- recode(trends_df$Metric, "Mean" = "Mean Conditions", "AR1" = "Weather Stability")
write.csv(
    trends_df,
    file = file.path(Dir, "trend_summary_per_peak.csv"),
    row.names = FALSE
)

#' 3. make plots of correlation estimate:
## main text plot, just overall trend with error bars, no individual peaks
Mean_mainBoxplot_gg <- ggplot(
    trends_df %>% filter(PeakID != "Overall"),
    aes(x = Variable, y = slope_per_decade, fill = Season)
) +
    geom_boxplot() +
    facet_grid2(factor(Metric, levels = c("Mean Conditions", "Weather Stability")) ~ factor(Variable, levels = climate_vars$display), scales = "free", independent = "y") +
    scale_fill_manual(values = c("Pre-Monsoon" = PreColour, "Post-Monsoon" = PostColour)) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
    labs(
        x = "Variable", y = "Slope per Decade",
        fill = "Season"
    ) +
    theme_bw() +
    theme(legend.position = "top")

ggsave(
    plot = Mean_mainBoxplot_gg,
    filename = file.path(Dir, "Mean_mainBoxplot.png"),
    width = 24, height = 16, units = "cm", dpi = 300
)

Mean_suppBoxplot_gg <- ggplot(
    trends_df %>% filter(PeakID != "Overall"),
    aes(x = PeakID, y = slope_per_decade, color = Season)
) +
    geom_point(size = 3) +
    facet_grid2(factor(Metric, levels = c("Mean Conditions", "Weather Stability")) ~ factor(Variable, levels = climate_vars$display), scales = "free", independent = "y") +
    scale_color_manual(values = c("Pre-Monsoon" = PreColour, "Post-Monsoon" = PostColour)) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
    labs(
        x = "Peaks", y = "Slope per Decade",
        color = "Season"
    ) +
    theme_bw() +
    theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

ggsave(
    plot = Mean_suppBoxplot_gg,
    filename = file.path(Dir, "Mean_suppBoxplot.png"),
    width = 24, height = 18, units = "cm", dpi = 300
)

# PANEL B: SD CHANGE ======================================================
sdA_peak_data <- peaks_ts_df %>%
    filter(Variable %in% climate_vars$name) %>%
    filter(PeakID %in% eightks_sf$PKNAME) %>%
    mutate(Year = as.numeric(format(Date, "%Y"))) %>%
    left_join(climate_vars, by = c("Variable" = "name")) %>%
    mutate(
        Variable = ifelse(is.na(display), Variable, display),
        Season = factor(Season,
            levels = c("Pre", "Post"),
            labels = c("Pre-Monsoon", "Post-Monsoon")
        )
    ) %>%
    group_by(Variable, Season, PeakID, Year) %>%
    summarise(Value = sd(Value, na.rm = TRUE), .groups = "drop")

sdB_peak_data <- pred_df %>%
    filter(Variable %in% climate_vars$name) %>%
    filter(PeakID %in% eightks_sf$PKNAME) %>%
    left_join(climate_vars, by = c("Variable" = "name")) %>%
    mutate(
        Season = factor(Season,
            levels = c("Pre", "Post"),
            labels = c("Pre-Monsoon", "Post-Monsoon")
        ),
        Variable = ifelse(is.na(display), Variable, display)
    ) %>%
    group_by(Variable, Season, PeakID, Year) %>%
    summarise(Value = mean(AR1, na.rm = TRUE), .groups = "drop")

sd_trend_summary <- dplyr::bind_rows(
    sdA_peak_data %>% mutate(Metric = "SD"),
    sdB_peak_data %>% mutate(Metric = "AR1")
) %>%
    group_by(Variable, Season, PeakID, Metric) %>%
    group_modify(~ trend_summary(.x)) %>%
    ungroup()

sd_trend_summaryOverall <- sd_trend_summary %>%
    mutate(PeakID = "Overall") %>%
    group_by(Variable, Season, PeakID, Metric) %>%
    summarise(
        n = sum(n),
        slope_per_year = mean(slope_per_year, na.rm = TRUE),
        slope_low_per_year = mean(slope_low_per_year, na.rm = TRUE),
        slope_high_per_year = mean(slope_high_per_year, na.rm = TRUE),
        slope_per_decade = mean(slope_per_decade, na.rm = TRUE),
        slope_low_per_decade = mean(slope_low_per_decade, na.rm = TRUE),
        slope_high_per_decade = mean(slope_high_per_decade, na.rm = TRUE),
        p_value = mean(p_value, na.rm = TRUE),
        pearson_r = mean(pearson_r, na.rm = TRUE),
        .groups = "drop"
    )

sd_trends_df <- rbind(sd_trend_summaryOverall, sd_trend_summary) %>%
    arrange(Variable, Season, Metric, PeakID)
sd_trends_df$Metric <- recode(sd_trends_df$Metric, "SD" = "Standard Deviation", "AR1" = "Weather Stability")

write.csv(
    sd_trends_df,
    file = file.path(Dir, "trend_summary_sd_per_peak.csv"),
    row.names = FALSE
)

SD_mainBoxplot_gg <- ggplot(
    sd_trends_df %>% filter(PeakID != "Overall"),
    aes(x = Variable, y = slope_per_decade, fill = Season)
) +
    geom_boxplot() +
    facet_grid2(factor(Metric, levels = c("Standard Deviation", "Weather Stability")) ~ factor(Variable, levels = climate_vars$display), scales = "free", independent = "y") +
    scale_fill_manual(values = c("Pre-Monsoon" = PreColour, "Post-Monsoon" = PostColour)) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
    labs(
        x = "Variable", y = "Slope per Decade",
        fill = "Season"
    ) +
    theme_bw() +
    theme(legend.position = "top")

SD_suppBoxplot_gg <- ggplot(
    sd_trends_df %>% filter(PeakID != "Overall"),
    aes(x = PeakID, y = slope_per_decade, color = Season)
) +
    geom_point(size = 3) +
    facet_grid2(factor(Metric, levels = c("Standard Deviation", "Weather Stability")) ~ factor(Variable, levels = climate_vars$display), scales = "free", independent = "y") +
    scale_color_manual(values = c("Pre-Monsoon" = PreColour, "Post-Monsoon" = PostColour)) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "grey") +
    labs(
        x = "Peaks", y = "Slope per Decade",
        color = "Season"
    ) +
    theme_bw() +
    theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )

ggsave(
    plot = plot_grid(
        SD_mainBoxplot_gg,
        SD_suppBoxplot_gg + theme(legend.position = "bottom"),
        ncol = 1,
        rel_heights = c(1, 1.2), labels = c("A", "B")
    ),
    filename = file.path(Dir, "SD_suppBoxplot.png"),
    width = 24, height = 30, units = "cm", dpi = 300
)

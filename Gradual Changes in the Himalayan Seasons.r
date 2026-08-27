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

# PREPARATION =============================================================
peaks_ts_df <- Peaks_ts_DAILY_df

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

# LINEAR TREND SUMMARY ====================================================
## lmerTest::lmer is used in place of lme4::lmer so that summary() includes
## Satterthwaite p-values in the coefficients table (the "Pr(>|t|)" column).
fit_lmer_trend <- function(peak_data, metric_label) {
    do.call(rbind, lapply(unique(peak_data$Variable), FUN = function(y) {
        do.call(rbind, lapply(unique(peak_data$Season), FUN = function(x) {
            df <- peak_data %>%
                filter(Variable == y, Season == x) %>%
                mutate(Year_c = Year - min(Year), PeakID = factor(PeakID))

            m1 <- lmerTest::lmer(Value ~ Year_c + (1 | PeakID), data = df)
            ci   <- confint(m1, parm = "Year_c", level = 0.95)
            slp  <- unname(fixef(m1)["Year_c"])

            tibble::tibble(
                Variable           = y,
                Season             = x,
                PeakID             = "Overall",
                Metric             = metric_label,
                n                  = nrow(df),
                slope_per_year     = slp,
                slope_low_per_year = ci[1],
                slope_high_per_year= ci[2],
                slope_per_decade   = slp * 10,
                slope_low_per_decade  = ci[1] * 10,
                slope_high_per_decade = ci[2] * 10,
                p_value            = summary(m1)$coefficients["Year_c", "Pr(>|t|)"],
                pearson_r          = NA_real_
            )
        }))
    }))
}

MeanTrendAll <- fit_lmer_trend(mainA_peak_data, "Mean")
ARTrendAll   <- fit_lmer_trend(mainB_peak_data, "AR1")

trends_df <- rbind(MeanTrendAll, ARTrendAll, mainAB_trend_summary)
## rename Mean to "Mean Conditions" and AR1 to "Weather Stability" in the Metric column
trends_df$Metric <- recode(trends_df$Metric, "Mean" = "Mean Conditions", "AR1" = "Weather Stability")
## saving file making sure it is properly delimited and that values are not rounded
write.csv2(
    trends_df,
    file = file.path(Dir, "trend_summary_per_peak.csv"),
    row.names = FALSE
)

# PLOTTING ================================================================
## shared margin geometry constants
seg_gap <- 1 # gap between last data year and start of the margin
seg_len <- 1.6 # length of each dashed segment

## helper: build the lm-endpoint annotation data for a peak data frame
make_trend_annotations <- function(peak_data) {
    max_year <- max(peak_data$Year, na.rm = TRUE)
    ## Pre- and Post-Monsoon get horizontally offset columns so their
    ## arrows and Δ labels never overlap even when y-values are close.
    ## Pre: left column,  Post: right column (each seg_len wide, seg_gap apart)
    x_start_pre  <- max_year + seg_gap
    x_end_pre    <- x_start_pre + seg_len
    x_mid_pre    <- (x_start_pre + x_end_pre) / 2
    x_start_post <- x_end_pre + 1.2
    x_end_post   <- x_start_post + seg_len
    x_mid_post   <- (x_start_post + x_end_post) / 2

    peak_data %>%
        group_by(Variable, Season) %>%
        group_modify(~ {
            fit <- lm(Value ~ Year, data = .x)
            pred <- predict(fit, newdata = data.frame(Year = c(min(.x$Year), max(.x$Year))))
            tibble::tibble(trend_start = pred[1], trend_end = pred[2], diff = pred[2] - pred[1])
        }) %>%
        ungroup() %>%
        mutate(
            x_start = ifelse(Season == "Pre-Monsoon", x_start_pre,  x_start_post),
            x_end   = ifelse(Season == "Pre-Monsoon", x_end_pre,    x_end_post),
            x_mid   = ifelse(Season == "Pre-Monsoon", x_mid_pre,    x_mid_post)
        )
}

## helper: produce the full annotated facet plot for one metric
## y_max:      optional upper y-axis cap (e.g. 1 for AR1)
## diff_digits: decimal places in the Δ label
plot_metric <- function(peak_data, y_label, y_max = NA, diff_digits = 2) {
    annotations <- make_trend_annotations(peak_data)
    max_year <- max(peak_data$Year)
    diff_fmt <- paste0("Δ %+.", diff_digits, "f")

    ggplot(peak_data, aes(x = Year, y = Value, fill = Season, color = Season)) +
        geom_smooth(se = TRUE, span = 0.25) +
        geom_smooth(method = "lm", linetype = "dashed", se = FALSE) +
        geom_segment(
            data = annotations,
            aes(x = x_start, xend = x_end, y = trend_start, yend = trend_start),
            linetype = "dashed", linewidth = 0.6, inherit.aes = TRUE
        ) +
        geom_segment(
            data = annotations,
            aes(x = x_start, xend = x_end, y = trend_end, yend = trend_end),
            linetype = "dashed", linewidth = 0.6, inherit.aes = TRUE
        ) +
        geom_segment(
            data = annotations,
            aes(x = x_mid, xend = x_mid, y = trend_start, yend = trend_end),
            arrow = arrow(length = unit(0.15, "cm"), ends = "both"), linewidth = 0.5,
            inherit.aes = TRUE
        ) +
        geom_text(
            data = annotations %>% filter(Season == "Post-Monsoon"),
            aes(x = x_mid, y = pmax(trend_start, trend_end), label = sprintf(diff_fmt, diff)),
            vjust = -0.6, size = 3, fontface = "bold", show.legend = FALSE, inherit.aes = TRUE
        ) +
        geom_text(
            data = annotations %>% filter(Season == "Pre-Monsoon"),
            aes(x = x_mid, y = pmin(trend_start, trend_end), label = sprintf(diff_fmt, diff)),
            vjust = 1.6, size = 3, fontface = "bold", show.legend = FALSE, inherit.aes = TRUE
        ) +
        scale_fill_manual(values = c("Pre-Monsoon" = PreColour, "Post-Monsoon" = PostColour)) +
        scale_color_manual(values = c("Pre-Monsoon" = PreColour, "Post-Monsoon" = PostColour)) +
        scale_x_continuous(
            expand = expansion(mult = c(0.02, 0), add = c(0, seg_gap + 2 * seg_len + 1.2 + seg_len + 1)),
            breaks = scales::breaks_pretty()(range(peak_data$Year)),
            labels = function(x) ifelse(x >= min(peak_data$Year) & x <= max_year, x, "")
        ) +
        scale_y_continuous(limits = c(NA, y_max)) +
        facet_wrap(~ factor(Variable, levels = climate_vars$display), scales = "free_y") +
        theme_bw() +
        labs(x = "Year", y = y_label) +
        theme(
            legend.position  = "bottom",
            legend.direction = "horizontal",
            legend.text      = element_text(size = 12)
        )
}

ggsave(
    plot = plot_metric(mainA_peak_data, y_label = "Mean Value"),
    filename = file.path(Dir, "Figure_MeanTrend.png"),
    width = 32, height = 19, units = "cm", dpi = 300
)

ggsave(
    plot = plot_metric(mainB_peak_data, y_label = "AR1", diff_digits = 3),
    filename = file.path(Dir, "SUPP_AR1Trend.png"),
    width = 32, height = 19, units = "cm", dpi = 300
)

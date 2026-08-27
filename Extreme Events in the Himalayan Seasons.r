#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  ANALYSIS & VISUALIZATIONS: Statistical analysis and visualization of
#'      - number of extreme events (HIGH / LOW) per season per year
#'      - average run-length of extreme events (HIGH / LOW) per season per year
#'  DEPENDENCIES:
#'  - "Himalayan Summit Success - DATA.R"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
source("Himalayan Summit Success - DATA.r")

Dir <- file.path(Dir.Exports, "Extreme Events in the Himalayan Seasons")
if (!dir.exists(Dir)) dir.create(Dir, recursive = TRUE)

# PREPARATION =============================================================
peaks_ts_df <- Peaks_ts_DAILY_df

## restrict to the same peaks and variables used in the main script
extreme_base <- peaks_ts_df %>%
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
    )

# METRIC 1: Number of extremes per season per year ========================
extreme_counts <- extreme_base %>%
    filter(Extreme %in% c("HIGH", "LOW")) %>%
    group_by(Variable, Season, PeakID, Year, Extreme) %>%
    summarise(Value = n(), .groups = "drop")

count_HIGH <- extreme_counts %>%
    filter(Extreme == "HIGH") %>%
    select(-Extreme)
count_LOW <- extreme_counts %>%
    filter(Extreme == "LOW") %>%
    select(-Extreme)

# METRIC 2: Average run-length of extremes per season per year ============
## Run-length = consecutive days with Extreme == target label within a
## Variable × Season × PeakID × Year group. We compute all runs, then
## average them per group × year.
avg_runlength <- function(df, target) {
    df %>%
        arrange(Variable, Season, PeakID, Year, Date) %>%
        group_by(Variable, Season, PeakID, Year) %>%
        group_modify(~ {
            is_extreme <- .x$Extreme == target
            ## rle gives lengths of consecutive runs; keep only extreme runs
            r <- rle(is_extreme)
            run_lengths <- r$lengths[r$values]
            tibble::tibble(Value = if (length(run_lengths) == 0) NA_real_ else mean(run_lengths))
        }) %>%
        ungroup()
}

runlen_HIGH <- avg_runlength(extreme_base, "HIGH")
runlen_LOW <- avg_runlength(extreme_base, "LOW")

# TREND SUMMARY (lm per peak, then lmerTest mixed model overall) ==========
trend_summary <- function(df) {
    df <- df %>% filter(!is.na(Year), !is.na(Value), !is.nan(Value))

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

## per-peak lm trends (mirrors mainAB_trend_summary in the main script)
extreme_trend_summary <- dplyr::bind_rows(
    count_HIGH %>% mutate(Metric = "N Extremes HIGH"),
    count_LOW %>% mutate(Metric = "N Extremes LOW"),
    runlen_HIGH %>% mutate(Metric = "Run Length HIGH"),
    runlen_LOW %>% mutate(Metric = "Run Length LOW")
) %>%
    group_by(Variable, Season, PeakID, Metric) %>%
    group_modify(~ trend_summary(.x)) %>%
    ungroup()

## mixed-model overall trend (mirrors fit_lmer_trend in the main script)
## lmerTest::lmer provides Satterthwaite p-values in summary()$coefficients
fit_lmer_extreme <- function(peak_data, metric_label) {
    do.call(rbind, lapply(unique(peak_data$Variable), FUN = function(y) {
        do.call(rbind, lapply(unique(peak_data$Season), FUN = function(x) {
            df <- peak_data %>%
                filter(Variable == y, Season == x) %>%
                filter(!is.na(Value), !is.nan(Value)) %>%
                mutate(Year_c = Year - min(Year, na.rm = TRUE), PeakID = factor(PeakID))

            ## need at least 2 peaks (groups) and sufficient observations
            if (dplyr::n_distinct(df$PeakID) < 2 || nrow(df) < 5 ||
                dplyr::n_distinct(df$Year) < 2) {
                return(tibble::tibble(
                    Variable = y, Season = x, PeakID = "Overall",
                    Metric = metric_label, n = nrow(df),
                    slope_per_year = NA_real_, slope_low_per_year = NA_real_,
                    slope_high_per_year = NA_real_, slope_per_decade = NA_real_,
                    slope_low_per_decade = NA_real_, slope_high_per_decade = NA_real_,
                    p_value = NA_real_, pearson_r = NA_real_
                ))
            }

            m1 <- lmerTest::lmer(Value ~ Year_c + (1 | PeakID), data = df)
            ci <- confint(m1, parm = "Year_c", level = 0.95)
            slp <- unname(fixef(m1)["Year_c"])

            tibble::tibble(
                Variable = y,
                Season = x,
                PeakID = "Overall",
                Metric = metric_label,
                n = nrow(df),
                slope_per_year = slp,
                slope_low_per_year = ci[1],
                slope_high_per_year = ci[2],
                slope_per_decade = slp * 10,
                slope_low_per_decade = ci[1] * 10,
                slope_high_per_decade = ci[2] * 10,
                p_value = summary(m1)$coefficients["Year_c", "Pr(>|t|)"],
                pearson_r = NA_real_
            )
        }))
    }))
}

NHighTrendAll <- fit_lmer_extreme(count_HIGH, "N Extremes HIGH")
NLowTrendAll <- fit_lmer_extreme(count_LOW, "N Extremes LOW")
RLHighTrendAll <- fit_lmer_extreme(runlen_HIGH, "Run Length HIGH")
RLLowTrendAll <- fit_lmer_extreme(runlen_LOW, "Run Length LOW")

trends_extreme_df <- rbind(
    NHighTrendAll, NLowTrendAll, RLHighTrendAll, RLLowTrendAll,
    extreme_trend_summary
)

write.csv2(
    trends_extreme_df,
    file = file.path(Dir, "trend_summary_extremes_per_peak.csv"),
    row.names = FALSE
)

# PLOTTING ================================================================
seg_gap <- 1 # gap between last data year and start of the margin
seg_len <- 1.6 # length of each dashed segment

make_trend_annotations <- function(peak_data) {
    max_year <- max(peak_data$Year, na.rm = TRUE)
    ## Pre- and Post-Monsoon get horizontally offset columns so their
    ## arrows and Δ labels never overlap even when y-values are close.
    ## Pre: left column,  Post: right column (each seg_len wide, seg_gap apart)
    x_start_pre <- max_year + seg_gap
    x_end_pre <- x_start_pre + seg_len
    x_mid_pre <- (x_start_pre + x_end_pre) / 2
    x_start_post <- x_end_pre + 1.2
    x_end_post <- x_start_post + seg_len
    x_mid_post <- (x_start_post + x_end_post) / 2

    peak_data %>%
        group_by(Variable, Season) %>%
        group_modify(~ {
            d <- .x %>% filter(!is.na(Value), !is.nan(Value))
            if (nrow(d) < 2 || dplyr::n_distinct(d$Year) < 2) {
                return(tibble::tibble(trend_start = NA_real_, trend_end = NA_real_, diff = NA_real_))
            }
            fit <- lm(Value ~ Year, data = d)
            pred <- predict(fit, newdata = data.frame(Year = c(min(d$Year), max(d$Year))))
            tibble::tibble(trend_start = pred[1], trend_end = pred[2], diff = pred[2] - pred[1])
        }) %>%
        ungroup() %>%
        mutate(
            x_start = ifelse(Season == "Pre-Monsoon", x_start_pre, x_start_post),
            x_end   = ifelse(Season == "Pre-Monsoon", x_end_pre, x_end_post),
            x_mid   = ifelse(Season == "Pre-Monsoon", x_mid_pre, x_mid_post)
        )
}

## y_max:       optional upper y-axis cap
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
    plot = plot_grid(
        plot_metric(count_HIGH, y_label = "N Extreme Days (HIGH)") +
            theme(legend.position = "none"),
        plot_metric(count_LOW, y_label = "N Extreme Days (LOW)"),
        labels = c("A", "B"), label_size = 16, ncol = 1, align = "v"
    ),
    filename = file.path(Dir, "SUPP_NExtremes.png"),
    width = 32, height = 26, units = "cm", dpi = 300
)

ggsave(
    plot = plot_grid(
        plot_metric(runlen_HIGH, y_label = "Mean Run Length (HIGH) [days]") +
            theme(legend.position = "none"),
        plot_metric(runlen_LOW, y_label = "Mean Run Length (LOW) [days]"),
        labels = c("A", "B"), label_size = 16, ncol = 1, align = "v"
    ),
    filename = file.path(Dir, "SUPP_RunLength.png"),
    width = 32, height = 26, units = "cm", dpi = 300
)


#' this here explains why there are no trendlines for pre-monsoon season high snow cover: the upper bound is simply almost the same as the maximum value of the data
var_df <- peaks_ts_df[peaks_ts_df$Variable == "snow_cover", ]
do.call(rbind, lapply(eightks_sf$PKNAME, FUN = function(peak) {
    data.frame(
        Variable = "snow_cover",
        Peak = peak,
        Season = rep(c("Pre", "Post"), each = 2),
        Bound = rep(c("Lower", "Upper"), 2),
        Value = c(
            quantile(var_df[var_df$Season == "Pre" & var_df$PeakID == peak, "Value"], probs = c(0.05, 0.95), na.rm = TRUE),
            quantile(var_df[var_df$Season == "Post" & var_df$PeakID == peak, "Value"], probs = c(0.05, 0.95), na.rm = TRUE)
        )
    )
}))

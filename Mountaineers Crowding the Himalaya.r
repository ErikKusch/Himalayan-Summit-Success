#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  ANALYSIS & VISUALIZATIONS:
#'      - Expedition number and expedition size over time
#'      - Summit bid windows over time
#'      - Change-point analysis of mortality and summit success
#'  DEPENDENCIES:
#'  - "Himalayan Summit Success - DATA.r"
#'  - "PlottingFunctions.r"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
# source("Himalayan Summit Success - DATA.r")
source("PlottingFunctions.r")
Dir <- file.path(Dir.Exports, "Crowding in Himalayan mountaineering")
if (!dir.exists(Dir)) dir.create(Dir, recursive = TRUE)

# BAYESIAN MODELS =========================================================
## Data -------------------------------------------------------------------
df_data <- Expeditions_df
df_data$year_0 <- df_data$YEAR - min(df_data$YEAR) # add year 0 as starting year
df_data <- df_data %>%
    filter(PEAKID %in% TargetIDs)
breaks <- seq(Years_vec[1], tail(Years_vec, 1), by = 10)
df_data$YearBin <- cut(
    df_data$YEAR,
    breaks = breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = paste0(breaks[-length(breaks)], ":", breaks[-1])
)

## Total Number of Expeditions --------------------------------------------
df_expeditions <- df_data %>%
    group_by(PEAKID, year_0) %>%
    summarise(
        n_expeditions = n(),
        .groups = "drop"
    )
if (file.exists(file.path(Dir.Exports, "model_TE.RData"))) {
    load(file.path(Dir.Exports, "model_TE.RData"))
} else {
    model_TE <- brms::brm(
        n_expeditions ~ year_0 + year_0:PEAKID + (1 | PEAKID),
        data = df_expeditions,
        family = poisson(),
        chains = 4,
        cores = 4,
        iter = 10000,
        warmup = 5000
    )
    save(model_TE, file = file.path(Dir.Exports, "model_TE.RData"))
}
# summary(model_TE)
# pp_check(model_TE, ndraws = 500)
conditions <- expand.grid(
    PEAKID = sort(unique(df_expeditions$PEAKID)),
    year_0 = sort(unique(df_expeditions$year_0))
)
rownames(conditions) <- NULL
TE_draws <- model_TE %>%
    add_epred_draws(newdata = conditions, re_formula = NULL) %>%
    rename(TE = .epred)

## Total Expedition Members -----------------------------------------------
df_members <- df_data %>%
    group_by(PEAKID, year_0) %>%
    summarise(
        total_members = sum(TOTMEMBERS, na.rm = TRUE),
        .groups = "drop"
    )
if (file.exists(file.path(Dir.Exports, "model_TM.RData"))) {
    load(file.path(Dir.Exports, "model_TM.RData"))
} else {
    model_TM <- brms::brm(
        total_members ~ year_0 + year_0:PEAKID + (1 | PEAKID),
        data = df_members,
        family = poisson(),
        chains = 4,
        cores = 4,
        iter = 10000,
        warmup = 5000
    )
    save(model_TM, file = file.path(Dir.Exports, "model_TM.RData"))
}
# summary(model_TM)
# pp_check(model_TM, ndraws = 500)
conditions <- expand.grid(
    PEAKID = sort(unique(df_members$PEAKID)),
    year_0 = sort(unique(df_members$year_0))
)
rownames(conditions) <- NULL
TM_draws <- model_TM %>%
    add_epred_draws(newdata = conditions, re_formula = NULL) %>%
    rename(TM = .epred)

## Expedition Size --------------------------------------------------------
df_filtered <- df_data
if (file.exists(file.path(Dir.Exports, "model_ES.RData"))) {
    load(file.path(Dir.Exports, "model_ES.RData"))
} else {
    model_ES <- brms::brm(
        TOTMEMBERS ~ year_0 +
            +year_0:PEAKID +
            (1 | PEAKID),
        data = df_filtered,
        family = poisson(),
        chains = 4,
        cores = 4,
        # cores = parallel::detectCores(),
        iter = 10000,
        warmup = 5000
    )
    save(model_ES, file = file.path(Dir.Exports, "model_ES.RData"))
}
# summary(model_ES)
# pp_check(model_ES, ndraws = 500)
conditions <- expand.grid(
    PEAKID = sort(unique(df_filtered$PEAKID)),
    year_0 = sort(unique(df_filtered$year_0))
)
rownames(conditions) <- NULL
ES_draws <- model_ES %>%
    add_epred_draws(newdata = conditions, re_formula = NULL) %>%
    rename(ES = .epred)

## Summit Bid Window ------------------------------------------------------
df_filtered <- df_data %>%
    filter(DIFF_BCtoSMT >= 0)
if (file.exists(file.path(Dir.Exports, "model_SB.RData"))) {
    load(file.path(Dir.Exports, "model_SB.RData"))
} else {
    model_SB <- brms::brm(
        DIFF_BCtoSMT ~ year_0 +
            +year_0:PEAKID +
            (1 | PEAKID),
        data = df_filtered,
        # family = "gaussian",
        family = poisson(),
        chains = 4,
        cores = 4,
        # cores = parallel::detectCores(),
        iter = 10000,
        warmup = 5000
    )
    save(model_SB, file = file.path(Dir.Exports, "model_SB.RData"))
}
# summary(model_SB)
# pp_check(model_SB, ndraws = 500)
conditions <- expand.grid(
    PEAKID = sort(unique(df_filtered$PEAKID)),
    year_0 = sort(unique(df_filtered$year_0))
)
rownames(conditions) <- NULL
SB_draws <- model_SB %>%
    add_epred_draws(newdata = conditions, re_formula = NULL) %>%
    rename(SB = .epred)

# CHANGEPOINT ANALYSIS ====================================================
## Data -------------------------------------------------------------------
SeasonCounts_df <- df_data %>%
    count(YEAR, SEASON) %>%
    pivot_wider(
        names_from = SEASON,
        values_from = n,
        names_prefix = "Season_",
        values_fill = 0
    )
SeasonCounts_df$TOTAL <- SeasonCounts_df$Season_1 + SeasonCounts_df$Season_3
SeasonCounts_df$PROPORTION <- SeasonCounts_df$Season_3 / SeasonCounts_df$TOTAL
clean_df <- SeasonCounts_df[!is.na(SeasonCounts_df$PROPORTION), ]
data_ranks <- rank(clean_df$PROPORTION)

## Analysis -------------------------------------------------------------------
res_amoc <- cpt.mean(data_ranks, method = "AMOC")
cp_index <- cpts(res_amoc)
if (length(cp_index) > 0) {
    # Translate the index back to the actual year from your dataset
    cp_year <- clean_df$YEAR[cp_index]

    # Calculate medians of the original proportions for reporting
    median_before <- median(clean_df$PROPORTION[1:cp_index])
    median_after <- median(clean_df$PROPORTION[(cp_index + 1):nrow(clean_df)])

    cat("--- Median-based Change Point Analysis ---\n")
    cat("Significant change point found in year:", cp_year, "\n")
    cat("Median up to (and including)", cp_year, ":", round(median_before, 4), "\n")
    cat("Median after", cp_year, ":", round(median_after, 4), "\n")
} else {
    cat("No significant change point found in the median.\n")
}

## Plotting ---------------------------------------------------------------
ShadeSteps <- 150
UpperBreaks <- seq(0.5, 1, length.out = ShadeSteps + 1)
LowerBreaks <- seq(0, 0.5, length.out = ShadeSteps + 1)

UpperShade_df <- data.frame(
    ymin = UpperBreaks[-(ShadeSteps + 1)],
    ymax = UpperBreaks[-1],
    alpha = seq(0, 0.28, length.out = ShadeSteps)
)

LowerShade_df <- data.frame(
    ymin = LowerBreaks[-(ShadeSteps + 1)],
    ymax = LowerBreaks[-1],
    alpha = rev(seq(0, 0.28, length.out = ShadeSteps))
)

Plots_Changepoint <- ggplot(
    SeasonCounts_df,
    aes(x = YEAR, y = PROPORTION)
) +
    stat_smooth(
        data = SeasonCounts_df %>% filter(YEAR <= cp_year),
        method = "lm",
        color = "black", # PreColour,
        linewidth = 1,
        level = 0.75
    ) +
    stat_smooth(
        data = SeasonCounts_df %>% filter(YEAR >= cp_year),
        method = "lm",
        color = "black", # PostColour,
        linewidth = 1,
        level = 0.75
    ) +
    geom_rect(
        data = LowerShade_df,
        aes(ymin = ymin, ymax = ymax, alpha = alpha),
        xmin = -Inf,
        xmax = Inf,
        fill = PreColour,
        inherit.aes = FALSE
    ) +
    geom_rect(
        data = UpperShade_df,
        aes(ymin = ymin, ymax = ymax, alpha = alpha),
        xmin = -Inf,
        xmax = Inf,
        fill = PostColour,
        inherit.aes = FALSE
    ) +
    scale_alpha_identity() +
    geom_point() +
    # scale_color_gradient2(
    #     low = PreColour,
    #     mid = "#858585",
    #     high = PostColour,
    #     midpoint = 0.5,
    #     name = "Proportion of Post-Monsoon Season Expeditions"
    # ) +
    # geom_line() +
    # Add horizontal segments for medians (equivalent to segments())
    # geom_segment(
    #     aes(x = min(clean_df$YEAR), xend = cp_year, y = median_before, yend = median_before),
    #     color = "black", linewidth = 1.5 # Note: linewidth instead of lwd for ggplot2
    # ) +
    # geom_segment(
    #     aes(x = cp_year, xend = max(clean_df$YEAR), y = median_after, yend = median_after),
    #     color = "black", linewidth = 1.5
    # ) +
    # Add vertical line at change point (equivalent to abline(v = cp_year))
    geom_vline(xintercept = cp_year, color = "#5a0000", linewidth = 1, linetype = "dashed") +
    # add horizontal, dashed line at 0.5
    geom_hline(yintercept = 0.5, color = "#858585", linewidth = 0.2, linetype = "dashed") +
    # Add label for change point (equivalent to text())
    annotate(
        "label",
        x = cp_year, y = 0.125,
        label = paste("Changepoint:", cp_year),
        hjust = -0.1,
        color = "#5a0000"
    ) +
    ggrepel::geom_label_repel(
        aes(label = TOTAL),
        box.padding = 0.35,
        point.padding = 0.5,
        segment.color = "grey50"
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(y = "Proportion of Post-Monsoon Season Expeditions", x = "Year") +
    lims(y = c(0, 1))
# Plots_Changepoint

# PLOT SAVING =============================================================
## Bayesian Plot Making ---------------------------------------------------
Plots_ExpedSizes <- FUN.BayesianPlot(
    plot_ls = list(
        df = left_join(
            ES_draws %>% select(PEAKID, year_0, .draw, ES),
            TE_draws %>% select(PEAKID, year_0, .draw, TE),
            by = c("PEAKID", "year_0", ".draw")
        ),
        cols = c("#003575", "#754000"),
        columns = c("ES", "TE"),
        names = c("Individual Expedition Sizes [#]", "Total Number of Expeditions [#]")
    ),
    TotalExped = TRUE
)

Plots_ExpedSummitBid <- FUN.BayesianPlot(
    plot_ls = list(
        df = left_join(
            SB_draws %>% select(PEAKID, year_0, .draw, SB),
            TM_draws %>% select(PEAKID, year_0, .draw, TM),
            by = c("PEAKID", "year_0", ".draw")
        ),
        cols = c("#6F0075", "#067500"),
        columns = c("SB", "TM"),
        names = c("Summit Bid Window [days]", "Total Expedition Sizes [#]")
    ),
    TotalExped = TRUE
)
Plots_ExpedSummitBid$Main

## Main Text --------------------------------------------------------------
Main_gg <- cowplot::plot_grid(
    plot_grid(
        Plots_ExpedSummitBid$Main,
        Plots_ExpedSizes$Main,
        nrow = 1, align = "h", axis = "tb", labels = c("A", "B")
    ),
    Plots_Changepoint,
    ncol = 1, labels = c("", "C"), rel_heights = c(0.85, 1)
)
# Main_gg
ggsave(
    Main_gg,
    filename = file.path(Dir, "Figure_Mountaineers_Crowding_the_Himalaya.png"),
    width = 36/1.25, height = 34/1.25, units = "cm", dpi = 600
)

## Supplement -------------------------------------------------------------
ggsave(
    Plots_ExpedSummitBid$Supp,
    filename = file.path(Dir, "SUPP_Mountaineers_Crowing_the_Himalaya_ExpedSummitBid.png"),
    width = 21, height = 21, units = "cm", dpi = 600
)

ggsave(
    Plots_ExpedSizes$Supp,
    filename = file.path(Dir, "SUPP_Mountaineers_Crowing_the_Himalaya_ExpedNumberAndSizes.png"),
    width = 21, height = 21, units = "cm", dpi = 600
)

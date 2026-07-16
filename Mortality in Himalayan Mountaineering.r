#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  ANALYSIS & VISUALIZATIONS:
#'      - Mortality over time (both in general and conitional on any death)
#'      - Visualisations of mortality by reason and altitude
#'      - Avalanche mortality by season and decade
#'  DEPENDENCIES:
#'  - "Himalayan Summit Success - DATA.r"
#'  - "PlottingFunctions.r"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
# source("Himalayan Summit Success - DATA.r")
source("PlottingFunctions.r")
Dir <- file.path(Dir.Exports, "Mortality in Himalayan Mountaineering")
if (!dir.exists(Dir)) dir.create(Dir, recursive = TRUE)

# BAYESIAN MODEL OF MORTALITY =============================================
## Data -------------------------------------------------------------------
df_data <- ModelData_df <- Expeditions_df
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
## Actual Model -----------------------------------------------------------
if (file.exists(file.path(Dir.Exports, "model_MT.RData"))) {
    load(file.path(Dir.Exports, "model_MT.RData"))
} else {
    model_MT <- brm(
        bf(
            Mortality_all ~ year_0 + (year_0 | PEAKID), # random intercepts & slopes per peak
            phi ~ 1, # constant precision
            zoi ~ year_0 + (year_0 | PEAKID) # model zero/one inflation similarly
        ),
        data = df_data,
        family = zero_one_inflated_beta(),
        chains = 4,
        cores = 4,
        iter = 10000,
        warmup = 5000,
        inits = "0",
        seed = 123,
        control = list(adapt_delta = 0.95, max_treedepth = 15)
    )
    save(model_MT, file = file.path(Dir.Exports, "model_MT.RData"))
}
# summary(model_MT)
# pp_check(model_MT, ndraws = 500)
conditions <- expand.grid(
    PEAKID = sort(unique(df_data$PEAKID)),
    year_0 = sort(unique(df_data$year_0))
)
rownames(conditions) <- NULL
# Conditional mortality (given >0)
mortality_draws <- model_MT %>%
    add_epred_draws(newdata = conditions, re_formula = NULL) %>%
    rename(mort_rate = .epred)
# Probability of any mortality
hurdle_draws <- model_MT %>%
    add_epred_draws(newdata = conditions, dpar = "zoi", re_formula = NULL) %>%
    rename(prob_mort = zoi)

## Preparing Plotting Data and Colours ------------------------------------
MortDeath <- list(
    df = left_join(
        mortality_draws %>% select(PEAKID, year_0, .draw, mort_rate),
        hurdle_draws %>% select(PEAKID, year_0, .draw, prob_mort),
        by = c("PEAKID", "year_0", ".draw")
    ),
    cols = c("#006D75", "#750800"),
    columns = c("mort_rate", "prob_mort"),
    names = c("Probability of Any Mortality [%]", "Conditional Mortality [%]")
)
MortDeath$df <- MortDeath$df %>%
    mutate(
        mort_rate = mort_rate * 100,
        prob_mort = prob_mort * 100
    )

## Plotting ---------------------------------------------------------------
Plots_Bayes <- FUN.BayesianPlot(
    plot_ls = MortDeath, ScaleFac = 5
)
Plots_Bayes$Main

# MORTALITY BY CAUSE ======================================================
## Data -------------------------------------------------------------------
d <- members_df %>%
    mutate(
        DEATH     = as.logical(DEATH),
        MYEAR     = as.numeric(as.character(MYEAR)),
        MSEASON   = as.numeric(as.character(MSEASON)),
        DEATHHGTM = as.numeric(as.character(DEATHHGTM)),
        DEATHTYPE = as.numeric(as.character(DEATHTYPE))
    ) %>%
    filter(
        DEATH == TRUE,
        PEAKID %in% TargetIDs,
        MYEAR >= Years_vec[1], MYEAR <= tail(Years_vec, 1),
        !is.na(DEATHHGTM)
    )
d <- d[d$DEATHHGTM != 0, ]
panel_a_data <- d %>%
    filter(DEATHTYPE %in% as.numeric(names(DEATHTYPE_LABELS))) %>%
    mutate(cause = factor(DEATHTYPE_LABELS[as.character(DEATHTYPE)], levels = CAUSE_ORDER))

n_by_cause <- panel_a_data %>% count(cause)
cause_labels_with_n <- setNames(
    sprintf("%s\n(n=%d)", n_by_cause$cause, n_by_cause$n),
    as.character(n_by_cause$cause)
)
panel_a_data$cause_label <- factor(cause_labels_with_n[as.character(panel_a_data$cause)],
    levels = cause_labels_with_n[CAUSE_ORDER]
)

## Plotting ---------------------------------------------------------------
Plots_Cause <- ggplot(panel_a_data, aes(x = cause_label, y = DEATHHGTM)) +
    geom_boxplot(
        outlier.shape = NA, fill = "grey", color = "#333333",
        linewidth = 0.6, width = 0.55
    ) +
    geom_hline(yintercept = 8000, color = "#999999", linewidth = 0.5, linetype = "dotted") +
    annotate("text",
        x = 5.4, y = 8000, label = "8000 m", angle = 90,
        vjust = -0.4, size = 2.6, color = "#777777"
    ) +
    coord_cartesian(ylim = c(3800, 8600)) +
    labs(x = NULL, y = "Death altitude (m a.s.l.)") +
    theme_bw(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0, size = 11),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8)
    )

# AVALANCHE-MORTALITY BY ELEVATION AND SEASON =============================
## Data -------------------------------------------------------------------
d <- members_df %>%
    mutate(
        MYEAR     = as.numeric(as.character(MYEAR)),
        MSEASON   = as.numeric(as.character(MSEASON)),
        DEATHHGTM = as.numeric(as.character(DEATHHGTM)),
        DEATHTYPE = as.numeric(as.character(DEATHTYPE)),
        DEATH     = as.logical(DEATH)
    ) %>%
    filter(
        DEATH == TRUE,
        PEAKID %in% TargetIDs,
        DEATHTYPE == 7, # avalanche
        MYEAR >= Years_vec[1], MYEAR <= tail(Years_vec, 1),
        MSEASON %in% c(1, 3), # 1 = pre-Monsoon, 3 = post-Monsoon
        !is.na(DEATHHGTM)
    ) %>%
    mutate(
        season_label = factor(MSEASON,
            levels = c(1, 3),
            labels = c("Pre-Monsoon", "Post-Monsoon")
        ),
        peak_label = factor(PEAKID,
            levels = TargetIDs,
            labels = eightks_sf$PKNAME[match(TargetIDs, eightks_sf$ID)]
        )
    )
d <- d[d$DEATHHGTM != 0, ]
# cat(sprintf("\nPART 1: %d avalanche deaths with valid altitude, season, and year.\n", nrow(d)))

# --- Summary table: n and median altitude per peak x season ---
summary_tbl <- d %>%
    group_by(peak_label, season_label) %>%
    summarise(n = n(), median_alt = median(DEATHHGTM), .groups = "drop") %>%
    arrange(peak_label, season_label)
# print(summary_tbl)
write.csv(summary_tbl, file.path(Dir, "deadly_avalanche_mean_altitude_summary.csv"), row.names = FALSE)

# --- Pooled (all peaks combined) comparison ---
pre_all <- d$DEATHHGTM[d$season_label == "Pre-Monsoon"]
post_all <- d$DEATHHGTM[d$season_label == "Post-Monsoon"]
pooled_test <- wilcox.test(pre_all, post_all)
sink(file.path(Dir, "deadly_avalanche_mean_altitude_wilcoxtest.txt"))
cat(sprintf(
    "\nPooled: Pre median=%.0fm (n=%d), Post median=%.0fm (n=%d), Wilcoxon p=%.4f\n",
    median(pre_all), length(pre_all), median(post_all), length(post_all), pooled_test$p.value
))

# --- Per-peak Wilcoxon test, only where both seasons have n >= 3 ---
n_by_peak_season <- d %>%
    count(peak_label, season_label) %>%
    pivot_wider(names_from = season_label, values_from = n, values_fill = 0)

testable_peaks <- n_by_peak_season %>%
    filter(`Pre-Monsoon` >= 3, `Post-Monsoon` >= 3) %>%
    pull(peak_label)

cat("\nPer-peak Wilcoxon tests (only peaks with n>=3 in both seasons):\n")
for (p in as.character(testable_peaks)) {
    sub <- d %>% filter(peak_label == p)
    test <- wilcox.test(DEATHHGTM ~ season_label, data = sub)
    diff <- median(sub$DEATHHGTM[sub$season_label == "Post-Monsoon"]) -
        median(sub$DEATHHGTM[sub$season_label == "Pre-Monsoon"])
    cat(sprintf("  %-20s diff=%+.0fm  p=%.4f\n", p, diff, test$p.value))
}
excluded_peak_labels <- setdiff(eightks_sf$PKNAME[match(TargetIDs, eightks_sf$ID)], as.character(testable_peaks))
cat(
    "Peaks excluded from testing (n<3 in at least one season): ",
    paste(excluded_peak_labels, collapse = ", "), "\n"
)
sink()

## Plotting ---------------------------------------------------------------
d_plot <- d %>% # add decade identity
    filter(MYEAR <= tail(Years_vec, 1) - 1) %>%
    mutate(decade_label = cut(
        MYEAR,
        breaks = seq(Years_vec[1] - 1, tail(Years_vec, 1) - 1 + 10, by = 10),
        include.lowest = TRUE,
        right = FALSE,
        labels = paste0(seq(Years_vec[1] - 1, tail(Years_vec, 1) - 1, by = 10), "s")
    ))

Plots_Avalanche <- ggplot(d_plot, aes(x = decade_label, y = DEATHHGTM, fill = season_label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55) +
    # geom_jitter(aes(color = season_label), width = 0.15, size = 1.3, alpha = 0.7) +
    # facet_wrap(~peak_label, nrow = 3) +
    scale_fill_manual(values = c("Pre-Monsoon" = PreColour, "Post-Monsoon" = PostColour)) +
    scale_color_manual(values = c("Pre-Monsoon" = PreColour, "Post-Monsoon" = PostColour)) +
    geom_hline(yintercept = 8000, color = "#999999", linewidth = 0.5, linetype = "dotted") +
    annotate("text",
        x = 5.4, y = 8000, label = "8000 m", angle = 90,
        vjust = -0.4, size = 2.6, color = "#777777"
    ) +
    coord_cartesian(ylim = c(3800, 8600)) +
    labs(x = NULL, y = "Death altitude (m a.s.l.)", fill = "") +
    theme_bw(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0, size = 11),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8),
        legend.position = # place legend inside plot area
            c(0.185, 0.9),
        legend.background = element_rect(fill = "transparent") # Transparent legend box
    )

# PLOT SAVING =============================================================
## Main Text --------------------------------------------------------------
Main_gg <- plot_grid(
    Plots_Bayes$Main,
    plot_grid(
        Plots_Cause,
        Plots_Avalanche,
        nrow = 1,
        labels = c("B", "C")
    ),
    ncol = 1, labels = c("A", ""), rel_heights = c(1, 1)
)

ggsave(
    Main_gg,
    filename = file.path(Dir, "Figure_Mortality_in_Himalayan_Mountaineering.png"),
    width = 36, height = 28, units = "cm", dpi = 600
)

## Supplement -------------------------------------------------------------
ggsave(
    Plots_Bayes$Supp,
    filename = file.path(Dir, "SUPP_Mortality_in_Himalayan_Mountaineering.png"),
    width = 21, height = 42, units = "cm", dpi = 600
)

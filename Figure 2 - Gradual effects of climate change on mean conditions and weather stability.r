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

message("Figure 2 Plotting")
FName <- "Figure2"

# PREPARATION =============================================================
peaks_ts_df <- peaks_ts_df %>%
    filter(Variable %in% MainVars$VarName) %>%
    filter(PeakID %in% eightks_sf$PKNAME) %>%
    mutate(YEAR = as.numeric(substr(Date, 1, 4))) %>%
    filter(Variable %in% MainVars$VarName) %>%
    group_by(PeakID, YEAR, Season, Variable) %>%
    summarise(
        Mean = mean(Value, na.rm = TRUE),
        SD = sd(Value, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    left_join(MainVars, by = c("Variable" = "VarName")) %>%
    mutate(Season = factor(Season,
        levels = c("Pre", "Post"),
        labels = c("Pre-Monsoon", "Post-Monsoon")
    )) %>%
    mutate(Variable = ifelse(is.na(ClearName), Variable, ClearName)) %>%
    select(-ClearName)

pred_df <- pred_df %>%
    filter(Variable %in% MainVars$VarName) %>%
    filter(PeakID %in% eightks_sf$PKNAME) %>%
    left_join(MainVars, by = c("Variable" = "VarName")) %>%
    mutate(Season = factor(Season,
        levels = c("Pre", "Post"),
        labels = c("Pre-Monsoon", "Post-Monsoon")
    ))

# PANEL A: MEAN CHANGE ====================================================
mainA_data <- peaks_ts_df %>%
    group_by(Variable, Season, YEAR) %>%
    summarise(Mean = mean(Mean, na.rm = TRUE), .groups = "drop")

mainB_data <- pred_df %>%
    group_by(ClearName, Season, Year) %>%
    summarise(AR1 = mean(AR1, na.rm = TRUE), .groups = "drop")
colnames(mainB_data)[1] <- "Variable"

mainAB_ls <- lapply(unique(mainA_data$Variable), FUN = function(v) {
    # v <- unique(mainA_data$Variable)[1]
    mean_iter <- mainA_data %>% filter(Variable == v)
    ar1_iter <- mainB_data %>% filter(Variable == v)

    mean_min <- min(mean_iter$Mean, na.rm = TRUE)
    mean_max <- max(mean_iter$Mean, na.rm = TRUE)
    ar1_min <- min(ar1_iter$AR1, na.rm = TRUE)
    ar1_max <- max(ar1_iter$AR1, na.rm = TRUE)

    if (isTRUE(all.equal(ar1_max, ar1_min))) {
        scale_factor <- 1
        offset <- mean_min - ar1_min
    } else {
        scale_factor <- (mean_max - mean_min) / (ar1_max - ar1_min)
        offset <- mean_min - ar1_min * scale_factor
    }

    ggplot() +
        stat_smooth(
            data = mean_iter,
            aes(x = YEAR, y = Mean, color = Season, fill = Season, group = Season, linetype = "Climate Trend"),
            method = "lm",
            se = TRUE
        ) +
        stat_smooth(
            data = ar1_iter,
            aes(x = Year, y = AR1 * scale_factor + offset, color = Season, fill = Season, group = Season, linetype = "Autocorrelation (AR1)"),
            method = "lm",
            se = TRUE
        ) +
        scale_color_manual(values = c("Pre-Monsoon" = "#6b74eb", "Post-Monsoon" = "#b45303")) +
        scale_fill_manual(values = c("Pre-Monsoon" = "#6b74eb", "Post-Monsoon" = "#b45303")) +
        scale_linetype_manual(
            values = c("Climate Trend" = "solid", "Autocorrelation (AR1)" = "dashed"),
            name = "Trendline"
        ) +
        scale_y_continuous(
            name = v,
            sec.axis = sec_axis(~ (. - offset) / scale_factor, name = "Day-to-Day Stability/Predictability (AR1)")
        ) +
        labs(x = "Year") +
        theme_bw() +
        guides(
            color = guide_legend(
                order = 1,
                keywidth = grid::unit(1.8, "cm"),
                override.aes = list(linewidth = 1.1)
            ),
            fill = "none",
            linetype = guide_legend(
                order = 2,
                keywidth = grid::unit(1.8, "cm"),
                override.aes = list(color = "black", linewidth = 1.1)
            )
        ) +
        theme(
            legend.position = "top",
            legend.key.width = grid::unit(1.8, "cm")
        )
})

get_legend <- function(myggplot) {
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

mainAB <- cowplot::plot_grid(
    get_legend(mainAB_ls[[1]]),
    cowplot::plot_grid(
    plotlist = lapply(mainAB_ls, FUN = function(x){x+theme(legend.position = "none")}), 
    ncol = 1, align = "v", labels = "AUTO"),
    ncol = 1, rel_heights = c(0.025, 1)
)

ggsave(
    mainAB,
    filename = file.path(Dir.Exports, paste0(FName, "_main.png")),
    width = 24, height = 32, units = "cm", dpi = 300
)


# SUPPLEMENT ==============================================================
## Tables -----------------------------------------------------------------
stop("Save numbers associated with the figure to CSV files for supplementary material")

## Individual peak plots --------------------------------------------------
suppA <- ggplot(mean_df, aes(x = YEAR, y = Mean, color = PeakID, group = PeakID)) +
    stat_smooth(method = "lm", se = TRUE) +
    facet_wrap(~ Season + Variable, scales = "free") +
    scale_color_viridis_d() +
    labs(
        x = "Year", y = "Mean Value"
    ) +
    theme_bw() +
    theme(legend.position = "top")

stop("also produce supp for pred_df")
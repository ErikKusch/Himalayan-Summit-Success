#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  ANALYSIS & VISUALIZATIONS: Statistical analysis and visualization of
#'      - Weather predictability
#'      - Extreme weather events
#'      - Avalanches and storms
#'  DEPENDENCIES:
#'  - "Himalayan Summit Success - DATA.R"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
# source("Himalayan Summit Success - DATA.r")
Dir <- file.path(Dir.Exports, "Objective hazards in Himalayan mountaineering exacerbated by climate change")
if (!dir.exists(Dir)) dir.create(Dir, recursive = TRUE)

# AVALANCHES AND STORMS ===================================================
## Data -------------------------------------------------------------------
Expeditions_df <- Expeditions_df %>%
    # Ensure dates are Date objects
    mutate(
        SMTDATE = as.Date(SMTDATE),
        BCDATE = as.Date(BCDATE),
        SMT_month = month(SMTDATE),
        BC_month = month(BCDATE)
    ) %>%
    # Define SEASON
    mutate(
        SEASON = case_when(
        SMT_month %in% 3:5 & BC_month %in% 3:5 ~ "Pre-Monsoon",
        SMT_month %in% 9:11 & BC_month %in% 9:11 ~ "Post-Monsoon",
        TRUE ~ NA_character_ # remove all other data
    )
    )

seasonal_df <- Expeditions_df %>%
    group_by(YEAR, SEASON) %>%
    summarise(
        storms = sum(STORM, na.rm = TRUE),
        avalanches = sum(AVALANCHE, na.rm = TRUE),
        deaths = sum(MDEATHS, na.rm = TRUE)
    ) %>%
    arrange(SEASON, YEAR) %>%
    group_by(SEASON) %>%
    mutate(
        cum_storms = cumsum(storms),
        cum_avalanches = cumsum(avalanches)
    ) %>%
    ungroup()

# Compute cumulative ratio per season
seasonal_df <- seasonal_df %>%
    mutate(ratio_avalanche_storm = ifelse(storms > 0, avalanches / storms, NA))
seasonal_df <- seasonal_df %>%
    group_by(SEASON) %>%
    mutate(cum_ratio = ifelse(cum_storms > 0, cum_avalanches / cum_storms, NA)) %>%
    ungroup()

# label placement
mid_x <- seasonal_df %>%
    group_by(SEASON) %>%
    summarise(mid_x = mean(cum_storms, na.rm = TRUE))
seasonal_df <- seasonal_df %>%
    left_join(mid_x, by = "SEASON") %>%
    mutate(
        label_hjust = ifelse(cum_storms > mid_x, 1, 0),
        label_nudge_x = ifelse(cum_storms > mid_x, -0.5, 0.5)
    )
seasonal_df <- seasonal_df[!is.na(seasonal_df$SEASON), ]

# secondary axis scaling
max_cum_avalanches <- max(seasonal_df$cum_avalanches, na.rm = TRUE)
max_cum_ratio <- max(seasonal_df$cum_ratio, na.rm = TRUE)
max_cum_storms <- max(seasonal_df$cum_storms, na.rm = TRUE)
scale_factor <- max_cum_avalanches / max_cum_ratio
time_scale_factor <- max_cum_avalanches / max_cum_storms

## Plotting ---------------------------------------------------------------
# (1) ratio of avalanches to storms over time by season
# cumulative avalanches vs storms with cumulative ratio, combined seasons
AvaStormCombined_gg <- ggplot(
    seasonal_df,
    aes(
        x = cum_storms,
        y = cum_avalanches,
        color = factor(SEASON, levels = c("Pre-Monsoon", "Post-Monsoon"))
        )
    ) +
    geom_line(aes(group = SEASON), size = 1.25) +
    geom_point(aes(size = deaths, shape = deaths == 0), alpha = 0.6, na.rm = TRUE) +
    # stat_smooth(
    #     aes(group = SEASON),
    #     method = "lm",
    #     se = TRUE,
    #     size = 0.8,
    #     na.rm = TRUE
    # ) +
    geom_text_repel(
        data = subset(seasonal_df, !is.na(deaths) & deaths > 0),
        aes(label = YEAR, hjust = label_hjust),
        nudge_x = subset(seasonal_df, !is.na(deaths) & deaths > 0)$label_nudge_x * 2,
        direction = "both",
        segment.color = "grey50",
        size = 3.8,
        box.padding = 0.8,
        point.padding = 0.5,
        show.legend = FALSE
    ) +
#   geom_line(
#     aes(y = cum_ratio * scale_factor, group = SEASON),
#     linetype = "dashed",
#     size = 1
#   ) +
    scale_color_manual(
        values = c("Pre-Monsoon" = PreColour, "Post-Monsoon" = PostColour),
        name = "Season"
    ) +
    # scale_y_continuous(
    #     name = "Cumulative Avalanche Reports",
    #     sec.axis = sec_axis(
    #     trans = ~ . / scale_factor,
    #     name = "Cumulative Avalanche / Cumulative Storm Ratio"
    #     )
    # ) +
    scale_size_continuous(
        name = "Deaths per Season",
        range = c(3, 9),
        breaks = function(x) {
            br <- pretty(x)
            br[br > 0]
        }
    ) +
    scale_shape_manual(
        values = c("TRUE" = 5, "FALSE" = 18),
        labels = c("TRUE" = "0 deaths", "FALSE" = ">0 deaths"),
        name = "Death Count"
    ) +
    labs(
        y = "Cumulative Avalanche Reports",
        x = "Cumulative Storm Reports"
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 12),
        axis.title.y.left = element_text(size = 12),
        axis.title.y.right = element_text(size = 12),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.justification = "center",
        legend.spacing.x = grid::unit(0.8, "cm"),
        legend.key.height = grid::unit(0.5, "cm"),
        legend.text = element_text(vjust = 0.5),
        legend.title = element_text(vjust = 0.5)
    ) +
    guides(
        color = guide_legend(order = 1, nrow = 1, byrow = TRUE),
        size = guide_legend(order = 2, nrow = 1, byrow = TRUE),
        shape = guide_legend(order = 3, nrow = 1, byrow = TRUE)
    )
AvaStormCombined_gg

# (2) cumulative storms over time by season
AvaTimeStorms_gg <- ggplot(
    seasonal_df,
    aes(
        x = YEAR,
        y = cum_storms,
        color = factor(SEASON, levels = c("Pre-Monsoon", "Post-Monsoon")),
        group = SEASON
    )
    ) +
    geom_line(size = 1.15) +
    geom_point(pch = 19, size = 4, alpha = 1, na.rm = TRUE) +
    scale_color_manual(
        values = c("Pre-Monsoon" = PreColour, "Post-Monsoon" = PostColour),
        name = "Season"
    ) +
    labs(
        x = "Year",
        y = "Cumulative Storm Reports"
    ) +
    theme_bw() +
    theme(
        legend.position = "top",
        legend.direction = "horizontal",
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold")
    ) +
    lims(x = c(1975, NA))
    AvaTimeStorms_gg

# (3) cumulative avalanches over time by season
AvaTimeAvalanches_gg <- ggplot(
    seasonal_df,
    aes(
        x = YEAR,
        y = cum_avalanches,
        color = factor(SEASON, levels = c("Pre-Monsoon", "Post-Monsoon")),
        group = SEASON
    )
    ) +
    geom_line(size = 1.15) +
    geom_point(pch = 18, size = 4, alpha = 1, na.rm = TRUE) +
    scale_color_manual(
        values = c("Pre-Monsoon" = PreColour, "Post-Monsoon" = PostColour),
        name = "Season"
    ) +
    labs(
        x = "Year",
        y = "Cumulative Avalanche Reports"
    ) +
    theme_bw() +
    theme(
        legend.position = "top",
        legend.direction = "horizontal",
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold")
    ) +
    lims(x = c(1975, NA))
AvaTimeAvalanches_gg

# convenience object containing both plots
plot_grid(
    AvaStormCombined_gg,
    plot_grid(
        AvaTimeStorms_gg+theme(legend.position = "none"), 
        AvaTimeAvalanches_gg+theme(legend.position = "none"),
        nrow = 1, align = "v", labels = c("(B) Cumulative Storms", "(C) Cumulative Avalanches"), label_size = 12),
    ncol = 1,
    align = "v",
    labels = c("(A) Cumulative Avalanche vs Storm Reports"),
    label_size = 12
    )





















# EXTREMES ================================================================
## Data -------------------------------------------------------------------


## Plotting ---------------------------------------------------------------








































































# PLOTS ===================================================================

## Extremes ---------------------------------------------------------------
# Step 1: Filter for main variables and relevant extremes
df_filtered <- peaks_ts_df %>%
  filter(
    Variable %in% MainVars$VarName, # limit to main variables
    PeakID %in% eightks_sf$PKNAME,
    Extreme %in% c("HIGH", "LOW")
  ) %>%
  mutate(
    Year = year(Date)
  ) %>%
  arrange(Variable, PeakID, Year, Season, Date)

# Step 2: Compute counts per Variable × PeakID × Year × Season × Extreme
extreme_counts <- df_filtered %>%
  group_by(Variable, PeakID, Year, Season, Extreme) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  )

# Step 3: Compute mean extreme ratio per group
df2 <- peaks_ts_df %>%
  filter(
    Variable %in% MainVars$VarName, # limit to main variables
    PeakID %in% eightks_sf$PKNAME
  ) %>%
  mutate(
    Year = year(Date)
  ) %>%
  arrange(Variable, PeakID, Year, Season, Date)
mean_ratio_df <- df2 %>%
  group_by(Variable, PeakID, Year, Season) %>%
  summarise(
    MeanExtremeRatioHigh = mean(ExtremeRatioHigh, na.rm = TRUE),
    MeanExtremeRatioLow = mean(ExtremeRatioLow, na.rm = TRUE),
    .groups = "drop"
  )
mean_ratio_long <- mean_ratio_df %>%
  pivot_longer(
    cols = c(MeanExtremeRatioHigh, MeanExtremeRatioLow),
    names_to = "Extreme",
    values_to = "MeanExtremeRatio"
  ) %>%
  mutate(
    Extreme = case_when(
      Extreme == "MeanExtremeRatioHigh" ~ "HIGH",
      Extreme == "MeanExtremeRatioLow" ~ "LOW"
    )
  ) %>%
  arrange(Variable, PeakID, Season, Extreme, Year)

# Step 4: Compute average run length per group
avg_run_df <- df_filtered %>%
  group_by(Variable, PeakID, Year, Season, Extreme) %>%
  arrange(Date) %>% # ensure chronological order
  summarise(
    AvgRun = {
      # convert Date to integer days (numeric)
      dates_int <- as.integer(as.Date(Date))
      # identify consecutive runs
      run_id <- cumsum(c(TRUE, diff(dates_int) != 1))
      # lengths of each run
      run_lengths <- table(run_id)
      # average run length
      mean(as.numeric(run_lengths))
    },
    .groups = "drop"
  )

# Step 5: Merge all summaries into one long-format summary
all_combinations <- expand_grid(
  Variable = MainVars$VarName,
  PeakID = unique(df_filtered$PeakID),
  Year = unique(df_filtered$Year),
  Season = c("Pre", "Post"),
  Extreme = c("HIGH", "LOW")
)
summary_long_full <- all_combinations %>%
  left_join(extreme_counts, by = c("Variable", "PeakID", "Year", "Season", "Extreme")) %>%
  left_join(mean_ratio_long, by = c("Variable", "PeakID", "Year", "Season", "Extreme")) %>%
  left_join(avg_run_df, by = c("Variable", "PeakID", "Year", "Season", "Extreme")) %>%
  # 3. Replace NAs for counts and run length with 0
  mutate(
    Count = replace_na(Count, 0),
    AvgRun = replace_na(AvgRun, 0)
  ) %>%
  arrange(Variable, PeakID, Season, Extreme, Year) %>%
  group_by(Variable, PeakID, Season, Extreme) %>%
  mutate(
    Cumulative_Count = cumsum(Count),
    Cumulative_MeanExtremeRatio = cumsum(MeanExtremeRatio),
    Cumulative_RunLength = cumsum(AvgRun)
  ) %>%
  ungroup() %>%
  mutate(
    Season = case_when(
      Season == "Pre" ~ "Pre-Monsoon",
      Season == "Post" ~ "Post-Monsoon",
      TRUE ~ Season
    ),
    Season = factor(Season, levels = c("Pre-Monsoon", "Post-Monsoon"))
  )

## non-cumulative lines
#### this is supplement!!!!
extremes_combined <- lapply(MainVars$VarName, FUN = function(var) {
  # var = "2m_temperature"
  # Filter data for the current variable
  plot_data <- summary_long_full %>%
    filter(Variable == var, Extreme %in% c("LOW", "HIGH"))

  # Compute scaling factor for HIGH relative to LOW
  low_max <- max(plot_data$MeanExtremeRatio[plot_data$Extreme == "LOW"], na.rm = TRUE)
  high_max <- max(plot_data$MeanExtremeRatio[plot_data$Extreme == "HIGH"], na.rm = TRUE)
  scale_factor <- low_max / high_max

  # Create the plot
  p <- ggplot(plot_data, aes(x = Year, y = MeanExtremeRatio, color = Extreme, shape = Extreme)) +
    geom_point(size = 2) +
    # geom_line(aes(group = Extreme)) +
    # Separate smoothing lines for LOW and HIGH
    geom_smooth(
      data = plot_data %>% filter(Extreme == "LOW"),
      aes(x = Year, y = MeanExtremeRatio),
      method = "lm", se = FALSE, color = "#003575"
    ) +
    geom_smooth(
      data = plot_data %>% filter(Extreme == "HIGH"),
      aes(x = Year, y = MeanExtremeRatio),
      method = "lm", se = FALSE, color = "#754000"
    ) +
    facet_grid(PeakID ~ Season) +
    scale_color_manual(values = c("LOW" = "#003575", "HIGH" = "#754000")) +
    scale_shape_manual(values = c("LOW" = 25, "HIGH" = 24)) +
    labs(x = "Year", y = "Ratio of Daily Values to Extreme Threshold (95% Quantile)") +
    theme_bw() +
    theme(legend.position = "top") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50")

  return(p)
})
names(extremes_combined) <- MainVars$VarName

## cumulative lines
### this is main text!!!!
### Compute linear model stats per Season ---------------------------
# Nest by Season
extremes_cumul <- lapply(MainVars$VarName, FUN = function(var) {
  extremes_ls <- lapply(c("LOW", "HIGH"), FUN = function(extreme) {
    plot_data <- summary_long_full %>%
      filter(Extreme == extreme, Variable == var)
    nested_lm <- plot_data %>%
      group_by(Season) %>%
      nest() %>%
      mutate(
        lm_count = map(data, ~ lm(Cumulative_Count ~ Year, data = .x)),
        lm_run   = map(data, ~ lm(Cumulative_RunLength ~ Year, data = .x))
      )

    # Extract slope and R² for Cumulative_Count safely
    lm_stats_count <- nested_lm %>%
      mutate(
        tidied = map(lm_count, tidy),
        glanced = map(lm_count, glance)
      ) %>%
      unnest(c(tidied, glanced), names_sep = "_") %>%
      filter(tidied_term == "Year") %>%
      select(Season, estimate = tidied_estimate, r.squared = glanced_r.squared)

    # Similarly for Cumulative_RunLength
    lm_stats_run <- nested_lm %>%
      mutate(
        tidied = map(lm_run, tidy),
        glanced = map(lm_run, glance)
      ) %>%
      unnest(c(tidied, glanced), names_sep = "_") %>%
      filter(tidied_term == "Year") %>%
      select(Season, estimate = tidied_estimate, r.squared = glanced_r.squared)


    # Determine annotation positions
    # X axis range
    x_min <- min(plot_data$Year, na.rm = TRUE)
    x_max <- max(plot_data$Year, na.rm = TRUE)
    x_range <- x_max - x_min

    # Y axis ranges
    y_count_min <- min(plot_data$Cumulative_Count, na.rm = TRUE)
    y_count_max <- max(plot_data$Cumulative_Count, na.rm = TRUE)
    y_count_range <- y_count_max - y_count_min

    y_run_min <- min(plot_data$Cumulative_RunLength, na.rm = TRUE)
    y_run_max <- max(plot_data$Cumulative_RunLength, na.rm = TRUE)
    y_run_range <- y_run_max - y_run_min

    # Pre- and Post-Monsoon positions as percentages of range
    annotation_positions_count <- tibble(
      Season = rev(c("Pre-Monsoon", "Post-Monsoon")),
      x = c(x_min + 0.8 * x_range, x_min + 0.2 * x_range),
      y = c(y_count_min + 0.2 * y_count_range, y_count_min + 0.8 * y_count_range)
    )

    annotation_positions_run <- tibble(
      Season = rev(c("Pre-Monsoon", "Post-Monsoon")),
      x = c(x_min + 0.8 * x_range, x_min + 0.2 * x_range),
      y = c(y_run_min + 0.2 * y_run_range, y_run_min + 0.8 * y_run_range)
    )

    # Merge annotation text
    annotation_count <- annotation_positions_count %>%
      left_join(lm_stats_count, by = "Season") %>%
      mutate(
        label = paste0(Season, "\nSlope = ", round(estimate, 3), "\nR² = ", round(r.squared, 3))
      )

    annotation_run <- annotation_positions_run %>%
      left_join(lm_stats_run, by = "Season") %>%
      mutate(
        label = paste0(Season, "\nSlope = ", round(estimate, 3), "\nR² = ", round(r.squared, 3))
      )

    # Plot 1: Cumulative_Count
    count_gg <- ggplot(plot_data, aes(x = Year, y = Cumulative_Count, group = PeakID, col = Season)) +
      geom_point(alpha = 0.5, size = 1.3) +
      stat_smooth(aes(group = Season, col = Season), method = "lm") +
      geom_text(
        data = annotation_count, # precomputed slope/R² labels
        aes(x = x, y = y, label = label, col = Season),
        inherit.aes = FALSE,
        size = 4
      ) +
      scale_color_manual(values = c("Pre-Monsoon" = "#003575", "Post-Monsoon" = "#754000")) +
      theme_bw() +
      labs(x = "Year", y = "Number of Extreme Events [#]") +
      theme(
        strip.text = element_text(face = "bold"),
        axis.title = element_text(size = 12),
        legend.position = "none"
      )

    # Second plot: Cumulative_RunLength
    run_gg <- ggplot(plot_data, aes(x = Year, y = Cumulative_RunLength, group = PeakID, col = Season)) +
      geom_point(alpha = 0.5, size = 1.3) +
      stat_smooth(aes(group = Season, col = Season), method = "lm") +
      geom_text(
        data = annotation_run, # precomputed slope/R² labels
        aes(x = x, y = y, label = label, col = Season),
        inherit.aes = FALSE,
        size = 4
      ) +
      scale_color_manual(values = c("Pre-Monsoon" = "#6F0075", "Post-Monsoon" = "#067500")) +
      theme_bw() +
      labs(x = "Year", y = "Cumulative Run Length (days)") +
      theme(
        strip.text = element_text(face = "bold"),
        axis.title = element_text(size = 12),
        legend.position = "none"
      )

    plot_grid(count_gg, run_gg, ncol = 2)
  })
  names(extremes_ls) <- c("LOW", "HIGH")
  extremes_ls
})
names(extremes_cumul) <- MainVars$VarName

## Saving Plots -----------------------------------------------------------
label_row <- function(text) {
  ggplot() +
    geom_rect(
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      fill = "white",
      color = NA
    ) +
    annotate(
      "text",
      x = 0, y = 0,
      label = text,
      hjust = 0,
      fontface = "bold",
      size = 5
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1), clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
}
### Main Text +++++++
label_A <- label_row("(A) Compound Extremes: Avalanches & Storms")
label_B <- label_row("(B) Day-to-Day Predictability")
label_C <- label_row("(C) Low-Temperature Extremes")
label_D <- label_row("(D) High-Temperature Extremes")

main_ggs <- plot_grid(
  label_A,
  AvaStorm_gg,
  label_B,
  AR_gg[[1]][[2]]$Density$All,
  label_C,
  extremes_cumul[["2m_temperature"]]$LOW,
  label_D,
  extremes_cumul[["2m_temperature"]]$HIGH,
  ncol = 1,
  rel_heights = c(
    0.08, 1, # A label + plot
    0.08, 1, # B label + plot
    0.08, 1, # C label + plot
    0.08, 1 # D label + plot
  )
)
ggsave(
  main_ggs,
  file = FName,
  width = 14, height = 22
)

### Supplement +++++++
#### AR Line trends per Variable
lapply(names(AR_gg), FUN = function(x) {
  # x <- "2m_temperature"
  plotlist <- lapply(AR_gg[[x]], "[[", "Lines")
  labellist <- lapply(names(plotlist), FUN = function(y) {
    label_row(paste0("(", LETTERS[which(names(AR_gg) == x) + 1], ") Day-to-Day Predictability: ", y))
  })
  AR_ggsupp <- plot_grid(
    labellist[[1]],
    plotlist[[1]],
    labellist[[2]],
    plotlist[[2]],
    labellist[[3]],
    plotlist[[3]],
    labellist[[4]],
    plotlist[[4]],
    labellist[[5]],
    plotlist[[5]],
    ncol = 1,
    rel_heights = c(
      0.08, 1, # A label + plot
      0.08, 1, # B label + plot
      0.08, 1, # C label + plot
      0.08, 1, # D label + plot
      0.08, 1 # E label + plot
    )
  )
  ggsave(
    AR_ggsupp,
    file = paste0(tools::file_path_sans_ext(FName), "_Supplement_AR_", x, ".png"),
    width = 20, height = 22 / 3 * 5
  )
})

#### Extreme Ratio Plots
lapply(names(extremes_combined), FUN = function(x) {
  # x <- "2m_temperature"
  ggsave(
    extremes_combined[[x]],
    file = paste0(tools::file_path_sans_ext(FName), "_Supplement_ExtremeRatios_", x, ".png"),
    width = 12, height = 16
  )
})

#### Cumulative Extreme Plots
lapply(names(extremes_cumul), FUN = function(x) {
  # x <- "windspeed"
  MainVars$ClearName[MainVars$VarName == x]
  label_A <- label_row(paste0("(A) Low-", MainVars$ClearName[MainVars$VarName == x], " Extremes"))
  label_B <- label_row(paste0("(B) High-", MainVars$ClearName[MainVars$VarName == x], " Extremes"))
  Cumuls_gg <- plot_grid(
    label_A,
    extremes_cumul[[x]][["LOW"]],
    label_B,
    extremes_cumul[[x]][["HIGH"]],
    ncol = 1,
    rel_heights = c(
      0.08, 1, # A label + plot
      0.08, 1 # B label + plot
    )
  )
  ggsave(
    Cumuls_gg,
    file = paste0(tools::file_path_sans_ext(FName), "_Supplement_CumulativeExtremes", x, ".png"),
    width = 14, height = 14
  )
})

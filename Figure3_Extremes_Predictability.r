#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  ANALYSIS & VISUALIZATIONS: Statistical analysis and visualization of
#'      - Weather predictability
#'      - Extreme weather events
#'  DEPENDENCIES:
#'  - "Himalayan Summit Success - DATA.R"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
source("Himalayan Summit Success - DATA.r")

message("Figure 3 Plotting")
FName <- file.path(Dir.Exports, "Figure3_.png")

# PLOTS ===================================================================
## Reported Avalanches and Storms -----------------------------------------
Expeditions_df <- read.csv(file.path(Dir.Data, "CleanedExpeditions.csv"))
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
  ) %>%
  filter(!is.na(SEASON)) # keep only pre/post-monsoon expeditions

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

# secondary axis scaling
max_cum_avalanches <- max(seasonal_df$cum_avalanches, na.rm = TRUE)
max_cum_ratio <- max(seasonal_df$cum_ratio, na.rm = TRUE)
scale_factor <- max_cum_avalanches / max_cum_ratio

# Define colors
primary_color <- "#006D75" # cumulative avalanche X storm line
secondary_color <- "#750800" # cumulative ratio line

AvaStorm_gg <- ggplot(seasonal_df, aes(x = cum_storms, y = cum_avalanches)) +
  # cumulative avalanche line
  geom_line(color = primary_color, size = 1) +
  # death points
  geom_point(aes(size = deaths), colour = primary_color, alpha = 0.6, na.rm = TRUE) +
  # linear smoother with multiple confidence intervals
  stat_smooth(method = "lm", colour = "#832c86", fill = "#832c86", alpha = 0.2, level = 0.95, na.rm = TRUE, size = 0.5) +
  stat_smooth(method = "lm", colour = "#832c86", fill = "#832c86", alpha = 0.3, level = 0.8, na.rm = TRUE, size = 0.5) +
  stat_smooth(method = "lm", colour = "#832c86", fill = "#832c86", alpha = 0.4, level = 0.5, na.rm = TRUE, size = 0.5) +
  # year labels
  geom_text_repel(
    data = subset(seasonal_df, !is.na(deaths) & deaths > 0),
    aes(label = YEAR, x = cum_storms, y = cum_avalanches, hjust = label_hjust),
    nudge_x = subset(seasonal_df, !is.na(deaths) & deaths > 0)$label_nudge_x * 2,
    direction = "both",
    segment.color = "grey50",
    size = 4,
    box.padding = 0.8,
    point.padding = 0.5
  ) +
  # cumulative ratio line, scaled
  geom_line(aes(y = cum_ratio * scale_factor), color = secondary_color, linetype = "dashed", size = 1) +
  # primary and secondary y-axis
  scale_y_continuous(
    name = "Cumulative Avalanche Reports",
    sec.axis = sec_axis(
      trans = ~ . / scale_factor,
      name = "Cumulative Avalanche / Cumulative Storm Ratio"
    )
  ) +
  scale_size_continuous(name = "Deaths per Season") +
  labs(
    x = "Cumulative Storm Reports"
  ) +
  theme_bw() +
  facet_wrap(~ factor(SEASON, levels = c("Pre-Monsoon", "Post-Monsoon")), scales = "free", ncol = 2) +
  # color the axis titles
  theme(
    axis.title.y.left = element_text(color = primary_color, size = 12),
    axis.title.y.right = element_text(color = secondary_color, size = 12),
    axis.title.x = element_text(size = 12),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  guides(size = guide_legend(nrow = 1, byrow = TRUE))

## Change in Autocorrelative Coefficient ----------------------------------
ar_cols <- grep("^AR", names(pred_df), value = TRUE)

AR_gg <- lapply(MainVars$VarName, FUN = function(i) {
  AR_gg <- lapply(ar_cols, function(col) {
    # col <- "AR10"

    ### Filter and prepare data -------------------------------------
    plot_df <- pred_df %>%
      filter(
        Variable == i,
        PeakID %in% eightks_sf$PKNAME
      ) %>%
      mutate(
        AR = .data[[col]],
        Year = as.numeric(Year),
        Season = recode(
          Season,
          "Pre"  = "Pre-Monsoon",
          "Post" = "Post-Monsoon"
        ),
        Season = factor(
          Season,
          levels = c("Pre-Monsoon", "Post-Monsoon") # enforce facet order
        )
      ) %>%
      left_join(
        summits_df %>% select(PKNAME, ID),
        by = c("PeakID" = "PKNAME")
      )

    ### Label data (per PeakID × Season trends) --------------------
    labelInfo <- lapply(
      split(plot_df, list(plot_df$ID, plot_df$Season)),
      function(dat) {
        if (nrow(dat) < 2) {
          return(NULL)
        }

        lmmod <- lm(AR ~ Year, data = dat)

        # Choose left or right label position based on renamed season
        s <- unique(dat$Season)
        label_year <- max(dat$Year) # if (s == "Pre-Monsoon") min(dat$Year) else max(dat$Year)

        LabelVal <- predict(lmmod, newdata = data.frame(Year = label_year))

        data.frame(
          ID        = unique(dat$ID),
          Season    = s,
          LabelPos  = LabelVal,
          LabelYear = label_year,
          Pval      = summary(lmmod)$coefficients["Year", "Pr(>|t|)"],
          Effect    = summary(lmmod)$coefficients["Year", "Estimate"]
        )
      }
    )
    labelInfo <- do.call(rbind, labelInfo)
    labelInfo$Fill <- as.character(labelInfo$Pval < 0.05)

    ### Overall seasonal trends -------------------------------------
    season_trends <- lapply(
      unique(plot_df$Season),
      function(season) {
        dat <- plot_df %>% filter(Season == season)
        if (nrow(dat) < 2) {
          return(NULL)
        }

        lmmod <- lm(AR ~ Year, data = dat)
        label_year <- max(dat$Year) # if (season == "Pre-Monsoon") min(dat$Year) else max(dat$Year)
        LabelVal <- predict(lmmod, newdata = data.frame(Year = label_year))

        data.frame(
          ID        = "Region",
          Season    = season,
          LabelPos  = LabelVal,
          LabelYear = label_year,
          Pval      = summary(lmmod)$coefficients["Year", "Pr(>|t|)"],
          Effect    = summary(lmmod)$coefficients["Year", "Estimate"],
          Fill      = "Region"
        )
      }
    )
    season_trends <- do.call(rbind, season_trends)
    labelInfo <- rbind(labelInfo, season_trends)

    ### Merge metadata ------------------------------------------------------
    plot_df <- plot_df %>%
      left_join(labelInfo, by = c("ID", "Season")) %>%
      mutate(StatSig = Pval < 0.05)

    ### Final plot (1-row facet, no bold) ----------------------------------
    AR_gg <- ggplot(plot_df, aes(x = Year, y = AR, group = interaction(ID, Season))) +
      geom_smooth(
        aes(linetype = StatSig, color = Effect),
        method = "lm",
        alpha = 0.3,
        linewidth = 0.8,
        show.legend = TRUE
      ) +
      geom_smooth(
        aes(group = Season),
        method = "lm",
        linewidth = 2,
        color = "#1a1b1b",
        show.legend = FALSE
      ) +
      geom_label_repel(
        data = labelInfo,
        aes(x = LabelYear, y = LabelPos, label = ID, fill = Effect),
        color = "#000000",
        size = 5,
        direction = "y",
        hjust = -0.4, # ifelse(labelInfo$Season == "Pre-Monsoon", 1.4, -0.4),
        nudge_x = 0,
        show.legend = FALSE
      ) +
      facet_wrap(~Season, ncol = 2, scales = "free_y") +
      scale_linetype_manual(
        name = "Statistical Significance of Trend",
        values = c(`FALSE` = "dashed", `TRUE` = "solid"),
        guide = guide_legend(
          title.position = "top",
          label.position = "bottom",
          nrow = 1,
          override.aes = list(size = 1, color = "black"),
          keywidth = 4,
          keyheight = 1,
          order = 2
        )
      ) +
      scale_color_viridis_c(
        direction = -1,
        name = "Mean Change per Year",
        guide = guide_colorbar(
          title.position = "top",
          barwidth = 15,
          barheight = 1.5,
          order = 1
        )
      ) +
      scale_fill_viridis_c(direction = -1) +
      labs(
        x = "Year",
        y = col
      ) +
      theme_bw() +
      lims(x = c(min(plot_df$Year) - 5, max(plot_df$Year) + 5)) +
      theme(
        strip.text = element_text(size = 12, face = "plain"), # no boldface
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.box.spacing = unit(5, "pt")
      ) +
      guides(fill = "none")

    ### Create decade bins ---------------------------------------------------
    breaks <- seq(1051, 2021, by = 10)
    plot_df$YearBin <- cut(
      plot_df$Year,
      breaks = breaks,
      include.lowest = TRUE,
      right = FALSE,
      labels = paste0(breaks[-length(breaks)], ":", breaks[-1])
    )

    ### Density Plot – Individual Peaks × Season ----------------------------
    Density_indivs <- ggplot(plot_df, aes(
      x = AR,
      color = YearBin,
      fill = YearBin
    )) +
      geom_density(alpha = 0.025, adjust = 1) +
      facet_grid(ID ~ Season, scales = "free_x") +
      scale_color_viridis_d(option = "C", direction = -1) +
      scale_fill_viridis_d(option = "C", direction = -1) +
      labs(
        x = col,
        y = "Density",
        color = "10-year bins",
        fill = "10-year bins"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.box.spacing = unit(0, "pt")
      ) +
      guides(
        color = guide_legend(nrow = 1, byrow = TRUE),
        fill  = guide_legend(nrow = 1, byrow = TRUE)
      )

    ### Density Plot – Entire Region by Season ------------------------------
    Density_all <- ggplot(plot_df, aes(
      x = AR,
      color = YearBin,
      fill = YearBin
    )) +
      geom_density(alpha = 0.05, adjust = 1) +
      scale_color_viridis_d(option = "C", direction = -1, name = "10-year bins") +
      scale_fill_viridis_d(option = "C", direction = -1, name = "10-year bins") +
      facet_wrap(~Season, ncol = 2, scales = "free_y") +
      labs(
        x = col,
        y = "Density"
      ) +
      theme_bw() +
      theme(
        legend.position = "top",
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.box.spacing = unit(0, "pt"),
        strip.text = element_text(size = 12, face = "plain") # matches your trend plots
      ) +
      guides(
        color = guide_legend(nrow = 1, byrow = TRUE),
        fill  = guide_legend(nrow = 1, byrow = TRUE)
      )

    ### Return list ----------------------------------------------------------
    list(Lines = AR_gg, Density = list(Indiv = Density_indivs, All = Density_all))
  })
  names(AR_gg) <- ar_cols
  AR_gg
})
names(AR_gg) <- MainVars$VarName

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

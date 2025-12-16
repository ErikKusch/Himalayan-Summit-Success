#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  ANALYSIS & VISUALIZATIONS: Statistical analysis and visualization of
#'      - mean climate change
#'  DEPENDENCIES:
#'  - "Himalayan Summit Success - DATA.R"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
source("Himalayan Summit Success - DATA.r")

message("Figure 2 Plotting")
FName <- file.path(Dir.Exports, "Figure2_.png")

# PLOTS ===================================================================

## Maps of Change ---------------------------------------------------------
symmetric_range <- function(x) {
    max_abs <- max(abs(range(x, na.rm = TRUE)))
    c(-max_abs, max_abs)
}

MapPlots <- lapply(1:nrow(MainVars), FUN = function(i) {
    # i <- 1
    ## extract data and get naming right
    RasterIter <- CDSData_ls[[MainVars$VarName[i]]]
    LegendTitle <- MainVars$ClearName[i]
    ClimateChangeTitle <- gsub("\\s*\\[.*?\\]", "", MainVars$ClearName[i])
    ## prepare data
    StartMean <- mean(RasterIter[[format(time(RasterIter), "%Y") %in% head(unique(format(time(RasterIter), "%Y")), 20)]])
    StopMean <- mean(RasterIter[[format(time(RasterIter), "%Y") %in% tail(unique(format(time(RasterIter), "%Y")), 20)]])
    ## cropping to tigher bounding box
    cropbox <- st_bbox(eightks_sf) + c(-0.5, -0.5, 0.5, 0.5)
    StartMean <- crop(StartMean, ext(c(cropbox$xmin, cropbox$xmax, cropbox$ymin, cropbox$ymax)))
    StopMean <- crop(StopMean, ext(c(cropbox$xmin, cropbox$xmax, cropbox$ymin, cropbox$ymax)))

    ## plotting
    Currentday_gg <- ggplot() +
        geom_spatraster(data = StopMean) +
        geom_sf(data = eightks_sf, shape = 2, colour = "white") +
        ggrepel::geom_text_repel(
            data = summits_df[summits_df$HEIGHTM >= 8000, ],
            aes(x = LON, y = LAT, label = PKNAME),
            max.overlaps = 30,
            colour = "white"
        ) +
        scale_fill_viridis_c(
            #  option = "C",
            name = LegendTitle,
            guide = guide_colourbar(title.vjust = 0.75)
        ) +
        labs(
            x = "Longitude",
            y = "Latitude",
            # title = paste("Average", LegendTitle, "During the Last 20 Years on Record"),
        ) +
        theme(
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.width = unit(2, "cm"),
            legend.key.height = unit(1, "cm"),
            panel.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c")
        )

    ClimateChange_gg <- ggplot() +
        geom_spatraster(data = StopMean - StartMean) +
        geom_sf(data = eightks_sf, shape = 2, colour = "white") +
        ggrepel::geom_text_repel(
            data = summits_df[summits_df$HEIGHTM >= 8000, ],
            aes(x = LON, y = LAT, label = PKNAME),
            max.overlaps = 30,
            colour = "white"
        ) +
        scale_fill_gradient2(
            low = "cyan",
            mid = "grey30",
            high = "red",
            midpoint = 0,
            limits = symmetric_range(values(StopMean - StartMean)), # This ensures symmetric color scale
            name = paste("Δ", LegendTitle),
            guide = guide_colourbar(title.vjust = 0.75)
        ) +
        labs(
            x = "Longitude",
            y = "Latitude",
            # title = paste(LegendTitle, "Change Between First and Last 20 Years on Record"),
        ) +
        theme(
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.width = unit(2, "cm"),
            legend.key.height = unit(1, "cm"),
            panel.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c"),
        )

    ## return
    plot_grid(Currentday_gg, ClimateChange_gg, ncol = 2)
})

## Density Plots per Season -----------------------------------------------
DensPlots <- lapply(1:nrow(MainVars), FUN = function(i) {
    # i = 1
    ### Prepare plot data -----
    plot_df <- peaks_ts_df[peaks_ts_df$Variable == MainVars$VarName[i], ] ## limit to variable
    plot_df <- plot_df[plot_df$PeakID %in% eightks_sf$PKNAME, ] ## limit to 8k peaks
    plot_df$Date <- as.POSIXct(plot_df$Date)
    plot_df$Year <- as.numeric(substr(plot_df$Date, 1, 4))
    plot_df <- plot_df %>%
        left_join(
            summits_df %>% select(PKNAME, ID),
            by = c("PeakID" = "PKNAME")
        )
    breaks <- seq(1051, 2021, by = 10)
    plot_df$YearBin <- cut(
        plot_df$Year,
        breaks = breaks,
        include.lowest = TRUE,
        right = FALSE,
        labels = paste0(breaks[-length(breaks)], ":", breaks[-1])
    )

    #### Individual peaks +++++
    Density_indivs <- ggplot(plot_df, aes(
        x = Value,
        color = YearBin,
        fill = YearBin
    )) +
        geom_density(alpha = 0.025, adjust = 1) +
        facet_grid(PeakID ~ Season, scales = "free_x") +
        scale_color_viridis_d(option = "C", direction = -1) +
        scale_fill_viridis_d(option = "C", direction = -1) +
        labs(
            x = MainVars$ClearName[i],
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

    #### Entire Region +++++
    plot_df$Season <- factor(plot_df$Season, levels = c("Pre", "Post"))
    plot_df$Season <- factor(
        plot_df$Season,
        levels = c("Pre", "Post"), # original values in your data
        labels = c("Pre-Monsoon", "Post-Monsoon") # desired facet labels
    )
    Density_all <- ggplot(plot_df, aes(x = Value, color = YearBin, fill = YearBin)) +
        geom_density(alpha = 0.05, adjust = 1) +
        # Continuous color scales for decades
        scale_color_viridis_d(option = "C", direction = -1, name = "10-year bins") +
        scale_fill_viridis_d(option = "C", direction = -1, name = "10-year bins") +
        # Facet by Season in one column
        facet_wrap(~Season, ncol = 2, scales = "free_y") +
        # Labels
        labs(
            x = MainVars$ClearName[i],
            y = "Density"
        ) +
        # Theme
        theme_bw() +
        theme(
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.direction = "horizontal",
            legend.justification = "center",
            legend.box.spacing = unit(0, "pt"),
            strip.text = element_text(size = 12, face = "bold")
        ) +
        guides(
            color = guide_legend(nrow = 1, byrow = TRUE),
            fill  = guide_legend(nrow = 1, byrow = TRUE)
        )

    ### Return plots -----
    list(Indiv = Density_indivs, All = Density_all)
})

## Linear Trendlines ------------------------------------------------------
LinePlots <- lapply(1:nrow(MainVars), FUN = function(i) {
    # i = 1
    ### Prepare plot data -----
    plot_df <- peaks_ts_df[peaks_ts_df$Variable == MainVars$VarName[i], ] ## limit to variable
    plot_df <- plot_df[plot_df$PeakID %in% eightks_sf$PKNAME, ] ## limit to 8k peaks
    plot_df$Date <- as.POSIXct(plot_df$Date)
    plot_df$Year <- as.numeric(substr(plot_df$Date, 1, 4))
    plot_df <- plot_df %>%
        left_join(
            summits_df %>% select(PKNAME, ID),
            by = c("PeakID" = "PKNAME")
        )
    breaks <- seq(1051, 2021, by = 10)
    plot_df$YearBin <- cut(
        plot_df$Year,
        breaks = breaks,
        include.lowest = TRUE,
        right = FALSE,
        labels = paste0(breaks[-length(breaks)], ":", breaks[-1])
    )

    ### Line plots of change over years -----
    #### Prepare Labels +++++
    # Prepare label information for each ID and Season
    labelInfo <- lapply(split(plot_df, list(plot_df$ID, plot_df$Season)), function(dat) {
        lmmod <- lm(Value ~ Year, data = dat)
        label_year <- if (unique(dat$Season) == "Pre") min(dat$Year) else max(dat$Year)
        Label <- predict(lmmod, newdata = data.frame(Year = label_year))
        data.frame(
            ID = unique(dat$ID),
            Season = unique(dat$Season),
            LabelPos = Label,
            LabelYear = label_year,
            Pval = summary(lmmod)$coefficients["Year", "Pr(>|t|)"],
            Effect = summary(lmmod)$coefficients["Year", "Estimate"]
        )
    })
    labelInfo <- do.call(rbind, labelInfo)
    labelInfo$Fill <- as.character(labelInfo$Pval < 0.5)

    season_trends <- lapply(unique(plot_df$Season), function(season) {
        dat <- plot_df %>% filter(Season == season)
        lmmod <- lm(Value ~ Year, data = dat)
        # Label positions: start for pre, end for post
        label_year <- if (season == "Pre") min(dat$Year) else max(dat$Year)
        Label <- predict(lmmod, newdata = data.frame(Year = label_year))
        data.frame(
            ID = "Region",
            Season = season,
            LabelPos = Label,
            LabelYear = label_year,
            Pval = summary(lmmod)$coefficients["Year", "Pr(>|t|)"],
            Fill = "Region",
            Effect = summary(lmmod)$coefficients["Year", "Estimate"]
        )
    })
    season_trends <- do.call(rbind, season_trends)
    labelInfo <- rbind(labelInfo, season_trends)

    # Merge labels back to main dataframe for significance
    plot_df <- plot_df %>%
        left_join(labelInfo, by = c("ID", "Season")) %>%
        mutate(StatSig = Pval < 0.05)
    plot_df$Effect <- as.numeric(plot_df$Effect)
    labelInfo$Effect <- as.numeric(labelInfo$Effect)

    #### Actual Plot +++++
    Return_gg <- ggplot(plot_df, aes(x = Year, y = Value, group = interaction(ID, Season))) +
        # Trendlines for each peak-season
        geom_smooth(aes(linetype = StatSig, col = Effect), method = "lm", alpha = 0.3, linewidth = 0.8, show.legend = TRUE) +
        # Overall trendlines per season
        geom_smooth(aes(group = Season), method = "lm", linewidth = 2, col = "#1a1b1b", show.legend = FALSE) +
        # Labels for each ID
        geom_label_repel(
            data = labelInfo,
            aes(x = LabelYear, y = LabelPos, label = ID, fill = Effect),
            color = "#000000",
            size = 8,
            hjust = ifelse(labelInfo$Season == "Pre", 1.5, -1),
            direction = "y",
            nudge_x = 0,
            show.legend = FALSE
        ) +
        # Linetype legend for StatSig
        scale_linetype_manual(
            name = "Statistical Significance of Trend",
            values = c(`FALSE` = "dashed", `TRUE` = "solid"),
            guide = guide_legend(
                title.position = "top",
                label.position = "bottom",
                nrow = 1,
                override.aes = list(
                    size = 1.5, # thicker lines
                    color = "black" # force legend lines to black
                ),
                keywidth = 4,
                keyheight = 1,
                order = 2 # shows second in legend
            )
        ) +
        # Continuous color bar for Effect
        scale_color_viridis_c(
            direction = -1,
            name = "Mean Change Year over Year",
            guide = guide_colorbar(
                title.position = "top",
                barwidth = 15,
                barheight = 1.5,
                order = 1 # shows first in legend
            ),
        ) +
        # Fill scale for labels (optional if you want label fill colors)
        scale_fill_viridis_c(direction = -1) +
        # Axis labels
        labs(
            x = "Year",
            y = MainVars$ClearName[i]
        ) +
        # Plot limits
        lims(x = c(min(plot_df$Year) - 5, max(plot_df$Year) + 5)) +
        # Theme and legend layout
        theme_bw() +
        theme(
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.direction = "horizontal",
            legend.justification = "center",
            legend.box.spacing = unit(5, "pt") # small spacing between legends
        ) +
        # Remove unwanted fill legend for labels
        guides(
            fill = "none"
        )

    ### Return plots -----
    Return_gg
})

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
label_A <- label_row("(A) Air Temperature Change")
label_B <- label_row("(B) Windspeed Change")
label_C <- label_row("(C) Snow Cover Change")

main_ggs <- plot_grid(
    label_A,
    plot_grid(
        MapPlots[[1]],
        DensPlots[[1]]$All,
        nrow = 2,
        rel_heights = c(1, 0.8)
    ),
    label_B,
    plot_grid(
        MapPlots[[2]],
        DensPlots[[2]]$All,
        nrow = 2,
        rel_heights = c(1, 0.8)
    ),
    label_C,
    plot_grid(
        MapPlots[[3]],
        DensPlots[[3]]$All,
        nrow = 2,
        rel_heights = c(1, 0.8)
    ),
    ncol = 1,
    rel_heights = c(
        0.08, 1, # A label + plot
        0.08, 1, # B label + plot
        0.08, 1 # C label + plot
    )
)
ggsave(
    main_ggs,
    file = FName,
    width = 14, height = 22
)

### Supplement +++++++
ggsave(
    plot_grid(plotlist = LinePlots, ncol = 1),
    file = paste0(tools::file_path_sans_ext(FName), "_Supplement_LineTrends.png"),
    width = 16, height = 22
)

lapply(1:length(DensPlots), FUN = function(x) {
    ggsave(
        DensPlots[[x]]$Indiv,
        file = paste0(tools::file_path_sans_ext(FName), "_Supplement_Density_", MainVars$VarName[[x]], ".png"),
        width = 16, height = 22
    )
})

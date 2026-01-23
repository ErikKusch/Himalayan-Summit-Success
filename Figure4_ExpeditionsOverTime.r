#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  ANALYSIS & VISUALIZATIONS: Statistical analysis and visualization of
#'      - Expeditions throughout time and seasons
#'  DEPENDENCIES:
#'  - "Himalayan Summit Success - DATA.R"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
source("Himalayan Summit Success - DATA.r")

message("Figure 4 Plotting")
FName <- file.path(Dir.Exports, "Figure4_.png")

# PLOTS ===================================================================
Data_ls <- list(
    Main = ModelData_df[ModelData_df$PKNAME %in% eightks_sf$PKNAME, ],
    Supplement_AllPeaks = ModelData_df
)

### Expeditions over Time -------------------------------------------------
# Calculate expedition duration and create year column
Plot_ls <- lapply(names(Data_ls), FUN = function(Name) {
    ModelData_df <- Data_ls[[Name]]
    ModelData_df$Duration <- as.numeric(difftime(as.Date(ModelData_df$TERMDATE),
        as.Date(ModelData_df$BCDATE),
        units = "days"
    ))
    ModelData_df$Year <- as.numeric(format(as.Date(ModelData_df$BCDATE), "%Y"))

    # Create daily presence data with member counts
    expedition_days <- do.call(rbind, lapply(1:nrow(ModelData_df), function(i) {
        dates <- seq(as.Date(ModelData_df$BCDATE[i]),
            as.Date(ModelData_df$TERMDATE[i]),
            by = "day"
        )
        data.frame(
            Date = dates,
            Year = as.numeric(format(dates, "%Y")),
            DayOfYear = as.numeric(format(dates, "%j")),
            Members = ModelData_df$TOTMEMBERS[i] # Add member count
        )
    }))

    # Create complete grid of years and days
    all_years <- 1951:2021
    all_days <- 1:365
    complete_grid <- expand.grid(Year = all_years, DayOfYear = all_days)

    # Count expeditions and total members per day
    expedition_counts <- aggregate(
        cbind(Expeditions = 1, Members = Members) ~ Year + DayOfYear,
        data = expedition_days,
        FUN = sum
    )
    expedition_counts$Year <- as.numeric(as.character(expedition_counts$Year))
    expedition_counts$DayOfYear <- as.numeric(as.character(expedition_counts$DayOfYear))

    # Merge with complete grid to include all years and days
    expedition_counts <- merge(
        complete_grid,
        expedition_counts,
        by = c("Year", "DayOfYear"),
        all.x = TRUE
    )

    # Replace NAs with zeros
    expedition_counts$Expeditions[is.na(expedition_counts$Expeditions)] <- 0
    expedition_counts$Members[is.na(expedition_counts$Members)] <- 0

    # Reshape data to long format for faceting
    expedition_counts <- tidyr::pivot_longer(
        expedition_counts,
        cols = c(Expeditions, Members),
        names_to = "type",
        values_to = "value"
    )
    expedition_counts <- expedition_counts[expedition_counts$value != 0, ]

    # Create separate plots for expeditions and members
    ExpeditionPlots_ls <- lapply(unique(expedition_counts$type), function(plot_type) {
        data_subset <- expedition_counts[expedition_counts$type == plot_type, ]

        # Create title based on type
        title <- if (plot_type == "Expeditions") {
            "Number of Active Expeditions Throughout the Year"
        } else {
            "Total Number of Expedition Members Throughout the Year"
        }

        ggplot(
            data_subset,
            aes(
                x = DayOfYear,
                y = as.factor(Year),
                fill = value
            )
        ) +
            # Add vertical lines for pre and post monsoon seasons
            geom_vline(xintercept = c(60, 151), color = "white", alpha = 0.3, linetype = "dashed") + # March (60) to May (151)
            geom_vline(xintercept = c(244, 304), color = "white", alpha = 0.3, linetype = "dashed") + # September (244) to November (304)
            geom_tile() +
            scale_fill_viridis_c(
                name = "", # plot_type,
                option = "C"
            ) +
            scale_y_discrete(
                name = "Year",
                limits = rev(as.character(1951:2021)), # Reverse to show earliest at bottom
                breaks = rev(as.character(seq(1951, 2021, by = 5))) # Show only every 5th year
            ) +
            scale_x_continuous(
                name = "Month",
                breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
                limits = c(1, 365), # Set limits to show full year
                expand = c(0, 0) # Remove spacing at edges
            ) +
            # labs(
            #     title = title,
            #     subtitle = "Vertical lines indicate pre-monsoon and post-monsoon seasons"
            # ) +
            theme_bw() +
            theme(
                legend.position = "right",
                legend.key.width = unit(0.5, "cm"),
                legend.key.height = unit(2, "cm"),
                panel.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c"),
                panel.grid = element_line(color = "grey30"),
                text = element_text(color = "white"),
                axis.text = element_text(color = "white"),
                plot.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c"),
                legend.background = element_rect(fill = "#2c2c2c"),
                legend.text = element_text(color = "white"),
                legend.title = element_text(color = "white"),
                panel.spacing = unit(0, "lines"), # Remove spacing between tiles
                axis.ticks = element_blank() # Remove axis ticks for cleaner look
            )
    })
    # Name the list elements
    names(ExpeditionPlots_ls) <- unique(expedition_counts$type)

    ### Proportions of Expeditions over Time ----------------------------------

    SeasonCounts_df <- ModelData_df %>%
        count(YEAR, SEASON) %>%
        pivot_wider(
            names_from = SEASON,
            values_from = n,
            names_prefix = "Season_",
            values_fill = 0
        )
    SeasonCounts_df$TOTAL <- SeasonCounts_df$Season_1 + SeasonCounts_df$Season_3
    SeasonCounts_df$PROPORTION <- SeasonCounts_df$Season_3 / SeasonCounts_df$TOTAL
    # saveRDS(SeasonCounts_df, file = paste0("SeasonExpeds_", Name, ".rds"))


    ## raw proportion plot
    SeasonsExpedsProp <- ggplot(
        SeasonCounts_df,
        aes(x = YEAR, y = PROPORTION)
    ) +
        geom_point() +
        ggrepel::geom_label_repel(aes(label = TOTAL),
            box.padding   = 0.35,
            point.padding = 0.5,
            segment.color = "grey50"
        ) +
        geom_smooth(method = "loess", fill = "#156082", col = "#156082") +
        theme_bw() +
        labs(y = paste("Proportion of Post-Monsoon Season Expeditions"), x = "Year")

    ## changepoint analysis
    clean_df <- SeasonCounts_df[!is.na(SeasonCounts_df$PROPORTION), ]
    data_ranks <- rank(clean_df$PROPORTION)
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

    SeasonsExpeds <- ggplot(
        SeasonCounts_df,
        aes(x = YEAR, y = PROPORTION)
    ) +
        geom_point() +
        # geom_line() +
        # Add horizontal segments for medians (equivalent to segments())
        geom_segment(
            aes(x = min(clean_df$YEAR), xend = cp_year, y = median_before, yend = median_before),
            color = "blue", linewidth = 1.5 # Note: linewidth instead of lwd for ggplot2
        ) +
        geom_segment(
            aes(x = cp_year, xend = max(clean_df$YEAR), y = median_after, yend = median_after),
            color = "blue", linewidth = 1.5
        ) +
        # Add vertical line at change point (equivalent to abline(v = cp_year))
        geom_vline(xintercept = cp_year, color = "red", linewidth = 1, linetype = "dashed") +
        # Add label for change point (equivalent to text())
        annotate(
            "text",
            x = cp_year, y = max(clean_df$PROPORTION),
            label = paste("Break:", cp_year),
            hjust = 1.2,
            color = "red"
        ) +
        ggrepel::geom_label_repel(
            aes(label = TOTAL),
            box.padding = 0.35,
            point.padding = 0.5,
            segment.color = "grey50"
        ) +
        theme_bw() +
        labs(y = "Proportion of Post-Monsoon Season Expeditions", x = "Year")

    ## return
    list(
        ExpeditionPlots_ls = ExpeditionPlots_ls,
        SeasonsExpeds = SeasonsExpeds,
        SeasonsExpedsProp = SeasonsExpedsProp
    )
})
names(Plot_ls) <- names(Data_ls)

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

label_A <- label_row("(A) Expeditions with established Basecamp")
label_B <- label_row("(B) Proportion of Pre- vs. Post-Monsoon Season Expeditions")
label_C <- label_row("(C) Number of Mountaineers with established Basecamp")

lapply(names(Plot_ls), FUN = function(Name) {
    # print(Name)
    ExpeditionPlots_ls <- Plot_ls[[Name]]$ExpeditionPlots_ls
    main_ggs <- plot_grid(
        label_A,
        ExpeditionPlots_ls$Expeditions,
        label_B,
        Plot_ls[[Name]]$SeasonsExpeds,
        label_C,
        ExpeditionPlots_ls$Members,
        ncol = 1,
        rel_heights = c(
            0.08, 1, # A label + plot
            0.08, 1, # B label + plot
            0.08, 1 # C label + plot
        )
    )
    Name <- paste0(Name, "_")
    FName <- paste0(
        tools::file_path_sans_ext(FName),
        paste0(
            gsub(pattern = "Main_", replacement = "", Name), ".png"
        )
    )
    ggsave(
        main_ggs,
        file = FName,
        width = 14, height = 22
    )
})

FUN.BayesianPlot <- function(plot_ls, TotalExped = FALSE, ScaleFac = 1) {
    df <- plot_ls$df
    cols <- plot_ls$cols
    columns <- plot_ls$columns
    Names <- plot_ls$names
    colnames(df)[match(columns, colnames(df))] <- c("Var1", "Var2")

    ### Supplement (for each peakid) -------
    Peaks_ls <- lapply(unique(df$PEAKID), FUN = function(k) {
        dfiter <- df[df$PEAKID == k, ]
        scale_factor <- max(dfiter$Var1, na.rm = TRUE) / max(dfiter$Var2, na.rm = TRUE) / ScaleFac
        Supp_gg <- ggplot(dfiter, aes(x = year_0 + 1951)) +
            # Expedition Size (ES) - primary axis
            stat_lineribbon(
                aes(
                    y = Var1,
                    fill = Names[1],
                    color = Names[1],
                    linetype = Names[1]
                ),
                .width = c(.95, 0.8, 0.5),
                alpha = 0.2,
                size = 1
            ) +
            # Summit Bid Window (SB) - secondary axis, scaled
            stat_lineribbon(
                aes(
                    y = Var2 * scale_factor,
                    fill = Names[2],
                    color = Names[2],
                    linetype = Names[2]
                ),
                .width = c(.95, 0.8, 0.5),
                alpha = 0.2,
                size = 1
            ) +
            # Dual y-axis
            scale_y_continuous(
                name = Names[1],
                sec.axis = sec_axis(~ . / scale_factor, name = Names[2])
            ) +
            # Manual legend (colors + linetype)
            scale_color_manual(
                name = "",
                values = setNames(cols[1:2], Names[1:2])
            ) +
            scale_fill_manual(
                name = "",
                values = setNames(cols[1:2], Names[1:2])
            ) +
            scale_linetype_manual(
                name = "",
                values = setNames(c("solid", "dashed"), Names[1:2])
            ) +
            # Facets per PEAKID
            facet_wrap(~PEAKID, scales = "free_y", ncol = 3) +
            # Labels
            labs(x = "Year") +
            # Theme
            theme_bw() +
            theme(
                axis.title.y = element_text(color = cols[1]),
                axis.title.y.right = element_text(color = cols[2]),
                legend.position = "top",
                legend.direction = "horizontal",
                legend.key.width = unit(2, "cm"),
                strip.text = element_text(face = "bold")
            )
        Supp_gg
    })
    # Fuse into one plot
    Peaks_clean <- lapply(Peaks_ls, function(p) {
        p + theme(legend.position = "none")
    })
    panel <- plot_grid(plotlist = Peaks_clean, ncol = 3, align = "v") # Vertical stack
    legend <- get_legend(Peaks_ls[[1]]) # Recover legend from first plot
    Supp_gg <- plot_grid(panel, legend, ncol = 1, rel_heights = c(1, 0.02)) # Add outer legend to stack
    # Supp_gg

    ### Main Text (for each peakid) -------
    if (TotalExped) { # TotalExped flag fixes some scaling for the ExpedSize panel
        df <- df %>%
            ungroup() %>%
            group_by(.draw, year_0) %>%
            summarise(
                Var1 = mean(Var1, na.rm = TRUE),
                Var2 = sum(Var2, na.rm = TRUE),
                .groups = "drop"
            )
    }
    scale_factor <- max(df$Var1, na.rm = TRUE) / max(df$Var2, na.rm = TRUE) / ScaleFac # scale factor for dual y-axis
    Main_gg <- ggplot(df, aes(x = year_0 + 1951)) +
        # Expedition Size (ES) - primary axis
        stat_lineribbon(
            aes(
                y = Var1,
                fill = Names[1],
                color = Names[1],
                linetype = Names[1]
            ),
            .width = c(.95, 0.8, 0.5),
            alpha = 0.2,
            size = 1
        ) +
        # Summit Bid Window (SB) - secondary axis, scaled
        stat_lineribbon(
            aes(
                y = Var2 * scale_factor,
                fill = Names[2],
                color = Names[2],
                linetype = Names[2]
            ),
            .width = c(.95, 0.8, 0.5),
            alpha = 0.2,
            size = 1
        ) +
        # Dual y-axis
        scale_y_continuous(
            name = Names[1],
            sec.axis = sec_axis(~ . / scale_factor, name = Names[2])
        ) +
        # Manual legend (colors + linetype)
        scale_color_manual(
            name = "",
            values = setNames(cols[1:2], Names[1:2])
        ) +
        scale_fill_manual(
            name = "",
            values = setNames(cols[1:2], Names[1:2])
        ) +
        scale_linetype_manual(
            name = "",
            values = setNames(c("solid", "dashed"), Names[1:2])
        ) +
        # Labels
        labs(x = "Year") +
        # Theme
        theme_bw() +
        guides(
            color = guide_legend(ncol = 2),
            fill = guide_legend(ncol = 2),
            linetype = guide_legend(ncol = 2)
        ) +
        theme(
            axis.title.y = element_text(color = cols[1]),
            axis.title.y.right = element_text(color = cols[2]),
            legend.position = "top",
            legend.direction = "vertical",
            legend.box = "vertical",
            legend.key.width = unit(2, "cm"),
            strip.text = element_text(face = "bold")
        )
    # Main_gg

    ### Return -------
    list(Main = Main_gg, Supp = Supp_gg)
}

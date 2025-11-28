#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  ANALYSIS & VISUALIZATIONS: Statistical analysis and visualization of
#'      - summit bid window size
#'      - expedition size
#'      - mortality and death rates
#'  DEPENDENCIES:
#'  - "Himalayan Summit Success - DATA.R"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
source("Himalayan Summit Success - DATA.r")

message("Figure 1 Plotting")
FName <- file.path(Dir.Exports, "Figure1_.png")

# PREPARATION =============================================================
df_data <- ModelData_df
df_data$year_0 <- df_data$YEAR - min(df_data$YEAR) # add year 0 as starting year
df_data <- df_data %>%
    filter(HEIGHTM >= 8000)
breaks <- seq(1051, 2021, by = 10)
df_data$YearBin <- cut(
    df_data$YEAR,
    breaks = breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = paste0(breaks[-length(breaks)], ":", breaks[-1])
)

# MODELS ==================================================================

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
summary(model_TE)
pp_check(model_TE, ndraws = 500)
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
summary(model_TM)
pp_check(model_TM, ndraws = 500)
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
summary(model_ES)
pp_check(model_ES, ndraws = 500)
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
summary(model_SB)
pp_check(model_SB, ndraws = 500)
conditions <- expand.grid(
    PEAKID = sort(unique(df_filtered$PEAKID)),
    year_0 = sort(unique(df_filtered$year_0))
)
rownames(conditions) <- NULL
SB_draws <- model_SB %>%
    add_epred_draws(newdata = conditions, re_formula = NULL) %>%
    rename(SB = .epred)


## Mortality & Death Rates ------------------------------------------------
df_filtered <- df_data
if (file.exists(file.path(Dir.Exports, "model_MT.RData"))) {
    load(file.path(Dir.Exports, "model_MT.RData"))
} else {
    model_MT <- brm(
        bf(
            Mortality_all ~ year_0 + (year_0 | PEAKID), # random intercepts & slopes per peak
            phi ~ 1, # constant precision
            zoi ~ year_0 + (year_0 | PEAKID) # model zero/one inflation similarly
        ),
        data = df_filtered,
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
summary(model_MT)
pp_check(model_MT, ndraws = 500)
conditions <- expand.grid(
    PEAKID = sort(unique(df_filtered$PEAKID)),
    year_0 = sort(unique(df_filtered$year_0))
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
# Combine draws by .draw, year_0, PEAKID
# combined_draws <- left_join(
#     mortality_draws %>% select(PEAKID, year_0, .draw, mort_rate),
#     hurdle_draws %>% select(PEAKID, year_0, .draw, prob_mort),
#     by = c("PEAKID", "year_0", ".draw")
# )


# PLOTTING ================================================================
## Preparing Plotting Data and Colours ------------------------------------
Dataframes <- list(
    ExpedSummitBid = list(
        df = left_join(
            ES_draws %>% select(PEAKID, year_0, .draw, ES),
            SB_draws %>% select(PEAKID, year_0, .draw, SB),
            by = c("PEAKID", "year_0", ".draw")
        ),
        cols = c("#003575", "#754000"),
        columns = c("ES", "SB"),
        names = c("Individual Expedition Sizes [#]", "Summit Bid Window [days]")
    ),
    ExpedSizes = list(
        df = left_join(
            TM_draws %>% select(PEAKID, year_0, .draw, TM),
            TE_draws %>% select(PEAKID, year_0, .draw, TE),
            by = c("PEAKID", "year_0", ".draw")
        ),
        cols = c("#6F0075", "#067500"),
        columns = c("TM", "TE"),
        names = c("Total Expedition Sizes [#]", "Total Number of Expeditions [#]")
    ),
    MortDeath = list(
        df = left_join(
            mortality_draws %>% select(PEAKID, year_0, .draw, mort_rate),
            hurdle_draws %>% select(PEAKID, year_0, .draw, prob_mort),
            by = c("PEAKID", "year_0", ".draw")
        ),
        cols = c("#006D75", "#750800"),
        columns = c("mort_rate", "prob_mort"),
        names = c("Probability of Any Mortality [%]", "Conditional Mortality [%]")
    )
)

## Plot Generation --------------------------------------------------------
plotlist <- lapply(1:length(Dataframes), FUN = function(i) {
    # i = 2
    df <- Dataframes[[i]]$df
    cols <- Dataframes[[i]]$cols
    columns <- Dataframes[[i]]$columns
    Names <- Dataframes[[i]]$names
    colnames(df)[match(columns, colnames(df))] <- c("Var1", "Var2")

    ### Supplement (for each peakid) -------
    Peaks_ls <- lapply(unique(df$PEAKID), FUN = function(k) {
        dfiter <- df[df$PEAKID == k, ]
        scale_factor <- max(dfiter$Var1, na.rm = TRUE) / max(dfiter$Var2, na.rm = TRUE)
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
                name = "Estimate Type",
                values = setNames(cols[1:2], Names[1:2])
            ) +
            scale_fill_manual(
                name = "Estimate Type",
                values = setNames(cols[1:2], Names[1:2])
            ) +
            scale_linetype_manual(
                name = "Estimate Type",
                values = setNames(c("solid", "dashed"), Names[1:2])
            ) +
            # Facets per PEAKID
            facet_wrap(~PEAKID, scales = "free_y", nrow = 4) +
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
    panel <- plot_grid(plotlist = Peaks_clean, ncol = 1, align = "v") # Vertical stack
    legend <- get_legend(Peaks_ls[[1]]) # Recover legend from first plot
    Supp_gg <- plot_grid(panel, legend, ncol = 1, rel_heights = c(1, 0.02)) # Add outer legend to stack
    # Supp_gg

    ### Main Text (for each peakid) -------
    if (i == 2) {
        df <- df %>%
            ungroup() %>%
            group_by(.draw, year_0) %>%
            summarise(
                Var1 = sum(Var1, na.rm = TRUE),
                Var2 = sum(Var2, na.rm = TRUE),
                .groups = "drop"
            )
    }
    scale_factor <- max(df$Var1, na.rm = TRUE) / max(df$Var2, na.rm = TRUE)
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
            name = "Estimate Type",
            values = setNames(cols[1:2], Names[1:2])
        ) +
        scale_fill_manual(
            name = "Estimate Type",
            values = setNames(cols[1:2], Names[1:2])
        ) +
        scale_linetype_manual(
            name = "Estimate Type",
            values = setNames(c("solid", "dashed"), Names[1:2])
        ) +
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
    # Main_gg

    ### Return -------
    list(Main = Main_gg, Supp = Supp_gg)
})

## Plot Fusing & Saving ---------------------------------------------------
### Supplement -------
combined <- plot_grid(plotlist = lapply(plotlist, "[[", "Supp"), nrow = 1)
ggsave(
    combined,
    file = paste0(tools::file_path_sans_ext(FName), "_Supplement.png"),
    width = 20, height = 45
)

### Main Text ------
img <- readPNG(file.path(Dir.Exports, "AreaMap.png")) # replace with your file path
Map_png <- rasterGrob(img, interpolate = TRUE) # convert to a grob
combined <- plot_grid(
    plot_grid(plotlist = lapply(plotlist, "[[", "Main"), nrow = 1),
    Map_png,
    ncol = 1
)
ggsave(
    combined,
    file = FName,
    width = 20, height = 16
)


















# ## Summit Bid Window & Number of Expeditions ------------------------------
# ### Main Text --------
# # entire region, totals have to be summed

# ### Supplement --------
# # individual peaks



# ## Size of Indvidual Expeditions & Total Size of all Expeditions ----------

# ### Main Text --------
# # entire region, totals have to be summed

# ### Supplement --------
# # individual peaks





# ## Mortality & Death Rate Panel -------------------------------------------
# ### Main Text --------
# # entire region, totals have to be summed

# ### Supplement --------
# # individual peaks














# ## Summit Bid Window & Expedition Size Panel ------------------------------
# ### Main Text Plot  ------------
# # Combine draws by .draw, year_0, PEAKID
# combined_draws2 <- left_join(
#     ES_draws %>% select(PEAKID, year_0, .draw, ES),
#     SB_draws %>% select(PEAKID, year_0, .draw, SB),
#     by = c("PEAKID", "year_0", ".draw")
# )

# # Compute scaling factor
# scale_factor <- max(combined_draws2$ES, na.rm = TRUE) / max(combined_draws2$SB, na.rm = TRUE) / 1.5

# BidMembers_gg <- ggplot(combined_draws2, aes(x = year_0 + 1951)) +
#     # Expedition Size (ES) - primary axis
#     stat_lineribbon(
#         aes(
#             y = ES,
#             fill = "Expedition Size [#]",
#             color = "Expedition Size [#]",
#             linetype = "Expedition Size [#]"
#         ),
#         .width = c(.95, 0.8, 0.5),
#         alpha = 0.2,
#         size = 1
#     ) +
#     # Summit Bid Window (SB) - secondary axis, scaled
#     stat_lineribbon(
#         aes(
#             y = SB / scale_factor,
#             fill = "Summit Bid Window [days]",
#             color = "Summit Bid Window [days]",
#             linetype = "Summit Bid Window [days]"
#         ),
#         .width = c(.95, 0.8, 0.5),
#         alpha = 0.2,
#         size = 1
#     ) +
#     # Dual y-axis
#     scale_y_continuous(
#         name = "Expedition Size [#]",
#         sec.axis = sec_axis(~ . * scale_factor, name = "Summit Bid Window [days]")
#     ) +
#     # Manual legend (colors + linetype)
#     scale_color_manual(
#         name = "Estimate Type",
#         values = c("Expedition Size [#]" = "#961765", "Summit Bid Window [days]" = "#8a8b21")
#     ) +
#     scale_fill_manual(
#         name = "Estimate Type",
#         values = c("Expedition Size [#]" = "#961765", "Summit Bid Window [days]" = "#8a8b21")
#     ) +
#     scale_linetype_manual(
#         name = "Estimate Type",
#         values = c("Expedition Size [#]" = "solid", "Summit Bid Window [days]" = "dashed")
#     ) +
#     # Facets per PEAKID
#     facet_wrap(~PEAKID, scales = "free_y", nrow = 4) +
#     # Labels
#     labs(x = "Year") +
#     # Theme
#     theme_bw() +
#     theme(
#         axis.title.y = element_text(color = "#961765"),
#         axis.title.y.right = element_text(color = "#8a8b21"),
#         legend.position = "top",
#         legend.direction = "horizontal",
#         legend.key.width = unit(2, "cm"),
#         strip.text = element_text(face = "bold")
#     )
# BidMembers_gg

# ### Supplement Plot  ------------

# ## Mortality & Death Rate Panel -------------------------------------------
# ### Main Text Plot  ------------
# # Scale factor for secondary axis
# scale_factor <- max(combined_draws$mort_rate / 2, na.rm = TRUE)

# MortDeaths_gg <- ggplot(combined_draws, aes(x = year_0 + 1951)) +
#     # Primary y-axis: probability of any mortality (blue)
#     stat_lineribbon(
#         aes(
#             y = prob_mort,
#             fill = "Probability of Any Mortality",
#             color = "Probability of Any Mortality",
#             linetype = "Probability of Any Mortality"
#         ),
#         .width = c(.95, 0.8, 0.5),
#         alpha = 0.2,
#         size = 1
#     ) +
#     # Secondary y-axis: conditional mortality (orange), scaled
#     stat_lineribbon(
#         aes(
#             y = mort_rate / scale_factor,
#             fill = "Conditional Mortality (given >0)",
#             color = "Conditional Mortality (given >0)",
#             linetype = "Conditional Mortality (given >0)"
#         ),
#         .width = c(.95, 0.8, 0.5),
#         alpha = 0.2,
#         size = 1
#     ) +
#     # Dual y-axis
#     scale_y_continuous(
#         name = "Probability of Any Mortality",
#         sec.axis = sec_axis(~ . * scale_factor, name = "Conditional Mortality (given >0)")
#     ) +
#     # Manual legend
#     scale_color_manual(
#         name = "Estimate Type",
#         values = c(
#             "Probability of Any Mortality" = "#000655",
#             "Conditional Mortality (given >0)" = "#6b3400"
#         )
#     ) +
#     scale_fill_manual(
#         name = "Estimate Type",
#         values = c(
#             "Probability of Any Mortality" = "#000655",
#             "Conditional Mortality (given >0)" = "#6b3400"
#         )
#     ) +
#     scale_linetype_manual(
#         name = "Estimate Type",
#         values = c(
#             "Probability of Any Mortality" = "solid",
#             "Conditional Mortality (given >0)" = "dashed"
#         )
#     ) +
#     # Facets
#     facet_wrap(~PEAKID, scales = "free_y") +
#     # Axis labels
#     labs(x = "Year") +
#     # Theme
#     theme_bw() +
#     theme(
#         axis.title.y = element_text(color = "#000655"),
#         axis.title.y.right = element_text(color = "#6b3400"),
#         legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.key.width = unit(2, "cm"),
#         strip.text = element_text(face = "bold")
#     )
# MortDeaths_gg

# ### Supplement Plot  ------------

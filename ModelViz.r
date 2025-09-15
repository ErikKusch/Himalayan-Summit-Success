# model__1 <- brms::brm(
#   STORM ~ YEAR +
#     +YEAR:PKNAME +
#     (1 | PKNAME),
#   data = df_filtered,
#   family = bernoulli(link = "logit"),
#   # family = poisson(),
#   chains = 4
#   # cores = parallel::detectCores(),
#   # iter = 10000,
#   # warmup = 5000
# )


# save(model__1, file = file.path(Dir.Exports, "Model__1.RData"))

# load(file.path(Dir.Exports, "Model__1.RData"))

model <- model__1


# View(model__1)



Summits <- unique(model$data$PKNAME)
Years <- min(model$data$YEAR):max(model$data$YEAR)

mod_pred <- add_fitted_draws(
  model = model,
  newdata = data.frame(
    PKNAME = rep(Summits, each = length(Years)),
    YEAR = rep(Years, length(Summits))
  ),
  n = 10
)

labelInfo <- split(mod_pred, mod_pred$PKNAME)
labelInfo <- pblapply(labelInfo, function(dat) {
  ELabel <- predict(loess(.value ~ YEAR, span = 0.8, data = dat), newdata = data.frame(YEAR = max(dat$YEAR)))
  summit <- unique(dat$PKNAME)
  data.frame(ELabel = ELabel, Summit = summit)
})
labelInfo <- do.call(rbind, labelInfo)

Line_gg <- ggplot(mod_pred, aes(x = YEAR, y = .value)) +
  geom_smooth(aes(group = PKNAME), col = "#535353", alpha = 0.2, method = "loess", span = 0.8) +
  geom_smooth(fill = "#156082", col = "#156082", method = "loess", span = 0.8) +
  geom_label_repel(
    data = labelInfo[labelInfo$ELabel != 0, ],
    aes(
      x = max(mod_pred$YEAR), y = ELabel,
      label = Summit
    ),
    color = "#535353",
    nudge_x = 7
  ) +
  theme_bw() +
  labs(x = "Year", y = "Porbability of Reported Storm") +
  lims(y = c(0, 1))


post_pars_lm1 <- get_variables(model)
mylist_lm1_b1 <- list(model, as.name("b_YEAR"))

# Obtain population-level parameter
pooled_lm1_b1 <- do.call(spread_draws, mylist_lm1_b1)
pooled_lm1_b1 <- pooled_lm1_b1 %>%
  mutate(
    PKNAME = "ALL",
    param = "Intercept",
    b_YEAR = b_YEAR
  )

b_plot_df <- model %>%
  spread_draws(b_YEAR, r_PKNAME[PKNAME, param]) %>%
  filter(param == "Intercept") %>%
  bind_rows(pooled_lm1_b1)

StatSig <- aggregate(b_YEAR ~ PKNAME, b_plot_df, FUN = quantile, c(0.05, 0.95))
StatSig$Direction <- sign(StatSig$b_YEAR[, 1]) + sign(StatSig$b_YEAR[, 2])
StatSig$Sig <- abs(StatSig$Direction) == 2
b_plot_df$StatSig <- NA
for (i in 1:nrow(b_plot_df)) {
  if (length(which(StatSig$Summit == b_plot_df$Summit[i])) == 0) {
    b_plot_df$StatSig[i] <- FALSE
  } else {
    b_plot_df$StatSig[i] <- StatSig$Sig[intersect(which(StatSig$Outcome == b_plot_df$Outcome[i]), which(StatSig$Summit == b_plot_df$Summit[i]))]
  }
}
b_plot_df$col <- ifelse(b_plot_df$StatSig, "green", "red")
b_plot_df$fill <- ifelse(b_plot_df$PKNAME == "ALL", "#156082", "#535353")

BExceeded_gg <-
  ggplot(
    b_plot_df,
    aes(
      y = factor(PKNAME, levels = rev(c("ALL", rev(labelInfo$Summit[order(labelInfo$ELabel)])))),
      x = b_YEAR
    )
  ) +
  stat_halfeye(aes(fill = fill)) +
  scale_fill_manual(values = c("#156082", "#535353")) +
  geom_boxplot(aes(col = col), width = 0.3, lwd = 1.1) +
  scale_color_manual(values = c("#003b05", "#5c0000"), breaks = c("green", "red")) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "BRMS Model Coefficient Posterior Samples", y = "")

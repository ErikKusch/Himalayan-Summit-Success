#' ####################################################################### #
#' PROJECT: [Himalayan Summit Success - Climate Data]
#' CONTENTS:
#'  ANALYSIS: Mortality caused by objective hazards (crevasse, icefall collapse, avalanche, falling rock/ice)
#'  DEPENDENCIES:
#'  - "Himalayan Summit Success - DATA.R"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
source("Himalayan Summit Success - DATA.r")

# PREPARATION =============================================================
memb <- read.dbf("Data/members.DBF", as.is = TRUE)
peaks <- read.dbf("Data/peaks.DBF", as.is = TRUE)

# ANALYSIS ================================================================
sel <- c("EVER", "CHOY", "MANA", "LHOT", "DHA1", "MAKA", "ANN1", "KANG", "LSHR")

# DEATHTYPE: 5 = crevasse, 6 = icefall collapse,
#            7 = avalanche, 8 = falling rock / ice
obj <- memb %>%
    filter(
        DEATH == TRUE,
        PEAKID %in% sel,
        as.numeric(MYEAR) >= 1950,
        DEATHTYPE %in% c(5, 6, 7, 8)
    ) %>%
    mutate(DEATHHGTM = as.numeric(DEATHHGTM))

# objective-hazard deaths per peak
obj %>%
    group_by(PEAKID) %>%
    summarise(n = n(), median_alt = median(DEATHHGTM, na.rm = TRUE)) %>%
    left_join(peaks %>% select(PEAKID, PKNAME), by = "PEAKID") %>%
    arrange(desc(n))

# single-day events with >= 3 deaths
obj %>%
    filter(!is.na(DEATHDATE)) %>%
    group_by(PEAKID, DEATHDATE) %>%
    summarise(n = n(), alt = median(DEATHHGTM, na.rm = TRUE), .groups = "drop") %>%
    filter(n >= 3) %>%
    left_join(peaks %>% select(PEAKID, PKNAME), by = "PEAKID") %>%
    arrange(desc(n))

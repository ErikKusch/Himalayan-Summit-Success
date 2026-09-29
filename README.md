# Shifting Risks in Himalayan Mountaineering Under a Changing Climate

This repository contains the data, R scripts, and exported results for the study of how climate change, commercialization, and increasing expedition efficiency are reshaping risk in high-altitude Himalayan mountaineering.

## Overview

The study combines two complementary sources:

- Daily ERA5-Land climate reanalysis fields for air temperature, snowfall, wind, and snow depth across the Himalaya from 1951 onward.
- Himalayan Database expedition records for Nepalese peaks, restricted to summits above 8,000 m and peaks with more than 30 recorded entries.

The analyses evaluate trends in summit-bid windows, expedition mortality, expedition size and crowding, seasonal climbing activity, and concurrent climate conditions during the pre-monsoon (March-May) and post-monsoon (September-November) seasons. Climate extremes are defined relative to season-specific 5th and 95th percentiles of the 1951-2021 baseline.

## Main findings

We demsontrate a shift toward a more fragile mountaineering system:

- Summit-bid windows have shortened by approximately 44% since the 1950s.
- Summit-level air temperatures have increased by approximately 0.7 degrees C.
- Mean wind speeds have increased, while snow depth has declined.
- High-wind days and low-snow conditions have become substantially more common.
- The pre-monsoon season is increasingly characterized by unstable wind and snowfall conditions associated with snowstorms and avalanche risk.
- The post-monsoon season is increasingly characterized by warm events, snowpack weakening, and greater icefall and rockfall exposure.
- Epedition ratio of Post- to Pre-monsoon expeditions increased from 50% prior to 1991 to 67% therafter.
- Per-expedition mortality has declined while multi-fatality events have become more frequent, indicating that individual safety improvements can coexist with increasing collective exposure.

Together, these results suggest that Himalayan mountaineering is becoming more efficient but also more temporally concentrated and collectively vulnerable.

## Repository structure

### Analysis scripts

R scripts ---------

Current analysis workflow

- `Himalayan Summit Success - DATA.r`
    Main data-preparation script. Loads expedition and peak data, loads or
    downloads ERA5-Land rasters, selects eligible peaks, extracts hourly and
    daily summit time series, calculates wind speed and seasonal percentile
    bounds, flags extreme conditions, calculates predictability metrics, and
    creates analysis-ready exports.

- `Mortality in Himalayan Mountaineering.r`
    Bayesian mortality analyses, including mortality trends, models controlling
    for expedition size, mortality by cause and altitude, and avalanche mortality
    summaries by season.

- `Mountaineers Crowding the Himalaya.r`
    Analyses of expedition counts, total mountaineer counts, expedition size,
    summit-bid windows, and change points in seasonal expedition activity.

- `Gradual Changes in the Himalayan Seasons.r`
    Peak-level and overall trend analyses of mean seasonal climate conditions and
    weather stability, represented by autocorrelation metrics.

- `Extreme Events in the Himalayan Seasons.r`
    Counts and run lengths of high and low seasonal extremes, with per-peak and
    mixed-model trend summaries.

- `Objective hazards in Himalayan mountaineering exacerbated by climate change.r`
    Seasonal analyses of storms, avalanches, deaths, avalanche-to-storm ratios,
    and avalanche-death altitude differences.

- `PlottingFunctions.r`
    Shared plotting functions, including Bayesian uncertainty ribbon plots used
    by the mortality and crowding analyses.

The file `X - PersonalSettings.R` is intentionally excluded from this deposit. It contains local settings and may contain climate-data-store credentials.

### Data

The `Data/` directory is not tracked through GitHub due to file size limitations. These files can be found at: https://doi.org/10.5281/zenodo.23042453

### Results and exports

The `Exports/` directory contains fitted model objects, analysis-ready model data, peak-level time series, predictability estimates, and summaries from the climate, mortality, and crowding analyses.

## Citation

TBD
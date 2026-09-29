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

- `Himalayan Summit Success - DATA.r`: data preparation and construction of the
	analysis datasets.
- `Mortality in Himalayan Mountaineering.r`: mortality and fatality analyses.
- `Mountaineers Crowding the Himalaya.r`: expedition size, participation, and
	crowding analyses.
- `Gradual Changes in the Himalayan Seasons.r`: long-term seasonal climate
	trends.
- `Extreme Events in the Himalayan Seasons.r`: seasonal extreme-event and
	predictability analyses.
- `Objective hazards in Himalayan mountaineering exacerbated by climate change.r`:
	climate-related objective-hazard analyses.
- `ANALYSIS_MortalityandObjectiveHazards.r`: combined mortality and objective
	hazard analysis.
- `PlottingFunctions.r`: shared plotting functions.
- `X - PersonalSettings.R`: local path and personal R settings.

### Data

The `Data/` directory is not tracked through GitHub due to file size limitations. These files can be found at: https://doi.org/10.5281/zenodo.23042453

### Results and exports

The `Exports/` directory contains fitted model objects, analysis-ready model data, peak-level time series, predictability estimates, and summaries from the climate, mortality, and crowding analyses.

## Citation

TBD
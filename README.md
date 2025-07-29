# NESP Building for the Future Project
Code repository for NESP 5.8 Building for the Future project

# Executive Summary 
This report proposes a dataset-agnostic methodology for generating location-specific future weather files using outputs from regional climate models tailored for building energy modelling in Australia. Traditional weather files used in energy simulations, typically based on historical meteorological observations, do not account for future climate change, limiting their effectiveness for resilience planning and climate-adapted design. To address this gap, we introduce a future-focused approach that integrates high-resolution reanalysis data with regional climate model projections, capturing both shifts in average climate and extremes for 2050. 

In this study we have applied the method to BARRA-2 reanalysis data to represent historical climate and the Australian Climate Service (ACS) CMIP6 future climate projections, covering both low (SSP1-2.6) and high (SSP3-7.0) emissions scenarios. However, the approach can be applied to other regional projections from the National Partnership for Climate Projections (NPCP) suite. Data selection prioritises model credibility, high spatial and temporal resolution, full variable coverage, and ensemble representation to capture climate uncertainty (see Section 2: Data Selection). 

The methodology consists of three key steps: 

- Variable Extraction: Relevant weather variables are extracted for each NatHERS climate zone location. These include daily variables for selecting representative months and hourly variables for constructing synthetic years. 
- Bias Adjustment: Given that climate model outputs contain biases, we apply the Quantile-Delta-Change (QDC) method to align future projections with the observed historical climate. This technique maintains the distributional characteristics of the modelled future climate while correcting for known biases. 
- Weather File Generation: Two types of synthetic weather files, Typical Meteorological Year (TMY) and Extreme Meteorological Year (XMY) are produced. TMY files represent average climate conditions. XMY files capture the most intense, the most severe and the longest projected heatwaves to carry out building performance resilience assessments. 
The resulting TMY and XMY files are formatted in the EnergyPlus Weather (EPW) format and include comprehensive metadata for use in building simulation platforms like EnergyPlus and NatHERS. Each location will have multiple files corresponding to different emissions scenarios, enabling scenario-based analysis of building performance under climate change.


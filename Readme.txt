# Prototype NatHERS TMY and XMY Weather Files from BARPA-R CMIP6 Simulations

## Disclaimer

These weather files are prototype datasets generated for testing purposes only.

They have been produced as part of a research and development activity and must not be used to inform building design, regulatory compliance, or energy efficiency assessments.

These files are not approved NatHERS weather files and are not suitable for operational use.

Distribution of these files is restricted to the designated testing group. They must not be shared outside this group without permission from the project team.

---

# Project

NESP Climate Systems Hub – Building for the Future

This work forms part of the Building for the Future project under the National Environmental Science Program (NESP).

The objective of this prototype dataset is to explore the use of regional climate model projections to generate future climate weather files compatible with the NatHERS simulation framework.

The prototype files support testing of:

* workflows for generating NatHERS-compatible weather files from climate model simulations
* bias-correction approaches suitable for building simulation variables
* potential future climate datasets for building performance modelling

---

# Input Climate Data

The prototype weather files are derived from simulations produced with the BARPA-R regional climate model, driven by CMIP6 global climate models.

## Regional climate model

Model: BARPA-R (Su et al. 2022)
Institution: Bureau of Meteorology
Domain: AUS-15 (~15 km resolution)
Temporal resolution: Hourly

BARPA-R dynamically downscales CMIP6 global climate simulations to produce high-resolution climate projections over Australia.

---

## Driving global climate models

The following CMIP6 models were used:

| Global Climate Model | Run         |
|----------------------|-------------|
| ACCESS-CM2           | r4i1p1f1    |
| ACCESS-ESM1-5        | r6i1p1f1    |
| CESM2                | r11i1p1f1   |
| CMCC-ESM2            | r1i1p1f1    |
| EC-Earth3            | r1i1p1f1    |
| MPI-ESM1-2-HR        | r1i1p1f1    |
| NorESM2-MM           | r1i1p1f1    |

The ensemble member identifiers follow the CMIP6 convention r*i*p*f*, representing the realisation (r), initialization method (i), physics configuration (p), and forcing dataset (f) used for the simulation.

---

## Emission scenarios

Two CMIP6 scenarios are included:

* SSP1-2.6 (ssp126) – low emissions scenario
* SSP3-7.0 (ssp370) – high emissions scenario

---

## Variables extracted from BARPA-R

The following hourly variables are extracted from BARPA-R simulations:

| Variable | Description                  |Units |
| -------- | ---------------------------- |------|
| tas      | Near-surface air temperature |K     |
| huss     | Specific humidity            |kg/kg |
| sfcWind  | 10 m wind speed              |m/s   |
| psl      | Mean sea level pressure      |kPa   |
| uas      | Zonal wind component         |m/s   |
| vas      | Meridional wind component    |m/s   |
| clt      | Total cloud cover            |%     |
| rsdsdir  | Direct shortwave radiation   |W m⁻² |
| rsdsdif  | Diffuse shortwave radiation  |W m⁻² |

### Selection of Representative BARPA-R Grid Cells

NatHERS station locations do not generally coincide exactly with the centres of the BARPA-R model grid (~15 km resolution). To obtain representative climate data for each site, station coordinates are mapped to the closest suitable BARPA-R grid cell using a set of physical selection criteria.

For each NatHERS location, the great-circle distance to all BARPA-R grid cells is calculated. Candidate grid cells are then filtered based on land fraction and elevation constraints to ensure that the selected grid cell is physically representative of the station location.

A grid cell is considered eligible if it satisfies the following conditions:

- land fraction ≥ 80 %
- surface elevation (orography) > 0 m
- grid-cell elevation within ±150 m of the NatHERS station elevation

Among all eligible grid cells, the grid cell with the minimum great-circle distance to the station location is selected.
The final station location is updated to the coordinates of the selected grid cell. The station elevation is kept to compute the correct station-level surface pressure (see next step).

---

# Pre-processing Adjustments

Prior to bias correction, several variables are transformed to ensure consistency with the NatHERS observational reference dataset. This step ensures that bias correction is performed on like-for-like variables between the climate model simulations and the NatHERS reference dataset.

### Surface pressure

BARPA-R provides mean sea level pressure (psl).
For NatHERS weather files, station-level surface pressure (ps) is required.

Surface pressure is derived from mean sea level pressure using the site elevation from the NatHERS dataset and standard atmospheric relationships. This adjustment ensures that pressure values represent the local station elevation rather than sea level conditions.

---

# Bias Correction Method

Selected variables are bias-corrected using Quantile Delta Change (QDC) scaling relative to a NatHERS observational reference dataset.

The following variables are bias-corrected:

* tas
* huss
* sfcWind
* psl
* rsds

The QDC approach preserves projected climate change signals while adjusting the distribution of model outputs to better match observed climatology (Irving & Macadam 2024). For radiation variables, bias-adjusted values at solar altitude angles below 10° are replaced with the corresponding observational values to avoid artefacts near sunrise and sunset.

Key characteristics of the implementation:

* applied at hourly resolution using a +/- 1 hour sliding temporal window for smooth transitioning between hours
* quantile mapping performed using 100 quantiles, allowing the bias correction to capture distributional differences between model and observations while maintaining robust quantile estimation given the available sample size (175320 hourly time steps)

---

# Variables Not Bias Corrected

The following variables are currently not bias corrected due to methodological limitations:

| Variable       | Reason                                                  |
| -------------- | ------------------------------------------------------- |
| clt            | NatHERS cloud cover is reported in octas (ordinal cate- |
|                | gorical data), which is not compatible with the contin- |
|                | uous quantile mapping approach used here.               |
| wind direction | Circular variable requiring specialised bias correction |

These variables are therefore retained directly from the NatHERS data.
This limitation should be considered when interpreting prototype results.

---

# Derived Variables

Additional variables required for NatHERS weather files are derived during processing.

## Shortwave radiation bias adjustment

To maintain physical consistency, the bias correction is applied only to the global horizontal irradiance (rsds). The direct and diffuse components are then reconstructed using the direct radiation fraction derived from the raw BARPA-R model output.

## Direct radiation fraction

The fraction of direct radiation is calculated from the unadjusted BARPA-R data as:

f_dir(month, hour) = mean(rsdsdir) / mean(rsds)

where the means are calculated for each combination of month and hour across the full time series.
This produces a monthly-hourly climatology of the direct radiation fraction.
For each timestep t, the corresponding fraction is assigned based on its month and hour.

## Reconstruction of radiation components

After bias correcting the global radiation (rsds_adj), the direct and diffuse components are reconstructed as:

rsdsdir_adj = rsds_adj * f_dir(time)
rsdsdif_adj = rsds_adj - rsdsdir_adj

This approach ensures:

* the physical constraint rsds = rsdsdir + rsdsdif
* stable direct/diffuse partitioning
* preservation of the model-simulated seasonal and diurnal structure of the direct fraction

## Direct Normal Irradiance (DNI)

The NatHERS reference dataset includes Direct Normal Irradiance (DNI), while BARPA-R provides direct horizontal irradiance (rsdsdir).
To ensure comparable variables for the NatHERS weather file format, DNI is derived using the solar altitude angle (slr_alt) from the NatHERS dataset:

DNI = rsdsdir / sin(slr_alt)

where:

* rsdsdir  = direct horizontal irradiance
* slr_alt  = solar altitude angle

For very low solar elevation angles, where the projection factor becomes unstable, DNI is set to zero to avoid unrealistic values.

---

# Unit Conversions for NatHERS Format

Variables are converted to units required by the NatHERS weather file format.

| Variable            | Conversion        |
| ------------------- | ----------------- |
| Temperature         | Kelvin → °C       |
| Specific humidity   | kg/kg → g/kg      |
| Radiation           | retained as W m⁻² |
| Wind speed          | retained as m/s   |
| Atmospheric pressure| kPa → Pa          |

---

# Typical Meteorological Year (TMY) and Extreme Meteorological Year (XMY) (!!!FOR SURENDRA TO COMPLETE!!!)

Prototype NatHERS weather files include TMY and XMY datasets constructed from the bias-corrected climate simulations.
The following methodology is currently under development and subject to refinement.

### Typical/Representative Meteorological Year (TMY/RMY)

* selection of representative months from the multi-year climate simulation
* evaluation of candidate months based on statistical similarity to long-term climate distributions
* assembly of representative months into a synthetic typical year

---

### Extreme Meteorological Year (XMY)

Extreme Meteorological Years (XMY) are constructed to represent climate conditions associated with severe heatwave events relevant for building performance.

Three XMY datasets are produced, each based on a different heatwave characteristic:

- the longest heatwave
- the most severe heatwave
- the most intense heatwave

These years are selected from the climate simulations based on heatwave metrics calculated from temperature data, allowing testing of building performance under a range of extreme heat conditions.

---

# Prototype File Structure

Each file represents:

* a single location (e.g. Melbourne)
* a single climate model (e.g. ACCESS-ESM1.5)
* a single emission scenario (e.g.ssp370)
* a future 20-year simulation period (e.g. 2041-2060)
* a weather file type identifier with method (e.g. Sandia RMY)

Example filename:

Melbourne_AUS-15_ACCESS-ESM1-5_ssp370_r6i1p1f1_BOM_BARPA-R_v1-r1_1hr_QDC-NatHERS_sandia_RMY_2041-2060

---

# Important Caveats

These prototype datasets have several limitations:

* Some variables are not bias corrected (cloud cover, wind direction)
* The TMY/XMY generation methodology is experimental
* Results depend on climate model uncertainty
* The files are not validated for building regulatory use

Further research and validation are required before such datasets could be considered for operational NatHERS applications.

---

# Project Contributors

This work was undertaken as part of the NESP Building for the Future project.
* Project lead: Ramona Dalla Pozza
* Contributors: Surendra Rauniyar, David Hoffmann

Contributors include researchers from:

* Bureau of Meteorology
* Australian Climate Service
* collaborating research partners within the NESP Climate Systems Hub

---

# Acknowledgements

Climate projections used in this dataset were produced using the BARPA-R regional climate modelling system developed by the Bureau of Meteorology. More details can be found in the Bureau Research Report #069.

# References

Su, C.-H., Stassen, C., Howard, E., Ye, H., Bell, S. S., Pepler, A., Dowdy, A. J., Tucker, S. O., Franklin, C. (2022), BARPA: New development of ACCESS-based regional climate modelling for Australian Climate Service, Bureau Research Report 069, accessed online http://www.bom.gov.au/research/publications/researchreports/BRR-069.pdf 

Irving, D. and Macadam, I. (2024) Application-Ready Climate Projections from CMIP6 using the Quantile Delta Change method. CSIRO Climate Innovation Hub Technical Note 5. https://doi.org/10.25919/03by-9y62 



































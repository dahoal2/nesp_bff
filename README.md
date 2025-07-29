# NESP Building for the Future Project
Code repository for NESP 5.8 Building for the Future project

---

## Table of Contents
- [Executive Summary](#executive-summary)
- [Explanation of scripts](#script-overview)
- [Step 1: Data Preparation](#step-1-data-extraction)
- [Step 2: QDC-scaling](#step-2-qdc-scaling)
- [Step 3: Calculate missing variables and convert units](#step-3-calculate-missing-vars)
- [Step 4: Create TMY and XMY files](#step-4-create-tmy-and-xmy)

---

## Executive Summary
This report proposes a dataset-agnostic methodology for generating location-specific future weather files using outputs from regional climate models tailored for building energy modelling in Australia. Traditional weather files used in energy simulations, typically based on historical meteorological observations, do not account for future climate change, limiting their effectiveness for resilience planning and climate-adapted design. To address this gap, we introduce a future-focused approach that integrates high-resolution reanalysis data with regional climate model projections, capturing both shifts in average climate and extremes for 2050. 

In this study we have applied the method to BARRA-2 reanalysis data to represent historical climate and the Australian Climate Service (ACS) CMIP6 future climate projections, covering both low (SSP1-2.6) and high (SSP3-7.0) emissions scenarios. However, the approach can be applied to other regional projections from the National Partnership for Climate Projections (NPCP) suite. Data selection prioritises model credibility, high spatial and temporal resolution, full variable coverage, and ensemble representation to capture climate uncertainty. 

The methodology consists of four key steps: 

- **Step 1** - Variable Extraction: Relevant weather variables are extracted for each NatHERS climate zone location. These include daily variables for selecting representative months and hourly variables for constructing synthetic years. 
- **Step 2** - Bias Adjustment: Given that climate model outputs contain biases, we apply the Quantile-Delta-Change (QDC) method to align future projections with the observed historical climate. This technique maintains the distributional characteristics of the modelled future climate while correcting for known biases.
- **Step 3**: Converting units and calculating missing variables required for NatHERS files that are directly available from the model outputs. Those are wet bulb temperature, dew point temperature, wind direction (compass index), and diffuse solar irradiance.
- **Step 4** - Weather File Generation: Two types of synthetic weather files, Typical Meteorological Year (TMY) and Extreme Meteorological Year (XMY) are produced. TMY files represent average climate conditions. XMY files capture the most intense, the most severe and the longest projected heatwaves to carry out building performance resilience assessments. 
The resulting TMY and XMY files are formatted in the EnergyPlus Weather (EPW) format and include comprehensive metadata for use in building simulation platforms like EnergyPlus and NatHERS. Each location will have multiple files corresponding to different emissions scenarios, enabling scenario-based analysis of building performance under climate change.

# Workflow
**Requirements**: Needs access to projects xp65 (hh5 module), eg3 (NESP project), hq89 (CCAM-v2203-SN (CSIRO)), py18 (BARPA-R), ob53 (BARRA-R2)
## Script overview
**Main scripts**
- `step1_extracting_variables.py`: Executable python script that pairs with `pbs_step1` in the `pbs_scripts` directory to extract needed variables for the 12 locations of interest.
- `step1_extracting_variables.ipynb`: Developer version of the above .py script. It can be run on ARE with the same outcome as the above script. Since the extraction of the variables can take a long time (~35 hours for CCAM daily data, 3-8 hours for BARPA-R daily/hourly data), it's best to sumbit the above script as a PBS job.
- `rechunk_ccam_parallel.sh`: Processes hourly CCAM-v2203-SN (CSIRO) data by rechunking the files from `time:1, lat:612, lon:929` to `time:2000, lat:15, lon:15` (emperically tested to be the fasted chunking, but only tested a few). Reads each annual file from `hq89` (CSIRO-CCAM directory), rechunks it with `ncks -O --cnk_dmn` and saves it to `/scratch/eg3/dh4185/rechunked/<GCM>/<scenario>/`. The script uses `GNU` with 10 parallel stream to run all 10 variables in parallel.
- `step2_qdc_scaling.ipynb`: Takes output files from Step 1 and applies QDC-scaling to them using BARRA-R2 data as reference in `g/data/eg3/nesp_bff/step1_raw_data_extraction/BARRA-R2/`. Adjust settings to RCM (`BARPA-R` or `CCAM-v2203-SN`) and timescale (`1hr` or `day`). Runs very quickly, and produces diagnostic plots (histograms, scaling factor and climatology) that are saved in `/g/data/eg3/nesp_bff/plots/qdc_adjustfactor/`.
- `step3_calc_missing_variables.ipynb`: Converts units to those required by the NatHERS files and calculates variables that are not available from RCM  output (wet bulb temperature, dew point temperature, wind direction (compass index), and diffuse solar irradiance).
- `step4_tmy_xmy`: Script to create TMY and XMY files. Currently empty.
- `utils.py`: Provides main functions for scripts above.

**Other scripts**
- `pbs_scripts/pbs_rechunk_ccam_parallel`: PBS job script to do the rechunking for hourly CCAM data. Uses ~6 hours with 256GB memory and 28 CPUs.
- `pbs_scripts/pbs_step1`: PBS job script to execute `step1_extracting_variables.py`.
- `pbs_scripts/pbs_step1_CCAM_day`: Same as above, but specified for *daily CCAM* where only one GCM is selected at a time (e.g. `--gcms=NorESM2-MM`), so that once can set up 7 PBS jobs in parallel for faster processing of such files. All daily CCAM data has been extracted so that this shouldn't be needed for now.
- `pbs_scripts/run_step1.ipynb`: runs `step1_extracting_variables.py` in an Jupyter notebook. Can be neglected.
- `plotting_checks.ipynb`: Jupyter notebook for plotting interim outputs for cross-checking.
- `future_heat.ipynb`: Analysis of changes in temperature threshold days (>=30,35,40) for 20-year time periods centred on 2050 and 2080 (ssp126, ssp370) relative to 2005 (historical). Produces outputs for AGCDv2, BARRA-R2 and bias-adjusted (QME) BARPA-R and CCAM-v2203-SN (CSIRO) (from project kj66). Output saved in `/home/565/dh4185/mn51-dh4185/repos_collab/nesp_bff/data/future_heat/` (note, not controlled by git).
  

## Step 1: Data extraction
First, the necessary variables (listed in table below) are extracted for each of the 12 point locations (Figure 1) from raw RCM data. BARPA-R is very efficiently chunked for our operation which favours little chunking across time and lots of chunking along latitude and longitude dimensions. CCAM is chunked for each time step but continuously along latitude and longitude dimensions which requires the dataset to be fully loaded, and is a s such very compute intense. It is currently not saved in a different format (comms with Marcus Thatcher, but they're dicussing options). This takes considerably more time to process all variables for ONE location and ONE GCM: BARPA-R day: ~2min, hourly: ~5min. CCAM daily: ~25min, hourly: >7.5h hours. Hence, CCAM hourly data is preprocessed (rechunked and stored on /scratch/eg3/dh4185/) to interim files per year, only loaded if all files for one GCM are present.

<table style="font-size: 12px;">
  <tr>
    <td>
      <b>Table 1: Daily and Hourly Variables</b>
      <table>
        <tr>
          <th>Daily Variable</th>
          <th>Hourly Variable</th>
        </tr>
        <tr><td>temperature_max (<code>tasmax</code>)</td><td>temperature (<code>tas</code>)</td></tr>
        <tr><td>temperature_min (<code>tasmin</code>)</td><td>cloud_cover (<code>clt</code>)</td></tr>
        <tr><td>humidity_specific_max (from hourly <code>huss</code>)</td><td>humidity_relative (<code>hurs</code>)</td></tr>
        <tr><td>humidity_specific_min (from hourly <code>huss</code>)</td><td>humidity_specific (<code>huss</code>)</td></tr>
        <tr><td>pressure (<code>psl</code>)</td><td>pressure (<code>psl</code>)</td></tr>
        <tr><td>wind_speed_10m (<code>sfcWind</code>)</td><td>wind_speed_10m (<code>sfcWind</code>)</td></tr>
        <tr><td>wind_speed_10m_max (<code>sfcWindmax</code>)</td><td>wind_direction_u (<code>uas</code>)</td></tr>
        <tr><td>solar_global (<code>rsds</code>)</td><td>wind_direction_v (<code>vas</code>)</td></tr>
        <tr><td>solar_direct (<code>rsdsdir</code>)</td><td>cloud_cover (<code>clt</code>)</td></tr>
        <tr><td>–</td><td>solar_global (<code>rsds</code>)</td></tr>
        <tr><td>–</td><td>solar_direct (<code>rsdsdir</code>)</td></tr>
      </table>
    </td>
    <td style="vertical-align: top; padding-left: 20px;">
      <b>Figure 1: 12 selected point locations with elevation.</b>
      <img src="locations_on_satmap.png" alt="Plot" width="700"/>
    </td>
  </tr>
</table>



**Key tasks**:
- When extracting daily and hourly data from BARPA-R or daily data from CSIRO-CCAM, run `qsub pbs_step1`. In the pbs script one can specify the RCM (BARPA-R *or* CCAM-v2203-SN), GCM(s) (if none specified, then all available for selected RCM are processed), scenario (historical, ssp370, ...), frequency (hourly or daily), start and end year (select 30-year periods, e.g. 2050 for the first offering), and number of dask workers depending on job size (use 20 when `ncpu=14` and `mem=64GB` in submitted PBS job). E.g.:
  ```bash
  python3 -W ignore ../step1_extracting_variables.py --rcm=CCAM-v2203-SN --scenario=ssp370 --freq=day --startYear=2035 --endYear=2064 --nworkers=20
  ```
  If no GCM is selected, all available for the selected RCM are processed. Output files are stored in `/g/data/eg3/nesp_bff/step1_raw_data_extraction/<RCM>/`. The file names are following the CORDEX naming convention, starting with the location name instead of the variable.
- When extracting hourly CSIRO-CCAM files, these need to be preprocessed first by rechunking the files from `time:1, lat:612, lon:929` to `time:2000, lat:15, lon:15`. The bash script `rechunk_ccam_parallel.sh` reads each file from `hq89` (CSIRO-CCAM directory), rechunks it with `ncks -O --cnk_dmn` and saves it to `/scratch/eg3/dh4185/rechunked/<GCM>/<scenario>/`. This is meant to be a temporary storage as hourly data of all variables for one GCM and 60 years (30 years histrorical and 30 years ssp370) takes up ~6.3TB of storage. Don't store more than two GCMs. Process the files on /scratch/ asap after completion with `qsub pbs_step1`! Running hourly CCAM data extraction with `qsub pbs_step1` without preprocessing yields a message (*"No files for GCM {_gcm} and {_var} exists. Run rechunk_ccam.sh first."*) and skips to the next gcm without computation. Run qsub pbs_rechunk_ccam_parallel
  It's best practice to run two PBS jobs simultaneously, one for processing historical data, one for processing future (ssp126 or ssp370) data.
- After rechunking and processing hourly CCAM data, remove interim files from `/scratch/` to free up space. Track `eg3` disk quota with
```bash
lquota | grep eg3
```


---

# Step 2: QDC-scaling
The bias adjustment process uses high-quality historical baseline data to correct climate model output. Instead of station observations, we use historical weather data from the BARRA-R2 reanalysis product. A multi-decade baseline, typically 30 years, is needed to capture regional climate variability. Currently, hourly climate projection data lack bias correction, which is essential for aligning them with observed climate statistics. To address this, we apply the Quantile-Delta-Change (QDC) method, a distribution-based statistical technique that adjusts each quantile of observed data using differences between historical and future model outputs. This approach preserves the observed distribution while incorporating the climate change signal. Unlike simpler methods that only shift the mean, QDC captures non-linear distributional changes and is better suited for applications needing realistic extremes (Irving & Macadam, 2024).


**Key tasks**:
- Run `process_data.py` for statistical analysis
- Apply filtering and aggregation methods
- Output is stored in `/outputs/intermediate/`

---

## Step 3: Visualization and Output

Describe how results are visualized and shared. Include screenshots or links if available.

**Key tasks**:
- Generate plots using `visualize_results.ipynb`
- Export maps and charts to `/outputs/figures/`
- Optional: run dashboard via Streamlit

---

## Step 4: Review and Validation

Explain how results are validated or interpreted, and any peer review or stakeholder feedback mechanisms.

**Key tasks**:
- Run validation checks in `validate_outputs.py`
- Compare outputs to baseline datasets
- Document insights and uncertainties

---

## License

[MIT](LICENSE) or another license type.

---

## Contact

For questions or collaboration requests, please contact [Your Name](mailto:your.email@example.com).

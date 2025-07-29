# NESP Building for the Future Project
Code repository for NESP 5.8 Building for the Future project

---

## Table of Contents
- [Executive Summary](#executive-summary)
- [Step 1: Data Preparation](#step-1-data-preparation)
- [Step 2: Processing and Analysis](#step-2-processing-and-analysis)
- [Step 3: Visualization and Output](#step-3-visualization-and-output)
- [Step 4: Review and Validation](#step-4-review-and-validation)

---

## Executive Summary
This report proposes a dataset-agnostic methodology for generating location-specific future weather files using outputs from regional climate models tailored for building energy modelling in Australia. Traditional weather files used in energy simulations, typically based on historical meteorological observations, do not account for future climate change, limiting their effectiveness for resilience planning and climate-adapted design. To address this gap, we introduce a future-focused approach that integrates high-resolution reanalysis data with regional climate model projections, capturing both shifts in average climate and extremes for 2050. 

In this study we have applied the method to BARRA-2 reanalysis data to represent historical climate and the Australian Climate Service (ACS) CMIP6 future climate projections, covering both low (SSP1-2.6) and high (SSP3-7.0) emissions scenarios. However, the approach can be applied to other regional projections from the National Partnership for Climate Projections (NPCP) suite. Data selection prioritises model credibility, high spatial and temporal resolution, full variable coverage, and ensemble representation to capture climate uncertainty (see Section 2: Data Selection). 

The methodology consists of three key steps: 

- Variable Extraction: Relevant weather variables are extracted for each NatHERS climate zone location. These include daily variables for selecting representative months and hourly variables for constructing synthetic years. 
- Bias Adjustment: Given that climate model outputs contain biases, we apply the Quantile-Delta-Change (QDC) method to align future projections with the observed historical climate. This technique maintains the distributional characteristics of the modelled future climate while correcting for known biases. 
- Weather File Generation: Two types of synthetic weather files, Typical Meteorological Year (TMY) and Extreme Meteorological Year (XMY) are produced. TMY files represent average climate conditions. XMY files capture the most intense, the most severe and the longest projected heatwaves to carry out building performance resilience assessments. 
The resulting TMY and XMY files are formatted in the EnergyPlus Weather (EPW) format and include comprehensive metadata for use in building simulation platforms like EnergyPlus and NatHERS. Each location will have multiple files corresponding to different emissions scenarios, enabling scenario-based analysis of building performance under climate change.

# Workflow
**Requirements**: Needs access to projects xp65 (hh5 module), eg3 (NESP project), hq89 (CSIRO-CCAM), py18 (BARPA-R), ob53 (BARRA-R2)

## Step 1: Data extraction
First, the necessary variables are extracted for each point location from raw RCM data. BARPA-R is very efficiently chunked four our operation which favours little chunking across time and lots of chunking along lat and lon. CCAM is chunked for each time step but not at all along lat and lon dimensions which requires the dataset to be fully loaded. It is currently not saved in a different format (comms with Marcus Thatcher). This takes considerable more time to process: BARPA-R day: ~2min, hourly: ~5min. CCAM daily: ~25min, hourly: >7.5h hours. Hence, CCAM hourly data is preprocessed (rechunked and stored on /scratch/eg3/dh4185/) to interim files per year, only loaded if all files for one GCM are present.

**Key tasks**:
- When extracting daily and hourly data from BARPA-R or daily data from CSIRO-CCAM, run `qsub pbs_step1`. In the pbs script one can specify the following options:
  `--rcm=CCAM-v2203-SN --gcms=NorESM2-MM --scenario=historical --freq=day --startYear=1985 --endYear=2014 --nworkers=10`. If no GCM is selected, all available for the selected RCM are processed. Output files are stored in `/g/data/eg3/nesp_bff/step1_raw_data_extraction/<RCM>/`. The file names are following the CORDEX naming convention, starting with the location name instead of the variable.
- When extracting hourly CSIRO-CCAM files, these need to be preprocessed first by rechunking the files from `time:1, lat:612, lon:929` to `time:2000, lat:15, lon:15` (emperically tested to be the fasted chunking, but only tested a few). The bash script `rechunk_ccam_parallel.sh` reads each file from `hq89` (CSIRO-CCAM directory), rechunks it with `ncks -O --cnk_dmn` and saves it to `/scratch/eg3/dh4185/rechunked/<GCM>/<scenario>/`. This is meant to be a temporary storage as hourly data for one GCM and 60 years (30 years histrorical and 30 years ssp370) takes up ~6.3TB of storage. Don't store more than two GCMs. Process the files on /scratch/ asap after completion with `qsub pbs_step1`!
-  without that step yields a message and nothing is computed.  
- Run `prepare_data.py` to clean and format the data

---

## Step 2: Processing and Analysis

Explain how data is processed and any key algorithms or scripts involved.

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

'''
######################################################################
## BUREAU OF METEOROLOGY
## National Environmental Science Program (NESP2)
## NESP CS5.8 - Building for the Future
##
## DATE:        June-2024
## SCRIPT:      utils.py
## AUTHOR:      david.hoffmann@bom.gov.au
##
## PURPOSE:     Core functions and dictionaries for scripts in this
                repository to derive TMY and XMY files from BARRA2 
                reanalysis and BARPA-R and CSIRO-CCAM CMIP6 climate data.
##
######################################################################
'''
import xarray as xr
import os
import sys
import glob
import numpy as np
from xclim import sdba

import warnings
import logging

warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')
logging.getLogger().setLevel(logging.CRITICAL)

####< Functions

def apply_hourly_qdc_sliding_window(da_obs, da_model_historical, da_model_future, var, kind='+', window=1):
    """
    Apply Quantile Delta Mapping for hourly data using a sliding window around each hour.
    
    Parameters:
    - da_obs: observed hourly data (e.g., reference)
    - da_model_historical: historical model hourly data
    - da_model_future: future model hourly data
    - kind: '+' or '*' for additive or multiplicative scaling
    - window: number of hours before/after to include (1 = 3hr window)
    
    Returns:
    - da_adjusted: xarray DataArray with adjusted values for all 24 hours
    - qdc_models: dict of trained QDC models keyed by hour
    """
    adjusted_chunks = []
    adjusted_slices = {}
    qdc_models = {}

    for hour in range(24):
        # Define sliding window: [hour - window, hour, hour + window]
        hours_to_use = [(hour + offset) % 24 for offset in range(-window, window + 1)]
        # print(hours_to_use)
        
        # Filter datasets for selected hours
        obs_slice = da_obs.where(da_obs.time.dt.hour.isin(hours_to_use), drop=True)
        hist_slice = da_model_historical.where(da_model_historical.time.dt.hour.isin(hours_to_use), drop=True)
        fut_slice = da_model_future.where(da_model_future.time.dt.hour.isin(hours_to_use), drop=True)

        # Train QDC
        QDC = sdba.QuantileDeltaMapping.train(
            fut_slice,
            hist_slice,
            nquantiles=100,
            group='time.month',
            kind=kind,
        )
        qdc_models[hour] = QDC
        
        # Adjust obs_slice using the QDC model
        adjusted = QDC.adjust(obs_slice, interp='linear')
        
        # Apply radiation-specific handling for multiplicative method
        if var in ['rsds', 'rsdsdir']:
            adjusted = radiation_0_adjustment(obs_slice, adjusted)
            
        # Clip variables that are bound by 0-100 (hurs and clt)
        if var in ["hurs","clt"]:
            adjusted = adjusted.clip(min=0, max=100) 

        # Extract only the target hour
        # For multiplicative adjustments (e.g., radiation), mask where adjusted is NaN but obs is not
        if kind == "*":
            ref = obs_slice.where(obs_slice.time.dt.hour == hour, drop=True)
            adjusted_hour = adjusted.where(adjusted.time.dt.hour == hour, drop=True)
            # Replace both NaN and inf (and -inf) with reference values
            invalid_mask = xr.ufuncs.isnan(adjusted_hour) | xr.ufuncs.isinf(adjusted_hour)
            adjusted_hour = xr.where(invalid_mask, ref, adjusted_hour)

        else:
            adjusted_hour = adjusted.where(adjusted.time.dt.hour == hour, drop=True)

        # adjusted_hour = adjusted.where(adjusted.time.dt.hour == hour, drop=True)
        adjusted_chunks.append(adjusted_hour)
        adjusted_slices[hour] = adjusted_hour

    # Concatenate all hourly chunks and sort
    da_adjusted = xr.concat(adjusted_chunks, dim='time').sortby('time')
    
    return da_adjusted, qdc_models, adjusted_slices

def radiation_0_adjustment(reference, adjusted, threshold=1.0):
    """
    Prevent NaNs when applying multiplicative QDC to radiation data.
    Keeps night-time values (where reference ≈ 0) unchanged.
    Radiation threshold for night is set to 1 Wm-2.
    """
    is_day = reference > threshold
    rad_adjusted = xr.where(is_day, adjusted, reference)
    return rad_adjusted

def align_future_to_historical(ds_model_hist, ds_model_future):
    """
    Aligns the future dataset's time dimension to match the reference dataset,
    dropping extra days if necessary (e.g., leap days in a different calendar).

    Parameters
    ----------
    ds_model_hist : xr.Dataset
        Reference dataset with desired time axis.
    ds_model_future : xr.Dataset
        Future dataset to be aligned.

    Returns
    -------
    ds_model_future_aligned : xr.Dataset
        Future dataset with time dimension aligned to ds_model_hist.
    """
    len_hist = ds_model_hist.time.size
    len_future = ds_model_future.time.size
    # print(f"Time length: Ref = {len_hist}, model = {len_future}")

    if len_future > len_hist:
        # Drop the last time steps to match length
        ds_model_future = ds_model_future.isel(time=slice(0, len_hist))
    elif len_future < len_hist:
        raise ValueError(f"Future dataset is too short (length {len_future}), but reference has {len_hist}.")

    # Assign the reference time axis
    ds_model_future = ds_model_future.assign_coords(time=ds_model_hist.time)

    return ds_model_future

def process_humidity(da,_var):
    print(f"Resample hourly data for {_var}.")
    if _var == 'humidity_specific_max':
        da = da.resample(time='1D').max().rename('hussmax')
    elif _var == 'humidity_specific_min':
        da = da.resample(time='1D').min().rename('hussmin')
    
    return da

def process_time(da,_var,_timescale):
    if _timescale == "day":
        da['time'] = da['time'].dt.floor('D')
        da = da.sel(time=~da.get_index("time").duplicated())
    elif _timescale == "1hr" and _var in ['clt','rsds','rsdsdir']:
        da['time'] = da['time'].dt.floor('h')
        da = da.sel(time=~da.get_index("time").duplicated())
    return da
    
def plot_qdc_hourly_diagnostics(da_obs_1hr,da_model_hist_1hr,da_model_fut_1hr,
                                QDC_dict,adjusted_slices,var,model,location):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    import matplotlib.gridspec as gridspec
    
    selected_hours = [0,6,12,18] #[0,2,4,6,8,10,12,14,16,18,20,22]
    
    for hour in selected_hours:
        qdc = QDC_dict[hour]
        adj = adjusted_slices[hour]
    
        obs_hour = da_obs_1hr.where(da_obs_1hr.time.dt.hour == hour, drop=True)
        hist_hour = da_model_hist_1hr.where(da_model_hist_1hr.time.dt.hour == hour, drop=True)
        fut_hour = da_model_fut_1hr.where(da_model_fut_1hr.time.dt.hour == hour, drop=True)
    
        # Set up GridSpec
        fig = plt.figure(figsize=(12, 8))
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.2])
    
        ax0 = fig.add_subplot(gs[0, 0])  # Monthly means
        ax1 = fig.add_subplot(gs[0, 1])  # Heatmap
        ax2 = fig.add_subplot(gs[1, :])  # Full-width histogram
    
        fig.suptitle(f"{location}: QDC Diagnostics for {model} {var} - Hour {hour:02d} UTC", fontsize=14, y=0.95)
    
        # 1a. Monthly means
        obs_hour.groupby("time.month").mean().plot(ax=ax0, label="'Obs' hist. (2000)", color="blue")
        hist_hour.groupby("time.month").mean().plot(ax=ax0, label="Model hist. (2000)", color="green")
        fut_hour.groupby("time.month").mean().plot(ax=ax0, label="Model future (2050-ssp370)", color="red")
        adj.groupby("time.month").mean().plot(ax=ax0, label="'Obs' adjusted", linestyle="--", color="purple")
        ax0.set_title("Monthly Mean by Hour")
        ax0.legend()
    
        # 1b. Adjustment factor heatmap
        af = qdc.ds.af
        vmin, vmax = float(af.min()), float(af.max())
        if vmin < 0 and vmax > 0:
            cmap = "RdBu_r"
            center = 0.0
        else:
            cmap = "viridis"
            center = None
    
        af.transpose("month", "quantiles").plot.pcolormesh(
            ax=ax1, cmap=cmap, center=center,
            x="quantiles", y="month",
            cbar_kwargs={"label": "Adjustment Factor"}
        )
        ax1.set_title("QDC Adjustment Factors")
    
        # 2. Histogram
        sns.histplot(obs_hour.values.flatten(), ax=ax2, label="'Obs' hist. (2000)", kde=True, stat="density", bins=50, color="blue")
        sns.histplot(hist_hour.values.flatten(), ax=ax2, label="Model hist. (2000)", kde=True, stat="density", bins=50, color="green")
        sns.histplot(fut_hour.values.flatten(), ax=ax2, label="Model future (2050-ssp370)", kde=True, stat="density", bins=50, color="red")
        sns.histplot(adj.values.flatten(), ax=ax2, label="'Obs' adjusted", kde=True, stat="density", bins=50, color="purple")
        ax2.set_title("Distributions")
        ax2.legend()
    
        # plt.tight_layout()
        # plt.savefig(fig_file, bbox_inches='tight')


# Function to compute great-circle distance using NumPy (vectorised)
def great_circle_distance(lat1, lon1, lat2, lon2):
    R = 6371.0  # Earth's radius in km
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    return R * c  # Distance in km

def update_locations(da_sftlf,dict_locations):
    ### Find nearest grid cell on land to lat/lon coordinates from stations
    # Extract grid information
    latitudes = da_sftlf['lat'].values
    longitudes = da_sftlf['lon'].values
    valid_mask = da_sftlf.values >= 90  # Mask where terrain fraction is at least 90
    
    # Get valid lat/lon pairs
    valid_lats, valid_lons = np.meshgrid(latitudes, longitudes, indexing='ij')  # Ensure correct shape
    valid_lats, valid_lons = valid_lats[valid_mask], valid_lons[valid_mask]  # Apply mask
    
    # Find the closest valid grid cell for each location & store in a new dictionary
    updated_locations = {}
    for name, info in dict_locations.items():
        lat, lon, elev = info['Lat'], info['Lon'], info['Elev']
        distances = great_circle_distance(lat, lon, valid_lats, valid_lons)
        
        min_idx = np.argmin(distances)  # Index of closest point
        new_lat, new_lon = valid_lats[min_idx], valid_lons[min_idx]
        
        updated_locations[name] = {'Lat': new_lat, 'Lon': new_lon, 'Elev': elev}  # Preserve elevation
        
    return updated_locations
             
def preprocess_location(ds,lat,lon):
    return ds.sel(lat=lat,lon=lon, method="nearest")
             # .sel(time=slice(f'{start_year}-01-01', f'{end_year}-12-31'))

def get_git_hash():
    """Returns the git hash for the working repository"""
    git_repo = git.Repo(sys.argv[0], search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    git_hash = git.Repo(git_root).heads[0].commit
    git_text = " (Git hash: %s)" %(str(git_hash)[0:7])
    return git_text
    
###< Metadata

locations = {
    'Darwin': {'Lat': -12.42, 'Lon': 130.89, 'Elev': 30},
    'Cairns': {'Lat': -16.87, 'Lon': 145.75, 'Elev': 2},
    'Brisbane': {'Lat': -27.48, 'Lon': 153.04, 'Elev': 8},
    'Longreach': {'Lat': -23.44, 'Lon': 144.28, 'Elev': 192},
    'Mildura': {'Lat': -34.23, 'Lon': 142.09, 'Elev': 50},
    'Adelaide': {'Lat': -34.95, 'Lon': 138.51, 'Elev': 2},
    'Perth': {'Lat': -31.92, 'Lon': 115.87, 'Elev': 25},
    'Sydney': {'Lat': -33.95, 'Lon': 151.17, 'Elev': 39},
    'Melbourne': {'Lat': -37.67, 'Lon': 144.83, 'Elev': 31},
    'Canberra': {'Lat': -35.31, 'Lon': 149.19, 'Elev': 577},
    'Hobart': {'Lat': -42.89, 'Lon': 147.33, 'Elev': 51},
    'Thredbo': {'Lat': -36.50, 'Lon': 148.29, 'Elev': 1380}
}

model_dict = {
    "BARRA-R2":{"root_dir":"/g/data/ob53/BARRA2/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/{}/",
                 "sftlf":"/g/data/ob53/BARRA2/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/fx/sftlf/latest/sftlf_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1.nc",
                 "grid": "AUS-11",
                 "org": "BOM",
                 "gcms":{
                     "ERA5":{"mdl_run":"hres","version": "v1","created":"latest"}}},
    "BARPA-R":{"root_dir":"/g/data/py18/BARPA/output/CMIP6/DD/AUS-15/BOM/{}/{}/{}/{}/{}/{}/",
               "sftlf":"/g/data/py18/BARPA/output/CMIP6/DD/AUS-15/BOM/ACCESS-ESM1-5/historical/r6i1p1f1/BARPA-R/v1-r1/fx/sftlf/latest/sftlf_AUS-15_ACCESS-ESM1-5_historical_r6i1p1f1_BOM_BARPA-R_v1-r1_fx.nc",
               "grid": "AUS-15",
               "org": "BOM",
               "gcms":{
                   "ACCESS-ESM1-5":{"mdl_run":"r6i1p1f1","version": "v1-r1","created": "latest"},
                   "ACCESS-CM2":{"mdl_run":"r4i1p1f1","version": "v1-r1","created": "latest"},
                   "CESM2":{"mdl_run":"r11i1p1f1","version": "v1-r1","created": "latest"},
                   "CMCC-ESM2":{"mdl_run":"r1i1p1f1","version": "v1-r1","created": "latest"},
                   "EC-Earth3":{"mdl_run":"r1i1p1f1","version": "v1-r1","created": "latest"},
                   "MPI-ESM1-2-HR":{"mdl_run":"r1i1p1f1","version": "v1-r1","created": "latest"},
                   "NorESM2-MM":{"mdl_run":"r1i1p1f1","version": "v1-r1","created": "latest"}}},
    "CCAM-v2203-SN":{"root_dir":"/g/data/hq89/CCAM/output/CMIP6/DD/AUS-10i/CSIRO/{}/{}/{}/{}/{}/{}/",
               "sftlf":"/g/data/hq89/CCAM/output/CMIP6/DD/AUS-10i/CSIRO/ACCESS-ESM1-5/historical/r6i1p1f1/CCAM-v2203-SN/v1-r1/fx/sftlf/v20231206/sftlf_AUS-10i_ACCESS-ESM1-5_historical_r6i1p1f1_CSIRO_CCAM-v2203-SN_v1-r1.nc",
               "grid": "AUS-10i",
               "org": "CSIRO",
               "gcms":{
                   "ACCESS-ESM1-5":{"mdl_run":"r6i1p1f1","version": "v1-r1","created": "v20240327"},
                   "ACCESS-CM2":{"mdl_run":"r4i1p1f1","version": "v1-r1","created": "v20231206"},
                   "CESM2":{"mdl_run":"r11i1p1f1","version": "v1-r1","created": "v20231206"},
                   "CMCC-ESM2":{"mdl_run":"r1i1p1f1","version": "v1-r1","created": "v20231206"},
                   "CNRM-ESM2-1":{"mdl_run":"r1i1p1f2","version": "v1-r1","created": "v20231206"},
                   "EC-Earth3":{"mdl_run":"r1i1p1f1","version": "v1-r1","created": "v20231206"},
                   "NorESM2-MM":{"mdl_run":"r1i1p1f1","version": "v1-r1","created": "v20231206"}}}
    }

vars_1hr = {
    'temperature': ['tas'],
    'humidity_relative': ['hurs'],
    'humidity_specific': ['huss'],
    'wind_speed_10m': ['sfcWind'],
    'pressure': ['psl'],
    'wind_direction_u': ['uas'],
    'wind_direction_v': ['vas'],
    'cloud_cover': ['clt'],
    'solar_global': ['rsds'],
    'solar_direct': ['rsdsdir']
}

vars_day = {
    'temperature_max': ['tasmax'],
    'temperature_min': ['tasmin'],
    'humidity_specific_max': ['huss'],
    'humidity_specific_min': ['huss'],
    'pressure': ['psl'],
    'wind_speed_10m': ['sfcWind'],
    'wind_speed_10m_max': ['sfcWindmax'],
    'solar_global': ['rsds'],
    'solar_direct': ['rsdsdir']
}


cmap_dict = {
    "tasmax": "OrRd",
    "tasmin": "OrRd",
    "hussmax": "PuBuGn",
    "hussmin": "PuBuGn",
    "psl": "cividis",
    "sfcWind": "YlGnBu",
    "sfcWindmax": "YlGnBu",
    "rsds": "YlOrBr",
    "rsdsdir": "YlOrRd"
}

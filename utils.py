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
from scipy.spatial import QhullError

import warnings
import logging

warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')
logging.getLogger().setLevel(logging.CRITICAL)

####< Functions
import numpy as np
import xarray as xr

def _proj_from_altitude(solar_altitude_deg: xr.DataArray) -> xr.DataArray:
    """
    Projection factor for converting between direct horizontal and DNI.
    cos(zenith) = sin(altitude)
    """
    alt_rad = np.deg2rad(solar_altitude_deg)
    return np.sin(alt_rad)


def direct_horizontal_to_dni(
    direct_h: xr.DataArray,
    solar_altitude_deg: xr.DataArray,
    *,
    min_proj: float = 0.05,     # ~3° solar elevation threshold
    max_dni: float = 1400.0,
    out_name: str = "DNI",
) -> xr.DataArray:
    """
    Convert direct horizontal irradiance (W/m²) to Direct Normal Irradiance (DNI, W/m²).
    DNI = direct_h / sin(altitude)

    Safeguards:
      - sun below horizon / very low elevation (proj <= min_proj) -> 0
      - clip DNI to [0, max_dni]
    """
    proj = _proj_from_altitude(solar_altitude_deg)

    dni = xr.where(
        proj > min_proj,
        direct_h / proj,
        0.0
    )

    dni = dni.clip(min=0.0, max=max_dni)
    dni.attrs.update({"units": "W m-2", "description": "Direct normal irradiance"})
    return dni.rename(out_name)


def dni_to_direct_horizontal(
    dni: xr.DataArray,
    solar_altitude_deg: xr.DataArray,
    *,
    min_proj: float = 0.05,      # ~3° solar elevation threshold
    max_direct_h: float = 1400.0,
    out_name: str = "rsdsdir",
) -> xr.DataArray:
    """
    Convert Direct Normal Irradiance (DNI, W/m²) to direct horizontal irradiance (W/m²).
    direct_h = DNI * sin(altitude)

    Safeguards:
      - sun below horizon / very low elevation (proj <= min_proj) -> 0
      - clip direct_h to [0, max_direct_h]
    """
    proj = _proj_from_altitude(solar_altitude_deg)

    direct_h = xr.where(
        proj > min_proj,
        dni * proj,
        0.0
    )

    direct_h = direct_h.clip(min=0.0, max=max_direct_h)
    direct_h.attrs.update({"units": "W m-2", "description": "Direct horizontal irradiance"})
    return direct_h.rename(out_name)

def fix_duplicate_times(ds):
    """Remove duplicate time entries, keeping the first occurrence."""
    if not ds.get_index("time").is_unique:
        print(f"Found {ds.get_index('time').duplicated().sum()} duplicate time entries.")
        ds = ds.sel(time=~ds.get_index("time").duplicated())
        print("Duplicates removed.")
    else:
        print("No duplicates found.")
    return ds
    
def convert_sea_level_pressure_to_station_pressure(mslp_kpa, h_m):
    """Return the atmospheric station pressure in Pa given sea level pressure in hPa*10 and station elevation in m.
    https://www.weather.gov/media/epz/wxcalc/stationPressure.pdf
    https://www.weather.gov/media/epz/wxcalc/pressureConversion.pdf
    From Surendra's script
    """

    # convert (or keep) pressure and elevation inputs as floats
    # Pa = float(Pa)
    # h_m = float(h_m)

    # # convert from kPa to inHg
    Pa_inHg = mslp_kpa * 0.29529983071445

    # # calculate station pressure according to formula from https://www.weather.gov/epz/wxcalc_stationpressure
    Pstn_inHg = Pa_inHg * ((288 - 0.0065*h_m)/288)**5.2561

    # # convert from inHg to Pa
    Pstn_kpa = 3.3863886666666714 * Pstn_inHg

    # Pstn_kpa = mslp_kpa *((288 - 0.0065*h_m)/288)**5.2561
    Pstn_kpa.attrs = mslp_kpa.attrs.copy()
    return Pstn_kpa

import numpy as np
import xarray as xr
from xclim import sdba

def radiation_0_adjustment_solarfix(reference: xr.DataArray,
                           adjusted: xr.DataArray,
                           *,
                           threshold: float = 1.0) -> xr.DataArray:
    """
    Prevent NaNs/inf when applying multiplicative QDC to radiation data by keeping
    night-time values (where reference ~ 0) unchanged.
    """
    is_day = reference > threshold
    return xr.where(is_day, adjusted, reference)

def apply_hourly_qdc_sliding_window_solarfix(
    da_obs: xr.DataArray,
    da_model_historical: xr.DataArray,
    da_model_future: xr.DataArray,
    var: str,
    *,
    nq: int = 100,
    kind: str = "+",
    window: int = 1,
    solar_alt: xr.DataArray | None = None,
    min_alt: float = 10.0,
    rad_night_threshold: float = 1.0,
    clip_bounds: dict | None = None,
    # what to use when masking low-sun / spikes / invalids for radiation
    # - "obs": use observational reference (old behaviour)
    # - "future_model": use raw BARPA-R future (recommended for future-signal fidelity)
    # - "historical_model": use raw BARPA-R historical
    radiation_mask_source: str = "future_model",
):
    """
    Hourly QDC with sliding window, plus radiation-specific post-adjustment fixes.

    Radiation handling:
      - If solar_alt is provided, values with solar altitude < min_alt are replaced with the
        chosen replacement source (default: raw future model).
      - Physically unrealistic spikes (>1300 W/m²) are replaced with the same replacement source.
      - For multiplicative kind ("*"), any NaN/Inf in the final hourly product are replaced; for
        radiation variables this replacement is also controlled by `radiation_mask_source`
        (default: raw future model), otherwise it falls back to the observational reference.
    """
    if kind not in {"+", "*"}:
        raise ValueError("kind must be '+' or '*'")

    is_radiation = var in {"rsds", "rsdsdir", "rsdsdif"}
    clip_bounds = clip_bounds or {}

    if radiation_mask_source not in {"obs", "future_model", "historical_model"}:
        raise ValueError("radiation_mask_source must be one of: 'obs', 'future_model', 'historical_model'")

    adjusted_chunks = []
    adjusted_slices = {}
    qdc_models = {}

    for hour in range(24):
        hours_to_use = [(hour + offset) % 24 for offset in range(-window, window + 1)]

        obs_slice = da_obs.where(da_obs.time.dt.hour.isin(hours_to_use), drop=True)
        hist_slice = da_model_historical.where(da_model_historical.time.dt.hour.isin(hours_to_use), drop=True)
        fut_slice = da_model_future.where(da_model_future.time.dt.hour.isin(hours_to_use), drop=True)

        # Train QDC
        QDC = sdba.QuantileDeltaMapping.train(
            fut_slice,
            hist_slice,
            nquantiles=nq,
            group="time.month",
            kind=kind,
        )
        qdc_models[hour] = QDC

        # Adjust - try linear interpolation first and fall back to 'nearest' if it fails for radiation
        try:
            adjusted = QDC.adjust(obs_slice, interp="linear")
        except (QhullError, ValueError):
            if is_radiation:
                print(f"⚠ Fallback to nearest for {var}, hour={hour}")
                adjusted = QDC.adjust(obs_slice, interp="nearest")
            else:
                raise

        # Radiation-specific post-adjustment handling
        if is_radiation:
            # Keep night-time values unchanged (multiplicative kind)
            if kind == "*":
                adjusted = radiation_0_adjustment(
                    obs_slice,
                    adjusted,
                    threshold=rad_night_threshold,
                )

            # Choose what we replace masked radiation with (obs vs raw model)
            if radiation_mask_source == "obs":
                repl = obs_slice.sel(time=adjusted.time)
            elif radiation_mask_source == "historical_model":
                repl = hist_slice.sel(time=adjusted.time)
            else:  # "future_model"
                repl = fut_slice.sel(time=adjusted.time)

            if solar_alt is not None:
                alt_slice = solar_alt.sel(time=adjusted.time)

                # 1) Low-solar-altitude handling
                adjusted = xr.where(
                    alt_slice < min_alt,
                    repl,
                    adjusted,
                )

            # 2) Physically unrealistic spikes
            adjusted = xr.where(
                adjusted > 1300,
                repl,
                adjusted,
            )

        # Clip bounded variables
        if var in {"hurs", "clt"}:
            adjusted = adjusted.clip(min=0, max=100)

        if var in clip_bounds:
            lo, hi = clip_bounds[var]
            adjusted = adjusted.clip(min=lo, max=hi)

        # Extract only the target hour
        adjusted_hour = adjusted.where(adjusted.time.dt.hour == hour, drop=True)

        # Replace invalid values (multiplicative kind)
        if kind == "*":
            invalid = xr.ufuncs.isnan(adjusted_hour) | xr.ufuncs.isinf(adjusted_hour)

            # Default replacement for invalids is observational reference at the target hour,
            # but for radiation variables we keep future-signal fidelity by using `repl`.
            if is_radiation:
                # repl is defined only in the radiation block above; recompute robustly here
                if radiation_mask_source == "obs":
                    ref_hour = obs_slice.where(obs_slice.time.dt.hour == hour, drop=True)
                elif radiation_mask_source == "historical_model":
                    ref_hour = hist_slice.where(hist_slice.time.dt.hour == hour, drop=True)
                else:  # "future_model"
                    ref_hour = fut_slice.where(fut_slice.time.dt.hour == hour, drop=True)
            else:
                ref_hour = obs_slice.where(obs_slice.time.dt.hour == hour, drop=True)

            adjusted_hour = xr.where(invalid, ref_hour, adjusted_hour)

        adjusted_chunks.append(adjusted_hour)
        adjusted_slices[hour] = adjusted_hour

    da_adjusted = xr.concat(adjusted_chunks, dim="time").sortby("time")
    return da_adjusted, qdc_models, adjusted_slices


def apply_qdc_sliding_window(da_obs, da_model_historical, da_model_future, var, nq=100, kind='+', window=1, nslots=None):
    """
    Apply Quantile Delta Mapping (QDC) with a sliding window across the diurnal cycle.
    Works with both hourly (24 slots/day) and half-hourly (48 slots/day) data.

    Parameters
    ----------
    da_obs : xarray.DataArray
        Observed data with time coordinate (hourly or half-hourly).
    da_model_historical : xarray.DataArray
        Historical model data (same resolution as da_obs).
    da_model_future : xarray.DataArray
        Future model data (same resolution as da_obs).
    var : str
        Variable name (needed for special cases: radiation, humidity, etc.).
    kind : str, optional
        '+' or '*' for additive or multiplicative scaling. Default='+'
    window : int, optional
        Sliding window size in slots (1 = ±1 slot). Default=1
    nslots : int or None, optional
        Force number of slots (24 or 48). If None, detect from data.
        - Hourly data → always 24
        - Half-hourly data → default 48, but you can override with 24.

    Returns
    -------
    da_adjusted : xarray.DataArray
        Bias-adjusted data.
    qdc_models : dict
        Dictionary of trained QDC models, keyed by slot index.
    adjusted_slices : dict
        Dictionary of adjusted slices, keyed by slot index.
    """

    # --- Detect whether data is hourly or half-hourly ---
    minutes = da_obs.time.dt.minute.values

    if (minutes % 60 == 0).all():
        # hourly → force 24 slots
        detected = 24
        slot_index = da_obs.time.dt.hour
        print("Detected hourly data → using 24 slots.")
        if nslots not in (None, 24):
            raise ValueError("Hourly data can only be used with 24 slots.")

    elif (minutes % 30 == 0).all():
        # half-hourly
        detected = 48
        if nslots is None:
            nslots = 48  # default for half-hourly
        elif nslots not in (24, 48):
            raise ValueError("Half-hourly data must use 24 or 48 slots.")
        if nslots == 48:
            slot_index = da_obs.time.dt.hour * 2 + da_obs.time.dt.minute // 30
            print("Detected half-hourly data → using 48 slots.")
        else:
            slot_index = da_obs.time.dt.hour
            print("Detected half-hourly data → aggregating into 24 slots.")

    else:
        raise ValueError("Data must be hourly or half-hourly aligned.")

    adjusted_chunks = []
    adjusted_slices = {}
    qdc_models = {}

    for slot in range(nslots):
        # Define sliding window across slots
        slots_to_use = [(slot + offset) % nslots for offset in range(-window, window + 1)]

        # Filter datasets for selected slots
        obs_slice  = da_obs.where(slot_index.isin(slots_to_use), drop=True)
        hist_slice = da_model_historical.where(slot_index.isin(slots_to_use), drop=True)
        fut_slice  = da_model_future.where(slot_index.isin(slots_to_use), drop=True)

        # Train QDC
        QDC = sdba.QuantileDeltaMapping.train(
            fut_slice,
            hist_slice,
            nquantiles=nq,
            group='time.month',
            kind=kind,
        )
        qdc_models[slot] = QDC

        # Adjust observations
        from scipy.spatial import QhullError
        import warnings
        
        try:
            adjusted = QDC.adjust(obs_slice, interp='linear')
        except (ValueError, QhullError) as e:
            # fallback to nearest — more robust for low number of points / degenerate data
            warnings.warn(f"Linear QDM interpolation failed (fallback to 'nearest') for var {var if 'var' in locals() else ''}: {e}")
            adjusted = QDC.adjust(obs_slice, interp='nearest')
        # adjusted = QDC.adjust(obs_slice, interp='linear')

        # Radiation-specific handling for multiplicative case
        if var in ['rsds', 'rsdsdir']:
            adjusted = radiation_0_adjustment(obs_slice, adjusted)

        # Clip bounded variables
        if var in ["hurs", "clt"]:
            adjusted = adjusted.clip(min=0, max=100)

        # Extract only the target slot
        if nslots == 24:
            target_mask = (adjusted.time.dt.hour == slot)
        else:  # 48-slot case
            target_hour = slot // 2
            target_minute = (slot % 2) * 30
            target_mask = (adjusted.time.dt.hour == target_hour) & (adjusted.time.dt.minute == target_minute)

        adjusted_slot = adjusted.where(target_mask, drop=True)

        # Handle NaNs/infs for multiplicative method
        if kind == "*":
            ref = obs_slice.where(target_mask, drop=True)
            invalid_mask = xr.ufuncs.isnan(adjusted_slot) | xr.ufuncs.isinf(adjusted_slot)
            adjusted_slot = xr.where(invalid_mask, ref, adjusted_slot)

        adjusted_chunks.append(adjusted_slot)
        adjusted_slices[slot] = adjusted_slot

    # Concatenate all slots and sort by time
    da_adjusted = xr.concat(adjusted_chunks, dim='time').sortby('time')

    return da_adjusted, qdc_models, adjusted_slices


def apply_hourly_qdc_sliding_window(da_obs, da_model_historical, da_model_future, var, nq=100, kind='+', window=1):
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
            nquantiles=nq,
            group='time.month',
            kind=kind,
        )
        qdc_models[hour] = QDC
        
        # Adjust obs_slice using the QDC model
        adjusted = QDC.adjust(obs_slice, interp='linear')
        
        # Apply radiation-specific handling for multiplicative method
        if var in ['rsds', 'rsdsdir','rsdsdif']:
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
    elif _timescale == "1hr" and _var in ['clt','rsds','rsdsdir','rsdsdif','pr']:
        da['time'] = da['time'].dt.floor('h')
        da = da.sel(time=~da.get_index("time").duplicated())
    return da
    
def plot_qdc_hourly_diagnostics(da_obs_1hr,da_model_hist_1hr,da_model_fut_1hr,
                                QDC_dict,adjusted_slices,var,model,location,hour,ssp,fut_period):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    import matplotlib.gridspec as gridspec
    
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
    
    fig.suptitle(f"{location}: QDC Diagnostics for {model} {var} - Hour {hour:02d} UTC", fontsize=14, y=0.97)
    
    # 1a. Monthly means
    obs_hour.groupby("time.month").mean().plot(ax=ax0, label="'NatHERS' hist. (1995-2014)", color="blue")
    hist_hour.groupby("time.month").mean().plot(ax=ax0, label=f"{model} hist. (1995-2014)", color="green")
    fut_hour.groupby("time.month").mean().plot(ax=ax0, label=f"{model} future ({fut_period}-{ssp})", color="red")
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
    sns.histplot(obs_hour.values.flatten(), ax=ax2, label="'NatHERS' hist. (1995-2014)", kde=True, stat="density", bins=50, color="blue")
    sns.histplot(hist_hour.values.flatten(), ax=ax2, label=f"{model} hist. (1995-2014)", kde=True, stat="density", bins=50, color="green")
    sns.histplot(fut_hour.values.flatten(), ax=ax2, label=f"{model} future ({fut_period}-{ssp})", kde=True, stat="density", bins=50, color="red")
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

import numpy as np

def update_locations_orog_sftlf(da_sftlf, da_orog, dict_locations, *,
                     land_fraction_threshold=80,
                     orog_threshold=150,
                     return_gridcell_elev=True):
    """
    Update station locations to the closest eligible grid cell, matching the
    plotting routine criteria.

    Eligibility:
      - land fraction >= land_fraction_threshold
      - orog > 0
      - orog within [max(0, elev - orog_threshold), elev + orog_threshold]
        (symmetric tolerance; adjust if you only want upper bound)

    Parameters
    ----------
    da_sftlf : xr.DataArray (lat, lon)
    da_orog  : xr.DataArray (lat, lon)
    dict_locations : dict
        {name: {'Lat':..., 'Lon':..., 'Elev':...}, ...}
    land_fraction_threshold : float
    orog_threshold : float (m)
    return_gridcell_elev : bool
        If True, update Elev to the gridcell orog. If False, preserve station Elev.

    Returns
    -------
    updated_locations : dict
    """
    # Basic coordinate sanity
    valid_lats = da_orog.lat
    valid_lons = da_orog.lon

    updated_locations = {}

    for name, info in dict_locations.items():
        lat0, lon0, elev0 = info["Lat"], info["Lon"], info["Elev"]

        # Distance field (lat, lon)
        distances = great_circle_distance(lat0, lon0, valid_lats, valid_lons)

        # Elevation eligibility window (symmetric; lower bound clipped to 0)
        lower = max(0.0, float(elev0) - float(orog_threshold))
        upper = float(elev0) + float(orog_threshold)

        elig = (
            (da_orog > 0)
            & (da_sftlf >= land_fraction_threshold)
            & (da_orog >= lower)
            & (da_orog <= upper)
        )

        distances = distances.where(elig)

        # If no eligible cells, fall back to land-only (or closest anywhere)
        if distances.isnull().all():
            # fallback 1: land fraction only (and orog > 0)
            elig2 = (da_orog > 0) & (da_sftlf >= land_fraction_threshold)
            distances2 = utils.great_circle_distance(lat0, lon0, valid_lats, valid_lons).where(elig2)

            if distances2.isnull().all():
                # fallback 2: closest anywhere on grid
                distances2 = utils.great_circle_distance(lat0, lon0, valid_lats, valid_lons)

            idx = distances2.argmin(dim=("lat", "lon"))
            new_lat = distances2.lat[idx["lat"]].item()
            new_lon = distances2.lon[idx["lon"]].item()
        else:
            idx = distances.argmin(dim=("lat", "lon"))
            new_lat = distances.lat[idx["lat"]].item()
            new_lon = distances.lon[idx["lon"]].item()

        # Choose elevation to store
        if return_gridcell_elev:
            new_elev = float(da_orog.sel(lat=new_lat, lon=new_lon).item())
        else:
            new_elev = elev0

        updated_locations[name] = {"Lat": new_lat, "Lon": new_lon, "Elev": np.round(new_elev,1)}

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
    'Melbourne': {'Lat': -37.666, 'Lon': 144.832, 'Elev': 118.8},
    'Canberra': {'Lat': -35.305, 'Lon': 149.201, 'Elev': 580.0},
    'Darwin': {'Lat': -12.424, 'Lon': 130.893, 'Elev': 35.0},
    'Cairns': {'Lat': -16.874, 'Lon': 145.746, 'Elev': 8.3},
    'Brisbane': {'Lat': -27.392, 'Lon': 153.129, 'Elev': 9.5},
    'Longreach': {'Lat': -23.437, 'Lon': 144.277, 'Elev': 192.5},
    'Mildura': {'Lat': -34.236, 'Lon': 142.087, 'Elev': 51.1},
    'Adelaide': {'Lat': -34.921, 'Lon': 138.622, 'Elev': 51.0},
    'Perth': {'Lat': -31.927, 'Lon': 115.976, 'Elev': 20.0},
    'Sydney': {'Lat': -33.941, 'Lon': 151.173, 'Elev': 5.0},
    'Hobart': {'Lat': -42.89, 'Lon': 147.328, 'Elev': 51.4}#,
    # 'Thredbo': {'Lat': -36.49, 'Lon': 148.30, 'Elev': 1380.0}
}

model_dict = {
    "BARRA-R2":{"root_dir":"/g/data/ob53/BARRA2/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/{}/",
                 "sftlf":"/g/data/ob53/BARRA2/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/fx/sftlf/latest/sftlf_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1.nc",
                 "orog":"/g/data/ob53/BARRA2/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/fx/orog/latest/orog_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1.nc",
                 "grid": "AUS-11",
                 "org": "BOM",
                 "gcms":{
                     "ERA5":{"mdl_run":"hres","version": "v1","created":"latest"}}},
    "BARRA-C2":{"root_dir":"/g/data/ob53/BARRA2/output/reanalysis/AUST-04/BOM/ERA5/historical/hres/BARRA-C2/v1/{}/",
                 "sftlf":"/g/data/ob53/BARRA2/output/reanalysis/AUST-04/BOM/ERA5/historical/hres/BARRA-C2/v1/fx/sftlf/latest/sftlf_AUST-04_ERA5_historical_hres_BOM_BARRA-C2_v1.nc",
                 "orog":"/g/data/ob53/BARRA2/output/reanalysis/AUST-04/BOM/ERA5/historical/hres/BARRA-C2/v1/fx/orog/latest/orog_AUST-04_ERA5_historical_hres_BOM_BARRA-C2_v1.nc",
                 "grid": "AUST-04",
                 "org": "BOM",
                 "gcms":{
                     "ERA5":{"mdl_run":"hres","version": "v1","created":"latest"}}},
    "BARPA-R":{"root_dir":"/g/data/py18/BARPA/output/CMIP6/DD/AUS-15/BOM/{}/{}/{}/{}/{}/{}/",
               "sftlf":"/g/data/py18/BARPA/output/CMIP6/DD/AUS-15/BOM/ACCESS-ESM1-5/historical/r6i1p1f1/BARPA-R/v1-r1/fx/sftlf/latest/sftlf_AUS-15_ACCESS-ESM1-5_historical_r6i1p1f1_BOM_BARPA-R_v1-r1_fx.nc",
               "orog":"/g/data/py18/BARPA/output/CMIP6/DD/AUS-15/BOM/ACCESS-ESM1-5/historical/r6i1p1f1/BARPA-R/v1-r1/fx/orog/latest/orog_AUS-15_ACCESS-ESM1-5_historical_r6i1p1f1_BOM_BARPA-R_v1-r1_fx.nc",
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
    "BARPA-C":{"root_dir":"/g/data/py18/BARPA/output/CMIP6/DD/AUST-04/BOM/{}/{}/{}/{}/{}/{}/",
               "sftlf":"/g/data/py18/BARPA/output/CMIP6/DD/AUST-04/BOM/ACCESS-ESM1-5/historical/r6i1p1f1/BARPA-C/v1-r1/fx/sftlf/latest/sftlf_AUST-04_ACCESS-ESM1-5_historical_r6i1p1f1_BOM_BARPA-C_v1-r1_fx.nc",
               "orog":"/g/data/py18/BARPA/output/CMIP6/DD/AUST-04/BOM/ACCESS-ESM1-5/historical/r6i1p1f1/BARPA-C/v1-r1/fx/orog/latest/orog_AUST-04_ACCESS-ESM1-5_historical_r6i1p1f1_BOM_BARPA-C_v1-r1_fx.nc",
               "grid": "AUST-04",
               "org": "BOM",
               "gcms":{
                   "ACCESS-ESM1-5":{"mdl_run":"r6i1p1f1","version": "v1-r1","created": "latest"},
                   "EC-Earth3":{"mdl_run":"r1i1p1f1","version": "v1-r1","created": "latest"}}},
    "CCAM-v2203-SN":{"root_dir":"/g/data/hq89/CCAM/output/CMIP6/DD/AUS-10i/CSIRO/{}/{}/{}/{}/{}/{}/",
               "sftlf":"/g/data/hq89/CCAM/output/CMIP6/DD/AUS-10i/CSIRO/ACCESS-ESM1-5/historical/r6i1p1f1/CCAM-v2203-SN/v1-r1/fx/sftlf/v20231206/sftlf_AUS-10i_ACCESS-ESM1-5_historical_r6i1p1f1_CSIRO_CCAM-v2203-SN_v1-r1.nc",
                "orog":"/g/data/hq89/CCAM/output/CMIP6/DD/AUS-10i/CSIRO/ACCESS-ESM1-5/historical/r6i1p1f1/CCAM-v2203-SN/v1-r1/fx/orog/v20231206/orog_AUS-10i_ACCESS-ESM1-5_historical_r6i1p1f1_CSIRO_CCAM-v2203-SN_v1-r1.nc",
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
    # 'solar_global': ['rsds'],
    'solar_direct': ['rsdsdir'],
    'solar_diffuse': ['rsdsdif']
    # 'precipitation': ['pr']
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

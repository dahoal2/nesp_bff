#!/usr/bin/env python3
import argparse
import os
import glob
import sys
import numpy as np
import xarray as xr

# Local atmos checkout
sys.path.append("/home/565/dh4185/mn51-dh4185/repos_collab/atmos/")
from atmos.thermo import wet_bulb_temperature

def convert_temp_unit(da: xr.DataArray) -> xr.DataArray:
    units = da.attrs.get("units", None)
    if units not in ("degC", "K"):
        raise ValueError(f"Unexpected temperature units {units!r} (expected 'K' or 'degC')")
    if units == "K":
        da = da - 273.16
        da.attrs["units"] = "degC"
    return da


def cloud_frac_to_oktas(clt_pct: xr.DataArray) -> xr.DataArray:
    clt = clt_pct.clip(0, 100)
    bins = np.array([-0.1, 1, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 100.1])

    okta = xr.apply_ufunc(
        np.digitize,
        clt,
        kwargs={"bins": bins},
        vectorize=True,
        output_dtypes=[np.int16],
    ) - 1

    # okta = okta.where(clt.notnull())
    okta.attrs["units"] = "oktas"
    okta.attrs["description"] = "cloud oktas 0-8"
    return okta.rename("total_cloud_cover")
# def cloud_frac_to_oktas(clt_frac):
#     clt_frac = clt_frac.clip(0, 100)
#     bins = [-0.1,1,12.5,25,37.5,50,62.5,75,87.5,100.1]
#     # oktas = np.arange(0,9) # 0 to 8
#     okta = xr.apply_ufunc(
#         lambda x: np.digitize(x, bins) -1,
#         clt_frac,
#         input_core_dims=[[]],
#         output_core_dims=[[]],
#         vectorize=True,
#         output_dtypes=[int]
#     ).rename('clt')
#     okta_ = okta.where(okta.notnull())
#     okta_.attrs["units"] = "oktas"
#     okta_.attrs["description"] = "cloud oktas 0-8"
#     return okta#.rename("clt")
    
def _proj_from_altitude(solar_altitude_deg: xr.DataArray) -> xr.DataArray:
    """
    Projection factor for converting between direct horizontal and DNI.
    cos(zenith) = sin(altitude)
    """
    alt_rad = np.deg2rad(solar_altitude_deg)
    return np.sin(alt_rad)


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

def wind_dir_degree_from_vector(uas: xr.DataArray, vas: xr.DataArray) -> xr.DataArray:
    from xclim.indicators import convert
    _, wind_dir_degrees = convert.wind_speed_from_vector(uas=uas, vas=vas)
    return wind_dir_degrees.rename("windDir_degrees")


def degree_to_compassIndex(degrees: xr.DataArray, sfcWind: xr.DataArray, calm_threshold=0.026) -> xr.DataArray:
    deg = degrees.where(np.isfinite(degrees))
    idx = (np.floor((deg - 11.25) / 22.5) % 16) + 1
    idx = xr.where((sfcWind > calm_threshold) & idx.notnull(), idx, 0)
    idx = idx.astype("int16")
    idx.attrs["units"] = "index"
    idx.attrs["description"] = "0-16, 0=calm, 1=NNE, 16=N"
    idx.attrs["long_name"] = "Near-Surface Wind Direction"
    idx.attrs["standard_name"] = "wind_direction"
    return idx.rename("wind_dir")#("Wind_direction")


def process_one(infile: str, outfile: str, ds_nathers: xr.Dataset, *, compress: bool = True) -> None:
    ds = xr.open_dataset(infile)

    # Temperature
    tas = convert_temp_unit(ds["tas"]).rename("tas")#("Dry_bulb_temperature")
    tas.attrs["cell_methods"] = "time: point (interval: 1H)"
    tas.attrs["long_name"] = "Near-Surface Air Temperature"
    tas.attrs["standard_name"] = "air_temperature" 

    # Wet bulb temperature
    twbt = wet_bulb_temperature(
        ds["psl"], ds["tas"], ds["huss"],
        saturation="isobaric", phase="liquid", polynomial=True
    )
    twbt.attrs["units"] = ds["tas"].attrs["units"] 
    # print(twbt)
    twbt = convert_temp_unit(twbt).rename("twbt")#("Wet_bulb_temperature")
    twbt.attrs["cell_methods"] = "time: point (interval: 1H)"
    twbt.attrs["long_name"] = "Near-Surface Wet Bulb Temperature"
    twbt.attrs["standard_name"] = "wet_bulb_temperature"
    
    # Moisture (note: this is specific humidity in g/kg, naming kept as you used)
    huss = (ds["huss"] * 1000).rename("huss")#("Absolute_moisture_content")
    huss.attrs["units"] = "g/kg"
    huss.attrs["cell_methods"] = "time: point (interval: 1H)"
    huss.attrs["long_name"] = "Near-Surface Specific Humidity"
    huss.attrs["standard_name"] = "specific_humidity"

    # Pressure
    psl = (ds["psl"] / 1000).rename("psl")#("Atmospheric_pressure")
    psl.attrs["units"] = "kPa"
    psl.attrs["cell_methods"] = "time: point (interval: 1H)"
    psl.attrs["long_name"] = "Surface Air Pressure"
    psl.attrs["standard_name"] = "surface_air_pressure"

    # Wind
    sfcWind = ds["sfcWind"].rename("wind_speed")#("Wind_speed")
    sfcWind.attrs["units"] = "m s-1"
    sfcWind.attrs["cell_methods"] = "time: point (interval: 1H) area: interpolation (method: bilinear)"
    sfcWind.attrs["long_name"] = "Near-Surface Wind Speed"
    sfcWind.attrs["standard_name"] = "wind_speed"
    wind_dir = wind_dir_degree_from_vector(ds["uas"], ds["vas"])
    wind_direction = degree_to_compassIndex(wind_dir, sfcWind, calm_threshold=0.026)
    wind_direction.attrs["cell_methods"] = "time: point (interval: 1H) area: interpolation (method: bilinear)"

    # Radiation 
    # rsdsdir is in fact DNI
    dni = ds['rsdsdir'].rename("rsdsdir")#('Direct_solar_irradiance')
    dni.attrs["cell_methods"] = "time: mean (interval: 1H)"
    dni.attrs["units"] = "W m-2"
    dni.attrs["long_name"] = "Surface Direct Downwelling Shortwave Radiation"
    dni.attrs["standard_name"] = "surface_direct_downwelling_shortwave_flux_in_air"
    
    rsdsdif = ds['rsdsdif'].rename("rsdsdif")#('Diffuse_solar_irradiance')
    rsdsdif.attrs["cell_methods"] = "time: mean (interval: 1H)"
    rsdsdif.attrs["units"] = "W m-2"
    rsdsdif.attrs["long_name"] = "Surface Diffuse Downwelling Shortwave Radiation"
    rsdsdif.attrs["standard_name"] = "surface_diffuse_downwelling_shortwave_flux_in_air"

    # Convert DNI to direct horizontal radiation to derive rsds
    proj = np.sin(np.deg2rad(ds_nathers["slr_alt"]))
    min_proj = np.sin(np.deg2rad(5))   # ~0.087
                        
    rsdsdir = xr.where(
        proj > min_proj,
        dni * proj,
        0.0
    ).rename("rsdsdir")
    # rsdsdir = dni_to_direct_horizontal(dni, proj_alt, out_name="rsdsdir")

    rsds = (rsdsdir+rsdsdif).rename("rsds")#('Global_solar_irradiance')
    rsds.attrs["units"] = "W m-2"
    rsds.attrs["cell_methods"] = "time: mean (interval: 1H)"
    rsds.attrs["long_name"] = "Surface Downwelling Shortwave Radiation"
    rsds.attrs["standard_name"] = "surface_downwelling_shortwave_flux_in_air"

    # Cloud
    clt_okta = cloud_frac_to_oktas(ds["clt"])
    clt_okta.attrs["cell_methods"] = "time: mean (interval: 1 hour)"
    clt_okta.attrs["long_name"] = "Total Cloud Cover"
    clt_okta.attrs["standard_name"] = "cloud_octas"

    ds_new = xr.merge(
        [tas, twbt, huss, psl, sfcWind, wind_direction, rsds, dni, rsdsdif, clt_okta],
        compat="no_conflicts",
    )

    ds_new_, ds_nathers_ = xr.align(ds_new, ds_nathers, join="inner")

    vars_nathers = ds_nathers_[["flag-tas","flag-huss","flag-psl","flag-wind_speed","flag-total_cloud_cover","flag-wind_dir","slr_alt","slr_azm",
                "flag-rsds","flag-rsdsdif","flag-rsdsdir"]]
    
    ds_out = xr.merge([ds_new_,vars_nathers],compat='override')
    ds_out = ds_out.assign_coords(round_method=ds_nathers.round_method)


    # Remove global attrs + drop known extras
    ds_out.attrs = {}
    ds_out = ds_out.drop_vars(["crs"], errors="ignore")

    # Kick out 29 Februarys if it exists
    ds_out_no_leap = ds_out.convert_calendar("noleap") 
    
    # mask_feb29 = (ds_out.time.dt.month == 2) & (ds_out.time.dt.day == 29)    
    # if mask_feb29.any():
    #     ds_out_no_leap = ds_out.sel(time=~mask_feb29)
    # else:
    #     ds_out_no_leap = ds_out
        
    # Write float vars as float32 (+ optional compression)
    encoding = {}
    for v in ds_out_no_leap.data_vars:
        if np.issubdtype(ds_out_no_leap[v].dtype, np.floating):
            encoding[v] = {"dtype": "float32"}
            if compress:
                encoding[v].update({"zlib": True, "complevel": 4})

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    ds_out_no_leap.to_netcdf(outfile, encoding=encoding)


def main():
    ap = argparse.ArgumentParser(description="Step3: calc missing NatHERS vars from step2 QDC files (per location).")
    ap.add_argument("--loc", required=True, help="Location prefix, e.g. Cairns")
    ap.add_argument("--input_dir", default="/g/data/eg3/nesp_bff/step2_qdc_scaling/BARPA-R/")
    ap.add_argument("--output_dir", default="/g/data/eg3/nesp_bff/step3_calc_missing_vars/")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--no_compress", action="store_true")
    args = ap.parse_args()

    nathers_dir = "/g/data/eg3/spr548/projects/nesp_bff/data/raw_data/NatHERS/historical/"
    ds_nathers = xr.open_dataset(glob.glob(f"{nathers_dir}{args.loc}_*_NatHERS_UTC_19891231-20151231.nc")[0]).sel(time=slice("1995","2014")).load()

    pattern = os.path.join(args.input_dir, f"{args.loc}_*1hr_*_QDC-*.nc")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files found for {args.loc} using pattern: {pattern}")

    print(f"{args.loc}: found {len(files)} files")

    for infile in files:
        outfile = os.path.join(
            args.output_dir,
            os.path.basename(infile).replace(".nc", "_NatHERSvars.nc")
        )

        if (not args.overwrite) and os.path.exists(outfile):
            print(f"Skip (exists): {os.path.basename(outfile)}")
            continue

        print(f"Process: {os.path.basename(infile)}")
        process_one(infile, outfile, ds_nathers, compress=(not args.no_compress))

    print("Done.")


if __name__ == "__main__":
    main()
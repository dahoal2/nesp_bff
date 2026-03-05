#!/usr/bin/env python3
import argparse
import os
import glob
import time
import numpy as np
import pandas as pd
import xarray as xr

# Make sure utils import works
import sys
sys.path.append("/home/565/dh4185/mn51-dh4185/repos_collab/nesp_bff/")
import utils
from utils import locations, model_dict


def to_pandas_time(ds: xr.Dataset) -> xr.Dataset:
    """Convert time coord to pandas datetime64[ns] so xarray .dt works reliably."""
    time_pd = pd.to_datetime(ds.time.values.astype(str))
    return ds.assign_coords(time=time_pd)


def build_paths(root_dir, rcm, loc, gcm, scenario, fut_period, starty, endy, hist_slice, obs_name):
    in_dir = f"{root_dir}step1_raw_data_extraction/{rcm}/"
    out_dir = f"{root_dir}step2_qdc_scaling/{rcm}/"

    out_file = (
        f"{out_dir}{loc}_"
        f"{model_dict[rcm]['grid']}_"
        f"{gcm}_{scenario}_"
        f"{model_dict[rcm]['gcms'][gcm]['mdl_run']}_"
        f"{model_dict[rcm]['org']}_"
        f"{rcm}_{model_dict[rcm]['gcms'][gcm]['version']}_"
        f"1hr_{starty+1}-{endy}_QDC-{obs_name}.nc"
    )

    model_historical_file = (
        f"{in_dir}{loc}_"
        f"{model_dict[rcm]['grid']}_"
        f"{gcm}_historical_"
        f"{model_dict[rcm]['gcms'][gcm]['mdl_run']}_"
        f"{model_dict[rcm]['org']}_"
        f"{rcm}_{model_dict[rcm]['gcms'][gcm]['version']}_"
        f"1hr_{hist_slice}.nc"
    )

    model_future_file = (
        f"{in_dir}{loc}_"
        f"{model_dict[rcm]['grid']}_"
        f"{gcm}_{scenario}_"
        f"{model_dict[rcm]['gcms'][gcm]['mdl_run']}_"
        f"{model_dict[rcm]['org']}_"
        f"{rcm}_{model_dict[rcm]['gcms'][gcm]['version']}_"
        f"1hr_{starty}-{endy}.nc"
    )

    return out_file, model_historical_file, model_future_file


def main():
    p = argparse.ArgumentParser(description="Step2 QDC scaling for one (loc, scenario, fut_period).")
    p.add_argument("--loc", required=True, help="Location key (must exist in utils.locations)")
    p.add_argument("--scenario", required=True, choices=["ssp126", "ssp370"], help="Scenario")
    # p.add_argument("--fut_period", required=True, choices=["2030", "2050", "2070"], help="Future period label")
    p.add_argument("--rcm", default="BARPA-R", help="RCM name (default: BARPA-R)")
    p.add_argument("--obs_name", default="NatHERS", help="Observation dataset name tag")
    p.add_argument("--nslots", type=int, default=24, help="Number of slots for sliding-window QDC")
    p.add_argument("--window", type=int, default=1, help="Window size for sliding-window QDC")
    p.add_argument("--nq", type=int, default=100, help="Number of quantiles for QDC")
    p.add_argument("--root_dir", default="/g/data/eg3/nesp_bff/", help="Root directory")
    p.add_argument("--nathers_dir", default="/g/data/eg3/spr548/projects/nesp_bff/data/raw_data/NatHERS/historical/", help="NatHERS directory")
    p.add_argument("--hist_slice", default="1994-2014", help="Historical file slice label (used in filename)")
    p.add_argument("--hist_years", default="1995-2014", help="Actual historical years to select (e.g. 1995-2014)")
    p.add_argument("--do_plots", action="store_true", help="Generate diagnostics plots (slow)")
    p.add_argument("--plot_dir", default="/g/data/eg3/nesp_bff/plots/qdc_adjustfactor/", help="Plot output dir")
    args = p.parse_args()

    fut_dict = {
        "2030": {"start": 2020, "end": 2040},
        "2050": {"start": 2040, "end": 2060},
        "2070": {"start": 2060, "end": 2080},
    }

    vars_qdc = ["tas", "huss", "sfcWind", "psl", "rsds"]  # "rsdsdir", "rsdsdif"
    plot_vars = vars_qdc #["tas", "huss", "sfcWind", "psl", "rsds", "rsdsdir", "rsdsdif"]
    vars_1hr = ["tas", "huss", "sfcWind", "psl", "rsds", "rsdsdir", "rsdsdif"] # "uas", "vas", "clt" -> cloud and wind dir from NatHERS in step3 

    loc = args.loc
    if loc not in locations:
        raise ValueError(f"Unknown loc '{loc}'. Must be one of: {list(locations.keys())}")
    
    if loc == "Thredbo":
        print("Skipping Thredbo (as per notebook).")
        return
    
    # Load NatHERS once per location (outside fut_period loop)
    # root_nathers = "/g/data/eg3/spr548/projects/nesp_bff/data/raw_data/NatHERS/historical/"
    nathers_glob = glob.glob(f"{args.nathers_dir}{loc}_*NatHERS_UTC_19891231-20151231.nc")
    if not nathers_glob:
        raise FileNotFoundError(f"NatHERS file not found for loc={loc} in {args.nathers_dir}/")
    ds_nathers = xr.open_dataset(nathers_glob[0]).sel(time=slice(*args.hist_years.split("-"))).load()
    ds_nathers = to_pandas_time(ds_nathers)
    
    for fut_period in ["2030", "2050", "2070"]:
        starty = fut_dict[fut_period]["start"]
        endy = fut_dict[fut_period]["end"]
    
        print(f"=== Running: loc={loc} scenario={args.scenario} fut_period={fut_period} rcm={args.rcm} ===")
    
        for gcm in model_dict[args.rcm]["gcms"]:
            t0 = time.time()
            print(f"\n--- GCM: {gcm} ---")
    
            out_file, hist_file, fut_file = build_paths(
                args.root_dir, args.rcm, loc, gcm, args.scenario, fut_period,
                starty, endy, args.hist_slice, args.obs_name
            )
    
            if os.path.exists(out_file):
                print(f"Output exists, skipping: {out_file}")
                continue
    
            print(f"Opening historical: {hist_file}")
            ds_hist = xr.open_dataset(hist_file, use_cftime=True).sel(
                time=slice(args.hist_years.split("-")[0], args.hist_years.split("-")[1])
            ).load()
            ds_hist = to_pandas_time(ds_hist)
    
            print(f"Opening future: {fut_file}")
            ds_fut = xr.open_dataset(fut_file, use_cftime=True).sel(
                time=slice(str(starty + 1), str(endy))
            ).load()
            ds_fut = to_pandas_time(ds_fut)
    
            ds_nathers_loc = ds_nathers.copy()
            ds_nathers_loc["lat"] = ds_hist["lat"]
            ds_nathers_loc["lon"] = ds_hist["lon"]
    
            ds_fut_aligned = utils.align_future_to_historical(ds_hist, ds_fut)
            
            # Keep a copy of the raw future for fraction calculation (before any QDC changes)
            ds_fut_raw = ds_fut_aligned  # this is raw model future, aligned in time
            f_dir_time = utils.compute_fdir_time(ds_fut_raw)

            def _drop_feb29(da: xr.DataArray) -> xr.DataArray:
                return da.where(~((da.time.dt.month == 2) & (da.time.dt.day == 29)), drop=True)
            
            def _has_feb29_time(time):
                return bool(((time.dt.month == 2) & (time.dt.day == 29)).any())
            
            model_has_feb29 = _has_feb29_time(ds_fut_aligned["time"]) or _has_feb29_time(ds_hist["time"])

            if not model_has_feb29:
                ds_nathers_loc = _drop_feb29(ds_nathers_loc)
    
            var_datasets = []

            for v in vars_1hr:
                print(f"Var: {v}")
            
                # Skip rsdsdir/rsdsdif here; we'll build them after rsds is adjusted
                if v in {"rsdsdir", "rsdsdif"}:
                    continue
            
                if v in vars_qdc:
                    print(f"  QDC: {v}")
                    kind = "+" if v == "tas" else "*"
            
                    if v == "psl":
                        print("  Converting sea-level pressure to station pressure...")
                        ds_hist = ds_hist.assign(
                            psl=utils.convert_sea_level_pressure_to_station_pressure(ds_hist.psl, locations[loc]["Elev"])
                        )
                        ds_fut_aligned = ds_fut_aligned.assign(
                            psl=utils.convert_sea_level_pressure_to_station_pressure(ds_fut_aligned.psl, locations[loc]["Elev"])
                        )
            
                    da_obs = ds_nathers_loc[v]
                    da_hist = ds_hist[v]
                    da_fut = ds_fut_aligned[v]

                    if v == "rsds":
                        da_hist.attrs["units"] = da_obs.attrs.get("units", da_hist.attrs.get("units", ""))
                        da_fut.attrs["units"] = da_obs.attrs.get("units", da_fut.attrs.get("units", ""))
                    
                 #=======================================================

                    # slr_alt=None

                    # #============== convert direct radiation ==================
                    # if v == "rsdsdir":
                    #     # Outname 'rsdsdir' even though it's technically 'DNI' but has to match var
                    #     # names in obs
                    #     slr_alt = ds_nathers_loc["slr_alt"]          # degrees
                    #     da_hist_al, slr_alt_al = xr.align(ds_hist[v], slr_alt, join="left")
                    #     da_fut_al, slr_alt_al = xr.align(ds_fut_aligned[v], slr_alt, join="left")

                    #     min_proj = np.sin(np.deg2rad(8))   # ~0.087
                    #     proj = np.sin(np.deg2rad(slr_alt_al))
                        
                    #     da_hist = xr.where(
                    #         proj > min_proj,
                    #         da_hist_al / proj,
                    #         0.0
                    #     ).rename("rsdsdir")
                    #     da_hist.attrs.update({"units": "W m-2", "description": "Direct normal irradiance"})

                    #     da_fut = xr.where(
                    #         proj > min_proj,
                    #         da_fut_al / proj,
                    #         0.0
                    #     ).rename("rsdsdir")
                    #     da_fut.attrs.update({"units": "W m-2", "description": "Direct normal irradiance"})

                 #=======================================================
                    
            
                    da_adjusted, qdc_models, adjusted_slices = utils.apply_hourly_qdc_sliding_window_solarfix(
                        da_obs=da_obs,
                        da_model_historical=da_hist,
                        da_model_future=da_fut,
                        var=v,
                        nq=args.nq,
                        kind=kind,
                        window=args.window,
                        solar_alt=ds_nathers_loc["slr_alt"] if v == "rsds" else None,  # optional
                        min_alt=10.0,
                        radiation_mask_source="obs",  # you can keep/adjust
                    )
            
                    # Add adjusted variable
                    var_datasets.append(da_adjusted.rename(v).to_dataset())
            
                    # If we just adjusted rsds, reconstruct rsdsdir/rsdsdif
                    if v == "rsds":
                        da_hist.attrs["units"] = da_obs.attrs.get("units", da_hist.attrs.get("units", ""))
                        da_fut.attrs["units"] = da_obs.attrs.get("units", da_fut.attrs.get("units", ""))
                        
                        # Align fraction to adjusted time (should already match, but be robust)
                        f_dir = f_dir_time.sel(time=da_adjusted.time)
            
                        rsdsdir_adj = (da_adjusted * f_dir).rename("rsdsdir")
                        rsdsdif_adj = (da_adjusted - rsdsdir_adj).rename("rsdsdif")
            
                        # keep metadata tidy
                        rsdsdir_adj.attrs.update({"units": da_adjusted.attrs.get("units", "W m-2"),
                                                  "description": "Direct component derived from rsds using future-model direct fraction"})
                        rsdsdif_adj.attrs.update({"units": da_adjusted.attrs.get("units", "W m-2"),
                                                  "description": "Diffuse component derived as rsds - rsdsdir"})
            
                        # Safety: enforce non-negative and rsdsdir<=rsds
                        rsdsdir_adj = rsdsdir_adj.clip(min=0)
                        rsdsdif_adj = rsdsdif_adj.clip(min=0)

                        # QC: should be ~0 aside from float error
                        closure = np.abs(da_adjusted - (rsdsdir_adj + rsdsdif_adj)).max()
                        print("Max rsds closure error:", float(closure))
            
                        var_datasets.append(rsdsdir_adj.to_dataset())
                        var_datasets.append(rsdsdif_adj.to_dataset())

                        if args.do_plots and v in plot_vars:
                            import matplotlib.pyplot as plt
                            for hour in [0, 12]:
                                fig_file = os.path.join(
                                    args.plot_dir,
                                    os.path.basename(out_file).replace(".nc", f"_AdjFact_{v}_{hour}UTC.pdf"),
                                )
                                if not os.path.exists(fig_file):
                                    utils.plot_qdc_hourly_diagnostics(
                                        da_obs, da_hist, da_fut,
                                        qdc_dict, adjusted_slices,
                                        v, f"{args.rcm}-{gcm}",
                                        loc, hour,
                                        args.scenario, fut_period
                                    )
                                    plt.tight_layout()
                                    plt.savefig(fig_file, bbox_inches="tight")
                                    plt.close()

            
                else:
                    # Fill with unadjusted variables from NatHERS
                    da_fut = ds_nathers[v] #ds_fut_aligned[v]
                    var_datasets.append(da_fut.rename(v).to_dataset())
    
    
            ds_out = xr.merge(var_datasets)
            print(ds_out)
            print(ds_out.to_dataframe()[vars_1hr].describe())
            os.makedirs(os.path.dirname(out_file), exist_ok=True)
            print(f"Writing: {out_file}")
            ds_out.to_netcdf(out_file)
    
            print(f"Done {gcm}, {fut_period} in {(time.time() - t0)/60:.2f} min")
    
    print("\nAll done.")


if __name__ == "__main__":
    main()
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

    vars_qdc = ["tas", "huss", "sfcWind", "psl", "rsdsdir", "rsdsdif"]
    plot_vars = ["tas", "huss", "sfcWind", "psl", "rsds", "rsdsdir", "rsdsdif"]
    vars_1hr = ["tas", "huss", "sfcWind", "psl", "uas", "vas", "clt", "rsds", "rsdsdir", "rsdsdif"]

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
    
            var_datasets = []
    
            for v in vars_1hr:
                print(f"Var: {v}")
    
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

                    slr_alt=None
    
                    if v == "rsds":
                        da_hist.attrs["units"] = da_obs.attrs.get("units", da_hist.attrs.get("units", ""))
                        da_fut.attrs["units"] = da_obs.attrs.get("units", da_fut.attrs.get("units", ""))

                    #============== convert direct radiation ==================
                    if v == "rsdsdir":
                        # Outname 'rsdsdir' even though it's technically 'DNI' but has to match var
                        # names in obs
                        slr_alt = ds_nathers_loc["slr_alt"]          # degrees
                        da_hist_al, slr_alt_al = xr.align(ds_hist[v], slr_alt, join="left")
                        da_fut_al, slr_alt_al = xr.align(ds_fut_aligned[v], slr_alt, join="left")

                        min_proj = np.sin(np.deg2rad(8))   # ~0.087
                        proj = np.sin(np.deg2rad(slr_alt_al))
                        
                        da_hist = xr.where(
                            proj > min_proj,
                            da_hist_al / proj,
                            0.0
                        ).rename("rsdsdir")
                        da_hist.attrs.update({"units": "W m-2", "description": "Direct normal irradiance"})

                        da_fut = xr.where(
                            proj > min_proj,
                            da_fut_al / proj,
                            0.0
                        ).rename("rsdsdir")
                        da_fut.attrs.update({"units": "W m-2", "description": "Direct normal irradiance"})

                        # print(xr.merge([da_obs.rename("rsdsdir_obs").to_dataset(),
                        #                 da_fut.rename("rsdsdir_fut").to_dataset(),
                        #                 da_hist.rename("rsdsdir_hist").to_dataset()]).to_dataframe()[['rsdsdir_obs',
                        #                                                                  'rsdsdir_fut',
                        #                                                                  'rsdsdir_hist']].describe())

                     #=======================================================
                    
                    # da_adj, qdc_dict, adjusted_slices = utils.apply_qdc_sliding_window(
                    #     da_obs=da_obs,
                    #     da_model_historical=da_hist,
                    #     da_model_future=da_fut,
                    #     var=v,
                    #     nq=args.nq,
                    #     kind=kind,
                    #     window=args.window,
                    #     nslots=args.nslots,
                    # )

                    da_adjusted, qdc_models, adjusted_slices = utils.apply_hourly_qdc_sliding_window_solarfix(
                        da_obs=da_obs,
                        da_model_historical=da_hist,
                        da_model_future=da_fut,
                        var=v,
                        nq=args.nq,
                        kind=kind,
                        window=args.window,
                        solar_alt=slr_alt,  # degrees
                        min_alt=10.0,
                        # clip_bounds={"rsdsdir": (0, 1400)},
                        radiation_mask_source="future_model"
                    )

                    var_datasets.append(da_adjusted.rename(v).to_dataset())
                        
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
                    da_fut = ds_fut_aligned[v]
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
    main()                                                                                                                                                                 
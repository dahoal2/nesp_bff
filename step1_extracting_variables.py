#!/usr/bin/env python3

# Ignore warnings as especially xclim/xsdba is sending out a warning about changing libraries.
import warnings
import logging
warnings.filterwarnings('ignore') 
logging.getLogger("distributed.worker.memory").setLevel(logging.ERROR)
logging.getLogger('flox').setLevel(logging.ERROR)

import xarray as xr
import os
import sys
import dask.distributed
import glob
from dask.distributed import Client
from dask.diagnostics import ProgressBar
import tempfile
import dask
import numpy as np
import time
import argparse

######### Import utils and dataset_finder
# Path to cloned NESP repo
sys.path.append('/home/565/dh4185/mn51-dh4185/repos_collab/nesp_bff/')
import utils

# Path to cloned ACS dataset finder repo (https://github.com/AusClimateService/dataset_finder)
sys.path.append('/home/565/dh4185/mn51-dh4185/repos_collab/dataset_finder/')
from dataset_finder import *

# Static metadata dictionaries from utils.py
from utils import locations, model_dict, vars_1hr, vars_day

# Set memory limit depending on job allocation to omit crashing. Generally a 14 CPU/64GB mem job for daily and hourly BARPA and BARRA works fine.
# Since hourly CCAM data processed per year a 2 CPU/9GB mem job works, and is more resource efficient as it needs to run a very long time.
def get_memory_limit():
    try:
        ncpus = int(os.environ.get("PBS_NCPUS", 1))
        # Conservative estimate: 4500 MB per CPU
        print(f"Using {ncpus * 4500}MB memory.")
        return f"{ncpus * 4500}mb"
    except Exception:
        return "32000mb"  # fallback default


##### Main loop 

def main(inargs):

    ##### Parameters
    _rcm = inargs.rcm
    _scenario = inargs.scenario
    start_y = inargs.startYear
    end_y = inargs.endYear
    _freq = inargs.freq
    root_dir = inargs.outputDir
    _nworkers = inargs.nworkers
    #####

    # Sets variable list depending on frequency
    _vars = vars_1hr if _freq == "1hr" else vars_day
    
    # Output dir
    # root_dir = "/g/data/eg3/nesp_bff/step1_raw_data_extraction/"
    if _rcm == "CCAM-v2203-SN":
        out_dir = f"{root_dir}CSIRO-CCAM/"
    else:
        out_dir = f"{root_dir}{_rcm}/"

    # Set up dask
    dask.config.set({
        "distributed.scheduler.worker-saturation": 1, #< This should use the new behaviour which helps with memory pile up
        })
    client = Client(n_workers=_nworkers,
                    threads_per_worker=1,
                    local_directory = tempfile.mkdtemp(),
                    memory_limit=get_memory_limit()) 

    print(f"---------- {_rcm} for '{_freq}' data ----------")
    updated_locations = utils.update_locations(xr.open_dataset(model_dict[_rcm]["sftlf"]).sftlf, locations)
    
    # Iterating though the 12 locations in the updated locations dictionary
    for loc in updated_locations:
        start_time_loc = time.time()
        print(f"========================== {loc} =======================")
        lat = updated_locations[loc]['Lat']
        lon = updated_locations[loc]['Lon']

        gcms_to_process = inargs.gcms if hasattr(inargs, "gcms") and inargs.gcms is not None else list(model_dict[_rcm]["gcms"].keys())

        # Iterating through GCMs for the selected RCM
        for _gcm in gcms_to_process:
            start_time_gcm = time.time()
            print(f"***** {_gcm} *****")
            
            # Boolen to specify if data for hourly CCAM has been rechunked
            should_continue = False

            # Specifying output file name in line was naming convention
            out_file = (
                f"{out_dir}{loc}_"
                f"{model_dict[_rcm]['grid']}_"
                f"{_gcm}_{_scenario}_"
                f"{model_dict[_rcm]['gcms'][_gcm]['mdl_run']}_"
                f"{model_dict[_rcm]['org']}_"
                f"{_rcm}_{model_dict[_rcm]['gcms'][_gcm]['version']}_"
                f"{_freq}_{start_y}-{end_y}.nc"
            )

            # Check if file already exists. If a file needs to be recomputed, it has to be deleted manually first
            if not os.path.exists(out_file):
                print(f"Processing: {out_file}.....")
                var_list = []

                # Iterating through the variables (daily or hourly var_list)
                for _var in _vars:
                    
                    start_time_var = time.time()
                    _timescale = _freq
                    print(f"{_var}: {_vars[_var]}")

                    # Maximum and minimum specific humidity (hussmax, hussmin) is not provided at daily timescale and needs to be 
                    # derived from hourly data.
                    if _timescale == "day" and _var in ['humidity_specific_max', 'humidity_specific_min']:
                        _timescale = "1hr"
                        
                    # BARPA-R is very efficiently chunked four our operation which favours little chunking across time
                    # and lots of chunking along lat and lon. CCAM is chunked for each time step but not at all
                    # along lat and lon dimensions which requires the dataset to be fully loaded. This takes con-
                    # siderable more time to process: BARPA-R day: ~2min, hourly: ~5min. CCAM daily: ~25min, hourly: >7.5h hours
                    # Hence, CCAM hourly data is preprocessed to interim files per year, and then loaded and concatenated.                   
                    if _rcm == "CCAM-v2203-SN" and _timescale == "1hr" and _var not in ['humidity_specific_max',
                                                                                        'humidity_specific_min']:
                        print(f"Use hourly data for {_var}.")
                        start_time_CCAM_1hr = time.time()
                        # Read proprocessd/rechunked hourly CCAM from /scratch/eg3
                        scratch_dir = f"/scratch/eg3/dh4185/rechunked/{_gcm}/{_scenario}/"
                        rechunk_files = sorted(glob.glob(
                            f"{scratch_dir}{_vars[_var][0]}_"
                            f"{model_dict[_rcm]['grid']}_"
                            f"{_gcm}_{_scenario}_"
                            f"{model_dict[_rcm]['gcms'][_gcm]['mdl_run']}_"
                            f"{model_dict[_rcm]['org']}_{_rcm}_"
                            f"{model_dict[_rcm]['gcms'][_gcm]['version']}_1hr_*.nc"))
                        
                        if len(rechunk_files) != 30 and len(rechunk_files) >= 1:
                            print(f"Files don't cover 30 years from {start_y} to {end_y}. Check files and "
                                  f"rerun rechunk_ccam.sh")
                            if rechunk_files:
                                for file in rechunk_files:
                                    print(file)
                                should_continue = True
                                break
                        elif len(rechunk_files) == 0:
                            print(f"No files for GCM {_gcm} and {_var} exists. Run "
                                  f"rechunk_ccam.sh first.")
                            should_continue = True
                            break
                        else:
                            # print(rechunk_files)
                            # Read all years and preprocessing lat/lon selection
                            da = xr.open_mfdataset(rechunk_files, parallel=True,
                                                                preprocess=lambda ds: utils.preprocess_location(ds, lat, lon))[_vars[_var][0]]
                            da = da.chunk({'time': -1}).sel(time=slice(str(start_y),str(end_y)))
                            # Aliging time coordinates (mix of variables at half hour and full hours)
                            da_all = utils.process_time(da,_vars[_var][0],_timescale)
                            var_list.append(da_all.to_dataset())
                                                    
                            print(f"Processing time for {_var}: {((time.time() - start_time_var)/60):.2f} minutes\n")

                    # If BARPA-R or BARRA-R2 at daily or hourly timescale, or CCAM at daily time scale selected process all years at once.
                    elif _rcm in ["BARPA-R","BARRA-R2"] or _rcm == "CCAM-v2203-SN" and _timescale == "day":
                        print(f"Use data from disk for {_var}.")
                        # Get file paths using the ACS dataset finder
                        all_data = get_datasets("ACS_DS",
                                            rcm=_rcm, gcm=_gcm, scenario=_scenario,
                                            grid=model_dict[_rcm]["grid"],
                                            org=model_dict[_rcm]["org"],
                                            mdl_run=model_dict[_rcm]["gcms"][_gcm]["mdl_run"],
                                            ver=model_dict[_rcm]["gcms"][_gcm]["version"],
                                            timescale=_timescale,
                                            year=year_range(start_y, end_y)).select(var=_vars[_var], exact_match=True)

                        # Read all years and preprocessing lat/lon selection
                        da = xr.open_mfdataset(all_data.get_files(), parallel=True,
                                                preprocess=lambda ds: utils.preprocess_location(ds, lat, lon))[_vars[_var][0]]
                        da = da.chunk({'time': -1})
                        
                        # Using hourly huss data to determine daily hussmax and hussmin
                        if _var == 'humidity_specific_max' or _var == 'humidity_specific_min':
                            da = utils.process_humidity(da,_var)
                        da_all = utils.process_time(da,_vars[_var][0],_timescale)
                        var_list.append(da_all.to_dataset())
                    
                    else:
                        print("Inappropriate RCM, GCM, timescale requested. Check Settings.")
                        break
               
                    print(f"Processing time for {_var}: {((time.time() - start_time_var)/60):.2f} minutes\n")

                if should_continue:
                    print("Move to the next GCM.")
                    continue  # Move to the next GCM


                # Remove unwanted variables
                cleaned_list = [da.drop_vars(["bnds", "height", "level_height",
                                              "model_level_number", "sigma"], 
                                              errors="ignore") for da in var_list]
    
                da_var = xr.merge(cleaned_list)
                print(da_var)

                saver = da_var.to_netcdf(out_file,compute=False)
                future = client.persist(saver)
                dask.distributed.progress(future)
                future.compute()
                print(f"Saved: {out_file}")
                                    
                print(f"Processing time for {_rcm}-{_gcm}: {((time.time() - start_time_gcm)/60):.2f} minutes\n")
                
            else:
                print(f"File already exists: {out_file}")
    
        print(f"Done with location: {loc}\n")
    
    print("All processing complete.")

    #< Close the client
    client.close()


if __name__ == '__main__':
    extra_info =""" 
author:
    David Hoffmann, david.hoffmann@bom.gov.au
    Created: 08/07/2025
    Modified: 18/07/2025
"""
    description = """
    Extract variables for NatHERS weather files from reanalysis and climate projections data.
        """    
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)


    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Extract location data from RCMs")
    parser.add_argument("--rcm", required=True, choices = ["BARRA-R2","BARPA-R","CCAM-v2203-SN"], help="RCM name, e.g., BARPA-R")
    parser.add_argument("--gcms", nargs="+", choices=["ACCESS-ESM1-5", "ACCESS-CM2", "CESM2", "CMCC-ESM2", "CNRM-ESM2-1", "EC-Earth3", "NorESM2-MM"],help="Optional list of GCMs to process. If omitted, all available GCMs will be processed.")
    parser.add_argument("--scenario", required=True, help="Scenario, e.g., ssp370 or historical")
    parser.add_argument("--startYear", type=int, required=True, help="Start year")
    parser.add_argument("--endYear", type=int, required=True, help="End year")
    parser.add_argument("--freq", required=True, choices = ['1hr','day'], help="Use hourly (1hr) or daily (day) frequency")
    parser.add_argument("--outputDir", type=str, default="/g/data/eg3/nesp_bff/step1_raw_data_extraction/", help="Output directory on Gadi. Default is '/g/data/eg3/nesp_bff/step1_raw_data_extraction/{rcm}/'")
    parser.add_argument("--nworkers", type=int, required=True, help="Number of wdask workers. 20-25 recommended for 32-64GB job size.")
    args = parser.parse_args()
    main(args)

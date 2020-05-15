"""
Module of analysis and processing functions for the ARNA project
"""

import os
import sys
import xarray as xr
import glob
import numpy as np
import AC_tools as AC
import pandas as pd
from netCDF4 import Dataset
from datetime import datetime as datetime_
#import matplotlib.dates as mdates
#import time
import datetime as datetime
import time
from time import gmtime, strftime
import matplotlib.pyplot as plt
import seaborn as sns
import gc
#import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
import seaborn as sns
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
#from cartopy.mpl.ticker import LatitudeLocator, LongitudeLocator
from shapely.geometry.polygon import LinearRing
import xesmf as xe
import os
#from bs4 import BeautifulSoup
import requests
from PIL import Image, ImageDraw
import PIL
from multiprocessing import Pool
from functools import partial
import matplotlib


def main():
    """
    Main driver for analysis for the ARNA campaign
    """
    # ---  Forecast operation steps
#    do_operational_ARNA_forecast_steps()

    # ---  Forecast (testing etc)
#    get_latest_GEOSCF_fcast_data_ALL()
    # Get lastest GEOS5 forecast
#    get_latest_GEOS5_fcast_data()
    #
#    get_latest_3D_fields_2plot()
    #
#    build_diagnostic_plots_4fcast()
    #
#    do_analysis_of_assimulation_output_JUST_GEOSCF()

    # --- Assimulation
    # Get the assimulation data
#    download_GEOSCF_assim_data()
    # Do analysis for GEOS-CF
#    do_analysis_of_assimulation_output_JUST_GEOSCF()
    # Do plots for analysis for GEOSCF+GEOS5
#    do_analysis_of_assim_GEOSCF_GEOS5()


    # - run the extraction of ARNA flights from GEOS-Cf data
#    extract_GEOS54all_ARNA_flights()
    plt_timeseries_comparisons4ARNA_flights()
    plt_alt_binned_comparisons4ARNA_flights()

    pass # This file is now treated as a module and not called directly



def download_GEOSCF_assim_data():
    """
    Get assimilation (replay) GEOS-CF product
    """
    # - Get cubes of example data around Cape Verde
    # Extract the obervational data?
    retrieve_OPenDAP_data = True
    if retrieve_OPenDAP_data:
        #
        years = (2018, 2019)
#       years = (2018, )    # For testing.
#        years = (2019, )    # For testing.
#        limit_lvls = True  # Just consider the 3 sample levels for now
        for year in years:
            print(year)
            # - Get 3D data
            # Just for five levels of interest
#            get_GEOSCF_data_cubes4campaign(year=year, limit_lvls=True)
            # For all levels
    #        get_GEOSCF_data_cubes4campaign(year=year, limit_lvls=False)
    #        get_GEOSCF_data_cubes4campaign(year=year, limit_lvls=False)
            # For vertical slices on longitudes
            get_GEOSCF_data_cubes4campaign(year=year, limit_lvls=False,
                                           limit_lons=True)

            # - Get 2D data
#            get_data_surface4campaign(year=year)

            # - Compile the downloaded species files into a single NetCDF
    #        compile_OPeNDAP_files_by_day(year=year)


def download_GEOS5_assim_data(mode='assim'):
    """
    Get assimilation (replay) GEOS-CF product
    """
    #
    # - Get cubes of example data around Cape Verde
    # Extract the obervational data?
    retrieve_OPenDAP_data = True
    if retrieve_OPenDAP_data:
        #
#        years = (2018, 2019)
#       years = (2018, )    # For testing.
        years = (2019, )    # For testing.
        limit_lvls = True  # Just consider the 3 sample levels for now
        for year in years:
            print(year)

            # - Get 3D data
#            AC.get_GEOS5_as_ds_via_OPeNDAP()
            # Call once to get the dust data
            get_GEOS5_data_cubes4campaign(year=year)
            # Call again to get the pressure and heights
            collection = 'inst3_3d_asm_Nv'
            vars2use = ['pl', 'h']
            get_GEOS5_data_cubes4campaign(year=year, collection=collection,
                                          vars2use=vars2use)


    # - Now post process files into single doy files.
    post_process_GEOS5_files = True
    if post_process_GEOS5_files:
        #
#        collection = 'inst3_3d_aer_Nv'
        data_folder = get_local_folder('NASA_data')
        folder2use = data_folder + '//GEOS_5/ARNA/assim/cuboids/'
        filestr = 'ARNA_GEO5_{}_Cape_Verde_{}_{:0>3}_*_ALL_lvls.nc'
#        years = (2019, )    # For testing.
        years = (2018, 2019)
        for year in years:
            print(year)
            # Get doys
#            doys2use = np.arange(32, 46)
            # Loop by doy
#            for doy in doys2use:
                # Select files for doy and open
#                files2use = glob.glob( folder2use+filestr.format( '*', year, doy) )
#                files2use = [i for i in files2use if 'REGRIDDED' not in i ]
#                print( files2use )
            #
            regrid_GEOS5_files_in_folder(folder=folder, doys2use=doys2use,
                                         year=year )



def do_analysis_of_assim_GEOSCF_GEOS5(do_core_lon_analysis=True,
                                      do_horizon_analysis=False,
                                      ):
    """
    Do analysis using dust from GEOS-5 and NOy from GEOS-CF
    """
    # Data locations
    NASA_data = get_local_folder('NASA_data')
    G5_folder = NASA_data + '/GEOS_5/ARNA/'
    GCF_folder = NASA_data + '/GEOS_CF/ARNA/'
    # - Do analysis on a altitude slice basis
    if do_horizon_analysis:
        # Do analysis on a year by year basis
        years = (2018, 2019)
#        years = (2018, )    # For testing.
#        years = (2019, )    # For testing.
        limit_lvls = True  # Just consider the 5 sample levels for now
        for year in years:
            print(year)
            # - Get GEOS-CF data
            folder2use = GCF_folder + '/assim/alt_slice/'
            files2use = glob.glob('{}/ARNA*{}_0*.nc'.format(folder2use, year))
            files2use = list(sorted(files2use))
            # Open all the files as a single dataset
            ds = xr.open_mfdataset(files2use)
            print(ds.time)
            # Update the variable names to be in GEOS-Chem format
            ds = update_GEOSCF_names2GEOSChem_names(ds=ds)
            # Update the units to GEOS-Chem/interpretable ones
            ds = convert_GEOSCF_units(ds=ds)
            # Add extra analysis species
            ds = add_extra_derived_species(ds=ds)
            # - Get GEOS-5 data
            # Get data
            folder2use = G5_folder + '/assim/cuboids/'
            files2use = glob.glob('{}/ARNA*{}_0*REGRID*.nc'.format(folder2use, year))
            files2use = list(sorted(files2use))
            # Open all the files as a single dataset
            ds5 = xr.open_mfdataset(files2use)
            print(ds5.time)
            # - Process the data
            # Add to the GEOS-CF dataset
            ds['Dust'] = ds5['Dust']
            ds['Dust'].attrs = ds5['Dust'].attrs
            # Only plot points where both data sources are present
            only_plot_where_GEOS5 = True
            if only_plot_where_GEOS5:
                time_boolean = [i in ds5.time.values for i in ds.time.values]
                ds = ds.isel(time=time_boolean)
            # Plot and then save together
            folder4plots_str = '{}/GEOS_CF/ARNA/assim/plots_with_GEOS5/alt_slice/'
            folder4plots = folder4plots_str.format(NASA_data)
            # - Plot the two layer plot (setting plotting to use the Dust settings.)
            plot_spatial_concs_2layer(ds, verbose=True, var2plot1='NOy',
                                      var2plot2='Dust', folder=folder4plots)

    # - Do analysis on a longitudinal slice basis
    if do_core_lon_analysis:
        years = (2018, 2019)
#        years = (2018, )    # For testing.
        limit_lvls = True  # Just consider the 3 sample levels for now
        for year in years:
            print('plotting for:',  year)
            # - Get data and process
            folder2use = GCF_folder + '/assim/lon_slice/'
            # Get the data for previous years
            files2use = glob.glob('{}/ARNA*{}_0*.nc'.format( folder2use, year ) )
            files2use = list(sorted(files2use))
            # Open all the files as a single dataset
            ds = xr.open_mfdataset(files2use)
            print(ds.time)
            # Update the variable names to be in GEOS-Chem format
            ds = update_GEOSCF_names2GEOSChem_names(ds=ds)
            # Update the units to GEOS-Chem/interpretable ones
            ds = convert_GEOSCF_units(ds=ds)
            # Add extra analysis species
            ds = add_extra_derived_species(ds=ds)
            # - Plot up spatial averages
            # Reset the environment to remove sea-born settings
            sns.reset_orig()
            # - Get GEOS-5 data
            # Get data
            folder2use = G5_folder + '/assim/cuboids/'
            fstr = '{}/ARNA*{}_0*REGRID*.nc'
            files2use = glob.glob(fstr.format( folder2use, year ) )
            files2use = list(sorted(files2use))
            # Open all the files as a single dataset
            ds5 = xr.open_mfdataset(files2use)
            print(ds5.time)
            # Only consider the slice for lon
            lons2use = [-18 ,-19.5, -21, -22.5, -24, -25.5]
            ds5 = ds5.isel(lon=[(i in lons2use) for i in ds5.lon]   )
            # Add to the GEOS-CF dataset
            ds['Dust'] = ds5['Dust']
            ds['Dust'].attrs = ds5['Dust'].attrs
            only_plot_where_GEOS5 = True
            if only_plot_where_GEOS5:
                time_boolean = [i in ds5.time.values for i in ds.time.values ]
                ds = ds.isel(time=time_boolean)
            # Attributes
            attrs = ds.lev.attrs
            attrs['units'] = 'millibar'
            ds.lev.attrs = attrs
            # Plot and then save together
            folder4plots_str = '{}/GEOS_CF/ARNA/assim/plots_with_GEOS5/lon_slice/'
            folder4plots = folder4plots_str.format(NASA_data)
            # GEOS-CF
            plt_spatial_2layer_vertical_lon(ds, folder=folder4plots,
                                                   var2plot1='NOy',
                                                   var2plot2='Dust',)



def do_analysis_of_assimulation_output_JUST_GEOSCF(do_core_lon_analysis=True,
                                                   do_horizon_analysis=False
                                                   ):
    """
    Do analysis of assimilation (replay) GEOS-CF product
    """
    # - Do analysis on a altitude slice basis
    # Where to look for data
    NASA_data = get_local_folder('NASA_data')
    folder = NASA_data + 'GEOS_CF/ARNA/assim/alt_slice/'
    if do_horizon_analysis:
        # Do analysis on a year by year basis
    #    years = (2018, 2019)
        years = (2018, )    # For testing.
        limit_lvls = True  # Just consider the 3 sample levels for now
        for year in years:
            print(year)
            # - Get data and process
            # Get the data for previous years
            files2use = glob.glob('{}/ARNA*{}_0*.nc'.format( folder, year ) )
            # Open all the files as a single dataset
            ds = xr.open_mfdataset(files2use)
            print(ds.time)
            # Update the variable names to be in GEOS-Chem format
            ds = update_GEOSCF_names2GEOSChem_names(ds=ds)
            # Update the units to GEOS-Chem/interpretable ones
            ds = convert_GEOSCF_units(ds=ds)
            # Add extra analysis species
            ds = add_extra_derived_species(ds=ds)
            # - Plot up spatial averages
            # Reset the environment to remove sea-born settings
            sns.reset_orig()
            #
            plt_avg_spatial_by_lvl(ds, year, verbose=True)
            # Plot the two layer plot
            plot_spatial_concs_2layer(ds, verbose=True)

            # - General general stats on NOy and dust
            # Now use seaborn settings
            sns.set(color_codes=True)
            sns.set_context("talk")
            # ... for the region
            region = 'Cape_Verde'
            summarise_stats_on_species_in_ds4lvls(ds, region=region, year=year)
            # ... for a smaller area around Cape Verde
            dsCV = select_smaller_area_around_Cape_Verde(ds.copy())
            region = 'Cape_Verde_LOCAL'
            summarise_stats_on_species_in_ds4lvls(dsCV, region=region,
                                                  year=year)
            # - Plot up PDFs of distribution
            # ... for the region
            region = 'Cape_Verde'
            PDF_on_species_in_ds4lvls(ds, region=region, year=year)
            # ... for a smaller area around Cape Verde
            region = 'Cape_Verde_LOCAL'
            PDF_on_species_in_ds4lvls(dsCV, region=region, year=year)
            # NOTE: region is defined as a smaller shape than the whole map
            # Plot up a map to show the subsampled region
            plt_smaller_area_around_CVAO(ds)

            # - Plot up timeseries plots for key species
            #
            region = 'Cape_Verde_LOCAL'
            plt_timeseries4ds(dsCV, region=region, year=year)


    # - Do analysis on a longitudinal basis
    # Where to look for data
    folder = NASA_data + '/GEOS_CF/ARNA//assim/lon_slice/'
    # Do vertical
    if do_core_lon_analysis:
        # Do analysis on a year by year basis
#        years = (2018, 2019)
        years = (2018, )    # For testing.
        limit_lvls = True  # Just consider the 3 sample levels for now
        for year in years:
            print('processing for:',  year)
            # - Get data and process
            # Get the data for previous years
            files2use = glob.glob('{}/ARNA*{}_0*.nc'.format( folder, year ) )
            # Open all the files as a single dataset
            ds = xr.open_mfdataset(files2use)
            print(ds.time)
            # Update the variable names to be in GEOS-Chem format
            ds = update_GEOSCF_names2GEOSChem_names(ds=ds)
            # Update the units to GEOS-Chem/interpretable ones
            ds = convert_GEOSCF_units(ds=ds)
            # Add extra analysis species
            ds = add_extra_derived_species(ds=ds)

            # - Plot up spatial averages
            # Reset the environment to remove sea-born settings
            sns.reset_orig()
    #        plt_avg_spatial_by_lvl(ds, year, verbose=True)
            # Plot NOy and dust separately
            # Plot together
            folder4plots = './'
            plt_spatial_2layer_vertical_lon(ds, folder=folder4plots)



def get_all_surface_data_for_CVAO4dates():
    """
    Get all the GEOS-CF surface data for CVAO
    """
    # TODO
    pass


def summarise_stats_on_species_in_ds4lvls(ds, region='Cape_Verde', extr_str='',
                                          vars2use=None, year=2018):
    """
    Summarise the stats on species at different heights
    """
    # Which variables to include in analysis?
    if not isinstance(vars2use, list):
        vars2use = [i for i in ds.data_vars]
    # Setup a dataframe to store data in
    df = pd.DataFrame()
    # Loop by variable
    vars2use = [i for i in ds.data_vars]
    for var in vars2use:
        # Loop by level
        lvls = np.array(ds.lev.values)
        for lvl in lvls:
            print(var, lvl)
            # Setup a name to save data to
            varname = '{}-{:.0f}hPa'.format(var, lvl)
            # Flatten the data
            vals = ds[var].sel(lev=lvl).values.flatten()
            S = pd.Series(vals).describe()
            # Save to the main DataFrame
            df[varname] = S
            del vals, S
    # Save summary statistics
    savename = 'ARNA_stats_species_{}_{}_{}.csv'.format(region, year, extr_str)
    df.to_csv(savename)


if __name__ == "__main__":
    main()

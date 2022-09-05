#!/usr/bin/python
"""
Module of driver plotting functions for during the ARNA campaign
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
import datetime as datetime
import time
from time import gmtime, strftime
import matplotlib.pyplot as plt
import seaborn as sns
import gc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
import seaborn as sns
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from shapely.geometry.polygon import LinearRing
import xesmf as xe
import os
import requests
from PIL import Image, ImageDraw
import PIL
from multiprocessing import Pool
from functools import partial
import matplotlib
# import ARNA analysis/campaign code as a module
import arna as ar
#from arna import mk_core_plts4fcast_GEOSCF_GEOS5
from arna import plot_spatial_concs_2layer, plt_spatial_2layer_vertical_lat
from arna import plt_spatial_2layer_vertical_lon


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
    pass


def do_operational_ARNA_forecast_steps(dt=None):
    """
    Do daily work pipeline for ARNA project - REDUNDENT: Now run externally
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Make the folder structure for plots and data
    ar.mk_folder_structure4fcast(dt=dt)

    # Get GEOS-CF data
#    get_latest_GEOSCF_fcast_data_ALL()
#    get_latest_GEOSCF_fcast_data_alt_slice()
#    get_latest_GEOSCF_fcast_data_lat_slice()
#    get_latest_GEOSCF_fcast_data_lon_slice()
    # Get GEOS-5 data
#    get_latest_GEOS5_fcast_data()

    # Check all the files successfully downloaded. if not get missing files.
    n_failed_doys = ar.check4failed_downloads(dt=dt)
    print(n_failed_doys)

    # Process GEOS-5 data - folder currently needs to be provided!
    folder = ar.get_GEOS_data_folder4dt(dt=dt, product='GEOS_5')
    ar.regrid_GEOS5_files_in_folder(folder=folder)

    # Make 2 layer plots for GEOS5+GEOSCF
#    ar.mk_core_plts4fcast_GEOSCF_GEOS5( None )

    # Move plots to <>
    # /GEOS_CF/ARNA/fcast/2019_12_15/plots.incGEOS5
    # Lon_slice / alt_slice (frames folder)

    #
#    ar.regrid_files_then_do_ind_plots()
#    ar.mk_ind_plsts4fcast_GEOSCF_GEOS5()

    # Animate forecast data to videos.

    # Get SDS-WAS forecast plots (download directly)

    # Get SDS-WAS forecast netcdf files?

    # Get GEOS5 diagnostic plots
#    ar.get_latest_GEOS5_diagnostics()

    # Move folders into synced core +24 and +48 folders

    # Sync these files to websfiles
    ar.mv_plots2webfiles(dt=dt)


def mk_missing_ARNA_plots4dts(dts=None, mk_plots=True):
    """
    Make plots for specific datetimes
    """
    # Use the last 18 days, unless specific dates are provided
    if isinstance(dts, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
        # Use the last 18 days
        dts = [AC.add_days(dt, i*-1) for i in range(1, 18)]
    # Loop a call to checking and plotting for a given date
    for dt in dts:
        ar.mk_missing_ARNA_plots4dt(dt=dt, mk_plots=mk_plots)


def mk_missing_ARNA_plots4dt(dt=None, mk_plots=True):
    """
    Make plots for a specific datetime
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Check if the folders are not complete
    dstr = dt.strftime('%Y/%m/%d %H:%M')
    print('Checking plots have been made for {}'.format(dstr))
    # Check if the plots are already present
    folders2get = ar.which_plot_folders_are_not_complete4dt(dt=dt)
    # Now only download those that are not complete
    if len(folders2get) == 0:
        print('All plots present and correct')
    else:
        print('WARNING: None/not all plots present! - ', folders2get)
        # Hardwire folder names to check
        subfolders2check = {
            'alt_slice': 0,
            'alt_slice.zoomed': 1,
            'lat_slice': 2,
            'lon_slice': 3,
        }
        # Regrid the files...
        folder = ar.get_GEOS_data_folder4dt(dt=dt, product='GEOS_5')
        ar.regrid_GEOS5_files_in_folder(folder=folder)
        # Make the folder structure for the plots
        mk_folder_structure4fcast(dt=dt)
        # Check sub folder by folder and do call the specific fucntion for given case
        for subfolder in subfolders2check.keys():
            if (subfolder in folders2get) and mk_plots:
                print("Making '{}' plots for '{}'".format(subfolder, dstr))
                plot_type = subfolders2check[subfolder]
                # What type of plot is the subfolder for?
                mk_core_plts4fcast_GEOSCF_GEOS5(plot_type, dt=dt,
                                                do_core_alt_analysis=False,
                                                do_zoomed_alt_analysis=False,
                                                do_core_lon_analysis=False,
                                                do_core_lat_analysis=False,
                                                )
        # Make the individual plots too
        subfolder = 'alt_slice.individual'
        if (subfolder in folders2get) and mk_plots:
            print("Making '{}' plots for '{}'".format(subfolder, dstr))
            ar.regrid_files_then_do_ind_plots(dt=dt)

        # Also check the GMAO folders
        # NOTE: this will not work, as the files on the GMAO site are latest ones
        # A new function could be written to GMAO 2D plots, but not datagrams
#        subfolder = 'plots.GMAO'
#        if subfolder in folders2get:
#            print("Making '{}' plots for '{}'".format(subfolder,dstr ))
#            GMAO
#            get_latest_GEOS5_diagnostics(dt=dt)


def regrid_files_then_do_core_plots(dt=None, plt_in_parallel=False):
    """
    Check files are present, then regrid GEOS5 files and plot core plots
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Check that the files are all present and correct...
    n_failed_doys = ar.check4failed_downloads(dt=dt, retry=False)
    print(n_failed_doys)
    plot_anyway = True
    if (n_failed_doys == 0) or plot_anyway:
        # Regrid the files...
        folder = ar.get_GEOS_data_folder4dt(dt=dt, product='GEOS_5')
        ar.regrid_GEOS5_files_in_folder(folder=folder)
        # Make the folder structure for the plots
        ar.mk_folder_structure4fcast(dt=dt)
        # Plot the data by species in parallel or series
        if plt_in_parallel:
            # Parallelise over variables
            plot_type = {
                'alt_slice': 0,
                'alt_slice.zoomed': 1,
                'lat_slice': 2,
                'lon_slice': 3,
            }
            plot_type.keys()
            pstr = 'Plotting variables in parallel across a pool of {}'
            print(pstr.format(len(plot_type)))
            # Setup this number of pools
            p = Pool(len(plot_type))
            # plot all of the plots types at the same time
            p.map(partial(mk_core_plts4fcast_GEOSCF_GEOS5, dt=dt,
                          do_core_alt_analysis=False,
                          do_zoomed_alt_analysis=False,
                          do_core_lon_analysis=False,
                          do_core_lat_analysis=False,
                          ), plot_type)
            # Close the pool
            p.close()
        else:
            mk_core_plts4fcast_GEOSCF_GEOS5(None, dt=dt)
    else:
        # Print a warning to screen if case failed.
        pstr = 'WARNING: NO REGRID/PLOTTING AS ALL FILES NOT PRESENT/CORRECT!'
        print(pstr)
        print(dt)


def regrid_files_then_do_ind_plots(dt=None, plt_in_parallel=False):
    """
    Check files are present, then regrid GEOS5 files and individual plots
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Check that the files are all present and correct...
    n_failed_doys = ar.check4failed_downloads(dt=dt, retry=False)
    print(n_failed_doys)
    plot_anyway = True
    if (n_failed_doys == 0) or plot_anyway:
        # Regrid the files...
        folder = ar.get_GEOS_data_folder4dt(dt=dt, product='GEOS_5')
        ar.regrid_GEOS5_files_in_folder(folder=folder)
        # Make the folder structure for the plots
        ar.mk_folder_structure4fcast(dt=dt)
        # Plot the data by species in parallel or series
        if plt_in_parallel:
            # Parallelise over variables
            vars2use = ['NOy', 'NO2', 'O3', 'CO', 'Dust']
            vars2use = [[i] for i in vars2use]  # make sure these are lists!
            pstr = 'Plotting variables in parallel across a pool of {}'
            print(pstr.format(len(vars2use)))
            # Setup this number of pools
            p = Pool(len(vars2use))
            # Now map partially across the pool
            p.map(partial(mk_ind_plsts4fcast_GEOSCF_GEOS5,
                          dt=dt,), vars2use)
            # Close the pool
            p.close()
        else:
            mk_ind_plsts4fcast_GEOSCF_GEOS5(None, dt=dt)
    else:
        #
        pstr = 'WARNING: STOPPED REGRID/PLOT - NOT ALL FILES PRESENT/CORRECT!'
        print(pstr)
        print(dt)


def mk_ind_plsts4fcast_GEOSCF_GEOS5(vars2use, dt=None,
                                    only_plot_where_GEOS5=True):
    """
    Make alt slice plots on a individual species basis from GEOS-CF output
    """
    # Get most recent fcast data as a dataset
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Use default variables if a list not provided
    if isinstance(vars2use, type(None)):
        vars2use = ['NOy', 'NO2', 'O3', 'CO', 'Dust']
    # Get GEOS-CF folder for datettime
    collection = 'chm_inst_1hr_g1440x721_p23'
    G5_data_4dt = ar.get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                             collection=collection)
    # Get GEOS-5 folder for datetime
    collection = 'inst3_3d_aer_Np'
    G5_data_4dt = ar.get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                             collection=collection)
    # Get plots folder
    G5_plots_4dt = ar.get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                              inc_collection=False)
    G5_plots_4dt += '/plots/'
    # - Get GEOS-CF data
    ds = ar.get_most_recent_GEOSCF_data(dt=dt)
    t0_CF = ds.time.values[0]
    t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
    # - Get GEOS-5 data
    fstr = '{}/ARNA*_{}_*REGRID*.nc'
    files2use = glob.glob(fstr.format(G5_data_4dt, dt.year))
    files2use = list(sorted(files2use))
    # Open all the files as a single dataset
    ds5 = xr.open_mfdataset(files2use)
    t0_G5 = ds5.time.values[0]
    t0_G5_str = AC.dt64_2_dt([t0_G5])[0].strftime('%Y/%m/%d %H:%M')
    print(ds5.time, t0_G5_str)
    # Limit alt to certain levels
    hPa_heights = [1000, 900, 800, 700, 600, 500]
    ds5 = ds5.isel(lev=[(i in hPa_heights) for i in ds5.lev])
    # Add GEOS5 Dust to the GEOS-CF dataset
    ds['Dust'] = ds5['Dust']
    ds['Dust'].attrs = ds5['Dust'].attrs
    # Only plot for points where both dust and noy data exist
    if only_plot_where_GEOS5:
        ds = ds.isel(time=[i in ds5.time.values for i in ds.time.values])
    # Add a string for when the GEOS runs were initiated
    assert t0_CF_str == t0_G5_str
    extr_title_str = ' (GEOS-CF/5 from {})'.format(t0_G5_str)
    # Folder to save plots
    folder = G5_plots_4dt + '/alt_slice.individual/'
    # Make single species plots for GEOSCF and GEOS5
    ar.plot_individual_spec_alt_slices(ds, vars2use=vars2use, folder=folder,
                                       extr_title_str=extr_title_str)


def mk_core_plts4fcast_GEOSCF_GEOS5(plot_type, dt=None,
                                    do_core_alt_analysis=True,
                                    do_zoomed_alt_analysis=True,
                                    do_core_lon_analysis=True,
                                    do_core_lat_analysis=True,
                                    only_plot_where_GEOS5=True,
                                    testing_mode=False,
                                    ):
    """
    Do analysis using dust from GEOS-5 and NOy from GEOS-CF

    TODO: split off plotting functions to be called separately.
    (e.g. horizontal, vertical, )
    """
    # Set a datetime to use, if not provided
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # GEOS-5 folder
    collection = 'inst3_3d_aer_Np'
    G5_data_4dt = ar.get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                             collection=collection)
    # Root plot saving folder
    G5_plot_dir = ar.get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                             inc_collection=False)
    G5_plot_dir += '/plots/'

    # - Do analysis on a altitude slice basis
    if do_core_alt_analysis or (plot_type == 0):
        # Get GEOS-CF data
        ds = ar.get_most_recent_GEOSCF_data(dt=dt)
        t0_CF = ds.time.values[0]
        t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
        # - Get GEOS-5 data
        fstr = '{}/ARNA*_{}_*REGRID*.nc'
        files2use = glob.glob(fstr.format(G5_data_4dt, dt.year))
        files2use = list(sorted(files2use))
        # Open all the files as a single dataset
        # , combine='nested', concat_dim='time')
        ds5 = xr.open_mfdataset(files2use)
        t0_G5 = ds5.time.values[0]
        t0_G5_str = AC.dt64_2_dt([t0_G5])[0].strftime('%Y/%m/%d %H:%M')
        print(ds5.time, t0_G5_str)
        # - Process GEOS-CF/5 data
        # Limit alt to certain levels
        hPa_heights = [1000, 900, 800, 700, 600, 500]
        ds5 = ds5.isel(lev=[(i in hPa_heights) for i in ds5.lev])
        # Add the GEOS-5 data to the GEOS-CF dataset
        ds['Dust'] = ds5['Dust']
        ds['Dust'].attrs = ds5['Dust'].attrs
        #  Only plot for points where both dust and noy data exist
        if only_plot_where_GEOS5:
            ds = ds.isel(time=[i in ds5.time.values for i in ds.time.values])
        # Add a string for when the GEOS runs were initiated
        assert t0_CF_str == t0_G5_str
        extr_title_str = ' (GEOS-CF/5 from {})'.format(t0_G5_str)
        # Set save location for plots
        if testing_mode:
            folder = './'  # unhash for testing
        else:
            folder = G5_plot_dir + '/alt_slice/'
        # Plot the two layer plot (setting plotting to use the Dust settings.)
        plot_spatial_concs_2layer(ds, verbose=True,
                                  # Plot NOy, then dust
                                  #                                  var2plot1='NOy', var2plot2='Dust',
                                  # Plot Dust, the NOy
                                  var2plot1='Dust', var2plot2='NOy',
                                  region='Cape_Verde', folder=folder,
                                  testing_mode=testing_mode,
                                  extr_title_str=extr_title_str)

    if do_zoomed_alt_analysis or (plot_type == 1):
        # Get GEOS-CF data
        ds = ar.get_most_recent_GEOSCF_data(dt=dt)
        t0_CF = ds.time.values[0]
        t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
        # - Get GEOS-5 data
        # Get data
        fstr = '{}/ARNA*_{}_*REGRID*.nc'
        files2use = glob.glob(fstr.format(G5_data_4dt, dt.year))
        files2use = list(sorted(files2use))
        # Open all the files as a single dataset
        # , combine='nested', concat_dim='time')
        ds5 = xr.open_mfdataset(files2use)
        t0_G5 = ds5.time.values[0]
        t0_G5_str = AC.dt64_2_dt([t0_G5])[0].strftime('%Y/%m/%d %H:%M')
        print(ds5.time, t0_G5_str)
        # Limit alt to certain levels
        hPa_heights = [1000, 900, 800, 700, 600, 500]
        ds5 = ds5.isel(lev=[(i in hPa_heights) for i in ds5.lev])
        # Add to the GEOS-CF dataset
        ds['Dust'] = ds5['Dust']
        ds['Dust'].attrs = ds5['Dust'].attrs
        #  Only plot for points where both dust and noy data exist
        if only_plot_where_GEOS5:
            ds = ds.isel(time=[i in ds5.time.values for i in ds.time.values])
        # AddNOy * Dust (GEOS5)
        ds['NOy*Dust'] = ds['Dust'] * ds['NOy']
        # Add a string for when the GEOS runs were initiated
#        extr_title_str = ' (G-CF:{}, G-5:{})'.format( t0_CF_str, t0_G5_str )
        assert t0_CF_str == t0_G5_str
        extr_title_str = ' (GEOS-CF/5 from {})'.format(t0_G5_str)
        # Also put the zoomed in plots and save to folder
        if testing_mode:
            folder = './'  # unhash for testing
        else:
            folder = G5_plot_dir + 'alt_slice.zoomed/'
        # Plot the two layer plot (setting plotting to use the Dust settings.)
        plot_spatial_concs_2layer(ds, verbose=True,
                                  # Plot NOy, then dust
                                  #                                  var2plot1='NOy', var2plot2='Dust',
                                  # Plot Dust, the NOy
                                  var2plot1='Dust', var2plot2='NOy',
                                  region='Cape_Verde_Flying',
                                  add_max_vals_as_txt=True,
                                  folder=folder,
                                  testing_mode=testing_mode,
                                  extr_title_str=extr_title_str)

    # - Do analysis on a longitudinal slice basis
    if do_core_lon_analysis or (plot_type == 2):
        # - Get GEOS-CF data
        filestr = 'ARNA*{}*-18.0_*.nc'
        ds = ar.get_most_recent_GEOSCF_data(dt=dt, filestr=filestr)
        t0_CF = ds.time.values[0]
        t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
        # - Get GEOS-5 data
        # Get data
        fstr = '{}/ARNA*_{}_*REGRID*.nc'
        files2use = glob.glob(fstr.format(G5_data_4dt, dt.year))
        files2use = list(sorted(files2use))
        # Open all the files as a single dataset
        # , combine='nested', concat_dim='time')
        ds5 = xr.open_mfdataset(files2use)
        t0_G5 = ds5.time.values[0]
        t0_G5_str = AC.dt64_2_dt([t0_G5])[0].strftime('%Y/%m/%d %H:%M')
        print(ds5.time, t0_G5_str)
        # Only consider the slices for lon
        lons2use = [-18, -19.5, -21, -22.5, -24, -25.5]
        ds5 = ds5.isel(lon=[(i in lons2use) for i in ds5.lon])
        # Add to the GEOS-CF dataset
        ds['Dust'] = ds5['Dust']
        ds['Dust'].attrs = ds5['Dust'].attrs
        # Only plot for the GEOS5 time steps.
        if only_plot_where_GEOS5:
            ds = ds.isel(time=[i in ds5.time.values for i in ds.time.values])
        # Attributes
        attrs = ds.lev.attrs
        attrs['units'] = 'millibar'
        ds.lev.attrs = attrs
        # Add a string for when the GEOS runs were initiated
        assert t0_CF_str == t0_G5_str
        extr_title_str = ' (GEOS-CF/5 from {})'.format(t0_G5_str)
        # Also put the lon slice in plots and save to folder
        if testing_mode:
            folder = './'  # unhash for testing
        else:
            folder = G5_plot_dir + 'lon_slice/'
        # GEOS-CF
        plt_spatial_2layer_vertical_lon(ds, folder=folder,
                                        # Plot NOy, then dust
                                        #                                               var2plot1='NOy', var2plot2='Dust',
                                        # Plot Dust, the NOy
                                        var2plot1='Dust',
                                        var2plot2='NOy',
                                        testing_mode=testing_mode,
                                        extr_title_str=extr_title_str)

    # - Do analysis on a latitude slice basis
    if do_core_lat_analysis or (plot_type == 3):
        # - Get GEOS-CF data
        filestr = 'ARNA*{}*ALL_lvls_lats_*.nc'
        ds = ar.get_most_recent_GEOSCF_data(dt=dt, filestr=filestr)
        t0_CF = ds.time.values[0]
        t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
        # - Get GEOS-5 data
        # Get data
        fstr = '{}/ARNA*_{}_*REGRID*.nc'
        files2use = glob.glob(fstr.format(G5_data_4dt, dt.year))
        files2use = list(sorted(files2use))
        # Open all the files as a single dataset
        # , combine='nested', concat_dim='time')
        ds5 = xr.open_mfdataset(files2use)
        t0_G5 = ds5.time.values[0]
        t0_G5_str = AC.dt64_2_dt([t0_G5])[0].strftime('%Y/%m/%d %H:%M')
        print(ds5.time, t0_G5_str)
        # Only consider the slices for lon
        lats2use = [12, 13, 14, 15, 16, 17]
        ds5 = ds5.isel(lat=[(i in lats2use) for i in ds5.lat])
        # Add to the GEOS-CF dataset
        ds['Dust'] = ds5['Dust']
        ds['Dust'].attrs = ds5['Dust'].attrs
        # Only plot for the GEOS5 time steps.
        if only_plot_where_GEOS5:
            ds = ds.isel(time=[i in ds5.time.values for i in ds.time.values])
        # Attributes
        attrs = ds.lev.attrs
        attrs['units'] = 'millibar'
        ds.lev.attrs = attrs
        # Add a string for when the GEOS runs were initiated
        assert t0_CF_str == t0_G5_str
        extr_title_str = ' (GEOS-CF/5 from {})'.format(t0_G5_str)
        # Also put the lat slice in plots and save to folder
        if testing_mode:
            folder = './'  # unhash for testing
        else:
            folder = G5_plot_dir + 'lat_slice/'
        # GEOS-CF
        plt_spatial_2layer_vertical_lat(ds, folder=folder,
                                        # Plot NOy, then dust
                                        #                                             var2plot1='NOy', var2plot2='Dust',
                                        # Plot Dust, the NOy
                                        var2plot1='Dust',
                                        var2plot2='NOy',
                                        testing_mode=testing_mode,
                                        extr_title_str=extr_title_str)


def get_latest_GEOSCF_fcast_data_ALL(dt=None, just_check_yesterday=True,
                                     rm_existing_file=False, debug=True):
    """
    Get the latest GEOS-CF forecast data (if new data available)
    """
    if isinstance(dt, type(None)):
        if just_check_yesterday:
            # Just use yesterday for Now
            TNow = AC.time2datetime([gmtime()])[0]
            dt = AC.add_days(TNow, -1)
            dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
        else:
            # What is the last data there is data locally for?
            last_local_data_date = ar.what_is_latest_data_locally(debug=debug)
            assert_err = 'Stopping as no data for last 5 days'
            assert not isinstance(last_local_data_date, type(None)), assert_err
            # When was the last forecast?
            last_fcast_start = ar.when_should_GEOSCF_have_last_run_from()
            # Use this as the datetime to download
            dt = last_fcast_start
            ass_str = 'The most recent data is already locally available'
            assert last_fcast_start > last_local_data_date, ass_str
            # Download data if available from NASA
#            lastest_GEOSCF_avail = check_if_latest_GEOSCF_is_available(last_fcast_start)
    else:
        # Using the provided datetime (dt)
        pass

    # - Retrieve the latest GEOS-CF data - sliced by altitude
    ar.download_GEOSCF_fcast_data4date(dt=dt, limit_lvls=True,
                                       rm_existing_file=rm_existing_file,
                                       limit_lons=False)
    # - Retrieve the latest GEOS-CF data - sliced by longitude
    ar.download_GEOSCF_fcast_data4date(dt=dt, limit_lvls=False,
                                       rm_existing_file=rm_existing_file,
                                       limit_lons=True)
    # - Retrieve the latest GEOS-CF data - sliced by latitude
    ar.download_GEOSCF_fcast_data4date(dt=dt, limit_lvls=False,
                                       rm_existing_file=rm_existing_file,
                                       limit_lons=False, limit_lats=True)


if __name__ == "__main__":
    main()

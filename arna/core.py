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

    # ---  assimulation
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


def do_operational_ARNA_forecast_steps(dt=None):
    """
    Do daily work pipeline for ARNA project - REDUNDENT: Now run externally
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Make the folder structure for plots and data
    mk_folder_structure4fcast(dt=dt)

    # Get GEOS-CF data
#    get_latest_GEOSCF_fcast_data_ALL()
#    get_latest_GEOSCF_fcast_data_alt_slice()
#    get_latest_GEOSCF_fcast_data_lat_slice()
#    get_latest_GEOSCF_fcast_data_lon_slice()
    # Get GEOS-5 data
#    get_latest_GEOS5_fcast_data()

    # Check all the files successfully downloaded. if not get missing files.
    n_failed_doys = check4failed_downloads(dt=dt)
    print(n_failed_doys)

    # Process GEOS-5 data - folder currently needs to be provided!
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5')
    regrid_GEOS5_files_in_folder(folder=folder)

    # Make 2 layer plots for GEOS5+GEOSCF
#    mk_core_plts4fcast_GEOSCF_GEOS5( None )

    # Move plots to <>
    # /GEOS_CF/ARNA/fcast/2019_12_15/plots.incGEOS5
    # Lon_slice / alt_slice (frames folder)

    #
#    regrid_files_then_do_ind_plots()
#    mk_ind_plsts4fcast_GEOSCF_GEOS5()


    # Animate forecast data to videos.

    # Get SDS-WAS forecast plots (download directly)

    # Get SDS-WAS forecast netcdf files?

    # Get GEOS5 diagnostic plots
#    get_latest_GEOS5_diagnostics()


    # Move folders into synced core +24 and +48 folders

    # Sync these files to websfiles
    mv_plots2webfiles(dt=dt)


def mk_core_plot_folders_then_mv2webfiles(dt=None, mv2webfiles=True,
                                          debug=True):
    """
    Make core folders (+?? hours), then move these to webfiles
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # - mv the +24/+48 files into the core folders
    dstr = dt.strftime('%Y/%m/%d %H:%M')
    copy_files2core_plot_folders(dt=dt)
    # - Now move the folder to webfiles
    if mv2webfiles:
        TNow = AC.time2datetime( [gmtime()] )[0]
        pstr = "Started moving files for {} to webfiles @ {}"
        print(  pstr.format( dstr, TNow.strftime('%Y/%m/%d %H:%M') ) )
        mv_plots2webfiles(dt=dt)
    # - Now move the files to google drive
    TNow = AC.time2datetime( [gmtime()] )[0]
    pstr = "Started moving files for {} to google drive @ {}"
    print(  pstr.format( dstr, TNow.strftime('%Y/%m/%d %H:%M') ) )
    # Move the files
    mv_plots2google_drive(dt=dt, debug=debug)
    # print that the job is finished.
    TNow = AC.time2datetime( [gmtime()] )[0]
    pstr = "Finished moving files for {} to google drive @ {}"
    print(  pstr.format( dstr, TNow.strftime('%Y/%m/%d %H:%M') ) )


def get_data_for_fcasts4datetimes(dts=None):
    """
    Download all forecast data in bulk for given datetimes (dts)
    """
    # Use the last 18 days, unless specific dates are provided
    if isinstance(dts, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
        # Use the last 18 days
        dts = [AC.add_days(dt, i*-1) for i in range(1, 18) ]
    # Loop by date and download data if available and not downloaded.
    for dt in dts:
        dstr = dt.strftime('%Y/%m/%d %H:%M')
        pstr = 'Checking historic GEOS5-CF forecasts are downloaded for {}'
        print(pstr.format(dstr))
        # Check whether files are available
        # TODO - Add check for GEOS5/CF for a specific date
        # Check if the files have been downloaded
        n_failed_doys = check4failed_downloads(dt=dt, retry=False,
                                               inc_GEOS5_Nv_collection=True)
        print(n_failed_doys)
        # If there files are not present, then download in series
        # TODO: Confirm why this call repeated - does this matter?
        n_failed_doys = check4failed_downloads(dt=dt, retry=True,
                                               inc_GEOS5_Nv_collection=True)
        print(n_failed_doys)


def get_expanded_variable_GMAO_files(verbose=True, debug=False):
    """
    Get the expanded variable NetCDF files from GMAO
    """
    from bs4 import BeautifulSoup
    import requests
    import re
    import wget
    # What is the root URL for the data?
    URL = 'https://gmao.gsfc.nasa.gov/gmaoftp/geoscf/ARNA/'
    # Access the html via requests
    r = requests.get(URL)
    html = r.text
    soup = BeautifulSoup(html, "html.parser")
    links = soup.findAll(href=re.compile("\.nc4$"))
    # Where to save the files?
    folder = get_local_folder('NASA_data')
    folder += 'GEOS_CF/ARNA/assim/expanded_variable_files/'
    AC.mk_folder(folder, verbose=verbose)
    # Check that the URL was correct
    assert len(links) > 0, 'WARNING: No links found @ URL ({})'.format(URL)
    # Get the one image
    for link in links:
        filename =  link['href']
        # Check it file present...
        if os.path.isfile(folder+filename):
            pstr = 'WARNING: Not downloading GMAO file as it exists! ({})'
            print(pstr.format(folder+filename) )
        # If not, download
        else:
            print( URL+filename )
            wget.download(URL+filename, folder+filename)


def mk_test_plots_from_GMAO_output(dt=None, load_all_dates=True):
    """
    Make a set of test plots for the expanded NASA GMAO plots
    """
    # Which day(s) to use for testing output?
    if isinstance(dt, type(None)) and (not load_all_dates):
    # Use the last 18 days, unless specific dates are provided
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
        # Setup the glob string for the files
        dstr = dt.strftime('%Y/%m/%d')
        gfstr = 'GEOS-CF.v01.rpl.ACSIS-ARNA.{}*.nc4'.format(dstr)

    elif load_all_dates:
        gfstr = 'GEOS-CF.v01.rpl.ACSIS-ARNA.*.nc4'
    else:
        # Setup the glob string for the files
        dstr = dt.strftime('%Y/%m/%d')
        gfstr = 'GEOS-CF.v01.rpl.ACSIS-ARNA.{}*.nc4'.format(dstr)
    # Where are the files?
    folder = get_local_folder('NASA_data')
    folder += 'GEOS_CF/ARNA/assim/expanded_variable_files/'
    # Load the NetCDF files
    files2use = glob.glob(folder+gfstr)
    ds = xr.open_mfdataset(files2use)
    # Update the units for lev
    HPa_l = get_GEOSCF_vertical_levels(native_levels=True)
    attrs = ds.lev.attrs.copy()
    attrs['units'] = 'hPa'
    attrs['standard_name'] = 'Pressure'
    attrs['long_name'] = 'Pressure'
    ds.lev.attrs = attrs
    ds.lev.values = [ HPa_l[int(i)] for i in ds.lev.values -1]

    # Update units
    ds = convert_GEOSCF_units(ds=ds, debug=True)
    # Now do a handful of quick plots
#    extra_str = 'GMAO_EXPANDED_TEST'
#    vars2plot = [i for i in ds.data_vars]
    vars2plot = ['O3', 'CO', 'Cl', 'U' , 'V', 'T', 'IO', 'BrO']
    title_date = 'avg 24/25th jan'
    folder = './'
#    for var2plot in vars2plot[:3]:
    for var2plot in vars2plot:
#        for lev2use in ds.lev.values:
#        for lev2use in [72, 51]:
        for lev2use in [525., 985.]:
#            for lev2use in ds.lev.values[:2]: # unhash if testing
            print( var2plot, lev2use, title_date)
#            ds_tmp = ds[[var2plot]].sel(lev=lev2use)
            bool1 = ds.lev.values == lev2use
            ds_tmp = ds[[var2plot]].isel(lev=bool1)
            ds_tmp.squeeze()
            #
            ds_tmp =  ds_tmp.mean(dim='time')
#            if not isinstance(extr_title_str, type(None)):
#                title += '\n '+ extr_title_str
            # Get the LateX for of the species name
            try:
                LaTeX_spec = AC.latex_spec_name(var2plot)
            except KeyError:
                pstr = 'WARNING: not converted {} to LaTeX form'
                print(pstr.format(var2plot))
                LaTeX_spec = var2plot
            # Force use of standard name as long name
            attrs = ds[var2plot].attrs
            attrs['long_name'] = LaTeX_spec
            ds_tmp[var2plot].attrs = attrs
            # Setup a string for the title
            title = '[{}] @ level {:0>2} on {}'
            title = title.format(LaTeX_spec, lev2use, title_date)
            # Set extra string for filename
            extra_str = 'ARNA_lev_{:0>2}_{}'.format(lev2use,
                                                    'GMAO_EXPANDED_TEST')
            # Now plot
            quick_map_plt_CV_1layer(ds_tmp, var2plot=var2plot,
                                             use_local_CVAO_area=True,
                                             extra_str=extra_str,
                                             extend='both',
                                             title=title,
                                             folder=folder,
                                             save_plot=True )

            # Do some garbage collecting
            gc.collect()


def mk_missing_ARNA_plots4dts(dts=None, mk_plots=True):
    """
    Make plots for specific datetimes
    """
    # Use the last 18 days, unless specific dates are provided
    if isinstance(dts, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
        # Use the last 18 days
        dts = [AC.add_days(dt, i*-1) for i in range(1, 18) ]
    # Loop a call to checking and plotting for a given date
    for dt in dts:
        mk_missing_ARNA_plots4dt(dt=dt, mk_plots=mk_plots)


def mk_missing_ARNA_plots4dt(dt=None, mk_plots=True):
    """
    Make plots for a specific datetime
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Check if the folders are not complete
    dstr = dt.strftime('%Y/%m/%d %H:%M')
    print('Checking plots have been made for {}'.format(dstr))
    # Check if the plots are already present
    folders2get = which_plot_folders_are_not_complete4dt(dt=dt)
    # Now only download those that are not complete
    if len(folders2get) == 0 :
        print('All plots present and correct')
    else:
        print('WARNING: None/not all plots present! - ', folders2get)
        # Hardwire folder names to check
        subfolders2check = {
        'alt_slice': 0,
        'alt_slice.zoomed': 1,
        'lat_slice' : 2,
        'lon_slice': 3,
        }
        # Regrid the files...
        folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5')
        regrid_GEOS5_files_in_folder(folder=folder)
        # Make the folder structure for the plots
        mk_folder_structure4fcast(dt=dt)
        # Check sub folder by folder and do call the specific fucntion for given case
        for subfolder in subfolders2check.keys():
            if (subfolder in folders2get) and  mk_plots:
                print("Making '{}' plots for '{}'".format(subfolder,dstr ))
                plot_type = subfolders2check[ subfolder ]
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
            print("Making '{}' plots for '{}'".format(subfolder,dstr ))
            regrid_files_then_do_ind_plots(dt=dt)

        # Also check the GMAO folders
        # NOTE: this will not work, as the files on the GMAO site are latest ones
        # A new function could be written to GMAO 2D plots, but not datagrams
#        subfolder = 'plots.GMAO'
#        if subfolder in folders2get:
#            print("Making '{}' plots for '{}'".format(subfolder,dstr ))
#            GMAO
#            get_latest_GEOS5_diagnostics(dt=dt)


def which_plot_folders_are_not_complete4dt(dt):
    """
    Check which plots have been made for a specific datetime (dt)
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Get the root plot folder
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5', inc_collection=False)
    folder += '/plots/'
    # Hardwire folder names to check
    subfolders2check =[
    'plots.GMAO',
    'alt_slice',
    'lon_slice',
    'lat_slice',
    'alt_slice.zoomed',
    'alt_slice.individual'
    ]
    if not os.path.isdir(folder):
        return subfolders2check
    else:
        # Setup a dataframe to store values
        df = pd.DataFrame()
        # Loop subfolders and get info on folders
        for subfolder in subfolders2check:
            # Evaluate full folder name
            subfolder_full = '{}/{}'.format(folder, subfolder)
            if os.path.isdir(subfolder_full):
                d = AC.get_stats_on_files_in_folder_as_dict(folder=subfolder_full)
            else:
                d = {'#': 0}
            df[subfolder] = pd.Series(d)
        # - Now decide if all the file are present
    #    df = df.T
        # Hardwire number of expected files
        nfiles4ARNA = {
        'plots.GMAO': 169,
        'alt_slice': 246,
        'lon_slice': 246,
        'lat_slice': 246,
        'alt_slice.zoomed': 246,
        'alt_slice.individual': 1230,
    #    'core.plus_024H.2020_01_22_12_00': 61,
    #    'core.plus_024H.2020_01_22_12_00': 61,
    #    'core.plus_030H.2020_01_22_18_00': 61,
    #    'core.plus_048H.2020_01_23_12_00': 61,
    #    'core.plus_054H.2020_01_23_18_00': 61,
    #    'core.plus_072H.2020_01_24_12_00': 61,
    #    'core.plus_078H.2020_01_24_18_00': 55, # Note this should be 61!
        }
        # Check the right number of files are present
        correct_num_col = 'Correct #'
        for col in df.columns:
            correct_num = df[col]['#'] == nfiles4ARNA[col]
            df.loc[ correct_num_col , col ] = correct_num
        # just select the name of the
        names = df.loc[ :, df.T[correct_num_col].values == False ].columns
        return list(names)




def copy_files2core_plot_folders(dt=None, verbose=True, debug=False):
    """
    Copy the files to core plot folders (+24, +48) for a given date
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Setup a string for the date of the forecast
    dt0_Str = dt.strftime('%Y/%m/%d %H:%M')
    # Also setup datetime objects for the +24 fcast and +48 fcasts
    dt24 = AC.add_days(dt, 1)
    dt30 = AC.add_hrs(dt, 24+6)
    dt48 = AC.add_days(dt, 2)
    dt54 = AC.add_hrs(dt, 48+6)
    dt72 = AC.add_days(dt, 3)
    dt78 = AC.add_hrs(dt, 72+6)
    # Get the root plot folder
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5', inc_collection=False)
    folder += '/plots/'

    # - Copy the files for the +24 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt.year, dt.month, dt24.day, dt24.hour, dt24.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(24, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt24.year, dt24.month, dt24.day, dt24.hour, dt24.minute)
    # Get a list of files to copy
    files = glob.glob( '{}/*/*{}.png'.format(folder, dstr) )
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_024*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # print some information about the files found
    pstr = 'Moving {:>3} files to {} from forecast started on {}'
    print( pstr.format(len(files), fstr, dt0_Str) )
    # Loop by files and move
    for file in files:
        filename  = file.split('/')[-1]
        if debug:
            print(filename, file )
        os.popen('cp {} {}'.format(file, dfolder+filename))

    # - Copy the files for the +30 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt30.year, dt30.month, dt30.day, dt30.hour, dt30.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(30, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    files = glob.glob( '{}/*/*{}.png'.format(folder, dstr) )
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_030*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # Loop by files and move
    for file in files:
        filename  = file.split('/')[-1]
        if debug:
            print(filename, file )
        os.popen('cp {} {}'.format(file, dfolder+filename))

    # - Copy the files for the +48 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt48.year, dt48.month, dt48.day, dt48.hour, dt48.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(48, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    files = glob.glob( '{}/*/*{}.png'.format(folder, dstr) )
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_048*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # Loop by files and move
    for file in files:
        filename  = file.split('/')[-1]
        if debug:
            print(filename, file )
        os.popen('cp {} {}'.format(file, dfolder+filename))

    # - Copy the files for the +54 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt54.year, dt54.month, dt54.day, dt54.hour, dt54.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(54, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    files = glob.glob( '{}/*/*{}.png'.format(folder, dstr) )
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_054*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # Loop by files and move
    for file in files:
        filename  = file.split('/')[-1]
        if debug:
            print(filename, file )
        os.popen('cp {} {}'.format(file, dfolder+filename))

    # - Copy the files for the +72 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt72.year, dt72.month, dt72.day, dt72.hour, dt72.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(72, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    files = glob.glob( '{}/*/*{}.png'.format(folder, dstr) )
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_072*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # Loop by files and move
    for file in files:
        filename  = file.split('/')[-1]
        if debug:
            print(filename, file )
        os.popen('cp {} {}'.format(file, dfolder+filename))

    # - Copy the files for the +78 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt78.year, dt78.month, dt78.day, dt78.hour, dt78.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(78, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    files = glob.glob( '{}/*/*{}.png'.format(folder, dstr) )
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_078*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob( '{}/*/*{}*'.format(folder, gstr) )
    # Loop by files and move
    for file in files:
        filename  = file.split('/')[-1]
        if debug:
            print(filename, file )
        os.popen('cp {} {}'.format(file, dfolder+filename))


def regrid_files_then_do_core_plots(dt=None, plt_in_parallel=False):
    """
    Check files are present, then regrid GEOS5 files and plot core plots
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Check that the files are all present and correct...
    n_failed_doys = check4failed_downloads(dt=dt, retry=False)
    print(n_failed_doys)
    plot_anyway=True
    if (n_failed_doys == 0) or plot_anyway:
        # Regrid the files...
        folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5')
        regrid_GEOS5_files_in_folder(folder=folder)
        # Make the folder structure for the plots
        mk_folder_structure4fcast(dt=dt)
        # Plot the data by species in parallel or series
        if plt_in_parallel:
            # Parallelise over variables
            plot_type = {
            'alt_slice': 0,
            'alt_slice.zoomed': 1,
            'lat_slice' : 2,
            'lon_slice': 3,
            }
            plot_type.keys()
            pstr = 'Plotting variables in parallel across a pool of {}'
            print( pstr.format(len(plot_type)) )
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
        # Print a waring to screen if case failed.
        pstr = 'WARNING: NO REGRID./PLOTTING AS ALL FILES NOT PRESENT/CORRECT!'
        print(pstr)
        print(dt)


def regrid_files_then_do_ind_plots(dt=None, plt_in_parallel=False):
    """
    Check files are present, then regrid GEOS5 files and individual plots
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Check that the files are all present and correct...
    n_failed_doys = check4failed_downloads(dt=dt, retry=False)
    print(n_failed_doys)
    plot_anyway=True
    if (n_failed_doys == 0) or plot_anyway:
        # Regrid the files...
        folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5')
        regrid_GEOS5_files_in_folder(folder=folder)
        # Make the folder structure for the plots
        mk_folder_structure4fcast(dt=dt)
        # Plot the data by species in parallel or series
        if plt_in_parallel:
            # Parallelise over variables
            vars2use =  [ 'NOy', 'NO2', 'O3', 'CO', 'Dust']
            vars2use = [[i] for i in vars2use ] # make sure these are lists!
            pstr = 'Plotting variables in parallel across a pool of {}'
            print( pstr.format(len(vars2use)) )
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


def get_Google_drive_ID4ARNA_forecast_plots():
    """
    Retrieve the google drive folder ID for the ARNA fcast folder
    """
    # Hardcoded for now...
    # TODO: retrieve this based on a walk through directories
    GD_ID = '1u4TAwsefzfFJ5OSa0PwiMgF_sanNR1QI'
    return GD_ID


def mv_plots2google_drive(dt=None, verbose=True, debug=False):
    """
    Move plots to a google drive account in an automated manner
    """
    from pydrive.auth import GoogleAuth
    from pydrive.drive import GoogleDrive
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Make the datetime string to be used for folders
    dt_str = '{}_{:0>2}_{:0>2}_{:0>2}z'
    dt_str = dt_str.format(dt.year, dt.month, dt.day, dt.hour)
    # What are the root folders for the data/plots
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                     inc_collection=False)
    folder += '/plots/'
    # Set location to store data
    gauth = GoogleAuth()
    drive = GoogleDrive(gauth)
    # Authenticate the connections
    gauth.CommandLineAuth()
    # Get the ID for the core
    fcast_GD_folder_ID = get_Google_drive_ID4ARNA_forecast_plots()
    # Create a folder for the date (dt)
    folder_metadata = {
#        'title' : dt_str+'_earth0_PyDrive_sync',
        'title' : dt_str,
        #
        "parents":  [{"id": fcast_GD_folder_ID}],
        # The mimetype defines this new file as a folder
        'mimeType' : 'application/vnd.google-apps.folder'
    }
    GD_folder_root = drive.CreateFile(folder_metadata)
    GD_folder_root.Upload()
    # Get info on the folder
    root_folder_title = GD_folder_root['title']
    root_folder_id = GD_folder_root['id']
    # Loop through files and directories to upload
    for path, directories, files in os.walk(folder):
        # What is the subroot for the folder?
        subfolder = path.split( dt_str )[-1].split('/plots/')[-1]
        if verbose:
            print( path, directories, files)
        # Create the subfolder for plots
        metadata = {
        'title': subfolder,
        "parents":  [{"id": root_folder_id}],
        "mimeType": "application/vnd.google-apps.folder"
        }
        GD_subfolder = drive.CreateFile(metadata)
        GD_subfolder.Upload()
        GD_subfolder_title = GD_subfolder['title']
        GD_subfolder_id = GD_subfolder['id']
        if verbose:
            print( subfolder, GD_subfolder_title, GD_subfolder_id )
        # Loop and upload by file
        for file in files:
            if debug:
                print(file)
            # Metadata for file
            metadata = {
            'title': file,
            "parents":  [{"id": GD_subfolder_id}],
            }
            # Set a file space
            GD_file = drive.CreateFile(metadata)
            # point at the required file
            GD_file.SetContentFile('{}/{}'.format(path,file))
            GD_file.Upload() # Upload the file.


def mv_GEOS_data_from_earth2viking(dt=None, user=None):
    """
    Move the GEOS-5/CF data from earth0 to viking
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Which user to use?
    if isinstance(user, type(None)):
        import getpass
        user = getpass.getuser()
    # - Get the locations of the GEOS-5 data
    # Get the locations on
    host = 'viking'
    collection = 'inst3_3d_aer_Np'
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                     collection=collection,
                                     host=host )
    # Make the local folder if it is not present
    AC.mk_folder(folder=folder)
    # Get the remote folder
    host = 'earth0'
    remote_folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                            collection=collection,
                                            host=host )
    # Check that the remote folder is present?  and the files are correct?
    # TODO...
    # Call scp via os
    host = 'earth0.york.ac.uk'
    call_str = "scp -r {}@{}:{}/* {}".format(user, host, remote_folder, folder)
    os.system(call_str)

    # - Get the locations of the GEOS-5 data
    # Get the locations on
    host = 'viking'
    collection = 'chm_inst_1hr_g1440x721_p23'
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                     collection=collection,
                                     host=host )
    # Make the local folder if it is not present
    AC.mk_folder(folder=folder)
    # Get the remote folder
    host = 'earth0'
    remote_folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                            collection=collection,
                                            host=host )
    # Check that the remote folder is present?  and the files are correct?
    # TODO...
    # Call scp via os
    host = 'earth0.york.ac.uk'
    call_str = "scp {}@{}:{}/* {}".format(user, host, remote_folder, folder)
    os.system(call_str)


def mv_plots2webfiles(dt=None, debug=True):
    """
    Move plots to York's 'webfiles' server in an automated manner
    """
    import paramiko
    import sys
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Make the datetime string to be used for folders
    dt_str = '{}_{:0>2}_{:0>2}_{:0>2}z'
    dt_str = dt_str.format(dt.year, dt.month, dt.day, dt.hour)
    # What are the root folders for the data/plots
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                     inc_collection=False)
    folder += '/plots/'
    # Define data location and ports etc
    host = "webfiles.york.ac.uk" # hard-coded
    port = 22
    transport = paramiko.Transport((host, port))
    password = "niangw?ndi" # hard-coded
    username = "chem631" # hard-coded
    # Connect to webfiles
    transport.connect(username=username, password=password)
    # Open a sftp link
    sftp = paramiko.SFTPClient.from_transport(transport)
    # Set location to store data
    dt_root_path = './SHERWEN_TOMAS/ARNA/fcast/{}/'.format(dt_str)
    # Create a folder for the date
    mkdir_if_not_present(sftp=sftp, path=dt_root_path)
    # Loop through files and directories to upload
    for path, directories, files in os.walk(folder):
        # What is the subroot for the folder?
#        dt_sub_path = dt_root_path + path.split( dt_str+'_test' )[-1] # TESTING!
        subfolder = path.split( dt_str )[-1].split('/plots/')[-1]
        dt_sub_path = '{}/{}/'.format(dt_root_path, subfolder)
        if debug:
            print( path, directories, files)
            print('')
            print(dt_sub_path)
        # open connection
        sftp = paramiko.SFTPClient.from_transport(transport)
        # Create the path if it doesn't exist
        mkdir_if_not_present(sftp=sftp, path=dt_sub_path)
        # Close connection
        sftp.close()
        # Loop and upload by file
        for file in files:
            # open connection
            sftp = paramiko.SFTPClient.from_transport(transport)
            # tranfer file
            sftp.put('{}/{}'.format(path,file), dt_sub_path+file)
            # Close connection
            sftp.close()


def mkdir_if_not_present(sftp=None, path=None, verbose=True):
    """
    Make a remote path via sftp if not present
    """
    try:
        sftp.chdir(path)  # Test if path exists
        if verbose:
            pstr = 'WARNING: remote path exists, so folder not created - {}'
            print( pstr.format(path) )
    except IOError:
        sftp.mkdir(path)  # Create path
        sftp.chdir(path)
        if verbose:
            print('CREATED remote path that did not exist - {}'.format(path) )


def mk_folder_structure4fcast(dt=None, mode='fcast', verbose=False):
    """
    Make folders required for a date forecast data & plots 4 ARNA
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Check that the forecast date is for the noon forecast
    assert dt.hour == 12, 'WARNING forecast date not starting at noon!'
    # Make the datetime string to be used for folders
    dt_str = '{}_{:0>2}_{:0>2}_{:0>2}z'
    dt_str = dt_str.format(dt.year, dt.month, dt.day, dt.hour)
    # Also setup datetime objects for the +24 fcast and +48 fcasts
    dt24 = AC.add_days(dt, 1)
    dt30 = AC.add_hrs(dt, 24+6)
    dt48 = AC.add_days(dt, 2)
    dt54 = AC.add_hrs(dt, 48+6)
    dt72 = AC.add_days(dt, 3)
    dt78 = AC.add_hrs(dt, 72+6)
    # What are the root folders for the data/plots
    NASA_data = get_local_folder('NASA_data')
    G5_folder = NASA_data + 'GEOS_5/ARNA/'
    GCF_folder = NASA_data + 'GEOS_CF/ARNA/'
    # Setup a list to populate with folders to make
    folders = []
    # - Setup folder structure for GEOS-CF structure
    # Setup folder for the data
    folders += ['{}/{}/{}'.format(GCF_folder, mode, dt_str)]
    # - Setup folder structure for GEOS-5 structure
    # Setup folder for the data
    folders += ['{}/{}/{}'.format(G5_folder, mode, dt_str)]
    # Setup folder for the plots
    ex_str = 'plots'
    plots_root = '{}/{}/{}/{}'.format(G5_folder, mode, dt_str, ex_str)
    folders += [plots_root]
    # those directly snapped from GMAO
    ex_str = 'plots.GMAO'
    folders += ['{}/{}'.format(plots_root, ex_str)]
    # Make strings for folders for +24H and +48H
    dstr = 'core.plus_{:0>3}H.{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr24 = dstr.format(24, dt24.year, dt24.month, dt24.day, dt24.hour, dt24.minute)
    dstr30 = dstr.format(30, dt30.year, dt30.month, dt30.day, dt30.hour, dt30.minute)
    dstr48 = dstr.format(48, dt48.year, dt48.month, dt48.day, dt48.hour, dt48.minute)
    dstr54 = dstr.format(54, dt54.year, dt54.month, dt54.day, dt54.hour, dt54.minute)
    dstr72 = dstr.format(72, dt72.year, dt72.month, dt72.day, dt72.hour, dt72.minute)
    dstr78 = dstr.format(78, dt78.year, dt78.month, dt78.day, dt78.hour, dt78.minute)
    # List of subfolders to make
    subfolders = [
    'alt_slice', 'lon_slice', 'lat_slice', 'alt_slice.zoomed',
    'alt_slice.individual', dstr24, dstr30, dstr48, dstr54, dstr72, dstr78,
    ]
    folders += ['{}/{}'.format(plots_root, i) for i in subfolders]
    # - Now loop and make folders
    for folder2mk in folders:
        AC.mk_folder(folder=folder2mk, verbose=verbose)


def mk_ind_plsts4fcast_GEOSCF_GEOS5(vars2use, dt=None,
                                    only_plot_where_GEOS5=True):
    """
    Make alt slice plots on a individual species basis from GEOS-CF output
    """
    # Get most recent fcast data as a dataset
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Use default variables if a list not provided
    if isinstance(vars2use, type(None)):
        vars2use =  [ 'NOy', 'NO2', 'O3', 'CO', 'Dust']
    # Get GEOS-CF folder for datettime
    collection = 'chm_inst_1hr_g1440x721_p23'
    G5_data_4dt = get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                          collection=collection)
    # Get GEOS-5 folder for datetime
    collection = 'inst3_3d_aer_Np'
    G5_data_4dt = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                          collection=collection)
    # Get plots folder
    G5_plots_4dt = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                           inc_collection=False)
    G5_plots_4dt += '/plots/'
    # - Get GEOS-CF data
    ds = get_most_recent_GEOSCF_data(dt=dt)
    t0_CF = ds.time.values[0]
    t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
    # - Get GEOS-5 data
    fstr = '{}/ARNA*_{}_*REGRID*.nc'
    files2use = glob.glob(fstr.format( G5_data_4dt, dt.year ) )
    files2use = list(sorted(files2use))
    # Open all the files as a single dataset
    ds5 = xr.open_mfdataset(files2use)
    t0_G5 = ds5.time.values[0]
    t0_G5_str = AC.dt64_2_dt([t0_G5])[0].strftime('%Y/%m/%d %H:%M')
    print(ds5.time, t0_G5_str)
    # Limit alt to certain levels
    hPa_heights = [1000, 900, 800, 700, 600, 500]
    ds5 = ds5.isel(lev=[(i in hPa_heights) for i in ds5.lev] )
    # Add GEOS5 Dust to the GEOS-CF dataset
    ds['Dust'] = ds5['Dust']
    ds['Dust'].attrs = ds5['Dust'].attrs
    # Only plot for points where both dust and noy data exist
    if only_plot_where_GEOS5:
        ds = ds.isel(time=[i in ds5.time.values for i in  ds.time.values ] )
    # Add a string for when the GEOS runs were initiated
    assert t0_CF_str == t0_G5_str
    extr_title_str = ' (GEOS-CF/5 from {})'.format( t0_G5_str )
    # Folder to save plots
    folder = G5_plots_4dt + '/alt_slice.individual/'
    # Make single species plots for GEOSCF and GEOS5
    plot_individual_spec_alt_slices(ds, vars2use=vars2use, folder=folder,
                                    extr_title_str=extr_title_str)


def convert_aircraft_locs2table():
    """
    Make a csv file with details on the airports linked to ARNA campaign
    """
    locs2use = ['Dakar', 'DSS', 'Sao Vicente Airport', 'VXE', 'Praia Airport', 'RAI', 'Gran Canaria Airport', 'LPA', 'Lisbon Airport', 'LIS', 'Paris (Charles de Gaulle) Airport', 'CDG']
    # Loop by location
    d = {}
    for loc in locs2use:
        lon, lat, alt = AC.get_loc(loc)
        # Add to dictionary
        d[loc] = {'Longitude':lon, 'Latitude': lat, 'Altitude' : alt}
    # Compile to dataframe and then save.
    pd.DataFrame(d).T.round(2).to_csv('ARNA_Airport_locs.csv')


def plot_up_longitudinally_sampled_locs(ds, var2use='noy', extents=None,
                                       add_detailed_map=True,
                                       add_flyable_range_as_circle=True,
                                       add_ARNA_locs=True,
                                       ):
    """
    plot up locations that data is sliced by lon
    """
    # Use example data
    if isinstance(ds, type(None)):
        folder = '/Users/tomassherwen/Google_Drive/Data/ARNA/GEOS_CF/'
        folder += '/data_GEOSCF_2019_12_14/'
        filename = 'ARNA_GEOSCF_chm_inst_1hr_g1440x721_p23_Cape_Verde_2019_353_noy_'
        filename += 'lvls_1000_900_800_700_600_500.nc'
        ds = xr.open_dataset( folder + filename )

    # Local area analysed as Cape Verde
    x0 = -30
    x1 =-10
    y0 = 0
    y1 = 25
    # - Select the data
    # Just get an example dataset
    ds = ds[[var2use]]
    # Select a single level and time
    ds = ds.sel(time=ds.time[0])
    ds = ds.sel(lev=ds.lev[0])
    # Set values region
    bool1 = ((ds.lon >= x0) & (ds.lon <= x1)).values
    bool2 = ((ds.lat >= y0) & (ds.lat <= y1)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    # Set all values to 1
    arr = ds[var2use].values
    arr[:] = np.NaN
    ds[var2use].values = arr
    # Plot the data
    projection = ccrs.PlateCarree()
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection=projection, aspect='auto', alpha=0.5)
    # Plot up the dummy data
    LatVar = 'lat'
    LonVar = 'lon'
    ds[var2use].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                             transform=ccrs.PlateCarree())
    # Now plot as a linear ring
    lats = (5, 35, 35, 5)
    lons2plot = [-18 ,-19.5, -21, -22.5, -24, -25.5]
    for lon in lons2plot:
        # Set lats
        lons = (lon, lon+0.125, lon+0.125, lon)
        #
        ring = LinearRing(list(zip(lons, lats)))
        ax.add_geometries([ring], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='green',
                          zorder=10, linestyle='-',
                          )
    # Add some grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=.5, color='gray', alpha=0.25, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    # Mark known places to help geo-locate viewers
    if add_ARNA_locs:
        locs2plot  = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
        for loc2plot in locs2plot:
            lon, lat, alt = AC.get_loc(loc2plot)
            # Now plot up locations
            ax.plot(lon, lat, 'bo', markersize=5, markerfacecolor='none',
                    markeredgewidth=2,
                    zorder=10,
                    markeredgecolor='black',
                    transform=ccrs.PlateCarree())
            # Add a label for the location?
#            ax.text(lon, lat+0.25, loc2plot, transform=ccrs.PlateCarree())
    # Add the flyable range of the FAAM BAE146
    if add_flyable_range_as_circle:
#        n_points = 1000
        # Approximate from James' max distance
        # ( 16.8331-13 ) *110667.45
        locs4circles = 'Dakar', 'Sao Vicente Airport',
        for loc in locs4circles:
            # Get locations to centre circle on
            lon, lat, alt = AC.get_loc(loc)
            # Radius in degrees
#            radius = 16.8331-13
            radius = 21 - 16.8331
            # Plot up circle
            ax.add_patch(mpatches.Circle(xy=[lon, lat],
                                         radius=radius,
                                         transform=projection,
                                         facecolor='none',
                                         edgecolor='grey',
                                         linestyle=':',
                                         zorder=10
                                         ))
    # Get limits of plotting data
    if isinstance(extents, type(None)):
        x0 = float(ds[LonVar].min())
        x1 = float(ds[LonVar].max())
        y0 = float(ds[LatVar].min())
        y1 = float(ds[LatVar].max())
        extents = (x0, x1, y0, y1)
    ax.set_extent(extents, crs=ccrs.PlateCarree())
    # Beautify the figure/plot
    add_detailed_map = True
    if add_detailed_map:
        # Add borders to map
        ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                       facecolor='none', zorder=50)
        # Also add minor islands (inc. Cape Verde)
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor=None,
                                                facecolor='none')
        ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)
    # Save the plot to disk
    savename = 'spatial_plot_Cape_Verde_plotted_data_lons'
    savename = AC.rm_spaces_and_chars_from_str(savename)
    plt.savefig(savename+'.png', dpi=dpi)


def plot_up_latitudinally_sampled_locs(ds, var2use='noy', extents=None,
                                       add_detailed_map=True,
                                       add_flyable_range_as_circle=True,
                                       add_ARNA_locs=True,
                                       ):
    """
    plot up locations that data is sliced by lon
    """
    # Use example data
    if isinstance(ds, type(None)):
        folder = '/Users/tomassherwen/Google_Drive/Data/ARNA/GEOS_CF/'
        folder += '/data_GEOSCF_2019_12_14/'
        filename = 'ARNA_GEOSCF_chm_inst_1hr_g1440x721_p23_Cape_Verde_2019_353_noy_'
        filename += 'lvls_1000_900_800_700_600_500.nc'
        ds = xr.open_dataset( folder + filename )
    # Local area analysed as Cape Verde
    x0 = -30
    x1 =-10
    y0 = 0
    y1 = 25
    # - Select the data
    # Just get an example dataset
    ds = ds[[var2use]]
    # Select a single level and time
    ds = ds.sel(time=ds.time[0])
    ds = ds.sel(lev=ds.lev[0])
    # Set values region
    bool1 = ((ds.lon >= x0) & (ds.lon <= x1)).values
    bool2 = ((ds.lat >= y0) & (ds.lat <= y1)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    # Set all values to 1
    arr = ds[var2use].values
    arr[:] = np.NaN
    ds[var2use].values = arr
    # Plot the data
    projection = ccrs.PlateCarree()
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection=projection, aspect='auto', alpha=0.5)
    # Plot up the dummy data
    LatVar = 'lat'
    LonVar = 'lon'
    ds[var2use].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                             transform=ccrs.PlateCarree())
    # Now plot as a linear ring
    lons = (-15, -30, -30, -15)
    lats2plot = [12, 13, 14, 15, 16, 17]
    for lat in lats2plot:
        # Set lats
        lats = (lat, lat+0.125, lat+0.125, lat)
        #
        ring = LinearRing(list(zip(lons, lats)))
        ax.add_geometries([ring], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='green',
                          zorder=10, linestyle='-',
                          )
    # Add some grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=.5, color='gray', alpha=0.25, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    # Mark known places to help geo-locate viewers
    if add_ARNA_locs:
        locs2plot  = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
        for loc2plot in locs2plot:
            lon, lat, alt = AC.get_loc(loc2plot)
            # Now plot up locations
            ax.plot(lon, lat, 'bo', markersize=5, markerfacecolor='none',
                    markeredgewidth=2,
                    zorder=10,
                    markeredgecolor='black',
                    transform=ccrs.PlateCarree())
            # Add a label for the location?
#            ax.text(lon, lat+0.25, loc2plot, transform=ccrs.PlateCarree())
    # Add the flyable range of the FAAM BAE146
    if add_flyable_range_as_circle:
        # Approximate from James' max distance
        # ( 16.8331-13 ) *110667.45
        locs4circles = 'Dakar', 'Sao Vicente Airport',
        for loc in locs4circles:
            # Get locations to centre circle on
            lon, lat, alt = AC.get_loc(loc)
            # Radius in degrees
#            radius = 16.8331-13
            radius = 21 - 16.8331
            # Plot up circle
            ax.add_patch(mpatches.Circle(xy=[lon, lat],
                                         radius=radius,
                                         transform=projection,
                                         facecolor='none',
                                         edgecolor='grey',
                                         linestyle=':',
                                         zorder=10
                                         ))
    # Get limits of plotting data
    if isinstance(extents, type(None)):
        x0 = float(ds[LonVar].min())
        x1 = float(ds[LonVar].max())
        y0 = float(ds[LatVar].min())
        y1 = float(ds[LatVar].max())
        extents = (x0, x1, y0, y1)
    ax.set_extent(extents, crs=ccrs.PlateCarree())
    # Beautify the figure/plot
    if add_detailed_map:
        # Add borders to map
        ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                       facecolor='none', zorder=50)
        # Also add minor islands (inc. Cape Verde)
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor=None,
                                                facecolor='none')
        ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)
    # Save plot
    savename = 'spatial_plot_Cape_Verde_plotted_data_lons'
    savename = AC.rm_spaces_and_chars_from_str(savename)
    plt.savefig(savename+'.png', dpi=dpi)


def plot_individual_spec_alt_slices(ds, folder='./',
                                    extr_title_str=None,
                                    vars2use=None):
    """
    Plot up other core species individually
    """
    # Use all variables in dataset if
    if isinstance(vars2use, type(None)):
        vars2use = list(ds.data_vars)
    # Loop by variables  and plot for each time and level
    for var2plot in vars2use:
        # Get the LateX for of the species name
        try:
            LaTeX_spec = AC.latex_spec_name(var2plot)
        except KeyError:
            print('WARNING: not converted {} to LaTeX form'.format(var2plot))
            LaTeX_spec = var2plot
        # Then loop by times in output
        dts = AC.dt64_2_dt( ds.time.values )
#        for n_time, t in enumerate( dts[:3] ): # unhash if testing
        for n_time, t in enumerate( dts ):
            # Sub select for time
            ds_tmp = ds[[var2plot]].sel(time=ds.time.values[n_time])
            # Create a date string to use to save file
            date_str = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
            date_str = date_str.format(t.year, t.month, t.day, t.hour,
                                       t.minute)
            title_date = t.strftime('%Y/%m/%d %H:%M')
            for lev2use in ds.lev.values:
#            for lev2use in ds.lev.values[:2]: # unhash if testing
                print( var2plot, lev2use, title_date)
                ds_tmpII = ds_tmp[[var2plot]].sel(lev=lev2use)
                # Setup a string for the title
                title = '[{}] @ {:.0f}hPa on {}'
                title = title.format(LaTeX_spec, lev2use, title_date)
                if not isinstance(extr_title_str, type(None)):
                    title += '\n '+ extr_title_str
                # Force use of standard name as long name
                attrs = ds_tmpII[var2plot].attrs
                attrs['long_name'] = LaTeX_spec
                ds_tmpII[var2plot].attrs = attrs
                # Set extra string for filename
                extra_str = 'ARNA_lev_{:.0f}_hPa_{}'.format(lev2use, date_str)
                # now plot
                quick_map_plt_CV_1layer(ds_tmpII, var2plot=var2plot,
                                                 use_local_CVAO_area=True,
                                                 extra_str=extra_str,
                                                 extend='both',
                                                 title=title,
                                                 folder=folder,
                                                 save_plot=True )
                # Do some clean up
                plt.close('all')
                del ds_tmpII
            del ds_tmp
            gc.collect()


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


def regrid_GEOS5_files_in_folder(folder=None, dt=None, doys2use=None,
                                 year=None,
                                 collection='inst3_3d_aer_Np',
                                 remake_files=False, debug=False ):
    """
    Regrid all GEOS-5 NetCDF files in folder to GEOS-CF format
    """
    #
    # Use the latest downloaded GEOS-5 folder if one is not set.
    if isinstance(folder, type(None)):
        # Use yesterday as the date if this is not provided
        if isinstance(dt, type(None)):
            # Just use yesterday for Now
            Tnow = AC.time2datetime( [gmtime()] )[0]
            # Get the 5-day forecast at noon...
            dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
            # Use yesterday
            dt =  AC.add_days(dt, -1)
        # Get the GEOS5 folder for a given date
        folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                         collection=collection)

    # Get current year if a year not specified
    if isinstance(year, type(None)):
        year = AC.time2datetime( [gmtime()] )[0].year
    # Which files are available
    files_in_folder = glob.glob(folder+'ARNA*_{}_*.nc'.format(year))
    # Make sure already regridded files are not used
    files_in_folder = [i for i in files_in_folder if 'REGRIDDED' not in i ]
    # Which doys to use
    if isinstance(doys2use, type(None)):
        doys2use = [i.split('/')[-1] for i in files_in_folder ]
        doys2use = [i.split(str(year))[-1][1:4] for i in doys2use ]
        doys2use = list(sorted((set(doys2use))))
    # Save name for the regridded file?
    writestr = 'ARNA_GEO5_Combined_Cape_Verde_{}_{:0>3}_REGRIDDED_ALL_lvls.nc'
    # Loop by doy and process
    for doy in doys2use:
        # Select files for a given doy
        doy_str = '_{:0>3}_'.format(doy)
        files2use = [i for i in files_in_folder if doy_str in i]
        if debug:
            print(len(files2use), [i.split('/')[-1] for i in files2use], '')
        # Set the name of the file to save
        filename2save = writestr.format( year, doy)
        if os.path.isfile(folder + filename2save) and not remake_files:
            if debug:
                pstr = 'WARNING: not regridding file as file already present'
                print(pstr)
        else:
            # Open all files as a single xarray dataset
            ds = xr.open_mfdataset(files2use)
            if debug:
                print([i for i in ds.time.values])
            # Convert files to a single file with variable for dust
            ds = convert_GEOS5_dust2_single_var_in_ug_m3(ds=ds,
                                                         rtn_just_new_var=True)
            # Regridd file to GEOS-CF template file
            # NOTE: this is hardcoded for now
            tplate_folder = get_local_folder('NASA_data')
            tplate_folder += 'GEOS_CF/ARNA/fcast/2019_11_06/data/'
            tplate_fname = 'ARNA_GEOSCF_chm_inst_1hr_g1440x721_p23_Cape_Verde_'
            tplate_fname += '2019_315_pm25du_rh35_gcc_lvls_1000_850_500.nc'
            # Regrid GEOS5 to GEOS-CF
            # TODO: add option to not vertically regrid
            ds = regrid_GEOS52GEOSCF_coordinates(ds, collection=collection,
                                                 tplate_folder=tplate_folder,
                                                 tplate_fname=tplate_fname)
            # Save as a new file
            ds.to_netcdf(folder+filename2save)


def regrid_GEOS52GEOSCF_coordinates(ds=None, collection='inst3_3d_aer_Np',
                                    tplate_folder=None, tplate_fname=None,
#                                    regrid_vertically=False,
                                    lat_var='lat', lon_var='lon',
                                    vars2regrid=None):
    """
    Regrid GEOS-5 output to be on the same grid as GEOS-CF
    """
    # Which variables to regrid?
    if isinstance(vars2regrid, type(None)):
        vars2regrid = list(ds.data_vars)

    # - Regrid horizontally (lat, lon)
    # Get lats and lons to regrid to
    dsT = xr.open_dataset( tplate_folder+tplate_fname )
    lat = dsT[lat_var]
    lon = dsT[lon_var]
    # Create a new dataset to template to from this
    ds_out = xr.Dataset({
        # 'time': ( ['time'], ds['time'] ),
        'lat': (['lat'], lat),
        'lon': (['lon'], lon),
    })
    # Create a regridder (to be reused )
    regridder = xe.Regridder(ds, ds_out, 'bilinear', reuse_weights=True)
    ds_l = []
    for var2use in vars2regrid:
        # Create a dataset to re-grid into
        ds_out = xr.Dataset({
#            'time': ( ['time'], ds['time'] ),
            'lat': (['lat'], lat),
            'lon': (['lon'], lon),
        })
        # Get a DataArray
        dr = ds[var2use]
        # Build regridder
        dr_out = regridder(dr)
        # Important note: Extra dimensions must be on the left, i.e. (time, lev, lat, lon) is correct but (lat, lon, time, lev) would not work. Most data sets should have (lat, lon) on the right (being the fastest changing dimension in the memory). If not, use DataArray.transpose or numpy.transpose to preprocess the data.
        # Exactly the same as input?
        xr.testing.assert_identical(dr_out['time'], ds['time'])
        # Save variable
        ds_l += [dr_out]

    # Combine variables
    dsN = xr.Dataset()
    for n, var2use in enumerate(vars2regrid):
        dsN[var2use] = ds_l[n]
        # use attributes of input dataset
        dsN[var2use].attrs = ds[var2use].attrs
    # Clean up
#    regridder.clean_weight_file()

    # - Regrid vertically (alt)
#    if regrid_vertically:
    if collection == 'inst3_3d_aer_Nv':
        # use bulk values for altitude
        lon, lat, alt = AC.get_latlonalt4res(res='4x5',
                                             full_vertical_grid=True)
        # NOTE: GEOS-5 output has its z axis flipped.
        dsN.lev.values = alt[::-1]
        # Which
        HPa_l = get_GEOSCF_vertical_levels( native_levels=True )
        hPa_as_km = [i for i in AC.hPa_to_Km(HPa_l) ]
        # Now interpolate
        dsN = dsN.interp(lev=hPa_as_km)
        # Put the units back as millibar
        dsN.lev.values = HPa_l
    # Return the dataset
    return dsN


def convert_GEOS5_dust2_single_var_in_ug_m3(ds=None, rtn_just_new_var=True):
    """
    Convert GEOS5 dust (in v/v) into ug/m3
    """
    # Name for combined data
    NewVar = 'Dust'
    # If more than 3 dust tracers, then sum up the dust tracers
    if len( [i for i in ds.data_vars if 'du' in i ] ) > 3:
        dust_tracers = ['du{:0>3}'.format(i) for i in range(1,6) ]
        ds[NewVar] = ds[dust_tracers[0]].copy()
        for tra in dust_tracers[1:]:
            ds[NewVar].values = ds[NewVar].values + ds[tra].values
    else:
        ds = ds.rename( {'du': NewVar} )
    # Now convert to ug/3m
    print(ds.time)
    data = ds[NewVar].values
    ds[NewVar].values = AC.convert_spec_v_v_2_ugm3( data=data, spec='DST1' )
    # Update the units
    attrs = ds[NewVar].attrs
    attrs['units'] = '$\mu$g m$^{-3}$'
    attrs['long_name'] = 'Dust ($\Sigma$ bins 1-5)'
    ds[NewVar].attrs = attrs
    if rtn_just_new_var:
        return ds[[NewVar]]
    else:
        return ds


def get_most_recent_GEOSCF_data(dt=None,
                                limit_lvls=True, limit_lons=False,
                                filestr='ARNA*{}*lvls_10*.nc',
                                collection='chm_inst_1hr_g1440x721_p23',
                                ):
    """
    Get most recent GEOS-CF output (downloaded and processed)
    """
    if isinstance(dt, type(None)):
        # Just use yesterday for Now
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # Get the date/location string for data
    filename = filestr.format( dt.year )
    folder2use = get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                         collection=collection)
    glob_str = '{}/{}'.format( folder2use, filename )
    files2use = glob.glob(glob_str)
    files2use = list(sorted(files2use))
    asstr = 'WARNING no files found matching {}'
    assert len(files2use), asstr.format(glob_str)
    # Open all the files as a single dataset
    ds = xr.open_mfdataset(files2use)
    print(ds.time)
    # Update the variable names to be in GEOS-Chem format
    ds = update_GEOSCF_names2GEOSChem_names(ds=ds)
    # Update the units to GEOS-Chem/interpretable ones
    ds = convert_GEOSCF_units(ds=ds)
    # Add extra analysis species
    ds = add_extra_derived_species(ds=ds)
    return ds


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
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt =  AC.add_days(dt, -1)
    # GEOS-5 folder
    collection = 'inst3_3d_aer_Np'
    G5_data_4dt = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                          collection=collection)
    # Root plot saving folder
    G5_plot_dir = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                          inc_collection=False)
    G5_plot_dir += '/plots/'

    # - Do analysis on a altitude slice basis
    if do_core_alt_analysis or (plot_type == 0):
        # Get GEOS-CF data
        ds = get_most_recent_GEOSCF_data(dt=dt)
        t0_CF = ds.time.values[0]
        t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
        # - Get GEOS-5 data
        fstr = '{}/ARNA*_{}_*REGRID*.nc'
        files2use = glob.glob(fstr.format( G5_data_4dt, dt.year ))
        files2use = list(sorted(files2use))
        # Open all the files as a single dataset
        ds5 = xr.open_mfdataset(files2use)#, combine='nested', concat_dim='time')
        t0_G5 = ds5.time.values[0]
        t0_G5_str = AC.dt64_2_dt([t0_G5])[0].strftime('%Y/%m/%d %H:%M')
        print(ds5.time, t0_G5_str)
        # - Process GEOS-CF/5 data
        # Limit alt to certain levels
        hPa_heights = [1000, 900, 800, 700, 600, 500]
        ds5 = ds5.isel(lev=[(i in hPa_heights) for i in ds5.lev] )
        # Add the GEOS-5 data to the GEOS-CF dataset
        ds['Dust'] = ds5['Dust']
        ds['Dust'].attrs = ds5['Dust'].attrs
        #  Only plot for points where both dust and noy data exist
        if only_plot_where_GEOS5:
            ds = ds.isel(time=[i in ds5.time.values for i in  ds.time.values ] )
        # Add a string for when the GEOS runs were initiated
        assert t0_CF_str == t0_G5_str
        extr_title_str = ' (GEOS-CF/5 from {})'.format( t0_G5_str )
        # Set save location for plots
        if testing_mode:
            folder = './' # unhash for testing
        else:
            folder = G5_plot_dir + '/alt_slice/'
        # Plot the two layer plot (setting plotting to use the Dust settings.)
        plot_spatial_concs_2layer(ds, verbose=True,
                                  # Plot NOy, then dust
#                                  var2plot1='NOy', var2plot2='Dust',
                                  # Plot Dust, the NOy
                                  var2plot1='Dust', var2plot2='NOy',
                                  region='Cape_Verde',folder=folder,
                                  testing_mode=testing_mode,
                                  extr_title_str=extr_title_str )

    if do_zoomed_alt_analysis or (plot_type == 1):
        # Get GEOS-CF data
        ds = get_most_recent_GEOSCF_data(dt=dt)
        t0_CF = ds.time.values[0]
        t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
        # - Get GEOS-5 data
        # Get data
        fstr = '{}/ARNA*_{}_*REGRID*.nc'
        files2use = glob.glob(fstr.format( G5_data_4dt, dt.year))
        files2use = list(sorted(files2use))
        # Open all the files as a single dataset
        ds5 = xr.open_mfdataset(files2use)#, combine='nested', concat_dim='time')
        t0_G5 = ds5.time.values[0]
        t0_G5_str = AC.dt64_2_dt([t0_G5])[0].strftime('%Y/%m/%d %H:%M')
        print(ds5.time, t0_G5_str)
        # Limit alt to certain levels
        hPa_heights = [1000, 900, 800, 700, 600, 500]
        ds5 = ds5.isel(lev=[(i in hPa_heights) for i in ds5.lev] )
        # Add to the GEOS-CF dataset
        ds['Dust'] = ds5['Dust']
        ds['Dust'].attrs = ds5['Dust'].attrs
        #  Only plot for points where both dust and noy data exist
        if only_plot_where_GEOS5:
            ds = ds.isel(time=[i in ds5.time.values for i in  ds.time.values ] )
        # AddNOy * Dust (GEOS5)
        ds['NOy*Dust'] = ds['Dust'] * ds['NOy']
        # Add a string for when the GEOS runs were initiated
#        extr_title_str = ' (G-CF:{}, G-5:{})'.format( t0_CF_str, t0_G5_str )
        assert t0_CF_str == t0_G5_str
        extr_title_str = ' (GEOS-CF/5 from {})'.format( t0_G5_str )
        # Also put the zoomed in plots and save to folder
        if testing_mode:
            folder = './' # unhash for testing
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
                                  extr_title_str=extr_title_str )

    # - Do analysis on a longitudinal slice basis
    if do_core_lon_analysis or (plot_type == 2):
        # - Get GEOS-CF data
        filestr = 'ARNA*{}*-18.0_*.nc'
        ds = get_most_recent_GEOSCF_data(dt=dt, filestr=filestr)
        t0_CF = ds.time.values[0]
        t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
        # - Get GEOS-5 data
        # Get data
        fstr = '{}/ARNA*_{}_*REGRID*.nc'
        files2use = glob.glob(fstr.format( G5_data_4dt, dt.year ) )
        files2use = list(sorted(files2use))
        # Open all the files as a single dataset
        ds5 = xr.open_mfdataset(files2use)#, combine='nested', concat_dim='time')
        t0_G5 = ds5.time.values[0]
        t0_G5_str = AC.dt64_2_dt([t0_G5])[0].strftime('%Y/%m/%d %H:%M')
        print(ds5.time, t0_G5_str)
        # Only consider the slices for lon
        lons2use = [-18 ,-19.5, -21, -22.5, -24, -25.5]
        ds5 = ds5.isel(lon=[(i in lons2use) for i in ds5.lon] )
        # Add to the GEOS-CF dataset
        ds['Dust'] = ds5['Dust']
        ds['Dust'].attrs = ds5['Dust'].attrs
        # Only plot for the GEOS5 time steps.
        if only_plot_where_GEOS5:
            ds = ds.isel(time=[i in ds5.time.values for i in  ds.time.values ] )
        # Attributes
        attrs = ds.lev.attrs
        attrs['units'] = 'millibar'
        ds.lev.attrs = attrs
        # Add a string for when the GEOS runs were initiated
        assert t0_CF_str == t0_G5_str
        extr_title_str = ' (GEOS-CF/5 from {})'.format( t0_G5_str )
        # Also put the lon slice in plots and save to folder
        if testing_mode:
            folder = './' # unhash for testing
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
                                               extr_title_str=extr_title_str )

    # - Do analysis on a latitude slice basis
    if do_core_lat_analysis or (plot_type == 3):
        # - Get GEOS-CF data
        filestr = 'ARNA*{}*ALL_lvls_lats_*.nc'
        ds = get_most_recent_GEOSCF_data(dt=dt, filestr=filestr)
        t0_CF = ds.time.values[0]
        t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
        # - Get GEOS-5 data
        # Get data
        fstr = '{}/ARNA*_{}_*REGRID*.nc'
        files2use = glob.glob(fstr.format( G5_data_4dt, dt.year ) )
        files2use = list(sorted(files2use))
        # Open all the files as a single dataset
        ds5 = xr.open_mfdataset(files2use)#, combine='nested', concat_dim='time')
        t0_G5 = ds5.time.values[0]
        t0_G5_str = AC.dt64_2_dt([t0_G5])[0].strftime('%Y/%m/%d %H:%M')
        print(ds5.time, t0_G5_str)
        # Only consider the slices for lon
        lats2use =  [12, 13, 14, 15, 16, 17]
        ds5 = ds5.isel(lat=[(i in lats2use) for i in ds5.lat] )
        # Add to the GEOS-CF dataset
        ds['Dust'] = ds5['Dust']
        ds['Dust'].attrs = ds5['Dust'].attrs
        # Only plot for the GEOS5 time steps.
        if only_plot_where_GEOS5:
            ds = ds.isel(time=[i in ds5.time.values for i in  ds.time.values ] )
        # Attributes
        attrs = ds.lev.attrs
        attrs['units'] = 'millibar'
        ds.lev.attrs = attrs
        # Add a string for when the GEOS runs were initiated
        assert t0_CF_str == t0_G5_str
        extr_title_str = ' (GEOS-CF/5 from {})'.format( t0_G5_str )
        # Also put the lat slice in plots and save to folder
        if testing_mode:
            folder = './' # unhash for testing
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
                                               extr_title_str=extr_title_str )


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


def get_GEOS_data_folder4dt(dt=None, product='GEOS_CF', host=None,
                            mode='fcast', collection='inst3_3d_aer_Np',
                            inc_collection=True):
    """
    Get the data folder location for a given GEOS product and datetime
    """
    # Where is the data for the product?
    NASA_data = get_local_folder('NASA_data', host=host)
    folder = '{}/{}/ARNA/'.format(NASA_data, product)
    # And for the datetime
    date_str = '{}_{:0>2}_{:0>2}_{:0>2}z/'
    date_str = date_str.format(dt.year, dt.month, dt.day, dt.hour)
    # Return the string including collection and production mode
    if inc_collection:
        return '{}/{}/{}/data.{}/'.format( folder, mode, date_str, collection )
    else:
        return '{}/{}/{}/'.format( folder, mode, date_str )


def get_latest_GEOSCF_fcast_data_ALL(dt=None, just_check_yesterday=True,
                                     rm_existing_file=False, debug=True):
    """
    Get the latest GEOS-CF forecast data (if new data available)
    """
    if isinstance(dt, type(None)):
        if just_check_yesterday:
            # Just use yesterday for Now
            TNow = AC.time2datetime( [gmtime()] )[0]
            dt = AC.add_days(TNow, -1)
            dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
        else:
            # What is the last data there is data locally for?
            last_local_data_date = what_is_latest_data_locally(debug=debug)
            assert_err = 'Stopping as no data for last 5 days'
            assert not isinstance(last_local_data_date, type(None)), assert_err
            # When was the last forecast?
            last_fcast_start = when_should_GEOSCF_have_last_run_from()
            # Use this as the datetime to download
            dt = last_fcast_start
            ass_str = 'The most recent data is already locally available'
            assert last_fcast_start > last_local_data_date, ass_str
            # Download data if available from nasa
#            lastest_GEOSCF_avail = check_if_latest_GEOSCF_is_available(last_fcast_start)
    else:
        # Using the provided datetime (dt)
        pass

    # - Retrieve the latest GEOS-CF data - sliced by altitude
    download_GEOSCF_fcast_data4date(dt=dt, limit_lvls=True,
                                    rm_existing_file=rm_existing_file,
                                    limit_lons=False)
    # - Retrieve the latest GEOS-CF data - sliced by longitude
    download_GEOSCF_fcast_data4date(dt=dt, limit_lvls=False,
                                    rm_existing_file=rm_existing_file,
                                    limit_lons=True)
    # - Retrieve the latest GEOS-CF data - sliced by latitude
    download_GEOSCF_fcast_data4date(dt=dt, limit_lvls=False,
                                    rm_existing_file=rm_existing_file,
                                    limit_lons=False, limit_lats=True)



def get_latest_GEOSCF_fcast_data(dt=None, just_check_yesterday=True,
                                 limit_lvls=False,
                                 vars2use=None, doys2use=None,
                                 retry_download_if_data_not_present=True,
                                 limit_lons=False, rm_existing_file=False,
                                 limit_lats=False, debug=True):
    """
    Get the latest GEOS-CF forecast data (if new data available)
    """
    if debug:
        pstr = 'limit_lons={}, limit_lvls={}, limit_lats={} - in {}'
        func = 'get_latest_GEOSCF_fcast_data'
        print( pstr.format( limit_lons, limit_lvls, limit_lats, func ) )
    # Set date to use if not provided
    if isinstance(dt, type(None)):
        if just_check_yesterday:
            # Just use yesterday for Now
            TNow = AC.time2datetime( [gmtime()] )[0]
            dt = AC.add_days(TNow, -1)
            dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
        else:
            # What is the last data there is data locally for?
            last_local_data_date = what_is_latest_data_locally(debug=debug)
            assert_err = 'Stopping as no data for last 5 days'
            assert not isinstance(last_local_data_date, type(None)), assert_err
            # When was the last forecast?
            last_fcast_start = when_should_GEOSCF_have_last_run_from()
            # Use this as the datetime to download
            dt = last_fcast_start
            ass_str = 'The most recent data is already locally available'
            assert last_fcast_start > last_local_data_date, ass_str
            # Download data if available from nasa
#            lastest_GEOSCF_avail = check_if_latest_GEOSCF_is_available(last_fcast_start)
    else:
        # Using the provided datetime (dt)
        pass
    # Setup a text string with for dt
    dstr = dt.strftime('%Y/%m/%d %H:%M')

    # -
    # Is the data available
    dt_is_lastest_run = is_dt_latest_GEOSCF(dt=dt)
    # If it is, just download it
    if dt_is_lastest_run:
        # - Retrieve the latest GEOS-CF data - sliced by alt, lat, or lon
        download_GEOSCF_fcast_data4date(dt=dt,
                                        rm_existing_file=rm_existing_file,
                                        vars2use=vars2use,
                                        doys2use=doys2use,
                                        limit_lvls=limit_lvls,
                                        limit_lons=limit_lons,
                                        limit_lats=limit_lats)
    elif retry_download_if_data_not_present:
        # Minutes to wait
        times2wait = 30
        min2wait = 4
        pstr = 'WARNING: data not availible ({}), so trying every {} min.({}x)'
        print( pstr.format(dstr, min2wait, times2wait) )
        for time2wait in range(times2wait):
            pstr = 'WARNING: Wait for {:>2} min @ {}, then retrying download'
            TNow = AC.time2datetime( [gmtime()] )[0]
            print( pstr.format(min2wait, TNow.strftime('%Y/%m/%d %H:%M') ) )
            time.sleep(min2wait*60)
            pstr = 'Finished waiting @ {}'
            print(pstr.format(TNow.strftime('%Y/%m/%d %H:%M')))
            # Check if the data is availible
            dt_is_lastest_run = is_dt_latest_GEOSCF(dt)
            if dt_is_lastest_run:
                print('Data now available! Attempting download.')
                # - Retrieve the latest GEOS-CF data - sliced by alt, lat, or lon
                download_GEOSCF_fcast_data4date(dt=dt,
                                             rm_existing_file=rm_existing_file,
                                                vars2use=vars2use,
                                                doys2use=doys2use,
                                                limit_lvls=limit_lvls,
                                                limit_lons=limit_lons,
                                                limit_lats=limit_lats)
                print('WARNING: Stopping. Completed attempt at data download.')
                sys.exit()
            else:
                print('Data still not available...')
    else:
        print('ERROR: Data not available for date so cannot attempt download')


def get_latest_GEOSCF_fcast_data_alt_slice(dt=None, debug=True):
    """
    Wrapper to get alt slices from the latest GEOS-CF forecast data
    """
    # Hardcode options here for now.
    limit_lvls=True
    limit_lons=False
    # Call the main fcast data retrieval function
    get_latest_GEOSCF_fcast_data(dt=dt,
                                 limit_lvls=limit_lvls,
                                 limit_lons=limit_lons,
#                                 limit_lats=limit_lats
                                 )


def check4GEOSCF_failed_downloads(dt=None, folder=None, n_doys=6,
                                  var2use='noy',
                                  limit_lons=False, limit_lvls=False,
                                  limit_lats=False,
                                  debug=False):
    """
    Check that all the GEOSCF files have downloaded correctly
    """
    import os
    import glob
    # Print
    if debug:
        pstr = 'limit_lons={}, limit_lvls={}, limit_lats={} - in {}'
        func = 'check4GEOSCF_failed_downloads'
        print( pstr.format( limit_lons, limit_lvls, limit_lats, func ) )
    # Sliced by longitude
    if (limit_lvls == True) and (limit_lons==False) and (limit_lats==False):
        ext_str = '*lvls_1000_900_800_700_600_500*'
    # Sliced by altitude
    if (limit_lvls==False) and (limit_lons==True) and (limit_lats==False):
        ext_str = '*_lons_*'
    # Sliced by latitude
    if (limit_lvls==False) and (limit_lons==False) and (limit_lats==True):
        ext_str = '*_lats_*'
    # Also only consider one variable at a fime
    file_str = '*{}*{}*'.format(var2use, ext_str)
    # - Work out what the expected files for a date are
    expected_dates = [AC.add_days(dt, i) for i in range(n_doys)]
    dfT = pd.DataFrame(index=expected_dates)
    expected_doys = list(dfT.index.dayofyear.values)
    # -  Now test for files
    # Get file(s)
    files = glob.glob(folder+file_str)
    files = list(sorted(files))
    # Check if any files were found?
    if len(files) > 0:
        # Get size(s) of files  and the doys they are for
        doys = [i.split(str(dt.year))[-1][1:4] for i in files ]
        sizes = [os.path.getsize(i) for i in files]
        # Compile into a dataframe
        df = pd.DataFrame({'doy': doys, 'sizes': sizes})
        if debug:
            print(df)
        # Get details on dimensions
        dims = []
        for file in files:
            try:
                dims2add = [dict(xr.open_dataset(file).dims)]
            except OSError:
                dims2add = [[]]
            # Add the dimensions to the save list
            if (len(dims2add[0]) != 4):
                dims2add = [{'lat': 0, 'lev': 0, 'lon': 0, 'time': 0}]
            dims += dims2add
        # Add the dimensions
        for n_file, file in enumerate(files):
            dims4file = dims[n_file]
            for key in dims4file.keys():
                df.loc[n_file,key ] = dims4file[key]
        # - Add checks that all the lat, lon and, lev dims are the same length?
        if debug:
            print(df)
        coords2test = ['lon', 'lat', 'lev']
        for coord in coords2test:
            mode_ = df[coord].mode().values[0]
            # Add a column with a test for the coordinates
            var2use = 'eq_{}_mode'.format(coord)
            df[var2use] = df[coord] == mode_
        # - Do tests on doys and save those that fail to a list
        failed_doys = []
        # - Check if all coords are modal
        # Check that all coords are the same
        cols2use = [i for i in df.columns if 'mode' in i ]
        if debug:
            print( df[cols2use].values.all() )
            print( '')
            print( df )
        if df[cols2use].values.all():
            pass
        else:
            for doy in df['doy'].values:
                tmp = df.loc[df['doy']==doy,:]
                if debug:
                    print(tmp)
                if tmp[cols2use].values.all():
                    pass
                else:
                    failed_doys += [int(tmp['doy'].values.astype(int))]

        # - Check if the expected number of timestamps are in each time
        expected_times = [12.0, 24.0, 24.0, 24.0, 24.0, 13.0][:df.shape[0]]
        var2use = 'expected_time'
        df[var2use] = df['time'].values == expected_times
        # Save a list of the doys for which the retrieval failed for
        tmp_df = df.loc[df['expected_time'] == False, :]['doy']
        failed_doys += list(tmp_df.values.astype(int))
        # - Check if all the expected doys are present
        retrieved_doys = [int(i) for i in df['doy'].values]
        failed_doys += [i for i in expected_doys if (i not in retrieved_doys)]
        # Print a warning if there are doys included that are not expected
        unexpected_doys = [i for i in retrieved_doys if i not in expected_doys]
        if len(unexpected_doys) > 0:
            pstr = 'WARNING: Unexpected doys have been downloaded -'
            print(pstr, unexpected_doys)
        # - Return a list of the failed doys
        if debug:
            print(df, failed_doys)
        return list(sorted(set(failed_doys)))
    else:
        pstr = "WARNING: No files for '{}'  - using expected doys (folder:{})"
        print(pstr.format(var2use, folder))
        return expected_doys


def check4GEOS5_failed_downloads(dt=None, folder=None, n_doys=6, var2use='du',
                                 debug=False):
    """
    Check that all the GEOS5 files have downloaded correctly
    """
    import os
    import glob
    # Get file(s), size(s), and the doys they are for
    # Also only consider one variable at a fime
    files = glob.glob(folder+'*_{}_*.nc*'.format(var2use))
    # Make sure only none regridded files are checked
    files = [i for i in files if 'REGRIDDED' not in i ]
    files = list(sorted(files))
    # - Work out what the expected files for a date are
    expected_dates = [AC.add_days(dt, i) for i in range(n_doys)]
    dfT = pd.DataFrame(index=expected_dates)
    expected_doys = list(dfT.index.dayofyear.values)
    # Check if any files were found?
    if len(files) > 0:
        # Get size(s) of files  and the doys they are for
        doys = [i.split(str(dt.year))[-1][1:4] for i in files ]
        sizes = [os.path.getsize(i) for i in files]
        # Compile into a dataframe
        df = pd.DataFrame({'doy': doys, 'sizes': sizes})
        # Get details on dimensions
        dims = []
        for file in files:
            try:
                dims2add = [dict(xr.open_dataset(file).dims)]
            except OSError:
                dims2add = [[]]
            # Add the dimensions to the save list
            if debug:
                print(len(dims2add), dims2add )
            if (len(dims2add[0]) != 4):
                dims2add = [{'lat': 0, 'lev': 0, 'lon': 0, 'time': 0}]
            dims += dims2add
        # Add the dimensions
        for n_file, file in enumerate(files):
            dims4file = dims[n_file]
            for key in dims4file.keys():
                df.loc[n_file,key ] = dims4file[key]
        # Do tests on other variables
    #    assert len(set(df['lat'])) == 1, 'WARNING: files do not have the same # of lats!'
    #    assert len(set(df['lon'])) == 1, 'WARNING: files do not have the same # of lats!'
    #    assert len(set(df['lev'])) == 1, 'WARNING: files do not have the same # of lats!'
        # - add checks that all the lat, lon and, lev dims are the same length?
        coords2test = ['lon', 'lat', 'lev']
        for coord in coords2test:
            mode_ = df[coord].mode().values[0]
            # Add a column with a test for the coordinates
            var2use = 'eq_{}_mode'.format(coord)
            df[var2use] = df[coord] == mode_
        # - Do tests on doys and save those that fail to a list
        failed_doys = []
        # - Check if all coords are modal
        # Check that all coords are the same
        cols2use = [i for i in df.columns if 'mode' in i ]
        if df[cols2use].values.all():
            pass
        else:
            for doy in df['doy'].values:
                tmp = df.loc[df['doy']==doy,:]
                if tmp[cols2use].values.all():
                    pass
                else:
                    failed_doys += [int(tmp['doy'].values)]

        # - Check if the expected number of timestamps are in each time
        expected_times = [4., 8., 8., 8., 8., 5.][:df.shape[0]]
        var2use = 'expected_time'
        df[var2use] = df['time'].values == expected_times
        # Save a list of the doys for which the retrieval failed for
        df_tmp = df.loc[df['expected_time'] == False, :]['doy']
        failed_doys += list( df_tmp.values.astype(int) )
        # - Check if all the expected doys are present
        expected_doys = list(dfT.index.dayofyear.values)
        retrieved_doys = [int(i) for i in df['doy'].values]
        failed_doys += [i for i in expected_doys if (i not in retrieved_doys)]
        # Print a warning if there are days
        unexpected_doys = [i for i in retrieved_doys if i not in expected_doys]
        if len(unexpected_doys) > 0:
            pstr = 'WARNING: Unexpected doys have been downloaded -'
            print(pstr, unexpected_doys)
        # - Return a list of the failed doys
        if debug:
            print(failed_doys)
        return list(sorted(set(failed_doys)))
    else:
        pstr = "WARNING: No files found for GEOS5 so using expected doys (folder:{})"
        print(pstr.format(folder))
        return expected_doys



def get_latest_GEOSCF_fcast_data_lon_slice(dt=None, debug=True):
    """
    Wrapper to get lon slices from the latest GEOS-CF forecast data
    """
    # Hardcode options here for now.
    limit_lvls = False
    limit_lons = True
    # Call retrieval of forecast data
    get_latest_GEOSCF_fcast_data(dt=dt,
                                 limit_lvls=limit_lvls,
                                 limit_lons=limit_lons,
#                                 limit_lats=limit_lats
                                 )


def get_latest_GEOSCF_fcast_data_lat_slice(dt=None, debug=True):
    """
    Wrapper to get lat slices from the latest GEOS-CF forecast data
    """
    # Hardcode options here for now.
    limit_lvls = False
    limit_lons = False
    limit_lats = True
    # Call retrieval of forecast data
    get_latest_GEOSCF_fcast_data(dt=dt,
                                 limit_lvls=limit_lvls,
                                 limit_lons=limit_lons,
                                 limit_lats=limit_lats
                                 )


def get_latest_GEOS5_fcast_data(dt=None, mode='fcast',
                                collection='inst3_3d_aer_Np', doys2use=None,
                                vars2use=None,
#                                collection='inst3_3d_aer_Nv',
#                                dt=None,
                                retry_download_if_data_not_present=True,
                                rm_existing_file=False,
                                debug=True):
    """
    Get the latest GEOS-5 forecast data (5-day started @ noon)
    """
    # Use yesterday at noon if no date provided
    if isinstance(dt, type(None)):
        TNow = AC.time2datetime( [gmtime()] )[0]
        dt = AC.add_days(TNow, -1)
        dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
    # Setup a text string with for dt
    dstr = dt.strftime('%Y/%m/%d %H:%M')
    # NOTE: this is current dealt with in get_GEOS5_as_ds_via_OPeNDAP (in fcast mode)
    try:
        ds = AC.get_GEOS5_as_ds_via_OPeNDAP(dt=dt, mode=mode,
                                            collection=collection )
        # Temporarily do not allow passing of date to OPeNDAP fetcher
#        ds = AC.get_GEOS5_as_ds_via_OPeNDAP(dt=None, mode=mode, collection=collection)
    # Check it has files...
    except OSError:
        print('GEOS5 - Failure loading data for ', dt )
        try:
            dt = AC.add_days(dt, -1)
            print('GEOS5 - Looking on the day before', dt )
            ds = AC.get_GEOS5_as_ds_via_OPeNDAP(dt=dt, mode=mode,
                                                collection=collection)
            # Temporarily do not allow passing of date to OPeNDAP fetcher
#            ds = AC.get_GEOS5_as_ds_via_OPeNDAP(dt=None, mode=mode, collection=collection)
        except:
            print( 'WARNING: FAILED to find GEOS5 latest dt... stopping now.')
    # Get initial date in forecast
    dt = AC.dt64_2_dt( [ds['time'].values[0]])[0]
    # What would the folder be called for this noon forecast?
    NASA_data = get_local_folder('NASA_data')
    folder = NASA_data + 'GEOS_5/ARNA/fcast/'
    folderstr = '{}/{}_{:0>2}_{:0>2}_{:0>2}z/data.{}/'
    folder = folderstr.format(folder, dt.year, dt.month, dt.day, dt.hour,
                              collection)
    # Check if the folder is present, if not create it
    if os.path.isdir(folder):
        print('WARNING: folder already exists')
    else:
        os.makedirs(folder)
        print('Created folder for data ({})'.format(folder))

    # Set the variables to extract
    if isinstance(vars2use, type(None)):
        # Make the default variable du for inst3_3d_aer_Np (only one variabel regardless)
        if collection == 'inst3_3d_aer_Np':
            vars2use = ['du']
    # Check if the requested date is present.
    dt_is_lastest_run = is_dt_latest_GEOS5(dt=dt, collection=collection,
                                           mode=mode)
    override_GEOS5_check = True
    if dt_is_lastest_run or override_GEOS5_check:
        # Download this data if not already present
        get_GEOS5_data_cubes4collection(ds=ds, mode=mode,
                                        collection=collection,
                                        folder=folder,
                                        vars2use=vars2use, doys2use=doys2use,
                                        rm_existing_file=rm_existing_file,
                                        )

    elif retry_download_if_data_not_present:
        # Minutes to wait
        times2wait = 5
        min2wait = 4
        pstr = 'WARNING: data not availible ({}), so trying again every {} minutes ({}x)'
        print( pstr.format(dstr, min2wait, times2wait) )
        for time2wait in range(times2wait):
            pstr = 'WARNING: Waiting for {:>2} min @ {}, then re-attempting download'
            TNow = AC.time2datetime( [gmtime()] )[0]
            print( pstr.format(min2wait, TNow.strftime('%Y/%m/%d %H:%M') ) )
            time.sleep(min2wait*60)
            print('Finished waiting @ {}'.format(TNow.strftime('%Y/%m/%d %H:%M')))
            # Check if the data is availible
            dt_is_lastest_run = is_dt_latest_GEOS5(dt)
            if dt_is_lastest_run:
                print('Data now available! Attempting download.')
                # Download this data if not already present
                get_GEOS5_data_cubes4collection(ds=ds, mode=mode,
                                                collection=collection,
                                                folder=folder,
                                                vars2use=vars2use,
                                                doys2use=doys2use,
                                                rm_existing_file=True,
                                                )

                print('WARNING: Stopping. Completed attempt at data download')
                sys.exit()
            else:
                print('Data still not available...')
    else:
        print('ERROR: Data not available for date so cannot attempt download')

    # Force a garbage clean up
    gc.collect()
    # - Regrid the files if all are present
    print( 'Checking if files sucessfully downloaded. If so, regridding them!')
    #
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                     collection=collection)
    failed_doys = check4GEOS5_failed_downloads(dt=dt, folder=folder)
    # If there are no failed days, then regrid
    if len(failed_doys) == 0:
        print( 'Attempting to regrid GEOS5 files')
        regrid_GEOS5_files_in_folder(folder=folder)
        print( 'Regridded GEOS5 files')
    else:
        print('Not all files present so, ')


def check4failed_downloads(dt=None, retry=True, n_doys=6,
                           inc_GEOS5_Nv_collection=False,
                           debug=False):
    """
    Check for failed downloads of GEOS-5/GEOS-CF files, then retry getting these
    """
    # Just use noon yesterday if no date tis provided
    if isinstance(dt, type(None)):
        TNow = AC.time2datetime( [gmtime()] )[0]
        dt = AC.add_days(TNow, -1)
        dt = datetime.datetime(dt.year, dt.month, dt.day, 12)

    # Local variables
    dt_fmt = dt.strftime('%Y/%m/%d %H:%M')
    sum_n_failed_doys = 0
    # What are the expected doys?
    expected_dates = [AC.add_days(dt, i) for i in range(n_doys)]
    dfT = pd.DataFrame(index=expected_dates)
    expected_doys = list(dfT.index.dayofyear.values)
    # - GEOS5 downloads
    collection = 'inst3_3d_aer_Np'
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                     collection=collection)
    failed_doys = check4GEOS5_failed_downloads(dt=dt, folder=folder)
    n_failed_doys = len(failed_doys)
    sum_n_failed_doys += n_failed_doys
    pstr = 'GEOS5 - there were {} failed dowloads on {}'
    print(pstr.format(n_failed_doys, dt_fmt))
    if (n_failed_doys > 0):
        print('The failed days were: ', failed_doys)
    if (n_failed_doys > 0) and retry:
        print('')
        get_latest_GEOS5_fcast_data(dt=dt, doys2use=failed_doys,
                                    rm_existing_file=True)

    # Also download the GEOS5 Nv collection
    if inc_GEOS5_Nv_collection:
        collection = 'inst3_3d_aer_Nv'
        folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                         collection=collection)
        vars2use = ['du{:0>3}'.format(i) for i in range(1,6)]
        for var2use in vars2use:
            failed_doys = check4GEOS5_failed_downloads(dt=dt, folder=folder,
                                                       var2use=var2use)
            n_failed_doys = len(failed_doys)
    #        sum_n_failed_doys += n_failed_doys
            pstr = 'GEOS5 - there were {} failed dowloads for {} on {}'
            print(pstr.format(n_failed_doys, var2use, dt_fmt))
            if (n_failed_doys > 0):
                print('The failed days were: ', failed_doys)
            if (n_failed_doys > 0) and retry:
                print('')
                get_latest_GEOS5_fcast_data(dt=dt, doys2use=failed_doys,
                                            collection=collection,
                                            vars2use=[var2use],
                                            rm_existing_file=True)

    # - GEOSCF downloads
    collection = 'chm_inst_1hr_g1440x721_p23'
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                     collection=collection)
    vars2use = list(convert_GEOSCF_var2GEOSChem_name(rtn_dict=True).keys())
    # First check for lon_slice
    limit_lats = False
    limit_lons = True
    limit_lvls = False
    # Loop by variable
    for var2use in vars2use:
        failed_doys = check4GEOSCF_failed_downloads(dt=dt, folder=folder,
                                                    var2use=var2use,
                                                    limit_lons=limit_lons,
                                                    limit_lats=limit_lats,
                                                    limit_lvls=limit_lvls,
                                                    debug=debug)
        n_failed_doys = len(failed_doys)
        sum_n_failed_doys += n_failed_doys
        pstr = "GEOSCF (lon slice): {} failed downloads for '{}' on {}"
        print(pstr.format(n_failed_doys, var2use, dt_fmt))
        if (n_failed_doys > 0):
            print('The failed days were: ', failed_doys)
        if (n_failed_doys > 0) and retry:
            print('')
            get_latest_GEOSCF_fcast_data(dt=dt,
                                         vars2use=[var2use],
                                         doys2use=failed_doys,
                                         limit_lvls=limit_lvls,
                                         limit_lats=limit_lats,
                                         limit_lons=limit_lons,
                                         rm_existing_file=True,
                                         debug=debug)

    # Then lat_slice
    limit_lats = True
    limit_lons = False
    limit_lvls = False
    # Loop by variable
    for var2use in vars2use:
        failed_doys = check4GEOSCF_failed_downloads(dt=dt, folder=folder,
                                                    var2use=var2use,
                                                    limit_lons=limit_lons,
                                                    limit_lats=limit_lats,
                                                    limit_lvls=limit_lvls,
                                                    debug=debug)
        n_failed_doys = len(failed_doys)
        sum_n_failed_doys += n_failed_doys
        pstr = "GEOSCF (lat slice) - {} failed downloads for '{}' on {}"
        print(pstr.format(n_failed_doys, var2use, dt_fmt))
        if (n_failed_doys > 0):
            print('The failed days were: ', failed_doys)
        # If the doys do not line up... delete the files and start again

        # Download missing files
        if (n_failed_doys > 0) and retry:
            print('')
            get_latest_GEOSCF_fcast_data(dt=dt,
                                         vars2use=[var2use],
                                         doys2use=failed_doys,
                                         limit_lvls=limit_lvls,
                                         limit_lats=limit_lats,
                                         limit_lons=limit_lons,
                                         rm_existing_file=True,
                                         debug=debug)

    # Then alt_slice
    limit_lats = False
    limit_lons = False
    limit_lvls = True
    for var2use in vars2use:
        failed_doys = check4GEOSCF_failed_downloads(dt=dt, folder=folder,
                                                    var2use=var2use,
                                                    limit_lons=limit_lons,
                                                    limit_lats=limit_lats,
                                                    limit_lvls=limit_lvls,
                                                    debug=debug)
        n_failed_doys = len(failed_doys)
        sum_n_failed_doys += n_failed_doys
        pstr = "GEOSCF (alt slice) - {} failed dowloads for '{}' on {}"
        print(pstr.format(n_failed_doys, var2use, dt_fmt))
        if (n_failed_doys > 0):
            print('The failed days were: ', failed_doys)
        if (n_failed_doys > 0) and retry:
            print('')
            get_latest_GEOSCF_fcast_data(dt=dt,
                                         vars2use=[var2use],
                                         doys2use=failed_doys,
                                         limit_lvls=limit_lvls,
                                         limit_lats=limit_lats,
                                         limit_lons=limit_lons,
                                         rm_existing_file=True,
                                         debug=debug)
    # Return the total number of failed download doys
    return sum_n_failed_doys


def is_dt_latest_GEOSCF(dt=None, mode='fcast',
                        collection='chm_inst_1hr_g1440x721_p23'):
    """
    Check if the latest GEOS-CF run is available on OPeNDAP
    """
    # Just use noon yesterday if no date tis provided
    if isinstance(dt, type(None)):
        TNow = AC.time2datetime( [gmtime()] )[0]
        dt = AC.add_days(TNow, -1)
        dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
    # Check  containers
#    ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode, date=dt)
    # Temporarily do not allow passing of date to OPeNDAP fetcher
    ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode,
                                         date=None)
    last_start_date = AC.dt64_2_dt( [ds.time.values[0]] )[0]
    # Return check of files
    if last_start_date == dt:
        return True
    else:
        return False


def is_dt_latest_GEOS5(dt=None, mode='fcast',
                       collection='chm_inst_1hr_g1440x721_p23'):
    """
    Check if the latest GEOS-5 run is available on OPeNDAP
    """
    # Just use noon yesterday if no date tis provided
    if isinstance(dt, type(None)):
        TNow = AC.time2datetime( [gmtime()] )[0]
        dt = AC.add_days(TNow, -1)
        dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
    # Check  containers
    # Temporarily do not allow passing of date to OPeNDAP fetcher
    ds = AC.get_GEOS5_as_ds_via_OPeNDAP(collection=collection, mode=mode,
                                        dt=dt)
    last_start_date = AC.dt64_2_dt( [ds.time.values[0]] )[0]
    # Return check of files
    if last_start_date == dt:
        return True
    else:
        return False


def download_GEOSCF_fcast_data4date(dt=None, ds=None, mode='fcast',
                                    collection='chm_inst_1hr_g1440x721_p23',
                                    vars2use=None, doys2use=None,
                                    rm_existing_file=False,
                                    limit_lvls=True, limit_lons=False,
                                    limit_lats=False):
    """
    Retrieve the latest GEOS-CF data
    """
    # Get a pointer at the latest forecast data
    if isinstance(ds, type(None)):
#        ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode, date=dt)
        ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode,
                                             date=None)
    # use the latest forecast data in the dataset if one not provided
    if isinstance(dt, type(None)):
        dt = AC.dt64_2_dt( [ds.time.values[0]] )[0]
    pstr = 'Attempting to download GEOS-CF data ({}) starting on {}'
    print( pstr.format( collection, dt.strftime('%Y/%m/%d %H:%M')) )
    # Where to save the data
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                     collection=collection)
    # Check if the folder is present, if not create it
    if os.path.isdir(folder):
        print('WARNING: folder already exists ({})'.format(folder))
        print('(IF DATA IS PRESENT IT WILL BE OVERWRITEN - CHECK THIS!)')
        pstr = 'SETTINGS: limit_lvls={}; limit_lons={}; limit_lats={}'
        print(pstr.format(limit_lvls, limit_lons, limit_lats))
    else:
        os.makedirs(folder)
        print('Created folder for data ({})'.format(folder))
    # - Get the 3D data
    # Now retrieve the data -
    get_GEOSCF_data_cubes4collection(ds=ds, mode=mode, collection=collection,
                                     folder=folder, vars2use=vars2use,
                                     doys2use=doys2use, dt=dt,
                                     limit_lvls=limit_lvls,
                                     limit_lons=limit_lons,
                                     rm_existing_file=rm_existing_file,
                                     limit_lats=limit_lats )

    # - Get the 2D data
#    collection = 'chm_inst_1hr_g1440x721_p23'
#    ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode)
#    get_GEOSCF_data_cubes4collection(ds=ds, mode=mode, collection=collection)


def when_should_GEOSCF_have_last_run_from():
    """
    GEOS-CF is run every midnight
    """
    # Time now (GMT)
    TNow = AC.time2datetime( [gmtime()] )[0]
    # Therefore most recent forecast started when?
    dt = datetime.datetime(TNow.year, TNow.month, TNow.day, 12, )
    # Return as datetime
    return dt


def what_is_latest_data_locally(only_check_last5dates=True, debug=True):
    """
    Check what the latest data is locally
    """
    # Where to look for data
    folder = get_local_folder('NASA_data')
    folder += '/GEOS_CF/ARNA/fcast/'
    # What is the latest data that has been?
    subfolders = list(sorted(glob.glob(folder+'/*')))
    # Only look at the data in the last 5 folders for expendiency
    if only_check_last5dates:
        subfolders = subfolders[-5:]
    # Convert the subfolder names into datetimes
    dates4folders = [i.split(folder)[-1].split('_') for i in subfolders]
    print(dates4folders)
#    dates4folders = [i for i in dates4folders if ('z' not in i) ]
    dates4folders = [[int(ii) for ii in i[:3]] for i in dates4folders ]
    dates4folders = [datetime.datetime(*i) for i in dates4folders]
    # Does the folder actually contain the data?
    CF_data_present_list = []
    for n_folder, folders in enumerate( subfolders ):
        date = dates4folders[n_folder]
        if debug:
            pstr = "Checking data in folder for date '{}': {}"
            print(pstr.format(date, folder))
        CF_data_present_list += [is_GEOSCF_data_in_folder(folder)]
    #
    set_of_bools = list(set(CF_data_present_list))
    if len(set_of_bools) ==1:
        if set_of_bools[0] == True:
            if debug:
                print('Data present locally for all dates checked!')
            return dates4folders[-1]
        else:
            pstr = 'WARNING: No local data present 4 the last 5 dates checked!'
            print(pstr)
            return None


def is_GEOSCF_data_in_folder(folder):
    """
    check if GEOS-CF data is in the folder
    """
    # TODO - write a checker for GEOS-CF
    return False


def set_limits4ar_plotted_range(var, verbose=False):
    """
    Retrieve plotting limits for a given species
    """
    # Default is that there are no plotting limits
    vmin = None
    vmax = None
    # Dust (GEOS-5 product)
    if 'Dust' in var:
        vmin = 5
        vmax = 135
    # All particulate matter
    elif 'PM' in var:
        vmin = 5
        vmax = 50
    # All particulate matter
    elif 'CO' in var:
        vmin = 80
        vmax = 150
    # All particulate matter
    elif 'O3' in var:
        vmin = 0
        vmax = 50
    # NO2
    elif 'NO2' in var:
        vmin = 0
        vmax = 30/1E3
    # NOy
    elif 'NOy' in var:
        vmin = 0
        vmax = 1
    else:
        if verbose:
            print('WARNING: No range set for {}'.format(var))
    return vmin, vmax


def select_smaller_area_around_Cape_Verde(ds):
    """
    Sub-select a smaller region just around Cape Verde
    """
    # What is the smaller area around Cape Verde defined as?
#    lat_min = 15 # used for initial work
    lat_min = 14
    lat_max = 18
    lon_max = -21
    lon_min = -26
    # Create booleans to sub select
    bool1 = ((ds.lon >= lon_min) & (ds.lon <= lon_max)).values
    bool2 = ((ds.lat >= lat_min) & (ds.lat <= lat_max)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    return ds


def plt_smaller_area_around_CVAO(ds, LonVar='lon',
                                 LatVar='lat',
                                 use_local_CVAO_area=False):
    """
    Plot up the area around Cape Verde
    """
    # Make a new variable using an old one as a template
    var2copy = [i for i in ds.data_vars][0]
    NewVar = 'Cape_Verde_LOCAL'
    dsCV = ds.copy()
    dsCV[NewVar] = ds[var2copy].copy()
    # Set all values to zero
    arr = dsCV[NewVar].values
    arr[:] = 1
    dsCV[NewVar].values = arr
    # Set extents - from main dataset
    x0 = float(ds[LonVar].min())
    x1 = float(ds[LonVar].max())
    y0 = float(ds[LatVar].min())
    y1 = float(ds[LatVar].max())
    extents = (x0, x1, y0, y1)
    # Sub-select data
    dsCV = select_smaller_area_around_Cape_Verde(dsCV)
    # Select surface and average over time
    dsCV = dsCV.sel(lev=dsCV.lev[0])
    dsCV = dsCV.mean(dim='time')
    # Now plot
    quick_map_plt_CV_1layer(dsCV, var2plot=NewVar, extents=extents,
                                     use_local_CVAO_area=use_local_CVAO_area,
                                     extra_str='ARNA', save_plot=True)


def plt_timeseries4ds(ds, region='Cape_Verde', extr_str='',
                      vars2use=None, year=2018, verbose=False,
                      show_plot=False, dpi=320):
    """
    Plot timeseries of data at different heights
    """
    # Which variables to include in analysis?
    if not isinstance(vars2use, list):
        vars2use = [i for i in ds.data_vars]
    # - Now plot up species as PDf based on level
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context("talk", font_scale=0.75)
    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_timeseries_{}_{}'.format(region, year)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    color = 'red'
    alpha = 0.3
    xtickrotation = 45
    # Now loop and plot
    for lev in ds.lev.values:
        # Get vars for levels
        #        suffix = '-{:.0f}hPa'.format(lev)
        #        vars2plot = [i for i in df.columns if i.endswith(suffix) ]
        for var in vars2use:
            #            var = df_var.split(suffix)[0]
            print(lev, var)
            # Get the LateX for of the species name
            try:
                LaTeX_spec = AC.latex_spec_name(var)
            except KeyError:
                print('WARNING: not converted {} to LaTeX form'.format(var))
                LaTeX_spec = var
            # Retrive units
            units = ds[var].units
            # Set title
            title_str = 'Timeseries of [{}] @ {:.0f}hPa in Feb {}'
            title = title_str.format(LaTeX_spec, lev, year)
            # Setup the plot
            fig = plt.figure(figsize=(10, 6))
            ax = fig.add_subplot(111)
            # Setup a dictionary of colours
#            color_dict = dict(zip(lvls, AC.get_CB_color_cycle() ) )
            # Select data for species and level
            da = ds.sel(lev=lev)[var]
            # For sorting by time
            da = da.sortby('time')
            # Get mean, min and max
            mean_vals = da.mean(dim=['lat', 'lon'], skipna=True).values
            min_vals = da.min(dim=['lat', 'lon'], skipna=True).values
            max_vals = da.max(dim=['lat', 'lon'], skipna=True).values
            time = AC.dt64_2_dt(da.time.values)
            # Now plot up as a time series
            if mean_vals[np.isfinite(mean_vals)].shape[0] == 0:
                prtstr = 'WARNING: Not plotting {} @ {:.0f}hPa because all NaNs!'
                print(prtstr.format(var, lev))
                pass
            else:
                ax.plot(time, mean_vals, color=color, label=var,
                        #                linewidth=lw, ls=ls
                        )
                ax.fill_between(time, mean_vals, min_vals,
                                alpha=alpha, facecolor=color)
                ax.fill_between(time, mean_vals, max_vals,
                                alpha=alpha, facecolor=color)
                # Add a title
                plt.title(title)
                # Add a legend
                plt.legend()
                # Rotate dates by 45 degrees
#                xticks = ax.get_xticks()
#                labels = ax.get_xticklabels()
#                ax.set_xticks(xticks)
#                ax.set_xticklabels(labels, rotation=xtickrotation)
#                ax.set_xticklabels(labels)
#                print(labels)
#                ax.set_xticklabels(xticks, rotation=xtickrotation)
                # Make sure that all titles and labels are readable
                try:
                    plt.tight_layout()
                except:
                    pass
                # Save to PDF
                AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
                if show_plot:
                    plt.show()
            plt.close()
            del da
    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


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


def PDF_on_species_in_ds4lvls(ds, region='Cape_Verde', extr_str='',
                              vars2use=None, year=2018, verbose=False,
                              show_plot=False, dpi=320):
    """
    plot the stats on species at different heights as a PDF
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
            S = pd.Series(vals)
            # Save to the main DataFrame
            df[varname] = S
            del vals, S

    # - Now plot up species as PDf based on level
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context("talk")
    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_PDF_of_concs_by_level_{}_{}'.format(region, year)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    # Now loop and plot
    for var in vars2use:
        #    for var in vars2use[:2]:
        # Get vars for levels
        vars2plot = [i for i in df.columns if i.startswith(var+'-')]
        # Get the LateX for of the species name
        try:
            LaTeX_spec = AC.latex_spec_name(var)
        except KeyError:
            print('WARNING: not converted {} to LaTeX form'.format(var))
            LaTeX_spec = var
        # Retrive units
        units = ds[var].units
        # Set title
        title = 'PDF of [{}] for hPa levels in Feb {}'.format(LaTeX_spec, year)
        # Setup the plot
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        # Setup a dictionary of colours
        lvls = np.array(ds.lev.values)
        color_dict = dict(zip(lvls, AC.get_CB_color_cycle()))
        # Plot up by level
        for lev2use in lvls:
            levstr = '{:.0f}hPa'.format(lev2use)
            var2plot = [i for i in vars2plot if levstr in i]
            assert len(
                var2plot) == 1, 'ERROR: There must only be one variable per level!'
            if verbose:
                print(var, lev2use, var2plot)
            # Retrieve the data
            arr = df[var2plot].values
            # only consider values that are not NaN
            arr = arr[~np.isnan(arr)]
            # Plot up and add title
            axlabel = '{} ({})'.format(LaTeX_spec, units)
            ax = sns.distplot(arr, axlabel=axlabel, label=levstr,
                              color=color_dict[lev2use], ax=ax)
            # Make sure the scaling is correct
            ax.autoscale()
        # Add a title
        plt.title(title)
        # Add a legend
        plt.legend()
        # Make sure that all titles and labels are readable
        plt.tight_layout()
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()
    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def add_extra_derived_species(ds, debug=False):
    """
    Add extra derived variables of interest
    """
    # Add NOy divided by PM2.5 concentration
#     NewVar = 'NOy/PM2.5(dust)'
#     ds[NewVar] = ds['NOy'] / ds[ 'PM2.5(dust)']
#     attrs = ds[NewVar].attrs
#     attrs['units'] = 'unitless'
#     ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration
    NewVar = 'NOy*PM2.5(dust)'
    ds[NewVar] = ds['NOy'] * ds['PM2.5(dust)']
    attrs = ds[NewVar].attrs
    attrs['units'] = 'unitless'
    ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration (only for where dust >50 ug/m3)
#     CopyVar =  'NOy'
#     NewVar = 'NOy/PM2.5(dust>50ug/m3)'
#     ds[NewVar] = ds[CopyVar].copy()
#     vals2use = ds[ 'PM2.5(dust)'].copy().values
#     vals2use[vals2use<50] = np.NaN
#     ds[NewVar].values = ds['NOy'].values / vals2use
#     attrs = ds[NewVar].attrs
#     attrs['units'] = 'unitless'
#     ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration (only for where dust >25 ug/m3)
    CopyVar = 'NOy'
    NewVar = 'NOy/PM2.5(dust>10ug/m3)'
    ds[NewVar] = ds[CopyVar].copy()
    vals2use = ds['PM2.5(dust)'].copy().values
    vals2use[vals2use < 10] = np.NaN
    ds[NewVar].values = ds['NOy'].values / vals2use
    attrs = ds[NewVar].attrs
    attrs['units'] = 'unitless'
    ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration (only for where dust >25 ug/m3)
    CopyVar = 'NOy'
    NewVar = 'NOy/PM2.5(dust>25ug/m3)'
    ds[NewVar] = ds[CopyVar].copy()
    vals2use = ds['PM2.5(dust)'].copy().values
    vals2use[vals2use < 25] = np.NaN
    ds[NewVar].values = ds['NOy'].values / vals2use
    attrs = ds[NewVar].attrs
    attrs['units'] = 'unitless'
    ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration (only for where dust >25 ug/m3)
    CopyVar = 'NOy'
    NewVar = 'NOy(where>10ug/m3 dust)'
    ds[NewVar] = ds[CopyVar].copy()
    NOy = ds[NewVar].values
    dust = ds['PM2.5(dust)'].values
    NOy[dust < 10] = np.NaN
    ds[NewVar].values = NOy
    attrs = ds[NewVar].attrs
    attrs['units'] = 'ppbv'
    ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration (only for where dust >0.5 ug/m3)
    CopyVar = 'NOy'
    NewVar = 'NOy/PM2.5(dust>0.5ug/m3)'
    ds[NewVar] = ds[CopyVar].copy()
    vals2use = ds['PM2.5(dust)'].copy().values
    vals2use[vals2use < 0.5] = np.NaN
    ds[NewVar].values = ds['NOy'].values / vals2use
    attrs = ds[NewVar].attrs
    attrs['units'] = 'unitless'
    ds[NewVar].attrs = attrs
    # print data variables if debugging...
    if debug:
        print('-'*10, [i for i in ds.data_vars], '-'*10)

    return ds


def get_species_units(input):
    """
    Retrieve units for GEOS-CF species and derived values
    """
    if 'NOy/' in input:
        units = 'unitless'
    elif 'pm25' in input:
        units = '$\mu$g m$^{-3}$'
    elif 'Dust' in input:
        units = '$\mu$g m$^{-3}$'
    elif 'NOy' in input:
        units = 'ppbv'
    else:
        units = 'v/v'
    return units


def get_vmin_value4var(input):
    """
    Set a vmin value for a variable (below which values are set to NaN)
    """
    #
    if 'Dust' in input:
        vmin = 15
    elif 'NOy' in input:
        vmin = 0.5
    else:
        print('WARNING: vmin case not set for variable! ({})'.format(input))
    return vmin


def get_cmap4var(input):
    """
    Get a colour-map for a specific variable.
    """
    #
    if 'Dust' in input:
        if input == 'Dust':
            nticks=9-1
#            nticks=10-1
#             nticks=11-1
#             nticks=12-1
#             ticks = np.linspace(vmin2, vmax2, nticks+1)
#             ticks = [15*(i+1) for i in range(9)]
            ticks = [15*(i) for i in range(11)]
        # for PM2.5 dust
        else:
            # Get min and max values to plot
            vmin2, vmax2 = set_limits4ar_plotted_range(input)
            vmin2 = get_vmin_value4var(input)
            # Now set the number of ticks and range
            nticks=8-1
            ticks = np.linspace(vmin2, vmax2, nticks+1)
        cmap = get_reduced_cmap(cmap='Reds', maxval=0.75)
#        cmap = get_reduced_cmap(cmap='Greens', maxval=0.75)
#        cmap = get_reduced_cmap(cmap='Greys', maxval=0.75)
        cmap = AC.mk_discrete_cmap( nticks=nticks, cmap=cmap)
    elif input == 'NOy':
        ticks = None
        nticks = 6-1
        cmap = AC.mk_discrete_cmap( nticks=nticks, cmap='Blues')
    else:
        print('WARNING: vmin case not set for variable! ({})'.format(input))
    return cmap, ticks, nticks


def convert_GEOSCF_units(ds, debug=False):
    """
    Convert from GEOS-CF units to GEOS-Chem units
    """
    vars2use = [i for i in ds.data_vars]
    for var in vars2use:
        # Get standard conversion units
        try:
            units, scaleby = AC.tra_unit(var, scale=True)
            if debug:
                print(var, units, scaleby)
        except KeyError:
            print('WARNING: setting {} to unitless!'.format(var))
            scaleby = 1
            units = 'unitless'
        # update the data by the scaling
        ds[var].values = ds[var].values.copy() * scaleby
        # updates units
        attrs = ds[var].attrs.copy()
        attrs['units'] = units
        ds[var].attrs = attrs
        # Do some garbage cleaning
        gc.collect()
    return ds


def update_GEOSCF_names2GEOSChem_names(ds):
    """
    Update the GEOS-CF names to be in GEOS-Chem form
    """
    # Get a dictionary of the GEOS-CF and GEOS-Chem names
    d = convert_GEOSCF_var2GEOSChem_name(rtn_dict=True)
    return ds.rename(name_dict=d)


def convert_GEOSCF_var2GEOSChem_name(input=None, rtn_dict=False):
    """
    Convert the GEOS-CG species name to GEOS-Chem name
    """
    d = {
        'noy': 'NOy', 'co': 'CO', 'no2': 'NO2', 'o3': 'O3',
        'pm25_rh35_gcc': 'PM2.5',
        'pm25du_rh35_gcc': 'PM2.5(dust)',
        #    'pm25su_rh35_gcc' : 'PM2.5(SO4)',
        #    'pm25ni_rh35_gcc' : 'PM2.5(NIT)',
        #    'pm25soa_rh35_gc' : 'PM2.5(SOA)',
        #    'pm25ss_rh35_gcc' : 'PM2.5(SSA)',
        #    'pm25bc_rh35_gcc' : 'PM2.5(BC)',
        #    'pm25oc_rh35_gcc' : 'PM2.5(OC)',
        'so2': 'SO2',
    }
    if rtn_dict:
        return d
    else:
        return d[input]


def plt_avg_spatial_by_lvl(ds, year=2018, vars2use=None,
                           verbose=False,
                           use_local_CVAO_area=False,
                           show_plot=False, dpi=320):
    """
    Make a PDF of average spatial concentrations by level
    """
    # Setup PDF to save plots to.
    savetitle = 'ARNA_avg_spatial_concs_{}'.format(year)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    # Which variables to plot
    if not isinstance(vars2use, list):
        vars2use = [i for i in ds.data_vars]
    vars2use = ['NOy', 'PM2.5(dust)', ]
    # Now loop and plot
    for var in vars2use:
        # Plot up by level
        for lev2use in ds.lev.values:
            if verbose:
                print(var, lev2use)
            # Get units for species
            units = ds[var].units
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lev=lev2use).mean(dim='time')
            # Get the LateX for of the species name
            try:
                LaTeX_spec = AC.latex_spec_name(var)
            except KeyError:
                print('WARNING: not converted {} to LaTeX form'.format(var))
                LaTeX_spec = var
            # Set title
            title = 'Average [{}] @ {:.0f}hPa in Feb {}'.format(
                LaTeX_spec, lev2use, year)
            # Plot up and add title
            quick_map_plt_CV_1layer(ds_tmp, var2plot=var, title=title,
                                    use_local_CVAO_area=use_local_CVAO_area,
                                    save_plot=False, units=units)
            del ds_tmp
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()
    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plot_average_spatial_concs4lon(ds, year=2018, vars2use=None,
                                   use_local_CVAO_area=False,
                                   show_plot=False, dpi=320,
                                   verbose=False,):
    """
    Make a PDF of average spatial concentrations by level
    """
    # Setup PDF to save plots to.
    savetitle = 'ARNA_avg_spatial_concs_{}'.format(year)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    # Which variables to plot
    if not isinstance(vars2use, list):
        vars2use = [i for i in ds.data_vars]
    vars2use = ['NOy', 'PM2.5(dust)', ]
    # Now loop and plot
    for var in vars2use:
#    for var in vars2use[:2]:  # For testing
        # Plot up by level
        #        for lev2use in ds.lev.values:
        for lon2use in list(ds.lon.values)  :
            if verbose:
                print(var, lev2use)
            # Get units for species
            units = ds[var].units
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lon=lon2use).mean(dim='time')
            # Get the LateX for of the species name
            try:
                LaTeX_spec = AC.latex_spec_name(var)
            except KeyError:
                print('WARNING: not converted {} to LaTeX form'.format(var))
                LaTeX_spec = var
            # Set title
            title = 'Average [{}] @ {:.0f}hPa in Feb {}'.format(
                LaTeX_spec, lev2use, year)
            # Plot up and add title
#            quick_map_plt_CV_1layer(ds_tmp, var2plot=var, title=title,
#                                             use_local_CVAO_area=use_local_CVAO_area,
#                                             save_plot=False, units=units)
            # vertical plot
            del ds_tmp
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()
    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def quick_map_plt_CV_1layer(ds, var2plot=None, extra_str='',
                            projection=ccrs.PlateCarree(),
                            save_plot=True, show_plot=False,
                            savename=None,
                            units=None, title=None,
                            LatVar='lat', LonVar='lon', fig=None,
                            ax=None,
                            extents=None,
                            region='Cape_Verde',
                            use_local_CVAO_area=True,
                            add_flyable_range_as_circle=True,
                            add_flyable_range=False,
                            add_detailed_map=True,
                            add_ARNA_locs=True,
                            extend='neither', folder='./',
                            dpi=320):
    """
    Plot up a quick spatial plot of data using cartopy

    Parameters
    -------
    ds (xr.Dataset): dataset object holding data to plot
    var2plot (str): variable to plot within the dataset
    LatVar, LonVar (str): variables to use for latitude and longitude
    save_plot (bool): save the plot as a .png ?
    show_plot (bool): show the plot on screen
    dpi (int): resolution to use for saved image (dots per square inch)
    savename (str): name to use for png of saved .png
    extra_str (str): extra string to append to save .png
    projection (cartopy.crs obj.):  projection to use
    fig (figure instance): matplotlib figure instance
    ax (axis instance): axis object to use

    Returns
    -------
    (None)
    """
    # Use the 1st data variable if not variable given
    if isinstance(var2plot, type(None)):
        pstr = 'WARNING: No variable to plot(var2plot), trying 1st data_var'
        print(pstr)
        var2plot = list(ds.data_vars)[0]

    # Setup figure and axis and plot
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=(10, 6))
#        fig = plt.figure()
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111, projection=projection, aspect='auto')
    # Setup plotted range
    vmin, vmax = set_limits4ar_plotted_range(var2plot)
    # Now plot up
    ds[var2plot].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                             transform=ccrs.PlateCarree(),
                             vmin=vmin, vmax=vmax, extend=extend)
    # Add some grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=.5, color='gray', alpha=0.25, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    # Limit plot to Cape Verde region
    if use_local_CVAO_area:
        x0 = -30
        x1 =-10
        y0 = 0
        y1 = 25
        extents = (x0, x1, y0, y1)
    # Mark known places to help geo-locate viewers
    if add_ARNA_locs:
#        colours = AC.get_CB_color_cycle()
        locs2plot  = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
        for loc2plot in locs2plot:
            lon, lat, alt = AC.get_loc(loc2plot)
            # Now plot up locations
            ax.plot(lon, lat, 'bo', markersize=5, markerfacecolor='none',
                    markeredgewidth=2,
                    zorder=10,
                    markeredgecolor='black',
                    transform=ccrs.PlateCarree())
            # Add a label for the location?
#            ax.text(lon, lat+0.25, loc2plot, transform=ccrs.PlateCarree())
    # Add a box to show the flyable range
    if add_flyable_range:
        # Get the minimum
        d = get_max_flying_range4BAE146()
        min_lon = d['min_lon']
        max_lon = d['max_lon']
        min_lat = d['min_lat']
        max_lat = d['max_lat']
        # Create new lists
        lons = [min_lon, min_lon, max_lon, max_lon]
        lats = [min_lat, max_lat, max_lat, min_lat]
        # Now plot as a linear ring
        ring = LinearRing(list(zip(lons, lats)))
        ax.add_geometries([ring], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='grey',
                          zorder=10, linestyle=':',
                          )
    if add_flyable_range_as_circle:
#        n_points = 1000
        # Approximate from James' max distance
        # ( 16.8331-13 ) *110667.45
        locs4circles = 'Dakar', 'Sao Vicente Airport',
        for loc in locs4circles:
            # Get locations to centre circle on
            lon, lat, alt = AC.get_loc(loc)
            # Radius in degrees
#            radius = 16.8331-13
            radius = 21 - 16.8331
            # Plot up circle
            ax.add_patch(mpatches.Circle(xy=[lon, lat],
                                         radius=radius,
                                         transform=projection,
                                         facecolor='none',
                                         edgecolor='grey',
                                         linestyle=':',
                                         zorder=10
                                         ))
    # Get limits of plotting data
    if isinstance(extents, type(None)):
        x0 = float(ds[LonVar].min())
        x1 = float(ds[LonVar].max())
        y0 = float(ds[LatVar].min())
        y1 = float(ds[LatVar].max())
        extents = (x0, x1, y0, y1)
    ax.set_extent(extents, crs=ccrs.PlateCarree())
    # Beautify the figure/plot
    if add_detailed_map:
        # Add borders for countries
        ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                       facecolor='none', zorder=50)
        # Also add minor islands (inc. Cape Verde)
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor=None,
                                                facecolor='none')
        ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)
    # Update the colour bar lanel
    if not isinstance(units, type(None)):
        im = ax.images
        cb = im[-1].colorbar
        cb.set_label('{}'.format(units))
    # Add a generic title if one is not provided
    if not isinstance(title, type(None)):
        plt.title(title)
    # Save the plot?
    if save_plot:
        if isinstance(savename, type(None)):
            savename = 'spatial_plot_{}_{}'.format(var2plot, extra_str)
        savename = AC.rm_spaces_and_chars_from_str(savename)
        plt.savefig(folder+savename+'.png', dpi=dpi)
    if show_plot:
        plt.show()


def get_max_flying_range4BAE146(calculate_online=False):
    """
    Get the maximum flying extent of the BAE146 aircraft
    """
    if calculate_online:
        min_lat, max_lat, min_lon, max_lon = 90, -90, 180, -180
        locs2plot  = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
        for loc2plot in locs2plot:
            lon, lat, alt = AC.get_loc(loc2plot)
            min_lat = min(min_lat, lat)
            max_lat = max(max_lat, lat)
            min_lon = min(min_lon, lon)
            max_lon = max(max_lon, lon)
        # What buffers to use for lat and lon
        # ~69 miles for latitude - so add buffer of 500/69 (7 degrees lat)
        lat_buffer = 7.0
        # ~55 miles for longitude - so add buffer of 9 degrees (500/55 = 9 degrees Lon)
        lon_buffer = 9.0
        # Update min and max lons
        min_lat -= lat_buffer
        max_lat += lat_buffer
        min_lon -= lon_buffer
        max_lon += lon_buffer
    else:
        # Use lats and lons.
        min_lat = 13
        max_lat = 21
        min_lon = -29
        max_lon = -21
    # Return data as a dictionary
    d = {
    'min_lat' : min_lat,
    'max_lat' : max_lat,
    'min_lon' : min_lon,
    'max_lon' : max_lon,
    'max_alt' : 8000
    }
    return d


def plot_spatial_concs_2layer(ds, show_plot=False, folder=None,
                              var2plot1='NOy', var2plot2='PM2.5(dust)',
                              extr_title_str='', region='Cape_Verde',
                              add_max_vals_as_txt=False,
                              testing_mode=False,
                              verbose=False, testing=False ):
    """
    Plot up a two layer plot on a single map for given levels
    """
    # Local variables
    try:
        LaTeX_spec1 = AC.latex_spec_name(var2plot1)
    except KeyError:
        LaTeX_spec1 = var2plot1
    try:
        LaTeX_spec2 = AC.latex_spec_name(var2plot2)
    except KeyError:
        LaTeX_spec2 = var2plot2
    # Set data below threshold to zero based on variable name
    ds = set_values_below_range2NaNs4spec(var=var2plot1, ds=ds)
    ds = set_values_below_range2NaNs4spec(var=var2plot2, ds=ds)
    # Set lists of variables to loop and plot
    if testing_mode:
        levs2use = [700]
        times2use = ds.time[:4]
    else:
        levs2use = ds.lev.values
        times2use =  ds.time.values
    # Plot up by level and time
    for time2use in times2use:
        for lev2use in levs2use:
            # Get time as human readable string
            dstr = AC.dt64_2_dt([time2use])[0].strftime('%Y/%m/%d %H:%M')
            # Print out status
            pstr = "'plotting 2layer @ {:.0f}hPa on {}"
            print(pstr.format(lev2use, dstr) )
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lev=lev2use, time=time2use)
            # Set title
            title = '[{}] & [{}] @ {:.0f}hPa on {}'.format(
                LaTeX_spec1, LaTeX_spec2, lev2use, dstr)
            # Add extra string to existing title string
            title += '\n '+ extr_title_str
            # Save plots
            extra_str = 'lev_{}_dt_{}'.format( lev2use, dstr )
            quick_map_plt_2layer(ds_tmp, var2plot1=var2plot1,
                                  folder=folder, region=region,
                                  var2plot2=var2plot2, title=title,
                                  add_max_vals_as_txt=add_max_vals_as_txt,
                                  save_plot=True, extra_str=extra_str)
            # Tidy up...
            plt.close('all')


def set_values_below_range2NaNs4spec(var=None, ds=None):
    """
    To improve aesthetics of plots, values below a certain threshold are removed
    """
    # Limit plotted NOy values to those above 0.5 pptv
    if var == 'NOy':
        arr = ds[var].values
        arr[arr<0.5] = np.NaN
        ds[var].values = arr
    # Limit Dust values to those about
    elif 'Dust' in var:
        arr = ds[var].values
        arr[arr<15] = np.NaN
        ds[var].values = arr
    else:
        pstr = "WARNING: No case set for variable '{}', so not restricting array values"
        print(pstr.format(var))
    return ds


def plt_spatial_2layer_vertical_lon(ds, show_plot=False, folder=None,
                                    var2plot1='NOy',
                                    var2plot2='PM2.5(dust)',
                                    extr_title_str=None,
                                    testing_mode=False,
                                    ):
    """
    Plot up a two layer plot on a single map for given levels
    """
    # Local variables
    try:
        LaTeX_spec1 = AC.latex_spec_name(var2plot1)
    except KeyError:
        LaTeX_spec1 = var2plot1
    try:
        LaTeX_spec2 = AC.latex_spec_name(var2plot2)
    except KeyError:
        LaTeX_spec2 = var2plot2
    # Set data below threshold to zero based on variable name
    ds = set_values_below_range2NaNs4spec(var=var2plot1, ds=ds)
    ds = set_values_below_range2NaNs4spec(var=var2plot2, ds=ds)
    # Set lists of variables to loop and plot
    if testing_mode:
        lons2use = [-24.]
        times2use = ds.time[:4]
    else:
        lons2use = ds.lon.values
        times2use = ds.time.values
    # Plot by time and lon
    for time2use in times2use:
        for lon2use in lons2use:
            # Get time as human readable string
            dstr = AC.dt64_2_dt([time2use])[0].strftime('%Y/%m/%d %H:%M')
            # Print out status
            pstr = "'plotting 2layer @ {:.1f}$^{}$W on {}"
            print( pstr.format( lon2use*-1, '{\circ}', dstr) )
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lon=lon2use, time=time2use)
            # Set title
            title_str = '[{}] & [{}] @ {:.1f}$^{}$W on {}'
            title = title_str.format( LaTeX_spec1, LaTeX_spec2, lon2use*-1,
                                    '{\circ}', dstr)
            if not isinstance(extr_title_str, type(None)):
                title += extr_title_str
            # Save plots
            extra_str = 'lon_{}E_dt_{}'.format( lon2use, dstr )
            # Update the long_names - var1
            attrs = ds_tmp[var2plot1].attrs
            attrs['long_name'] = LaTeX_spec1
            ds_tmp[var2plot1].attrs = attrs
            # Update the long_names - var2
            attrs = ds_tmp[var2plot2].attrs
            attrs['long_name'] = LaTeX_spec2
            ds_tmp[var2plot2].attrs = attrs
            # Now call plotter
            quick_lon_plot_2layer(ds_tmp, var2plot1=var2plot1,
                                  folder=folder,
                                  var2plot2=var2plot2, title=title,
                                  save_plot=True, extra_str=extra_str
                                  )
            plt.close('all')


def plt_spatial_2layer_vertical_lat(ds, show_plot=False, folder=None,
                                    var2plot1='NOy',
                                    var2plot2='PM2.5(dust)',
                                    extr_title_str=None,
                                    testing_mode=False,
                                    ):
    """
    Plot up a two layer plot on a single map for given levels
    """
    # Local variables
    try:
        LaTeX_spec1 = AC.latex_spec_name(var2plot1)
    except KeyError:
        LaTeX_spec1 = var2plot1
    try:
        LaTeX_spec2 = AC.latex_spec_name(var2plot2)
    except KeyError:
        LaTeX_spec2 = var2plot2
    # Set data below threshold to zero based on variable name
    ds = set_values_below_range2NaNs4spec(var=var2plot1, ds=ds)
    ds = set_values_below_range2NaNs4spec(var=var2plot2, ds=ds)
    # Set lists of variables to loop and plot
    if testing_mode:
        lats2use = [16.]
        times2use = ds.time[:4]
    else:
        lats2use = ds.lat.values
        times2use = ds.time.values
    # Plot up by level and time
    for time2use in times2use:
        for lat2use in lats2use:
            # Get time as human readable string
            dstr = AC.dt64_2_dt([time2use])[0].strftime('%Y/%m/%d %H:%M')
            # Print out status
            pstr = "'plotting 2layer @ {:.1f}$^{}$N on {}"
            print( pstr.format( lat2use, '{\circ}', dstr) )
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lat=lat2use, time=time2use)
            # Set title
            title_str = '[{}] & [{}] @ {:.1f}$^{}$N on {}'
            title = title_str.format( LaTeX_spec1, LaTeX_spec2, lat2use,
                                    '{\circ}', dstr)
            if not isinstance(extr_title_str, type(None)):
                title += extr_title_str
            # Save plots
            extra_str = 'lat_{}N_dt_{}'.format( lat2use, dstr )
            # Update the long_names - var1
            attrs = ds_tmp[var2plot1].attrs
            attrs['long_name'] = LaTeX_spec1
            ds_tmp[var2plot1].attrs = attrs
            # Update the long_names - var2
            attrs = ds_tmp[var2plot2].attrs
            attrs['long_name'] = LaTeX_spec2
            ds_tmp[var2plot2].attrs = attrs
            # Now call plotter
            quick_lat_plt_2layer(ds_tmp, var2plot1=var2plot1,
                                  folder=folder,
                                  var2plot2=var2plot2, title=title,
                                  save_plot=True, extra_str=extra_str
                                  )
            plt.close('all')


def mk_video_from_plts(FileStr=None):
    """
    Convert images frames into videos using ffmpeg
    """
    FileStr = 'spatial_plot_NOy_PM2_5_dust__lev_1000_0_dt_*.png'
    #
    print('WARNING: ffmpeg calls here have been switched off')
#    ffmpegstr = "ffmpeg -framerate 3 -pattern_type glob -i '{FileStr}' -c:v libx264 -pix_fmt yuv420p out.mp4"

#ffmpeg -framerate 3 -pattern_type glob -i 'spatial_plot_NOy_PM2_5_dust__lev_500_0_dt_2019_11_*.png' -c:v libx264 -pix_fmt yuv420p out_2019_11_06_5day_fcast_500hPa.mp4
    os.system(" ")


def plot_CVAO_region_on_global_map(ds, var2use='NOy'):
    """
    Plot a global map to show ARNA campaign region
    """
    # - extents to use
    # extracted data from OPeNDAP
#     x0 = -30
#     x1 =-10
#     y0 = 0
#     y1 = 25
    # Local area analysed as Cape Verde
    x0 = -30
    x1 =-10
    y0 = 0
    y1 = 25
    # - Select the data
    # Just get an example dataset
    ds = ds[[var2use]]
    # Select a single level and time
    ds = ds.sel(time=ds.time[0])
    ds = ds.sel(lev=ds.lev[0])
    # Set values region
    bool1 = ((ds.lon >= x0) & (ds.lon <= x1)).values
    bool2 = ((ds.lat >= y0) & (ds.lat <= y1)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    # Set all values to 1
    arr = ds[var2use].values
    arr[:] = 1
    ds[var2use].values = arr
    # Plot the data
    projection = ccrs.Robinson()
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection=projection, aspect='auto', alpha=0.5)
    LatVar = 'lat'
    LonVar = 'lon'
    ds[var2use].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                             transform=ccrs.PlateCarree())
    # Beautify the figure/plot
    ax.coastlines()
    ax.set_global()
    # Force global perspective
    ax.set_global() # this will force a global perspective
    # Save
    savename = 'spatial_plot_Cape_Verde_flying_area'
    savename = AC.rm_spaces_and_chars_from_str(savename)
    plt.savefig(savename+'.png', dpi=dpi)


def get_reduced_cmap(cmap='Reds', minval=0.0, maxval=0.75, npoints=100):
    """
    Get a reduced colormap object (cmap)
    """
    import matplotlib
    # Load color map
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    # Make a truncated colour map
    Nstr = 'trunc({n},{a:.2f},{b:.2f})'
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                Nstr.format(n=cmap.name, a=minval, b=maxval),
                cmap(np.linspace(minval, maxval,  npoints))
                                                               )
    return cmap


def quick_map_plt_2layer(ds, var2plot1=None, var2plot2=None, extra_str='',
                         projection=ccrs.PlateCarree(), folder=None,
                         save_plot=True, show_plot=False, savename=None,
                         units=None, title=None,
                         LatVar='lat', LonVar='lon', fig=None, ax=None,
                         extents=None,
                         add_ARNA_locs=True,
                         use_local_CVAO_area=True, region='Cape_Verde',
                         extend='both',
                         add_flyable_range_as_box=False,
                         add_flyable_range_as_circle=True,
                         add_detailed_map=True,
                         add_max_vals_as_txt=False,
                         dpi=320):
    """
    Plot up a quick spatial plot of data using cartopy

    Parameters
    -------
    ds (xr.Dataset): dataset object holding data to plot
    var2plot (str): variable to plot within the dataset
    LatVar, LonVar (str): variables to use for latitude and longitude
    save_plot (bool): save the plot as a .png ?
    show_plot (bool): show the plot on screen
    dpi (int): resolution to use for saved image (dots per square inch)
    savename (str): name to use for png of saved .png
    extra_str (str): extra string to append to save .png
    projection (cartopy.crs obj.):  projection to use
    fig (figure instance): matplotlib figure instance
    ax (axis instance): axis object to use
    add_ARNA_locs (bool):
    use_local_CVAO_area (bool):

    Returns
    -------
    (None)
    """
    # Use the 1st data variable if not variable given
    if isinstance(var2plot1, type(None)):
        pstr = 'WARNING: No variable to plot (var2plot), trying 1st data_var'
        print(pstr)
        var2plot1 = list(ds.data_vars)[0]
#        var2plot1 = 'NOy'
    if isinstance(var2plot2, type(None)):
        pstr = 'WARNING: No variable to plot (var2plot), trying 1st data_var'
        print(pstr)
        var2plot2 = list(ds.data_vars)[0]
#        var2plot2 = 'PM2.5(dust)'
    # Setup figure and axis and plot
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=(10, 6))
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111, projection=projection, aspect='auto')

    # - Plot first var
    alpha = 0.5
    # Setup plotted range
    vmin1, vmax1 = set_limits4ar_plotted_range(var2plot1)
    vmin1 = get_vmin_value4var(var2plot1)
    units1 = get_species_units(var2plot1)
    cmap, ticks, nticks = get_cmap4var(var2plot1)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap':cmap, 'extend':extend, }
    else:
        cbar_kwargs  = {'ticks':ticks, 'cmap':cmap, 'extend':extend,  }
    # Now plot up var1
    ds[var2plot1].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                              transform=ccrs.PlateCarree(),
                              vmin=vmin1, vmax=vmax1,
                              zorder=1, alpha=alpha,
                              cmap=cmap,
                              cbar_kwargs=cbar_kwargs,
#                             extend=extend,
                              )
    # Update the units on the colour bar panel
    im = ax.images
    cb = im[-1].colorbar
    try:
        LaTeX_spec1 = AC.latex_spec_name(var2plot1)
    except KeyError:
        LaTeX_spec1 = var2plot1
    cb.set_label('{} ({})'.format(LaTeX_spec1, units1))
    # - Plot second var
    alpha = 0.4
    # Now plot up var 2
    vmin2, vmax2 = set_limits4ar_plotted_range(var2plot2)
    vmin2 = get_vmin_value4var(var2plot2)
    units2 = get_species_units(var2plot2)
    cmap, ticks, nticks = get_cmap4var(var2plot2)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap': cmap, 'extend' : extend, }
    else:
        cbar_kwargs  = {'ticks': ticks, 'cmap': cmap, 'extend' : extend, }
    # Now plot up var2
    ds[var2plot2].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                              transform=ccrs.PlateCarree(),
                              vmin=vmin2, vmax=vmax2,
                              zorder=1, alpha=alpha, cmap=cmap,
                              cbar_kwargs=cbar_kwargs,
#                             extend=extend,
                              )
    # Update the units on the colour bar panel
    im = ax.images
    cb = im[-1].colorbar
    try:
        LaTeX_spec2 = AC.latex_spec_name(var2plot2)
    except KeyError:
        LaTeX_spec2 = var2plot2
    cb.set_label('{} ({})'.format(LaTeX_spec2, units2))
    # Update the xaxis ticks (lon) to deg W
#     update_lon_units = True
#     if update_lon_units:
#         xticks = ax.get_xticks()
#         xticklabels = ax.get_xticklabels()
#         ax.set_xticklabels([str(i)*-1 for i in xticks])

    # - Update plot aesthetics
    # Add some grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=.5, color='gray', alpha=0.25, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    # Just plot over the CVAO region?
    if use_local_CVAO_area and (region != 'Cape_Verde_Flying'):
        x0 = -30
        x1 =-10
        y0 = 0
        y1 = 25
        extents = (x0, x1, y0, y1)
    elif (region == 'Cape_Verde_Flying'):
        x0 = -29.1
        x1 =-15.9
        y0 = 11.9
        y1 = 21.1
        extents = (x0, x1, y0, y1)
    # Add extra lat and lon grid libnes
    if (region == 'Cape_Verde_Flying'):
        # Which X tickst to use?
        # x axis
        xticks = np.arange(ds.lon.values.min(), ds.lon.values.max(), 0.2   )
        # y axis
        yticks = np.arange(ds.lat.values.min(), ds.lat.values.max(), 0.2   )
        # --- New approach
#        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#        LON_FORMATTER = LONGITUDE_FORMATTER(dms=True)#, auto_hide=True)
#        LON_FORMATTER.set_locs(xticks)
#        LAT_FORMATTER = LATITUDE_FORMATTER(dms=True)#, auto_hide=True)
#        LAT_FORMATTER.set_locs(yticks)
#        gl.xformatter = LON_FORMATTER
#        gl.yformatter = LAT_FORMATTER
        # --- Old approach
#        xmajor = np.arange(ds.lon.values.min(), ds.lon.values.max(), 1   )
#        xminor = [i for i in xticks if i not in xmajor]
#        xminor_locator = mticker.FixedLocator(xminor)
#        xmajor_locator = mticker.FixedLocator(xmajor)
#        ax.xaxis.set_major_locator(xmajor_locator)
#        ax.xaxis.set_minor_locator(xminor_locator)
        gl.xlocator = mticker.FixedLocator(xticks) # last working setting...
#        ymajor = np.arange(ds.lat.values.min(), ds.lat.values.max(), 1 )
#        yminor = [i for i in yticks if i not in ymajor]
        gl.ylocator = mticker.FixedLocator(yticks) # last working setting...
#        get_labels
#        ymajor_locator = mticker.FixedLocator(ymajor)
#        yminor_locator = mticker.FixedLocator(yminor)
#        ax.yaxis.set_major_locator(ymajor_locator)
#        ax.yaxis.set_minor_locator(yminor_locator)
        # tight off the main labels.
#        gl.xlabels_bottom = False
#        gl.ylabels_left = False
        gl.xlabels_bottom = True # last working setting...
        gl.ylabels_left = True # last working setting...
        # Add main labels
#        for lon in xmajor:
#            buffer = 0.25
#            ax.text( y0, lon, '{}'.format(lon), fontsize=10, alpha=0.5,
#                    horizontalalignment='center' )
#        for lat in ymajor:
#            buffer = 0.25
#            ax.text( lat, x0, '{}'.format(lon), fontsize=10, alpha=0.5,
#                    horizontalalignment='center' )
        # Make axis label text smaller
        gl.xlabel_style = {'size': 6, 'rotation' :90 }
        gl.ylabel_style = {'size': 6, }

    # Mark a known place to help us geo-locate ourselves
    if add_ARNA_locs:
        colours = AC.get_CB_color_cycle()
        locs2plot  = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
        for loc2plot in locs2plot:
            lon, lat, alt = AC.get_loc(loc2plot)
            ax.plot(lon, lat, 'bo', markersize=5, markerfacecolor='none',
                    markeredgewidth=2,
#                    markeredgecolor=colours[0],
                    markeredgecolor='black',
                    transform=ccrs.PlateCarree())
#            ax.text(lon, lat+0.25, loc2plot, transform=ccrs.PlateCarree())

    # Add a box to show the flyable range
    if add_flyable_range_as_box:
        # Get the minimum
        d = get_max_flying_range4BAE146()
        min_lon = d['min_lon']
        max_lon = d['max_lon']
        min_lat = d['min_lat']
        max_lat = d['max_lat']
        # Create new lists
        lons = [min_lon, min_lon, max_lon, max_lon]
        lats = [min_lat, max_lat, max_lat, min_lat]
        # Now plot as a linear ring
        ring = LinearRing(list(zip(lons, lats)))
        ax.add_geometries([ring], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='grey',
                          zorder=10, linestyle=':',
                          )
    if add_flyable_range_as_circle:
#        n_points = 1000
        # Approximate from James' max distance
        # ( 16.8331-13 ) *110667.45
        locs4circles = 'Dakar', 'Sao Vicente Airport',
        for loc in locs4circles:
            # Get locations to centre circle on
            lon, lat, alt = AC.get_loc(loc)
            # Radius in degrees
#            radius = 16.8331-13
            radius = 21 - 16.8331
            # Plot up circle
            ax.add_patch(mpatches.Circle(xy=[lon, lat],
                                         radius=radius,
#                                         color='red',
#                                         alpha=0.3,
                                         transform=projection,
                                         facecolor='none',
#                                         edgecolor='grey',
                                         edgecolor='black',
                                         linestyle=':',
                                         linewidth=3.0,
                                         zorder=100
                                         ))
    # Get limits of plotting data
    if isinstance(extents, type(None)):
        x0 = float(ds[LonVar].min())
        x1 = float(ds[LonVar].max())
        y0 = float(ds[LatVar].min())
        y1 = float(ds[LatVar].max())
        extents = (x0, x1, y0, y1)
    ax.set_extent(extents, crs=ccrs.PlateCarree())
    # Beautify the figure/plot
    if add_detailed_map:
        # Add borders for countries
        ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                       facecolor='none', zorder=50)
        # Also add minor islands (inc. Cape Verde)
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor=None,
                                                facecolor='none')
        ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)

    # Update the xaxis ticks (lon) to deg W
#    update_lon_units = True
#    if update_lon_units:
#         xticks = ax.get_xticks()
#         xticklabels = ax.get_xticklabels()
#         ax.set_xticklabels([str(i)*-1 for i in xticks])

        # Trun of cartopy axis
#        gl.xlabels_bottom = False
        #
        # Add main labels
#         for lon in np.arange(-27, -12, 3):
#             buffer = 0.25
#             ax.text(x0-buffer, lon, '{}'.format(lon), fontsize=10, alpha=0.5,
#                     horizontalalignment='center', transform=ccrs.PlateCarree())


    # Add the grid box with maximum NOy*Dust (print value)
    if (region == 'Cape_Verde_Flying') and add_max_vals_as_txt:
        # Sub select data of interest
        try:
            ds_tmp = ds[['NOy', 'Dust', 'NOy*Dust']]
            # Get extents
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            # Reduce the dataset to this size.
            print(xlim, ylim)
            bool1 = ((ds_tmp.lon >= xlim[0]) & (ds_tmp.lon <= xlim[1])).values
            bool2 = ((ds_tmp.lat >= ylim[0]) & (ds_tmp.lat <= ylim[1])).values
            # Cut by lon, then lat
            ds_tmp = ds_tmp.isel(lon=bool1)
            ds_tmp = ds_tmp.isel(lat=bool2)
            # Find the
            da = ds_tmp[ 'NOy*Dust']
            print( da.where(da==da.max(), drop=True).squeeze() )
            lon = da.where(da==da.max(), drop=True).squeeze().lon.values
            lat = da.where(da==da.max(), drop=True).squeeze().lat.values
            # Add a cross on the map.
            radius = ( 21 - 16.8331 ) / 4
            # Plot up circle
            projection = ccrs.PlateCarree()
            ax.add_patch(mpatches.Circle(xy=[lon, lat],
                                         radius=radius,
                                         transform=projection,
                                         facecolor='none',
                                         edgecolor='black',
                                         linestyle='-',
                                         linewidth=3.0,
                                         zorder=100
                                         ))
            # Add a label saying the location and values for NOy and Dust
            # Add label for the airport
#            xtext = (xlim[1]-xlim[0])/2
#            ytext = (ylim[1]-ylim[0])/2
#            print(xtext, ytext)
#            xtext, ytext = -22.2, 20.
            xtext, ytext = -16.2, 20.
#            print(xtext, ytext)
            #
            ds_tmp = ds_tmp.sel(lon=lon)
            ds_tmp = ds_tmp.sel(lat=lat)
            NOy_at_loc = float(ds_tmp['NOy'].values)
            Dust_at_loc = float(ds_tmp['Dust'].values)
            lon_nav_units = convert_decimcal_degress2nav_format([lon*-1])[0]
            lat_nav_units = convert_decimcal_degress2nav_format([lat])[0]
#            plot_txt = 'NOy*Dust max. @ {}N {}E'.format(lat_nav_units, lon_nav_units)
            plot_txt = 'NOy*Dust max. @ {}N {}W'.format(lat_nav_units,
                                                        lon_nav_units)
            # Also add the values for NOy+Dust
            units1 = 'ppbv'
            units2 = '$\mu$g m$^{-3}$'
            pstr = '\n (NOy={:.2f} {}, Dust={:.1f} {})'
            plot_txt += pstr.format(NOy_at_loc, units1, Dust_at_loc, units2)

            # Now add to plot
            ax.text(xtext, ytext, plot_txt,
#                     '{}'.format(loc_),
                    fontsize=10,
#                    alpha=0.5,
#                    horizontalalignment='center',
                    horizontalalignment='right',
                    transform=projection,
                    zorder=100,
                    )
        except:
            print(ds.data_vars)
            pstr = 'WARNING: not adding waypoint as NOy*Dust var not present'
            print(pstr)
    # Add a generic title if one is not provided
    if not isinstance(title, type(None)):
        plt.title(title)
    # Save the plot?
    if save_plot:
        if isinstance(savename, type(None)):
            savename = 'spatial_plot_{}_{}_{}_{}'
            savename = savename.format(region, var2plot1, var2plot2, extra_str)
        savename = AC.rm_spaces_and_chars_from_str(savename)
        if isinstance(folder, type(None)):
            folder = './'
        plt.savefig(folder+savename+'.png', dpi=dpi)
    if show_plot:
        plt.show()


def convert_decimcal_degress2nav_format(degrees_dec, just_degrees_minutes=True,
                                        LaTeX_format=True, debug=False ):
    """
    convert decimal degree to degrees, minutes and then decimal seconds
    """
    # Setup a list to store values
    vals2rtn = []
    # Loop by degrees provided
    for deg in degrees_dec:
        degrees = int(deg)
        degrees_frac = abs(AC.myround(deg, base=0.01, integer=False) - degrees)
        minutes = AC.myround(degrees_frac * 60,  base=1, integer=False)
        if just_degrees_minutes:
            pass
        else:
            print('TODO - calculation not setup!')
        if debug:
            print( 'Conversions for: {:.7f}'.format(deg))
            # Into decimal location (the same)
            pstr = 'Decimal format for location: {:>3}{:0>5}'
            print( pstr.format(degrees, degrees_frac) )
            # Into degrees, minutes and then decimal seconds
            if just_degrees_minutes:
                pstr = "Pilot format for location: {:>3}o {:0>2}'"
                print( pstr.format(degrees, minutes) )
                print('')
            else:
                pstr = "Pilot format for location: {:>3}o {:0>2}.{:0>2}'"
                print( pstr.format(degrees, minutes, seconds_dec) )
                print('')
        # Return in LaTeX form?
        if LaTeX_format:
            if just_degrees_minutes:
                degrees_str = str(degrees)
                minutes = '{:0>2}'.format(minutes)
                vals2rtn += [degrees_str+'$^{\circ}$'+"{}'".format(minutes)]
            else:
                print('TODO - calculation not setup!')
        else:
            if just_degrees_minutes:
                vals2rtn += ['{}{}'.forma(degrees, minutes)]
            else:
                print('TODO - calculation not setup!')
    return vals2rtn


def get_visibility_reports( dts=None, folder='./', debug=False  ):
    """
    Get the visibility reports from SDS-WAS
    """
    import wget
    # Which dates to use
    if isinstance(dts, type(None)):
        Tnow = AC.time2datetime( [gmtime()] )[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day,)
        # Use yesterday
        dt =  AC.add_days(dt, -1)
        # Use the last 18 days
        dts = [AC.add_days(dt, i*-1) for i in range(0, 5) ]
    # URL for SDS-WAS address
    URL_str = 'https://sds-was.aemet.es/archive/images/visibility/'
    URL_str += '{}/{:0>2}/images/{}{:0>2}{:0>2}_visibility.png'
    # for dt in dts
    for dt in dts:
        if debug:
            print(dt)
        URL = URL_str.format(dt.year, dt.month, dt.year, dt.month, dt.day )
        if debug:
            print(URL)
        filename = URL.split('/')[-1]
        if debug:
            print(filename)
        wget.download(URL, folder+filename)


def get_latest_GEOS5_diagnostics(dt=None,
                                 add_circles2GMAO_spatial_maps=True,
                                 verbose=False, debug=False):
    """
    Get the latest GEOS5 plots (datagrams, weather maps etc)
    """
    # Local variables
    prefix = 'ARNA'
    # If no day given, use previous day
    if isinstance(dt, type(None)):
        # Just use yesterday for Now
        TNow = AC.time2datetime( [gmtime()] )[0]
        dt = AC.add_days(TNow, -1)
        dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
    # Folder to dave plots
    G5_folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                        inc_collection=False)
    folder2save = G5_folder + '/plots/plots.GMAO/'
    # Check the folder exists and create if not
    # Check if the folder is present, if not create it
    if os.path.isdir(folder2save):
        if verbose:
            print('WARNING: folder already exists ({})'.format(folder2save))
    else:
        os.makedirs(folder2save)
        if verbose:
            print('Created folder for data ({})'.format(folder2save))

    #  - Get the latest datagram plots
    plts2get = [
    # Dust at Cape Verde
    'du_16.7_-23.0', 'dumass_16.7_-23.0',
    # Nitrate, CO, BC at Cape Verde
    'nimass_16.7_-23.0', 'bc_16.7_-23.0', 'bcmass_16.7_-23.0',
    # Total aerosol at Cape Verde
    'total_16.7_-23.0', 'co_16.7_-23.0',
    # Dust and total at Dakar
    'du_14.4_-17.0', 'dumass_14.4_-17.0', 'total_14.4_-17.0',
    ]
    AC.get_GEOS5_datagram_plots(folder=folder2save, dt=dt, plts2get=plts2get,
                                prefix=prefix)
    # - Get Weather maps
    ptype = 'wxmaps'
    taus = [i*24 for i in range(6)]
    taus += [0+(3*i) for i in range(3*9)] # Get 3 hourly values for the 1st 72 hours
    taus = list(sorted(set(taus)))
    field = 'precip'
    AC.get_GEOS5_online_diagnostic_plots(folder=folder2save, ptype=ptype,
                                         field=field, dt=dt, taus=taus,
                                         prefix=prefix)

    # - Get composition maps
    ptype = 'chem2d'
    taus = [i*24 for i in range(6)]
    taus += [0+(3*i) for i in range(3*9)] # Get 3 hourly values for the 1st 72 hours
    taus = list(sorted(set(taus)))
    fields = ['duaot', 'cobbaf', 'bcsmass', 'niaot', 'nismass', ]
    for field in fields:
        AC.get_GEOS5_online_diagnostic_plots(folder=folder2save, ptype=ptype,
                                             dt=dt, field=field, taus=taus,
                                             prefix=prefix)

    # - Add a circle over Dakar and Sao Vicente
    if add_circles2GMAO_spatial_maps:
        files = glob.glob( folder2save+'*.png')
        files = [i for i in files if ('chem2d' in i) or ('wxmaps' in i) ]
        for file in files:
            if debug:
                print('Adding Circles to file: {}'.format( file ) )
            # Open the image
            image = PIL.Image.open(file)
            draw = ImageDraw.Draw(image)
            # Add circles
            r = 39
            coords_list = [(774, 410), (840, 434) ]
            for coords in coords_list:
                x,y = coords
                # Calculate the locations
                leftUpPoint = (x-r, y-r)
                rightDownPoint = (x+r, y+r)
                twoPointList = [leftUpPoint, rightDownPoint]
                draw.ellipse(twoPointList, fill=None,
                             outline='Grey',
                             width=3 )
            # Save resulting image
            image.save( file )
    dt_str = dt.strftime('%Y/%m/%d %H:%M')
    print( 'Finished get_latest_GEOS5_diagnostics for {}'.format(dt_str) )


def quick_lon_plot_2layer(ds, var2plot1=None, var2plot2=None, extra_str='',
                          save_plot=True, show_plot=False, savename=None,
                          units=None, title=None, folder=None,
                          LatVar='lat', LonVar='lon', LevVar='lev',
                          fig=None, ax=None, extents=None,
                          add_ARNA_locs=True, use_local_CVAO_area=True,
                          extend='both', ylim=(0, 10), xlim=(5, 30),
                          dpi=320):
    """
    Plot up a quick longitude-altitude plot of data using cartopy

    Parameters
    -------
    ds (xr.Dataset): dataset object holding data to plot
    var2plot (str): variable to plot within the dataset
    LatVar, LonVar (str): variables to use for latitude and longitude
    save_plot (bool): save the plot as a .png ?
    show_plot (bool): show the plot on screen
    dpi (int): resolution to use for saved image (dots per square inch)
    savename (str): name to use for png of saved .png
    extra_str (str): extra string to append to save .png
    projection (cartopy.crs obj.):  projection to use
    fig (figure instance): matplotlib figure instance
    ax (axis instance): axis object to use
    add_ARNA_locs (bool):
    use_local_CVAO_area (bool):

    Returns
    -------
    (None)
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # Use the 1st data variable if not variable given
    if isinstance(var2plot1, type(None)):
        print('WARNING: No variable to plot (var2plot), trying 1st data_var')
#        var2plot1 = list(ds.data_vars)[0]
        var2plot1 = 'NOy'
    if isinstance(var2plot2, type(None)):
        print('WARNING: No variable to plot (var2plot), trying 1st data_var')
#        var2plot2 = list(ds.data_vars)[0]
        var2plot2 = 'PM2.5(dust)'
    # Setup figure and axis and plot
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=(12, 6))
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111)
    # Local variables (metres to kft)
    m2kft = 3.281
    # Convert lev to metres using hydrostatic aprox
    convert2km = True
    if convert2km:
        attrs = ds.lev.attrs
        print(attrs)
        if attrs['units'] == 'millibar':
            hPa_l = list(ds.lev.values).copy()
            ds.lev.values = AC.hPa_to_Km(hPa_l)
            attrs['units'] = 'km'
        else:
            pass
    # Then convert into kft?
    convert2kft = True
    if convert2kft:
        attrs = ds.lev.attrs
        print(attrs)
        if attrs['units'] == 'km':
#            hPa_l = list(ds.lev.values).copy()
            ds.lev.values = ds.lev.values * m2kft
            attrs['units'] = 'kft'
        else:
            pass
        #
        ylim = (ylim[0]*m2kft, ylim[1]*m2kft)

    # Setup devider axis for colourbar
    divider = make_axes_locatable(ax)
    # - Plot first var
    # Setup plotted range
    alpha = 0.5
    vmin1, vmax1 = set_limits4ar_plotted_range(var2plot1)
    vmin1 = get_vmin_value4var(var2plot1)
    units1 = get_species_units(var2plot1)
    cmap, ticks, nticks = get_cmap4var(var2plot1)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap': cmap, 'extend' : extend, 'pad':0.075, }
    else:
        cbar_kwargs  = {
        'ticks': ticks, 'cmap': cmap, 'extend' : extend, 'pad':0.075,
        }
    # Now plot up var1 - using pcolormesh
    cbar_ax = divider.append_axes("right", "2%", pad="1%")
    # Now plot
    ds[var2plot1].plot.pcolormesh(x=LatVar, y=LevVar, ax=ax,
                                 vmin=vmin1, vmax=vmax1,
                                 zorder=1, alpha=alpha,
#                                 origin='lower',
                                 yincrease=True,
                                 cmap=cmap,
#                                 extend=extend,
                                 cbar_ax=cbar_ax,
                                 cbar_kwargs=cbar_kwargs,
                                 )
    # Remove the title
    ax.set_title("")
    # - Plot second var
    # Setup a new axis for the colour bar
    cbar_ax = divider.append_axes("right", "2%", pad="7%")
    # Now plot up var 2
    alpha = 0.4
    vmin2, vmax2 = set_limits4ar_plotted_range(var2plot2)
    vmin2 = get_vmin_value4var(var2plot2)
    units2 = get_species_units(var2plot2)
    cmap, ticks, nticks = get_cmap4var(var2plot2)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap': cmap, 'extend' : 'both', }
    else:
        cbar_kwargs  = {'ticks': ticks, 'cmap': cmap, 'extend' : 'both', }
    # using pcolormesh
    ds[var2plot2].plot.pcolormesh(x=LatVar, y=LevVar, ax=ax,
                                 vmin=vmin2, vmax=vmax2,
                                 zorder=1, alpha=alpha,
                                 cmap=cmap,
#                                 origin='lower',
                                 yincrease=True,
                                 cbar_kwargs=cbar_kwargs,
                                 cbar_ax=cbar_ax,
                                 )
    # Remove the title
    ax.set_title("")
    # Add a box to show flyable area
    # 8km alt (and lon range)
    d = get_max_flying_range4BAE146()
    max_lat = d['max_lat']
    min_lat = d['min_lat']
    if convert2km and (not convert2kft):
        max_alt = d['max_alt'] /1E3
    elif convert2kft:
        max_alt = (d['max_alt'] / 1E3 * m2kft)
    else:
        pass
    max_alt_axis_coord = float(max_alt)/float(ylim[-1])
    print(max_alt_axis_coord)
    xrange = float(xlim[-1] - xlim[0])
    # Plot this up
    ax.axvline(x=max_lat, linestyle=':', alpha=0.5, color='grey', zorder=100,
               ymin=0, ymax=max_alt_axis_coord, linewidth=3.0)
    ax.axvline(x=min_lat, linestyle=':', alpha=0.5, color='grey', zorder=100,
               ymin=0, ymax=max_alt_axis_coord, linewidth=3.0)
    ax.axhline(y=max_alt, linestyle=':', alpha=0.5, color='grey', zorder=100,
               xmin=(min_lat - xlim[0] ) / xrange,
               xmax=(max_lat - xlim[0] ) / xrange,
               linewidth=3.0)

    # Add locations for airports
    locs2plot = [
    'Praia Airport', 'Dakar',  'Gran Canaria Airport', 'Sao Vicente Airport',
    'Lisbon Airport',  'Paris (Charles de Gaulle) Airport'
    ]
    for n, loc_ in enumerate( locs2plot ):
        lon_, lat_, NIU = AC.get_loc(loc_)
        # If lat in plotted range, then plot
        if (lat_ > xlim[0]) and (lat_ < xlim[-1]):
            ax.axvline(x=lat_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100, linewidth=3.0)
            # Add label for airports
            if n % 2 == 0:
                buffer = 0
            else:
                buffer = -0.5
            # Make a special case for Sao Vicente
            if loc_ == 'Sao Vicente Airport':
                buffer = -0.25
            # Set a base height for
            base = 9
            if convert2kft:
                base = base * m2kft
                buffer = buffer *3*1.5
            # Add label for the airport
            ax.text(lat_, base+buffer, '{}'.format(loc_), fontsize=10,
                    alpha=0.5,
                    horizontalalignment='center' )

    # Add lines for kft heights
    if convert2kft:
        hPa_heights = [1000, 900, 800, 700, 600, 500]
        km_heights = AC.hPa_to_Km(hPa_heights)
        kft_heights = [i*m2kft for i in km_heights]
        for n, height_ in enumerate( kft_heights ):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-2.5, height_, '{:.0f} hPa'.format(hPa_heights[n]),
                    fontsize=10,
                    alpha=0.5 )
    else:
        # Add lines for kft heights
        kft_heights = [20000, 15000, 10000, 5000]
        m_heights = [i/m2kft/1E3 for i in kft_heights]
        for n, height_ in enumerate( m_heights ):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-5, height_,
                    '{:.0f} kft'.format(kft_heights[n]/1E3),
                    fontsize=10, alpha=0.5 )

    # Beautify the figure/plot
    # Limit the yaxis
    ax.set_ylim(ylim)
    # Limit the xaxis
    ax.set_xlim(xlim)

    # Add hPA labels. - temporally remove .
#     press2plot = [1000, 900, 800, 700, 600, 500, 400]
#     ind = [ hPa_l.index(i) for i in press2plot ]
#     ax2.set_yticks(ds.lev[ind])
#     ax2.set_yticklabels([str(int(i)) for i in press2plot])
#     ax2.set_ylabel('Pressure [hPa]')

    # Add minor ticks for the x axis
#    ax.xaxis.grid(True, which='minor')

    # Add a generic title if one is not provided
    if not isinstance(title, type(None)):
#        ax.set_title(title)
        fig.suptitle(title)
    # Force tight layout
    plt.tight_layout(rect=(0, 0, 1, 0.95),)
#    plt.tight_layout()
    # Save the plot?
    if save_plot:
        if isinstance(savename, type(None)):
            sstr = 'spatial_plot_{}_{}_{}'
            savename = sstr.format(var2plot1, var2plot2, extra_str)
        savename = AC.rm_spaces_and_chars_from_str(savename)
        if isinstance(folder, type(None)):
            folder = './'
        plt.savefig(folder+savename+'.png', dpi=dpi)
    if show_plot:
        plt.show()


def get_coordinates_from_NetCDF_file(ds=None, folder=None, filename=None,
                                    falt_var='ALT_GIN',
                                    flat_var='LAT_GIN', flon_var='LON_GIN',
                                    AltVar='hPa', LonVar='lon', LatVar='lat',
                                    ftime_var='Time', TimeVar='time',
                                    convert_m2hPa=False, drop_NaNs=True):
    """
    Get locations (lat, lon, alt) from NetCDF files
    """
    import pandas as pd
    # Make a dataframne of locations from NEtCDF file
    if isinstance(ds, type(None)):
        ds = xr.open_dataset(folder+filename)
    df = pd.DataFrame()
    df[AltVar] = ds[falt_var].values
    df[LonVar] = ds[flon_var].values
    df[LatVar] = ds[flat_var].values
    df.index = ds[ftime_var].values
    # Convert metres of height to hPa
    if convert_m2hPa:
        df.loc[:,AltVar] = AC.hPa_to_Km(df[AltVar].values/1E3, reverse=True, )
    # Drop where there are not values for all coordinates
    if drop_NaNs:
        df = df.dropna()
    return df


def get_ARNA_flights_as_dfs():
    """
    Retrieve the ARNA flights as a list of dataframes
    """
    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    dfs = {}
    for flight_ID in flight_IDs:
        print( flight_ID )
        try:
            df = AC.get_FAAM_locations_as_df(flight_ID=flight_ID )
            dfs[flight_ID] = df
        except:
            print( 'WARNING: failed for {}'.format( flight_ID ) )
    return dfs


def plt_alt_binned_comparisons4ARNA_flights(dpi=320, show_plot=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context("paper", font_scale=0.75)
    # Which flights to plot
    flights_nums = [ 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # plot the altitude as a shadow on top of the plots
    plt_alt_as_shadow =  True

    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod = {}
    for flight_ID in flight_IDs:
        dfs_mod[flight_ID] = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        dfs_obs[flight_ID] = get_FAAM_core4flightnum(flight_ID=flight_ID )

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        df_mod = dfs_mod[flight_ID]

        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_altitude_binned_{}'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

        # - Plot up location of flights
        # Reset sns for spatial plot
#        sns.reset_orig()
        # New figure
        fig = plt.figure()
        # Get lat and lons
        lons = df_obs['LON_GIN'].values
        lats = df_obs['LAT_GIN'].values
        # Get dates of flight
        sdate_str = df_obs.index.min().strftime('%x %H:%M').strip()
        edate_str = df_obs.index.max().strftime('%x %H:%M').strip()
        # Make title
        title_str = 'Flight track for ARNA flight {} ({}-{})'
        title4plot = title_str.format(flight_ID, sdate_str, edate_str)
        #
        projection=ccrs.PlateCarree
        central_longitude = 0
        fig = plt.figure(dpi=dpi, facecolor='w', edgecolor='k')
        # Setup a cartopy projection for plotting
        ax = fig.add_subplot(111,
                             projection=projection(
                                 central_longitude=central_longitude)
                             )
        # Beautify
        ax.gridlines()
        # Add borders for countries
        ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                       facecolor='none', zorder=50)
        # Also add minor islands (inc. Cape Verde)
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor=None,
                                                facecolor='none')
        ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)
        # Plot settings
        marker='o'
        s=2
        alpha = 0.75
        cmap_list = AC.get_CB_color_cycle()
        # Now plot locations as scatter points on plot
        ax.scatter(lons, lats, color=cmap_list[0], s=s, marker=marker,
                   alpha=alpha,
                   label=flight_ID,
                   transform=projection(), zorder=999
                   )
        # Local area analysed as Cape Verde
        x0 = -30
        x1 =-10
        y0 = 0
        y1 = 25
        extents = (x0, x1, y0, y1)
        # Get limits of plotting data
        if isinstance(extents, type(None)):
            x0 = float(ds[LonVar].min())
            x1 = float(ds[LonVar].max())
            y0 = float(ds[LatVar].min())
            y1 = float(ds[LatVar].max())
            extents = (x0, x1, y0, y1)
        ax.set_extent(extents, crs=ccrs.PlateCarree())
        # Add a title to the plot
        plt.title(title4plot)
        plt.tight_layout()
        # Save to PDF
        if show_plot:
            plt.show()
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - put observations and vars to plot into a dictionary
        # Force alt to be in units of km
        ALT_var = 'Altitude (km)'
        Y_unit = ALT_var
        df_mod[ALT_var] = AC.hPa_to_Km( df_mod['model-lev'].values )
        df_obs[ALT_var] = df_obs['ALT_GIN'].values / 1E3
        #
        data_d = {'GEOS-CF': df_mod, 'Obs.':df_obs}

        # - Now plot up flight time series plots by variable
        title_str =  "Altitude binned '{}' ({}) during flight '{}'"
        # Setup color dictinoary
        color_dict = {'GEOS-CF': 'red', 'Obs.':'k'}
        unit_d = {}
        mod2obs_varnames = {
        'CO':'CO_AERO', 'O3':'O3_TECO', 'NO2':'no2_mr', 'NO':'no_mr',
        'HNO2':'hono_mr',
        'NOx':'NOx'
        }
        units_d = {
        'CO':'ppbv', 'O3':'ppbv', 'NO2':'pptv', 'NO':'pptv', 'NOx':'pptv',
        'HNO2':'pptv', 'HONO':'pptv',
        }
        range_d = {
        'CO':(50, 400), 'O3':(-10, 100), 'NO2':(-50, 500), 'NO':(-50, 500),
        'NOx':(-50, 500),
        'HNO2':(-50, 500), 'HONO':(-50, 500),
        }
        # - by variable
        runs = list(sorted(data_d.keys()))
		# Which variables to use?
        vars2plot = mod2obs_varnames.keys()
        print(vars2plot)
        print(df_obs.columns)
        vars2plot = [
        i for i in vars2plot if mod2obs_varnames[i] in df_obs.columns
        ]
        # What bins should be used?
        bins = [0.5*i for i in np.arange(15)]
        for var2plot in vars2plot:
            fig = plt.figure()
            ax = plt.gca()
            # Now loop data
            for n_key, key_ in enumerate(runs):
                print(n_key, key_, var2plot )
                #
                if key_ == 'Obs.':
                    varname = mod2obs_varnames[var2plot]
                else:
                    varname = var2plot
                # Setup an axis label
                units = units_d[var2plot]
                xlabel = '{} ({})'.format( var2plot, units )
                # Add alt to DataFrame
                df = pd.DataFrame({
                var2plot: data_d[key_][varname], ALT_var: data_d[key_][ALT_var]
                })
                #
                if key_ != 'Obs.':
                    scaleby = AC.get_unit_scaling(units)
                    df[var2plot] = df[var2plot].values * scaleby

                # drop any NaNs from the DataFrame
                s_shape = df.shape
                df.dropna(axis=0, how='any', inplace=True)
                if s_shape != df.shape:
                    pcent = (float(df.shape[0]) - s_shape[0])/s_shape[0] * 100.
                    pstr_dtr = 'WANRING dropped values - shape {}=>{} ({:.2f})'
                    print(pstr_dtr.format(s_shape, df.shape, pcent))
                # Plot up as binned boxplots using existing function
                try:
                    AC.binned_boxplots_by_altitude(df=df, fig=fig, ax=ax,
                                                   var2bin_by=ALT_var,
                                                   label=key_, xlabel=xlabel,
                                                   binned_var=var2plot,
                                                   num_of_datasets=len(runs),
                                                   bins=bins,
                                                   widths = 0.15,
                                                   dataset_num=n_key,
                                                   color=color_dict[key_])
                # Make NOx species be on a log scale
#                 if spec in NOx_specs:
#                     ax.set_xscale('log')
#                     ax.set_xlim( (1E-5, 1E3) )
#                 else:
                    ax.set_xscale('linear')
                except:
                    pass

                # Beautify plot
                plt.legend()
                plt.title(title_str.format(var2plot, units, flight_ID ))
                plt.xlim(range_d[var2plot])

            # Save to PDF
        #        fig.legend(loc='best', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
        #        plt.legend()
        #        plt.tight_layout()
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()

        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_timeseries_comparisons4ARNA_flights(dpi=320, show_plot=False):
    """
    Plot up timeseries comparisons between core observations and model data
    """
    import seaborn as sns
    # Now use Seaborn settings
    sns.set(color_codes=True)
    sns.set_context("paper", font_scale=0.75)
    # Which flights to plot
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
    flights_nums = [ 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # plot the altitude as a shadow on top of the plots
    plt_alt_as_shadow =  True
    aspect_ratio = 0.25
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod = {}
    for flight_ID in flight_IDs:
        dfs_mod[flight_ID] = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        dfs_obs[flight_ID] = get_FAAM_core4flightnum(flight_ID=flight_ID )

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        df_mod = dfs_mod[flight_ID]
        # get the begining and end of the flight from the extracted model times
        xylim_min = AC.add_minutes( df_mod.index.min(), -15)
        xylim_max = AC.add_minutes( df_mod.index.max(), 15 )
        xticks = df_mod.resample('15T' ).mean().index.values
        xticks = AC.dt64_2_dt( xticks )
        xticks_labels = [ i.strftime('%Y/%m/%d %H:%M') for i in xticks]

        # Setup a file
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_timeseries_flighttrack_{}'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

        # - Plot up location of flights
        # Reset sns for spatial plot
#        sns.reset_orig()
        # New figure
        fig = plt.figure()
        # Get lat and lons
        lons = df_obs['LON_GIN'].values
        lats = df_obs['LAT_GIN'].values
        # Get dates of flight
        sdate_str = df_obs.index.min().strftime('%x %H:%M').strip()
        edate_str = df_obs.index.max().strftime('%x %H:%M').strip()
        # Make title
        title_str = 'Flight track for ARNA flight {} ({}-{})'
        title4plot = title_str.format(flight_ID, sdate_str, edate_str)
        #
        projection=ccrs.PlateCarree
        central_longitude = 0
        fig = plt.figure(dpi=dpi, facecolor='w', edgecolor='k')
        # Setup a cartopy projection for plotting
        ax = fig.add_subplot(111,
                             projection=projection(
                                 central_longitude=central_longitude)
                             )
        # Beautify
        ax.gridlines()
        # Add borders for countries
        ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                       facecolor='none', zorder=50)
        # Also add minor islands (inc. Cape Verde)
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor=None,
                                                facecolor='none')
        ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)
        # Plot settings
        marker='o'
        s=2
        alpha = 0.75
        cmap_list = AC.get_CB_color_cycle()
        # Now plot locations as scatter points on plot
        ax.scatter(lons, lats, color=cmap_list[0], s=s, marker=marker,
                   alpha=alpha,
                   label=flight_ID,
                   transform=projection(), zorder=999
                   )
        # Local area analysed as Cape Verde
        x0 = -30
        x1 =-10
        y0 = 0
        y1 = 25
        extents = (x0, x1, y0, y1)
        # Get limits of plotting data
        if isinstance(extents, type(None)):
            x0 = float(ds[LonVar].min())
            x1 = float(ds[LonVar].max())
            y0 = float(ds[LatVar].min())
            y1 = float(ds[LatVar].max())
            extents = (x0, x1, y0, y1)
        ax.set_extent(extents, crs=ccrs.PlateCarree())
        # Add a title to the plot
        plt.title(title4plot)
        plt.tight_layout()
        # Save to PDF
        if show_plot:
            plt.show()
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Now plot up flight time series plots by variable
        title_str =  "Timeseries of '{}' ({}) during flight '{}'"
        # Now use Seaborn settings
        sns.set(color_codes=True)
        sns.set_context("paper", font_scale=0.75)

        # - Plot up carbon monoxide
        w, h = matplotlib.figure.figaspect(aspect_ratio)
        fig = plt.figure(figsize=(w, h))
        ax = fig.add_subplot(111)
        units = 'ppbv'
        var2plot = 'CO'
        obs_var2plot = 'CO_AERO'
        plt.plot(df_obs.index, df_obs[obs_var2plot].values, label='Obs.',
                 color='k' )
        mod_var2plot = 'CO'
        plt.plot(df_mod.index, df_mod[ mod_var2plot ].values*1E9,
                 label='GEOS-CF',
                 color='red' )
        # Beautify plot
        plt.title(title_str.format(var2plot, units, flight_ID ))
        plt.ylim(50, 400)
        plt.ylabel( '{} ({})'.format( var2plot, units) )
        plt.xlim(xylim_min, xylim_max)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels, rotation=45)
        print(xticks_labels)
        # Add a shadow of the altitude
        ax2 = ax.twinx()
        mod_var2plot = 'model-lev'
        # Invert the second y-axis
        if plt_alt_as_shadow:
            ax2.plot(df_mod.index, df_mod[ mod_var2plot ].values,
                     label='Altitude',
                     color='grey', zorder=100, alpha=0.25  )
            ax2.set_ylabel('Altitude (hPa)')
            ax2.grid(None)
            ax2.invert_yaxis()
            # Force use of the same ticks
            ax2.set_xticks(xticks)
            ax2.set_xticklabels(xticks_labels, rotation=45)

        # Save to PDF
        fig.legend(loc='best', bbox_to_anchor=(1,1),
                   bbox_transform=ax.transAxes)
        plt.tight_layout()
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()

        # - Plot up ozone
        w, h = matplotlib.figure.figaspect(aspect_ratio)
        fig = plt.figure(figsize=(w, h))
        #adjustFigAspect(fig, aspect=7.0)
        ax = fig.add_subplot(111)
#         xleft, xright = ax.get_xlim()
#         ybottom, ytop = ax.get_ylim()
#         ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*aspect_ratio)
        units = 'ppbv'
        var2plot = 'Ozone'
        obs_var2plot = 'O3_TECO'
        ln1 = plt.plot(df_obs.index, df_obs[obs_var2plot].values, label='Obs.',
                       color='k'  )
        mod_var2plot = 'O3'
        ln2 = plt.plot(df_mod.index, df_mod[ mod_var2plot ].values*1E9,
                       label='GEOS-CF', color='red'
                       )
        # Beautify plot
        title_str = "Timeseries of '{}' ({}) during flight '{}'"
        plt.title(title_str.format(var2plot, units, flight_ID ))
        plt.ylim(-10, 100)
        plt.ylabel( '{} ({})'.format( var2plot, units) )
        plt.xlim(xylim_min, xylim_max)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels, rotation=45)
        # Add a shadow of the altitude
        ax2 = ax.twinx()
        mod_var2plot = 'model-lev'
        # Invert the second y-axis
        if plt_alt_as_shadow:
            ax2.plot(df_mod.index, df_mod[ mod_var2plot ].values,
                     label='Altitude',
                      color='grey', zorder=100, alpha=0.25  )
            ax2.set_ylabel('Altitude (hPa)')
            ax2.grid(None)
            ax2.invert_yaxis()
            # Force use of the same ticks
            ax2.set_xticks(xticks)
            ax2.set_xticklabels(xticks_labels, rotation=45)

        # Save to PDF
        fig.legend(loc='best', bbox_to_anchor=(1,1),
                   bbox_transform=ax.transAxes)
        plt.tight_layout()
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()

        # - Plot up NO2
        try:
            # Setup for plot
            fig = plt.figure()
            w, h = matplotlib.figure.figaspect(aspect_ratio)
            fig = plt.figure(figsize=(w, h))
            #adjustFigAspect(fig, aspect=7.0)
            ax = fig.add_subplot(111)
            # Setup for specific variable
            units = 'pptv'
            var2plot = 'NO2'
            obs_var2plot = 'no2_mr'
            plt.plot(df_obs.index, df_obs[obs_var2plot].values, label='Obs.',
                     color='k' )
            mod_var2plot = 'NO2'
            plt.plot(df_mod.index, df_mod[ mod_var2plot ].values*1E12,
                     label='GEOS-CF',
                     color='red' )
            # Beautify plot
            plt.title(title_str.format(var2plot, units, flight_ID ))
            plt.ylim(-50, 500)
            plt.xlim(xylim_min, xylim_max)
            plt.ylabel( '{} ({})'.format( var2plot, units) )
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks_labels, rotation=45)
            # Add a shadow of the altitude
            ax2 = ax.twinx()
            mod_var2plot = 'model-lev'
            # Invert the second y-axis
            if plt_alt_as_shadow:
                ax2.plot(df_mod.index, df_mod[mod_var2plot].values,
                         label='Altitude',
                         color='grey', zorder=100, alpha=0.25  )
                ax2.set_ylabel('Altitude (hPa)')
                ax2.grid(None)
                ax2.invert_yaxis()
                # Force use of the same ticks
                ax2.set_xticks(xticks)
                ax2.set_xticklabels(xticks_labels, rotation=45)
            # Save to PDF
            fig.legend(loc='best', bbox_to_anchor=(1,1),
                       bbox_transform=ax.transAxes)
            plt.tight_layout()
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()
        except:
            print('Failed to plot NO2')

        # - Plot up NO
        try:
            # Setup for plot
            w, h = matplotlib.figure.figaspect(aspect_ratio)
            fig = plt.figure(figsize=(w, h))
            ax = fig.add_subplot(111)
            # Setup for specific variable
            units = 'pptv'
            var2plot = 'NO'
            obs_var2plot = 'no_mr'
            plt.plot(df_obs.index, df_obs[obs_var2plot].values, label='Obs.',
                     color='k' )
            mod_var2plot = 'NO'
            plt.plot(df_mod.index, df_mod[ mod_var2plot ].values*1E12,
                     label='GEOS-CF',
                     color='red' )
            # Beautify plot
            plt.title(title_str.format(var2plot, units, flight_ID ))
            plt.ylim(-50, 500)
            plt.xlim(xylim_min, xylim_max)
            plt.ylabel( '{} ({})'.format( var2plot, units) )
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks_labels, rotation=45)
            # Add a shadow of the altitude
            ax2 = ax.twinx()
            mod_var2plot = 'model-lev'
            # Invert the second y-axis
            if plt_alt_as_shadow:
                ax2.plot(df_mod.index, df_mod[mod_var2plot].values,
                         label='Altitude',
                         color='grey', zorder=100, alpha=0.25  )
                ax2.set_ylabel('Altitude (hPa)')
                ax2.grid(None)
                ax2.invert_yaxis()
                # Force use of the same ticks
                ax2.set_xticks(xticks)
                ax2.set_xticklabels(xticks_labels, rotation=45)
            # Save to PDF
            fig.legend(loc='best', bbox_to_anchor=(1,1),
                       bbox_transform=ax.transAxes)
            plt.tight_layout()
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()
        except:
            print('Failed to plot NO')

        # - Plot up NOx
        try:
            # Setup for plot
            w, h = matplotlib.figure.figaspect(aspect_ratio)
            fig = plt.figure(figsize=(w, h))
            ax = fig.add_subplot(111)
            # Setup for specific variable
            units = 'pptv'
            var2plot = 'NOx'
            obs_var2plot = 'NOx'
            plt.plot(df_obs.index, df_obs[obs_var2plot].values, label='Obs.',
                     color='k' )
            mod_var2plot = 'NOx'
            plt.plot(df_mod.index, df_mod[ mod_var2plot ].values*1E12,
                     label='GEOS-CF', color='red' )
            # Beautify plot
            plt.title(title_str.format(var2plot, units, flight_ID ))
            plt.ylim(-50, 500)
            plt.xlim(xylim_min, xylim_max)
            plt.ylabel( '{} ({})'.format( var2plot, units) )
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks_labels, rotation=45)
            # Add a shadow of the altitude
            ax2 = ax.twinx()
            mod_var2plot = 'model-lev'
            # Invert the second y-axis
            if plt_alt_as_shadow:
                ax2.plot(df_mod.index, df_mod[mod_var2plot].values,
                         label='Altitude',
                         color='grey', zorder=100, alpha=0.25  )
                ax2.set_ylabel('Altitude (hPa)')
                ax2.grid(None)
                ax2.invert_yaxis()
                # Force use of the same ticks
                ax2.set_xticks(xticks)
                ax2.set_xticklabels(xticks_labels, rotation=45)
            # Save to PDF
            fig.legend(loc='best', bbox_to_anchor=(1,1),
                       bbox_transform=ax.transAxes)
            plt.tight_layout()
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()
        except:
            print('Failed to plot NOx')

        # - Plot up HNO2
        try:
            # Setup for plot
            w, h = matplotlib.figure.figaspect(aspect_ratio)
            fig = plt.figure(figsize=(w, h))
            ax = fig.add_subplot(111)
            # Setup for specific variable
            units = 'pptv'
            var2plot = 'HONO'
            obs_var2plot = 'hono_mr'
            plt.plot(df_obs.index, df_obs[obs_var2plot].values, label='Obs.',
                     color='k' )
            mod_var2plot = 'HNO2'
            plt.plot(df_mod.index, df_mod[ mod_var2plot ].values*1E12,
                     label='GEOS-CF', color='red' )
            # Beautify plot
            plt.title(title_str.format(var2plot, units, flight_ID ))
            plt.ylim(-50, 500)
            plt.xlim(xylim_min, xylim_max)
            plt.ylabel( '{} ({})'.format( var2plot, units) )
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks_labels, rotation=45)
            # Add a shadow of the altitude
            ax2 = ax.twinx()
            mod_var2plot = 'model-lev'
            # Invert the second y-axis
            if plt_alt_as_shadow:
                ax2.plot(df_mod.index, df_mod[mod_var2plot].values,
                         label='Altitude', color='grey', zorder=100,
                         alpha=0.25  )
                ax2.set_ylabel('Altitude (hPa)')
                ax2.grid(None)
                ax2.invert_yaxis()
                # Force use of the same ticks
                ax2.set_xticks(xticks)
                ax2.set_xticklabels(xticks_labels, rotation=45)
            # Save to PDF
            fig.legend(loc='best', bbox_to_anchor=(1,1),
                       bbox_transform=ax.transAxes)
            plt.tight_layout()
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()
        except:
            print('Failed to plot NOx')

        # - Plot up temperature
        # Setup for plot
        w, h = matplotlib.figure.figaspect(aspect_ratio)
        fig = plt.figure(figsize=(w, h))
        #adjustFigAspect(fig, aspect=7.0)
        ax = fig.add_subplot(111)
        # Setup for specific variable
        units = '$^{\circ}$C'
        var2plot = 'Temperature'
        obs_var2plot = 'TAT_DI_R'
        plt.plot(df_obs.index, df_obs[ obs_var2plot ].values -273.15,
                 label='Obs.',
                 color='k' )
        mod_var2plot = 'T'
        plt.plot(df_mod.index, df_mod[mod_var2plot].values, label='GEOS-CF',
                 color='red'  )
        # Beautify plot
        plt.title(title_str.format(var2plot, units, flight_ID ))
        plt.ylim(-30, 30)
        plt.ylabel( '{} ({})'.format( var2plot, units) )
        plt.xlim(xylim_min, xylim_max)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels, rotation=45)
        # Add a shadow of the altitude
        ax2 = ax.twinx()
        mod_var2plot = 'model-lev'
        # Invert the second y-axis
        if plt_alt_as_shadow:
            ax2.plot(df_mod.index, df_mod[mod_var2plot].values,
                     label='Altitude',
                      color='grey', zorder=100, alpha=0.25  )
            ax2.set_ylabel('Altitude (hPa)')
            ax2.grid(None)
            ax2.invert_yaxis()
            # Force use of the same ticks
            ax2.set_xticks(xticks)
            ax2.set_xticklabels(xticks_labels, rotation=45)
        # Save to PDF
        fig.legend(loc='best', bbox_to_anchor=(1,1),
                   bbox_transform=ax.transAxes)
        plt.tight_layout()
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()

        # - Plot up Eastward wind
        # Setup for plot
        w, h = matplotlib.figure.figaspect(aspect_ratio)
        fig = plt.figure(figsize=(w, h))
        ax = fig.add_subplot(111)
        # Setup for specific variable
        units = 'm s$^{-1}$'
        var2plot = 'Eastward wind'
        obs_var2plot = 'U_C'
        plt.plot(df_obs.index, df_obs[obs_var2plot].values, label='Obs.',
                 color='k' )
        mod_var2plot = 'U'
        plt.plot(df_mod.index, df_mod[mod_var2plot].values, label='GEOS-CF',
                 color='red'  )
        # Beautify plot
        plt.title(title_str.format(var2plot, units, flight_ID ))
        plt.ylim(-25, 25)
        plt.ylabel( '{} ({})'.format( var2plot, units) )
        plt.xlim(xylim_min, xylim_max)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels, rotation=45)
        # Add a shadow of the altitude
        ax2 = ax.twinx()
        mod_var2plot = 'model-lev'
        # Invert the second y-axis
        if plt_alt_as_shadow:
            ax2.plot(df_mod.index, df_mod[mod_var2plot].values,
                     label='Altitude',
                      color='grey', zorder=100, alpha=0.25  )
            ax2.set_ylabel('Altitude (hPa)')
            ax2.grid(None)
            ax2.invert_yaxis()
            # Force use of the same ticks
            ax2.set_xticks(xticks)
            ax2.set_xticklabels(xticks_labels, rotation=45)
        # Save to PDF
        fig.legend(loc='best', bbox_to_anchor=(1,1),
                   bbox_transform=ax.transAxes)
        plt.tight_layout()
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()

        # - Plot up Eastward wind
        # Setup for plot
        w, h = matplotlib.figure.figaspect(aspect_ratio)
        fig = plt.figure(figsize=(w, h))
        ax = fig.add_subplot(111)
        # Setup for specific variable
        units = 'm s$^{-1}$'
        var2plot = 'Northward wind'
        obs_var2plot = 'V_C'
        mod_var2plot = 'V'
        plt.plot(df_obs.index, df_obs[obs_var2plot].values, label='Obs.',
                 color='k' )
        plt.plot(df_mod.index, df_mod[mod_var2plot].values, label='GEOS-CF',
                 color='red'  )
        # Beautify plot
        plt.legend()
        plt.title(title_str.format(var2plot, units, flight_ID ))
        plt.ylim(-25, 25)
        plt.ylabel( '{} ({})'.format( var2plot, units) )
        plt.xlim(xylim_min, xylim_max)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks_labels, rotation=45)
        # Add a shadow of the altitude
        ax2 = ax.twinx()
        mod_var2plot = 'model-lev'
        # Invert the second y-axis
        if plt_alt_as_shadow:
            ax2.plot(df_mod.index, df_mod[mod_var2plot].values,
                     label='Altitude',
                     color='grey', zorder=100, alpha=0.25  )
            ax2.set_ylabel('Altitude (hPa)')
            ax2.grid(None)
            ax2.invert_yaxis()
            # Force use of the same ticks
            ax2.set_xticks(xticks)
            ax2.set_xticklabels(xticks_labels, rotation=45)
        # Save to PDF
        fig.legend(loc='best', bbox_to_anchor=(1,1),
                   bbox_transform=ax.transAxes)
        plt.tight_layout()
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()

        # - Plot up Latitude
        try:
            # Setup for plot
            w, h = matplotlib.figure.figaspect(aspect_ratio)
            fig = plt.figure(figsize=(w, h))
            ax = fig.add_subplot(111)
            # Setup for specific variable
            units = '$^{\circ}$N'
            var2plot = 'Latitude'
            obs_var2plot = 'LAT_GIN'
            mod_var2plot = 'model-lat'
            plt.plot(df_obs.index, df_obs[obs_var2plot].values, label='Obs.',
                     color='k' )
            plt.plot(df_mod.index, df_mod[mod_var2plot].values,
                     label='GEOS-CF',
                     color='red'  )
            # Beautify plot
            plt.ylabel( '{} ({})'.format( var2plot, units) )
            plt.title(title_str.format(var2plot, units, flight_ID ))
            plt.xlim(xylim_min, xylim_max)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks_labels, rotation=45)
            # Add a shadow of the altitude
            ax2 = ax.twinx()
            mod_var2plot = 'model-lev'
            # Invert the second y-axis
            if plt_alt_as_shadow:
                ax2.plot(df_mod.index, df_mod[mod_var2plot].values,
                         label='Altitude',
                         color='grey', zorder=100, alpha=0.25  )
                ax2.set_ylabel('Altitude (hPa)')
                ax2.grid(None)
                ax2.invert_yaxis()
                # Force use of the same ticks
                ax2.set_xticks(xticks)
                ax2.set_xticklabels(xticks_labels, rotation=45)
            # Save to PDF
            fig.legend(loc='best', bbox_to_anchor=(1,1),
                       bbox_transform=ax.transAxes)
            plt.tight_layout()
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()
        except:
            print('Failed to plot Latitude')

        # - Plot up Longitude
        try:
            # Setup for plot
            w, h = matplotlib.figure.figaspect(aspect_ratio)
            fig = plt.figure(figsize=(w, h))
            ax = fig.add_subplot(111)
            # Setup for specific variable
            units = '$^{\circ}$E'
            var2plot = 'Longitude'
            obs_var2plot = 'LON_GIN'
            mod_var2plot = 'model-lon'
            plt.plot(df_obs.index, df_obs[obs_var2plot].values, label='Obs.',
                     color='k' )
            plt.plot(df_mod.index, df_mod[mod_var2plot].values,
                     label='GEOS-CF',
                     color='red'  )
            # Beautify plot
            plt.ylabel('{} ({})'.format( var2plot, units) )
            plt.title(title_str.format(var2plot, units, flight_ID ))
            plt.xlim(xylim_min, xylim_max)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks_labels, rotation=45)
            # Add a shadow of the altitude
            ax2 = ax.twinx()
            mod_var2plot = 'model-lev'
            # Invert the second y-axis
            if plt_alt_as_shadow:
                ax2.plot(df_mod.index, df_mod[mod_var2plot].values,
                         label='Altitude',
                         color='grey', zorder=100, alpha=0.25  )
                ax2.set_ylabel('Altitude (hPa)')
                ax2.grid(None)
                ax2.invert_yaxis()
                # Force use of the same ticks
                ax2.set_xticks(xticks)
                ax2.set_xticklabels(xticks_labels, rotation=45)
            # Save to PDF
            fig.legend(loc='best', bbox_to_anchor=(1,1),
                       bbox_transform=ax.transAxes)
            plt.tight_layout()
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()
        except:
            print('Failed to plot Longitude')

        # - Plot up altitude
        try:
            # Setup for plot
            w, h = matplotlib.figure.figaspect(aspect_ratio)
            fig = plt.figure(figsize=(w, h))
            ax = fig.add_subplot(111)
            # Setup for specific variable
            units = 'hPa'
            var2plot = 'Altitude'
            obs_var2plot = 'ALT_GIN'
            mod_var2plot = 'model-lev'
            # Local variables (metres to kft)
            vals = AC.hPa_to_Km( df_obs['ALT_GIN'].values/1E3, reverse=True )
            plt.plot(df_obs.index, vals, label='Obs.', color='k' )
            plt.plot(df_mod.index, df_mod[mod_var2plot].values,
                     label='GEOS-CF',
                     color='red'  )
            # Beautify plot
            plt.legend()
            plt.ylabel( '{} ({})'.format( var2plot, units) )
            plt.title(title_str.format(var2plot, units, flight_ID ))
            plt.xlim(xylim_min, xylim_max)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks_labels, rotation=45)
            # Invert the y-axis
            plt.gca().invert_yaxis()
            # Save to PDF
            plt.tight_layout()
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()
        except:
            print('Failed to plot altitude')

        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def adjustFigAspect(fig,aspect=1):
    '''
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    '''
    xsize,ysize = fig.get_size_inches()
    minsize = min(xsize,ysize)
    xlim = .4*minsize/xsize
    ylim = .4*minsize/ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(left=.5-xlim,
                        right=.5+xlim,
                        bottom=.5-ylim,
                        top=.5+ylim)


def get_FAAM_core4flightnum(flight_ID='C216' ):
    """
    Get the core FAAM flight data for a specific flight
    """
	# Where are the files? Which files to use?
    folder = '/users/ts551/scratch/data/ARNA/CEDA/v1/'
    files2use = glob.glob( folder + '*_{}*'.format(flight_ID.lower()) )
    # Open all the files together
    # Cannot as time is not capitalised in both
#    ds = xr.open_mfdataset(files2use)
    #
    FAAM_filename = [ i for i in files2use if 'core_faam' in i ]
    asstr = 'More than one core FAAM file found for flight!'
    assert len(FAAM_filename) == 1, asstr
    ds = xr.open_dataset(FAAM_filename[0])
    ds = ds.rename({'Time':'time'})
    # Check for the core NOx file
    Nitrates_filename = [ i for i in files2use if 'core-nitrates' in i ]
    if len(Nitrates_filename) >=1:
        ass_str = 'More than one nitrates file found for flight!'
        assert len(Nitrates_filename) == 1, ass_str
        ds2 = xr.open_dataset(Nitrates_filename[0])
        ds2 = ds2.mean(dim='sps10')
        # Merge the files
        ds = xr.merge([ds, ds2])
		# Get the HONO data too
        try:
            ds3 = xr.open_dataset(Nitrates_filename[0],
                                 group='non_core_group')
            ds3 = ds3.mean(dim='sps10')
            # Merge the files
            ds = xr.merge([ds, ds3])
        except OSError:
            print('WARNING: no HONO data found for {}'.format(flight_ID))
    # Convert to a dataframe
    df = ds.to_dataframe()
    # Add NOx as combined NO and NO2
    try:
        df['NOx'] = df['no_mr'].values + df['no2_mr'].values
    except:
        print('WARNING: failed to add NOx to flight {}'.format(flight_ID))
	# Resample the data?
    resample_data = True
    if resample_data:
        df = df.resample('1T' ).mean()
    return df


def get_GEOSCF_output4flightnum( flight_ID='C216' ):
    """
    Get the extracted GEOS-CF flight data for a specific flight
    """
    # Where is the extract GEOS-CF data?
    folder = '/users/ts551/scratch/data/ARNA/GEOS-CF/extracted_planeflight/'
    # Extract the data for a specific flight
    file2use = glob.glob( folder + '*_{}.csv'.format(flight_ID) )[0]
    df = pd.read_csv( file2use )
    # NOTE: this is a kludge for an early version of the file
#    df = df.T
#    new_header = df.iloc[0] #grab the first row for the header
#    df = df[1:] #take the data less the header row
#    df.columns = new_header #set the header row as the df header
    # Make the datatime the index and remove and unneeded columns
    df.index = df['Datetime'].values
    df.index = pd.to_datetime( df.index.values )
    # Add temperature in deg C
    df['TempK'] = df['T'].copy()
    df['T'] = df['T'].values - 273.15
    # Add NOx as combined NO and NO2
    df['NOx']  = df['NO'].values + df['NO2'].values
	# Resample the data?
    resample_data = True
    if resample_data:
        df = df.resample('1T' ).mean()
    return df


def extract_GEOS54all_ARNA_flights():
    """
    Extract GEOS-CF model data for all ARNA flights
    """
    # Which flights to plot
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
	# Just use non-transit ARNA flights
    flights_nums = [ 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # Loop by flight and extract the files
    for flight_ID in flight_IDs:
        # Extract data for flight
        df = extract_GEOS54ARNA_flight(flight_ID=flight_ID)
        # Save to csv.
        filename = 'ARNA_flightpath_extracted_from_GEOSCF_{}.csv'.format(flight_ID)
        folder = './'
        df.to_csv(folder+filename)


def extract_GEOS54ARNA_flight(flight_ID='C216'):
    """
    Extract ARNA flightpath from GEOS-CF assimilation data
    """
    # -  Get the measurement flight tracks
    # Manually set FAAM flight file to use for now...
    filename = 'core_faam_*_*_r*_{}_1hz.nc'.format(flight_ID.lower())
    folder = '/users/ts551/scratch/data/ARNA/CEDA/v1/'
    file2use = glob.glob(folder+filename)
    assert len(file2use) == 1, 'WARNING: more that one file found!'
    ds = xr.open_dataset( file2use[0] )
    #
    df = get_coordinates_from_NetCDF_file(ds=ds, convert_m2hPa=True)
    dt = AC.dt64_2_dt([df.index.values[0]])[0]
    # - Extract the flight tracks from the model
    #  Get the model data for the days near the flight
    dsGCF = get_GEOS_assim_expanded_dataset4ARNA(dt=dt)
    # Extract for known locations
    dfE = extract_ds4df_locs(ds=dsGCF, df=df,)
    return dfE


def extract_ds4df_locs(ds=None, df=None, LonVar='lon', LatVar='lat',
                       TimeVar='time',
                       AltVar='hPa', dsAltVar='hPa',
                       dsLonVar='lon', dsLatVar='lat', dsTimeVar='time',
                       vars2extract=None, debug=False, testing_mode=False):
    """
    Extract a xr.Dataset as for given locations (aka offline planeflight)
    """
    # Check the DataFrame index is datetime.datetime
    ass_str = 'WARNING: DataFrame index is not of datetime.datetime dtype'
#    assert df.index.dtype == datetime.datetime, ass_str
    # Extract all of the data variables unless a specific list is provided
    if isinstance(vars2extract, type(None)):
        vars2extract = list(ds.data_vars)
    # - Create a data frame for values
    dfN = pd.DataFrame()
    # - Loop via time
    loop_via_time = False
    if loop_via_time:
        if testing_mode:
            times2use = df.index.values[:10]
        else:
            times2use = df.index.values
        for dftime in times2use:
            # get the times for a specific data
            df_tmp = df.loc[df.index == dftime, :]
            if debug:
                print( df_tmp.shape )
            lat_ = df_tmp[LatVar].values
            hPa_ = df_tmp[AltVar].values
            lon_ = df_tmp[LonVar].values
            #
            ds_tmp = ds.sel(lon=lon_, lat=lat_, lev=hPa_, time=dftime,
                            method='nearest')

            #
            for data_var in vars2extract:
                dfN.loc[dftime, data_var ] = ds_tmp[data_var].values
            #
            del ds_tmp, df_tmp
    else:
        # get indexes =en masse then extract with these
        d = calculate_idx2extract(ds=ds, df=df)
        #
        if testing_mode:
            times2use = df.index.values[:10]
        else:
            times2use =  df.index.values
        # Loop by timestamp
        for n, time in enumerate( times2use ):
            # get the times for a specific data
            lat_idx = d['lat'][n]
            lon_idx = d['lon'][n]
            lev_idx = d['hPa'][n]
            time_idx =  d['time'][n]
            #
            ds_tmp = ds.isel(lat=lat_idx, lon=lon_idx, time=time_idx,
                             lev=lev_idx)
#            d_tmp = ds_tmp.to_dict()
#            vals = ds_tmp.to_array().values
            vals = [ ds_tmp[i].data for i in vars2extract ]
            vals = np.array(vals)
            for nval, val in enumerate(vals):
                dfN.loc[vars2extract[nval], time] = vals[nval]
            # Add the model position coordinates...
            dfN.loc['model-lat', time] = float(ds_tmp['lat'].values)
            dfN.loc['model-lon', time] = float(ds_tmp['lon'].values)
            dfN.loc['model-lev', time] = float(ds_tmp['lev'].values)
            dfN.loc['model-time', time] = float(ds_tmp['time'].values)
            del ds_tmp, vals
#        # Add the tracer names to the index
#        dfN.index = vars2extract
        # Make datetime the index
        dfN = dfN.transpose()
        # Save the datetime as a column too
        dfN['Datetime'] = dfN.index.values
        # Update the model datetime to be in datetime units
        dfN['model-time'] = pd.to_datetime(dfN['model-time'].values)
    return dfN


def calculate_idx2extract(ds=None, df=None, LonVar='lon', LatVar='lat',
                          TimeVar='time',
                          AltVar='hPa', dsAltVar='lev',
                          dsLonVar='lon', dsLatVar='lat', dsTimeVar='time',
                          debug=False, testing_mode=False):
    """
    Calculated the indexes to extract of a dataset from a dataframe.
    """
    # Get arrays of the coordinate variables in the dataset
    ds_lat = ds[dsLatVar].values
    ds_lon = ds[dsLonVar].values
    ds_hPa = ds[dsAltVar].values
    ds_time = ds[dsTimeVar].values
    # Calculate the index individually by coordinate
    lat_idx = [AC.find_nearest(ds_lat, i) for i in df[LatVar].values]
    lon_idx = [AC.find_nearest(ds_lon, i) for i in df[LonVar].values]
    hPa_idx = [AC.find_nearest(ds_hPa, i) for i in df[AltVar].values]
    time_idx = [AC.find_nearest(ds_time, i) for i in df.index.values]
    # Return a dictionary of the values
    return {LatVar:lat_idx, LonVar:lon_idx, TimeVar:time_idx, AltVar:hPa_idx}


def get_GEOS_assim_expanded_dataset4ARNA( dt=None, update_lvl_units=True ):
    """
    Get GEOS-CF assimilation data (expanded) for the ARNA campaign
    """
    # Where is the expanded GEOS-CF data
    NASA_data = get_local_folder('NASA_data')
    folder = NASA_data + 'GEOS_CF/ARNA/assim/expanded_variable_files/'
    file_str = '*.nc4'
    # Make a dataframe of available files
    files = list(sorted(glob.glob(folder+file_str)))
    dates = [i.split('.nc4')[0][-14:] for i in files ]
    dts = [datetime.datetime.strptime(i, '%Y%m%d_%H%Mz') for i in dates]
    df = pd.DataFrame({'files':files}, index=dts)
    # If a date is given, only open dates a day before and after
    if not isinstance(dt, type(None)):
        # Use one dat before and after
        sdate = AC.add_days(dt, -1)
        edate = AC.add_days(dt, 1)
        # Limit data to one day before or after
        bool1 = (df.index >=sdate) & (df.index <= edate)
        df = df.loc[bool1, :]
    # get a list of those to extract
    files = df['files'].values.astype(list)
    #
    ds = xr.open_mfdataset(files)
    #
    # Update the units for lev
    update_lvl_units = True
    if update_lvl_units:
        HPa_l = get_GEOSCF_vertical_levels(native_levels=True)
        attrs = ds.lev.attrs.copy()
        attrs['units'] = 'hPa'
        attrs['standard_name'] = 'Pressure'
        attrs['long_name'] = 'Pressure'
        ds.lev.attrs = attrs
        ds.lev.values = [ HPa_l[int(i)] for i in ds.lev.values -1]
    # Return the updated ds
    return ds


def plot_up_flight_locations_from_FAAM_website(d=None,
                                               LatVar='latitude',
                                               LonVar='longitude',
                                               folder='./'):
    """
    Plot up flight locations for all ARNA flights on a map
    """
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    # Get flight locations as a dataframe
    if isinstance(d, type(None)):
        d = get_ARNA_flights_as_dfs()
    # Setup plot
    projection=ccrs.PlateCarree
    central_longitude = 0
    fig = plt.figure(dpi=dpi, facecolor='w', edgecolor='k')
    # Setup a cartopy projection for plotting
    ax = fig.add_subplot(111,
                         projection=projection(
                             central_longitude=central_longitude)
                         )
    # Beautify
    ax.gridlines()
    # Add borders for countries
    ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                   facecolor='none', zorder=50)
    # Also add minor islands (inc. Cape Verde)
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                            edgecolor=None,
                                            facecolor='none')
    ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)
    # Plot settings
    marker='o'
    s=2
    alpha = 0.75
    cmap_list = AC.get_CB_color_cycle()
    #
    flight_IDs = list(sorted(d.keys()))
    for flight_ID_n, flight_ID in enumerate( flight_IDs ):
        df = d[flight_ID]
        # Save the data as a csv file?
        save2csv = True
        if save2csv:
            df.to_csv('TEMP_{}.csv'.format(flight_ID) )
        # Get lats and lons for flight
        lats = df[LatVar]
        lons = df[LonVar]
        # Now plot locations as scatter points on plot
        ax.scatter(lons, lats, color=cmap_list[flight_ID_n], s=s,
                   marker=marker,
                   alpha=alpha,
                   label=flight_ID,
                   transform=projection(), zorder=999
                   )
    # Local area analysed as Cape Verde
    x0 = -30
    x1 =-10
    y0 = 0
    y1 = 25
    extents = (x0, x1, y0, y1)
    # Get limits of plotting data
    if isinstance(extents, type(None)):
        x0 = float(ds[LonVar].min())
        x1 = float(ds[LonVar].max())
        y0 = float(ds[LatVar].min())
        y1 = float(ds[LatVar].max())
        extents = (x0, x1, y0, y1)
    ax.set_extent(extents, crs=ccrs.PlateCarree())
    # Include a legend
    plt.legend(ncol=3, loc="lower left", fontsize='small', handletextpad=0.05)
    # Add a tiutle
    plt.suptitle( u'ARNA campaign flights (excluding transits)')
    # Save the plot
    filename = 'ARNA_campaign_flight_locations_ALL'
    plt.savefig(folder+filename, pad_inches=0.25,
                bbox_inches='tight'
                )
    plt.close()


def quick_lat_plt_2layer(ds, var2plot1=None, var2plot2=None, extra_str='',
                         save_plot=True, show_plot=False, savename=None,
                         units=None, title=None, folder=None,
                         LatVar='lat', LonVar='lon', LevVar='lev',
                         fig=None, ax=None, extents=None,
                         add_ARNA_locs=True, use_local_CVAO_area=True,
                         extend='both', ylim=(0, 10), xlim=(-30, -16.5),
                         dpi=320):
    """
    Plot up a quick latitude-altitude plot of data using cartopy

    Parameters
    -------
    ds (xr.Dataset): dataset object holding data to plot
    var2plot (str): variable to plot within the dataset
    LatVar, LonVar (str): variables to use for latitude and longitude
    save_plot (bool): save the plot as a .png ?
    show_plot (bool): show the plot on screen
    dpi (int): resolution to use for saved image (dots per square inch)
    savename (str): name to use for png of saved .png
    extra_str (str): extra string to append to save .png
    projection (cartopy.crs obj.):  projection to use
    fig (figure instance): matplotlib figure instance
    ax (axis instance): axis object to use
    add_ARNA_locs (bool):
    use_local_CVAO_area (bool):

    Returns
    -------
    (None)
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # Use the 1st data variable if not variable given
    if isinstance(var2plot1, type(None)):
        print('WARNING: No variable to plot (var2plot), trying 1st data_var')
        var2plot1 = list(ds.data_vars)[0]
        var2plot1 = 'NOy'
    if isinstance(var2plot2, type(None)):
        print('WARNING: No variable to plot (var2plot), trying 1st data_var')
#        var2plot2 = list(ds.data_vars)[0]
        var2plot2 = 'PM2.5(dust)'
    # Setup figure and axis and plot
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=(12, 6))
#        fig = plt.figure()
    if isinstance(ax, type(None)):
#        ax = fig.add_subplot(111, projection=projection, aspect='auto')
        ax = fig.add_subplot(111)
#    if isinstance(xlim, type(None)):
#        xlim = (ds.lat.values.min(), ds.lat.values.max())
    # Local variables (metres to kft)
    m2kft = 3.281
    # Add a second y axis - tempoaroiliy hash out
#    ax2 = ax.twinx()
#        ax.set_yscale('log')
#        ax.set_xscale('linear')
    # Convert lev to metres using hydrostatic aprox
    convert2km = True
    if convert2km:
        attrs = ds.lev.attrs
        print(attrs)
        if attrs['units'] == 'millibar':
            hPa_l = list(ds.lev.values).copy()
            ds.lev.values = AC.hPa_to_Km(hPa_l)
            attrs['units'] = 'km'
        else:
            pass
    # Then convert into kft?
    convert2kft = True
    if convert2kft:
        attrs = ds.lev.attrs
        print(attrs)
        if attrs['units'] == 'km':
#            hPa_l = list(ds.lev.values).copy()
            ds.lev.values = ds.lev.values * m2kft
            attrs['units'] = 'kft'
        else:
            pass
        #
        ylim = (ylim[0]*m2kft, ylim[1]*m2kft)

    # Divider axis for colourbar
    divider = make_axes_locatable(ax)
    # - Plot first var
    #
#    cbar_kwargs = {  'extend': extend }
    # Setup plotted range
    alpha = 0.5
    vmin1, vmax1 = set_limits4ar_plotted_range(var2plot1)
    vmin1 = get_vmin_value4var(var2plot1)
    units1 = get_species_units(var2plot1)
    cmap, ticks, nticks = get_cmap4var(var2plot1)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap': cmap, 'extend' : 'both', 'pad':0.075,}
    else:
        cbar_kwargs  = {
        'ticks': ticks, 'cmap': cmap, 'extend' : 'both','pad':0.075,
        }
    # Now plot up var1 - using pcolormesh
    cbar_ax = divider.append_axes("right", "2%", pad="1%")
    # Now plot
    ds[var2plot1].plot.pcolormesh(x=LonVar, y=LevVar, ax=ax,
                                 vmin=vmin1, vmax=vmax1,
                                 zorder=1, alpha=alpha,
#                                 origin='lower',
                                 yincrease=True,
                                 cmap=cmap,
#                                 extend=extend,
                                 cbar_ax=cbar_ax,
                                 cbar_kwargs=cbar_kwargs,
                                 )
    # Remove the title
    ax.set_title("")

    # - Plot second var
    # Now plot up var 2
    alpha = 0.4
    vmin2, vmax2 = set_limits4ar_plotted_range(var2plot2)
    vmin2 = get_vmin_value4var(var2plot2)
    units2 = get_species_units(var2plot2)
    cmap, ticks, nticks = get_cmap4var(var2plot2)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap': cmap, 'extend' : 'both', }
    else:
        cbar_kwargs  = {'ticks': ticks, 'cmap': cmap, 'extend' : 'both', }
    # Setup a colour bar axis
    cbar_ax = divider.append_axes("right", "2%", pad="7%")
    # using pcolormesh
    im = ds[var2plot2].plot.pcolormesh(x=LonVar, y=LevVar, ax=ax,
                                 vmin=vmin2, vmax=vmax2,
                                 zorder=1, alpha=alpha, cmap=cmap,
#                                 origin='lower',
                                 yincrease=True,
                                 cbar_kwargs=cbar_kwargs,
                                 cbar_ax=cbar_ax,
                                 )
    # Remove the title
    ax.set_title("")

    # Add a box to show flyable area
    # 8km alt (and lon range)
    d = get_max_flying_range4BAE146()
#    max_lat = d['max_lat']
#    min_lat = d['min_lat']
    max_lon = d['max_lon']
    min_lon = d['min_lon']
    if convert2km and (not convert2kft):
        max_alt = d['max_alt'] /1E3
    elif convert2kft:
        max_alt = (d['max_alt'] / 1E3 * m2kft)
    else:
        pass
#    print(max_lat, min_lat, max_alt  )
    max_alt_axis_coord = float(max_alt)/float(ylim[-1])
    print(max_alt_axis_coord)
    xrange = float(xlim[-1] - xlim[0])
    # Plot this up
    ax.axvline(x=max_lon, linestyle=':', alpha=0.5, color='grey', zorder=100,
               ymin=0, ymax=max_alt_axis_coord, linewidth=3.0)
    ax.axvline(x=min_lon, linestyle=':', alpha=0.5, color='grey', zorder=100,
               ymin=0, ymax=max_alt_axis_coord, linewidth=3.0)
    ax.axhline(y=max_alt, linestyle=':', alpha=0.5, color='grey', zorder=100,
#               xmin=min_lon/xrange, xmax=max_lon/xrange,
               xmin=(min_lon - xlim[0] ) / xrange,
               xmax=(max_lon - xlim[0] ) / xrange,
               linewidth=3.0)

    # Add locations for airports
#    locs2plot = ['DSS', 'RAI', 'VXE', 'LPA', 'LIS', 'CDG']
    locs2plot = [
    'Praia Airport',  'Sao Vicente Airport', 'Dakar',
#    'Gran Canaria Airport', 'Lisbon Airport',  'Paris (Charles de Gaulle) Airport'
    ]
    for n, loc_ in enumerate( locs2plot ):
        lon_, lat_, NIU = AC.get_loc(loc_)
        # If lat in plotted range, then plot
        if (lon_ > xlim[0]) and (lon_ < xlim[-1]):
    #        print(lat_)
            ax.axvline(x=lon_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=3.0)
            # Add label for airports
    #        if loc_ in ['RAI','DSS' ]:
            if n % 2 == 0:
                buffer = 0
            else:
                buffer = -0.5
            # Make a special case for Sao Vicente
            if loc_ == 'Sao Vicente Airport':
                buffer = -0.25

            # Set a base height for
            base = 9
            if convert2kft:
                base = base * m2kft
                buffer = buffer *3*1.5

            # Add label for the airport
            ax.text(lon_, base+buffer, '{}'.format(loc_), fontsize=10,
                    alpha=0.5,
                    horizontalalignment='center' )

    #        ax.annotate(loc_, xy=(lat_, 5), xytext=(lat_+buffer, 5+2),
    #            arrowprops=dict(facecolor='black', shrink=0.05))

    # Add lines for kft heights
    if convert2kft:
        hPa_heights = [1000, 900, 800, 700, 600, 500]
        km_heights = AC.hPa_to_Km(hPa_heights)
        kft_heights = [i*m2kft for i in km_heights]
        for n, height_ in enumerate( kft_heights ):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-3.5, height_, '{:.0f} hPa'.format(hPa_heights[n]),
                    fontsize=10,
                    alpha=0.5 )
    else:
        # Add lines for kft heights
        kft_heights = [20000, 15000, 10000, 5000]
        m_heights = [i/m2kft/1E3 for i in kft_heights]
        for n, height_ in enumerate( m_heights ):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-5, height_,
                    '{:.0f} kft'.format(kft_heights[n]/1E3),
                    fontsize=10, alpha=0.5 )

    # Beautify the figure/plot
    # Limit the yaxis
    ax.set_ylim(ylim)
    # Limit the xaxis
    ax.set_xlim(xlim)

    # Add hPA labels. - temporally remove .
#     press2plot = [1000, 900, 800, 700, 600, 500, 400]
#     ind = [ hPa_l.index(i) for i in press2plot ]
#     ax2.set_yticks(ds.lev[ind])
#     ax2.set_yticklabels([str(int(i)) for i in press2plot])
#     ax2.set_ylabel('Pressure [hPa]')

    # Add minor ticks for the x axis
#    ax.xaxis.grid(True, which='minor')

    # Update the lon axis labels
#    update_lon_units = True
#    if update_lon_units:
#        ax.tick_params(axis='x', which='both', labelbottom='off')
#        im.tick_params(axis='x', which='both', labelbottom='off')
#        pass

    # Add a generic title if one is not provided
    if not isinstance(title, type(None)):
#        ax.set_title(title)
        fig.suptitle(title)
    # Force tight layout
    plt.tight_layout(rect=(0, 0, 1, 0.95),)
#    plt.tight_layout()
    # Save the plot?
    if save_plot:
        if isinstance(savename, type(None)):
            sstr = 'spatial_plot_{}_{}_{}'
            savename = sstr.format(var2plot1, var2plot2, extra_str)
        savename = AC.rm_spaces_and_chars_from_str(savename)
        if isinstance(folder, type(None)):
            folder = './'
        plt.savefig(folder+savename+'.png', dpi=dpi)
    if show_plot:
        plt.show()


def get_GEOSCF_data_cubes4collection(ds=None, region='Cape_Verde',
                                     doys2use=None,
                                     limit_lvls=True, limit_lons=False,
                                     limit_lats=True,
                                     vars2use=None, dt=None,
                                     collection = 'chm_inst_1hr_g1440x721_p23',
                                     rm_existing_file=False,
                                     dates2use=None, mode='fcast',
                                     folder=None):
    """
    Extract cubes of data for region of interest during campaign period
    """
    # Which variables to use
    if isinstance(vars2use, type(None)):
    #    vars2use = ("pm25du_rh35_gcc", "pm25_rh35_gcc", "noy", "co", "no2", "o3",)
        vars2use = list(convert_GEOSCF_var2GEOSChem_name(rtn_dict=True).keys())
    #    vars2use = ('pm25su_rh35_gcc', 'pm25ni_rh35_gcc', 'pm25soa_rh35_gc', 'pm25ss_rh35_gcc', 'so2')

    # - Get the data as dataset via OPeNDAP
    if not isinstance(ds, type(None)):
#        ds = C2BR2F4CF_as_ds_via_OPeNDAP(collection=collection, mode=mode, date=dt)
        # Temporarily do not allow passing of date to OPeNDAP fetcher
        ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode,
                                             date=None)

    # - Only consider dates that are in dates2use list (if provided)
    # Use the datetime components via pandas as not yet fully in xarray.
    if not isinstance(dates2use, type(None)):
        df = pd.DataFrame(index=ds['time'].values)
        df['date'] = df.index.date
        date_bool = df['date'].isin(dates2use).values
        ds = ds.isel(time=date_bool)
    # Check that the initial retrieved data is the same as the requested on
    dt_str = dt.strftime('%Y/%m/%d %H:%M')
    t0_CF = ds.time.values[0]
    t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
    t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
    if dt_str != t0_CF_str:
        pstr = 'WARNING: provided date ({}) is not 1st date ({}) in dataset!'
        print( pstr.format( dt_str, t0_CF_str ) )
        print( 'WARNING: this data is being saved here: {}'.format(folder) )

    # - Subset the dataset by region
    if region == 'Cape_Verde':
        # only consider the region the the plane will fly through +2 days transport
        # 20W-55E, 40S-40N,
        bool1 = ((ds.lon >= -35) & (ds.lon <= 10)).values
        bool2 = ((ds.lat >= 0) & (ds.lat <= 60)).values
        # Cut by lon, then lat
        ds = ds.isel(lon=bool1)
        ds = ds.isel(lat=bool2)
    elif region == 'global':
        pass
    else:
        print('WARNING: exiting as no region set')
        sys.exit()

    # - Limit the vertical levels
    if limit_lvls:
#        lvls = [1000, 850, 500]
        lvls = [1000, 900, 800, 700, 600, 500]
        bool_lvls = [i in lvls for i in ds['lev']]
        ds = ds.isel(lev=bool_lvls)
        extr_str = 'lvls_'+'_'.join([str(i) for i in lvls])
    else:
        extr_str = 'ALL_lvls'

    # - Limit the number of longitudes?
    if limit_lons:
        lons2check = [-18 ,-19.5, -21, -22.5, -24, -25.5]
        idx = [ AC.find_nearest(ds.lon.values, i) for i in  lons2check]
        lons2use = ds.lon.values[ idx ]
        bool_lons = [i in lons2use for i in ds['lon'].values]
        ds = ds.isel(lon=bool_lons)
        extr_str += '_lons_'+'_'.join([str(i) for i in lons2use])
    else:
        pass
    # - Limit the number of latitudes?
    if limit_lats:
        lats2check = [12, 13, 14, 15, 16, 17]
        idx = [ AC.find_nearest(ds.lat.values, i) for i in lats2check]
        lats2use = ds.lat.values[ idx ]
        bool_lats = [i in lats2use for i in ds['lat'].values]
        ds = ds.isel(lat=bool_lats)
        extr_str += '_lats_'+'_'.join([str(i) for i in lats2use])

    # Setup a print strin
    if limit_lons and limit_lvls:
        pstr = 'sliced by lvls+lons'
    elif limit_lvls:
        pstr = 'sliced by lvls'
    elif limit_lons:
        pstr = 'sliced by lons'
    elif limit_lats:
        pstr = 'sliced by lats'
    else:
        pstr =''

    # -  Save the data by variable and day of year reduce packet size of transfers
    # Where to save?
    if isinstance(folder, type(None)):
        folder = './'
    # Save format for data (region, year, doy, variable)
    filestr = 'ARNA_GEOSCF_{}_{}_{}_{:0>3}_{}_{}.nc'
    if isinstance(doys2use, type(None)):
        doys2use = list(sorted(set(ds['time.dayofyear'].values)))
    # Loop by var and save
    for var in vars2use:
        # Now loop and save files in single day units
        print("Getting data via OPeNDAP for var '{}' - {}".format(var, pstr))
        for doy in doys2use:
            doy_int = int(float(doy))
            pstr = "Getting data via OPeNDAP for doy '{}' ({})"
            print(pstr.format(doy_int, var))
            # Subset for a single day
            ds_tmp = ds[var].sel(time=ds['time.dayofyear'] == doy)
            # If there is data in the dataset, then download then download this
            if len(ds_tmp['time'].values) > 0:
                # Get the year for the doy
                year = list(set(ds_tmp['time.year'].values))[0]
                # Set the name of the file to save
                fname = filestr.format(collection, region, year, doy_int, var,
                                       extr_str)
                # now save this single file
                if rm_existing_file:
                    AC.rm_file(folder=folder, filename=fname)
                ds_tmp.to_netcdf(folder+fname)
                # Do some memory management
                del ds_tmp
                gc.collect()
            else:
                pstr = "WARNING: no data for '{}' on doy ({}) - CALLING STOP!"
                print(pstr.format(var, doy_int))
                pstr = 'Please check the correct doys are downloaded for {}'
                print( pstr.format( dt_str ) )
                print( 'doys downloading: ', doys2use )
                sys.exit()


def get_GEOS5_data_cubes4collection(ds=None, region='Cape_Verde',
                                    doys2use=None, mode='fcast',
                                    vars2use=None,
                                    collection='inst3_3d_aer_Nv',
                                    folder=None, dates2use=None,
#                                    limit_lvls=False, limit_lons=False
                                    rm_existing_file=False,
                                    ):
    """
    Extract cubes of data for region of interest during campaign period
    """
    # Which variables to use
    if isinstance(vars2use, type(None)):
        # Below are the dust variables for inst3_3d_aer_Nv
        vars2use = ['du{:0>3}'.format(i) for i in range(1,6)]
    # Get the data as dataset via OPeNDAP
    if isinstance(ds, type(None)):
        ds = AC.get_GEOS5_as_ds_via_OPeNDAP(collection=collection, mode=mode,
                                            date=dt)
        # Temporarily do not allow passing of date to OPeNDAP fetcher
#        ds = AC.get_GEOS5_as_ds_via_OPeNDAP(collection=collection, mode=mode, date=None)

    # - Only consider dates that are in dates2use list (if provided)
    # Use the datetime components via pandas as not yet fully in xarray.
    if not isinstance(dates2use, type(None)):
        df = pd.DataFrame(index=ds['time'].values)
        df['date'] = df.index.date
        date_bool = df['date'].isin(dates2use).values
        ds = ds.isel(time=date_bool)

    # - Subset the dataset by region
    if region == 'Cape_Verde':
        # only consider the region the the plane will fly through +2 days transport
        # 20W-55E, 40S-40N,
        # -35E to 10E
        # 0N to 60N
        bool1 = ((ds.lon >= -35) & (ds.lon <= 10)).values
        bool2 = ((ds.lat >= 0) & (ds.lat <= 60)).values
        # Cut by lon, then lat
        ds = ds.isel(lon=bool1)
        ds = ds.isel(lat=bool2)
    elif region == 'global':
        pass
    else:
        print('WARNING: exiting as no region set')
        sys.exit()

    # - Limit the vertical levels
    extr_str = 'ALL_lvls'
    # -  Save the data by variable and day of year reduce packet size of transfers
    # Where to save?
    if isinstance(folder, type(None)):
        folder = './'
    # Save format for data (region, year, doy, variable)
    filestr = 'ARNA_GEOS5_{}_{}_{}_{:0>3}_{}_{}.nc'
    if isinstance(doys2use, type(None)):
        doys2use = list(set(ds['time.dayofyear'].values))
    # Make an extra string to print what slicing the data has had
    pstr = '- getting GEOS5 data as cuboids'
    # Loop by var and save
    for var in vars2use:
        # Now loop and save files in single day units
        print("Getting data via OPeNDAP for var '{}' {}".format(var, pstr))
        for doy in doys2use:
            doy_int = int(float(doy))
            pstr = "Getting data via OPeNDAP for doy '{}' ({})"
            print(pstr.format(doy_int, var))
            # Subset for a single day
            ds_tmp = ds[var].sel(time=ds['time.dayofyear'] == doy)
            # Get the year for the day in the model output
            year = list(set(ds_tmp['time.year'].values))[0]
            # now save this single file
            savename = filestr.format(collection, region, year, doy_int, var,
                                      extr_str)
            print(folder+savename)
            # Remove the file with this name, if it already there
            if rm_existing_file:
                AC.rm_file(folder=folder, filename=savename)
            # Now save the file
            ds_tmp.to_netcdf(folder+savename)
            # Do some memory management
            del ds_tmp
            gc.collect()


def get_GEOSCF_data_cubes4campaign(year=2018, region='Cape_Verde',
                                   doys2use=None,
                                   vars2use=None, folder=None,
                                   limit_lvls=True, limit_lons=False):
    """
    Extract cubes of data for region of interest during campaign period
    """
    # Which variables to use
    if isinstance(vars2use, type(None)):
    #    vars2use = ("pm25du_rh35_gcc", "pm25_rh35_gcc", "noy", "co", "no2", "o3",)
        vars2use = list(convert_GEOSCF_var2GEOSChem_name(rtn_dict=True).keys())
    #    vars2use = ('pm25su_rh35_gcc', 'pm25ni_rh35_gcc', 'pm25soa_rh35_gc', 'pm25ss_rh35_gcc', 'so2')

    # Setup lists of days to use
    dates2use = get_dates4campaigns(year=year)
    print('Using dates:', dates2use)

    # - Get the data as dataset via OPeNDAP
    collection = 'chm_inst_1hr_g1440x721_p23'
    ds = get_GEOSCF_assimlation_ds(collection=collection)

    # - Subset the dataset by date - Only consider dates that are in dates2use list
    # Use the datetime components via pandas as not yet fully in xarray.
    df = pd.DataFrame(index=ds['time'].values)
    df['date'] = df.index.date
    date_bool = df['date'].isin(dates2use).values
    ds = ds.isel(time=date_bool)

    # - Subset the dataset by region
    if region == 'Cape_Verde':
        # only consider the region the the plane will fly through +2 days transport
        # 20W-55E, 40S-40N,
        bool1 = ((ds.lon >= -35) & (ds.lon <= 10)).values
        bool2 = ((ds.lat >= 0) & (ds.lat <= 60)).values
        # Cut by lon, then lat
        ds = ds.isel(lon=bool1)
        ds = ds.isel(lat=bool2)
    elif region == 'global':
        pass
    else:
        print('WARNING: exiting as no region set')
        sys.exit()

    # - Limit the vertical levels
    extr_str = ''
    if limit_lvls:
#        lvls = [1000, 850, 500]
        lvls = [1000, 900, 800, 700, 600, 500]
        bool_lvls = [i in lvls for i in ds['lev']]
        ds = ds.isel(lev=bool_lvls)
        extr_str += 'lvls_'+'_'.join([str(i) for i in lvls])
    else:
        extr_str += 'ALL_lvls'
    # - Limit the longitudes levels
    if limit_lons:
        lons2check = [-18 ,-19.5, -21, -22.5, -24, -25.5]
        idx = [ AC.find_nearest(ds.lon.values, i) for i in  lons2check]
        lons2use = ds.lon.values[ idx ]
        bool_lons = [i in lons2use for i in ds['lon'].values]
        ds = ds.isel(lon=bool_lons)
        extr_str += '_lons_'+'_'.join([str(i) for i in lons2use])
    else:
        pass

    # -  Save the data by variable and day of year reduce packet size of transfers
    # Where to save?
    if isinstance(folder, type(None)):
        folder = './'
    # Save format for data (region, year, doy, variable)
    filestr = 'ARNA_GEOSCF_chm_inst_1hr_g1440x721_p23_{}_{}_{:0>3}_{}_{}.nc'
    if isinstance(doys2use, type(None)):
        doys2use = list(set(ds['time.dayofyear'].values))
    # Make an extra string to print what slicing the data has had
    if limit_lons and limit_lvls:
        pstr = 'sliced by lvls+lons'
    elif limit_lvls:
        pstr = 'sliced by lvls'
    elif limit_lons:
        pstr = 'sliced by lons'
    else:
        pstr =''

    # Loop by var and save
    for var in vars2use:
        # Now loop and save files in single day units
        print("Getting data via OPeNDAP for var '{}' {}".format(var, pstr))
        for doy in doys2use:
            doy_int = int(float(doy))
            pstr = "Getting data via OPeNDAP for doy '{}' ({})"
            print(pstr.format(doy_int, var))
            # Subset for a single day
            ds_tmp = ds[var].sel(time=ds['time.dayofyear'] == doy)
            # now copy locally the data
#            .copy()
            # now save this single file
            savename = filestr.format(region, year, doy_int, var, extr_str)
            print(folder+savename)
            ds_tmp.to_netcdf(folder+savename)
            # Do some memory management
            del ds_tmp
            gc.collect()


def get_GEOS5_data_cubes4campaign(ds=None, year=2018, region='Cape_Verde',
                                  doys2use=None,
                                  vars2use=None, collection='inst3_3d_aer_Nv',
                                  limit_lvls=False, limit_lons=False):
    """
    Extract cubes of data for region of interest during campaign period
    """
    # Which variables to use
    if isinstance(vars2use, type(None)):
        vars2use = ['du{:0>3}'.format(i) for i in range(1,6)]
    # Setup lists of days to use
    dates2use = get_dates4campaigns(year=year)
    print('Using dates:', dates2use)
    # - Get the data as dataset via OPeNDAP
    if isinstance(ds, type(None)):
        ds = AC.get_GEOS5_as_ds_via_OPeNDAP(collection=collection)
    # - Subset the dataset by date - Only consider dates that are in dates2use list
    # Use the datetime components via pandas as not yet fully in xarray.
    df = pd.DataFrame(index=ds['time'].values)
    df['date'] = df.index.date
    date_bool = df['date'].isin(dates2use).values
    ds = ds.isel(time=date_bool)

    # - Subset the dataset by region
    if region == 'Cape_Verde':
        # only consider the region the the plane will fly through +2 days transport
        # 20W-55E, 40S-40N,
        bool1 = ((ds.lon >= -35) & (ds.lon <= 10)).values
        bool2 = ((ds.lat >= 0) & (ds.lat <= 60)).values
        # Cut by lon, then lat
        ds = ds.isel(lon=bool1)
        ds = ds.isel(lat=bool2)
    elif region == 'global':
        pass
    else:
        print('WARNING: exiting as no region set')
        sys.exit()
    # State that all data (cuboids are being extracted)
    extr_str = 'ALL_lvls'

    # -  Save the data by variable and day of year reduce packet size of transfers
    # Where to save?
    if isinstance(folder, type(None)):
        folder = './'
    # Save format for data (region, year, doy, variable)
    filestr = 'ARNA_GEOS5_{}_{}_{}_{:0>3}_{}_{}.nc'
    if isinstance(doys2use, type(None)):
        doys2use = list(set(ds['time.dayofyear'].values))
    # Make an extra string to print what slicing the data has had
    if limit_lons and limit_lvls:
        pstr = 'sliced by lvls+lons'
    elif limit_lvls:
        pstr = 'sliced by lvls'
    elif limit_lons:
        pstr = 'sliced by lons'
    else:
        pstr =''

    # Loop by var and save
    for var in vars2use:
        # Now loop and save files in single day units
        print("Getting data via OPeNDAP for var '{}' {}".format(var, pstr))
        for doy in doys2use:
            doy_int = int(float(doy))
            pstr = "Getting data via OPeNDAP for doy '{}' ({})"
            print(pstr.format(doy_int, var))
            # Subset for a single day
            ds_tmp = ds[var].sel(time=ds['time.dayofyear'] == doy)
            # Get the year for the day in the model output
            year = list(set(ds_tmp['time.year'].values))[0]
            # now save this single file
            savename = filestr.format(collection, region, year, doy_int, var,
                                      extr_str)
            print(folder+savename)
            ds_tmp.to_netcdf(folder+savename)
            # Do some memory management
            del ds_tmp
            gc.collect()


def get_data_surface4campaign(year=2018, region='Cape_Verde', doys=None):
    """
    Extract surface data for region of interest during campaign period
    """
    # - Local variables
    # Get useful variables from the  3D cubes
#    vars2use = "pm25du_rh35_gcc", "pm25_rh35_gcc", "noy", "co",  "no2", "o3"
    vars2use = list(convert_GEOSCF_var2GEOSChem_name(rtn_dict=True).keys())
    # Add extra surface species
#    vars2use += [ 'hno3', 'nit', 'no', 'n2o5' 'dst4', 'dst3', 'dst2', 'dst1',]
    extra_vars = [
        'xyle', 'dst2', 'hno4', 'pm25su_rh35_gcc', 'pm25ni_rh35_gcc', 'ocpi', 'eoh',
        'benz', 'so2', 'rcho', 'h2o2', 'ald2', 'pm25soa_rh35_gc', 'n2o5', 'tolu',
        'pm25_rh35_gocar', 'dst4', 'no', 'ocpo', 'alk4', 'nh3', 'nh4', 'prpe', 'c2h6',
        'pm25bc_rh35_gcc', 'salc', 'hno3', 'ch4', 'c3h8', 'dst3', 'sala', 'bcpo',
        'dst1', 'macr', 'hcho', 'soap', 'acet', 'isop', 'soas', 'pm25oc_rh35_gcc',
        'mvk', 'nit', 'bcpi', 'pan', 'mek', 'pm25ss_rh35_gcc'
    ]
    vars2use += extra_vars

    # Also get other detail from the 2D fields
    vars2use = vars2use[::-1]
    dates2use = get_dates4campaigns(year=year)
    # Add to filename that the data is just the model surface
    extr_str = 'surface'
    # - Get the data as dataset via OPeNDAP
    collection = 'chm_tavg_1hr_g1440x721_v1'
    ds = get_GEOSCF_assimlation_ds(collection=collection)
    # - Subset the dataset by date - Only consider dates that are in dates2use list
    # Use the datetime components via pandas as not yet full yin xarray.
    df = pd.DataFrame(index=ds['time'].values)
    df['date'] = df.index.date
    date_bool = df['date'].isin(dates2use).values
    ds = ds.isel(time=date_bool)

    # - Subset the dataset by region
    if region == 'Cape_Verde':
        # only consider the region the the plane will fly through +2 days transport
        # 20W-55E, 40S-40N,
        bool1 = ((ds.lon >= -35) & (ds.lon <= 10)).values
        bool2 = ((ds.lat >= 0) & (ds.lat <= 60)).values
        # Cut by lon, then lat
        ds = ds.isel(lon=bool1)
        ds = ds.isel(lat=bool2)
    elif region == 'global':
        pass
    else:
        print('WARNING: exiting as no region set')
        sys.exit()

    # -  Save the date in the earth0
    # Where to save the extracted data?
    if isinstance(folder, type(None)):
        folder = './'
    # Save format for data (region, year, doy, variable)
    filestr = 'ARNA_GEOSCF_{}_{}_{}_{:0>3}_{}_{}.nc'
    # Loop by var and save
#    for var in vars2use:
    # Now loop and save files in single day units
#    print("Getting data via OPeNDAP for var '{}'".format(var))
    if isinstance(doys, type(None)):
        doys = list(set(ds['time.dayofyear'].values))
    for doy in doys:
        doy_int = int(float(doy))
        print("Getting data via OPeNDAP for doy '{}'".format(doy_int))
        # Subset for a single day
        ds_tmp = ds[vars2use].sel(time=ds['time.dayofyear'] == doy)
        # Now copy locally the data
        ds_tmp = ds_tmp.copy()
        # Now save this single file
        savename = filestr.format(
            collection, region, year, doy_int, 'ALL', extr_str)
        ds_tmp.to_netcdf(folder+savename)
        # Do some memory management
        del ds_tmp
        gc.collect()


def get_dates4campaigns(year=2018, all_dates=False,
                        ground_CVAO_campaign=False):
    """
    Get a list of dates for campaigns
    """
    # Setup lists of days to use
    days2use = np.arange(1, 15)
    if isinstance(year, int):
        dates2use = [datetime.date(year, 2, i) for i in days2use]
    elif all_dates:
        dates2use = [datetime.date(2018, 2, i) for i in days2use]
        dates2use += [datetime.date(2019, 2, i) for i in days2use]
    else:
        print('WARNING: A integer year must be provided')
        sys.exit()
    # The 2019 August CVAO ground campaign (+2 days flights)
    if ground_CVAO_campaign:
        days2use = np.arange(1, 28)
        dates2use = [datetime.date(2019, 8, i) for i in days2use]
    return dates2use


def get_GEOSCF_assimlation_ds(collection='chm_inst_1hr_g1440x721_p23',
                              mode='assim'):
    """
    Wrapper to get the GEOS-CF assimilation data
    """
    # Get entire archive as a xr.Dataset
    ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(mode=mode, collection=collection)
#    ds = get_GEOSCF_as_ds_via_OPeNDAP_LOCAL(mode=mode, collection=collection)
    return ds


def get_GEOSCF_vertical_levels(print_equivalents=False, native_levels=False):
    """
    get a dictionary of GEOS Composition Forecast (GEOS-CF) vertical levels
    """
    # Get a list of the pressure levels in GEOS-CF
    if native_levels:
        HPa_l  = [
        0.01, 0.02, 0.0327, 0.0476, 0.066, 0.0893, 0.1197, 0.1595, 0.2113, 0.2785, 0.365, 0.4758, 0.6168, 0.7951, 1.0194, 1.3005, 1.6508, 2.085, 2.6202, 3.2764, 4.0766, 5.0468, 6.2168, 7.6198, 9.2929, 11.2769, 13.6434, 16.4571, 19.7916, 23.7304, 28.3678, 33.81, 40.1754, 47.6439, 56.3879, 66.6034, 78.5123, 92.3657, 108.663, 127.837, 150.393, 176.93, 208.152, 244.875, 288.083, 337.5, 375.0, 412.5, 450.0, 487.5, 525.0, 562.5, 600.0, 637.5, 675.0, 700.0, 725.0, 750.0, 775.0, 800.0, 820.0, 835.0, 850.0, 865.0, 880.0, 895.0, 910.0, 925.0, 940.0, 955.0, 970.0, 985.0
        ]
    else:
        HPa_l = [
            1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400,
            350, 300, 250, 200, 150, 100, 50
        ]

    # Get altitudes in km, then convert to metres
    Alt = AC.hPa_to_Km(HPa_l)
    Alt = np.array(Alt) / 1E3
    # Print out a summary
    if print_equivalents:
        # For just HPa and km
        pstr = 'A pressure {:>4}HPa of is equiv to {:>4,.3f}km'
        Alt_dict = dict(zip(HPa_l, Alt*1E3))
        for HPa in HPa_l:
            print(pstr.format(HPa, Alt_dict[HPa]))
        # Also for kft
        pstr = 'A press. {:>4} HPa of is equiv to {:>4,.3f} km ({:>4,.3f} kft)'
        for HPa in HPa_l:
            print(pstr.format(HPa, Alt_dict[HPa], Alt_dict[HPa]/304.8))
    return HPa_l


def get_local_folder(key, host=None, rtn_dict=False):
    """
    Hold folders in a dictionary and return specific variables or a dictionary
    """
    import platform
    # Get the host
    host = platform.node()
    # - Set locations of York HPC
    if ('viking' in host):
        NASA_data = '/mnt/lustre/groups/chem-acm-2018/earth0_data/NASA/'
    elif ('earth0' in host):
        NASA_data = '/work/data/NASA/'
    else:
        print( 'NASA folder loction not known' )
    # - Setup a dictionary for the variables
    d = {
    'NASA_data': NASA_data,
#    'folder4plots': folder4plots,
    }
    return d[key]


if __name__ == "__main__":
    main()

"""
Functions for processing/analysing ARNA campaign observations
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

# Import from elsewhere in ARNA module
from . core import *
from . utils import *


def get_coordinates_from_NetCDF_file(ds=None, folder=None, filename=None,
                                    falt_var='PS_RVSM',
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
    # NOTE: The below conversion is not advised.
    #       Use the external pressure variable instead (PS_RVSM).
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


def set_flagged_data2NaNs(ds, VarName='no_mr', flag2use=0,
                          FlagName='no_flag'):
    """
    Set the flagged data in a dataset to be NaNs
    """
    # What do the flags mean? (copied from FAAM NetCDF metadata)
    # ( file: 'core-nitrates_faam_20200211_v002_r1_c224.nc')
    # Flag=0 indicates good data.
    # Flag=1 indicates reduced quality data.
    # Flag=2 indicates suspect data.
    # Flag=3 indicates missing or invalid data.
    # Create a boolean
    bool = ds[FlagName].values != flag2use
    # Make values without the flagged value into NaNs
    ds[VarName].loc[ bool ] = np.NaN
    return ds


def get_CIMS_data4flight(flight_ID='C216', resample_data=True):
    """
    Retrieve ToF-CIMS data from ARNA flights
    """
    # Where is the data?
    folder = '{}/{}/'.format(get_local_folder('ARNA_data'), 'CIMS')
    # What is the name of the sheet in the excel file?
    sheet_name = flight_ID
    # Retrive the core CIMS observations (use the time coordinate for index)
    filename = 'ACSIS6_0.1hz_Bromine.xlsx'
    df = pd.read_excel( folder + filename, sheet_name=sheet_name)
    dt_var = 'Date:Time'
    df.index = pd.to_datetime( df[dt_var].values )
    del df[dt_var]
    # Also get HNO3 / HONO (use the time coordinate for index)
    filename = 'ACSIS6_ARNA_1hz_HNO3_HONO_ToSend.xlsx'
    df2 = pd.read_excel( folder + filename, sheet_name=sheet_name)
    dt_var = 'date'
    df2.index = pd.to_datetime( df2[dt_var].values )
    del df2[dt_var]
    # Merge the files
    df = pd.concat([df,df2], axis="index")
	# Resample the data?
    if resample_data:
        df = df.resample('1T' ).mean()
    return df


def get_FAAM_core4flightnum(flight_ID='C216', version='v2020_06',
                            resample_data=True):
    """
    Get the core FAAM flight data for a specific flight
    """
	# Where are the files? Which files to use?
    folder = '{}/CEDA/{}/'.format( get_local_folder('ARNA_data'), version)
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
    # Set flagged data to NaNs
    ds = set_flagged_data2NaNs(ds, VarName='O3_TECO', FlagName='O3_TECO_FLAG')
    ds = set_flagged_data2NaNs(ds, VarName='CO_AERO', FlagName='CO_AERO_FLAG')
    # Check for the core NOx file
    Nitrates_filename = [ i for i in files2use if 'core-nitrates' in i ]
    if len(Nitrates_filename) >=1:
        ass_str = 'More than one nitrates file found for flight!'
        assert len(Nitrates_filename) == 1, ass_str
        ds2 = xr.open_dataset(Nitrates_filename[0])
        ds2 = ds2.mean(dim='sps10')
        # Remove flagged values for NO and NO2
        ds2 = set_flagged_data2NaNs(ds2, VarName='no_mr', FlagName='no_flag')
        ds2 = set_flagged_data2NaNs(ds2, VarName='no2_mr', FlagName='no2_flag')
        # Merge the files
        ds = xr.merge([ds, ds2])
		# Get the HONO data too (if it is present)
        try:
            ds3 = xr.open_dataset(Nitrates_filename[0],
                                 group='non_core_group')
            ds3 = ds3.mean(dim='sps10')
            # Remove flagged values for NO and NO2
            ds3 = set_flagged_data2NaNs(ds3, VarName='hono_mr',
                                        FlagName='hono_flag')
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
    if resample_data:
        df = df.resample('1T' ).mean()
    return df


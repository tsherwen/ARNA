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
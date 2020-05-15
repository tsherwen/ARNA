"""
Miscellaneous functions for ARNA campaign work (pre/during/after)
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

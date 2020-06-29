"""
Functions for processing/analysing ARNA campaign observations
"""

import os
import sys
import gc
import glob
import xarray as xr
import numpy as np
import AC_tools as AC
import pandas as pd
from netCDF4 import Dataset
from datetime import datetime as datetime_
import datetime as datetime
import time
from time import gmtime, strftime

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


def get_SWAS_data4flight(flight_ID=None):
    """
    Retrieve SWAS data for ARNA flights
    """
    # Where is the data?
    folder = '{}/{}/'.format(get_local_folder('ARNA_data'), 'SWAS')
    filename = 'ARNA-FIRSTLOOK-SWAS_JDL_typo_fix.csv'
    var2use = 'comments4'
    format = '%d/%m/%Y %H:%M:%S'
    # Or use latest file (NOTE: issue with column formating)
#    filename = 'ARNA-SECONDLOOK-SWAS.csv'
    df = pd.read_csv(folder+filename)
    # Get a
    print('TODO: Speak to Jimmy to confirm the datetime columns!')
#    var2use = 'SAMPLE START TIME'
#    format = '%d/%m/%Y %H:%M:%S'
    df.index =  pd.to_datetime(df[var2use].values, format=format)
    # If a flight ID stated, then only return points for that flight
    if isinstance(flight_ID, type(None)):
        pass
    else:
        # Get the beginning and end of the flight
        dfS = get_summary4flight(flight_ID=flight_ID)
        sdate = dfS.index.min()
        edate = dfS.index.max()
        # Only consider samples within this time
        df = df.loc[df.index >= sdate, :]
        df = df.loc[df.index <= edate, :]
    return df


def map_SWAS_var2GEOS_var(var, invert=False):
    """
    Map variables names from SWAS to GEOS variable names
    """
    d = {
#    '1_3_butadiene':,
#    '1_butene':,
#    '2_3_methylpentane':,
#    '224_tmp':,
    'acetaldehyde': 'ALD2',
    'acetone': 'ACET',
#    'acetylene':,
    'benzene': 'BENZ', # GEOSChem, but not GEOS-CF output
#    'benzenechb':,
    'cis_2_butene': 'PRPE', #  NOTE: lumped tracer for >= C3 alkenes
    'cyclo_pentane': 'ALK4', #  NOTE: lumped tracer for >= C4 Alkanes
    'dms': 'DMS', # GEOSChem, but not GEOS-CF output
#    'dmschb':,
    'ethane': 'C2H6',
#    'ethene':,
#    'ethylbenzene':,
#    'extra_1':,
#    'extra_2':,
#    'extra_3':,
#    'extra_4':,
#    'extra_5':,
#    'extra_l2':,
    'iso_butane': 'ALK4', #  NOTE: lumped tracer for >= C4 Alkanes
    'iso_butene': 'PRPE', #  NOTE: lumped tracer for >= C3 alkenes
    'iso_pentane': 'ALK4', #  NOTE: lumped tracer for >= C4 Alkanes
    'isoprene': 'PRPE', #  NOTE: lumped tracer for >= C3 alkenes
    'methanol': 'MOH',
    'mp_xylene': 'XYLE',
    'n_butane': 'ALK4', #  NOTE: lumped tracer for >= C4 Alkanes
    'n_heptane': 'ALK4', #  NOTE: lumped tracer for >= C4 Alkanes
    'n_hexane': 'ALK4', #  NOTE: lumped tracer for >= C4 Alkanes
    'n_octane': 'ALK4', #  NOTE: lumped tracer for >= C4 Alkanes
    'n_pentane': 'ALK4', #  NOTE: lumped tracer for >= C4 Alkanes
    'o_xylene': 'XYLE',
    'pent_1_ene':  'PRPE', #  NOTE: lumped tracer for >= C3 alkenes
    'propane': 'C3H8',
    'propene': 'PRPE', #  NOTE: lumped tracer for >= C3 alkenes
    'toluene': 'TOLU',
#    'toluenechb':,
    'trans_2_butene': 'PRPE', #  NOTE: lumped tracer for >= C3 alkenes
    'trans_2_pentene':'PRPE', #  NOTE: lumped tracer for >= C3 alkenes
    }
    # Invert the dictionary?
    if invert:
        d = {v: k for k, v in list(d.items())}
    return d[var]


def get_ARNA_flight_log_as_df():
    """
    Make a single pd.DataFrame with all flight summaries
    """
    flights_nums = [
#    216,
    217, 218, 219, 220, 221, 222, 223, 224, 225
    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    dfs = []
    for flight_ID in flight_IDs:
        dfs += [get_summary4flight( flight_ID=flight_ID )]
    # Combine and return as a single dataframe sorted by time
    df = pd.concat(dfs)
    df = df.sort_index()
    return df


def get_summary4flight(flight_ID='C217'):
    """
    retrieve a FAAM flight summary as a dataframe
    """
    folder = '{}/{}/'.format(get_local_folder('ARNA_data'), 'CEDA/v2020_06')
    filename = 'flight-sum_faam_*_*_{}.csv'.format(flight_ID.lower())
    file2use = glob.glob( folder+filename )
    ass_str = 'WARNING: {} flight summaries found present for flight {}!'
    assert len(file2use) == 1, ass_str.format(file2use, flight_ID)
    # Add Gotcha for missing header in FAAM file archived at CEDA
    columns = [
    'Event', 'Start', 'Start Hdg / °', 'Start Hgt / kft', 'Start Lat / °',
    'Start Long / °', 'Stop', 'Stop Hdg / °', 'Stop Hgt / kft',
    'Stop Lat / °', ' Stop Long / °', 'Comment',
    ]
    if flight_ID=='C217':
        header = None
        names = columns
    else:
        header = 0
        names = None
    # read file
    df = pd.read_csv(file2use[0], header=header, names=names)
    # Add a flight ID column
    df['flight_ID'] = flight_ID
    # Update start column to be in datatime format
    var2use = 'Start'
    format = '%Y-%m-%d %H:%M:%S'
    df.index = pd.to_datetime( df[var2use].values, format=format)
    return df


def get_CIMS_data4flight(flight_ID='C216', resample_data=True, debug=False):
    """
    Retrieve ToF-CIMS data from ARNA flights
    """
    # Where is the data?
    folder = '{}/{}/'.format(get_local_folder('ARNA_data'), 'CIMS')
    # What is the name of the sheet in the excel file?
    sheet_name = flight_ID
    # Retrive the core CIMS observations (use the time coordinate for index)
    filename = 'ACSIS6_0.1hz_Bromine.xlsx'
    format='%m/%d/%Y %H:%M:%S'
    df = pd.read_excel(folder + filename, sheet_name=sheet_name,
                       date_parser=format)
    xl = pd.ExcelFile( folder + filename)
    if debug:
        print(xl.sheet_names)
    dt_var = 'Date:Time'
    dates = AC.dt64_2_dt(df[dt_var].values )
    dates = [i.strftime(format) for i in dates ]
    dates = [datetime.datetime.strptime(i,'%d/%m/%Y %H:%M:%S') for i in dates]
    df.index = dates
    del df[dt_var]
    # Also get HNO3 / HONO (use the time coordinate for index)
    filename = 'ACSIS6_ARNA_1hz_HNO3_HONO_ToSend.xlsx'
    df2 = pd.read_excel( folder + filename, sheet_name=sheet_name)
    xl = pd.ExcelFile( folder + filename)
    if debug:
        print(xl.sheet_names)
    dt_var = 'date'
#    df2.index = pd.to_datetime( df2[dt_var].values, format='%m/%d/%Y %H:%M:%S')
    dates = AC.dt64_2_dt(df2[dt_var].values )
    dates = [i.strftime(format) for i in dates ]
    dates = [datetime.datetime.strptime(i,'%d/%m/%Y %H:%M:%S') for i in dates]
    df2.index = dates
    del df2[dt_var]
    # Merge the files
    df = pd.concat([df,df2], axis="index")
    # Update the variable names
    VarNameDict = {
    'BrO (ppt*)':'BrO', 'Br2 + HOBr (ppt*)':'Br2+HOBr', 'HONO (ppt*)':'HONO',
    'HNO3 (ppt*) ':'HNO3',  'HNO3 (ppt*)' : 'HNO3',
    }
    df = df.rename(columns=VarNameDict)
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


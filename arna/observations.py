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
#    filename = 'ARNA-FIRSTLOOK-SWAS_JDL_typo_fix.csv'
#    var2use = 'comments4'
#    format = '%d/%m/%Y %H:%M:%S'
    # Or use latest file (NOTE: issue with column formating)
#    filename = 'ARNA-SECONDLOOK-SWAS.csv'
    # Use the updated second look file
    filename = 'ARNA-SECONDLOOK-SWAS_v2.csv'
    df = pd.read_csv(folder+filename)
    print(filename)
    # Update the index to use the SWAS fire fime
    var2use = 'SAMPLE START TIME'
    format = '%d/%m/%Y %H:%M:%S'
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


def get_filters_data4flight(flight_ID='C217', all_flights=True):
    """
    Retrieve filters data from ARNA flights
    """
    # Where is the data?
    folder = '{}/{}/'.format(get_local_folder('ARNA_data'), 'Filters')
    # What is the name of the sheet in the excel file?
#    filename = 'Overview_all_filters_ACSIS_5_and_ARNA-1.xlsx'
#    filename = 'Overview_filters_ARNA_2.xlsx'
    filename = 'Overview_filters_ARNA_2_TMS_edits.xlsx'
    sheet_name = 'Sheet1'
    dfFULL = pd.read_excel(folder + filename, sheet_name=sheet_name)
    # Now Just look at core data of interest
    CoreCols = [
    'Day', 'Flight', 'Filter', 'height', 'Time on', 'Time off',
    'Airflow (stL)',
    ]
    # - Select nitrate data
    # Yes, GEOS-CF has sulfate variables output - 'NIT', 'NITs'
    NO3cols = [i for i in dfFULL.columns if 'NO3' in i]
    dfN = dfFULL[CoreCols + NO3cols]
    #
#    NO3_var = 'NO3.total'
    NO3_var2use = ['NO3.2', 'NO3.5']
    units = 'nanomoles/m3'
#    NO3_idx = [list(dfFULL.columns).index(i) for i in var2inc]

    # - Also save sulfate?
    # Yes, GEOS-CF has sulfate variables output - 'SO4', 'SO4s'
#    SO4cols = [i for i in dfFULL.columns if 'SO4' in i]
#    dfS = dfFULL[CoreCols + SO4cols]
    SO4_var2use = ['SO4.2', 'SO4.5']
    units = 'nanomoles/m3'
#    SO4_idx = [list(dfFULL.columns).index(i) for i in SO4_var2use]

    # - Now chop off excess headers and make sure formats are correct
    df = dfFULL.loc[dfFULL.index[2:],:]
    #
#    idx2use = [list(dfFULL.columns).index(i) for i in CoreCols]
#    idx2use += NO3_idx + SO4_idx
#    cols2use = [list(dfFULL.columns)[i] for i in idx2use ]
    df = df[ CoreCols + NO3_var2use + SO4_var2use ]
    # Replace values less than black/NaNs with np.NaN
    df = df.replace('lower than blank', np.NaN)
    df = df.replace('NaN', np.NaN)
    # Remove blanks (as these are NaNs)
    df = df.rename_axis(None)
    df = df.loc[ (df['height'] != 'blank').values,:]
    # Update sampling times to date times
    # Start time
    TimeOnVar = 'Time on'
    sdate_var = 'Sample Start'
    df[sdate_var] = df['Day'].astype(str) + ' ' +df[TimeOnVar].astype(str)
    format = '%Y-%m-%d %H:%M:%S'
    df[sdate_var] = pd.to_datetime( df[sdate_var].values, format=format)
    del df[TimeOnVar]
    # End time
    TimeOffVar = 'Time off'
    edate_var = 'Sample End'
    df[edate_var] = df['Day'].astype(str) + ' ' +df[TimeOffVar].astype(str)
    format = '%Y-%m-%d %H:%M:%S'
    df[edate_var] = pd.to_datetime( df[edate_var].values, format=format)
    del df[TimeOffVar]
    # calculate mid point of sampling and set this as index
    mdate_var = 'Sample Mid'
#     df[edate_var] - df[sdate_var]
#     def calc_mid_point_as_dt(sdate, edate):
#         print(sdate, edate)
# #        sdate = sdate.astype(numpy.datetime64)
# #        edate = edate.astype(numpy.datetime64)
#         sdate = sdate.datetime()
#         edate = edate.datetime()
#         print([type(i) for i in (sdate, edate)])
#         sdate, edate = AC.dt64_2_dt([sdate, edate])
#         print(sdate, edate)
#         interval = (edate - sdate).total_seconds()
#         assert interval > 0, 'WARNING: filter internval <0!'
#         half_interval = interval/2
#         return AC.add_secs(sdate, half_interval )
#     df[mdate_var] = df.apply(lambda x: calc_mid_point_as_dt(
#                              sdate=x[sdate_var],
#                              edate=x[edate_var],
#                              ), axis=1)
    # TODO - check dates (as some sample times are implausible)
    df[mdate_var] = df[edate_var] - df[sdate_var]

    # Just use the start time for now
    df.index = df[sdate_var]
    df = df.rename_axis(None)
    # - Just consider totals for species of interest
    NO3_var = 'NO3.total'
    df[NO3_var] = df[ NO3_var2use ].sum(axis=1)
    SO4_var = 'SO4.total'
    df[SO4_var] = df[ SO4_var2use ].sum(axis=1)
    del dfFULL
    # Convert to ug/m3 from 'nanomoles/m3'
    df[NO3_var] = df[NO3_var].values / 1E9 * AC.species_mass('NIT') * 1E6
    df[SO4_var] = df[SO4_var].values / 1E9 * AC.species_mass('SO4') * 1E6
    # Return all flights unless a specific flight requested
    if all_flights:
        return df
    else:
        return df.loc[ df['Flight']==flight_ID, :]


def get_CIMS_data4flight(flight_ID='C216', resample_data=True, debug=False):
    """
    Retrieve ToF-CIMS data from ARNA flights
    """
    # Where is the data?
    folder = '{}/{}/'.format(get_local_folder('ARNA_data'), 'CIMS')
    # What is the name of the sheet in the excel file?
    sheet_name = flight_ID
    # - Setup a shared (1Hz) time axis to merge onto
    # Use time index for timings of flight, then add a 2 hour buffer
    dfS = get_summary4flight(flight_ID=flight_ID)
    sdate = dfS.index.values.min()
    edate = dfS.index.values.max()
    sdate, edate = AC.dt64_2_dt([sdate, edate])
    sdate = AC.add_hrs(sdate, -2)
    edate = AC.add_hrs(sdate, 2)
    index = pd.date_range(start=sdate, end=edate, freq='1S')
    dfM = pd.DataFrame(index=index)
    # - Retrive the core CIMS observations (use the time coordinate for index)
    try:
        filename = 'ACSIS6_0.1hz_Bromine.xlsx'
        format = '%m/%d/%Y %H:%M:%S'
        df = pd.read_excel(folder + filename, sheet_name=sheet_name,
                           date_parser=format)
        xl = pd.ExcelFile( folder + filename)
        if debug:
            print(xl.sheet_names)
        dt_var = 'Date:Time'
        dates = AC.dt64_2_dt(df[dt_var].values )
        dates = [i.strftime(format) for i in dates ]
        dates = [datetime_.strptime(i,'%d/%m/%Y %H:%M:%S') for i in dates]
        df.index = dates
        del df[dt_var]
        # Merge the files
        dfM = pd.concat([dfM, df], axis="index")
    except:
        pstr = "WARNING: failed to include CIMS halogen data for '{}' in df"
        print(pstr.format(flight_ID))
    # - Also get HNO3 / HONO (use the time coordinate for index)
    try:
        filename = 'ACSIS6_ARNA_1hz_HNO3_HONO_ToSend.xlsx'
        df2 = pd.read_excel( folder + filename, sheet_name=sheet_name)
        xl = pd.ExcelFile( folder + filename)
        if debug:
            print(xl.sheet_names)
        dt_var = 'date'
    #    df2.index = pd.to_datetime( df2[dt_var].values, format='%m/%d/%Y %H:%M:%S')
        dates = AC.dt64_2_dt(df2[dt_var].values )
        dates = [i.strftime(format) for i in dates ]
        dates = [datetime_.strptime(i,'%d/%m/%Y %H:%M:%S') for i in dates]
        df2.index = dates
        del df2[dt_var]
        # Merge the files
        dfM = pd.concat([dfM, df2], axis="index")
    except:
        pstr = "WARNING: failed to include CIMS NOy data for '{}' in df"
        print(pstr.format(flight_ID))

    # - Also include HCN data
    try:
        filename = 'ACSIS6_ARNA2_HCN_James_TMS_Update.xlsx'
        df = pd.read_excel( folder + filename, sheet_name=sheet_name)
        xl = pd.ExcelFile( folder + filename)
        if debug:
            print(xl.sheet_names)
        dt_var = 'date_time'
    #    df2.index = pd.to_datetime( df2[dt_var].values, format='%m/%d/%Y %H:%M:%S')
        dates = AC.dt64_2_dt(df[dt_var].values )
        dates = [i.strftime(format) for i in dates ]
        dates = [datetime_.strptime(i,'%d/%m/%Y %H:%M:%S') for i in dates]
        df.index = dates
        del df[dt_var]
        # Merge the files
        dfM = pd.concat([dfM, df], axis="index")
    except:
        pstr = "WARNING: failed to include CIMS HCN data for '{}' in df"
        print(pstr.format(flight_ID))

    # Update the variable names
    VarNameDict = {
    'BrO (ppt*)':'BrO', 'Br2 + HOBr (ppt*)':'Br2+HOBr', 'HONO (ppt*)':'HONO',
    'HNO3 (ppt*) ':'HNO3',  'HNO3 (ppt*)' : 'HNO3',
    'HCN (ppt*)': 'HCN',
    }
    dfM = dfM.rename(columns=VarNameDict)
    # Include a flag for flight ID
    dfM['flight_ID'] = flight_ID
	# Resample the data?
    if resample_data:
        dfM = dfM.resample('1T' ).mean()
    return dfM


def add_derived_variables2FAAM_data(df):
    """
    Add variables derived from other variables
    """
    # A metric for dust
    try:
        VarName = 'IS_DUST'
        df[VarName] = False
        # PCASP total above
        PCASP_threshold = 200
        bool = df['PCAS2CON'] >= PCASP_threshold
        df.loc[bool, VarName] = True
        # Above the boundary layer - just set as 900 hPa for now
        hPa_threshold = 900
        bool = df['PS_RVSM'] >= hPa_threshold
        df.loc[bool ,VarName] = False
    except KeyError:
        print("Derived variable not added to dataframe ({})".format(VarName))
    # A indicator for straight and level runs (SLR)
    try:
        VarName = 'IS_SLR'
        df[VarName] = False
        # PCASP total above
        var2use = 'VELD_GIN'
        threshold = 1
        bool1 = (df[var2use] <= threshold) & (df[var2use] >= -threshold)
        var2use = 'ROLL_GIN'
        threshold = 2.5
        bool2 = (df[var2use] <= threshold) & (df[var2use] >= -threshold)
        df.loc[bool1 & bool2, VarName] = True
    except KeyError:
        print("Derived variable not added to dataframe ({})".format(VarName))
    return df


def add_biomass_flag2df(df, CIMSdf=None, flight_ID='C216', threshold=None):
    """
    Add a flag for biomass burning to dataframe using CIMs HCN data
    """
    # Retrieve the CIMS HCN data if not provided
    if isinstance(CIMSdf, type(None)):
        df = get_CIMS_data4flight(flight_ID=flight_ID)
    # What Criterion to use for biomass burning?
    # 3 sigma enhancement over the background, and 10 sigma for a direct plume.
    var2use = 'HCN'
    NewVar = 'IS_BB'
    if isinstance(threshold, type(None)):
        background = 10
        stats = df[var2use].describe()
        sigma = stats.loc[ stats.index=='std'].values[0]
        threshold = background + (sigma*3)
    # Add a boolean for biomass burning (BB)
    df[NewVar] = False
    df.loc[df[var2use]>=threshold,  NewVar] = True
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
    # do the same for vertical velocity and roll
    # 'Roll angle from POS AV 510 GPS-aided Inertial Nav. unit (positive for left wing up)'
    VarName = 'ROLL_GIN'
    FlagName = 'ROLL_GIN_FLAG'
    ds = set_flagged_data2NaNs(ds, VarName=VarName, FlagName=FlagName)
    # 'Aircraft velocity down from POS AV 510 GPS-aided Inertial Navigation     # 'Aircraft velocity down from POS AV 510 GPS-aided Inertial Navigation unit'
    VarName = 'VELD_GIN'
    FlagName = 'VELD_GIN_FLAG'
    ds = set_flagged_data2NaNs(ds, VarName=VarName, FlagName=FlagName)
    # Convert to a dataframe
    df = ds.to_dataframe()
    df.index.rename(None)
    del ds
    # Check for the core NOx file
    Nitrates_filename = [ i for i in files2use if 'core-nitrates' in i ]
    if len(Nitrates_filename) >=1:
        ass_str = 'More than one nitrates file found for flight!'
        assert len(Nitrates_filename) == 1, ass_str
        ds = xr.open_dataset(Nitrates_filename[0])
        ds = ds.mean(dim='sps10')
        # Remove flagged values for NO and NO2
        ds = set_flagged_data2NaNs(ds, VarName='no_mr', FlagName='no_flag')
        ds = set_flagged_data2NaNs(ds, VarName='no2_mr', FlagName='no2_flag')
        # Merge the files
        df2 = ds.to_dataframe()
        vars2use = [i for i in df2.columns if 'no' in i]
        df = pd.concat([df, df2[vars2use]], axis=1)
        del ds, df2
		# Get the HONO data too (if it is present)
        try:
            ds = xr.open_dataset(Nitrates_filename[0], group='non_core_group')
            ds = ds.mean(dim='sps10')
            # Remove flagged values for NO and NO2
            ds = set_flagged_data2NaNs(ds, VarName='hono_mr',
                                       FlagName='hono_flag')
            # Merge the files
            df2 = ds.to_dataframe()
            # Use the same main time index as rest of file
            df2.index = df.index.values
            df2.index.rename(None)
            df = pd.concat([df, df2], axis=1)
            del ds, df2
        except OSError:
            print('WARNING: no HONO data found for {}'.format(flight_ID))
    # Check for cloud physics file
    cloud_phy_filenames = [ i for i in files2use if 'cloud-phy' in i ]
    # Just include the main file that ends "<flight_ID>.nc"
    suffix = '{}.nc'.format( flight_ID.lower() )
    cloud_phy_filename = [i for i in cloud_phy_filenames if i.endswith(suffix)]
    if len(cloud_phy_filename) >=1:
        ass_str = 'More than one main cloud phys file found for flight!'
        assert len(cloud_phy_filename) == 1, ass_str
        try:
            ds = xr.open_dataset(cloud_phy_filename[0])
            # Remove flagged values for PCASP
            VarName = 'PCAS2CON'
            FlagName = 'PCAS2_FLAG'
            ds = set_flagged_data2NaNs(ds, VarName=VarName, FlagName=FlagName)
            # Just include PCASP concentration and a flag for this.
            df2 = ds[ [VarName, FlagName] ].to_dataframe()
            df2.index = df2.index.floor('S')
            df2.index = df2.index.values # rm naming of time coord 'PCAS2TSPM'
            # The index is not immediately equivalent - index duplicates
            # Just use the first value if index duplicated for now.
            df2 = df2.loc[~df2.index.duplicated(keep='first')]
            # Merge the files
            df = pd.concat([df, df2], axis=1)
        except KeyError:
            print('WARNING: no PCAS data found for {}'.format(flight_ID))

    # Add NOx as combined NO and NO2
    try:
        df['NOx'] = df['no_mr'].values + df['no2_mr'].values
    except:
        print('WARNING: failed to add NOx to flight {}'.format(flight_ID))
	# Resample the data?
    if resample_data:
        df = df.resample('1T' ).mean()
    # Add derived variables
    df = add_derived_variables2FAAM_data(df)
    return df


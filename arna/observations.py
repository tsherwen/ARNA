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
    # Read file
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
    interval_var = 'Sample interval'
    df[interval_var] = df[edate_var] - df[sdate_var]
    # Just use the middle of the timestep
    df.index = df[sdate_var] + (df[interval_var]/2)
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


def get_CIMS_data4flight(flight_ID='C225', resample_data=True, debug=False):
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


def add_biomass_flag2df(df, CIMSdf=None, flight_ID='C225', threshold=None):
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


def get_FAAM_core4flightnum(flight_ID='C225', version='v2020_06',
                            folder=None, resample_data=True):
    """
    Get the core FAAM flight data for a specific flight
    """
	# Where are the files? Which files to use?
    if isinstance(folder, type(None)):
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
            # Get values for surface area too (CDP/PCASP)
            ds = get_surface_area4flight(flight_ID=flight_ID, ds=ds,
                                         instrument='PCASP' )
            PCASPvar = 'PCASP-total-surface'
            ds = get_surface_area4flight(flight_ID=flight_ID, ds=ds,
                                         instrument='CDP' )
            CDPvar = 'CDP-total-surface'
            # Just include PCASP/CDP variables and a flag for PCASP
            df2 = ds[ [VarName, FlagName, PCASPvar] ].to_dataframe()
            df2.index = df2.index.floor('S')
            df2.index = df2.index.values # rm naming of time coord 'PCAS2TSPM'
            # The index is not immediately equivalent - index duplicates
            # Just use the first value if index duplicated for now.
            df2 = df2.loc[~df2.index.duplicated(keep='first')]
            # Repeat the process for the CDP time index
            # TODO: Find why 3 time indexes are used that are meant to be the same
            df3 = ds[ [CDPvar] ].to_dataframe()
            df3.index = df3.index.floor('S')
            df3.index = df3.index.values # rm naming of time coord 'PCAS2TSPM'
            # The index is not immediately equivalent - index duplicates
            # Just use the first value if index duplicated for now.
            df3 = df3.loc[~df3.index.duplicated(keep='first')]
            # Merge the files
            df = pd.concat([df, df2, df3], axis=1)
        except KeyError:
            print('WARNING: no PCASP data found for {}'.format(flight_ID))

    # Add NOx as combined NO and NO2
    try:
        df['NOx'] = df['no_mr'].values + df['no2_mr'].values
    except:
        print('WARNING: failed to add NOx to flight {}'.format(flight_ID))
    # Include flight_ID
    df['flight_ID'] =  flight_ID
	# Resample the data?
    if resample_data:
        df = df.resample('1T' ).mean()
    # Add derived variables
    df = add_derived_variables2FAAM_data(df)
    return df


def explore_FAAM_aerosol_data():
    """
    Explore the FAAM aerosol data (CDP, PCASP) from ARNA-2 campaign
    """
    # -- PCASP
    dsPCASP = get_FAAM_mineral_dust_calibration(instrument='PCASP',
                                                rtn_values=False)
    # -- CDP
    dsCDP = get_FAAM_mineral_dust_calibration(instrument='CDP',
                                              rtn_values=False)
    # only consider "potential dust" above a certain size?
    # Use 100 um for now


def get_FAAM_mineral_dust_calibration(instrument='PCASP', rtn_values=True):
    """
    Retrieve FAAM mineral dust calibration
    """
    # Location and name of calibration files?
    folder = '{}/FAAM/'.format( get_local_folder('ARNA_data'))
    if instrument == 'PCASP':
        #  NOTE: range ~0.1-4 microns
        filename = 'PCASP1_faam_20200128_v001_r000_cal.nc'
        # NOTE: dust values are a nc subgroup!
    #    group = 'bin_cal'
        group = 'bin_cal/mineral_dust'
    #    group = 'flow_cal'
       # The real part of the refractive index was taken as 1.53 which is a common value and is in the OPAC database. It is quite a bit smaller than the 1.547 that was reported by Weinzierl et al. [2011] but has been shown to have a relatively weak effect on the instrument response. The values of the imaginary part were based on references in  Ryder et al. [2019] along with the frequency distribution of k(550nm) presented in fig 9 of Ryder et al. [2013]. So the minimum value was extended from 0.0015i to 0.001i. Calculating the bin boundaries with these multiple Mie curves was done with Gaussian centre-weighted averaging with 0.001i and 0.0024i being +/-2 sigma extreme values.
    elif instrument == 'CDP':
        #  NOTE: range ~4-120 microns
        filename = 'CDP1_faam_20200208_v001_r000_cal.nc'
        # NOTE: dust values are a nc subgroup!
        group = 'master_cal/mineral_dust'
    # Open and return the widths and
    ds = xr.open_dataset( folder+filename, group=group )
    # Get values for bin centres and widths in microns (1E-6 metres)
    BinWidths = ds['dia_width'].values.flatten()
    BinCentres = ds['dia_centre'].values.flatten()
    d = {'BinWidths': BinWidths, 'BinCentres': BinCentres}
    if rtn_values:
        return d
    else:
        return ds


def get_surface_area4flight(flight_ID='C225', instrument='PCASP',
                            plt_up_values=False, ds=None):
    """
    Assuming spherical shape, calculate surface area
    """
    # - Consider the surface area Exclude the first channel
    # Retrieve calibration for mineral dust
    d = get_FAAM_mineral_dust_calibration(instrument=instrument)
    # Get the bin widths and centres (and convert to in metres)
    BinWidths = d['BinWidths'] *1E-6
    BinCentres = d['BinCentres'] *1E-6
    # What is the (max) radius of something in a given bin
#    R = BinCentres + BinWidths/2 # Assume all max radius of bin?
    R = BinCentres # assume all have same radius as middle of bin
    # Surface area (units microns^2 / binned particule)
    S = 4*np.pi*R**2

    # Get the FAAM dataset for specific flight if not provided
    if isinstance(ds, type(None)):
        folder = '{}/CEDA/{}/'.format( get_local_folder('ARNA_data'), version)
        files2use = glob.glob( folder + '*_{}*'.format(flight_ID.lower()) )
        cloud_phy_fnames = [ i for i in files2use if 'cloud-phy' in i ]
        # Just include the main file that ends "<flight_ID>.nc"
        suffix = '{}.nc'.format( flight_ID.lower() )
        cloud_phy_fname = [i for i in cloud_phy_fnames if i.endswith(suffix)]
        if len(cloud_phy_fname) >=1:
            ass_str = 'More than one main cloud phys file found for flight!'
            assert len(cloud_phy_filename) == 1, ass_str
        try:
            ds = xr.open_dataset(cloud_phy_filename[0])
        except KeyError:
            pstr = 'WARNING: no {} data found for {}'
            print(pstr.format(instrument, flight_ID))
    # Now extract PCASP/CDP data
    # Variable prefix?
    if instrument == 'PCASP':
        # ‘corrected’ for aspiration efficiency of particles entering the inlet based on some broad assumptions
        VarStr = 'PCAS2_{:0>2}'
        # not corrected for aspiration efficiency
#        VarStr = 'PCAS_{}_u'
        FlagName = 'PCAS2_FLAG'
        TimeVar = 'PCAS2TSPM'
    elif instrument == 'CDP':
        VarStr = 'CDP_{:0>2}'
        TimeVar = 'CDP_TSPM'
        FlagName = 'CDP_FLAG'
    # "It is traditional to skip bin 1 as the lower boundary is ‘undefined’, it is also where all the electronic noise ends up."
    # Values to use?
    range2use = np.arange(2, 31)
    vars2use = [VarStr.format(i) for i in range2use]
    # Remove flagged values for PCASP
    for var in vars2use:
        ds = set_flagged_data2NaNs(ds, VarName=var, FlagName=FlagName)
    # Now apply to get surface area by bin
    suffix = 'surface'
    for n_channel, channel in enumerate(range2use):
        Python_idx = channel-1
        ChannelName = vars2use[n_channel]
        VarName = '{}_{}'.format(vars2use[n_channel], suffix)
        surface = S[Python_idx]
        ds[VarName] = ds[ChannelName]*surface
        attrs = {
        'units': 'm3/cm-3',
        'long_name': 'Surface area of {}'.format(ChannelName),
        }
        ds[VarName].attrs = attrs
    # Plot this up to sanity check
    if plt_up_values:
        vars2plt = ['{}_{}'.format(i, suffix) for i in vars2use]
        vals2plt = ds[vars2plt].copy().mean(dim=TimeVar).to_array().values
        plt.bar(BinCentres[1:], vals2plt, width=BinWidths[1:] )
        ax = plt.gca()
        plt.yscale('log')
        units = 'm${^3}$/cm$^{-3}$'
        plt.ylabel( 'Binned surface area ({})'.format(units) )
        plt.title('Surface Area during ARNA-2')
        AC.save_plot('ARNA_aerosol_area_{}'.format(instrument))
        plt.close()
    # Then give a total and bin this value by <2.5um and >2.5um
    VarName = '{}-total-surface'.format(instrument)
    vars2use = ['{}_{}'.format(i, suffix) for i in vars2use]
    ds[VarName] = ds[vars2use[0]].copy()
    for var2use in vars2use[1:]:
        ds[VarName].values = ds[VarName].values + ds[var2use].values

    #
#    vars2rtn = ['{}-total-surface'.format(instrument)]
    return ds


def mk_file_of_flags():
    """
    Make csv files of flagging (e.g. SLR, dust) for ARNA flights
    """
    # Which flights to plot?
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
	# Just use non-transit ARNA flights
    flights_nums = [
    217,
    218, 219, 220, 221, 222, 223, 224, 225,
    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # Loop by flight and retrieve flights
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_FAAM_core4flightnum(flight_ID=flight_ID )
        df['flight_ID'] = flight_ID
        dfs_obs[flight_ID] = df

    # Only include a few variables for spaces
    vars2use = [
    'IS_DUST', 'IS_SLR', 'PS_RVSM', 'PS_RVSM_FLAG',
    'LAT_GIN', 'LAT_GIN_FLAG',
    'LON_GIN', 'LON_GIN_FLAG',
    'ALT_GIN', 'ALT_GIN_FLAG', 'flight_ID',
    ]
    version = 'v1_0'
    filename = 'ARNA_FAAM_flightpath_flagging_{}_{}.csv'
    for flight_ID in flight_IDs:
        df = dfs_obs[flight_ID]
        df[vars2use].to_csv(filename.format(flight_ID, version))
    # Save all the values together
    df = pd.concat([dfs_obs[i] for i in dfs_obs.keys()], axis=1)
    df[vars2use].to_csv(filename.format('ALL', version))


def get_FAAM_flights_df():
    """
    Retrieve DataFrame of FAAM BAe146 flights
    """
    # Flights to use...
    DataRoot = get_local_folder('DataRoot')
    folder = '/{}/FAAM/core_faam_NetCDFs/'.format(DataRoot)
    filename = 'FAAM_BAe146_Biomass_burning_and_ACISIS_flights_tabulated.csv'
    df = pd.read_csv(folder+filename)
    # Only consider the dates after Jan 2018
    DateVar = 'Date'
    df[DateVar] = pd.to_datetime( df[DateVar] )
    return df


def get_flighttracks4campaign(campaign='ARNA-2', PressVar="PS_RVSM",
                              LonVar='LON_GIN', TimeVar='Time',
                              LatVar='LAT_GIN', resample_data=True ):
    """
    Flight tracks for campaign
    """
    # Get dataframe of all flights and select those for a given campaign
    DataRoot = get_local_folder('DataRoot')
    folder = '/{}/FAAM/core_faam_NetCDFs/'.format(DataRoot)
    df = get_FAAM_flights_df()
    flight_IDs = df.loc[df['Campaign']==campaign, :]['Flight ID']
    # For flight in campaign flights
    dfs = []
    for flight_ID in flight_IDs:
        try:
            # Retrieve FAAM BAe146 Core NetCDF files
            filename = 'core_faam_*_{}_1hz.nc'.format(flight_ID.lower())
            file2use = glob.glob(folder+filename)
            if len(file2use) > 1:
                print('WARNING: more that one file found! (so using latest file)')
                print(file2use)
            ds = xr.open_dataset( file2use[0] )
            # Only select the variable of intereest and drop where these are NaNs
            df = ds[ [PressVar, LatVar, LonVar, TimeVar] ].to_dataframe()
            df = df.dropna()
            # Remove the index name and add index values to a column
            df.index.name = None
            dfs += [df.copy()]
            del ds, df
        except IndexError:
            pstr = 'WARNING: failed to extract flight data for {} ({})'
            print(pstr.format(flight_ID, campaign))
    # Return
    df = pd.concat(dfs, axis=0)
    if resample_data:
        df = df.resample('1T' ).mean()
    return df



def mk_planeflight_files4FAAM_campaigns(testing_mode=False):
    """
    Make plane-flight input files for various FAAM campaigns
    """
    # Location of flight data
    DataRoot = get_local_folder('DataRoot')
    folder = '/{}/FAAM/core_faam_NetCDFs/'.format(DataRoot)
    folder4csv = '/{}/FAAM/GEOSChem_planeflight_inputs/'.format(DataRoot)
    df = get_FAAM_flights_df()
#    if testing_mode:
    # Only consider flights in 2020
    DateVar = 'Date'
    df = df.loc[df[DateVar]  > datetime.datetime(2020,1,1), :]
    # flights to use?
    flight_IDs2use = [
   'C227',
   'C225', # Just re-done
   'C223', 'C224', 'C222', 'C221', 'C219', 'C220', 'C218', 'C217',
   'C212', 'C211', 'C210', 'C209', 'C208',
   'C207',
   'C206', 'C205', 'C204', 'C203', 'C202',
   'C201',
   'C200', 'C199',
   'C190', 'C189', 'C188', 'C187', 'C186', 'C185',
   'C184', 'C183',
   'C182',
   'C181',
   'C180', 'C179', 'C178', 'C145', 'C144',
   'C143',
   'C142', 'C141',
   'C140',
   'C139', 'C134', 'C133', 'C132', 'C129',
   'C106', 'C105', 'C104', 'C103'
    ]
    print(flight_IDs2use)
#    df = df.loc[ df['Flight ID'].isin(flight_IDs2use),:]
    # Extract variables of interest
    flight_IDs = df['Flight ID'].values
    #
    num_tracers = 203
#    rxn_nums = [
    # Just extract, HNO3, NO2, O3, NIT(s) for now
    # Numbers in *.eqn file for
#    1, 2, 3, 11, 16, 130, 131, 132, 133,
    # Index numbers for RCONST array
#    699, 700, 701, 702, 707, 803, 804, 805, 806
#    ]
    # Mannually setup slist
    met_vars = [
        'GMAO_ABSH', 'GMAO_PSFC', 'GMAO_SURF', 'GMAO_TEMP', 'GMAO_UWND',
        'GMAO_VWND', 'GMAO_PRES'
    ]
    species = ['OH', 'HO2']
    assert isinstance(num_tracers, int), 'num_tracers must be an integer'
    slist = ['TRA_{:0>3}'.format(i) for i in np.arange(1, num_tracers+1)]
    # Add Jvals
    JVN2use = np.arange(1,139)
    JVN2drop = [4, 5, 35, 52, 57, 58, 102]
    JVN2use = [i for i in JVN2use if i not in JVN2drop]
    JVAL_list = ['JVL_{:0>3}'.format(i) for i in JVN2use]
    slist = slist + species + met_vars + JVAL_list

    # Loop and extract FAAM BAe146 flights
    for flight_ID in flight_IDs:
        print(flight_ID)
        AC.mk_planeflight_input4FAAM_flight(folder=folder,
                                            folder4csv=folder4csv,
                                            testing_mode=testing_mode,
                                            num_tracers=num_tracers,
#                                            rxn_nums=rxn_nums,
                                            slist=slist,
                                            flight_ID=flight_ID,)
        gc.collect()


    # - Re-process the dates with two flights on one day
    # 1st double flight
    flight_IDs2use_l = [
        ['C199', 'C200'],
        ['C202', 'C203'],
        ['C204', 'C205'],
        ['C206', 'C207'],
        ['C208', 'C209'],
        ['C211', 'C212'],
        ['C223', 'C224'],
        ['C219', 'C220']
    ]
    for flight_IDs2use in flight_IDs2use_l:
        ds_l = []
        for flight_ID in flight_IDs2use:
            filename = 'core_faam_*_{}_1hz.nc'.format(flight_ID.lower())
            file2use = glob.glob(folder+filename)
            print(flight_ID, file2use)
            ds_l += [ xr.open_dataset( file2use[0] ) ]
        # Now re-make file but with values for both flights
        ds = xr.concat(ds_l, dim='Time')
        AC.mk_planeflight_input4FAAM_flight(ds=ds,
                                            folder=folder,
                                            folder4csv=folder4csv,
                                            testing_mode=testing_mode,
                                            num_tracers=num_tracers,
#                                            rxn_nums=rxn_nums,
                                            slist=slist,
                                            flight_ID=flight_ID,)
        gc.collect()
        del ds, ds_l





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
from datetime import timedelta
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
        df.loc[:, AltVar] = AC.hPa_to_Km(df[AltVar].values/1E3, reverse=True, )
    # Drop where there are not values for all coordinates
    if drop_NaNs:
        df = df.dropna()
    return df


def get_ARNA_flights_as_dfs():
    """
    Retrieve the ARNA flights as a list of dataframes
    """
    flight_nums = [216, 217, 218, 219, 220, 221, 222, 223, 224, 225]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    dfs = {}
    for flight_ID in flight_IDs:
        print(flight_ID)
        try:
            df = AC.get_FAAM_locations_as_df(flight_ID=flight_ID)
            dfs[flight_ID] = df
        except:
            print('WARNING: failed for {}'.format(flight_ID))
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
    ds[VarName].loc[bool] = np.NaN
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
    df.index = pd.to_datetime(df[var2use].values, format=format)
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
        'benzene': 'BENZ',  # GEOSChem, but not GEOS-CF output
        #    'benzenechb':,
        'cis_2_butene': 'PRPE',  # NOTE: lumped tracer for >= C3 alkenes
        'cyclo_pentane': 'ALK4',  # NOTE: lumped tracer for >= C4 Alkanes
        'dms': 'DMS',  # GEOSChem, but not GEOS-CF output
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
        'iso_butane': 'ALK4',  # NOTE: lumped tracer for >= C4 Alkanes
        'iso_butene': 'PRPE',  # NOTE: lumped tracer for >= C3 alkenes
        'iso_pentane': 'ALK4',  # NOTE: lumped tracer for >= C4 Alkanes
        'isoprene': 'PRPE',  # NOTE: lumped tracer for >= C3 alkenes
        'methanol': 'MOH',
        'mp_xylene': 'XYLE',
        'n_butane': 'ALK4',  # NOTE: lumped tracer for >= C4 Alkanes
        'n_heptane': 'ALK4',  # NOTE: lumped tracer for >= C4 Alkanes
        'n_hexane': 'ALK4',  # NOTE: lumped tracer for >= C4 Alkanes
        'n_octane': 'ALK4',  # NOTE: lumped tracer for >= C4 Alkanes
        'n_pentane': 'ALK4',  # NOTE: lumped tracer for >= C4 Alkanes
        'o_xylene': 'XYLE',
        'pent_1_ene':  'PRPE',  # NOTE: lumped tracer for >= C3 alkenes
        'propane': 'C3H8',
        'propene': 'PRPE',  # NOTE: lumped tracer for >= C3 alkenes
        'toluene': 'TOLU',
        #    'toluenechb':,
        'trans_2_butene': 'PRPE',  # NOTE: lumped tracer for >= C3 alkenes
        'trans_2_pentene': 'PRPE',  # NOTE: lumped tracer for >= C3 alkenes
    }
    # Invert the dictionary?
    if invert:
        d = {v: k for k, v in list(d.items())}
    return d[var]


def get_ARNA_flight_log_as_df():
    """
    Make a single pd.DataFrame with all flight summaries
    """
    flight_nums = [
        #    216,
        217, 218, 219, 220, 221, 222, 223, 224, 225
    ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    dfs = []
    for flight_ID in flight_IDs:
        dfs += [get_summary4flight(flight_ID=flight_ID)]
    # Combine and return as a single dataframe sorted by time
    df = pd.concat(dfs)
    df = df.sort_index()
    return df


def get_summary4flight(flight_ID='C217'):
    """
    Retrieve a FAAM flight summary as a dataframe
    """
    folder = '{}/{}/'.format(get_local_folder('ARNA_data'), 'CEDA/v2020_06')
    filename = 'flight-sum_faam_*_*_{}.csv'.format(flight_ID.lower())
    file2use = glob.glob(folder+filename)
    ass_str = 'WARNING: {} flight summaries found present for flight {}!'
    assert len(file2use) == 1, ass_str.format(file2use, flight_ID)
    # Add Gotcha for missing header in FAAM file archived at CEDA
    columns = [
        'Event', 'Start', 'Start Hdg / °', 'Start Hgt / kft', 'Start Lat / °',
        'Start Long / °', 'Stop', 'Stop Hdg / °', 'Stop Hgt / kft',
        'Stop Lat / °', ' Stop Long / °', 'Comment',
    ]
    if flight_ID == 'C217':
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
    df.index = pd.to_datetime(df[var2use].values, format=format)
    return df


def get_filters_data4flight(flight_ID='C217', all_flights=True):
    """
    Retrieve filters data from ARNA flights
    """
    # Where is the data?
    folder = '{}/{}/'.format(get_local_folder('ARNA_data'), 'Filters')
    # What is the name of the sheet in the excel file?
    filename = 'Aerosol_composition_for_submission.xlsx'
    xl = pd.ExcelFile(folder + filename)
    species = xl.sheet_names  # see all sheet names (species)
    dfs = [pd.read_excel(folder+filename, sheet_name=i) for i in species]
    # Now Just look at core data of interest
    TimeOnVar = 'Start_time'
    TimeOffVar = 'End_time'
    # Setup the main Dataframe, then add additional dataframes to this
    df = dfs[0]
    labels2drop = [
        TimeOnVar, TimeOffVar, 'Filter', 'Campaign', 'Flight', 'Airflow_stL',
        'Average_altitude_m', 'moles_m-3_in_air',
    ]
    for __df in dfs[1:]:
        __df = __df.drop(labels=labels2drop, axis=1)
        df = pd.concat([df, __df], axis=1)
    # Update sampling times to date times
    df[TimeOnVar] = pd.to_datetime(df[TimeOnVar].values)
    # End time
    edate_var = TimeOffVar
    df[TimeOffVar] = pd.to_datetime(df[TimeOffVar].values)
    # calculate mid point of sampling and set this as index
    interval_var = 'Sample interval'
    df[interval_var] = df[edate_var] - df[TimeOnVar]
    # Just use the middle of the timestep
    df.index = df[TimeOnVar] + (df[interval_var]/2)
    df = df.rename_axis(None)
    # Setup total values in ug/m3 (from nM/m3)
    units = 'nmoles_m-3'
    NewUnits = 'ug_m-3'
    UncertaintyStr = 'Total_{}_uncertainty_{}'
    SpeciesStr = 'Total_{}_{}'
    for spec in species:
        ObsVar = SpeciesStr.format(spec, units)
        UncertVar = UncertaintyStr.format(spec, units)
        SpecMass = AC.species_mass(spec)
        #  Add in new variables
        NewVar = SpeciesStr.format(spec, NewUnits)
        df[NewVar] = df[ObsVar].values / 1E9 * SpecMass * 1E6
        NewUncert = SpeciesStr.format(spec, NewUnits)
        df[NewUncert] = df[UncertVar].values / 1E9 * SpecMass * 1E6
    return df


def get_filters_data4flight_REDUNDANT(flight_ID='C217', all_flights=True):
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
    # NOTE: GEOS-CF has only sulfate + nitrate variables output
    # ('NIT', 'NITs', 'SO4', 'SO4s')
    # Observational variables (columns) to use
    NO3_var2use = ['NO3.2', 'NO3.5']
    SO4_var2use = ['SO4.2', 'SO4.5']
    Cl_var2use = ['Cl.2', 'Cl.5']
    NO2_var2use = ['NO2.2', 'NO2.5']
    C2O4_var2use = ['C2O4.2', 'C2O4.5']
    units = 'nanomoles/m3'

    # - Now chop off excess headers and make sure formats are correct
    df = dfFULL.loc[dfFULL.index[2:], :]
    # Just consider the available data in comparisons
    vars2use = NO3_var2use + SO4_var2use + Cl_var2use + NO2_var2use
    vars2use += C2O4_var2use
    df = df[CoreCols + vars2use]
    # Replace values less than black/NaNs with np.NaN
    df = df.replace('lower than blank', np.NaN)
    df = df.replace('NaN', np.NaN)
    # Remove blanks (as these are NaNs)
    df = df.rename_axis(None)
    df = df.loc[(df['height'] != 'blank').values, :]
    # Update sampling times to date times
    # Start time
    TimeOnVar = 'Time on'
    sdate_var = 'Sample Start'
    df[sdate_var] = df['Day'].astype(str) + ' ' + df[TimeOnVar].astype(str)
    format = '%Y-%m-%d %H:%M:%S'
    df[sdate_var] = pd.to_datetime(df[sdate_var].values, format=format)
    del df[TimeOnVar]
    # End time
    TimeOffVar = 'Time off'
    edate_var = 'Sample End'
    df[edate_var] = df['Day'].astype(str) + ' ' + df[TimeOffVar].astype(str)
    format = '%Y-%m-%d %H:%M:%S'
    df[edate_var] = pd.to_datetime(df[edate_var].values, format=format)
    del df[TimeOffVar]
    # calculate mid point of sampling and set this as index
    interval_var = 'Sample interval'
    df[interval_var] = df[edate_var] - df[sdate_var]
    # Just use the middle of the timestep
    df.index = df[sdate_var] + (df[interval_var]/2)
    df = df.rename_axis(None)
    # - Just consider totals for species of interest
    NO3_var = 'NO3.total'
    df[NO3_var] = df[NO3_var2use].sum(axis=1)
    SO4_var = 'SO4.total'
    df[SO4_var] = df[SO4_var2use].sum(axis=1)
    Cl_var = 'Cl.total'
    df[Cl_var] = df[Cl_var2use].sum(axis=1)
    NO2_var = 'NO2.total'
    df[NO2_var] = df[NO2_var2use].sum(axis=1)
    C2O4_var = 'C2O4.total'
    df[C2O4_var] = df[C2O4_var2use].sum(axis=1)
    del dfFULL
    # Convert to ug/m3 from 'nanomoles/m3'
    df[NO3_var] = df[NO3_var].values / 1E9 * AC.species_mass('NIT') * 1E6
    df[SO4_var] = df[SO4_var].values / 1E9 * AC.species_mass('SO4') * 1E6
    df[Cl_var] = df[Cl_var].values / 1E9 * AC.species_mass('Cl') * 1E6
    df[NO2_var] = df[NO2_var].values / 1E9 * AC.species_mass('NO2') * 1E6
    df[C2O4_var] = df[C2O4_var].values / 1E9 * AC.species_mass('C2O4') * 1E6
    # Return all flights unless a specific flight requested
    if all_flights:
        return df
    else:
        return df.loc[df['Flight'] == flight_ID, :]


def get_CIMS_data4flight(flight_ID='C225', resample_data=True,
                         verbose=False, debug=False):
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
        xl = pd.ExcelFile(folder + filename)
        if debug:
            print(xl.sheet_names)
        dt_var = 'Date:Time'
        dates = AC.dt64_2_dt(df[dt_var].values)
        dates = [i.strftime(format) for i in dates]
        dates = [datetime_.strptime(i, '%d/%m/%Y %H:%M:%S') for i in dates]
        df.index = dates
        del df[dt_var]
        # Merge the files
        dfM = pd.concat([dfM, df], axis="index")
    except:
        pstr = "WARNING: failed to include CIMS halogen data for '{}' in df"
        if verbose:
            print(pstr.format(flight_ID))
    # - Also get HNO3 / HONO (use the time coordinate for index)
    try:
        filename = 'ACSIS6_ARNA_1hz_HNO3_HONO_ToSend.xlsx'
        df = pd.read_excel(folder + filename, sheet_name=sheet_name)
        xl = pd.ExcelFile(folder + filename)
        if debug:
            print(xl.sheet_names)
        dt_var = 'date'
        dates = AC.dt64_2_dt(df[dt_var].values)
        dates = [i.strftime(format) for i in dates]
        dates = [datetime_.strptime(i, '%d/%m/%Y %H:%M:%S') for i in dates]
        df.index = dates
        del df[dt_var]
        # Merge the files
        dfM = pd.concat([dfM, df], axis="index")
    except:
        pstr = "WARNING: failed to include CIMS NOy data for '{}' in df"
        if verbose:
            print(pstr.format(flight_ID))

    # Remove the existing HNO3 variables and then include the latest data
    redundant_HNO3_vars = ['HNO3 (ppt*) ', 'HNO3 (ppt*)']
    pstr = "NOTE: Removed redundant HNO3 data ('{}') for flight '{}'"
    try:
        for var in redundant_HNO3_vars:
            if var in dfM.columns:
                del dfM[var]
                if verbose:
                    print(pstr.format(var, flight_ID))
    except KeyError:
        pass
    try:
        filename = 'ARNA2_HNO3_YORK_1Hz.csv'
        df = pd.read_csv(folder + filename)
        dt_var = 'date_time'
        format = '%Y-%m-%d %H:%M:%S'
        df.index = pd.to_datetime(df[dt_var].values, format=format)
        del df[dt_var]
        # Merge the files
        bool = df['flight_id'].values == flight_ID
        Nvals = df.loc[bool, :].shape[0]
        AssStr = "WARNING: No values found for CIMS HNO3 on flight '({})"
        assert Nvals > 0, AssStr.format(flight_ID)
        dfM = pd.concat([dfM, df.loc[bool, :]], axis="index")
    except:
        pstr = "WARNING: failed to include CIMS HNO3 data for '{}' in df"
        print(pstr.format(flight_ID))

    # - Also include HCN data
    try:
        filename = 'ACSIS6_ARNA2_HCN_James_TMS_Update.xlsx'
        df = pd.read_excel(folder + filename, sheet_name=sheet_name)
        xl = pd.ExcelFile(folder + filename)
        if debug:
            print(xl.sheet_names)
        dt_var = 'date_time'
    #    df2.index = pd.to_datetime( df2[dt_var].values, format='%m/%d/%Y %H:%M:%S')
        dates = AC.dt64_2_dt(df[dt_var].values)
        dates = [i.strftime(format) for i in dates]
        format = '%Y-%m-%d %H:%M:%S'
        dates = [datetime_.strptime(i, format) for i in dates]
        df.index = dates
        del df[dt_var]
        # Merge the files
        dfM = pd.concat([dfM, df], axis="index")
    except:
        pstr = "WARNING: failed to include CIMS HCN data for '{}' in df"
        print(pstr.format(flight_ID))

    # Update the variable names
    VarNameDict = {
        'BrO (ppt*)': 'BrO',
        'Br2 + HOBr (ppt*)': 'Br2+HOBr',
        'HONO (ppt*)': 'HONO',
        'HNO3 (ppt*) ': 'HNO3', 'HNO3 (ppt*)': 'HNO3', 'HNO3_ppt': 'HNO3',
        'HCN (ppt*)': 'HCN',
    }
    dfM = dfM.rename(columns=VarNameDict)
    # Include a flag for flight ID
    dfM['flight_ID'] = flight_ID
    # Resample the data?
    if resample_data:
        dfM = dfM.resample('1T').mean()
    return dfM


def get_CIMS_data4flight_TEMP_LOWER():
    """
    Temporality return the lower-upper limits for the HNO3 data
    """
    folder = '{}/{}/'.format(get_local_folder('ARNA_data'), 'CIMS')
    filename = 'ARNA2_HNO3_UPPER_LOWER_PPT.csv'
    df = pd.read_csv(folder+filename)
    # Add a datetime object index
    dt_var = 'date_time'
    format = '%d/%m/%Y %H:%M'
    df.index = pd.to_datetime(df[dt_var].values, format=format)
    del df[dt_var]
    return df


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
        df.loc[bool, VarName] = False
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
        sigma = stats.loc[stats.index == 'std'].values[0]
        threshold = background + (sigma*3)
    # Add a boolean for biomass burning (BB)
    df[NewVar] = False
    df.loc[df[var2use] >= threshold,  NewVar] = True
    return df


def get_FAAM_core4flightnum(flight_ID='C225', version='v2020_06',
                            folder=None, resample_data=True):
    """
    Get the core FAAM flight data for a specific flight
    """
    # Where are the files? Which files to use?
    if isinstance(folder, type(None)):
        folder = '{}/CEDA/{}/'.format(get_local_folder('ARNA_data'), version)
    files2use = glob.glob(folder + '*_{}*'.format(flight_ID.lower()))
    # Open all the files together
    # Cannot as time is not capitalised in both
#    ds = xr.open_mfdataset(files2use)
    #
    FAAM_filename = [i for i in files2use if 'core_faam' in i]
    asstr = 'More than one core FAAM file found for flight!'
    assert len(FAAM_filename) == 1, asstr
    ds = xr.open_dataset(FAAM_filename[0])
    ds = ds.rename({'Time': 'time'})
    # Set flagged data to NaNs
    ds = set_flagged_data2NaNs(ds, VarName='O3_TECO', FlagName='O3_TECO_FLAG')
    ds = set_flagged_data2NaNs(ds, VarName='CO_AERO', FlagName='CO_AERO_FLAG')
    # do the same for vertical velocity and roll
    # 'Roll angle from POS AV 510 GPS-aided Inertial Nav. unit (positive for left wing up)'
    VarName = 'ROLL_GIN'
    FlagName = 'ROLL_GIN_FLAG'
    ds = set_flagged_data2NaNs(ds, VarName=VarName, FlagName=FlagName)
    # 'Aircraft velocity down from POS AV 510 GPS-aided Inertial Navigation
    VarName = 'VELD_GIN'
    FlagName = 'VELD_GIN_FLAG'
    ds = set_flagged_data2NaNs(ds, VarName=VarName, FlagName=FlagName)
    # Convert to a dataframe
    df = ds.to_dataframe()
    df.index.rename(None)
    del ds
    # Check for the core NOx file
    Nitrates_filename = [i for i in files2use if 'core-nitrates' in i]
    if len(Nitrates_filename) >= 1:
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
    cloud_phy_filenames = [i for i in files2use if 'cloud-phy' in i]
    # Just include the main file that ends "<flight_ID>.nc"
    suffix = '{}.nc'.format(flight_ID.lower())
    cloud_phy_filename = [i for i in cloud_phy_filenames if i.endswith(suffix)]
    if len(cloud_phy_filename) >= 1:
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
                                         instrument='PCASP')
            PCASPvar = 'PCASP-total-surface'
            ds = get_surface_area4flight(flight_ID=flight_ID, ds=ds,
                                         instrument='CDP')
            CDPvar = 'CDP-total-surface'
            # Just include PCASP/CDP variables and a flag for PCASP
            df2 = ds[[VarName, FlagName, PCASPvar]].to_dataframe()
            df2.index = df2.index.floor('S')
            df2.index = df2.index.values  # rm naming of time coord 'PCAS2TSPM'
            # The index is not immediately equivalent - index duplicates
            # Just use the first value if index duplicated for now.
            df2 = df2.loc[~df2.index.duplicated(keep='first')]
            # Repeat the process for the CDP time index
            # TODO: Find why 3 time indexes are used that are meant to be the same
            df3 = ds[[CDPvar]].to_dataframe()
            df3.index = df3.index.floor('S')
            df3.index = df3.index.values  # rm naming of time coord 'PCAS2TSPM'
            # The index is not immediately equivalent - index duplicates
            # Just use the first value if index duplicated for now.
            df3 = df3.loc[~df3.index.duplicated(keep='first')]
            # Merge the files
            df = pd.concat([df, df2, df3], axis=1)
        except KeyError:
            print('WARNING: no PCASP data found for {}'.format(flight_ID))
    # Include the latest NOx data from Simone and Chris
    dfHONO = get_latest_NOx_HONO_data(flight_ID=flight_ID)
#    if not isinstance(dfHONO, type(None)):
    df = pd.concat([df, dfHONO])

    # Add NOx as combined NO and NO2
    try:
        df['NOx_mr'] = df['no_mr'].values + df['no2_mr'].values
    except:
        print('WARNING: failed to add NOx to flight {}'.format(flight_ID))
    try:
        df['NOx'] = df['NO_pptV'].values + df['NO2_pptV'].values
    except:
        print('WARNING: failed to add NOx to flight {}'.format(flight_ID))

    # Add a values for H2O
    try:
        ds = set_flagged_data2NaNs(ds, VarName='NV_LWC1_C',
                                   FlagName='NV_LWC1_C_FLAG')
        ds = set_flagged_data2NaNs(ds, VarName='NV_LWC2_C',
                                   FlagName='NV_LWC2_C_FLAG')

        ds['H2O'] = ds['NV_LWC1_C'].copy()
        ds['H2O'] += ds['NV_LWC2_C']
        ds['H2O_U'] = ds['NV_LWC1_C'].copy()
        ds['H2O_U'] += ds['NV_LWC2_C']
        # Air density  ( ρ =  P / RT )
#        ds = set_flagged_data2NaNs(ds, VarName='TAT_DI_R',
#                                       FlagName='TAT_DI_R_FLAG')
#        ds = set_flagged_data2NaNs(ds, VarName='PS_RVSM',
#                                       FlagName='PS_RVSM_FLAG')
#        GasConst = 8.314 4621 # m3 Pa K−1 mol−1
#        GasConst *= 0.01 # m3 hPa K−1 mol−1
#        GasConst /= AC.constants('AVG') # m3 hPa K−1
        GasConst = 0.167226  # J/ kg K
        ds['AIR_DEN'] = ds["PS_RVSM"]*10 / (GasConst * ds['TAT_DI_R'])
        # (kg/m3) => g/m3
        ds['AIR_DEN'] = ds['AIR_DEN']*1E3

        # convert g/M3 to ppbv
        ds['H2O'] = ds['H2O'] / ds['AIR_DEN']

    except:
        print('WARNING: failed to add H2O to flight {}'.format(flight_ID))

    # Include flight_ID
    df['flight_ID'] = flight_ID
    # Resample the data?
    if resample_data:
        df = df.resample('1T').mean()
    # Add derived variables
    df = add_derived_variables2FAAM_data(df)
    return df


def get_latest_NOx_HONO_data(flight_ID=None, version='v1'):
    """
    Use the latest NOx/HONO data from Simone and Chris
    """
    # Locations of data
    folder = '{}/{}/'.format(get_local_folder('ARNA_data'), 'FAAM')
    folder += 'ARNA_NOx_HONO_data/'
    FileStr = 'NOx_and_HONO'
    files2use = glob.glob('{}*{}*{}*'.format(folder, FileStr, version))
    files2use = list(sorted(files2use))
    # Process entire dataset or just one flight
    if isinstance(flight_ID, type(None)):
        dfs = []
        for file2use in files2use:
            df = pd.read_csv(file2use)
            flight_ID = file2use.split('/')[-1][:4]
            df['flight_ID'] = flight_ID
            dfs += [df]
        # Combine
        df = pd.concat(dfs)
    else:
        file2use = [i for i in files2use if flight_ID in i]
        if len(file2use) == 0:
            print('WARNING: No NOx/HONO data for {}'.format(flight_ID))
            return None
        else:
            assert len(file2use) == 1, 'STOPPING: more than one file!'
            file2use = file2use[0]
            df = pd.read_csv(file2use)
            flight_ID = file2use.split('/')[-1][:4]
            df['flight_ID'] = flight_ID
    # ... and make index the datetime
    df = df.sort_index()
    df.index.name = None
    df.index = pd.to_datetime(df[df.columns[0]])
    df.index = pd.to_datetime(df.index.values)
    del df[df.columns[0]]
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
    folder = '{}/FAAM/'.format(get_local_folder('ARNA_data'))
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
    ds = xr.open_dataset(folder+filename, group=group)
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
    BinWidths = d['BinWidths'] * 1E-6
    BinCentres = d['BinCentres'] * 1E-6
    # What is the (max) radius of something in a given bin
#    R = BinCentres + BinWidths/2 # Assume all max radius of bin?
    R = BinCentres  # assume all have same radius as middle of bin
    # Surface area (units microns^2 / binned particule)
    S = 4*np.pi*R**2

    # Get the FAAM dataset for specific flight if not provided
    if isinstance(ds, type(None)):
        folder = '{}/CEDA/{}/'.format(get_local_folder('ARNA_data'), version)
        files2use = glob.glob(folder + '*_{}*'.format(flight_ID.lower()))
        cloud_phy_fnames = [i for i in files2use if 'cloud-phy' in i]
        # Just include the main file that ends "<flight_ID>.nc"
        suffix = '{}.nc'.format(flight_ID.lower())
        cloud_phy_fname = [i for i in cloud_phy_fnames if i.endswith(suffix)]
        if len(cloud_phy_fname) >= 1:
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
        plt.bar(BinCentres[1:], vals2plt, width=BinWidths[1:])
        ax = plt.gca()
        plt.yscale('log')
        units = 'm${^3}$/cm$^{-3}$'
        plt.ylabel('Binned surface area ({})'.format(units))
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


def mk_file_of_flags(flights_nums=[]):
    """
    Make csv files of flagging (e.g. SLR, dust) for ARNA flights
    """
    # Which flights to plot?
    # Just use non-transit ARNA flights
    if len(flights_nums) == 0:
        flights_nums = [
            217, 218, 219, 220, 221, 222, 223, 224, 225,
        ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # Loop by flight and retrieve flights
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_FAAM_core4flightnum(flight_ID=flight_ID)
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
    df[DateVar] = pd.to_datetime(df[DateVar])
    return df


def get_flighttracks4campaign(campaign='ARNA-2', PressVar="PS_RVSM",
                              LonVar='LON_GIN', TimeVar='Time',
                              LatVar='LAT_GIN', resample_data=True):
    """
    Flight tracks for campaign
    """
    # Get dataframe of all flights and select those for a given campaign
    DataRoot = get_local_folder('DataRoot')
    folder = '/{}/FAAM/core_faam_NetCDFs/'.format(DataRoot)
    df = get_FAAM_flights_df()
    flight_IDs = df.loc[df['Campaign'] == campaign, :]['Flight ID']
    # For flight in campaign flights
    dfs = []
    for flight_ID in flight_IDs:
        try:
            # Retrieve FAAM BAe146 Core NetCDF files
            filename = 'core_faam_*_{}_1hz.nc'.format(flight_ID.lower())
            file2use = glob.glob(folder+filename)
            pstr = 'WARNING: more that one file found! (so using latest file)'
            if len(file2use) > 1:
                print(pstr)
                print(file2use)
            ds = xr.open_dataset(file2use[0])
            # Only select the variable of intereest and drop where these are NaNs
            df = ds[[PressVar, LatVar, LonVar, TimeVar]].to_dataframe()
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
        df = df.resample('1T').mean()
    return df


def mk_planeflight_files4sites(testing_mode=False):
    """
    Make plane-flight input files for various ground sites
    """
    # Location of flight data
    TYPE = ['CVO{}'.format(i) for i in range(1, 8)]
    #
    sdate = datetime.datetime(2015, 1, 1,)
    edate = datetime.datetime(2021, 1, 1,)
    dates = pd.date_range(sdate, edate, freq='T')
    # Get list of species
    num_tracers = 203
    slist = get_planeflight_slist2output(num_tracers=num_tracers)

    # for each location make a DataFrame, then conbime
    dfs = []
    for n, type_ in enumerate(TYPE):
        # Get locations
        LON, LAT, ALT = AC.get_loc(type_)
        PRESS = 1013.25  # AC.hPa_to_Km([ALT/1E3], reverse=True, )
        print(n, type_, LON, LAT, ALT)
        # dictionary of data
        nvar = len(dates)
        d = {
            'datetime': dates, 'LAT': [LAT]*nvar, 'LON': [LON]*nvar,
            'TYPE': [type_]*nvar, 'PRESS': [PRESS]*nvar}
        dfs += [pd.DataFrame(d, index=np.arange(nvar)+(n*1E6))]
    # combine all TYPE (sites) and sort by date
    df = pd.concat(dfs).sort_values('datetime', ascending=True)

    # Now print as files
    AC.prt_PlaneFlight_files_v12_plus(df=df, slist=slist,
                                      Extra_spacings=Extra_spacings)


def get_planeflight_slist2output(num_tracers=None, folder=None,
                                 OutputJVALs=False):
    """
    Store of planeflight slist to request outputs for
    """
    # Use the number of tracers in the input.geos file, unless specified
    if isinstance(num_tracers, type(None)):
        if isinstance(folder, type(None)):
            TRAs = AC.get_specieslist_from_input_geos(folder=folder)
            num_tracers = len(TRAs)
        else:
            num_tracers = 203  # Use value for GEOS-Chem v12.9.0
    # Mannually setup slist
    met_vars = AC.get_MetAer_vars2use4PF()
    species = AC.get_Species_vars2use4PF()
    assert isinstance(num_tracers, int), 'num_tracers must be an integer'
    slist = ['TRA_{:0>3}'.format(i) for i in np.arange(1, num_tracers+1)]
    slist = slist + species + met_vars
    # Add Jvals
    if OutputJVALs:
        JVN2use = np.arange(1, 139)
        JVN2drop = [4, 5, 35, 52, 57, 58, 102]
        JVN2use = [i for i in JVN2use if i not in JVN2drop]
        JVAL_list = ['JVL_{:0>3}'.format(i) for i in JVN2use]
        slist += JVAL_list
    return slist


def mk_planeflight_files4FAAM_campaigns(folder=None, testing_mode=False,
                                        num_tracers=None, OutputJVALs=True):
    """
    Make plane-flight input files for various FAAM campaigns
    """
    # Location of flight data
    DataRoot = get_local_folder('DataRoot')
    folderFAAMnetCDF = '/{}/FAAM/core_faam_NetCDFs/'.format(DataRoot)
#    folder4csv = '/{}/FAAM/GEOSChem_planeflight_inputs.ACID/'.format(DataRoot)
    df = get_FAAM_flights_df()
#    if testing_mode:
    # Only consider flights in 2020
    DateVar = 'Date'
    df = df.loc[df[DateVar] > datetime.datetime(2020, 1, 1), :]
    # flights to use?
    flight_IDs2use = [
        'C227',
        'C225',  # Just re-done
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
    # Get list of species to output from plane flight diagnostic
    slist = get_planeflight_slist2output(folder=folder,
                                         num_tracers=num_tracers,
                                         OutputJVALs=OutputJVALs)
    # Loop and extract FAAM BAe146 flights
    for flight_ID in flight_IDs:
        print(flight_ID)
        AC.mk_planeflight_input4FAAM_flight(folder=folderFAAMnetCDF,
#                                            folder4csv=folder4csv,
                                            testing_mode=testing_mode,
                                            num_tracers=num_tracers,
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
            file2use = glob.glob(folderFAAMnetCDF+filename)
            print(flight_ID, file2use)
            ds_l += [xr.open_dataset(file2use[0])]
        # Now re-make file but with values for both flights
        ds = xr.concat(ds_l, dim='Time')
        AC.mk_planeflight_input4FAAM_flight(ds=ds,
                                            folder=folderFAAMnetCDF,
#                                            folder4csv=folder4csv,
                                            testing_mode=testing_mode,
                                            num_tracers=num_tracers,
                                            slist=slist,
                                            flight_ID=flight_ID,)
        gc.collect()
        del ds, ds_l
    gc.collect()


def read_FIREXAQ_files(path, folder='merge', var=''):
    """
    Read files of observational data from FIREX-AQ campaign

    Notes
    ----
     - Original function credit: Matt Rowlinson
    """
    df_list=[]
    flag_list=[]
    files2use = sorted(glob.glob(f'{path}/{folder}/*{var}*.ict'))
    for infileN, infile in enumerate( files2use ):
        with open(infile) as thefile:
            try:
                header= np.array([next(thefile) for x in range(90) ])
            except:
                continue
            start = header[6].replace(',',' ').split()
            start_date = datetime_( int( start[0] ),
                                    int( start[1] ),
                                    int( start[2] ))
        # Find where the header ends and values begin - manually narrowed down
        for nskip in range(675,680):
            try:
                fh = np.loadtxt(infile, skiprows=nskip, delimiter=',')
                break
            except:
                continue
        thefile = open(infile,'r')
        c = thefile.readlines()
        column_names = c[nskip-1].replace(' ','').split(',')
        df = pd.DataFrame(fh, index=fh[:,0], columns=column_names)

        # Use a different approach for
        if (var=='thru'):
            df = find_FIREXAQ_times(df, start_date, UseTimeStart=True)
        else:
            df = find_FIREXAQ_times(df, start_date, UseTimeStart=False)
        # Include the RF from the file name #
        # NOTE: research flight (RF) ID not included in filename or files,
        #       so using filenumber instead. This will not work if reading
        #       merge file (var = 'thru').
        df['FileNumber'] = infileN

        df_list.append(df)
    df = pd.concat(df_list)
    return df


def find_FIREXAQ_times(df, t0, UseTimeStart=False):
    """
    Find FIREX-AQ times

    Notes
    ----
     - Original function credit: Matt Rowlinson
    """
    timex = []
    # Only will work for unmerged files.
    if UseTimeStart:
#        tstamp = df[df.columns[1]]
        tstamp = df['Time_Start']
        for nDate in range(len(tstamp)):
            # If using 'Time_Start' (in seconds)
            timex.append(t0 + timedelta(seconds=tstamp.values[nDate]))

    else:
        # NOTE: units are fractional days, not seconds.
        # Subtract int of first day
        tstamp = df['Fractional_Day']
        day0 = int(tstamp.values[0] )
        days_since_t0 = tstamp.values - day0
        # If using fractional day
        timex = [ AC.add_days(t0, i) for i in days_since_t0 ]
#        for nDate in range(len(tstamp)):
#            timex.append(AC.add_days(t0, days_since_t0[nDate]))
    df.index = timex
    return df


def get_FIREX_AQ_from_ICT_files(UseMergeFile=False):
    """
    Extract the merge file for FIREX-AQ
    """
    # Local variables
#    Folder = '/users/ts551/scratch/data/ARNA/TOMAS_FIREX/merge/'
#    FileName = 'firexaq-mrg60-dc8_merge_20190722_R1_thru20190905.ict'
    path = '/mnt/lustre/groups/chem-acm-2018/shared_data/FIREX-AQ/merge/'
    #
    if UseMergeFile:
        FileName = 'firexaq-mrg60-dc8_merge_20190722_R1_thru20190905.ict'
        # NOTE: below needs to be updated to use merged file
        dfM = read_FIREX_ICT_file(path, FileName)
    # Load files induvidually and merge these into a single file to use
    else:
        files2use = sorted(glob.glob('{}/*.ict*'.format(path)))
        # Skip the merge file
        files2use = [i for i in files2use if 'thru' not in i]

        # Loop files and combine into a single dataframe
        for nfile2use, file2use in enumerate(files2use):
            FileName = file2use.split(path)[-1]
            print(file2use, FileName)
            df = read_FIREX_ICT_file(path, FileName,)
            # Save the filename as a proxy for research flight number
            df['FileNumber'] = nfile2use + 1
            # Merge with the other dataframes
            if nfile2use == 0:
                dfM = df
            else:
                dfM = pd.concat([dfM,df])#,ignore_index=True)
    return dfM


def read_FIREX_ICT_file(path, FileName):
    """
    Read and return ICT files as DataFrames
    """
    # Setup a manual file reader for the ICT files.
    file2use = '{}/{}'.format(path, FileName)
    # Local variables
    HeaderLineStarts = 'Time_Start, Time_Stop, Day_Of_Year_YANG, Latitude_YANG'
    Year = 2019
    FirstDayOfYear = datetime.datetime(Year, 1, 1)
    DOYvar = 'Day_Of_Year_YANG'
    StartVar = 'Time_Start'
    # Extract file by reading line by line
    with open( file2use, 'r') as OpenedFile:

        # Read data after the head line has been read
        ReadDataHereOnwards = False
        data = []
        for line in OpenedFile:
            line = line.strip()
            # Extract data after header
            if ReadDataHereOnwards:
                data += [line.split(',')]
            # skip lines until header for data found
            if line.startswith(HeaderLineStarts):
                header = line.split(',')
                header = [i.strip() for i in header]
                ReadDataHereOnwards = True

        # Compile data and header into a pd.DataFrame
        df = pd.DataFrame(data, columns=header)
        # convert columns to floats where possible
        for col in df.columns:
            df.loc[:, col] = pd.to_numeric(df[col])

        # Update the index to be in datetime
        dates = []
        days = df[DOYvar].values
        for idx in df.index:
            day = df.loc[idx, DOYvar]
            seconds = df.loc[idx, StartVar]
            date = FirstDayOfYear + datetime.timedelta(int(day) - 1.0)
            date = AC.add_secs(date, seconds)
            dates += [date]
        df.index = dates
    return df


def get_FIREX_AQ_data(debug=False, RtnAllData=True,
                      FilterPollutedAirMasses=True,
                      RmObsBelowGround=True,
                      UpdateTimeeZone2LocalTime=True,
                      FilterByTimeOfDay=True,
                      SetFlaggedDataToNaN=True,
                      stime='10:00', etime='15:00'):
    """
    Retrieve FIREX-AQ data as a pandas DataFrame

    Notes
    ----
     - Original function credit: Matt Rowlinson
    """
    firex_vars = Get_FIREXAQ_variable_dict()
    keys = firex_vars.keys()
    # Read FIREX-AQ data
    path = '/mnt/lustre/groups/chem-acm-2018/shared_data/FIREX-AQ'
    # NOTE: if the merge file will be read if var='thru'
#    df0 = read_FIREXAQ_files(path, var='thru')
    df0 = get_FIREX_AQ_from_ICT_files(UseMergeFile=False)

    # Convert timezone and apply restrictions on data
    if UpdateTimeeZone2LocalTime:
        df0.index = df0.index.tz_localize('UTC').tz_convert('US/Pacific')
    if FilterByTimeOfDay:
        if not UpdateTimeeZone2LocalTime:
            print('WARNING: Selecting time of day in UTC, not local time')
        df0 = df0.between_time(stime, etime)
    # Filter out polluted air mass
    if FilterPollutedAirMasses:
        df0 = df0[df0['CO_DACOM_DISKIN'] < 100. ]
    # Filter out flagged data?
    if RmObsBelowGround:
        df0 = df0[df0['MSL_GPS_Altitude_YANG'] > 0. ]
    # Flag the data here that
    if SetFlaggedDataToNaN:
        FlagValue = -999999.000000
        for col in df0.columns:
            df0.loc[ df0[col] == FlagValue, col] = np.NaN

    # Return entire dataset or just a single species?
    if RtnAllData:
        return df0
    else:
        for var in keys:
            print( var )
            df = pd.concat([ df0[firex_vars[var]['firex']],
                             df0['Latitude_YANG'],
                             df0['Longitude_YANG'],
                             df0['MSL_GPS_Altitude_YANG'] ],
                             axis=1 )
            df.columns = [var,'Latitude','Longitude','Altitude']
            df = df[ df[var] > 0. ]
            if debug:
                print( df )

            # Save species to dictionary
            dfs[var] = df.copy()
        print(df)



def Get_FIREXAQ_variable_dict():
    """
    Function to store variables for FIREX-AQ campaign

    Notes
    ----
     - Original function credit: Matt Rowlinson
    """
    firex_vars = { 'EOH'  : {
                            'firex' : 'C2H5OH_TOGA_APEL',
                            'gc'    : 'SpeciesConc_EOH',
                            'conv'  : False,
                            'scale' : 1e12},
                   'CH4'  : {
                            'firex' : 'CH4_DACOM_DISKIN',
                            'gc'    : 'SpeciesConc_CH4',
                            'conv'  : False,
                            'scale' : 1e9},
                   #'HNO3-NO3' : {
                   #         'firex' : 'HNO3+submicron-NO3_SAGA_DIBB',
                   #         'gc'    : ['HNO3','NIT','NITD1','NITD2'],
                    #        'conv'  : False,
                    #        'scale' : 1e12},
                   'C2H6' : {
                            'firex' : 'Ethane_WAS_BLAKE',
                            'gc'    : 'SpeciesConc_C2H6',
                            'conv'  : False,
                            'scale' : 1e12},
                   'C3H8' : {
                            'firex' : 'Propane_WAS_BLAKE',
                            'gc'    : 'SpeciesConc_C3H8',
                            'conv'  : False,
                            'scale' : 1e12},
                   'BENZ' : {
                            'firex' : 'Benzene_WAS_BLAKE',
                            'gc'    : 'SpeciesConc_BENZ',
                            'conv'  : False,
                            'scale' : 1e12},
                   'TOLU' : {
                            'firex' : 'Toluene_WAS_BLAKE',
                            'gc'    : 'SpeciesConc_TOLU',
                            'conv'  : False,
                            'scale' : 1e12},
                   'HNO2' : {
                            'firex' : 'HNO2_NOAACIMS_VERES',
                            'gc'    : 'SpeciesConc_HNO2',
                            'conv'  : False,
                            'scale' : 1e12 },
                   'NH4' : {
                            'firex' : 'NH4_ug/m3_DIBB',
                            'gc'    : 'SpeciesConc_NH4',
                            'conv'  : True,
                            'scale' : 1e12,
                            'mm'    : 18.04},
                   'SO4' : {
                            'firex' : 'SO4_ug/m3_DIBB',
                            'gc'    : 'SpeciesConc_SO4',
                            'conv'  : True,
                            'scale' : 1e9,
                            'mm'    : 96.06},
                   'NIT-all' : {
                            'firex' : 'NO3_ug/m3_DIBB',
                            'gc'    : ['NIT','NITs','NITD1','NITD2','NITD3','NITD4'],
                            'conv'    : True,
                            'scale'  : 1e9,
                            'mm' : 62.0049 },
                   'NITa' : {
                            'firex' : 'NO3_ug/m3_DIBB',
                            'gc'    : ['NIT','NITD1','NITD2'],
                            'conv'    : True,
                            'scale'  : 1e9,
                            'mm'  : 62.0049},
    #               'HNO2' : {
    #                        'firex' : 'HNO2_NOAACIMS_VERES',
    #                        'gc'    : 'SpeciesConc_HNO2',
    #                        'conv'    : False,
    #                        'scale'  : 1e12 },
                   'NH3' : {
                            'firex' : 'NH3_UIOPTR_ppbV_WISTHALER',
                            'gc'    : 'SpeciesConc_NH3',
                            'conv'    : False,
                            'scale'  : 1e9 },
                    'HNO3' : {
                            'firex' : 'HNO3-1Hz_CIT_WENNBERG',
                            'gc'    : 'SpeciesConc_HNO3' ,
                            'conv'    : False,
                            'scale'  : 1e12 },
                    'O3' : {
                            'firex' : 'O3_CL_RYERSON',
                            'gc'    : 'SpeciesConc_O3' ,
                            'conv'    : False,
                            'scale'  : 1e9 },
                    'NO' : {
                            'firex' : 'NO_CL_RYERSON',
                            'gc'    : 'SpeciesConc_NO' ,
                            'conv'    : False,
                            'scale'  : 1e9 },
                    'NO2' : {
                            'firex' : 'NO2_CL_RYERSON',
                            'gc'    : 'SpeciesConc_NO2' ,
                            'conv'    : False,
                            'scale'  : 1e9 },
                    'NOx' : {
                            'firex' : 'NOx',
                            'gc'    : 'NOx' ,
                            'conv'    : False,
                            'scale'  : 1e9 },
                    }
    '''
    firex_vars = { 'NOy' : {
                            'firex' : 'NOy_CL_RYERSON',
                            'gc'    : ['NO', 'NO2', 'PAN', 'HNO3',  'PPN', 'R4N2', 'N2O5', 'HNO4', 'BrNO2',
                                       'BrNO3', 'MPN', 'PROPNN', 'NO3', 'HNO2', 'IONO', 'IONO2',
                                       'INO', 'ClNO2', 'ClNO3'] ,
                            'conv'    : False,
                            'scale'  : 1e9 },
                    'CO' : {
                            'firex' : 'CO_DACOM_DISKIN',
                            'gc'    : 'SpeciesConc_CO' ,
                            'conv'    : False,
                            'scale'  : 1e9 }
            }
    '''
    return firex_vars
"""
Functions to interact with GEOS-Chem output after the ARNA campaign
"""
import os
import sys
import glob
import gc
import requests
import re
import wget
import numpy as np
import pandas as pd
import xarray as xr
import xesmf as xe
import AC_tools as AC
import datetime as datetime
import time
from time import gmtime, strftime
# Import from elsewhere in ARNA module
from . core import *
from . utils import *
from . observations import *


def v12_9_TRA_XX_2_name(TRA_XX, folder=None, RTN_dict=False):
    """
    Convert tracer number to tracer name in GEOS-Chem v12.9
    """
    TRAs = AC.get_specieslist_from_input_geos(folder=folder)
    nums = np.arange(1, len(TRAs)+1)
    if RTN_dict:
        return dict(zip(nums, TRAs))
    else:
        return dict(zip(nums, TRAs))[TRA_XX]


def get_dict_of_GEOSChem_model_output(res='0.5x0.625',
                                      RunSet='MERRA2-0.5-initial'):
    """
    Retrieve dictionary of model run names and their full location paths
    """
    RunRoot = get_local_folder('RunRoot')
    if res == '4x5':
        CoreRunStr = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.'
        AcidRunStr = 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.'
        d = {}
        # Boundary condition resolution runs
#        Run = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.PF'
#        folder = '{}/{}/'.format(RunRoot, Run)
#        d = {'BC-BASE': folder}
        # Boundary condition (BC) resolution runs
#        Run = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.PF'
#        folder = '{}/{}/'.format(RunRoot, Run)
#        d = {'BC-BASE-I': folder}
        # Current Boundary condition (BC) resolution runs
        RunStr = 'BCs.repeat/'
        folder = '{}/{}{}/'.format(RunRoot, CoreRunStr, RunStr)
        d['BC-BASE-II'] = folder
        RunStr = 'DustUptake.II'
        folder = '{}/{}{}/'.format(RunRoot, CoreRunStr, RunStr)
        d['AcidUptake-4x5-II'] = folder
        # re-run boundary conditions?
        RunStr = 'DustUptake.BCs'
        folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
        d['AcidUptake-4x5-III'] = folder
        # JNITs
        RunStr = 'DustUptake.JNIT.BCs'
        folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
        d['AcidUptake-4x5-JNIT'] = folder
        # JNITx25
        RunStr = 'DustUptake.JNITx25.BCs'
        folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
        d['AcidUptake-4x5-JNITx25'] = folder
        # Isotherm
#        RunStr = 'DustUptake.JNIT.Isotherm.BCs'
        RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON'
        folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
        d['AcidUptake-4x5-Isotherm'] = folder
        # CVAO campaign in Reed et al 2017
        RunStr = 'geosfp_4x5_standard.v12.9.0.BASE.2015.Aug.ARNA.BCs.repeat'
        folder = '{}/{}/'.format(RunRoot, RunStr)
        d['BC-BASE-2015'] = folder
    elif res == '0.25x0.3125' and (RunSet == 'FP-MOYA-Nest'):
        # GEOS-FP 0.25 nested run
        Run = 'geosfp_4x5_standard.v12.9.0.BASE.2019.2020.MOYA1.Nest'
        folder = '{}/{}/'.format(RunRoot, Run)
        d = {'FP-Nest': folder}
    elif res == '0.25x0.3125':  # and (RunSet=='GEOS-FP-Nest'):
        CoreRunStr = 'geosfp_4x5_standard.v12.9.0.BASE.2019.2020.ARNA.Nest.'
        # GEOS-FP 0.25 nested run
        RunStr = 'repeat'
        folder = '{}/{}{}/'.format(RunRoot, CoreRunStr, RunStr)
        d = {'FP-Nest': folder}
        # GFAS x2
        RunStr = 'repeat.GFASx2'
        folder = '{}/{}{}/'.format(RunRoot, CoreRunStr, RunStr)
        d['FP-Nest-BBx2'] = folder
#       d = {'FP-Nest-BBx2': folder}
        # JNITs x25
        RunStr = 'repeat.JNITs.x25'
        folder = '{}/{}{}/'.format(RunRoot, CoreRunStr, RunStr)
        d['FP-Nest-JNITx25'] = folder
#        d = {'FP-Nest-JNITx25': folder}
        # With Dust active
#        RunStr = 'repeat.IV.JNIT.x25.Dust'
#        folder = '{}/{}{}/'.format(RunRoot, CoreRunStr, RunStr)
#        d['FP-Nest-Dust'] folder}

    elif res == '0.5x0.625':  # and (RunSet=='MERRA2-0.5-initial'):
        Run = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.1x1.PF'
        folder = '{}/{}/'.format(RunRoot, Run)
        d = {'BASE-0.5': folder}
        Run = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.DustUptake.0.5x0.625'
        folder = '{}/{}/'.format(RunRoot, Run)
        d['AcidUptake-0.5'] = folder
        # Add BASE 4x5 run for testing
        Run = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.PF'
        folder = '{}/{}/'.format(RunRoot, Run)
        d['BC-BASE'] = folder
    else:
        pass
    return d


def save_model_output2csv(RunSet='FP-MOYA-Nest', res='0.25x0.3125',
                          folder='./'):
    """
    Save model output as csv file by flight
    """
    import seaborn as sns
    # Which flights to plot?
    if (RunSet == 'FP-MOYA-Nest') and (res == '0.25x0.3125'):
        # Local settings/variables
        flight_IDs = ['C006', 'C007']
        sdate_d = {
            'C006': datetime.datetime(2017, 3, 1),
            'C007': datetime.datetime(2017, 3, 2),
        }
        # Loop by flight and retrieve the files as dataframes
        dfs_mod = {}
        for flight_ID in flight_IDs:
            # Get data
            sdate = sdate_d[flight_ID]
            dfs_mod_GC = get_GEOSChem4flightnum(flight_ID=flight_ID,
                                                res=res,
                                                RunSet=RunSet,
                                                sdate=sdate,
                                                )
            # Save to csv
            df = dfs_mod_GC[list(dfs_mod_GC.keys())[0]]
            filename_str = 'GC_planeflight_data_{}_{}'
            filename = filename_str.format(RunSet, flight_ID)
#            filename = AC.rm_spaces_and_chars_from_str(filename)
            df.to_csv(os.path.join(folder+filename+'.csv'))

    elif (res == '0.25x0.3125') and (RunSet == 'FP-Nest'):
        flight_nums = [
            #    217,
            218, 219, 220, 221, 222, 223, 224, 225,
        ]
        flight_IDs = ['C{}'.format(i) for i in flight_nums]
        # - Loop by flight and retrieve the files as dataframes (mod + obs)
        # Model
        dfs_mod_GC = {}
        for flight_ID in flight_IDs:
            dfs = get_GEOSChem4flightnum(flight_ID=flight_ID, res=res,
                                         RunSet=RunSet,)
            df = dfs[RunSet]
            # Add the derived variables to the dataframe
            df = add_deriv_vars2df(df=df)
#            dfs_mod[flight_ID] = df
            # Save to csv
#            df = dfs_mod_GC[ list(dfs_mod_GC.keys())[0] ]
            filename_str = 'GC_planeflight_data_{}_{}'
            filename = filename_str.format(RunSet, flight_ID)
#            filename = AC.rm_spaces_and_chars_from_str(filename)
            df.to_csv(os.path.join(folder+filename+'.csv'))


def get_GEOSChem4flightnum(flight_ID='C225', res='0.5x0.625', sdate=None,
                           RunSet='MERRA2-0.5-initial', resample_data=True,
                           debug=False):
    """
    Retrieve GEOS-Chem output for FAAM flight
    """
    # Where is the extract GEOS-CF data?
    RunDict = get_dict_of_GEOSChem_model_output(res=res, RunSet=RunSet,)
    # Asume just one run for now...
#    folder = RunDict[ list(RunDict.keys())[0] ]
    dfs = {}
    for Run in list(RunDict.keys()):
        # Extract the data for a specific flight
        folder = RunDict[Run]
        files2use = glob.glob(os.path.join(folder, '*plane.log*'))
        files2use = list(sorted(files2use))
        # Get start date of flight
        # (use the plane-flight from same day as sdate)
        if isinstance(sdate, type(None)):
            dfS = get_summary4flight(flight_ID=flight_ID)
            sdate = dfS.index.values.min()
            edate = dfS.index.values.max()
            sdate, edate = AC.dt64_2_dt([sdate, edate])
        sdate_str = sdate.strftime('%Y%m%d')
        file2use = [i for i in files2use if sdate_str in i]
        AssStr = 'WARNING: more than one ({}) planeflight file found! - {}'
        assert len(file2use) == 1, AssStr.format(len(file2use),file2use)
        file2use = file2use[0]
        # - Instead do this manually for now
        # as cannot as output issue in v12.9 (fixed in runs > initial 4x5)
        try:
            # - Use standard AC_tool extraction
            # Get Header information from first file
            vars, sites = AC.get_pf_headers(file2use, debug=debug)
            # Extract all points from file
            df, vars = AC.pf_csv2pandas(file=file2use, vars=vars, epoch=True,
                                        r_vars=True)
        except ValueError:
            # Open file and save the data into
            with open(file2use, 'rb') as file:
                lines = [i for i in file]
            # Extract as raw data in chunks
            lines_1 = lines[0::2]
            header_1 = lines_1[0].decode('utf-8').split()
            data_1 = [i.decode('utf-8').split() for i in lines_1[1:]]
            df = pd.DataFrame(data_1, columns=header_1)
            lines_2 = lines[1::2]
            header_2 = lines_2[0].decode('utf-8').split()
            data_2 = [i.decode('utf-8').split() for i in lines_2[1:]]
            df2 = pd.DataFrame(data_2, columns=header_2)
            # Now combine
            df = pd.concat([df, df2], axis=1)
            # Now process the meta data/update type formats
            # TODO: below could be faster...
            # Use infer obects? - df.infer_objects
            dtypes = {'POINT': object, 'TYPE': str,
                      'YYYYMMDD': str, 'HHMM': str}
            cols2use = [i for i in df.columns if i not in dtypes.keys()]
            df[cols2use] = df[cols2use].apply(pd.to_numeric)
        # Add a datetime index
        df = AC.DF_YYYYMMDD_HHMM_2_dt(df, rmvars=None, epoch=False)
        df.index.name = None
        # Update the variable names
        d = v12_9_TRA_XX_2_name(None, folder=folder, RTN_dict=True)
        d = dict([('TRA_{:0>3}'.format(i), d[i]) for i in d.keys()])
        df = df.rename(columns=d)
        # Add derived (GEOSchem) variables to df
        df = add_derived_GEOSChem_specs2df(df)
        # Resample the data?
        if resample_data:
            df = df.resample('1T').mean()
        # save df
        dfs[Run] = df.copy()
        del df
    return dfs


def add_derived_GEOSChem_specs2ds(ds, prefix=''):
    """
    Add derived GEOS-Chem variables to xr.Dataset
    """
    # Add NOx as combined NO and NO2
    ds = AC.AddChemicalFamily2Dataset(ds, fam='NOx', prefix=prefix)
    # NOy = no_no2_hno3_hno4_hono_2xn2o5_pan_organicnitrates_aerosolnitrates
    ds = AC.AddChemicalFamily2Dataset(ds, fam='NOy', prefix=prefix)
    ds = AC.AddChemicalFamily2Dataset(ds, fam='NOy-gas', prefix=prefix)
    # Include a variable of NOy where HNO3 is removed
    ds[prefix+'NOy-HNO3'] = ds[prefix+'NOy'] - ds[prefix+'HNO3']
    # Include a variable of NOy where HNO3 is removed
    ds[prefix+'NOy-HNO3-PAN'] = ds[prefix+'NOy'] -  \
                                ds[prefix+'HNO3'] -  \
                                ds[prefix+'PAN']
    # gas-phase (exc. PAN, HNO3, HNO4, Org-NIT, N2O5)
    ds[prefix+'NOy-Limited'] = ds[prefix+'NO'] + \
                               ds[prefix+'NO2'] + \
                               ds[prefix+'HNO2'] + \
                               ds[prefix+'NIT'] + \
                               ds[prefix+'NITs']
    ds = AC.AddChemicalFamily2Dataset(ds, fam='NIT-all', prefix=prefix)
    return ds


def add_derived_GEOSChem_specs2df(df):
    """
    Add derived GEOS-Chem variables to pd.DataFrame
    """
    # Add temperature in deg C
    df['T'] = df['GMAO_TEMP'].copy()
    df['T'] = df['GMAO_TEMP'].values - 273.15
    # Inc. V nd U with same variable names as GEOS-CF
    df['V'] = df['GMAO_VWND'].copy()
    df['U'] = df['GMAO_UWND'].copy()
    # Add NOx as combined NO and NO2
    df['NOx'] = df['NO'].values + df['NO2'].values
    # Add NOy as defined in GEOS-CF
    # NOy = no_no2_hno3_hno4_hono_2xn2o5_pan_organicnitrates_aerosolnitrates
    vars2use = AC.GC_var('NOy-all')
    df['NOy'] = df['N2O5'].copy()  #  2 N2O5 in NOy, so 2x via template
    for var in vars2use:
        try:
            df.loc[:, 'NOy'] = df['NOy'].values + df[var].values
        except KeyError:
            pass
    # Add a variable for gas-phase NOy (by subtracting aerosol nitrate)
    vars2use = AC.GC_var('NOy-gas')
    df['NOy-gas'] = df['N2O5'].copy()  #  2 N2O5 in NOy, so 2x via template
    for var in vars2use:
        try:
            df.loc[:, 'NOy-gas'] = df['NOy-gas'].values + df[var].values
        except KeyError:
            pass
    # Include a variable of NOy where HNO3 is removed
    # NOy = no_no2_hno3_hno4_hono_2xn2o5_pan_organicnitrates_aerosolnitrates
    df['NOy-HNO3'] = df['NOy'].values - df['HNO3'].values
    # Include a variable of NOy where HNO3 is removed
    df['NOy-HNO3-PAN'] = df['NOy'].values - \
        df['HNO3'].values - df['PAN'].values
    # gas-phase (exc. PAN, HNO3, HNO4, Org-NIT, N2O5)
    df['NOy-Limited'] = df['NO'].values + df['NO2'].values + \
        df['HNO2'].values + df['NIT'].values + df['NITs'].values
    # Uset the P-I variable as a model level variable
    df['model-lev'] = df['P-I'].copy()
    return df


def get_whole_related_campaign_data(save2csv=True, campaign='ARNA-1',
                                    resample_data=False, debug=False):
    """
    Get model data from additional CVAO campaigns (surface+airborne)
    """
    # Which data to use?
    RunDir = '/users/ts551/scratch/GC/rundirs/'
    BASE_str = 'geosfp_4x5_standard.v12.9.0.BASE'
    ARNA1Var = BASE_str+'.2019.2020.ARNA1.Nest.repeat.JVALS/'
    ARNA1SurVar = BASE_str+'.2019.2020.ARNA1.Nest.repeat.JVALS.CVAO.PF/'
    CVAO2015 = BASE_str+'.2015.Aug.Nest.repeat.JVALS.CVAO.PF/'
    run_dict = {
    'ARNA-1' : RunDir+ARNA1Var,
    'ARNA-1-surface' : RunDir+ARNA1SurVar,
    'CVAO-2015-surface' : RunDir+CVAO2015,
    }
    # Name of file to save?
    SaveName = 'GC_model_output_{}'.format(campaign)
    SaveName = AC.rm_spaces_and_chars_from_str(SaveName)
    # NetCDF directory
    sanity_check_model_runs = False
    if sanity_check_model_runs:
        for key in run_dict.keys():
            run_dict[key] = run_dict[key]+'/OutputDir/'
        # check generic stats
        RunStr = '/geosfp_4x5_standard.v12.9.0.BASE.2019.2020.ARNA.BCs.repeat/'
        RunStr += 'OutputDir/'
        REF_wd = RunDir + RunStr
        use_REF_wd4Met = True
        df = AC.get_general_stats4run_dict_as_df(run_dict=run_dict,
                                                 use_REF_wd4Met=use_REF_wd4Met,
                                                 REF_wd=REF_wd,
                                                 )
    # Extract the planeflight files
    folder = run_dict[campaign]
    files2use = list(sorted(glob.glob(os.path.join(folder, '*plane.log*'))))
    file2use = files2use[0]
    # Get Header information from first file
    vars, sites = AC.get_pf_headers(file2use, debug=debug)
    # Extract all points from file
    dfs = []
    for file2use in files2use:
        df, vars = AC.pf_csv2pandas(file=file2use, vars=vars, epoch=True,
                                    r_vars=True)

        # Add a datetime index
        df = AC.DF_YYYYMMDD_HHMM_2_dt(df, rmvars=None, epoch=False)
        df.index.name = None
        # Set the out of (nested) box values to NaNs
        OutOfBoxValue =  -1000.0
        for col in df.columns:
            try:
                df.loc[ df[col].values == OutOfBoxValue, col] = np.NaN
            except TypeError:
                print('{} - TypeError found for {} '.format(campaign, col))
        # Update the variable names
        d = v12_9_TRA_XX_2_name(None, folder=folder, RTN_dict=True)
        d = dict([('TRA_{:0>3}'.format(i), d[i]) for i in d.keys()])
        df = df.rename(columns=d)
        # Add derived (GEOSchem) variables to df
        df = add_derived_GEOSChem_specs2df(df)
        # Resample the data?
        if ('surface' not in campaign) or resample_data:
            df = df.resample('1T').mean()
        dfs += [df]
    # Concat the data frames
    df2 = pd.concat(dfs, axis=0)
    # Save to disk
    if save2csv:
        if ('surface' not in campaign):
            df2.to_csv(SaveName+'.csv')
        else:
            TYPES = list(set(df2['TYPE'].values))
            print(TYPES)
            for TYPE in TYPES:
                df2save = df2.loc[ df2['TYPE'] == TYPE, : ]
                df2save.to_csv('{}_{}.csv'.format(SaveName, TYPE))
    return df2


def regrid_restart4ARNA_highres_grid(folder=None, filename=None, res='1x1'):
    """
    Regrid a restart file for the ARNA grid
    """
    # File and location?
    if isinstance(folder, type(None)):
        RunRoot = get_local_folder('RunRoot')
        RunDir = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.1x1.PF/'
        folder = RunRoot + RunDir
    if isinstance(filename, type(None)):
        filename = 'GEOSChem.Restart.20200201_0000z.nc4.FILE2REGRID'
    ds = xr.open_dataset(folder+filename)
    # Get the lons and lats for a requested resolution
    lons, lats, alt = AC.get_latlonalt4res(res=res)
#    lons = np.arange(-180, 180, 1)
#    lats = np.arange(-90, 90, 1)
    lons = np.arange(-180, 180, 0.3125)
    lats = np.arange(-90, 90, 0.25)

#    OutFile = 'GEOSChem.Restart.20200120_0000z.REGRIDED.nc4'
    OutFile = '{}{}.nc4'.format(filename.split('nc4')[0], 'REGRIDED')
    ds = AC.regrid_restart_file4flexgrid(ds, OutFile=OutFile, folder=folder,
                                         lons=lons, lats=lats)
    # Now chop out the region for the restart file
    LonMax = 15.
    LonMin = -35.
    LatMin = 0.
    LatMax = 34.
    bool1 = ((ds.lon >= LonMin) & (ds.lon <= LonMax)).values
    bool2 = ((ds.lat >= LatMin) & (ds.lat <= LatMax)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    # Save the re-gridded file
    OutFile = '{}{}.nc4'.format(filename.split('nc4')[0], 'REGRIDED_CROPPED')
    ds.to_netcdf(os.path.join(folder, OutFile))

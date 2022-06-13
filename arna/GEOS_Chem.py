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


def PF_TRAXXX_2TracerName(TRA_XX, folder=None, RTN_dict=False):
    """
    Convert tracer number to tracer name in GEOS-Chem v12.9
    """
    TRAs = AC.get_specieslist_from_input_geos(folder=folder)
    nums = np.arange(1, len(TRAs)+1)
    if RTN_dict:
        return dict(zip(nums, TRAs))
    else:
        return dict(zip(nums, TRAs))[TRA_XX]


def get_dict_of_GEOSChem_model_output(res='0.5x0.625', folder4netCDF=False,
                                      RunSet='MERRA2-0.5-initial',
                                      GC_version='v12.9',
                                      CoreRunsOnly=False):
    """
    Retrieve dictionary of model run names and their full location paths

    Parameters
    -------
    folder4netCDF (bool): append subfolder for NetCDFs (OutputDir) to folder

    Notes
    -----
    """
    # Get the root directory for model run directories
    RunRoot = get_local_folder('RunRoot')
    # Which GEOS-Chem version is being used?
#     if 'v12.9' in GC_version:
#         Is_v12_9 = True
#     elif 'v13.4' in GC_version:
#         Is_v12_9 = False
#     else:
#         print("No runs in funciton for GC_version ('{}')".format(GC_version))
    if (res == '4x5') and (RunSet == 'PostHONO'):
        RunDict = get_RunDict_of_HONO_runs(RunSet=RunSet,
                                           GC_version=GC_version)
        #
        if CoreRunsOnly:
            runs2use = ['Base', 'NIThv', 'NIThvCap']
            dNew = {}
            for run in runs2use:
                dNew[run] = RunDict[run]
            RunDict = dNew

    elif (res == '4x5') and (RunSet == 'MERRA2-0.5-initial'):
        CoreRunStrMERRA2 = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.'
        CoreRunStrGEOSFP = 'geosfp_4x5_standard.v12.9.0.BASE.2019.2020.'
        AcidRunStr = 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.'
        RunDict = {}
        # JNIT and extra acid runs
        if (RunSet == 'ACID'):
            # re-run boundary conditions?
            #             RunStr = 'DustUptake.BCs'
            #             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
            #             RunDict['Acid-4x5'] = folder
            #            RunDict['Acid-4x5-III'] = folder
            #             RunStr = 'DustUptake.JNIT.BCs'
            #             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
            #             RunDict['Acid-4x5-JNIT'] = folder
            # ACID + J50
            RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags.v2.J00'
            folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
            RunDict['Acid-4x5-J00'] = folder
            # JNITx25
#             RunStr = 'DustUptake.JNITx25.BCs'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-J25'] = folder
            # Isotherm
    #        RunStr = 'DustUptake.JNIT.Isotherm.BCs'
#            RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II'
#             RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-Isotherm'] = folder
            # Isotherm
    #        RunStr = 'DustUptake.JNIT.Isotherm.BCs'
#            RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II'
            IsoRunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags.'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, IsoRunStr+'v2')
#             RunDict['Acid-4x5-Isotherm.v2.0'] = folder
            # ... Other Isotherm versions
#             RunStr = IsoRunStr + 'v2.very_low'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-Isotherm.v2.3'] = folder
#             RunStr = IsoRunStr + 'v2.low'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-Isotherm.v2.2'] = folder
#             RunStr = IsoRunStr + 'v2.medium'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-Isotherm.v2.1'] = folder
#             RunStr = IsoRunStr + 'v2.deli'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-Isotherm.v2.0.deli'] = folder
#             RunStr = IsoRunStr + 'v2.4'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-Isotherm.v2.4'] = folder
#             RunStr = IsoRunStr + 'v2.4.deli'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-Isotherm.v2.4.deli'] = folder
#             RunStr = IsoRunStr + 'v2.4.deli.H2O'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-Isotherm.v2.4.deli.H2O'] = folder
            # Version 3
            RunStr = IsoRunStr + 'v3.0.H2O'
            folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
            RunDict['Acid-4x5-Isotherm.v3.0.H2O'] = folder
            # Version 3 plus v3.0
            RunStr = IsoRunStr + 'v3.0.H2O.Acid'
            folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
            RunDict['Acid-4x5-Isotherm.v3.0.H2O.Acid'] = folder
            # Re-run acid exploration?

            # Isotherm + HONO 100%
#             RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags.HONO100'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-Isotherm-HONO100'] = folder
            # Isotherm + BBx3
#             RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags.BBx3'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-Isotherm-BBx3'] = folder
            # ACID + J50
            RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags.J50'
            folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
            RunDict['Acid-4x5-J50'] = folder
            # ACID + J50 (no BB)
#             RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags.J50.BBx0'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-J50-BBx0'] = folder
            # Acid plus v3.0

            # Isotherm + BBx3 + NH3x3
#             RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags'
#             RunStr += '.J50.BBx3.NH3x3'
#             folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
#             RunDict['Acid-4x5-J50-BBx3-NH3x3'] = folder
            # Isotherm + African BBx3 + African NH3x3
            RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags'
            RunStr += '.J50.BBx3AFRICA.NH3x3/'
            folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
            RunDict['Acid-4x5-J50-AfBBx3-NH3x3'] = folder
            # Temporary runs - MERRA-2 met
#             RunStr = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.PF'
#             folder = '{}/{}{}/'.format(RunRoot, RunStr, '')
#             RunDict['BASE-4x5-MERRA-2'] = folder

            # Isotherm, deliquescent limit, Cap Jscale 50, HONO channel 100%
            RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags'
            RunStr += '.v3.0.H2O.cap2J50.HONO100'
            folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
            RunDict['Acid-4x5-Isotherm.v3aqCap50H100'] = folder

            # 4 pptv HONO everywhere
            RunStr = 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags'
            RunStr += '.v2.J00.HourlyOutput.HONO4pptv'
            folder = '{}/{}{}/'.format(RunRoot, AcidRunStr, RunStr)
            RunDict['Acid-4x5-4pptHONO'] = folder

            # Only return the core runs
            if CoreRunsOnly:
                runs2use = [
                    #                    'Acid-4x5',
                    #                    'Acid-4x5-J25',
                    'Acid-4x5-J00',
                    'Acid-4x5-J50',
                    #                    'Acid-4x5-Isotherm.v2.4',
                    #                    'Acid-4x5-J50-AfBBx3-NH3x3',
                    #                    'Acid-4x5-Isotherm.v2.4',
                    #                    'Acid-4x5-Isotherm.v2',
                    #                    'Acid-4x5-Isotherm.v3.0.H2O',
                    #                    'Acid-4x5-Isotherm.v3.0.H2O.Acid',
                    # core runs as of March 2022
                    'Acid-4x5-Isotherm.v3aqCap50H100',
                    #                    'Acid-4x5-4pptHONO',
                    # Temp
                    #                    'BASE-4x5-MERRA-2',
                    #                    'Acid-4x5-J50-BBx0',
                ]
                dNew = {}
                for run in runs2use:
                    dNew[run] = RunDict[run]
                RunDict = dNew

        else:  # consider the variants on base runs (e.g. Biomass burning)
            # Boundary condition resolution runs
            #        Run = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.PF'
            #        folder = '{}/{}/'.format(RunRoot, Run)
            #        RunDict = {'BC-BASE': folder}
            # Boundary condition (BC) resolution runs
            #        Run = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.PF'
            #        folder = '{}/{}/'.format(RunRoot, Run)
            #        RunDict = {'BC-BASE-I': folder}

            # BASE runs, but with increases in JNITs/BB
            #         RunStr = 'BCs.repeat.JNITx100.GFASx3/'
            #         folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
            #         RunDict['BC-BASE-BBx3-J100'] = folder
            #         RunStr = 'BCs.repeat.JNITx50.GFASx3/'
            #         folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
            #         RunDict['BC-BASE-BBx3-J50'] = folder
            #         RunStr = 'BCs.repeat.JNITx25.GFASx3/'
            #         folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
            #         RunDict['BC-BASE-BBx3-J25'] = folder
            # Using the dust uptake setup
            #         RunStr = 'DustUptake.II'
            #         folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
            #         RunDict['Acid-4x5-II'] = folder
            # Current Boundary condition (BC) resolution runs
            RunStr = 'ARNA.BCs.repeat'
            folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
            RunDict['BC-BASE'] = folder
            # J25
#             RunStr = 'ARNA.BCs.repeat.JNITx25'
#             folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
#             RunDict['BC-BASE-J25'] = folder
            # J50
            RunStr = 'ARNA.BCs.repeat.JNITx50'
            folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
            RunDict['BC-BASE-J50'] = folder
            # J50 + HONO-100
            RunStr = 'ARNA.BCs.repeat.JNITx50.HONO100'
            folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
            RunDict['BC-BASE-J50-HONO100'] = folder
            # BBx2
#             RunStr = 'ARNA.GFASx2.BCs'
#             folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
#             RunDict['BC-BASE-BBx2'] = folder
            # BBx3
            RunStr = 'ARNA.GFASx3.BCs'
            folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
            RunDict['BC-BASE-BBx3'] = folder
            # BBx3 + J100
            RunStr = 'ARNA.BCs.repeat.JNITx100.GFASx3'
            folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
            RunDict['BC-BASE-BBx3-J100'] = folder
            # BBx3 + J25
#             RunStr = 'ARNA.BCs.repeat.JNITx25.GFASx3'
#             folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
#             RunDict['BC-BASE-BBx3-J25'] = folder
            # BBx3 + J50
#             RunStr = 'ARNA.BCs.repeat.JNITx50.GFASx3'
#             folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
#             RunDict['BC-BASE-BBx3-J50'] = folder
            # BBx4
#             RunStr = 'ARNA.GFASx4.BCs'
#             folder = '{}/{}{}/'.format(RunRoot, CoreRunStrGEOSFP, RunStr)
#             RunDict['BC-BASE-BBx4'] = folder
            # Only return the core runs
            if CoreRunsOnly:
                runs2use = ['BC-BASE', 'BC-BASE-BBx3', 'BC-BASE-J50']
                dNew = {}
                for run in runs2use:
                    dNew[run] = RunDict[run]
                RunDict = dNew
        # Extra runs...
        # CVAO campaign during 2015 (Reed et al 2017)
#         RunStr = 'geosfp_4x5_standard.v12.9.0.BASE.2015.Aug.ARNA.BCs.repeat'
#         folder = '{}/{}/'.format(RunRoot, RunStr)
#         RunDict['BC-BASE-2015'] = folder
    elif (res == '0.25x0.3125') and (RunSet == 'FP-MOYA-Nest'):
        # GEOS-FP 0.25 nested run
        Run = 'geosfp_4x5_standard.v12.9.0.BASE.2019.2020.MOYA1.Nest'
        folder = '{}/{}/'.format(RunRoot, Run)
        RunDict = {'FP-Nest': folder}
    elif (res == '0.25x0.3125'):  # and (RunSet=='GEOS-FP-Nest'):
        CoreRunStr = 'geosfp_4x5_standard.v12.9.0.BASE.2019.2020.ARNA.Nest.'
        # GEOS-FP 0.25 nested run
        RunStr = 'repeat'
        folder = '{}/{}{}/'.format(RunRoot, CoreRunStr, RunStr)
        RunDict = {'FP-Nest': folder}
        # GFAS x2
        RunStr = 'repeat.GFASx2'
        folder = '{}/{}{}/'.format(RunRoot, CoreRunStr, RunStr)
        RunDict['FP-Nest-BBx2'] = folder
#       RunDict = {'FP-Nest-BBx2': folder}
        # JNITs x25
        RunStr = 'repeat.JNITs.x25'
        folder = '{}/{}{}/'.format(RunRoot, CoreRunStr, RunStr)
        RunDict['FP-Nest-JNITx25'] = folder
#        RunDict = {'FP-Nest-JNITx25': folder}
        # With Dust active
#        RunStr = 'repeat.IV.JNIT.x25.Dust'
#        folder = '{}/{}{}/'.format(RunRoot, CoreRunStr, RunStr)
#        RunDict['FP-Nest-Dust'] folder}

    elif (res == '0.5x0.625'):  # and (RunSet=='MERRA2-0.5-initial'):
        Run = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.1x1.PF'
        folder = '{}/{}/'.format(RunRoot, Run)
        RunDict = {'BASE-0.5': folder}
        Run = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.DustUptake.0.5x0.625'
        folder = '{}/{}/'.format(RunRoot, Run)
        RunDict['Acid-0.5'] = folder
        # Add BASE 4x5 run for testing
        Run = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.PF'
        folder = '{}/{}/'.format(RunRoot, Run)
        RunDict['BC-BASE'] = folder
    else:
        pass
    # Include NetCDF subfolder in run directory strings
    if folder4netCDF:
        for key in RunDict.keys():
            RunDict[key] = RunDict[key] + '/OutputDir/'
    return RunDict


def get_RunDict_of_HONO_runs(RunSet='PostHONO', GC_version='v13.4',
                             folder4netCDF=False):
    """
    Retrieve the model runs investigating HONO
    """
    # Which runs to use
    RunRoot = get_local_folder('RunRoot')
    if (RunSet == 'PostHONO') and (GC_version == 'v13.4'):
        RunStr = 'gc_4x5_47L_geosfp_fullchem.v13.4.0-rc.2'
        RunDict = {
            'Base': '{}{}{}/'.format(RunRoot, RunStr, '.orig.1monthTest'),
            'min4pptHONO': '{}{}{}/'.format(RunRoot, RunStr,
                                             '.orig.ARNA.HONO.4pptv'),
            #        '4pptHONO': '{}{}{}/'.format(RunRoot, RunStr,
            #                                     '.orig.ARNA.HONO.4pptv.all'),
            #        'NOxSink': '{}{}{}/'.format(RunRoot, RunStr, ''),
            #    'HalNitratesx10': '{}{}{}/'.format(RunRoot, RunStr, ''),
            #        'NIThv': '{}{}{}/'.format(RunRoot, RunStr, ''),
            #        'OH+NO2': '{}{}{}/'.format(RunRoot, RunStr, ''),
            #         'HO2+NO':  '{}{}{}/'.format(RunRoot, RunStr, ''),
            #         'N2O5': '{}{}{}/'.format(RunRoot, RunStr, ''),
        }
    elif (RunSet == 'PostHONO') and (GC_version == 'v12.9'):
        RunStr = 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.DustUptake'
        RunStr += '.JNIT.Isotherm.BCs.repeat.ON.II.diags'
        RunStr2 = RunStr + '.v2.J00.HourlyOutput'
        RunStr3 = RunStr + '.v5.0.H2O.AcidII.100HONO.2018'

        RunDict = {
        'Base': '{}{}{}/'.format(RunRoot, RunStr,
                                 '.v2.J00.HourlyOutput.2018'),
        'min4pptHONO': '{}{}{}/'.format(RunRoot, RunStr2,
                                         '.HONO4pptv'),
        'min4pptHONOday': '{}{}{}/'.format(RunRoot, RunStr2,
                                         '.HONO4pptv.night'),
        '4pptHONO': '{}{}{}/'.format(RunRoot, RunStr2,
                                     '.HONO4pptv.all'),
        'HalNitratesx10': '{}{}{}/'.format(RunRoot, RunStr2,
                                           '.2018.HalNitrate'),
        'NIThv': '{}{}{}/'.format(RunRoot, RunStr,
                                    '.v4.0.H2O.AcidII.100HONO.2018.pH7'),
        'NIThvCap': '{}{}{}/'.format(RunRoot, RunStr,
                                    '.v4.0.H2O.pH7.100HONO.Cap50'),
        'NIThvpH': '{}{}{}/'.format(RunRoot, RunStr,
                                   '.v4.0.H2O.pH4.100HONO.Cap50'),
        'N2O5': '{}{}{}/'.format(RunRoot, RunStr2, '.2018.N2O5.repeat'),
        # The below runs do not yet have enough output to analysis
        # But should do before 9am 17th March
        'NOxSink': '{}{}{}/'.format(RunRoot, RunStr2,
                                    '.2018.NOxSinkv1.1'),
        'HO2+NOv1.Cap':  '{}{}{}/'.format(RunRoot, RunStr2,
                                     '.2018.HO2andNO.repeat'
                                     ),
        'HO2+NOv1':  '{}{}{}/'.format(RunRoot, RunStr2,
                                     '.2018.HO2andNO.repeat.NoCap'
                                     ),
        'HO2+NOv2noH2O':  '{}{}{}/'.format(RunRoot, RunStr2,
                                     '.2018.HO2andNOv2'
                                     ),
        'HO2+NOv2.Cap':  '{}{}{}/'.format(RunRoot, RunStr2,
                                     '.2018.HO2andNOv2.H2O.Capped'
                                     ),
        'HO2+NOv2':  '{}{}{}/'.format(RunRoot, RunStr2,
                                     '.2018.HO2andNOv2.H2O'
                                     ),

        'NIThv.HO2+NOv2.Cap':  '{}{}{}/'.format(RunRoot, RunStr3,
                                     '.pH7.HO2andNOv2.H2O.Capped'
                                     ),

        'OH+NO2': '{}{}{}/'.format(RunRoot, RunStr2, '.2018.NO2andOH'),
        # The below probably will not have enough useful data by then
        #        'Amedro2020': '{}{}{}/'.format(RunRoot, RunStr2,
        #                                      '.2018.NO2andOH.Amedro'),
        # Redundent
#         'NIThv': '{}{}{}/'.format(RunRoot, RunStr,
#                                   '.v3.0.H2O.AcidII.100HONO.2018/'),
#         'NIThvCap': '{}{}{}/'.format(RunRoot, RunStr,
#                                      '.v3.0.H2O.cap2J50.HONO100.2018'),
        }
    else:
        PrtStr = "Unkonwn RunSet ('{}') and GC verions ('{}')"
        print(PrtStr.format(RunSet, GC_version))
        sys.exit(0)
    # Include the output directory folder in directory strings
    if folder4netCDF:
        for key in RunDict.keys():
            RunDict[key] = '{}{}/'.format(RunDict[key], 'OutputDir')
    return RunDict


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
                           CoreRunsOnly=False,
                           debug=False):
    """
    Retrieve GEOS-Chem output for FAAM flight
    """
    # Where is the extract GEOS-CF data?
    RunDict = get_dict_of_GEOSChem_model_output(res=res, RunSet=RunSet,
                                                CoreRunsOnly=CoreRunsOnly)
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
        AssStr = 'WARNING: More than one ({}) planeflight file found! - {}'
        assert len(file2use) <= 1, AssStr.format(len(file2use), file2use)
        AssStr = 'WARNING: No planeflight files found in folder: {}'
        assert len(file2use) != 0, AssStr.format(folder)
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
        d = PF_TRAXXX_2TracerName(None, folder=folder, RTN_dict=True)
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
    # Add an all sulfate tracer
    NewVar = 'SO4-all'
    vars2use = AC.GC_var(NewVar)
    df[NewVar] = df['NIT'].values
    for var2use in vars2use:
        try:
            df[NewVar] = df[NewVar].values + df[var2use].values
        except KeyError:
            pass
    # And a all nitrate tracer
    NewVar = 'NIT-all'
    vars2use = AC.GC_var(NewVar)
    df[NewVar] = df['NIT'].values
    for var2use in vars2use:
        try:
            df[NewVar] = df[NewVar].values + df[var2use].values
        except KeyError:
            pass
    # Uset the P-I variable as a model level variable
    df['model-lev'] = df['P-I'].copy()
    return df


def get_whole_related_campaign_data(save2csv=True, campaign='ARNA-1',
                                    resample_data=False, debug=False):
    """
    Get model data from additional CVAO campaigns (surface+airborne)
    """
    # Which data to use?
    RunRoot = get_local_folder('RunRoot')
    BASE_str = 'geosfp_4x5_standard.v12.9.0.BASE'
    ARNA1Var = BASE_str+'.2019.2020.ARNA1.Nest.repeat.JVALS/'
    ARNA1SurVar = BASE_str+'.2019.2020.ARNA1.Nest.repeat.JVALS.CVAO.PF/'
    CVAO2015 = BASE_str+'.2015.Aug.Nest.repeat.JVALS.CVAO.PF/'
    run_dict = {
        'ARNA-1': '{}{}'.format(RunRoot, ARNA1Var),
        'ARNA-1-surface': '{}{}'.format(RunRoot, ARNA1SurVar),
        'CVAO-2015-surface': '{}{}'.format(RunRoot, CVAO2015),
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
        REF_wd = '{}{}'.format(RunRoot, RunStr)
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
        OutOfBoxValue = -1000.0
        for col in df.columns:
            try:
                df.loc[df[col].values == OutOfBoxValue, col] = np.NaN
            except TypeError:
                print('{} - TypeError found for {} '.format(campaign, col))
        # Update the variable names
        d = PF_TRAXXX_2TracerName(None, folder=folder, RTN_dict=True)
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
                df2save = df2.loc[df2['TYPE'] == TYPE, :]
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


def mkNetCDF_mask4CVAO_African_nest(res='2x2.5', invert_mask=False,
                                    regrid21x1=True):
    """
    Create a mask of sub-equatorial africa as a NetCDF to be read by HEMCO
    """
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    import xesmf as xe
    ds = AC.get_LWI_map(res=res, rtn_ds=True)
    ds = ds.transpose('lon', 'lat')
    var2use = 'LWI'
    # Get mask
    m = AC.get_CVAO_Africa_nest_Masked(res=res)
    # Name to save as
    if regrid21x1:
        res = '1x1'
    savename = 'ARNA_CVAO_African_nest_{}'.format(res)
    vals = ds[var2use].values
    # Change mask values
    if invert_mask:
        vals = np.where(m, vals, 0)
        vals = np.where(~m, vals, 1)
        ds[var2use].values = vals
        savename += '_inverted'
        ExtraStr = ' inverted'
    else:
        vals = np.where(~m, vals, 0)
        vals = np.where(m, vals, 1)
        ds[var2use].values = vals
        ExtraStr = ''
    # Only consider variable of interest
    ds = ds[[var2use]]

    # Plot up without mask present
    AC.quick_map_plot(ds, var2plot=var2use, savename=savename)
    plt.title(savename)

    # --- Add time to ds
    # Include a time dimension to meet COARDS
    try:
        ds = ds.expand_dims('time', axis=0)
    except ValueError:
        pass
    # Update Namings
    Var2SaveAs = 'MASK'
    ds = ds.rename({var2use: Var2SaveAs})
    var2use = Var2SaveAs
    # make sure ordering meets COARDS
    ds = ds.transpose('time', 'lat', 'lon')
    # --- Regrid
    if regrid21x1:
        # Regrid to use the same grid as the EMEP data (if the 2x2.5 fails)
        Fstr = get_local_folder('HEMCO_data')
        Fstr += '/MASKS/v2018-09/EMEP_mask.geos.1x1.20151222.nc'
        dsREF = xr.open_dataset(Fstr)
        # Create a dataset to re-grid into
        ds_out = xr.Dataset({
            'lat': (['lat'], dsREF['lat']),
            'lon': (['lon'], dsREF['lon']),
        })
        # Create a regidder (to be reused )
        regridder = xe.Regridder(ds, ds_out, 'bilinear', reuse_weights=True)
        # Loop and regrid variables
        ds_l = []
        for var_ in [var2use]:
            # Create a dataset to re-grid into
            ds_out = xr.Dataset({
                'lat': (['lat'], dsREF['lat']),
                'lon': (['lon'], dsREF['lon']),
            })
            # Get a DataArray
            dr = ds[var_]  # .to_array()
            # build regridder
            dr_out = regridder(dr)
            # Important note: Extra dimensions must be on the left, i.e. (time, lev, lat, lon) is correct but (lat, lon, time, lev) would not work. Most data sets should have (lat, lon) on the right (being the fastest changing dimension in the memory). If not, use DataArray.transpose or numpy.transpose to preprocess the data.
            # exactly the same as input
#            xr.testing.assert_identical(dr_out['time'], dr['time'])
            # Save variable
            ds_l += [dr_out]
            # Combine variables
            dsN = xr.Dataset()
            for n, var2use in enumerate([var2use]):
                dsN[var2use] = ds_l[n]
                dsN[var2use].attrs = ds[var2use].attrs
            # Clean up
            regridder.clean_weight_file()
        # Parse the lat and lon attrs from original array
        dsN['lat'].attrs = ds['lat'].attrs
        dsN['lon'].attrs = ds['lon'].attrs
        # Now use re-gridded file
        del ds
        ds = dsN

    # - Plot up
    # Select data to plot
    a = ds[var2use].mean(dim='time')
    # Plot up map with mask present
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection=ccrs.Robinson(), aspect='auto')
    a.plot.imshow(x='lon', y='lat', ax=ax, transform=ccrs.PlateCarree())
    plt.title(savename)
    plt.savefig(savename+'.png')
    plt.close()

    # - Update attrs etc
    try:
        ds = ds.expand_dims('time', axis=0)
    except ValueError:
        pass
    ds = ds.assign_coords(time=[0])
    attrs = ds['time'].attrs
    attrs['standard_name'] = 'Time'
    attrs['long_name'] = attrs['standard_name']
    attrs['axis'] = 'T'
    attrs['units'] = 'hours since 2000-01-01 00:00:00'
    attrs['calendar'] = 'standard'
    ds['time'].attrs = attrs
    # Save array as float32 instead of float 64
    ds[Var2SaveAs].values = ds[Var2SaveAs].astype(np.float32).values
    # Update MASK attributes
    MainVarTitle = 'MASK of African high-res model nest '+ExtraStr
    attrs = ds[var2use].attrs
    attrs['standard_name'] = Var2SaveAs
    attrs['long_name'] = MainVarTitle
#    attrs['units'] = int(1)
    attrs['units'] = 'unitless'
    attrs['add_offset'] = int(0)
    attrs['scale_factor'] = int(1)
    attrs['missing_value'] = float(-1e-32)
    attrs['_FillValue'] = float(-1e-32)
    # Remove the vestigial ctm_units attribute
    try:
        del attrs['ctm_units']
    except:
        pass
    ds[var2use].attrs = attrs
    # Global attributes
    attrs = ds.attrs
    attrs['conventions'] = 'COARDS'
    attrs['contact'] = 'tomas.sherwen@york.ac.uk'
    attrs['title'] = MainVarTitle
    attrs['history'] = 'Created {}'.format(time.ctime(time.time()))
#    attrs['format'] = 'NetCDF-4'
    attrs['format'] = 'NETCDF3_CLASSIC'
    ds.attrs = attrs
    # Save as NetCDF
#    ds.to_netcdf( savename+'.nc', engine='scipy' )
    ds.to_netcdf(savename+'.nc', format='NETCDF3_CLASSIC')


def tag_GC_simulations():
    """
    Setup tagging of NOx/HNO2 prod etc in KPP
    """
    # Diagnostics to use?
    d = get_tags_for_NOx_HONO()
    tags = d.values()
    for key in d.keys():
        print('{} : {};'.format(key, d[key]))
    # Also print out just using "P" as the prefix.
    for key in d.keys():
        print('P{} : {};'.format(d[key], d[key]))
    # prepare other output for GEOS-Chem input files
    extr_str = 'ARNA_Standard'
    AC.print_out_lines_for_gckpp_file(tags=tags, extr_str=extr_str)
    AC.prt_lines4species_database_yml(tags=tags, extr_str=extr_str)

    ptr_str = '{:<11}= IGNORE; {}'
    d = dict(zip(diags, tags))
    for key in d.keys():
        print(ptr_str.format(d[key], '{'+key+'}'))

    ptr_str = 'P{:<10}= IGNORE; {}'
    d = dict(zip(diags, tags))
    for key in d.keys():
        print(ptr_str.format(d[key], '{'+key+'}'))


def get_tags_for_NOx_HONO(AllTags=False):
    """
    Function to store tags for NOx/HONO
    """
    diags = [
        # Version 6 tags
        'ProdHNO2fromHvNIT', 'ProdHNO2fromHvNITs', 'ProdHNO2fromHvNITD1',
        'ProdHNO2fromHvNITD2', 'ProdHNO2fromHvNITD3', 'ProdHNO2fromHvNITD4',
        'ProdNO2fromHvNIT', 'ProdNO2fromHvNITs', 'ProdNO2fromHvNITD1',
        'ProdNO2fromHvNITD2', 'ProdNO2fromHvNITD3', 'ProdNO2fromHvNITD4',
        'ProdNO2fromHONO', 'ProdHNO2fromOHandNO', 'ProdHNO2fromHET',
        'ProdNOnHO2ChannelA', 'ProdNOnHO2ChannelB',
        # Version 7 tags
        'ProdHNO3fromNO2nOH','ProdNO3fromHNO3nOH',
        'PhotNO2', 'PhotHNO3', 'PhotHNO2',
        'ProdHNO3fromHetNO3', 'ProdNITfromHetNO3','ProdNITsfromHetNO3',
    ]
    prefix = 'TN{:0>3}'
    tags = [prefix.format(i+1) for i in range(len(diags))]
    # pair up numbering (so that runs with different diagnostics have same #s)?
    d = dict(zip(diags, tags))
    # Include the automatic tagging of NOx
    def mk_KPP_tag_from_rxn_str(rxn_str=None, search_str=None,
                                prefix='ProdfromRXN', ):
        """
        Create a variable for reaction
        """
        reactants = rxn_str.split('=')[0]
        reactants = reactants.replace(' + ', '_n_')
        reactants = reactants.replace(' {+M} ', '_M_').strip()
        products = rxn_str.split('=')[-1]
        products = products.replace(' + ', '_n_')
        products = products.replace(' {+M} ', '_M_').strip()
        products = products.replace(' {+M}', '_M').strip()
        products = products[:10]
        # Return a new reaction string
        return'{}_{}_{}_to_{}'.format(prefix, search_str, reactants, products)

    if AllTags:
        DataRoot = get_local_folder('DataRoot')
        folder = '{}{}'.format(DataRoot, '/ARNA/Misc/')
#        FName = 'Tagged_reactions_in_Standard_v12.9.1_ARNA_v8_POx_tagged.csv'
        FName = 'Tagged_reactions_in_Standard_v12.9_ARNA_v9_PL_NOx_tagged.csv'
        df = pd.read_csv(folder+FName)
#        df['RxnName'] = df['rxn_str'].map(mk_KPP_tag_from_rxn_str)
        df['RxnName'] = df.apply(lambda x:
                                 mk_KPP_tag_from_rxn_str(rxn_str=x['rxn_str'],
                                 search_str = x['search_str'], ),
                                 axis=1)

        # combine into main dictionary
        d2 = dict(zip( df['RxnName'], df['tag'].values ) )
        d = AC.merge_two_dicts(d, d2)
    return d


def get_DryDepAndWetDep_ds(wd=None, Specs=None, dates2use=None,
                             DDprefix='DryDep_',
                             DDVprefix='DryDepVel_',
                             WLprefix='WetLossConv_',
                             LSprefix='WetLossLS_',
                             NewDepPrefix = 'DryAndWetDep_',
                             verbose=False, debug=False,
                             ):
    """
    Get a combined wet and dry loss value for species in kg/m2/s
    """
    # Avogadro constant (Mol^-1)
    AVG = AC.constants('AVG')
    # Species to process to a combined depositional loss
    if isinstance(Specs, type(None)):
        Specs = ['NO2', 'HNO3', 'NIT', 'NITs', 'NIT-all' ]
        Specs += ['NITD{}'.format(i) for i in np.arange(1, 5)]

    # Get dry Deposition
    try:
        Dds = AC.get_DryDep_ds(wd=wd, dates2use=dates2use)
        # Add all all-NIT to dataset
        if any([('NIT' in i) for i in Specs]):
            Dds = AC.AddChemicalFamily2Dataset(Dds, fam='NIT-all',
                                                  prefix=DDprefix)
        No_DryDep = False

    except AssertionError:
        print('WARNING: No dry dep diagnostics found for {}'.format(wd))
        No_DryDep = True

    # Get convective scale wet deposition
    try:
        Cds = AC.get_WetLossConv_ds(wd=wd,
                                   dates2use=dates2use)

        # Add all all-NIT to dataset
        if any([('NIT' in i) for i in Specs]):
            Cds = AC.AddChemicalFamily2Dataset(Cds, fam='NIT-all',
                                                  prefix=WLprefix)
        No_C_WetDep = False
    except AssertionError:
        print('WARNING: No ConvWetDep diags. found for {}'.format(wd))
        No_C_WetDep = True

    # Get large scale wet deposition
    try:
        LSds = AC.get_WetLossLS_ds(wd=wd, dates2use=dates2use)

        # Add all all-NIT to dataset
        if any([('NIT' in i) for i in Specs]):
            LSds = AC.AddChemicalFamily2Dataset(LSds, fam='NIT-all',
                                              prefix=LSprefix)
        No_LS_WetDep = False
    except AssertionError:
        print('WARNING: No LSWetDep diagnostics found for {}'.format(wd))
        No_LS_WetDep = True

    # Use a copy of Dry Dep xr.Dataset as template to add new data too
    try:
        Dds
        ds = Dds.copy()
    except NameError:
        PrtStr = 'WARNING: Dep. output not found (Dry:{}, C Wet:{}, LS Wet:{})'
        print(PrtStr.format(No_DryDep, No_C_WetDep, No_LS_WetDep))
        return
    LongNameStr = 'Total (Dry and Wet) deposition flux of species {}'
    # Loop by requested species
    for Spec in Specs:
        # RMM
        SpecRMM = AC.species_mass(Spec)
        NewVar = '{}{}'.format(NewDepPrefix, Spec)
        long_name = LongNameStr.format(Spec)
        try:

            # Get dry dep deposition
            Var = '{}{}'.format(DDprefix, Spec)
            if verbose:
                print(Dds[Var])
            # Convert units to kg/m2/s (from )
            ExpectedUnits = 'molec cm-2 s-1'
            attrs = Dds[Var].attrs
            if (attrs['units'] == ExpectedUnits):
                if verbose:
                    PrtStr = "NOTE: updating '{}' DryDep units for {}"
                    print(PrtStr.format(Spec, wd))
                ds[NewVar] = Dds[Var].copy() / AVG * 1E-4 * SpecRMM
                units = 'kg m-2 s-1'
            else:
                PrtStr = "WARNING: Expected units of {}. Got '{}' for '{}'"
                print( PrtStr.format(ExpectedUnits, attrs['units'], wd) )
                sys.exit()

            # Get LS Wet dep
            Var = '{}{}'.format(LSprefix, Spec)
            try:
                if verbose:
                    print(LSds[Var])
                # Save to total dry/wet dep
                ds[NewVar] += LSds[Var].sum(dim='lev') / LSds['AREA']
                units = 'kg m-2 s-1'
            except KeyError:
                PrtStr = "WARNING: Skiping '{}' ({}) as not in dataset (wd:{})"
                print(PrtStr.format(Spec, Var, wd))

            # Get convective scale wet deposition
            Var = '{}{}'.format(WLprefix, Spec)
            try:
                if verbose:
                    print(Cds[Var])
                # Save to total dry/wet dep
                ds[NewVar] += Cds[Var].sum(dim='lev') / Cds['AREA']
                units = 'kg m-2 s-1'
            except KeyError:
                PrtStr = "WARNING: Skiping '{}' ({}) as not in dataset (wd:{})"
                print(PrtStr.format(Spec, Var, wd))

            # Update units
            attrs['units'] = units
            attrs['long_name'] = long_name
            ds[NewVar].attrs = attrs

        except KeyError:
            PrtStr = "WARNING: Dry and/or Wet dep '{}' not present for: '{}'"
            print( PrtStr.format(Spec, wd) )

    # Select a limited set of variables to return
    vars2rtn = [i for i in ds.data_vars if (DDVprefix not in i )]
    vars2rtn = [i for i in vars2rtn if (DDprefix not in i )]
    ds = ds[vars2rtn]

    # Check numbers if in verbose or debug mode
    if (verbose or debug):
        for Spec in Specs:
            NewVar = '{}{}'.format(NewDepPrefix, Spec)

            print(NewVar, (ds[NewVar] * ds['AREA']).values.sum() )

    # Return total deposition to the dictionary
    return ds[vars2rtn]


def get_NOx_budget_ds_dict_for_runs(ApplyMask=False,
                                    RunDict=None,
                                    RunSet='ACID', res='4x5',
                                    CoreRunsOnly=False,
                                    dates2use=None,
                                    trop_limit=False,
                                    ExtraConcVars=['NOx', 'NOy', 'NIT-all'],
                                    MaskName='inflow_CVAO_area',
                                    IncAllProdLossDiags=False,
                                    ConvertProdLossUnits=True,
                                    debug=False):
    """
    Get a single dictionary of xr.Datasets for models runs (NOx budget vars)
    """
    # Local variables
    SCprefix = 'SpeciesConc_'
    debug = True
    # Set runs to use and dates to of NetCDF files to use
    if isinstance(RunDict, type(None)):
        RunDict = get_dict_of_GEOSChem_model_output(res=res, RunSet=RunSet,
                                                    CoreRunsOnly=CoreRunsOnly,
                                                    folder4netCDF=True)
    # Retrieve tags for NOx budget and make an inverted dictionary
    TagD = get_tags_for_NOx_HONO(AllTags=True)
    TagDr = {v: k for k, v in list(TagD.items())}

    # Get core species concentrations
    ConcD = {}
    for key in RunDict.keys():
        ConcD[key] = AC.GetSpeciesConcDataset(wd=RunDict[key],
                                              dates2use=dates2use)
    # Get Loss/Sources vis Prod/Loss diagnostic
    ProdD = {}
    for key in RunDict.keys():
        try:
            ProdD[key] = AC.get_ProdLoss_ds(wd=RunDict[key],
                                            dates2use=dates2use)
        except AssertionError:
            print('WARNING: No Prod/Loss diagnostics found for {}'.format(key))

    # Get StateMet collection - for unit conversions... (And humidity)
    StateD = {}
    for key in RunDict.keys():
        StateD[key] = AC.get_StateMet_ds(wd=RunDict[key],
                                         dates2use=dates2use)
    # Get Photolysis rates
    JValD = {}
    for key in RunDict.keys():
        try:
            JValD[key] = AC.GetJValuesDataset(wd=RunDict[key],
                                              dates2use=dates2use)
        except AssertionError:
            print('WARNING: No JVal diagnostics found for {}'.format(key))

    # Get Emissions
    Specs = ['NO2', 'NO', ]
    EmissD = {}
    for key in RunDict.keys():
        try:
            EmissD[key] = AC.get_HEMCO_diags_as_ds(wd=RunDict[key],
                                                   dates2use=dates2use)
        except AssertionError:
            print('WARNING: No HEMCO diagnostics found for {}'.format(key))

    # Get dry and wet Deposition
    Specs = ['NO2', 'HNO3', 'NIT', 'NITs', 'NIT-all']
    DepD = {}
    for key in RunDict.keys():
        try:
#            DepD[key] = AC.get_DryDepAndWetDep_ds(wd=RunDict[key],
            DepD[key] = get_DryDepAndWetDep_ds(wd=RunDict[key], Specs=Specs,
                                                    dates2use=dates2use)
        except AssertionError:
            PrtStr = 'WARNING: No dry/wet dep diagnostics found for {}'
            print(PrtStr.format(key))

    # - Calculate the various quantities per run.
    Specs = ['NO2', 'NO', 'HNO2', 'HNO3',  'OH', 'TN', 'NIT']
    JNITvars = ['ProdHNO2fromHvNIT', 'ProdHNO2fromHvNITs']
    JNITvars += ['ProdHNO2fromHvNITD{}'.format(i) for i in np.arange(1, 5)]
    RMM_air = AC.constants('RMM_air')
    AVG = AC.constants('AVG')  # Avogadro constant (Mol^-1)
    NOxD = {}
    for key in RunDict.keys():
        # Use StateMet collection for physical variables (use a reference sim?)
        StateMet = StateD[key]
        prefix = 'Jval_'
        dsVars = StateMet.data_vars
        # Add moles as variable here to aid future calculations
        StateMet['Met_MOLES'] = StateMet['Met_AD'] / 1E3 / RMM_air
        vars2use = ['Met_AD', 'Met_AIRVOL', 'Met_AIRDEN', 'Met_PMID',
                    'FracOfTimeInTrop', 'Met_MOLES',
                    ]
        NOxD[key] = StateMet[vars2use + ['AREA']].copy()

        # Add species concentration
        ds = ConcD[key]
        Vars = ds.data_vars
        vars2use = [i for i in Vars if any(i.endswith(ii) for ii in Specs)]
        # Include any extra species / families?
        if (len(ExtraConcVars) > 0):
            if debug:
                print(ExtraConcVars)
            for ExtraVar in ExtraConcVars:
                if debug:
                    print(ExtraVar)
                # Check if family already present
                try:
                    ds[ExtraVar]
                except KeyError:
                    try:
                        ds = AC.AddChemicalFamily2Dataset(ds, fam=ExtraVar,
                                                          prefix=SCprefix)
                        if debug:
                            print(ExtraVar)
                            print(ds[ExtraVar])

    #                except Exception:
                    except:
                        # If extra concentration not setup as a family, then
                        # assume ExtraVar is a species variable
                        pass
                    vars2use += ['{}{}'.format(SCprefix, ExtraVar)]
                    if debug:
                        print(vars2use)
#            print('TODO: setup including additional species/families')
        vars2use = list(set(vars2use))
        ds = xr.merge([NOxD[key], ds[vars2use]])
        NOxD[key] = ds

        # Jvals from JVal diagnostic
        try:
            ds = JValD[key]
            Vars = ds.data_vars
            vars2use = [i for i in Vars if any(ii in i for ii in Specs)]
            ds = xr.merge([NOxD[key], ds[vars2use]])
            JScaleVar = 'Jscale'
            ds[JScaleVar] = ds['Jval_NIT'].copy()
            ds[JScaleVar] = ds[JScaleVar] / ds['Jval_HNO3']
            NOxD[key] = ds
        except:
            PrtStr = "WARNING: Failed to add Jvals for '{}': {}"
            print(PrtStr.format(key, RunDict[key]))

        # Add total wet/dry deposition fluxes
        try:
            ds = DepD[key]
            Vars = ds.data_vars
            vars2use = [i for i in Vars if any(ii in i for ii in Specs)]
            ds = xr.merge([NOxD[key], ds[vars2use]])
            NOxD[key] = ds

        except:
            PrtStr = "WARNING: Failed to add Dep values for '{}': {}"
            print(PrtStr.format(key, RunDict[key]))

        # HONO/NO2 production and general tagging
        try:
            ds = ProdD[key]
            Vars = list(ds.data_vars)
            if IncAllProdLossDiags:
                PLvars = ['Loss', 'Prod']
                vars2use = [i for i in Vars if any(ii in i for ii in PLvars)]
            else:
                vars2use = [i for i in Vars if any(ii in i for ii in Specs)]
            ds = xr.merge([NOxD[key], ds[vars2use]])
            # Rename the tags to be more descriptive.
            prefix = 'Prod_'
            OldVarnames = [i for i in vars2use if prefix+'T' in i]
            NewVarnames = [TagDr[i.split(prefix)[-1]] for i in OldVarnames]
            ds = ds.rename(dict(zip(OldVarnames, NewVarnames)))
            # Add JNIT-all as a variable
            NewVar = 'ProdHNO2fromHvNIT-all'
            ds[NewVar] = ds[JNITvars[0]].copy()
            attrs = ds[JNITvars[0]].attrs.copy()
            for var in JNITvars[1:]:
                ds[NewVar] = ds[NewVar] + ds[var]
            attrs['long_name'] = "Chemical production of all NIT"
            ds[NewVar].attrs = attrs
            NOxD[key] = ds
        except:
            PrtStr = "WARNING: Failed to add Prod/Loss for '{}': {}"
            print(PrtStr.format(key, RunDict[key]))

        # Get humidity and key variables from StateMet collection
#        'Met_SPHU'

    del ds
    gc.collect()

    # Convert units or prod/loss
    if ConvertProdLossUnits:
        AVG = AC.constants('AVG') # Avogadros constant
        PrtStr = "WARNING: No conversion for var ('{}') units ('{}')"
        for key in RunDict.keys():
            ds = NOxD[key]
            # Select variables to use
            vars2use = [i for i in ds.data_vars if 'Prod' in i]
            vars2use += [i for i in ds.data_vars if 'Loss' in i]
            for var in vars2use:
                if debug:
                    print(var)
                attrs = ds[var].attrs
                units = attrs['units']
                data = ds[var].copy()
                if debug:
                    print(units)
                if units == 'molec cm-3 s-1':
                    # convert to kg (N) / s-1
                    data = data / AVG * ds['Met_AIRVOL'] * 1E6 * 14 * 1E3
                    units = 'kg N s-1'
                    attrs['units'] = units
                    ds[var] = data
                    ds[var].attrs = attrs
                elif (units == 'kg N s-1'):
                    pass
                elif (units == 'kg s-1'):
                    # This is HNO3 loss oan sea-salt (are the units correct)
                    #                if sum_data:
                    pass
                else:
                    print(PrtStr.format(var, units))
            NOxD[key] = ds


    # Reduce the spatial focus of data to requested MaskName
    if ApplyMask:
        # Set spatial extents
        MaskDict = AC.GetMaskExtents(MaskName, )
        lowerlat = MaskDict['lowerlat']
        higherlat = MaskDict['higherlat']
        lowerlon = MaskDict['lowerlon']
        higherlon = MaskDict['higherlon']
        # Loop and reduce dataset spatially extent
        for key in NOxD.keys():
            ds = NOxD[key]
            # reduce area to CVAO...
            bool1 = ((ds.lon >= lowerlon) & (ds.lon <= higherlon)).values
            bool2 = ((ds.lat >= lowerlat) & (ds.lat <= higherlat)).values
            # Cut by lon, then lat
            ds = ds.isel(lon=bool1)
            ds = ds.isel(lat=bool2)
            NOxD[key] = ds

    # Limit the data to the troposphere
    if trop_limit:
        # Loop and reduce dataset scale
        for key in NOxD.keys():
            ds = NOxD[key]
            var2use = 'FracOfTimeInTrop'
            # Get variables with lev coordinate
            vars2use = [i for i in ds.data_vars if i != var2use]
            LevCoord4var = ['lev' in ds[i].coords for i in vars2use]
            vars2use = np.array(vars2use)[LevCoord4var]
            if debug:
                print(vars2use)
            # Mask subset of variables with boolean list
            __bool = ~(ds[var2use].values < 1)
            for __var in vars2use:
                ds[__var] = ds[__var].where(__bool.astype(bool))
            NOxD[key] = ds

    return NOxD

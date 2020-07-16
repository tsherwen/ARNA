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
from bs4 import BeautifulSoup
#
import matplotlib.pyplot as plt
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
import seaborn as sns
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from shapely.geometry.polygon import LinearRing
import matplotlib

# Import from elsewhere in ARNA module
from . core import *
from . utils import *
from . observations import *
from . plotting import plt_highres_modelling_region, add_scatter_points2cartopy_ax


def v12_9_TRA_XX_2_name(TRA_XX, RTN_dict=False):
    """
    Convert tracer number to tracer name in GEOS-Chem v12.9
    """
    TRAs = ['ACET', 'ACTA', 'AERI', 'ALD2', 'ALK4', 'ATOOH', 'BCPI', 'BCPO', 'BENZ', 'Br', 'Br2', 'BrCl', 'BrNO2', 'BrNO3', 'BrO', 'BrSALA', 'BrSALC', 'C2H6', 'C3H8', 'CCl4', 'CFC11', 'CFC113', 'CFC114', 'CFC115', 'CFC12', 'CH2Br2', 'CH2Cl2', 'CH2I2', 'CH2IBr', 'CH2ICl', 'CH2O', 'CH3Br', 'CH3CCl3', 'CH3Cl', 'CH3I', 'CH4', 'CHBr3', 'CHCl3', 'Cl', 'Cl2', 'Cl2O2', 'ClNO2', 'ClNO3', 'ClO', 'ClOO', 'CO', 'DMS', 'DST1', 'DST2', 'DST3', 'DST4', 'EOH', 'ETHLN', 'ETNO3', 'ETP', 'GLYC', 'GLYX', 'H1211', 'H1301', 'H2402', 'H2O', 'H2O2', 'HAC', 'HBr', 'HC5A', 'HCFC123', 'HCFC141b', 'HCFC142b', 'HCFC22', 'HCl', 'HCOOH', 'HI', 'HMHP', 'HMML', 'HNO2', 'HNO3', 'HNO4', 'HOBr', 'HOCl', 'HOI', 'HONIT', 'HPALD1', 'HPALD2', 'HPALD3', 'HPALD4', 'HPETHNL', 'I', 'I2', 'I2O2', 'I2O3', 'I2O4', 'IBr', 'ICHE', 'ICl', 'ICN', 'ICPDH', 'IDC', 'IDCHP', 'IDHDP', 'IDHPE', 'IDN', 'IEPOXA', 'IEPOXB', 'IEPOXD', 'IHN1', 'IHN2', 'IHN3', 'IHN4', 'INDIOL', 'INO', 'INPB', 'INPD', 'IO', 'IONITA', 'IONO', 'IONO2', 'IPRNO3', 'ISALA', 'ISALC', 'ISOP', 'ITCN', 'ITHN', 'LIMO', 'LVOC', 'LVOCOA', 'MACR', 'MACR1OOH', 'MAP', 'MCRDH', 'MCRENOL', 'MCRHN', 'MCRHNB', 'MCRHP', 'MEK', 'MENO3', 'MGLY', 'MOH', 'MONITA', 'MONITS', 'MONITU', 'MP', 'MPAN', 'MPN', 'MSA', 'MTPA', 'MTPO', 'MVK', 'MVKDH', 'MVKHC', 'MVKHCB', 'MVKHP', 'MVKN', 'MVKPC', 'N2O', 'N2O5', 'NH3', 'NH4', 'NIT', 'NITs', 'NO', 'NO2', 'NO3', 'NPRNO3', 'O3', 'OClO', 'OCPI', 'OCPO', 'OCS', 'OIO', 'PAN', 'pFe', 'PIP', 'PP', 'PPN', 'PROPNN', 'PRPE', 'PRPN', 'PYAC', 'R4N2', 'R4P', 'RA3P', 'RB3P', 'RCHO', 'RIPA', 'RIPB', 'RIPC', 'RIPD', 'RP', 'SALA', 'SALAAL', 'SALACL', 'SALC', 'SALCAL', 'SALCCL', 'SO2', 'SO4', 'SO4s', 'SOAGX', 'SOAIE', 'SOAP', 'SOAS', 'TOLU', 'XYLE']
    nums = np.arange(1, len(TRAs)+1)
    if RTN_dict:
        return dict(zip(nums, TRAs))
    else:
        return dict(zip(nums, TRAs))[TRA_XX]


def get_GEOSChem4flightnum(flight_ID='C225', resample_data=True):
    """
    Retrieve GEOS-Chem output for FAAM flight
    """
    # Where is the extract GEOS-CF data?
    RunDir = '/users/ts551/scratch/GC/rundirs/'
    SpecificRunStr = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.PF'
    folder = '{}/{}/'.format( RunDir, SpecificRunStr )
    # Extract the data for a specific flight
    files2use = glob.glob( os.path.join(folder,'*plane.log*') )
    # Get start date of flight - use the plane-flight from same day as sdate
    dfS = get_summary4flight(flight_ID=flight_ID)
    sdate = dfS.index.values.min()
    edate = dfS.index.values.max()
    sdate, edate = AC.dt64_2_dt([sdate, edate])
    sdate_str = sdate.strftime('%Y%m%d')
    file2use = [i for i in files2use if sdate_str in i ]
    assert len(file2use)==1, 'WARNING: more than one planeflight found!'
    file2use = file2use[0]
    # - Use standard AC_tool extraction - cannot as output issue in v12.9
    # Get Header infomation from first file
#    vars, sites = AC.get_pf_headers(file2use, debug=debug)
    # Extract all points from file
#    df, vars = AC.pf_csv2pandas(file=file2use, vars=vars, epoch=True,
#                                r_vars=True)
    # - Instead do this manually for now
    # Open file and save the data into
    with open(file2use, 'rb') as file:
        lines = [i for i in file]
    # Extract as raw data in chunks
    lines_1 = lines[0::2]
    header_1 = lines_1[0].decode('utf-8').split()
    data_1 = [i.decode('utf-8').split() for i in lines_1[1:] ]
    df = pd.DataFrame(data_1, columns=header_1)
    lines_2 = lines[1::2]
    header_2 = lines_2[0].decode('utf-8').split()
    data_2 = [i.decode('utf-8').split() for i in lines_2[1:] ]
    df2 = pd.DataFrame(data_2, columns=header_2)
    # Now combine
    df = pd.concat([df, df2], axis=1)
    # Now process the meta data/update type formats
    # TODO: below could be faster...
    # Use infer obects? - df.infer_objects
    dtypes = {'POINT':object, 'TYPE':str, 'YYYYMMDD':str, 'HHMM':str}
    cols2use = [i for i in df.columns if i not in dtypes.keys()]
    df[cols2use] = df[cols2use].apply(pd.to_numeric)
#     for col in df.columns:
#         try:
#             dtype = dtypes[col]
#         except KeyError:
#             dtype = np.float
#         df.loc[:, col] = df[col].astype(dtype)
    # Add a datetime index
    df = AC.DF_YYYYMMDD_HHMM_2_dt(df, rmvars=None, epoch=False)
    df.index.name = None
    # Add temperature in deg C
    df['T'] = df['GMAO_TEMP'].copy()
    df['T'] = df['GMAO_TEMP'].values - 273.15
    # Update the variable names
    d = v12_9_TRA_XX_2_name(None, RTN_dict=True)
    d = dict( [('TRA_{:0>3}'.format(i), d[i]) for i in d.keys()] )
    df = df.rename(columns=d)
    # Add NOx as combined NO and NO2
    df['NOx'] = df['NO'].values + df['NO2'].values
    # Add NOy as defined in GEOS-CF
    # NOy = no_no2_hno3_hno4_hono_2xn2o5_pan_organicnitrates_aerosolnitrates
    vars2use = [
    'BrNO3', 'ClNO3', 'ETHLN', 'ETNO3', 'HNO2', 'HNO3', 'HNO4', 'HONIT',
    'ICN', 'IDN', 'IHN1', 'IHN2', 'IHN3', 'IHN4', 'INDIOL',
    'INPB', 'INPD', 'IONITA', 'IONO', 'IONO2', 'IPRNO3', 'ITCN', 'ITHN',
    'MCRHN', 'MCRHNB', 'MENO3', 'MONITA', 'MONITS', 'MONITU', 'MPAN', 'MPN',
    'MVKN', 'N2O5', 'NIT', 'NITs', 'NO', 'NO2', 'NO3', 'NPRNO3', 'PAN',
    'PPN', 'PROPNN', 'R4N2',
    ]
    df['NOy'] = df[ 'N2O5'].copy() # 2 N2O5 in NOy, so 2x via template
    for var in vars2use:
        df.loc[:,'NOy'] = df['NOy'].values + df[var].values
    # Include a variable of NOy where HNO3 is removed
    # NOy = no_no2_hno3_hno4_hono_2xn2o5_pan_organicnitrates_aerosolnitrates
    df['NOy-HNO3'] = df['NOy'].values - df['HNO3'].values
    # Include a variable of NOy where HNO3 is removed
    df['NOy-HNO3-PAN'] = df['NOy'].values - df['HNO3'].values - df['PAN'].values
    # gas-phase (exc. PAN, HNO3, HNO4, Org-NIT, N2O5)
    df['NOy-Limited'] = df['NO'].values + df['NO2'].values + df['HNO2'].values + df['NIT'].values + df['NITs'].values
    # Uset the P-I variable as a model level variable
    df['model-lev'] = df['P-I'].copy()
	# Resample the data?
    if resample_data:
        df = df.resample('1T' ).mean()
    return df


def evaluate_regional_grid4GEOSChem(show_plot=False, dpi=320):
    """
    Evaluate the regional variable high(er) resolution (flex)grid for GEOS-chem
    """
    # - Load FAAM flighttrack data
    # Import all of the flights by campaign and plot on grid
    # All FAAM campaigns - included in James and Freya's analysis
    # NOTE: this will only use the downloaded files for now.
    campaign = 'ARNA-2'
    dfARNA = get_flighttracks4campaign(campaign)
    campaigns = [
    'ARNA-2', 'ACISIS-5', 'ACRUISE', 'ACISIS-4', 'MOYA-2', 'ACISIS-3',
    'ACISIS-2',
#    'Clarify', # FAAM files not downloaded
    'MOYA-1',
#    'ACISIS-1' # FAAM files not downloaded
    ]
    dfs = [get_flighttracks4campaign(i) for i in campaigns]

    # - Plot up high resolution modelling region around ARNA-2 flights
    savetitle = 'ARNA_high_resolution_model_grid'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    # Plot up a blank global background
    fig, ax = plt_highres_modelling_region(plot_blank_data=True)
    # Add the campaign flights to this
    campaign = 'ARNA-2'
    df = dfARNA
    LatVar = 'LAT_GIN'
    LonVar = 'LON_GIN'
    lats = df[LatVar].values
    lons = df[LonVar].values
    color_list = AC.get_CB_color_cycle()
    projection=ccrs.PlateCarree
    # Now scatter points on plot
    ax = add_scatter_points2cartopy_ax(ax=ax, lons=lons, lats=lats,
                                       color=color_list[0],
                                       label=campaign)

    fig.legend(loc=7)
    fig.suptitle('Flight-tracks during ARNA-2 campaign')
    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
    if show_plot:
        plt.show()
    plt.close()

    # - Plot up high resolution modelling region around BB FAAM flights
    fig, ax = plt_highres_modelling_region(plot_blank_data=True)
    for n_campaign, campaign in enumerate( campaigns ):
        df = dfs[campaigns.index(campaign)]
        print(campaign)
        lats = df[LatVar].values
        lons = df[LonVar].values
        ax = add_scatter_points2cartopy_ax(ax=ax, lons=lons, lats=lats,
                                           color=color_list[n_campaign],
                                           label=campaign)
    title = "Flight-tracks during FAAM campaigns in biomass burning analysis"
    fig.suptitle(title)
    fig.legend(loc=7)
    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
    if show_plot:
        plt.show()
    plt.close()

    # - Now evaluate biomass burning emissions for grid
    earth0_folder = '/mnt/lustre/groups/chem-acm-2018/earth0_data/'
    HEMCO_folder = earth0_folder + 'GEOS//ExtData/HEMCO/'
    # Get GFAS emissions for the year

    # Plot by seasons - TODO

    # Plot for feb 2020
    folder = HEMCO_folder+'/GFAS/v2018-09/2020/'
    filename = 'GFAS_202002.nc'
    ds = xr.open_dataset(folder+filename)
    # Update lon to be in degrees West -
    var2use = 'cofire'
    ds = ds[[var2use]]
    ds = ds.assign_coords({'lon':ds.lon.values -180})
    # Update name and scaling
    Uvar2use = '{} (1E-9 {})'.format(ds[var2use].long_name, ds[var2use].units)
    ds = ds.rename({var2use:Uvar2use})
    var2use = Uvar2use
    ds = ds[[var2use]].mean(dim='time') *1E9
    # Remove zero data
    arr = ds[var2use].values
    arr[arr <= 0] = np.NaN
    ds[var2use].values = arr
    # And Roll the variables too
#    ds = ds.roll(lon=-int(len(ds.lon)/2))
    arr = ds[var2use].values
    arr = np.roll(arr, -int(len(ds.lon)/2), axis=1)
    ds[var2use].values = arr
    # Plot up the data
    fig, ax = plt_highres_modelling_region(ds=ds,var2use=var2use,
                                           plot_blank_data=False,
                                           rm_colourbar=False )

    fig.suptitle('Biomass burning emissions (GFAS) - Feb 2020 (ARNA-2)')
    del ds
    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
    if show_plot:
        plt.show()
    plt.close()

    # - Now evaluate dust emission for grid
    # Plot by seasons - TODO

    # Plot for February
    folder = HEMCO_folder+'/OFFLINE_DUST/v2019-01/0.5x0.625/2019/02/'
    files2use = glob.glob(folder+'*nc')
    ds = xr.open_mfdataset(files2use)
    # Combine all dust emissions
    var2use = 'Total dust emission (kg/m2/s)'
    ds[var2use] = ds['EMIS_DST1'].copy()
    ds[var2use] = ds[var2use].values +ds['EMIS_DST2']
    ds[var2use] = ds[var2use].values +ds['EMIS_DST3']
    ds[var2use] = ds[var2use].values +ds['EMIS_DST4']
    ds = ds[[var2use]].mean(dim='time')
    # Remove zero data
    arr = ds[var2use].values
    arr[arr <= 0] = np.NaN
    ds[var2use].values = arr
    # Plot up the data
    fig, ax = plt_highres_modelling_region(ds=ds,var2use=var2use,
                                           plot_blank_data=False,
                                           rm_colourbar=False )

    fig.suptitle('Dust emissions (online) - Feb *2019* (ARNA-2)')
    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
    if show_plot:
        plt.show()
    plt.close()

    # - Now evaluate NOx emission for grid?

    # - Others variables to plot / consider?
    # Night lights?

    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def regrid_restart4ARNA_highres_grid():
    """
    Regrid a restart file for the ARNA grid
    """
    # File and location?
    RunRoot = '/users/ts551/scratch/GC/rundirs/'
    RunDir = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020.1x1.PF/'
    folder = RunRoot + RunDir
    InFile = 'GEOSChem.Restart.20200201_0000z.nc4.FILE2REGRID'
    ds = xr.open_dataset(folder+InFile)
    #
    lons = np.arange(-180, 180, 1)
    lats = np.arange(-90, 90, 1)
    OutFile = 'GEOSChem.Restart.20200120_0000z.REGRIDED.nc4'
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
    OutFile = 'GEOSChem.Restart.20200120_0000z.REGRIDED_CROPPED.nc4'
    ds.to_netcdf(os.path.join(folder,OutFile))


#!/usr/bin/python
"""
Driver for generic plotting functions following the ARNA campaign
"""
import arna as ar



import logging
import numpy as np
import pandas as pd
import time
import calendar
import datetime as datetime
from datetime import datetime as datetime_



def main():
    """
    Main driver function
    """
#    from arna import get_FAAM_core4flightnum, get_filters_data4flight, get_GEOSCF4flightnum, add_derived_FAAM_flags2df4flight, get_GEOSChem4flightnum, plt_flightpath_spatially_over_CVAO, get_CIMS_data4flight

    # Seaborn context for plots?
#    context = 'talk'
    context = 'paper'

    # - Plot up all comparisons by altitude together
#     ar.plt_comp_by_alt_4ARNA_together(context=context,
#                                       res='0.25x0.3125',
#                                       just_SLR=True,)
#     ar.plt_comp_by_alt_4ARNA_together(context=context,
#                                       res='0.25x0.3125',
#                                       just_SLR=False)
#     ar.plt_comp_by_alt_4ARNA_together(context=context,
#                                       res='4x5', RunSet=None,
#                                       just_SLR=False,)
    # temp for testing
#    flight_nums = []
#    RunSet = None
#    res = '4x5'
#    NOxAsLog = True
#    CoreRunsOnly = False
#    CoreRunsOnly = True
#    savetitle = 'ARNA_altitude_binned_combined_file_{}'.format(res)
#     ar.plt_comp_by_alt_4ARNA_together(context=context,
#                                       res=res, RunSet=RunSet,
#                                       flight_nums=flight_nums,
#                                       savetitle=savetitle,
#                                       just_SLR=False,
#                                       NOxAsLog=NOxAsLog,
#                                       CoreRunsOnly=CoreRunsOnly,
#                                       debug=True)

    # Output the same bulk plots for the ACID runs
    flight_nums = []
    RunSet = 'PostHONO'
    res = '4x5'
    NOxAsLog = True
    CoreRunsOnly = True
    savetitle = 'ARNA_altitude_binned_combined_file_{}_{}'
    savetitle = savetitle.format(res, RunSet)
    ar.plt_comp_by_alt_4ARNA_together(context=context,
                                      res=res, RunSet=RunSet,
                                      flight_nums=flight_nums,
                                      savetitle=savetitle,
                                      just_SLR=False,
                                      NOxAsLog=NOxAsLog,
                                      CoreRunsOnly=CoreRunsOnly,
                                      debug=True)

    # The same plots as above, but split by their own PDF file..
    # NOTE: below function fails with a ValueError
#     ar.plt_comp_by_alt_4ARNA_all(just_SLR=True, context=context,
#                                  RunSet='FP-Nest', inc_GEOSChem=True,
#                                  just_plot_GEOS_Chem=True,
#                                  res='0.25x0.3125', close_pdf=True)

    # Plot up data for SLRs with and without dust
#    ar.plt_comp_by_alt_4ARNA_all(just_SLR=False, context=context)
#    ar.plt_comp_by_alt_4ARNA_all(just_SLR=True, context=context)
#     ar.plt_comp_by_alt_4ARNA_all(just_SLR=False, context=context,
#                                  RunSet='FP-Nest', inc_GEOSChem=True,
#                                  just_plot_GEOS_Chem=True,
#                                  res='0.25x0.3125', close_pdf=True)
#
#     ar.plt_comp_by_alt_4ARNA_flights_CIMS(context=context,
#                                           RunSet='FP-Nest',
#                                           res='0.25x0.3125',
#                                           inc_GEOSChem=True,
#                                           just_SLR=True)
#
#     ar.plt_comp_by_alt_4ARNA_all_DUST(plt_model=False, context=context)
#     ar.plt_comp_by_alt_4ARNA_all_DUST(plt_model=True, context=context)
#     ar.plt_comp_by_alt_4ARNA_CIMS_all_DUST(context=context)

    # - Plot up core comparisons by flight as
    # As timeseries ...
#     ar.plt_ts_comp4ARNA_flights(inc_GEOSChem=False, context=context)
#     ar.plt_ts_comp4ARNA_flights(inc_GEOSChem=True, context=context,
#                                 just_plot_GEOS_Chem=True,
#                                 RunSet='FP-Nest', res='0.25x0.3125')
    # By altitude (and by flight)
#     ar.plt_comp_by_alt_4ARNA_flights(context=context)
#     ar.plt_comp_by_alt_4ARNA_flights(inc_GEOSChem=True, context=context,
#                                      just_plot_GEOS_Chem=True,
#                                      RunSet='FP-Nest', res='0.25x0.3125')

    # - Plot up ToF-CIMS data by flight as
    # As timeseries ...
#     ar.plt_ts_comp4ARNA_flights_CIMS(context=context)
#     ar.plt_ts_comp4ARNA_flights_CIMS(context=context,
#                                      RunSet='FP-Nest',
#                                      res='0.25x0.3125',
#                                      inc_GEOSChem=True,
#                                      flight_nums=flight_nums,
#                                      LatVar='LAT',
#                                      LonVar='LON',)
    # By altitude (and by flight)
#     ar.plt_comp_by_alt_4ARNA_flights_CIMS(context=context)
#     ar.plt_comp_by_alt_4ARNA_flights_CIMS(context=context,
#                                           RunSet='FP-Nest',
#                                           res='0.25x0.3125',
#                                           inc_GEOSChem=True,
#                                           )

    # - Plot up nitrate aerosol data by flight as
    # As timeseries ...
#     ar.plt_ts_comp4ARNA_flights_filters(context=context)
#     ar.plt_ts_comp4ARNA_flights_filters(context=context,
#                                         RunSet='FP-Nest',
#                                         res='0.25x0.3125',
#                                         inc_GEOSChem=True,
#                                         LatVar='LAT',
#                                         LonVar='LON',)
#     ar.plt_ts_comp4ARNA_flights_filters(context=context,
#                                         res='4x5',
#                                         inc_GEOSChem=True,
#                                         LatVar='LAT',
#                                         LonVar='LON',
#                                         )

    # Plot up nitrate, JNIT, and their project
#    AC.mk_tri_NO3_JNIT_combination_plt()

    # - Plot up SWAS data by flight
    # Plot up SWAS data
#     ar.plt_ts_comp4ARNA_flights_SWAS(context=context)

    # - Plot up PCASP/CDP data by flight as
    # NOTE: CAS data being ignored currently due to issue with mirror window
    # As timeseries ...
#    ar.plt_ts_comp4ARNA_flights_PCASP()

    # - Plot up velocity and Roll, amongst other core physical vars by flight
    # As timeseries ...
#     ar.plt_ts_comp4ARNA_flights_PHYSICAL_VARS(context=context)
    ar.plt_ts_comp4ARNA_flights_PHYSICAL_VARS(context=context,
                                              just_plot_GEOS_Chem=True,
                                              inc_GEOSChem=True,
                                              res='0..25x0.3125',
                                              RunSet='FP-Nest')
    # Plot up the temperature data from Hannah Price
    # N/A? this is only for 2019. Data to be worked up for 2020.

    # - Plot up SWAS data by flight
    # Plot a comparison of NOy
#     ar.plt_ts_comp4ARNA_flights_NOy_ALL(context=context)
#     ar.plt_ts_comp4ARNA_flights_NOy_ALL(context=context,
#                                         RunSet='FP-Nest',
#                                         res='0.25x0.3125',
#                                         inc_GEOSChem=True,
#                                         LatVar='LAT',
#                                         LonVar='LON',)

    # - Other misc. plotting tasks
#    explore_high_ozone_near_CVAO()
#    extract_GEOS54all_ARNA_flights()

    # Evaluate the high resolution modelling region
#     ar.evaluate_regional_grid4GEOSChem()

    # Also plot up for related biomass-burning flights in MOYA campaign
#     ar.plt_ts_comp4MOYA_flights()
#     ar.plt_ts_comp4MOYA_flights_PHYSICAL_VARS()

    # Plot seasonal and vertical comparisons of nitrate (CVAO)
    ar.plt_seasonal_comparisons_of_nitrate()
    ar.mk_vertical_comparisons_with_nitrate()

    # Do the planeflight Jscale analysis
    do_planeflight_campaign_Janalysis()

    # Summarise the relevant deposition and emissions numbers
    summarise_emissions4RunDict(RunDict=RunDict, dates2use=dates2use)
    summarise_deposition4RunDict(RunDict=RunDict, dates2use=dates2use)

    # - Plot up observational comparisons
    plt_comp_with_NASA_Atom()
    plt_comp_with_FIREX_AQ()
    plt_obs_based_FIREX_analysis()


def do_planeflight_campaign_Janalysis(flight_nums=[], run2use=None):
    """
    Analysis nitrate photolysis within model planeflight output
    """
    # local variables
    CoreRunsOnly = True
    RunSet = 'IGACset.tagged'
    # Which run to do the full analysis for?
    if isinstance(run2use, type(None)):
#        run2use = 'Iso.UnlimitedAll'
        run2use = 'Iso.Unlimited'

    # Which flights to plot? - Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [218, 219, 220, 221, 222, 223, 224, 225, ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]

    # Retrieve GEOS-Chem planeflight data
    dfs_mod_GC = {}
    for flight_ID in flight_IDs:
        dfs = ar.get_GEOSChem4flightnum(flight_ID=flight_ID,
                                     res=res,
                                     CoreRunsOnly=CoreRunsOnly,
                                     RunSet=RunSet,)
        for key in dfs.keys():
            df = dfs[key]
            df = ar.add_derived_FAAM_flags2df4flight(df=df,
                                                  flight_ID=flight_ID)
            df['flight_ID'] = flight_ID

            dfs[key] = df
        dfs_mod_GC[flight_ID] = dfs
    dfs_mod = dfs_mod_GC

    # Combine to a single dataframe (dictionary lists are by flight )
    df_mod = pd.concat([dfs_mod[i][run2use] for i in dfs_mod.keys()], axis=0)
    dfs_mod_ALL = {}
    if len(dfs.keys()) > 1:
        for key in dfs.keys():
            ModByFlight = [dfs_mod[i][key] for i in dfs_mod.keys()]
            dfs_mod_ALL[key] = pd.concat(ModByFlight, axis=0)

    # Only consider data during SLRs?
    if just_SLR:
#        df_obs = df_obs.loc[df_obs['IS_SLR'] == True, :]
        for key in dfs_mod_ALL.keys():
            df_mod = dfs_mod_ALL[key]
            df_mod = df_mod.loc[df_mod['IS_SLR'] == True, :]
            dfs_mod_ALL[key] = df_mod
        extr_str = '_JUST_SLR'
    else:
        extr_str = ''

    # Update Name for JHNO2, JHNO3 and add JScale
    # (Calculate the enhancement factor)
    JHNO3var = 'JHNO3'
    JHNO2var = 'JHNO2'
    JNITvar = 'JNIT'
    NITvar = 'NIT-all'
    JScale = 'JScale'
    NameDict = {'JVL_015': JHNO2var, 'JVL_016': JHNO3var,'JVL_132':JNITvar }
    for key in list(dfs_mod_ALL.keys()):
        df = dfs_mod_ALL[key]
        df = df.rename(mapper=NameDict, axis=1)
        df[JScale] = df[JNITvar] / df[JHNO3var]
        dfs_mod_ALL[key] = df


    # - Plot enhancement factor against pNO3
    fig, ax = plt.subplots()
    SaveName = 'ARNA_modelled_JScale_versus_NIT'

    X = df[NITvar].values * 1E12
    Y = df[JScale].values
    plt.scatter(X, Y,)
    ax.set_xlabel("[pNO$_3^-$]$_\mathrm{bulk}$ (pptv)")
    ax.set_ylabel('JScale')#, fontsize=fontsize)

    # Include a line for the isotherm?
    no3 = np.arange(0.01, 10000, 0.01) # nmoles/m3
    f11 = (100 *  19.48  * 0.43)   / (1. + ( 0.43 * no3 ) )
    # convert no3 in nmoles/m3 to ppt (assuming air den)
#    AIRDEN = AC.constants( 'AIRDEN' )  # g/cm3
    AIRDEN = 0.001225  # g/cm3
    RMM_air = AC.constants( 'RMM_air' )
    #  (1/(g/mol)) = (mol/g) ; (mol/g) * (g/cm3) = mol/cm3
    MOLS = (1/RMM_air) * AIRDEN
    MOLS = MOLS * 1E6
    no3_pptv = (no3 / 1E12) / MOLS *1E12
#    plt.plot(no3_pptv, f11, label='Isotherm (v5)')
#    plt.legend()

    AC.save_plot(dpi=320, title='{}{}'.format(SaveName, extr_str))
    plt.close()

    # plot enhancement factor against JHNO3 * pNO3
    fig, ax = plt.subplots()
    SaveName = 'ARNA_modelled_JScale_versus_NITxJHNO3'

    X = (df[JHNO3var].values*60.*60.) * (df[NITvar].values * 1E12)
    Y =  df[JScale].values
    plt.scatter(X, Y,)
    label = r"$\it{j}_\mathrm{HNO_3}$ $\times$ "
    label += r"[pNO$_3^-$]$_\mathrm{bulk}$ (ppt h$^{-1}$)"
    ax.set_xlabel(label)#, fontsize=fontsize)
    ax.set_ylabel('JScale')#, fontsize=fontsize)

    # Include a line for the isotherm?

    AC.save_plot(dpi=320, title='{}{}'.format(SaveName, extr_str))
    plt.close()


def summarise_emissions4RunDict(RunDict=None, dates2use=None):
    """
    Summarise the modelled emissions for provided runs
    """
    # Local variables and checking of variables provided
    if isinstance(RunDict, type(None)):
        print('WARNING: Dictionary of model output locations required')
        sys.exit()

    # Extract the emissions NetCDFs to xr.datasets
    EmissD = {}
    for key in RunDict.keys():
        EmissD[key] = AC.get_HEMCO_diags_as_ds(wd=RunDict[key],
                                               dates2use=dates2use)

    # Convert values for NOx and return
    df = pd.DataFrame()
    for key in RunDict.keys():
        ds = EmissD[key]

        # Loop and store all NOx linked variables
        vars2use = [i for i in ds.data_vars if 'NO' in i]
        for __var in vars2use:

            # Get species
            Species = __var.split('_')[0].split('Emis')[-1]
            dsMonths = ds['time.month'].values
#            if dsMonths != np.arange(1,13):
            if len(dsMonths) != 12:
                Pstr = 'WARNING: 12 months not found in dataset ({}): {}'
                print(Pstr.format(key, dsMonths))

            # Sum over altitude if in code
            if 'lev' in ds[var].coords:
                ds = ds.sum(dim='lev')

            # Convert numbers to Tg (N) / yr
            data = ds['AREA'] * ds[__var] / AC.species_mass(Species) * 14.0

            # Convert to annual numbers
            data = data * 60.*60.*24.*365 / len(dsMonths)

            # Save to DataFrame
            df.loc[__var, key] = data.sum() * 1E3 / 1E12

    # Save to csv
    ExtStr = RunSet
    SaveName = 'ARNA_Emissions_summary_{}.csv'.format(ExtStr)
    df.to_csv(SaveName)


def summarise_deposition4RunDict(RunDict=None, dates2use=None):
    """
    Summarise the modelled emissions for provided runs
    """
    from arna import get_DryDepAndWetDep_ds
    # Local variables and checking of variables provided
    if isinstance(RunDict, type(None)):
        print('WARNING: Dictionary of model output locations required')
        sys.exit()

    # Extract the emissions NetCDFs to xr.datasets
    Specs = ['NO2', 'HNO3', 'NIT', 'NITs', 'NIT-all']
    DepD = {}
    for key in RunDict.keys():
        try:
#            DepD[key] = AC.get_DryDepAndWetDep_ds(wd=RunDict[key],
            DepD[key] = get_DryDepAndWetDep_ds(wd=RunDict[key], Specs=Specs,
                                                    dates2use=dates2use,
                                                    debug=True)
        except AssertionError:
            PrtStr = 'WARNING: No dry/wet dep diagnostics found for {}'
            print(PrtStr.format(key))


    # Convert values for NOx and return
    df = pd.DataFrame()
    prefix = 'DryAndWetDep_'
    for key in DepD.keys():
        ds = DepD[key]
        #
        if isinstance(ds, type(None)):
            continue
        # Loop and store all NOx linked variables
        vars2use = [i for i in ds.data_vars if prefix in i]
        print(vars2use)
        for __var in vars2use:
            # Get species
            Species = __var.split(prefix)[-1]
            dsMonths = ds['time.month'].values
#            if dsMonths != np.arange(1,13):
            if len(dsMonths) != 12:
                Pstr = 'WARNING: 12 months not found in dataset ({}): {}'
                print(Pstr.format(key, dsMonths))

            # Sum over altitude if in code
            if 'lev' in ds[__var].coords:
                ds[__var] = ds[__var].sum(dim='lev')

            # Convert numbers to Tg (N) / yr
            data = ds['AREA'] * ds[__var] / AC.species_mass(Species) * 14.0

            # Convert to annual numbers
            data = data * 60.*60.*24.*365 / len(dsMonths)

            # Save to DataFrame
            df.loc[__var, key] = data.sum() * 1E3 / 1E12


    # Save to csv
    ExtStr = RunSet
    SaveName = 'ARNA_deposition_summary_{}.csv'.format(ExtStr)
    df.to_csv(SaveName)



def explore_high_ozone_near_CVAO():
    """
    High ozone seen in some locations just in nested GEOS output
    """
    import seaborn as sns
    import gc
    # - plot up locations of high ozone by flight.
    # Use the high res model
    RunSet = 'FP-Nest'
    res = '0.25x0.3125'
    # Which flights to plot? - Just use non-transit ARNA flights
    flights_nums = [218, 219, 220, 221, 222, 223, 224, 225, ]
    flight_IDs = ['C{}'.format(i) for i in flights_nums]
    # Get model data
    dfs_mod_GC = {}
    for flight_ID in flight_IDs:
        dfs_mod_GC[flight_ID] = ar.get_GEOSChem4flightnum(flight_ID=flight_ID,
                                                          res=res,
                                                          RunSet=RunSet,)
    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        dfs_obs[flight_ID] = ar.get_FAAM_core4flightnum(flight_ID=flight_ID)

    # Now plot high ozone locations by flight
    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_flighttrack_high_ozone'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    for flight_ID in flight_IDs:
        #    for flight_ID in flight_IDs[-1:]:
        # Get dataframes for flight
        df_obs = dfs_obs[flight_ID]
        df_mod = dfs_mod_GC[flight_ID]
        # subselect locations where ozone > 70 ppbv
        df2plot = df_mod['FP-Nest']
        LatVar = 'LAT'
        LonVar = 'LON'
        cut_off = 70
        df2plot = df2plot.loc[df2plot['O3'] > cut_off*1E-09, :]
        title_str = 'Locations with modelled ozone >{} ppbv during {}'
        if (df2plot.shape[0] != 0):
            prt_str = 'Plotting {} for values >{} ppbv'
            print(prt_str.format(flight_ID, cut_off))

            title = title_str.format(cut_off, flight_ID)
            ar.plt_flightpath_spatially_over_CVAO(df=df2plot,
                                                  flight_ID=flight_ID,
                                                  title=title,
                                                  LatVar=LatVar, LonVar=LonVar)

            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        else:
            prt_str = 'WARNING: Not plotting for {} as no O3 ppbv values >{}'
            print(prt_str.format(flight_ID, cut_off))
    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')

    # - Check the average modelled ozone vertical profile at Cape verde

    #
    RunSet = 'FP-Nest'
    res = '0.25x0.3125'
    RunDict = ar.get_dict_of_GEOSChem_model_output(res=res, RunSet=RunSet)
    Folder = '{}/OutputDir/'.format(RunDict[RunSet])

    #
    ds = AC.get_GEOSChem_files_as_ds(wd=Folder)
#    StateMet = AC.get_StateMet_ds(wd=Folder)
    ModelAlt = AC.gchemgrid('c_km_geos5')
    prefix = 'SpeciesConc_'
    specs2plot = [i for i in ds.data_vars if prefix in i]
    specs2plot = [i.split(prefix)[-1] for i in specs2plot][::-1]
#    specs2plot = ['O3', 'CO', 'NO2'][::-1]

    # Setup PDF to save PDF plots to
    import seaborn as sns
    sns.set(color_codes=True)
    sns.color_palette('colorblind')
    sns.set_context(context)
    savetitle = 'ARNA_vertical_above_CVAO_GEOSChem_campaign_01'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    save2png = False
    show_plot = False
    font_scale = 1

    from funcs4obs import gaw_2_loc
    for spec in specs2plot[:100]:
        #    for spec in specs2plot[-100:]:
        #    for spec in specs2plot[-100:-90]:
        site = 'CVO'
        lat, lon, alt, TZ = gaw_2_loc(site)
        ds_tmp = ds[prefix+spec].sel(lat=lat, lon=lon, method='nearest')
        ds_tmp = ds_tmp.mean(dim='time')
        ds_tmp['lev'] = ModelAlt
        units, scalby = AC.tra_unit(spec, scale=True)
        ds_tmp *= scalby
        #  plot up...
        ds_tmp.plot(y='lev')
        plt.ylabel('{} ({})'.format('Altitude', 'km'))
        plt.xlabel('{} ({})'.format(spec, units))
        plt.title('Vertical Profile at CVAO during ARNA-2')
        plt.ylim(0, 15)
        if spec == 'O3':
            plt.xlim(-20, 200)
        elif spec == 'NO2':
            plt.xlim(-0.05, 0.2)
        elif spec == 'CO':
            plt.xlim(40, 160)
        if save2png:
            save_str = 'ARNA_vertical_above_CVAO_GEOSChem_{}'
            AC.save_plot(save_str.format(spec))
            AC.close_plot()
        else:
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            if show_plot:
                plt.show()
            plt.close()
        # Do some memory management...
        gc.collect()

    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')

    # - Explore high ozone in vertical column in coarse output
    site = 'CVO'
    from funcs4obs import gaw_2_loc
    lat, lon, alt, TZ = gaw_2_loc(site)
    NIU, NIU, ModelAlt = AC.get_latlonalt4res('4x5', full_vert_grid=True)
    ModelhPa = AC.hPa_to_Km(ModelAlt, reverse=True)

    def hPa_to_Km_reverse(value, reverse=True):
        return AC.hPa_to_Km(value, reverse=reverse)

    def hPa_to_Km_local(value, reverse=False):
        return AC.hPa_to_Km(value, reverse=reverse)

    # dates2use
    dates2use = [datetime.datetime(2020, 1, 1+i) for i in range(31)]
    dates2use += [datetime.datetime(2020, 2, 1+i) for i in range(29)]

    # Model runs to use?
    res = '4x5'
    RunSet = 'ACID'
    CoreRunsOnly = True
    folder4netCDF = True
    RunDict = ar.get_dict_of_GEOSChem_model_output(res=res, RunSet=RunSet,
                                                   CoreRunsOnly=CoreRunsOnly,
                                                   folder4netCDF=folder4netCDF)
    # Add other runs?
    # merra
    RunRoot = ar.get_local_folder('RunRoot')
    RunStr = '/merra2_4x5_standard.v12.9.0.BASE.2019.2020.diags/OutputDir/'
    RunDict['MERRA.4x5'] = RunRoot + RunStr

    # no biomass burning
    RunStr = '/geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.DustUptake'
    RunStr += '.JNIT.Isotherm.BCs.repeat.ON.II.diags.J50.BBx0/OutputDir/'
    RunDict['Acid-4x5-J50-NoBB'] = RunRoot + RunStr

    # limit plotted values to within the month before the campaign end
    sdate = datetime.datetime(2020, 1, 15)
    edate = datetime.datetime(2020, 2, 15)
    limit_time2campaign = True

    # Get data
    prefix = 'SpeciesConc_'
    dsD = {}
    for key in RunDict.keys():
        ds = AC.GetSpeciesConcDataset(wd=RunDict[key],
                                      dates2use=dates2use)
        ds = AC.AddChemicalFamily2Dataset(ds, fam='NOx', prefix=prefix)
        if limit_time2campaign:
            bool1 = AC.dt64_2_dt(ds.time) >= sdate
            bool2 = AC.dt64_2_dt(ds.time) <= edate
            ds = ds.isel(time=bool1)
            ds = ds.isel(time=bool2)
        dsD[key] = ds

    # Get GEOS-CF data and...
#    dates2use = [datetime.datetime(2020, 2, 1+i) for i in range(29)]
    dsCF = ar.get_GEOS_assim_expanded_dataset4ARNA(dts=dates2use)
#    dsCF['NOx'] = dsCF['NO'] + dsCF['NO2']
#    dsCF = ar.add_extra_derived_species(dsCF)
    dsCF = AC.AddChemicalFamily2Dataset(dsCF, fam='NOx', prefix='')
    # Update altitude
#    HPa_l = AC.get_GEOSCF_vertical_levels(native_levels=True)
#    hPa_as_km = [i for i in AC.hPa_to_Km(HPa_l)]

    # ... save via disk
    VarNames = ['O3', 'CO', 'NOx']
    ds2plotCF = dsCF[VarNames].sel(lat=lat, lon=lon, method='nearest')
    if limit_time2campaign:
        bool1 = AC.dt64_2_dt(ds2plotCF.time) >= sdate
        bool2 = AC.dt64_2_dt(ds2plotCF.time) <= edate
        ds2plotCF = ds2plotCF.isel(time=bool1)
        ds2plotCF = ds2plotCF.isel(time=bool2)
    ds2plotCF = ds2plotCF.mean(dim='time')
    savename = 'TEMP_NetCDF_GEOSCF_VI.nc'
    ds2plotCF = AC.save_ds2disk_then_reload(ds2plotCF, savename=savename)

    # plot up
    savetitle = 'ARNA_vertical_above_CVAO_GEOSChem_campaign_period2020'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # plot up GEOS-Chem runs
    prefix = 'SpeciesConc_'
    for VarName in VarNames:
        for key in RunDict.keys():
            ds = dsD[key][[prefix+VarName]]
            ds = ds.sel(lat=lat, lon=lon, method='nearest')
            X = ds[prefix+VarName].values[0] * 1E9
            Y = np.array(ModelhPa[:len(X)])
#            Y = ds.lev.values
            plt.plot(X, Y, label=key)

        # Plot up GEOS-CG
        X = ds2plotCF[VarName].values * 1E9
        Y = ds2plotCF.lev.values
        plt.plot(X, Y, label='GEOS-CF')

        # Beautify and save
        ax = plt.gca()
        if VarName == 'O3':
            plt.xlim(0, 100)
        elif VarName == 'NOx':
            plt.xlim(0, 1)
            ax.set_xscale('log')

#        ax.invert_yaxis()
        # Add a twin y axis
#         secax = ax.secondary_yaxis('right',
#                                    functions=(hPa_to_Km_local,
#                                    hPa_to_Km_reverse)
#                                    )
        plt.ylim()
        ax.invert_yaxis()
#        ax.set_yscale('log')

        ax.set_ylabel('Pressure altitude (hPa)')
        print(ax.get_ylim())
        ax.set_ylim(1000, 100)

        secax = ax.secondary_yaxis('right', functions=(AC.hPa2Km, AC.km2hPa))
        secax.set_ylabel('Altitude (km)')
#        secax.set_yscale('linear')
        print(secax.get_ylim())
        secax.set_ylim(0, 20)

        plt.title(VarName)
        plt.legend()
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)  # , tight=True)
        plt.close()

    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def extract_GC_data4CVAO():
    """
    Extract model data for CVAO
    """
    # - Get photolysis surface data for Simone
    RunRoot = ar.get_local_folder('RunRoot')
    folder = 'geosfp_4x5_standard.v12.9.0.BASE.2019.2020.ARNA.BCs.'
    folder += 'TEST.PF_Jrates.JVALS.GLOBAL/OutputDir/'
    folder = '{}{}'.format(RunRoot, folder)
    # Open the dataset
    ds = AC.GetJValuesDataset(wd=folder)
    # Extract for CVAO
    site = 'CVO'
    lat, lon, alt, TZ = gaw_2_loc(site)
    ds_tmp = ds.sel(lat=lat, lon=lon, method='nearest')
    ds_tmp = ds_tmp.sel(lev=ds_tmp.lev.values[0])
    # save out to csv.
    vars2use = [i for i in ds.data_vars if 'UV' in i]
    vars2use += [i for i in ds.data_vars if 'Jval_' in i]
    vars2use = list(set(vars2use))
    # Reduce to the variables of intereest and then save to disk
    ds_tmp = ds_tmp[vars2use].squeeze()
    del ds_tmp['lat']
    del ds_tmp['lon']
    del ds_tmp['lev']
    ds_tmp = AC.save_ds2disk_then_reload(ds_tmp)
    df = ds_tmp.to_dataframe()
    df.to_csv('GC_JValues_Collection_4x5_FP-GLOBAL_CVAO.csv')


def test_new_planeflight_Jrate_output():
    """
    Test online rxn. output via plane-flight diagnostic from GEOS-Chem
    """
    from funcs4obs import gaw_2_loc
    # - Setup sites to use
    df = pd.DataFrame()
    GAW_sites = [
        'ASK', 'BRW', 'CGO', 'CMN', 'CPT', 'CVO', 'JFJ', 'LAU', 'MHD', 'MLO',
        'MNM', 'NMY', 'SMO', 'SPO', 'THD'
    ]

    # - Get plane flight output
    RunRoot = ar.get_local_folder('RunRoot')
    folder = RunRoot +'/geosfp_4x5_standard.v12.9.0.BASE.2019.2020.ARNA.BCs'
    folder += '.TEST.PF_Jrates.REA.VI/'
    files2use = list(sorted(glob.glob(folder + '/TEST_1day/*plane*')))
    file2use = files2use[0]
    # Get Header infomation from first file
    vars, sites = AC.get_pf_headers(file2use, debug=debug)
    # Extract all points from file
    dfs = []
    for file2use in files2use:
        print(file2use)
        df, NIU = AC.pf_csv2pandas(file=file2use, vars=vars, epoch=True,
                                   r_vars=True)
        dfs += [df]
    # Append the dataframes together
    df = dfs[0].append(dfs[1:])
    df = AC.DF_YYYYMMDD_HHMM_2_dt(df, rmvars=None, epoch=False)
    df.index.name = None

    # Process and save csv files by date
    filename = 'GC_planeflight_data_FP-GLOBAL_JVALS_ARNA1_{}.csv'
    for file2use in files2use:
        date = file2use.split('plane.log.')[-1]
        print(file2use)
        vars, sites = AC.get_pf_headers(file2use, debug=debug)
        df, NIU = AC.pf_csv2pandas(file=file2use, vars=vars, epoch=True,
                                   r_vars=True)
        # Update the datetime index
        df = AC.DF_YYYYMMDD_HHMM_2_dt(df, rmvars=None, epoch=False)
        df.index.name = None
        df.to_csv(filename.format(date))

    # - Get output from NetCDF diagnostics
    # Get Get J-values
    folder = RunRoot+'/geosfp_4x5_standard.v12.9.0.BASE.2019.2020.ARNA.BCs.TEST.PF_Jrates.REA.III/OutputDir/'
    ds = AC.GetJValuesDataset(wd=folder)
    # And get StateMet
    FileStr = 'GEOSChem.StateMet.*'
    ds2 = xr.open_mfdataset(folder+FileStr)
    ds2 = ds2['Met_PMID']
    ds = xr.merge([ds2, ds])

    s_d = {}
    for site in GAW_sites:
        # Get: lat, lon, alt (press), timezone (UTC)
        lat, lon, alt, TZ = gaw_2_loc(site)
        # Use nearest grid box
#                lev2use = AC.find_nearest_value(HPa_r, alt)

        # Get the closest location
        ds_tmp = ds.sel(lat=lat, lon=lon, method='nearest')
        # then select height
#                lev=ds.lev.values[lev2use]
        lev2use = ds_tmp['Met_PMID'].mean(dim='time')
        first_level = int(lev2use.values[0])
        lev2use = AC.find_nearest(lev2use.values, alt)
#                int(lev2use.sel(lev=alt, method='nearest').values)
        ds_tmp = ds_tmp.isel(lev=lev2use)
        print(site, alt, lev2use, lev2use == 0)
        # extract data
        prefix = 'Jval_'
        vars2use = [i for i in ds_tmp.data_vars if prefix in i]
        ds_tmp = ds_tmp[vars2use]
        name_dict = zip(vars2use, [i.split(prefix)[-1] for i in vars2use])
        name_dict = dict(name_dict)
        ds_tmp = ds_tmp.rename(name_dict=name_dict)
        #
        for coord in [i for i in ds_tmp.coords if i != 'time']:
            del ds_tmp[coord]
        S = ds_tmp.to_dataframe()
        s_d[site] = S

    # - plot up
    specs2plot = ['NO2', 'HNO3', 'HNO2', 'BrO', 'IO', 'CH2I2', 'O3', 'PAN']
    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_test_Jvals'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    for site in GAW_sites:
        print(site)
        #
        nc_data = s_d[site]
        pf_data = df.loc[df['TYPE'] == site, :]
        for spec in specs2plot:
            print(spec)

            # Plot up PF
            pf_spec = None
            if spec == 'BrO':
                pf_spec = 'JVL_028'
            elif spec == 'O3':
                pf_spec = 'JVL_002'
            elif spec == 'NO2':
                pf_spec = 'JVL_011'
            elif spec == 'IO':
                pf_spec = 'JVL_116'
            elif spec == 'CH2I2':
                pf_spec = 'JVL_123'
            elif spec == 'PAN':
                pf_spec = 'JVL_059'
            elif spec == 'HNO3':
                pf_spec = 'JVL_016'
            elif spec == 'HNO2':
                pf_spec = 'JVL_015'
            elif spec == 'O1D':
                pf_spec = 'JVL_002'  # Also ???
            else:
                print('case not setup for species: {}'.format(spec))
            if isinstance(pf_spec, str):
                print(spec, pf_spec)
                data = pf_data[pf_spec]
                plt.plot(data.index, data.values, label='PF')
            # Plot up NC
            data = nc_data[spec]
            plt.plot(data.index, data.values, label='NetCDF')
            plt.legend()
            plt.title('{} @ {}'.format(spec, site))
            # save
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()

    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')

    # - Plot up differences between J rates in nested and global run
    savetitle = 'ARNA_Jvals_Global_vs_nest_model'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    RunRoot = ar.get_local_folder('RunRoot')
    RunStr = 'geosfp_4x5_standard.v12.9.0.BASE.2019.2020.'
    folderNest = RunRoot + RunStr + 'ARNA.Nest.repeat.JVALS/'
    folderGlobal = RunRoot + RunStr + 'ARNA.BCs.TEST.PF_Jrates.JVALS.GLOBAL/'
    files2useNest = list(sorted(glob.glob(folderNest + '/*plane*log*')))
    files2useGlobal = list(sorted(glob.glob(folderGlobal + '/*plane*log*')))
    for nfile2use, file2useNest in enumerate(files2useNest):
        date = file2use.split('plane.log.')[-1]
        print(file2use)
        file2useGlobal = files2useGlobal[nfile2use]
        file2use_dict = {'Nest': file2useNest, 'Global': file2useGlobal}
        dfs = {}
        for key in file2use_dict.keys():
            file2use = file2use_dict[key]
            vars, sites = AC.get_pf_headers(file2use, debug=debug)
            df, NIU = AC.pf_csv2pandas(file=file2use, vars=vars, epoch=True,
                                       r_vars=True)
            # Update the datetime index
            df = AC.DF_YYYYMMDD_HHMM_2_dt(df, rmvars=None, epoch=False)
            df.index.name = None
            dfs[key] = df
        # Now plot
        for nkey, key in enumerate(list(dfs.keys())):
            dfs[key]['JVL_134'].plot(label=key, color=['blue', 'red'][nkey])
            plt.title('JHNO3 for flighttrack on {}'.format(date))
        # save
        plt.ylabel('J-rate (s$^{-1}$)')
        plt.legend()
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()
    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


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
    projection = ccrs.PlateCarree
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
    for n_campaign, campaign in enumerate(campaigns):
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
    ds = ds.assign_coords({'lon': ds.lon.values - 180})
    # Update name and scaling
    Uvar2use = '{} (1E-9 {})'.format(ds[var2use].long_name, ds[var2use].units)
    ds = ds.rename({var2use: Uvar2use})
    var2use = Uvar2use
    ds = ds[[var2use]].mean(dim='time') * 1E9
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
    fig, ax = plt_highres_modelling_region(ds=ds, var2use=var2use,
                                           plot_blank_data=False,
                                           rm_colourbar=False)

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
    ds[var2use] = ds[var2use].values + ds['EMIS_DST2']
    ds[var2use] = ds[var2use].values + ds['EMIS_DST3']
    ds[var2use] = ds[var2use].values + ds['EMIS_DST4']
    ds = ds[[var2use]].mean(dim='time')
    # Remove zero data
    arr = ds[var2use].values
    arr[arr <= 0] = np.NaN
    ds[var2use].values = arr
    # Plot up the data
    fig, ax = plt_highres_modelling_region(ds=ds, var2use=var2use,
                                           plot_blank_data=False,
                                           rm_colourbar=False)

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


def mk_comparisons_of_humidty():
    """
    Make comparisons between observed and modelled humidity
    """
    # --- GEOS-Chem
    # Get model data
    dfs_mod_GC
    run2use = 'Acid-4x5-J00'
    dfs = [dfs_mod_GC[i][run2use] for i in flight_IDs]
    df = pd.concat(dfs)

    # Add satuation pressure to df
    # 𝑒𝑠0 saturation vapor pressure at 𝑇0 (Pa)
    T0 = 273.16
    T = df['GMAO_TEMP'].values  # model temp (already in kelvin)
    CC_partial_solution = np.exp((17.67 * (T - T0)) / (T - 29.65))
    df['Es'] = 611.0 * CC_partial_solution

    # NOTE: the model expoerts absolute humidty, not specific humidity
    # 𝑞  specific humidity or the mass mixing ratio of water vapor to total air (dimensionless)
    q = df['GMAO_ABSH']  # unitless (which is ≈ 𝑤 )
    p = df['GMAO_PRES']  # HPa

    # And then calculate Ws ...
    # where "𝑝 pressure (Pa)"
    df['Ws'] = 0.622 * df['Es'] / p

    # Complete calculation
    df['RH'] = 0.263 * p * q * (CC_partial_solution**-1)

    # --- GEOS-CF
    df = pd.concat([dfs_mod_CF[i] for i in flight_IDs])
    df['Alt'] = AC.hPa_to_Km(df['model-lev'].values)

    # plot
    import seaborn as sns
    sns.set(color_codes=True)
    fig, ax = plt.subplots()

    plt.title('Modelled Relative Humidity for ARNA-2 flights')
    # Plot up model data
    plt.scatter(df['RH'].values, df['Alt'].values,
                label='Relative Humidity')
#    plt.hlines( 0.753 )
#    plt.vlines(x=0.753, ymin=1000, ymax=150 )
#    ax.invert_yaxis()

    # Add a second axis
    plt.legend()
    AC.save_plot(dpi=720)
    plt.close('all')


def mk_vertical_comparisons_with_nitrate():
    """
    """

    # folder

    # Now plot by location

    #

    pass


def plt_seasonal_species_at_sites():
    """
    Plot up seasonal comparisons at campaign sites (Bermuda and CVAO)
    """
    # Get data
    RunRoot = ar.get_local_folder('RunRoot')
    FolderStr = 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.DustUptake.'
    FolderStr += 'JNIT.Isotherm.BCs.repeat.ON.II.diags.v2.J00.HourlyOutput/'
    FolderStr = RunRoot + FolderStr + 'OutputDir/'
    RunDict = {'J00': FolderStr}
    # Dates to use?
    sdate = datetime.datetime(2019, 1, 1)
    edate = datetime.datetime(2019, 12, 31)
    dates2use = pd.date_range(sdate, edate, freq='1H')

    # Put data into a dictionary
    dsD = {}
    file_str = 'GEOSChem.SpeciesConcSubset.*.nc4'
    for key in RunDict.keys():
        ds = AC.get_GEOSChem_files_as_ds(wd=RunDict[key],
                                         file_str=file_str,
                                         dates2use=dates2use)
        # Select year and surface data
        bool1 = AC.dt64_2_dt(ds.time) >= sdate
        bool2 = AC.dt64_2_dt(ds.time) <= edate
        ds = ds.isel(time=bool1)
        ds = ds.isel(time=bool2)
        ds = ds.sel(lev=ds.lev.values[0])
        # Drop excess variables and rename speices
        drop_vars = ['hyam', 'hybm', 'hyai', 'hybi', 'P0', 'AREA']
        for var in drop_vars:
            try:
                del ds[var]
            except KeyError:
                print('Not deleting: {}'.format(var))
        prefix = 'SpeciesConc_'
        VarNames = [i.split(prefix)[-1] for i in ds.data_vars]
        name_dict = dict(zip(ds.data_vars, VarNames))
        ds = ds.rename(name_dict=name_dict)
        # Save to dict
        dsD[key] = ds

    # Other runs to plot

    # Sites to plot
    # Bermuda, CVAO
    sites = ['CVO', 'BMW']
    for site in sites:

        # Sub-select data for site
        lon, lat, alt = AC.get_loc(site)

        # Plot up
        dfs = {}
        for key in RunDict.keys():

            # Subselect data
            ds2plot = ds.sel(lat=lat, lon=lon, method='nearest')
            del ds2plot['lev']
            del ds2plot['ilev']
            del ds2plot['lat']
            del ds2plot['lon']

            # Save output to csv file
            savename = 'TEMP_NetCDF_{}.nc'.format(site)
            AC.save_ds2disk_then_reload(ds2plot, savename=savename)
#            df = ds2plot.to_dataframe()
#            df.to_csv('ARNA_GEOSChem_v12_9_0_{}_{}'.format(site, key))
#            dfs[key] = df
        # Now loop to plot the species

    # Now plot
    specs2plot = ['O3', 'CO', 'HCl', 'ClNO2', 'IO', 'BrO', ]

    # plot up as diurnal
    savetitle = 'GEOSChem_v12_9_0_seasonal_diel_at_{}'.format(site)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    for spec in specs2plot:

        units, scaleby = AC.tra_unit(spec, scale=True)

        AC.plot_up_diel_by_season(dfs={'model': df.copy()*scaleby}, spec=spec,
                                  sub_str=site,
                                  units=units)

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plt_bermuda_obs(debug=False):
    """
    Plot up Bermuda observations
    """
    # Switches
    ReadFromCSV = True
    UseV13output = True

    # Get observational data
    ExcelFileName = 'Bermuda_data_hourly.xlsx'
    Folder = '/users/ts551/scratch/data/ARNA/Bermunda/'
    format = '%Y-%m-%d %H:%M:%S'
    # Spring Data
    if ReadFromCSV:
        CSVFileName = 'Bermuda_data_hourly_spring.csv'
        dfSpring = pd.read_csv(Folder+CSVFileName)
        UnitsSpring = dfSpring.loc[0, :].to_dict()
        dfSpring = dfSpring.drop(0)
        TimeVar = 'Time (utc-3)'
        dfSpring.index = pd.to_datetime(dfSpring[TimeVar].values)
    else:
        SheetName = 'Bermuda_Spring_data_hourly'
        dfSpring = pd.read_excel(Folder+FileName, SheetName=SheetName,
                                 date_parser=format)
        # Save the units as attributes
        UnitsSpring = dfSpring.loc[0, :].to_dict()
        dfSpring = dfSpring.drop([0])  # Drop the unit row
        dfSpring.index = dfSpring['Time (utc-3)']

    # Summer Data
    if ReadFromCSV:
        CSVFileName = 'Bermuda_data_hourly_summer.csv'
        dfSummer = pd.read_csv(Folder+CSVFileName)
        UnitsSummer = dfSummer.loc[0, :].to_dict()
        dfSummer = dfSummer.drop(0)
        TimeVar = 'Time'
        dfSummer.index = pd.to_datetime(dfSummer[TimeVar].values)
    else:
        SheetName = 'Bermuda_Summer_data_hourly'
        dfSummer = pd.read_excel(Folder + FileName, SheetName=SheetName,
                                 date_parser=format)
        # Save the units as attributes
        UnitsSummer = dfSummer.loc[0, :].to_dict()
        dfSummer = dfSummer.drop([0])  # Drop the unit row
        dfSummer.index = dfSummer['Time (utc-3)']

    # Store the start and end dates of the observational period
    # With a one dat buffer to improve mundging with model
    seasons = ('Spring', 'Summer', )
    dPeriods = {
        'Summer': (datetime.datetime(2019, 8, 10),
                   datetime.datetime(2019, 9, 12)
                   ),
        'Spring': (datetime.datetime(2019, 4, 16),
                   datetime.datetime(2019, 5, 14)
                   ),
    }
    # Ensure the units are the same in obs. between spring and summer
    for key in UnitsSummer.keys():
        if debug:
            print(key)
        NewUnits = UnitsSummer[key]
        if (key in list(UnitsSpring.keys())):
            CurrentUnits = UnitsSpring[key]
            SameUnits = CurrentUnits == NewUnits
            if not (SameUnits):
                PrtStr = "'Units in list for: '{}', as: '{}', ({}, same?:{})"
                print(PrtStr.format(key, NewUnits, CurrentUnits, SameUnits))
                print('Why the units different?')

    # Combine data into a single dataframe
    dfObs = pd.concat([dfSpring, dfSummer], axis=0)
    # convert times to UTC
    index = AC.dt64_2_dt(dfObs.index.values)
    dfObs.index = AC.add_hrs(index, 3)

    # Combine columns for HONO (due to inconsistent naming)
    Var1 = '[HONO]'
    Var2 = '[HONO] '
    VarHONO_Obs = 'HONO_processed'
    dfHONO = pd.concat([dfObs[Var1].dropna(), dfObs[Var2].dropna()])
    dfObs[VarHONO_Obs] = dfHONO

    # Clean up other columns and
    vars2del = Var1, Var2, 'Time (utc-3)', 'Time', 'Time.1',
    for var in vars2del:
        try:
            del dfObs[var]
        except:
            print("Error: failed to delete var: '{}'".format(var))
    # force columns to be numeric
    for col in dfObs.columns:
        dfObs[col] = pd.to_numeric(dfObs[col], errors='coerce')

    # Also add model to the comparisons
    # Use the generic year of Bermuda obs run for Pete/
    if UseV13output:
        FileName = 'v13-4-0_bermua.csv'
        dfMod = pd.read_csv(Folder + FileName)
        dfMod.index = pd.to_datetime(dfMod['datetime'].values)
        # Just use 2018 data for now.
        dfMod = dfMod.loc[dfMod.index.year == 2018, :]
        # But kludge the year to be 2019
        index = AC.dt64_2_dt(dfMod.index.values)
        dfMod.index = [AC.update_year(i, year=2019) for i in index]

    else:
        FileName = 'NSFB_ARNA_GEOSChem_v12_9_0_BMW_J00.csv'
        dfMod = pd.read_csv(Folder + FileName)
        dfMod.index = pd.to_datetime(dfMod['time'].values)

    # Add NOx to obs and model
    var1 = '[NO] (NOx system)'
    var2 = '[NO2] (NOx system)'
    dfObs['NOx'] = dfObs[var1] + dfObs[var1]
    dfMod['NOx'] = dfMod['NO'] + dfMod['NO2']
    NIT_all = ['NITD4', 'NITD3', 'NITD2', 'NITD1', 'NITs', 'NIT', ]
    dfMod['NITs-all'] = dfMod[NIT_all].sum(axis=1)

    # Map obs. species names to model ones
    dObs2Mod = {
        'Time (utc-3)': np.NaN,
        'Time': np.NaN,
        '[HNO3]': 'HNO3',
        '[pNO3] corrected': 'NITs-all',
        '[NO] (NOx system)': 'NO',
        '[NO2] (NOx system)': 'NO2',
        '[O3]': 'O3',
        'WS': 'NOT IN CURRENT DATASET',
        'WD': 'NOT IN CURRENT DATASET',
        'AirTempC_Avg': 'GMAO_TEMP',
        'RH_Avg': 'NOT IN CURRENT DATASET',
        'BP_mmHg_Avg': 'NOT IN CURRENT DATASET',
        'TSP': 'NOT IN CURRENT DATASET',
        'J(HONO)': 'NOT IN CURRENT DATASET',
        #    '[HONO] ': 'HNO2',
        # add derived variables
        'NOx': 'NOx',
        VarHONO_Obs: 'HNO2',
    }
    dObs2Mod_r = {v: k for k, v in list(dObs2Mod.items())}

    # Which vars to plot
    # TODO
    vars2plot = [
        'O3', 'NO', 'NO2', 'NOx', 'HNO3', 'NITs-all', 'HNO2',
        # Add J-rates and GMAO values
        #    'GMAO_TEMP',
    ]

    # - Plot up whole data as a generic time series
    context = 'paper'
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context(context)
    savetitle = 'ARNA_Bermunda_comp_v13'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    for season in seasons:
        # Sub-select dates for season
        sdate, edate = dPeriods[season]

        for var in vars2plot:
            print(season, var)
            units, scaleby = AC.tra_unit(var, scale=True)
            # Force units to be the same as observations in pptv
            pptv_units = 'HNO3', 'HNO2', 'NO', 'NO2', 'NOx'
            if (var in pptv_units):
                units, scaleby = 'pptv', 1E12
            # V13 alreaedy scaled
            if UseV13output:
                scaleby = 1
            # Plot obs
            bool1 = dfObs.index >= sdate
            bool2 = dfObs.index <= edate
            df2plot = dfObs.loc[(bool1 & bool2), :]
            plt.plot(df2plot[dObs2Mod_r[var]].index,
                     df2plot[dObs2Mod_r[var]].values,
                     label='Obs.',
                     color='Black',
                     )

            # Try to plot model
            bool1 = dfMod.index >= sdate
            bool2 = dfMod.index <= edate
            df2plot = dfMod.loc[(bool1 & bool2), :]
            plt.plot(df2plot[var].index,
                     df2plot[var].values * scaleby,
                     label='Model',
                     color='Red',
                     )

            plt.legend()

            # Add a title
            PrtStr = "Timeseries {} @ Bermuda during *{}* campaign ({})"
            plt.title(PrtStr.format(var, season.lower(), units))

            # Update x axis label rotation
            plt.xticks(rotation=45)

            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()

        # Plot up daily cycles
        for var in vars2plot:
            fig, ax = plt.subplots()
            print(season, var)
            units, scaleby = AC.tra_unit(var, scale=True)
            # Force units to be the same as
            pptv_units = 'HNO3', 'HNO2', 'NO', 'NO2', 'NOx'
            if (var in pptv_units):
                units, scaleby = 'pptv', 1E12

            # Plot obs
            bool1 = dfObs.index >= sdate
            bool2 = dfObs.index <= edate
            df2plot = dfObs.loc[(bool1 & bool2), :]

            AC.BASIC_diel_plot(dates=df2plot[dObs2Mod_r[var]].index,
                               data=df2plot[dObs2Mod_r[var]].values,
                               label='Obs.',
                               color='Black',
                               spec=var, units=units,
                               fig=fig, ax=ax)

            # Try to plot model
            bool1 = dfMod.index >= sdate
            bool2 = dfMod.index <= edate
            df2plot = dfMod.loc[(bool1 & bool2), :]

            AC.BASIC_diel_plot(dates=df2plot[var].index,
                               data=df2plot[var].values * scaleby,
                               label='Obs.',
                               color='Red',
                               spec=var, units=units,
                               fig=fig, ax=ax)

            # Add a title
            PrtStr = "Diel cycle {} @ Bermuda during {} campaign ({})"
            plt.title(PrtStr.format(var, season.lower(), units))

            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()

    # Save PDF
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plt_comp_with_NASA_Atom():
    """
    Plot up comparisons with NASA ATom data

    NOTE: Analysis now in Prj_ATom_analysis.py
    """
    # Get lastest NASA ATom data and plot up

    #

    #

    pass



def plt_comp_with_FIREX_AQ(debug=False):
    """
    Plot up a comparison of of various species
    """
    # Retrieve observational data
#     dfObs = ar.get_FIREX_AQ_data(RtnAllData=True,
#                                  FilterByTimeOfDay=False,
#                                  UpdateTimeeZone2LocalTime=False,)
    dfObs = ar.get_FIREX_AQ_data(RtnAllData=True,
                                 FilterByTimeOfDay=False,
                                 UpdateTimeeZone2LocalTime=False,
                                 FilterPollutedAirMasses=True,
                                 RmFlaggedData=True,
                                 )

    LatVar = 'Latitude_YANG'
    LonVar = 'Longitude_YANG'
    AltVar = 'MSL_GPS_Altitude_YANG'
    TimeVar = 'datetime'
    dfObs[TimeVar] = dfObs.index.values
    # Add in pressure and time columns to the DataFrame
    PressVar = 'hPa'
    dfObs[PressVar] = dfObs[AltVar] / 1E3
    dfObs[PressVar] = dfObs[PressVar].map(AC.km2hPa)
    dfObs_BACKUP = dfObs.copy()
    # Check the comparability of the altitude variables
    if debug:
        vars2plot = [ PressVar, 'Pressure_Altitude_YANG']
        for __var in vars2plot:
            plt.plot( dfObs[AltVar], dfObs[__var])
            plt.xlabel(AltVar)
            plt.ylabel(__var)
            #
            SaveName = 'ARNA_FIREXAQ_Alt_Press_dim_{}_vs_{}'
            AC.save_plot(SaveName.format(AltVar, __var))
            plt.close('all')

    # which model runs to use?
    if isinstance(RunDict, type(None)):
        RunSet = 'IGAC.ARNAv14'
        res = '4x5'
        GC_version='v12.9'
        RunDict = ar.get_dict_of_GEOSChem_model_output(res=res, RunSet=RunSet,
                                                       GC_version=GC_version,
                                                       folder4netCDF=True)

    # Extract nearest point in observations
    DataRoot = '/users/ts551/scratch/data/ARNA/FIREX_AQ/'
    dfsMod = {}
    keys2use = list(sorted(RunDict.keys()))
#    keys2use =  'Andersen22b.TEMPII' ,  'J00'
    for key in keys2use:
        folder = RunDict[key]

        # Use the backed up version of dfObs
        dfObs = dfObs_BACKUP.copy()

        # Extract on a per day basis?
        # This will reduce overheads as fewer files need to be opened.
        # This isn't possible due to the approach used of merging flights
        # without retaining the Reserach Flight (RF) number
        # Instead select the broad area of the observations and the specific
        # dates to use. Then save this locally and then re-extract from this.
        # Or remerge the files and included a FileNumber variable.

        #
        path4GC = '/{}/GeosChem_Extracted/'.format(DataRoot)
        #
        FileNumbers = list(set(dfObs['FileNumber'].values))
        FileNumbers = list(sorted([i for i in FileNumbers if np.isfinite(i) ]))
        FileNumbers = [int(i) for i in FileNumbers]
#        FileNumbers = np.arange(1, 20)
        dfModAll = pd.DataFrame()
        for FileNumber in FileNumbers:
            # Select the values to extract
            __bool = dfObs['FileNumber'] == FileNumber
            df2extract = dfObs.loc[__bool,:].copy()
            ObsShape4FileNumber =  df2extract.shape
            print(FileNumber, ObsShape4FileNumber, key)

            # Try and find the values saved offline
            csv_savenameStr = 'FIREX_AQ_Extracted_FileNumber_{:0>3}_{}.csv'
            csv_savename = csv_savenameStr.format(FileNumber, key)
            try:
                dfMod = pd.read_csv(path4GC+csv_savename)
                # Apply the offline saved index
                dfMod.index = dfMod[ 'SavedIndex' ]
                PrtStr = "{} FileNumber Extracted ('{}'). Same shape? {} ({})"
                SameNumRows = dfMod.shape[0] == ObsShape4FileNumber[0]
                print(PrtStr.format(FileNumber, key, dfMod.shape, SameNumRows))
                PrtStr = 'NOTE: Using saved csv file already extracted: {}'
                print(PrtStr.format(csv_savename))
#                del dfMod.index.name
            except FileNotFoundError:
                Pstr = 'WARNING: data being extracted for FileNumber: {}'
                print(Pstr.format(FileNumber))
                # Now extract the data for the specific filename (research flight)
                dfMod = Get_GEOSChem4flighttracks(df=df2extract,
                                                  LonVar=LonVar,
                                                  LatVar=LatVar,
                                                  PressVar=PressVar,
                                                  TimeVar=TimeVar,
                                                  TempSaveDir=path4GC,
                                                  folder=folder)
                # Print a check that the same size data has been extracted
                PrtStr = "{} FileNumber Extracted ('{}'). Same shape? {} ({})"
                SameNumRows = dfMod.shape[0] == ObsShape4FileNumber[0]
                print(PrtStr.format(FileNumber, key, dfMod.shape, SameNumRows))
                # Save the file to csv
                save2csv=True
                if save2csv:
                    dfMod['SavedIndex'] = dfMod.index.values
                    dfMod.to_csv(path4GC+csv_savename)

            #
            dfModAll = pd.concat([dfModAll, dfMod])

            #
            dfsMod[key] = dfModAll

    # Plot up vertical plots etc to compare on a per species basis

    # Which variables to plot?
    AltVar = 'MSL_GPS_Altitude_YANG'
    HNO2_CIMS = 'HNO2_NOAACIMS_VERES'
#    HNO2_ACES = 'HNO2_ACES_WOMACK'
#    HONO_SAGA = 'HONO_SAGA_DIBB'
    # Shared plotting variables
    SCprefix = 'SpeciesConc_'
    xlabel = 'HONO (pptv)'
    bins = np.arange(1, 11)
#    num_of_datasets = 3
    num_of_datasets = len(dfsMod.keys()) + 1

    # --- Quick plot of the data binned into
    vars2plot = ['HNO2', ]
#    [HNO2_CIMS, ]
    for var2plot in vars2plot:
        df2plot = dfObs
#        var2plot = HONO_SAGA
        # Select observation variable
        ObsVar = HNO2_CIMS
        color = 'k'
#        label = 'HNO2 ({})'.format(ObsVar)
        label = 'NOAA CIMS'
        # drop NaNs / flagged data
        FlagValue = -999999.000000
        df2plot.loc[ df2plot[ObsVar] == FlagValue, : ] = np.NaN
        df2plot.loc[ df2plot[ObsVar] < 0.0, : ] = np.NaN

        df2plot = df2plot[ [ObsVar, AltVar] ].dropna()
        # Update metres to kilometres
        df2plot[ AltVar ] = df2plot[ AltVar ].values / 1E3

        # Setup plot and title
        fig, ax = plt.subplots()
        title = 'FIREX_AQ_{}_quick_plot_v9_extra'.format(var2plot)

        AC.binned_boxplots_by_altitude(df=df2plot, fig=fig, ax=ax,
                                       var2bin_by=AltVar,
                                       label=label, xlabel=xlabel,
                                       binned_var=ObsVar,
                                       num_of_datasets=num_of_datasets,
                                       bins=bins,
                                       widths=0.15,
                                       dataset_num=1,
                                       color=color)

        # Extra update of labels
        label_dict = {
        'Andersen22b.TEMPII': 'Andersen22b',
        'Andersen22b': 'Andersen22b',
        'Shah22': 'Shah22',
        'J00': 'Base',
        'Ye17': 'Ye17',
        'Kas18':'Kas18',
        }

        # Add model values to the plot
#        keys2use = list(sorted(dfsMod.keys()))
        keys2use = [ 'Andersen22b', 'J00' ]
        colors2use = AC.get_CB_color_cycle()
        color_dict = dict(zip(keys2use, colors2use))
        for nkey, key in enumerate( keys2use ):
            df2plot = dfsMod[key]
            ModVar = '{}{}'.format(SCprefix, var2plot)
            # Use the observed altitude that was extracted for now
            df2plot[ AltVar ] = dfObs[AltVar].copy()
            # Update metres to kilometres
            df2plot[ AltVar ] = df2plot[ AltVar ].values / 1E3

            # drop NaNs / flagged data
            df2plot = df2plot[ [ModVar, AltVar] ].dropna()

            # Update units for HNO2
            df2plot.loc[:, ModVar] = df2plot.loc[:, ModVar] * 1E12

            label = '{}'.format(label_dict[key])
            #
            AC.binned_boxplots_by_altitude(df=df2plot, fig=fig, ax=ax,
                                           var2bin_by=AltVar,
                                           label=label, xlabel=xlabel,
                                           binned_var=ModVar,
                                           num_of_datasets=num_of_datasets,
                                           bins=bins,
                                           widths=0.15,
                                           dataset_num=1+(nkey+1),
                                           color=color_dict[key])

        # Beautify plot
        plt.legend()

        # Save the plot
        AC.save_plot(title)
        plt.close('all')


def add_dsLev_idx_from_dsLev_value(df=None, folder=None,
                                   LevVar='ds-lev',
                                   LevIdxVar='ds-lev-idx'):
    """
    Add the vertical index from the ds-lev value
    """


    #



    # Get a value for offline pressure
    OfflinePress = AC.gchemgrid('c_hPa_geos5_bounds')







def Get_GEOSChem4flighttracks(df, folder=None, TempSaveDir='./',
                              LatVar='lat', LonVar='lon',
                              TimeVar='datetime',
                              PressVar='hPa'):
    """
    Extract 3D data to retrieve modelled flight tracks for obs. locations
    """
    # Ext the extents of the observational data
    LatMin = df[LatVar].values.min()
    LatMax = df[LatVar].values.max()
    LonMin = df[LonVar].values.min()
    LonMax = df[LonVar].values.max()
    # Select dates to use to extract hourly 3D output from
    def __convert_datetime2days(input):
        return datetime.datetime(*input.timetuple()[:3])
    dates2use = AC.dt64_2_dt(df[TimeVar].values)
    dates2use = [__convert_datetime2days( i) for i in dates2use ]
    dates2use = list(set(dates2use))
    Months = [i.month for i in dates2use]
    # Which dates to use for pressure values
    def __convert_datetime2months(input):
        """ Convert datetimes to Year, Month, 1st of Day of Month """
        return datetime.datetime(*list(input.timetuple()[:2])+[1])
    StateMetDates2Use = [__convert_datetime2months(i) for i in dates2use]
    StateMetDates2Use = list(set(StateMetDates2Use))

    # Get 3D data from model (SpeciesConcSubset)
    file_str = 'GEOSChem.SpeciesConcSubset.*.nc4'
    ds = AC.get_GEOSChem_files_as_ds(wd=folder, file_str=file_str,
                                     dates2use=dates2use)
    # Get stateMet and add to core dataset (NOTE: dimensions  must be constant)
    StateMetfolder = folder
#     print('WARNING: Using temporary folder for StateMet!')
#     StateMetfolder ='/users/ts551/scratch/GC/rundirs/P_ARNA/geosfp_4x5_aciduptake.v12.9.0.ARNA.Isotherm.Diags.v9.Base/OutputDir/'
    StateMet = AC.get_StateMet_ds(wd=StateMetfolder,
                                  dates2use=StateMetDates2Use)
    StateMetPress = 'Met_PMID'
    CopyVar = 'SpeciesConc_HNO2'
    ds[StateMetPress] = ds[CopyVar].mean(dim='time')
    ds[StateMetPress] = StateMet[StateMetPress].copy().mean(dim='time')

    # Reduce the file size and temporary
    bool1 = ((ds.lon >= LonMin) & (ds.lon <= LonMax)).values
    bool2 = ((ds.lat >= LatMin) & (ds.lat <= LatMax)).values
    # If not a single lat or lon between extents, use closest
    if (all(bool1) == False):
        idx = AC.find_nearest(ds.lon, LonMin )
        bool1[idx] = True
    if (all(bool2) == False):
        idx = AC.find_nearest(ds.lat, LatMin )
        bool2[idx] = True

    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    # Save reduced NetCDF to disk, then reload
    savename = 'TEMP_NetCDF_ARNA_VIII.nc'
    ds = AC.save_ds2disk_then_reload(ds, folder=TempSaveDir, savename=savename)

    # delete variables that are not coordinates


    # Calculate pressure
#    PS = 1013. / 10 # in Pa (so hPa/10)
    CalculateOnlinePressCoord = False
    if CalculateOnlinePressCoord:
        PS = ds['Met_PS1DRY']
        pressure = ds['hyam'] * ds['P0'] + ds['hybm'] * PS
    else:
        # Temp. Use the offline pressure coordinate
        # Calculate the press idx offline here
        # Use the average pressure vertical for the region of analysis
        ds_hPa = ds[StateMetPress].mean(dim=['lat', 'lon']).values

        # Then calculate the indexs
        d = AC.calc_4D_idx_in_ds(ds=ds, df=df,
                                 LonVar=LonVar,
                                 LatVar=LatVar,
                                 TimeVar=TimeVar,
                                 AltVar=PressVar,
                                 ds_hPa=ds_hPa,
                                 )
        # Then extract 4D netCDF en masse
        dfN = pd.DataFrame()
        vars2extract = list(ds.data_vars)
        SCprefix = 'SpeciesConc_'
        vars2extract = [ i for i in vars2extract if SCprefix in i ]
        times2use = df.index.values
        for n, time in enumerate(times2use):
            # get the times for a specific data
            lat_idx = d[LatVar][n]
            lon_idx = d[LonVar][n]
#            lev_idx = d[AltVar][n]
            lev_idx = d[PressVar][n]
            time_idx = d[TimeVar][n]
            # en masse extract indexes
            ds_tmp = ds.isel(lat=lat_idx, lon=lon_idx, time=time_idx,
                             lev=lev_idx)
            vals = [ds_tmp[i].data for i in vars2extract]
            vals = np.array(vals)
            for nval, val in enumerate(vals):
                dfN.loc[vars2extract[nval], time] = vals[nval]
            # Add the model position coordinates...
            dfN.loc['ds-lat', time] = float(ds_tmp['lat'].values)
            dfN.loc['ds-lon', time] = float(ds_tmp['lon'].values)
            dfN.loc['ds-lev', time] = float(ds_tmp['lev'].values)
            dfN.loc['ds-time', time] = float(ds_tmp['time'].values)
            del ds_tmp, vals
        # Make datetime the index
        dfN = dfN.transpose()
        # Save the datetime as a column too
        dfN['Datetime'] = dfN.index.values
        # Update the model datetime to be in datetime units
        dfN['ds-time'] = pd.to_datetime(dfN['ds-time'].values)

    # Delete the temporary NetCDF
    os.remove(TempSaveDir+savename)

    # Extract the nearest locations to the observations
    return dfN


def plt_obs_based_FIREX_analysis():
    """
    Plot a quick HNO2 plot with all campaign data
    """
    # Retrieve observational data
    dfObs = ar.get_FIREX_AQ_data(RtnAllData=True,
                                 FilterByTimeOfDay=False,
                                 UpdateTimeeZone2LocalTime=False,
                                 FilterPollutedAirMasses=True,
                                 RmFlaggedData=True,
                                 )

    # Which variables to plot?
    AltVar = 'MSL_GPS_Altitude_YANG'
    HNO2_CIMS = 'HNO2_NOAACIMS_VERES'
    HNO2_ACES = 'HNO2_ACES_WOMACK'
    HONO_SAGA = 'HONO_SAGA_DIBB'
    # Shared plotting variables
    xlabel = 'HONO (pptv)'
    color = 'k'
    bins = np.arange(1, 11)

    # --- Quick plot of the data binned into
    for var2plot in [HNO2_CIMS, HNO2_ACES, HONO_SAGA]:
        df2plot = dfObs
#        var2plot = HONO_SAGA
        label = 'HNO2 ({})'.format(var2plot)
        # drop NaNs / flagged data
        FlagValue = -999999.000000
        df2plot.loc[ df2plot[var2plot] == FlagValue, : ] = np.NaN
#        df2plot.loc[ df2plot[var2plot] < 0.0, : ] = np.NaN

        df2plot = df2plot[ [var2plot, AltVar] ].dropna()
        # Update metres to kilometres
        df2plot[ AltVar ] = df2plot[ AltVar ].values / 1E3

        # Setup plot and title
        title = 'FIREX_AQ_HNO2_quick_plot_{}'.format(var2plot)
        fig, ax = plt.subplots()

        AC.binned_boxplots_by_altitude(df=df2plot, fig=fig, ax=ax,
                                       var2bin_by=AltVar,
                                       label=label, xlabel=xlabel,
                                       binned_var=var2plot,
                                       num_of_datasets=1,
                                       bins=bins,
                                       widths=0.15,
                                       dataset_num=1,
                                       color=color)
        plt.legend()
        AC.save_plot(title)




def testing_v13_implimentation():
    """
    Do
    """
    # Local settings
    RunSet = 'v13.4.1month'
    res = '4x5'
    GC_version = 'v13.4'
    #
    RunDict = get_dict_of_GEOSChem_model_output(RunSet=RunSet, res=res,
                                                GC_version=GC_version,
                                                folder4netCDF=True)

    # Dates to use for the test run? - assume 1 month output
    if isinstance(dates2use, type(None)):
        sdate = datetime.datetime(2018, 6, 1) # Beginning for spin up (6months)
        edate = datetime.datetime(2018, 6, 30) # 3 months into analysis year
        dates2use = pd.date_range(sdate, edate, freq='1D')

    # Get basic stats
    ExtraSpecs = ['HNO3', 'HNO2', 'NIT', 'NITs',
#                  'NITD1', 'NITD2', 'NITD3', 'NITD4'
                  ]
    df = AC.get_stats4RunDict_as_df(RunDict=RunDict,
                                        extra_burden_specs=ExtraSpecs,
                                        extra_surface_specs=ExtraSpecs,
                                        dates2use=dates2use,
                                        REF_wd=RunDict['4pptHONO'],
                                        use_REF_wd4Met=True )

    # Check the Jscale values and related checks
    # Use exiting in do_post_campaign_NOx_analysis.py
#    explore_JVALS_in_JNITs_runs()

    #


def Add_NITD_to_restart_file(folder='./',
                             FileName='GEOSChem.Restart.20180701_0000z.nc4'):
    """
    Manually add dust uptake variables to the restart file

    Notes
    ---
     - The shipped restart file doens't have NITD in it. So add this by copying
     another file
    """
    # Template NetCDF with values to transfer to restart file
    RunRoot = '/mnt/lustre/users/ts551/GC/rundirs/P_ARNA/'
    folder = 'geosfp_4x5_aciduptake.v12.9.0.ARNA.Isotherm.Diags.v9.Base'
    FileName = 'GEOSChem.SpeciesConc.20190601_0000z.nc4'
    PathStr = '{}/{}/OutputDir/{}'
    dsT = xr.open_dataset(PathStr.format(RunRoot, folder, FileName) )
    SCprefix = 'SpeciesConc_'

    # Restart file to update values in
    folder = 'gc_4x5_47L_merra2_fullchem_aciduptake.v13.4.1.ARNA.Test/'
    FileName = 'GEOSChem.Restart.20180601_0000z.nc4.BACKUP_UPDATED_DATES'
    dsR = xr.open_dataset('{}/{}/{}'.format(RunRoot, folder, FileName) )
    prefix = 'SpeciesRst_'
    dsRdata_varsList = dsR.data_vars
    dsRdata_varsList = [i for i in dsRdata_varsList if prefix in i ]
    dsR_SpeciesList = [i.split(prefix)[-1] for i in ds2data_varsList ]

    # Update Restart value to "well spin up" values
    TemplateVar = 'SpeciesRst_NIT'
    Pstr = "NOTE: Updating species ('{}') in dsR with values from template ds"
    for data_var in dsT.data_vars:
        if SCprefix in data_var:
            Species = data_var.split(SCprefix)[-1]

            if Species not in dsR_SpeciesList:
                print(Pstr.format(Species))

                # Setup new species
                NewVar = '{}{}'.format(prefix, Species)
                dsR[NewVar] = dsR[TemplateVar].copy()
                # Copy across attributes within restart file and update
                attrs = dsR[TemplateVar].attrs.copy()
                long_name = 'Dry mixing ratio of species {}'.format(Species)
                attrs['long_name'] = long_name
                dsR[NewVar].attrs = attrs

                # Update the values in the restart object with the spun up ones
                SpunUpVar = '{}{}'.format(SCprefix, Species)
                dsR[NewVar].values = dsT[SpunUpVar].values

    # Save to netCDF
    FileName = 'GEOSChem.Restart.20180601_0000z.nc4.UPDATED'
    dsR.to_netcdf('{}/{}/{}'.format(RunRoot, folder, FileName) )


def check_values_in_GC_Restart_file(ds=None,  prefix='SpeciesRst_',
                                    vars2plot=[], dpi=320):
    """
    Check values in GEOS-Chem restart file
    """
    # Have a check
    if isinstance(ds, type(None)):
        FileName = 'GEOSChem.Restart.20180601_0000z.nc4'
        folder = 'gc_4x5_47L_merra2_fullchem_aciduptake.v13.4.1.ARNA.Test/'
        FileName = 'GEOSChem.Restart.20180601_0000z.nc4.BACKUP_UPDATED_DATES'
        ds = xr.open_dataset('{}/{}/{}'.format(RunRoot, folder, FileName) )

    # Which variables to check(var2check)?
    # Which species are these (species2check)?
    if isinstance(vars2check, type(None)):
        vars2check = ds.data_vars
        vars2check = [i for i in vars2check if prefix in i ]
        species2check = [i.split(prefix)[-1] for i in vars2check ]

    # Check if max values are within 4 orders of magnitude or background values
    # Use the yaml file for this.
    dfG = pd.DataFrame() # Global values
    dfS = pd.DataFrame() # Surface values
    # loop first and check species values
    for n, __var in enumerate( vars2check ):
        Species = species2check[n]
        print(n, __var, Species)

        # Get stats on surface concentrations
        data = ds[__var].isel(lev=(ds.lev==ds.lev[0]))
        dfS[Species] = pd.Series(data.values.flatten()).describe()

        # Get Global statistics
        data = ds[__var]
        dfG[Species] = pd.Series(data.values.flatten()).describe()

    SaveName = 'ARNA_Restartfile_check_global_values.csv'
    dfG.to_csv(SaveName)
    SaveName = 'ARNA_Restartfile_check_surface_values.csv'
    dfS.to_csv(SaveName)


    # Check surface concentrations against background values?
    # Check for large differences between surface and global (3D) values?

    # -
    # Species to plot?
    if len(vars2plot) == 0:
        vars2plot = ['O3', 'NO2', 'NO', 'CO', 'NIT', 'NITs']
        vars2plot += ['SO4D{}'.format(i) for i in range(1,5)]
        vars2plot += ['NITD{}'.format(i) for i in range(1,5)]
        vars2plot += ['DSTAL{}'.format(i) for i in range(1,5)]
        vars2plot += ['DST{}'.format(i) for i in range(1,5)]
    # Setup PDF to save values too
    savetitle = 'RestartFile_plotted_value_check_v13_4_1'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # Loop and plot requested species
    ds2plot = ds.copy().mean(dim='time').isel(lev=(ds.lev==ds.lev[0]))
    ds2plot = ds2plot.squeeze()
    for __var in vars2plot:
        var2plot = '{}{}'.format(prefix, __var)

        # Apply any unit conversions etc if possible

        # plot surface
        kwargs = {}
        AC.quick_map_plot(ds2plot, var2plot=var2plot, verbose=verbose,
                          **kwargs)

        TitleStr = "Surface [{}] in Restart file"
        plt.title(TitleStr.format(var2plot))

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

    # plot zonal
    ds2plot = ds.copy().mean(dim='time')
#    SMfolder = '/mnt/lustre/users/ts551/GC/rundirs/P_ARNA/'
#    SMfolder += 'gc_4x5_47L_merra2_fullchem.v13.4.1.ARNA.Test/OutputDir/'
#    StateMet = AC.get_StateMet_ds(wd=SMfolder, dates2use=dates2use)
#    StateMet = StateMet.mean(dim='time')
    for __var in vars2plot:
        var2plot = '{}{}'.format(prefix, __var)

        fig, ax = plt.subplots()
        im = AC.ds2zonal_plot(ds2plot, var2plot=var2plot, StateMet=StateMet,
                              fig=fig, ax=ax, **kwargs)
        TitleStr = "Zonal [{}] in Restart file"
        plt.title(TitleStr.format(var2plot))

        # Add a colourbar
#        kwargs = {'extend':'both'}
        fig.colorbar(im, orientation="horizontal", pad=0.2, extend='both',
                     **kwargs)
#                     format=format, label=units)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()


    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def check_Ye17_output_log_file_values():
    """
    Check the range of JScale values
    """

    # Location and name of log file to check
    folder = '/users/ts551/scratch/GC/rundirs/P_ARNA/'
    folder += 'gc_4x5_47L_merra2_fullchem_aciduptake.v13.4.1.ARNAv11/'
    filename = 'GC.geos.log'

    # Compile data from prefixed lines in log file.
#    DebugStr = '@@@ Ye Jscale + NO3_CONC III:'
#    DebugStr = '@@@ Applied Jscale (for NO3_CONC):'
    DebugStr = '@@@ Used Jscale (+ [NO3]):'
#    DebugStr = '@@@ NO3_CONC III (nM/m^3) - NIT(s)+NITD1-4'
    with open(folder+filename, 'r') as OpenedFile:
        lines = [i for i in OpenedFile if DebugStr in i ]
#    lines = [i for i in lines if (DebugStr in i)]
    data = [i.split(DebugStr)[-1].strip().split() for i in lines]
    Jscale = np.array(data)[:, 0].astype(np.float)
    NO3_CONC = np.array(data)[:, 1].astype(np.float)
#    NO3_CONC = np.array(data)[:, 0].astype(np.float)
    df = pd.DataFrame()
    df['NO3_CONC'] = NO3_CONC

    #
    df['Jscale'] = Jscale

    df.describe()


    # - plot up
    # Convert nitrate from nmoles/m3 => ug/M3 (hash lines out if not)
    def ug_m3_2nmoles_m3(x):
        return x / 1e-9 / (14.+16.+16.+16.) / 1e6
    def nmoles_m3_2ug_m3(x):
        return x * 1e-9 * (14.+16.+16.+16.) * 1e6


    fig, ax1 = plt.subplots(dpi=320)
#    no3 = np.arange(0.01, 10000, 0.01) # nmoles/m3

    plt.scatter( NO3_CONC, Jscale, )
    units = 'nmoles m$^{-3}$'
    plt.xlabel('bulk [NO$_{3}^{-}$]'+' ({})'.format(units))
    plt.ylabel('JScale')

    # Set scaling to loglog
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    TitleSuffix = '2CPU'
    TitleStr = 'Jscale calculated online with Ye17 param on {}'
    plt.title( TitleStr.format(TitleSuffix) )
    title = 'ARNA_online_Ye17_check_{}'.format(TitleSuffix)
    AC.save_plot(title=title)


if __name__ == "__main__":
    main()

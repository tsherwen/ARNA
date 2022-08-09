"""
Driver script to do NOx/HONO (post campaign) analysis
"""
#!/usr/bin/python
# - Packages
import numpy as np
from time import gmtime, strftime
import time
import glob
import AC_tools as AC
import sys
import pandas as pd
# import ARNA analysis/campaign code as a module
import arna as ar
from matplotlib.colors import LogNorm


def main():
    """
    Main driver function
    """
    # --- NOx Budget
    # Perform core NOx budget analysis
    plt_key_NOx_budget_terms()
    analyse_NOx_budget()

    # Plot up the NOx budget against latitude
    plt_NOx_routes_by_lat()

    #

    # --- Misc NOx work during campaign
    # Get the core FAAM data
#    ar.get_FAAM_core4flightnum()

    # Get the ToF-CIMS data
    ds2 = ar.get_CIMS_data4flight()

    # --- Genreal run testing and acid uptake run exploration
    # Explore NOy/NIT after inclusion of acid uptake
    explore_NOy_with_acid_uptake()
    explore_ARNA_period_with_acid_uptake()

    # Explore lightning sources within ARNA runs
#    plt_lightning_by_month()

    # Explore the output of Jvalues and PF code
    explore_JVALS_in_JNITs_runs()

    # Do debbuging of Isotherm runs
#    explore_print_statement_debug_of_aciduptake_photolysis()


def explore_NOy_with_acid_uptake():
    """
    Explore the global NOy budgets with acidic uptake to dust included
    """
    from AC_tools import species_mass
    # - Local variables
    # use a v12.9.1 compatible definition of NOx
#    NOySpecs = AC.GC_var('NOy')
    NOySpecs = ['NO', 'NO2', 'PAN', 'HNO3', 'HNO2', 'NOy', 'NOy-gas'
                #    'PMN', 'PPN', 'R4N2',
                #        'N2O5', 'HNO4',\
                #        'BrNO2', 'BrNO3',
                #    'MPN',
                #    'ISOPN', 'PROPNN',
                # 'MMN',\
                #        'NO3', 'HNO2', 'IONO', 'IONO2', 'INO', 'ClNO2', 'ClNO3'
                ]
    # include NITs
    NOySpecs += ['NIT', 'NITs']
    # Also consider dust uptake on NOx
    NOySpecsA = NOySpecs + ['NITD1', 'NITD2', 'NITD3', 'NITD4']
    # Set runs to use
    RunSet = 'ACID'
    res = '4x5'
    folder4netCDF = True
    RunDict = ar.get_dict_of_GEOSChem_model_output(res=res, RunSet=RunSet,
                                                   folder4netCDF=folder4netCDF)
#     RunRoot = '/users/ts551/scratch/GC/rundirs/'
#     RunStr = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020{}/'
#     run_dict = {
#         'BASE': RunRoot+RunStr.format('.BCs.repeat.III/spin_up/'),
#         #    'ACID.III': RunRoot + RunStr.format('.DustUptake.III/spin_up/'),
#         'ACID.IV': RunRoot+RunStr.format('.DustUptake.IV/OutputDir/'),
#         'JNIT': RunRoot+RunStr.format('.DustUptake.IV.JNIT/OutputDir/'),
#         'JNITx25': RunRoot+RunStr.format('.DustUptake.IV.JNIT.x25/OutputDir/'),
#     }
#    run_dict = d
    # - Analysis
    # Get generic stats on runs
#    dates2use = None
    dates2use = [datetime.datetime(2019, 1+i, 1) for i in range(12)]

    extra_specs = NOySpecs + ['SO2', 'SO4']
    df = AC.get_general_stats4run_dict_as_df(run_dict=RunDict,
                                             dates2use=dates2use,
                                             extra_burden_specs=extra_specs)

    # Just get NOy species in N Tg equilivents
    avg_over_time = True  # Note: burdens area averaged overtime
    prefix = 'SpeciesConc_'
#    dates2use = None
    dates2use = [datetime.datetime(2018, 3+i, 1) for i in range(6)]
    ref_spec = 'N'
    rm_strat = True
    use_time_in_trop = True
    # Mass unit scaling
    mass_scale = 1E3
    mass_unit = 'Tg'

    # Get all of the speciesConcs for runs as list of datasets
    dsD = {}
    for key in run_dict.keys():
        dsD[key] = AC.GetSpeciesConcDataset(wd=run_dict[key],
                                            dates2use=dates2use)
    # Core dataframe for storing calculated stats on runs
    df = pd.DataFrame()
    #
    keys2use = ['ACID', 'BASE']
    keys2use = ['BASE', 'ACID.III', 'ACID.IV', ][::-1]
    keys2use = ['BASE', 'ACID.IV', 'JNIT', 'JNITx25'][::-1]
    for key in keys2use:
        print(key, run_dict[key])
        # Use the appropriate list of species
        if ('ACID' in key) or ('NIT' in key):
            specs2use = NOySpecsA
        else:
            specs2use = NOySpecs
        vars2use = [prefix+i for i in specs2use]
        # Get StateMet object
        StateMet = AC.get_StateMet_ds(wd=run_dict[key], dates2use=dates2use)
        # Average burden over time
        ds = dsD[key].copy()  # .mean(dim='time', keep_attrs=True)

        # reduce datasets to region of analysis
        region = 'model_nest_area'
#        region = None
        if not isinstance(region, type(None)):
            ds_l = [StateMet, ds]
            for _n, _ds in enumerate(ds_l):
                #
                # reduce lat and lon
                d = ar.get_analysis_region(region)
                x0, x1, y0, y1 = d['x0'], d['x1'], d['y0'], d['y1']
                # Set values region
                bool1 = ((_ds.lon >= x0) & (_ds.lon <= x1)).values
                bool2 = ((_ds.lat >= y0) & (_ds.lat <= y1)).values
                # Cut by lon, then lat
                _ds = _ds.isel(lon=bool1)
                _ds = _ds.isel(lat=bool2)
                # reduce  alt to bottom ~8km (350 hPa)
                HPa_l = AC.get_GEOSCF_vertical_levels(native_levels=True)[::-1]
                bool = [i > 350 for i in HPa_l]
                _ds = _ds.isel(lev=bool)
                #
                ds_l[_n] = _ds
            StateMet, ds = ds_l

        S = AC.get_Gg_trop_burden(ds, vars2use=vars2use, StateMet=StateMet,
                                  use_time_in_trop=use_time_in_trop,
                                  avg_over_time=avg_over_time,
                                  rm_strat=rm_strat)
        # convert to ref spec equivalent (e.g. N for NO2, C for ACET)
        for spec in specs2use:
            #                ref_spec = get_ref_spec(spec)
            val = S[prefix+spec]
            S[prefix+spec] = val/species_mass(spec)*species_mass(ref_spec)
        # Upate varnames
        varnames = ['{} burden ({})'.format(i, mass_unit) for i in specs2use]
        S = S.rename(index=dict(zip(list(S.index.values), varnames)))
        # Save the values for run to central DataFrame
        df[key] = S

    # Sum up NOy
    df = df.T
    df.T.sum(axis=0)
    df['Sum'] = df.T.sum(axis=0)
    # Add sum of NIT(S) burden
    cols2use = [i for i in df.columns if 'NIT' in i]
    df['Sum (NITS-All)'] = df[cols2use].T.sum(axis=0)
    # Add Sum of NIT(S) burden-NITD
    cols2use = [i for i in cols2use if 'NIT ' not in i]
    cols2use = [i for i in cols2use if 'NITs ' not in i]
    df['Sum (NITS-D)'] = df[cols2use].T.sum(axis=0)
    df = df.T
    print(df[keys2use])

    # Calc percent difference
    df_pcent = df.copy()
    REF = 'BASE'
    for col in [i for i in df_pcent.columns if i != REF]:
        df_pcent[col] = (df_pcent[col] - df_pcent[REF]) / df_pcent[REF]*100
    print(df_pcent[keys2use])

    # Calc percent of each species/family of total NOy.
    df_pcent_NOy = df.copy()
    REF = 'BASE'
    for col in [i for i in df_pcent_NOy.columns if i != REF]:
        df_pcent_NOy[col] = df_pcent_NOy[col] / df_pcent_NOy[col]['Sum']*100
    print(df_pcent_NOy[keys2use])

    # - Get the data and plot up spatially
    dates2use = [datetime.datetime(2018, 8, 1)]
    ref_spec = 'N'
    rm_strat = True
    use_time_in_trop = True
    # Mass unit scaling
    mass_scale = 1E3
    mass_unit = 'Tg'

    # Get all of the speciesConcs for runs as list of datasets
    dsD = {}
    for key in run_dict.keys():
        dsD[key] = AC.GetSpeciesConcDataset(wd=run_dict[key],
                                            dates2use=dates2use)
    # Look runs and extract
#    keys2use = ['ACID', 'BASE']
#    keys2use = [ 'BASE', 'ACID.III', 'ACID.IV',][::-1]
    keys2use = ['BASE', 'ACID.IV', 'JNIT', 'JNITx25'][::-1]
    dsD_NOy = {}
    for key in keys2use:
        print(key, run_dict[key])
        # Use the appropriate list of species
        if ('ACID' in key) or ('NIT' in key):
            specs2use = NOySpecsA
        else:
            specs2use = NOySpecs
        vars2use = [prefix+i for i in specs2use]
        StateMet = AC.get_StateMet_ds(wd=run_dict[key], dates2use=dates2use)
        # Average burden over time
        ds = dsD[key].copy()  # .mean(dim='time', keep_attrs=True)
        ds = AC.get_Gg_trop_burden(ds, vars2use=vars2use, StateMet=StateMet,
                                   use_time_in_trop=use_time_in_trop,
                                   avg_over_time=avg_over_time,
                                   sum_spatially=False,
                                   rm_strat=rm_strat)
        # convert to ref spec equivalent (e.g. N for NO2, C for ACET)
        for spec in specs2use:
            #                ref_spec = get_ref_spec(spec)
            val = ds[prefix+spec]
            ds[prefix+spec] = val/species_mass(spec)*species_mass(ref_spec)
            # Upate varnames
        dsD_NOy[key] = ds[vars2use]

    #
    for key in dsD_NOy.keys():
        ds = dsD_NOy[key]
        # Make a total NOy value
        VarName = 'NOy'
        ds[VarName] = ds[list(ds.data_vars)[0]].copy()
        for var in list(ds.data_vars)[1:]:
            ds[VarName] += ds[var]
        # Make a total aerosol nitrate variable
        VarName = 'NIT_all'
        vars2use = [i for i in ds.data_vars if 'NIT' in i]
        print(vars2use)
        if len(vars2use) > 1:
            ds[VarName] = ds[vars2use[0]].copy()
            for var in vars2use[1:]:
                ds[VarName] + ds[var]
        # Make a dust nitrate variable
        vars2use = [i for i in ds.data_vars if 'NITD' in i]
        VarName = 'NITD'
        print(vars2use)
        if len(vars2use) > 1:
            ds[VarName] = ds[vars2use[0]].copy()
            for var in vars2use[1:]:
                ds[VarName] += ds[var]
        # add % NITS
        VarName = 'aer. NIT (% of NOy)'
        ds[VarName] = ds[list(ds.data_vars)[0]].copy()
        ds[VarName] = ds['NIT_all'] / ds['NOy'] * 100
        # add % NITD
        try:
            ds['NITD']
            VarName = 'aer. NIT(D) (% of NOy)'
            ds[VarName] = ds[list(ds.data_vars)[0]].copy()
            ds[VarName] = ds['NITD'] / ds['NOy'] * 100
        except KeyError:
            pass

        # add back into dictionary
        dsD_NOy[key] = ds

    # Plot sapatial plots for model runs
#    run2use = 'JNIT'
#    run2use = 'JNITx25'
    run2use = 'BASE'
    for run2use in dsD_NOy.keys():
        ds = dsD_NOy[run2use]
        kwargs = {'vmin': 0, 'vmax': 25, 'extend': 'max'}
        plt.close('all')

        # plot up aerosol nitrate %
        var2plot = 'aer. NIT (% of NOy)'
        AC.quick_map_plot(ds.sel(lev=ds.lev[0]), var2plot=var2plot, **kwargs)
        SaveName = 'ARNA_spatial_NIT_all_pcent_August_{}'
        AC.save_plot(SaveName.format(run2use))
        plt.close('all')

        # plot up aerosol nitrate %
        try:
            var2plot = 'aer. NIT(D) (% of NOy)'
            AC.quick_map_plot(ds.sel(lev=ds.lev[0]), var2plot=var2plot,
                              **kwargs)
            SaveName = 'ARNA_spatial_NIT_dust_pcent_August_{}'
            AC.save_plot(SaveName.format(run2use))
            plt.close('all')
        except KeyError:
            pass


def explore_ARNA_period_with_acid_uptake():
    """
    Explore the 3rd set of acid uptake runs over the ARNA time/space
    """
    from arna import add_derived_GEOSChem_specs2ds
    from funcs4obs import gaw_2_loc
#    from matplotlib.ticker import FormatStrFormatter
    RunDir = ar.get_local_folder('RunRoot')
    REF1 = 'Acid-4x5-J00'
    # Just use the runs from the core code
    RunSet = 'ACID'
    res = '4x5'
#    RunSet = 'FP-Nest'
    RunDict = ar.get_dict_of_GEOSChem_model_output(res=res, RunSet=RunSet)
    # Use the spun up or un spun up output?
    for key in RunDict.keys():
        RunDict[key] = RunDict[key]+'OutputDir/'
#        RunDict[key] = RunDict[key]+'spin_up/'

    # Use spin up year - just for testing
#    sdate = datetime.datetime(2018, 1, 1, )
#    edate = datetime.datetime(2019, 1, 1, )
    # Use spun up year before the campaign (2019)
    sdate = datetime.datetime(2019, 1, 1, )
    edate = datetime.datetime(2019, 12, 31, )
    # Just use the January 2020, to compare against sonde data
    sdate = datetime.datetime(2019, 1, 1, )
    edate = datetime.datetime(2020, 1, 29, )
    # Just use the dates of the ARNA campaign
    sdate = datetime.datetime(2019, 2, 1, )
    edate = datetime.datetime(2020, 2, 29, )
    # Just use the dates of the ARNA campaign - but in 2019
#    sdate = datetime.datetime(2019, 1, 1, )
#    edate = datetime.datetime(2019, 3, 1, )
    dates2use = pd.date_range(sdate, edate, freq='1H')

    # Get Key statistics for run
    ExSpecs = [
        'HNO2', 'HNO3', 'NIT', 'NITs',
        'NIT-all', 'NOy-gas', 'NOy', 'SO2', 'SOx', 'SO4-all', 'SO4', 'SO4s',
    ]
    df = AC.get_general_stats4run_dict_as_df(run_dict=RunDict,
                                             REF_wd=RunDict[REF1],
                                             extra_surface_specs=ExSpecs,
                                             extra_burden_specs=ExSpecs,
                                             dates2use=dates2use,
                                             use_REF_wd4Met=True,
                                             debug=debug)

    # Calculate lightning source and add to pd.DataFrame
    var2use = 'EmisNO_Lightning'
    varName = 'Lightning (Tg N/yr)'
    dsD = {}
    for key in RunDict.keys():
        dsD[key] = AC.get_HEMCO_diags_as_ds(wd=RunDict[key],
                                            dates2use=dates2use)
    df2 = pd.DataFrame()
    for key in RunDict.keys():
        ds = dsD[key]
        val = (ds[var2use].mean(dim='time').sum(dim='lev') * ds['AREA'])
        val2 = val.values.sum() * 60 * 60 * 24 * 365  # => /yr
        df2.loc[key, varName] = val2*1E3/1E12
    # add into main DataFrame
    df = df.T
    df = pd.concat([df, df2], axis=1)
    df = df.T

    # Update the units (from NO/kg/yr to N/kg/yr)
#    df = df / AC.species_mass('NO') * AC.species_mass('N')

    # Rename the keys to more readable names
#     rename_dict = {
#     'BASE.BC': 'BASE',
#     'BASE.GFASx2.BC': 'BASE.BBx2',
#     'ACID.BC': 'ACID',
#     'ACID.JNITx25.BC': 'ACID.JNITx25'
#     }
#     for key in list(RunDict.keys()):
#         print(key, (key in rename_dict.keys()))
#         if (key in rename_dict.keys()):
#             NewKey = rename_dict[key]
#             RunDict[NewKey] = RunDict[key]
#             RunDict.pop(key)

    # Get the model for all species by run
    prefix = 'SpeciesConc_'
    dsD = {}
    for key in RunDict.keys():
        ds = AC.get_GEOSChem_files_as_ds(wd=RunDict[key], dates2use=dates2use)
        ds = add_derived_GEOSChem_specs2ds(ds, prefix=prefix)
        dsD[key] = ds
    dsS = {}
    for key in RunDict.keys():
        dsS[key] = AC.get_StateMet_ds(wd=RunDict[key])

    ModelAlt = AC.gchemgrid('c_km_geos5')
    specs2plot = [i for i in ds.data_vars if prefix in i]
    specs2plot = [i.split(prefix)[-1] for i in specs2plot][::-1]

    import seaborn as sns
    sns.set(color_codes=True)
    sns.color_palette('colorblind')
    sns.set_context(context)
    save2png = False
    show_plot = False
    font_scale = 1
    site = 'CVO'
    lat, lon, alt, TZ = gaw_2_loc(site)

    def __get_data4plotting(ds, spec='O3', lat=lat, lon=lon,
                            ModelAlt=ModelAlt):
        """
        Helper function to get data for plotting at a given location
        """
        ds_tmp = ds[prefix+spec].sel(lat=lat, lon=lon, method='nearest')
        ds_tmp = ds_tmp.mean(dim='time')
        ds_tmp['lev'] = ModelAlt
        units, scalby = AC.tra_unit(spec, scale=True)
        # Override default scaling for HONO (HNO2)
        force_units_as_pptv = [
            'NO', 'NO2', 'NOx', 'HNO2'
            'NOy', 'HNO3', 'NIT-all', 'NOy-gas',
            'NOy-HNO3', 'NOy-HNO3-PAN', 'NOy-Limited',
        ]
        if (spec in force_units_as_pptv):
            units = 'pptv'
            scalby = 1E12
        ds_tmp *= scalby
        return ds_tmp, units

    def __beautify_plot(spec='O3', units='ppbv', ylim=(0, 15), ax=None,
                        yscale='linear', xscale='linear'):
        """
        Helper function to beautify the plot
        """
        if isinstance(ax, type(None)):
            ax = plt.gca()
        plt.ylabel('{} ({})'.format('Altitude', 'km'))
        plt.xlabel('{} ({})'.format(spec, units))
        plt.title('Annual mean vertical profile of {} at CVAO'.format(spec))
        plt.legend()
        plt.ylim(ylim)
        ax.set_yscale(yscale)
        ax.set_xscale(xscale)
        # Set ranges for species of interest
        if spec == 'O3':
            plt.xlim(-5, 125)
        elif spec == 'NO' and (xscale != 'log'):
            plt.xlim(-0.1, 100)
        elif spec == 'NO2' and (xscale != 'log'):
            plt.xlim(-0.05, 0.2*1E3)
        elif (spec == 'NOx') and (xscale != 'log'):
            plt.xlim(-0.05, 0.4*1E3)
#         elif spec == 'HNO2' and (xscale != 'log'):
#             plt.xlim(-0.001, 0.003*1E3)
#         elif spec == 'HNO3':
#             plt.xlim(-0.05, 0.5)
        elif spec == 'NOy':
            plt.xlim(-0.05, 1.0)
        elif spec == 'NOy-gas':
            plt.xlim(-0.05, 1.0)
        elif spec == 'NOy-HNO3':
            plt.xlim(-0.05, 1000)
        elif spec == 'NOy-HNO3-PAN':
            plt.xlim(-0.05, 500)
        elif spec == 'NOy-Limited':
            plt.xlim(-0.05, 500)
        elif spec == 'CO':
            plt.xlim(25, 120)
        else:
            pass
        # update the range for log scale plots
        if xscale == 'log':
            plt.xlim(0.1, 300)
#            plt.tick_params(axis='x', which='minor')
#            ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))

    # -- Plot and save to single PDF file
    specs2plot = [
        'O3', 'CO', 'NOx', 'NO2', 'NO', 'HNO2', 'NOy', 'HNO3', 'NIT-all',
        'NOy-gas',
        #    'NOy-HNO3', 'NOy-HNO3-PAN', 'NOy-Limited',
    ]
    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_vertical_above_CVAO_GEOSChem_campaign_global'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    # Species list, then loop dictionary of datasets
    for spec in specs2plot:
        for key in dsD.keys():
            ds = dsD[key].copy()
            ds, units = __get_data4plotting(ds, spec=spec)
            #  plot up...
            ds.plot(y='lev', label=key)
            __beautify_plot(spec=spec, units=units)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        if show_plot:
            plt.show()
        plt.close()
        # Do some memory management...
        gc.collect()

    # - Also add specific ratio plots

    # NO / NOx
    spec = 'NO:NOx'
    for key in dsD.keys():
        ds = dsD[key].copy()
        dsI, units = __get_data4plotting(ds, spec='NO')
        dsII, units = __get_data4plotting(ds, spec='NOx')
        ds = dsI / dsII
        #  plot up...
        ds.plot(y='lev', label=key)
        __beautify_plot(spec=spec, units=units)
    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
    if show_plot:
        plt.show()
    plt.close()

    # HONO / NOx
    spec = 'HONO:NOx'
    for key in dsD.keys():
        ds = dsD[key].copy()
        dsI, units = __get_data4plotting(ds, spec='HNO2')
        dsII, units = __get_data4plotting(ds, spec='NOx')
        ds = dsI / dsII
        #  plot up...
        ds.plot(y='lev', label=key)
        __beautify_plot(spec=spec, units=units)
    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
    if show_plot:
        plt.show()
    plt.close()

    # HONO / NOy
    spec = 'HONO:NOy'
    for key in dsD.keys():
        ds = dsD[key].copy()
        dsI, units = __get_data4plotting(ds, spec='HNO2')
        dsII, units = __get_data4plotting(ds, spec='NOy')
        ds = dsI / dsII
        #  plot up...
        ds.plot(y='lev', label=key)
        __beautify_plot(spec=spec, units=units)
    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
    if show_plot:
        plt.show()
    plt.close()

    # Add a log NO, NO2 and NOx plot
    specs2plot_AsLog = ['NO', 'NO2', 'NOx', ]  # 'HNO2']
    # Species list, then loop dictionary of datasets
    for spec in specs2plot_AsLog:
        for key in dsD.keys():
            ds = dsD[key].copy()
            ds, units = __get_data4plotting(ds, spec=spec)
            #  plot up...
            ds.plot(y='lev', label=key)
            __beautify_plot(spec=spec, units=units, xscale='log')
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        if show_plot:
            plt.show()
        plt.close()
        # Do some memory management...
        gc.collect()

    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')

    # - Plot up the surface difference between the runs
    matplotlib.rc_file_defaults()
    sns.reset_orig()
    REF = 'Acid-4x5-J00'
#    DIFF = 'Acid-4x5-J50'
    DIFF = 'Acid-4x5-Isotherm.v2'
    prefix = 'SpeciesConc_'
    vars2plot = ['O3', 'CO', 'NOx', 'NO2', 'NO', 'HNO2',
                 'NOy', 'HNO3', 'NIT-all', 'NOy-gas',
                 ]
    vars2plot = [prefix+i for i in vars2plot]
    pcent = True
    for DIFF in ['Acid-4x5-J50', 'Acid-4x5-Isotherm.v2']:
        savetitle = 'ARNA_surface_plots_{}_vs_{}'.format(REF, DIFF)
        AC.plt_spatial_diff_between_runs_at_lvl(dsD, REF=REF, DIFF=DIFF,
                                                pcent=pcent,
                                                savetitle=savetitle,
                                                vars2plot=vars2plot,
                                                prefix=prefix)

    # - Plot up the zonal difference in concentration between the runs
    REF = 'Acid-4x5-J00'
#    DIFF = 'Acid-4x5-J50'
    DIFF = 'Acid-4x5-Isotherm.v2'
    prefix = 'SpeciesConc_'
    vars2plot = ['O3', 'CO', 'NOx', 'NO2', 'NO', 'HNO2',
                 'NOy', 'HNO3', 'NIT-all', 'NOy-gas',
                 ]
    vars2plot = [prefix+i for i in vars2plot]
    pcent = True
#    for key in dsD.keys():
#        dsD[key] =  dsD[key].mean(dim='time')
    for DIFF in ['Acid-4x5-J50', 'Acid-4x5-Isotherm.v2']:
        savetitle = 'ARNA_zonal_plots_{}_vs_{}'.format(REF, DIFF)
        AC.plt_zonal_diff_between_runs(dsD, REF=REF, DIFF=DIFF,
                                       StateMet=dsS[DIFF],
                                       pcent=pcent,
                                       savetitle=savetitle,
                                       vars2plot=vars2plot,
                                       prefix=prefix)


def plt_key_NOx_budget_terms(RunDict=None):
    """
    Make spatial plots of key tagged routes
    """
    # Local switches
    PltAsLog = True
    verbose = True
    # Which runs to use?
    dates2use = [datetime.datetime(2019, i+1, 1) for i in range(12)]
#    dates2use = [datetime.datetime(2018, i+1, 1) for i in range()]
    trop_limit = False
    RunSet = 'PostHONO'
#    RunSet = 'IGACset'
#    GC_version = 'v13.4'
    GC_version = 'v12.9'
    res = '4x5'
    if instance(RunDict, type(None)):
        RunDict = ar.get_dict_of_GEOSChem_model_output(RunSet=RunSet,
                                                       GC_version=GC_version,
                                                       res=res)
    # Which runs to use?
    if RunSet == 'PostHONO':
        runs2use = ['Base', 'min4pptHONO', 'NIThv']
#        DIFFrun = 'min4pptHONO'
        DIFFrun = 'NIThv'
        REF_list = ['Base', 'min4pptHONO',]
        KeyList = list(sorted(RunDict.keys()))
        for key in KeyList:
            if key in runs2use:
                pass
            else:
                del RunDict[key]
    elif RunSet == 'IGACset':
        runs2use = NOxD.keys()
        DIFFrun = 'Iso.UnlimitedAll'
        REF_list = ['Base']

    else:
    #    for run in ['Acid-4x5-J50']
    #    run2use = 'Acid-4x5-J50'
    #    run2use = 'Acid-4x5-Isotherm.v2'
        runs2use = ['Acid-4x5-Isotherm.v2', 'Acid-4x5-J50', 'Acid-4x5-J00', ]
        REF_list = ['Acid-4x5-J50', 'Acid-4x5-J00']
        DIFFrun = 'Acid-4x5-Isotherm.v2'
    # Get NOx budget dictionary
    NOxD = ar.get_NOx_budget_ds_dict_for_runs(RunDict=RunDict,
                                              trop_limit=trop_limit,
                                              dates2use=dates2use)
    # Get the tags for the Prod/Loss tagging as a dictionary
    tags = ar.get_tags_for_NOx_HONO()
    ProdVars = ['ProdHNO2fromHvNIT-all',
                'ProdHNO2fromOHandNO',
                'ProdHNO2fromHET']

    # - Plot up the annual mean surface values
    for run2use in runs2use:
        ds = NOxD[run2use].copy()
        ds = ds.mean(dim='time').sel(lev=ds.lev[0])
        norm = LogNorm(vmin=0.5, vmax=2000)
#        norm = LogNorm(vmin=10, vmax=1000)
#        norm = LogNorm(vmin=10, vmax=1000)
#        norm = LogNorm(vmin=50, vmax=1000)
        norm = LogNorm(vmin=100, vmax=1000)
        savetitle = 'ARNA_spatial_NOx_budget_{}'.format(run2use)
        if PltAsLog:
            kwargs = {'norm': norm}
            savetitle += '_Log'
        else:
            kwargs = {}
        # Setup a PDF
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # Plot key route totals
        for var2plot in ['Prod_OH', 'Prod_HNO2', 'Jscale']:
            #    for var2plot in [ 'Prod_HNO2']:
            AC.quick_map_plot(ds, var2plot=var2plot, verbose=verbose, **kwargs)
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
    #        del ds2plot

        # Plot routes
        ds2plot = ds[ProdVars].copy()
        for var2plot in ProdVars:
            AC.quick_map_plot(ds2plot, var2plot=var2plot,
                              verbose=verbose, **kwargs)
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
    #        del ds2plot

        # Now plot routes as % of total HNO2 production
        ds2plot = ds[ProdVars].copy() / ds['Prod_HNO2']
        for var2plot in ProdVars:
            AC.quick_map_plot(ds2plot, var2plot=var2plot, verbose=verbose,
                              **kwargs)
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
    #        del ds2plot

        # Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')

        # - Zonally plot up the NOx budget values
        savetitle = 'ARNA_zonal_NOx_budget_{}'.format(run2use)
        if PltAsLog:
            kwargs = {'norm': norm}
            savetitle += '_Log'
        else:
            kwargs = {}
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # Plot routes
        ds = NOxD[run2use].copy().mean(dim='time')
        vars2plot = ['Prod_OH', 'Prod_HNO2', 'Jscale'] + ProdVars
        ds2plot = ds[vars2plot]
        for var2plot in vars2plot:
            fig, ax = plt.subplots()
            im = AC.ds2zonal_plot(ds2plot, var2plot=var2plot, StateMet=ds,
                                  fig=fig, ax=ax, **kwargs)
            TitleStr = "{}"
            plt.title(TitleStr.format(var2plot))

            # Add a colourbar
    #        kwargs = {'extend':'both'}
            fig.colorbar(im, orientation="horizontal", pad=0.2, extend='both',
                         **kwargs)
    #                     format=format, label=units)
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
    #        del ds2plot

        # Now plot routes as % of total HNO2 production
        ds = NOxD[run2use].copy().mean(dim='time')
        ds2plot = ds[ProdVars] / ds['Prod_HNO2']
    #    kwargs = {'vmin': 0, 'vmax': 1, }
        for var2plot in ProdVars:
            fig, ax = plt.subplots()
            im = AC.ds2zonal_plot(ds2plot, var2plot=var2plot, StateMet=ds,
                                  fig=fig, ax=ax, **kwargs)
            TitleStr = "'{}' / total HONO production"
            plt.title(TitleStr.format(var2plot))

            # Add a colourbar
            fig.colorbar(im, orientation="horizontal", pad=0.2, extend='both',
                         **kwargs)
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
    #        del ds2plot

        # Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')

    # --- Plot up source of HONO in isotherm runs
    # (relative to J50 and J00 if before "PostHONO" set)
    run2use = DIFFrun
    for REFrun in REF_list:
        savetitle = 'ARNA_NOx_budget_REF_{}_vs_{}'.format(REFrun, run2use)
        # Plot routes
        REF = NOxD[REFrun].copy().mean(dim='time')
        ds = NOxD[run2use].copy().mean(dim='time')
        vars2plot = ['Prod_OH', 'Prod_HNO2', 'Jscale'] + ProdVars
        ds2plot = ds[vars2plot] / REF[vars2plot]
        # Log plots?
        if PltAsLog:
            kwargs = {'norm': norm}
            savetitle += '_Log'
        else:
            kwargs = {}
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        for var2plot in vars2plot:
            fig, ax = plt.subplots()

            # - Spatial plot
            AC.quick_map_plot(ds2plot.sel(lev=ds.lev[0]), var2plot=var2plot,
                              verbose=verbose, **kwargs)
            TitleStr = "{} ({} vs. {})"
            plt.title(TitleStr.format(var2plot, REFrun, run2use))

            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()

            # - Zonal plot
            im = AC.ds2zonal_plot(ds2plot, var2plot=var2plot, StateMet=ds,
                                  fig=fig, ax=ax, **kwargs)
            TitleStr = "{} ({} vs. {})"
            plt.title(TitleStr.format(var2plot, REFrun, run2use))
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


# def calculate_NIT_production_totals(RunDict=None, trop_limit=True,):
#     """
#     Calculate the global nitrate production routes.
#     """
#
#     pass



def analyse_NOx_budget(RunDict=None, trop_limit=True,):
    """
    Analyse NOx budget in tagged model runs
    """
    # Use the testing runs for func unless others provided
    if isinstance(RunDict, type(None)):
        DataRoot = ar.get_local_folder('RunRoot')
        RunPrefix = 'geosfp_4x5_aciduptake.v12.9.0.ARNA.Isotherm.Diags.'
        FullRootStr = '{}/{}{}/OutputDir/'
        RunDict = {
        'Iso.Unlimited': FullRootStr.format(DataRoot, RunPrefix,
                                           'v9.Iso.UnlimitedAll'),
        'Base': FullRootStr.format(DataRoot, RunPrefix, 'v9.Base'),
        }

    #ÊLocal variables and switches
    AVG = AC.constants('AVG') # Avogadros constant
    if isinstance(dates2use, type(None)):
#        dates2use = [datetime.datetime(2019, i+1, 1) for i in range(12)]
        sdate = datetime.datetime(2019, 1, 1) # Beginning for analysis year
#        edate = datetime.datetime(2019, 3, 1) # 3 months into analysis year
        edate = datetime.datetime(2019, 12, 31) # 3 months into analysis year
        dates2use = pd.date_range(sdate, edate, freq='1D')

    trop_limit = True
#    MaskName = None
    MaskName = 'inflow_CVAO_area'
#    ApplyMask = True
    ApplyMask = False
    LoadSavedNetCDFs = True
    AvgOverTime = True

    # Get NOx budget dictionary of objects (saved offline on extracted online)
    if LoadSavedNetCDFs:
        REF_str = 'ARNA_Testing'
        for key in NOxD.keys():
            print( key )
            try:
                pass
            except:
                print('WARNING')

    else:
        NOxD = ar.get_NOx_budget_ds_dict_for_runs(trop_limit=trop_limit,
                                                  dates2use=dates2use,
                                                  RunDict=RunDict,
                                                  IncAllProdLossDiags=True,
                                                  ApplyMask=ApplyMask,
                                                  MaskName=MaskName,
#                                              CoreRunsOnly=True,
                                                  )


    # - - Now do analysis of the data - -
    # Remove the prod/loss variables with non standard output
    var2skip = [
       'ProdCOfromNMVOC', 'ProdCOfromCH4',
       'ProdSO4fromOxidationOnDust',
       'ProdNITfromHNO3uptakeOnDust',
       'ProdSO4fromUptakeOfH2SO4g',
       'ProdSO4fromO3s', 'ProdSO4fromSRHObr', 'ProdSO4fromSRO3',
       'ProdSO4fromHOBrInCloud', 'ProdSO4fromO3inSeaSalt',
       'ProdSO4fromO3inCloud', 'ProdSO4fromO2inCloudMetal',
       'ProdSO4fromH2O2inCloud', 'ProdOCPIfromOCPO', 'ProdBCPIfromBCPO',
        ]

    # - Global (values in bottom X km) totals (or averages for s^-1 variables)
    df = pd.DataFrame()
    for key in NOxD.keys():
        # Get dataset and setup Series for data
        S = pd.Series()
        ds = NOxD[key]
        #
#        NIU, NIU, Alt = AC.get_latlonalt4res(res='4x5', full_vert_grid=True)
#        AltBool = Alt < 5

        # Get the total (or average?) production via route
        # kg N s-1
        vars2use = [i for i in ds.data_vars if 'Prod' in i]
        vars2use += [i for i in ds.data_vars if 'Phot' in i]
        vars2use += [i for i in ds.data_vars if 'Loss' in i]
        # Exclude the non standard units
        vars2use = [ i for i in vars2use if i not in var2skip]

        for var in vars2use:
            print(var)
            units = ds[var].units
            data = ds[var].copy()
#            data = data.isel(lev=AltBool) # masking troposphere, so not needed
            if debug:
                print(units)
            if units == 'molec cm-3 s-1':
                # convert to kg (N) / s-1 (=>mol => *m2 (*1e6) => kg (*1E3) )
                data = data / AVG * ds['Met_AIRVOL'] * 1E6 * 14 * 1E3
#                if sum_data:
                data = data.sum(dim=['lat', 'lon', 'lev'])
            elif (units == 'kg N s-1'):
                #                if sum_data:
                data = data.sum(dim=['lat', 'lon', 'lev'])
            elif (units == 'kg s-1'):
                # This is HNO3 loss oan sea-salt (are the units correct?)
                #                if sum_data:
                HNO3seasaltVar = 'Loss of HNO3 on sea salt aerosols'
                if HNO3seasaltVar == ds[var].long_name:
                    data = data / species_mass('HNO3') * 14.0
                data = data.sum(dim=['lat', 'lon', 'lev'])
            else:
                PrtStr = 'WARNING: No conversion setup for Units ({})'
                print(PrtStr.format(units))
            if AvgOverTime:
                data = data.mean(dim='time')
            S[var] = data.values

        # Also include calculations for NIT production
        # (and convert to average flux... )
        var2use = 'ProdNITfromHNO3uptakeOnDust'
        data = ds[var2use].copy()
        data = data.sum(dim=['lat', 'lon', 'lev'])
        if AvgOverTime:
            data = data.mean(dim='time')

        # Already in dataframe
#        var2use = 'LossHNO3onSeaSalt'
#        ds[var2use]

        # Store in the dataframe
        df[key] = S


#    for var in var2del:
#        del df.T[var]

    # Get deposition sinks too
    vars2use = [
    'DryAndWetDep_NIT-all',
    'DryAndWetDep_HNO3',
    'DryAndWetDep_NO2',
    'DryAndWetDep_NIT',
    'DryAndWetDep_NITs',
    ]
    prefix = 'DryAndWetDep_'

    for key in NOxD.keys():
        ds = NOxD[key]
        for __var in vars2use:
            # Species RMM
            Species = __var.split(prefix)[-1]
            SpeciesRMM = AC.species_mass( Species )
            # Convert kg (NO) /m2/s => kg (NO) /s
            data = ds[__var] * ds['AREA']
            # & convert kg X => Kg  N
            data = data / SpeciesRMM * 14.
            # Convert to Tg from kg?
#            data = data * 1E3 / 1E12
            # convert from /s => /year
#            data = data *60.*60.*24*.365
            # all 12 months present?
#            DsMonths = ds['time.month'].values
#            if DsMonths == np.arange(1, 13):
#            data = data.values.sum() / 12
#            else:
            data = data.mean(dim='time').sum().values
            df.loc[ __var, key ] = data
#             if debug:
#                 print('WARNING: Assuming ')

    # Get Emission sources
    var2use = 'EmisNO_Total'
    EmissD = {}
    for key in NOxD.keys():
        ds = AC.get_HEMCO_diags_as_ds(wd=RunDict[key],
                                      dates2use=dates2use)
        EmissD[key] = ds
        # Convert kg (NO) /m2/s => kg (NO) /s
        data = ds[var2use] * ds['AREA']
        # & convert kg NO => Kg N
        data = data / 30.0 * 14.
        # convert kg (N) => Tg (N)
#        data = data * 1E3 / 1E12
        # convert from /s => /year
#        data = data *60.*60.*24.*365.
        # Archive an average over time
        data = data.mean(dim='time').sum().values
        df.loc[var2use, key ] = data

    # Print DataFrame and save full numbers to csv
    print(df)
    SaveName4CSV = 'ARNA_NOx_Budget_'
    if ApplyMask:
        SaveName4CSV += 'Masked_{}'.format(MaskName)
    SaveName2Use = '{}.{}'.format(SaveName4CSV, 'csv')
    df.to_csv(SaveName2Use)

    # Consider the values relative to the J00 run
    df2 = df.copy()
    REF = 'Base'
    for col in [i for i in df2.columns if i != REF]:
        df2.loc[:, col] = df2.loc[:, col] / df2[REF]
    print(df2)
    SaveName2Use = '{}_REF_{}.{}'.format(SaveName4CSV, REF, 'csv')
    df2.to_csv(SaveName2Use)

    # Consider the values relative to the J50 run
    df3 = df.copy()
    REF = 'Acid-4x5-J50'
    for col in [i for i in df3.columns if i != REF]:
        df3.loc[:, col] = df3.loc[:, col] / df3[REF]
    print(df3)
    SaveName2Use = '{}_REF_{}.{}'.format(SaveName4CSV, REF, 'csv')
    df3.to_csv(SaveName2Use)

    # Consider ...

    # - Just consider the tagged PNOx reaction routes.
    REF = 'Base'
    DIFF = 'Iso.Unlimited'
    prefix = 'ProdfromRXN_'
    # select core routes
    vars2use = [i for i in df.T.columns if prefix in i]
    df4 = df.T[vars2use].T
    # Add ratio as a column
    RatioVar = '{}/{}'.format(DIFF, REF)
    df4[RatioVar] = df4[DIFF]/df4[REF]
    df4 = df4.sort_values(REF, ascending=False )
    SaveName2Use = '{}_REF_{}_Descending.{}'.format(SaveName4CSV, REF, 'csv')
    df4.to_csv(SaveName2Use)
    # As above, but a version based on which ones have changed
    df4.sort_values(RatioVar, ascending=False)


    # - Print out core values for NOx Budget diagram.
    df_BACKUP = df.copy()
    vars4plot = [
    'EmisNO_Total',
    'DryAndWetDep_NIT-all',
    'DryAndWetDep_HNO3',
    ]

    # And chemistry routes?
    vars4plot += [
     'ProdHNO2fromHvNIT-all', # TN001-TN007
     'ProdNITsfromHetNO3', # 'TN024',
     'ProdNITfromHetNO3', # 'TN025',
     'ProdHNO3fromHetNO3', # 'TN023',
     'ProdNO3fromHNO3nOH' ,  # 'TN019',
     'ProdHNO3fromNO2nOH',#  # 'TN018',
     'ProdNOnHO2ChannelA', #  # 'TN016',
     'ProdHNO2fromHET',  # 'TN015',
     'ProdHNO2fromOHandNO',  # 'TN014',
     'ProdNO2fromHONO',  # 'TN013',
     'PhotNO2', #: 'TN020', # 'TN020',
     'PhotHNO3', #: 'TN021', # 'TN021',
     'PhotHNO2', #: 'TN022', # 'TN022',
    ]

    # Any extra reactions useful
    vars4plot += [
    'LossHNO3onSeaSalt',
    'ProdNITfromHNO3uptakeOnDust',
    ]

    #
    df5 = df.T
#    df5.T[vars4plot].T
    df5[vars4plot].T.to_csv('ARNA_NOx_diagram_numbers.csv')





def plt_NOx_routes_by_lat(NOxD=None, ):
    """
    Plot up NOx production routes by lat
    """
    # Function expects to be provided with the NOx dictionary of datasets
    if isinstance(NOxD, type(None)):
        NOxD = ar.get_NOx_budget_ds_dict_for_runs(trop_limit=trop_limit,
                                                  dates2use=dates2use,
                                                  RunDict=RunDict,
                                                  IncAllProdLossDiags=True,
                                                  ApplyMask=ApplyMask,
                                                  MaskName=MaskName,
#                                              CoreRunsOnly=True,
                                                  )

    # The unit conversions to Kg (N)/s have been done already
    #  Convert to zonal values too
    for key in NOxD.keys():
        # Select the variables from the dataset
        ds = NOxD[key]
        vars2use = [i for i in ds.data_vars if 'Prod' in i]
        vars2use += [i for i in ds.data_vars if 'Loss' in i]

        # TODO: Add molecule weighting...
        # temp average lon, time
        ds = ds.mean(dim=('time',), keep_attrs=True)

        # weight by molecules
        for __var in vars2use:
            attrs = ds[__var].attrs.copy()
            ds[__var] = (ds[__var] * ds['Met_MOLES']).sum(dim='lon',)
            ds[__var] = ds[__var] / ds['Met_MOLES'].sum(dim='lon',)
            ds[__var].attrs = attrs

        # sum over lev and update in dictionary of datasets
        ds = ds.sum(dim=('lev'), keep_attrs=True)
        NOxD[key] = ds

    # - Plot up PNOx variables by lat
    ColorList = AC.get_CB_color_cycle()
    fig, ax = plt.subplots()
#    keys2plot = list( NOxD.keys())
    keys2plot = ['Iso.Unlimited']
    for key in keys2plot:
        ds = NOxD[key]

        # plot total PNOx
        PNOxVar = 'Prod_NOx'
        var2plot = PNOxVar
        units = 'molec cm-3 s-1'
        Y = ds[var2plot].values
        X = ds['lat'].values
        plt.plot(X, Y, label=var2plot, color=ColorList[0])
        plt.ylabel('{} ({})'.format(PNOxVar, units))
        plt.xlabel('Latitude ($^{\circ}$N)')

#        plt.legend()

        # Underlay % NOx production from nitrate photolysis
        # Twin y axis and show route contribution to PNOx
        ax2 = ax.twinx()
        var2plot =  'ProdHNO2fromHvNIT-all'
        units = '%'
        Y = ds[var2plot].values  / ds[PNOxVar].values *100
        X = ds['lat'].values
        plt.plot(X, Y, label=var2plot, color=ColorList[1])
        plt.ylabel( '{} {}'.format(var2plot, units))
#
#        plt.legend()
    # Beatify plot and save
    fig.legend()
    SaveTitle = 'ARNA_PNOx_vs_lat'
    AC.save_plot(title=SaveTitle)
    plt.close('all')


    # - plot HNO2 key routes and total production
    keys2plot = ['Iso.Unlimited']
    for key in keys2plot:
        ds = NOxD[key]

        #
        PHNO2var = 'Prod_HNO2'
        var2plot = PHNO2var
        units = ds[var2plot].units
        Y = ds[var2plot].values
        X = ds['lat'].values
        plt.plot(X, Y, label=var2plot, color=ColorList[0])
        plt.ylabel('{} ({})'.format(var2plot, units))
        plt.xlabel('Latitude ($^{\circ}$N)')

        # then plot a stack plot of the other routes?
        vars2plot = ['ProdHNO2fromHvNIT-all', 'ProdHNO2fromHET',
                     'ProdHNO2fromOHandNO'
                     ]
        ax2 = ax.twinx()
        # Quick plot of lines individually
        for var2plot in vars2plot:
            print(var2plot)
            units = '%'
            Y = ds[var2plot].values  / ds[PHNO2var].values *100
            print(Y)
            X = ds['lat'].values
            plt.plot(X, Y, label=var2plot, color=ColorList[1])
            plt.ylabel( '{} {}'.format(var2plot, units))

        #



#        var2plot = 'ProdHNO2fromHET'


#        var2plot = 'ProdHNO2fromOHandNO'

    # Beatify plot and save
    fig.legend()
    SaveTitle = 'ARNA_PNO2_vs_lat'
    AC.save_plot(title=SaveTitle)
    plt.close('all')

    # - Stacked area plot using plot
    ColorList = AC.get_CB_color_cycle()
    fig, ax = plt.subplots()

    # Explore the data
    ds_tmp = ds[vars2plot+[PHNO2var]].squeeze()
    df = ds_tmp.to_dataframe()
    RouteSumVar = 'Sum of routes'
    df[RouteSumVar] = df[vars2plot].sum(axis=1)
    df['% of total'] = df[RouteSumVar] / df[PHNO2var] *100

    #
    col2use = [i for i in df.columns if i != PHNO2var]
    for col in col2use:
        df[col] = df[col] /df[PHNO2var] *100

    #
    df[vars2plot].plot.area()
    fig.legend()
    plt.ylabel('{} ({})'.format('HONO via route', units))

    # Add a x axis over the top to show the total
    ax2 = ax.twinx()
    #
    PHNO2var = 'Prod_HNO2'
    var2plot = PHNO2var
    units = ds[var2plot].units
    Y = ds[var2plot].values
    X = ds['lat'].values
    ax2.plot(X, Y, label=var2plot, color=ColorList[0])
    plt.ylabel('{} ({})'.format(var2plot, units))
    plt.xlabel('Latitude ($^{\circ}$N)')


    SaveTitle = 'ARNA_PNO2_vs_lat_stacked_pcent'
    AC.save_plot(title=SaveTitle)
    plt.close('all')



def plt_lightning_by_month(context='paper'):
    """
    Plot lightning seasonally and explore high values
    """
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context(context)
    folder = ar.get_local_folder('RunRoot')
    folder += 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.'
    folder1 = folder + 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags'
    folder1 = folder1 + '.ToGetNetCDFOutput/OutputAndSpinUpSymLink/'
    folder2 = folder + 'DustUptake.JNITx50.BCs.Light30/OutputDir/'
    d = {'Acid.J50': folder1, 'Acid.J50.LightScale': folder2}
    var2use = 'EmisNO_Lightning'
    varName = 'Lightning (Tg N/yr)'
    MetFields = 'GEOSFP'
#    MetFields = 'MERRA2'
    dfs = {}
    for key in d.keys():
        ds = AC.get_HEMCO_diags_as_ds(wd=d[key])
        val = (ds[var2use].sum(dim='lev') * ds['AREA'])
        val = val.sum('lat').sum('lon') * 60 * 60 * 24 * 365  # => /month
        val = val*1E3/1E12
        # Convert NO => N
        val = val / AC.species_mass('NO') * AC.species_mass('N')
        dfs[key] = val.to_pandas()
    # plot up
    units = 'Lightning (Tg N/year (equiv.))'
    linestyles = AC.get_ls(len(dfs))
    for nKey, key in enumerate(d.keys()):
        dfs[key].plot(label=key, ls=linestyles[nKey])
    plt.title('Global Lightning NOx source in {}'.format(units))
    plt.ylabel(units)
    AC.save_plot(title='ARNA_Global_lightning_source_{}'.format(MetFields))
    plt.close()
    # Plot up with an assumed *1/3 reduction
    df = pd.DataFrame()
    df['v12.9 (GEOS-FP)'] = dfs['Acid.J50']
    df['v12.9 (GEOS-FP) - reduced by 1/3'] = val2*(1-0.33)
    df['v12.9 (GEOS-FP) - reduced by 1/4'] = val2*(1-0.25)
    df['v12.9 (GEOS-FP) - reduced by 1/2'] = val2*(1-0.5)
    df.index.name = None
    # Add means to plot
    df['dates'] = pd.to_datetime(df.index.values)
    means = df.groupby(df.dates.dt.year).mean().round(1)
    del df['dates']
    print(means.T)
    df.plot()
    titleStr = 'ARNA_Global_lightning_source_{}_SCALED'
    AC.save_plot(title=titleStr.format(MetFields))
    plt.close()


def explore_JVALS_in_JNITs_runs():
    """
    Plotting of J-rates for nitrate and sulfate species
    """
    folder = ar.get_local_folder('RunRoot')
    folder += 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.'
    folder += 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags/'
    folder += 'OutputDir/'
    # Species to use?
    TRAs = ['HNO3', 'NIT', 'NITs', 'NITD1', 'NITD2', 'NITD3', 'NITD4']
    Jprefix = 'Jval_'
    vars2use = [Jprefix+i for i in TRAs]
#    dates2use = None
    # Get photolysis data
    dsP = AC.GetJValuesDataset(wd=folder, dates2use=dates2use)
    dsP = dsP[vars2use]
    # Get surface photolysis data
    TRAs += ['SO2', 'SO4', 'SO4s', 'SO4D1', 'SO4D2', 'SO4D3', 'SO4D4']
    dsS = AC.GetSpeciesConcDataset(wd=folder, dates2use=dates2use)
    # Work out which idx to use
    NIU, NIU, alt = AC.get_latlonalt4res(res='4x5')
    lvl_dict = {1: 0.071, 13: 1.999, 21: 4.886, }

    # Setup a PDF
    savetitle = 'ARNA_spatial_photolysis_analysis'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # Loop and plot photolysis rates
    for var in vars2use:
        for lvl_idx in [1, 13, 21]:
            ds2plot = dsP[[var]].mean(dim='time').isel(lev=lvl_idx)

            # Plot up the photolysis rates
            title = "Average '{}' @ {:.1f}km".format(var, lvl_dict[lvl_idx])
            AC.quick_map_plot(ds2plot, var2plot=var, title=title,
                              save_plot=False)
#            plt.title(title)
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()

            # Plot up the photolysis rates (as a fraction of JHNO3)
            REF = dsP[[Jprefix+'HNO3']].mean(dim='time').isel(lev=lvl_idx)
            NewVarName = 'JScale'
            ds2plot[NewVarName] = ds2plot[var].copy()
            ds2plot[NewVarName] = ds2plot[NewVarName] / REF[Jprefix+'HNO3']
            TitleStr = "Ratio of '{}':JHNO3 @ {:.1f}km"
            title = TitleStr.format(var, lvl_dict[lvl_idx])
            AC.quick_map_plot(ds2plot, var2plot=NewVarName, title=title,
                              save_plot=False)

            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
            del ds2plot

    # Also plot up the average surface concentrations of related species
    SCprefix = 'SpeciesConc_'
    vars2use = [SCprefix+i for i in TRAs]
    for var in vars2use:
        #        for lvl_idx in [1, 13, 21]:
        for lvl_idx in [13]:
            ds2plot = dsS[[var]].mean(dim='time').isel(lev=lvl_idx)

            # Plot up the photolysis rates
            AC.quick_map_plot(ds2plot, var2plot=var, title=title,
                              save_plot=False)
            title = "Average '{}' @ {:.1f}km".format(var, lvl_dict[lvl_idx])
            plt.title(title)
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()

    # Save entire pdf
#    if close_pdf:
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')

    # Make a quick summary of the Jscale by NIT species
    REF = 'Jval_HNO3'
    df =  pd.DataFrame()
    for __var in dsP.data_vars:
        S = pd.Series( (dsP[__var] / dsP[REF]).values.flatten() ).describe()
        df[__var] = S
    df.to_csv( '{}.csv'.format(savetitle) )

    # Make a quick summary of the dust update
    df2 = pd.DataFrame()
    dfs = {}
    for dsMonth in dsS['time.month'].values:
        ds = dsS.isel(time=dsS['time.month'] == dsMonth).mean(dim='time')
        for __var in TRAs:
            dsVarName = '{}{}'.format(SCprefix, __var)
            S = pd.Series(ds[ dsVarName ].values.flatten()).describe()
            df2[__var] = S
            dfs[dsMonth] = df2


def explore_Jscale_offline():
    """
    Explore resulting Jscale values from Andersen2022 parameterisation
    """

    # Get data to explore parameter space with
    sdate = datetime.datetime(2018, 6, 1) # Beginning for spin up (6months)
    edate = datetime.datetime(2019, 6, 1) # 12 months into model run
    dates2use = pd.date_range(sdate, edate, freq='1D')
    folder = '/users/ts551/scratch/GC/rundirs/P_ARNA/'
    folder += 'gc_4x5_47L_merra2_fullchem_aciduptake.v13.4.1.ARNAv10/'
    folder += 'OutputDir/'
    ds = AC.GetSpeciesConcDataset(wd=folder, dates2use=dates2use )

    # Setup ancillary data for conversions
    MolecVar = 'Met_MOLECS'
    StateMet = AC.get_StateMet_ds(wd=folder, dates2use=dates2use)
    StateMet = AC.add_molec_den2ds(StateMet, MolecVar=MolecVar)
    AVG = AC.constants( 'AVG' )

    # Add a total nitrate variable
    NO3Var = 'NIT-all'
    ds = AC.AddChemicalFamily2Dataset(ds, fam=NO3Var)
    SALVar = 'SAL-all'
    ds = AC.AddChemicalFamily2Dataset(ds, fam=SALVar)

    # Update names
    SCprefix = 'SpeciesConc_'
    VarNames = ds.data_vars
    NewNames = [i.split(SCprefix)[-1] for i in VarNames]
    ds = ds.rename( dict(zip(VarNames, NewNames)) )

    # Convert values to moles /m3
    # v/v => molecs/cm3 => mol/cm3 =>  mol/m3 => nmol/m3
    for __var in [NO3Var, SALVar]:
        data = ds[__var] * StateMet[MolecVar] / AVG
        data = data * 1E6 * 1E9
        ds[__var].values = data


    # Calculate the resultant Jscale for these values
    # Andersen22a
    Andersen22aVar = 'Andersen22a'
    no3 = ds[NO3Var]
    ds[Andersen22aVar] = ((8.6E2 * np.log(1. + (0.94*no3) ) ) / (no3) ) - 38.1

    # Andersen22b
    Andersen22bVar = 'Andersen22b'
    FACvar = 'NO3/Cl'
    ds[FACvar] =  ds[NO3Var] / ds[SALVar]
    ds[Andersen22bVar] = ( 31.65 * np.log(1. + 6.34 * ds[FACvar])) / ds[FACvar]

    # Also include the inverse chloride to nitrate ratio.
    FACvar_r = 'Cl/NO3'
    ds[FACvar_r] =   ds[SALVar] / ds[NO3Var]


    # print out statistics


#    NITspecies = ['NIT', 'NITs', 'NITD1', 'NITD2', 'NITD3', 'NITD4']
#    ds =


    # Loop and plot values of interest
    vars2plot = NO3Var, SALVar, FACvar, FACvar_r, Andersen22aVar, Andersen22bVar

    ds2plot = ds.copy().mean(dim='time')
    ds2plot = ds2plot.isel(lev=(ds.lev==ds.lev[0])).squeeze()


    savetitle = 'ARNA_chloride_nitrate_ratios'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # plot surface
    for __var in vars2plot:

        kwargs = {}
        AC.quick_map_plot(ds2plot, var2plot=__var, verbose=verbose,
                          **kwargs)
        TitleStr = "Average surface '{}' ".format(__var)
        plt.title(TitleStr.format(__var))

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

    # Also plot up a limited version of the Cl/NO3 ratio
    __var = FACvar
    kwargs = {'vmin':0, 'vmax':25, 'extend': 'both'}
    AC.quick_map_plot(ds2plot, var2plot=__var, verbose=verbose,
                      **kwargs)
    TitleStr = "Average surface '{}' limited (0-25)".format(__var)
    plt.title(TitleStr.format(__var))

    # Also plot up a limited version of the NO3/Cl ratio
    __var = FACvar_r
    kwargs = {'vmin':0, 'vmax':37, 'extend': 'both'}
    AC.quick_map_plot(ds2plot, var2plot=__var, verbose=verbose,
                      **kwargs)
    TitleStr = "Average surface '{}' limited (0-37)".format(__var)
    plt.title(TitleStr.format(__var))

    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
    plt.close()

    # zonal plot of these values too.
    ds2plot = ds.copy().mean(dim='time')
    StateMet2use = StateMet.copy().mean(dim='time')

    for __var in vars2plot:
        var2plot = __var

        fig, ax = plt.subplots()
        kwargs = {}
        im = AC.ds2zonal_plot(ds2plot, var2plot=var2plot,
                              StateMet=StateMet2use,
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


def plot_weighted_nitrate_Jscale()
    """
    Plot the nitrate weighted Jscale
    """
    # - TEMP - setup a local dictionary for data to use
    # Update below as available.
    RunDict = {}
    RunDict['Kas18'] = '/users/ts551/scratch/GC/rundirs/P_ARNA/gc_4x5_47L_merra2_fullchem_aciduptake.v13.4.1.ARNAv10.ARNA.Test_NITD_Jscale.remakeKPPnModel.II.Kas18.FullRun/'
    RunDict['Shah22'] = '/users/ts551/scratch/GC/rundirs/P_ARNA/gc_4x5_47L_merra2_fullchem_aciduptake.v13.4.1.ARNAv10.ARNA.Test_NITD_Jscale.remakeKPPnModel.II.Shah22.FullRun/'

    RunDict['Ye17'] = '/users/ts551/scratch/GC/rundirs/P_ARNA/gc_4x5_47L_merra2_fullchem_aciduptake.v13.4.1.ARNAv12.Ye17.FullRun/'
    #'/users/ts551/scratch/GC/rundirs/P_ARNA/gc_4x5_47L_merra2_fullchem_aciduptake.v13.4.1.ARNAv11.Ye17.1Hour'
    RunDict['Andersen22a'] = '/users/ts551/scratch/GC/rundirs/P_ARNA/gc_4x5_47L_merra2_fullchem_aciduptake.v13.4.1.ARNAv12.Andersen22a.FullRun/'
    RunDict['Andersen22b'] = '/users/ts551/scratch/GC/rundirs/P_ARNA/gc_4x5_47L_merra2_fullchem_aciduptake.v13.4.1.ARNAv12.Andersen22b.FullRun/'

    # Point to the OutputDir
    for key in RunDict.keys():
        RunDict[key] = RunDict[key] + '/OutputDir/'

    # - local variables
    #
    TotalJNITeffVar = 'JscaleEffNIT-all'
    JscaleNITvar = 'JscaleNIT'
    JscaleNITsvar = 'JscaleNITs'
    JscaleNITD1var = 'JscaleNITD1'
    JscaleNITD2var = 'JscaleNITD2'
    JscaleNITD3var = 'JscaleNITD3'
    JscaleNITD4var = 'JscaleNITD3'
    JNITonDust =  ['Andersen22a', 'Andersen22b']
    #
    NO3Var = 'NIT-all'
    SCprefix = 'SpeciesConc_'
    dsD = {}
    StateMetD = {}

    # - Loop and extract weighted Jscale for all runs at once.
#    key = 'Kas18'
    for key in RunDict.keys():

        # Retrieve folder with model output
        folder = RunDict[key]

        # Set dates2use on a per key basis.
#         if (key == 'Ye17'):
#             # Ye17 runs not progressing past ~3 months
#             sdate = datetime.datetime(2019, 1, 30)
#             edate = datetime.datetime(2019, 10, 30)
# #        elif (key == 'Andersen22a') or (key == 'Andersen22b'):
#         elif (key == 'Andersen22b'):
#             # Andersen22a/b failed ~6 months post spin up (use period for all)
#             sdate = datetime.datetime(2019, 1, 1)
#             edate = datetime.datetime(2019, 6, 30)
#         else:
#             sdate = datetime.datetime(2019, 1, 1)
#             edate = datetime.datetime(2019, 12, 30)
        sdate = datetime.datetime(2019, 1, 1)
        edate = datetime.datetime(2019, 12, 30)
        dates2use = pd.date_range(sdate, edate, freq='1D')
#        dates2use = None

        # Get data and use average value over time for now
        ds = AC.GetSpeciesConcDataset(wd=folder, dates2use=dates2use)
        ds = ds.mean(dim='time')
        # Add a total nitrate variable
        ds = AC.AddChemicalFamily2Dataset(ds, fam=NO3Var)
        # Get the Jrates
        dsJ = AC.GetJValuesDataset(wd=folder, dates2use=dates2use)
        dsJ = dsJ.mean(dim='time')

        # Setup ancillary data for conversions
        MolecVar = 'Met_MOLECS'
        StateMet = AC.get_StateMet_ds(wd=folder, dates2use=dates2use)
        StateMet = AC.add_molec_den2ds(StateMet, MolecVar=MolecVar)
        StateMetD[key] = StateMet
#        AVG = AC.constants( 'AVG' )

        # Jscale on NIT? * [NIT?]
        dsJ[JscaleNITvar] = dsJ['Jval_NIT'] / dsJ['Jval_HNO3']
        dsJ[JscaleNITsvar] = dsJ['Jval_NITs'] / dsJ['Jval_HNO3']
        if (key in JNITonDust):
            dsJ[JscaleNITD1var] = dsJ['Jval_NITD1'] / dsJ['Jval_HNO3']
            dsJ[JscaleNITD2var] = dsJ['Jval_NITD2'] / dsJ['Jval_HNO3']
            dsJ[JscaleNITD3var] = dsJ['Jval_NITD3'] / dsJ['Jval_HNO3']
            dsJ[JscaleNITD4var] = dsJ['Jval_NITD4'] / dsJ['Jval_HNO3']

        # Sum the NIT? Jscales * [NIT], then devide by total
        dsJ[TotalJNITeffVar] = dsJ[JscaleNITvar] * ds[SCprefix+'NIT']
        dsJ[TotalJNITeffVar] += dsJ[JscaleNITsvar] * ds[SCprefix+'NITs']
        if (key in JNITonDust):
            dsJ[TotalJNITeffVar] += dsJ[JscaleNITD1var] * ds[SCprefix+'NITD1']
            dsJ[TotalJNITeffVar] += dsJ[JscaleNITD2var] * ds[SCprefix+'NITD2']
            dsJ[TotalJNITeffVar] += dsJ[JscaleNITD3var] * ds[SCprefix+'NITD3']
            dsJ[TotalJNITeffVar] += dsJ[JscaleNITD4var] * ds[SCprefix+'NITD4']
        dsJ[TotalJNITeffVar] = dsJ[TotalJNITeffVar] / ds[SCprefix+NO3Var]

        # Save to a dictionary of datasets to plot
        dsD[key] = dsJ[[TotalJNITeffVar, JscaleNITvar, JscaleNITsvar]]

        del ds, dsJ

    # --- Plot up values weighted JScale values
    # elephants
    PltAsLog = True
    # Setup a PDF
    savetitle = 'ARNA_spatial_weighted_Jscale'
    if PltAsLog:
        savetitle += '_log'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # Plot up values at the surface
#    keys = ['Kas18', 'Shah22', 'Ye17', 'Andersen22a', 'Andersen22b']
    for key in dsD.keys():

        # Select data to use in plotting
        var2plot = TotalJNITeffVar
        ds2plot = dsD[key].copy().sel(lev=dsD[key].lev[0])

        # normalise the plot between vmin and vmax
        if PltAsLog:
            vmin = 1.0
            vmax = 200
            norm = LogNorm(vmin=vmin, vmax=vmax)
#            kwargs = {'vmin':vmin, 'vmax': vmax, 'norm': norm}
            kwargs = {'norm': norm}
        else:
            kwargs = {}

        AC.quick_map_plot(ds2plot, var2plot=var2plot, verbose=verbose,
                          **kwargs)
        TitleStr = "Average surface weighted JScale ('{}')"
        plt.title(TitleStr.format(key))

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
#        if show_plot:
#            plt.show()
        plt.close()
    # Do some memory management...
    gc.collect()


    # --- Plot up values zonally
    # Zonal plot of these values too.
    for key in dsD.keys():

        # Select data to use in plotting
        var2plot = TotalJNITeffVar
        ds2plot = dsD[key].copy()
        StateMet2use = StateMetD[key].copy().mean(dim='time')

        fig, ax = plt.subplots()
        if PltAsLog:
            vmin = 1.0
            vmax = 200
            norm = LogNorm(vmin=vmin, vmax=vmax)
#            kwargs = {'vmin':vmin, 'vmax': vmax,
            kwargs = {'norm': norm}
        else:
            kwargs = {}
        im = AC.ds2zonal_plot(ds2plot, var2plot=var2plot,
                              StateMet=StateMet2use,
                              fig=fig, ax=ax, **kwargs)
        TitleStr = "Zonal weighted JScale  ('{}')"
        plt.title(TitleStr.format(key))

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



def Check_surface_values4keyJscale_inputs(RunDict=None):
    """
    Check the surface values for NIT species just to striation seen in output
    """
    #
    for key in RunDict.keys():

        # Retrieve folder with model output
        folder = RunDict[key]

        # Set dates2use on a per key basis.
        if (key == 'Ye17'):
            # Ye17 runs not progressing past ~3 months
            sdate = datetime.datetime(2018, 9, 30)
            edate = datetime.datetime(2019, 10, 30)
        elif (key == 'Andersen22a') or (key == 'Andersen22b'):
            # Andersen22a/b failed ~6 months post spin up (use period for all)
            sdate = datetime.datetime(2019, 1, 1)
            edate = datetime.datetime(2019, 6, 30)
        else:
            sdate = datetime.datetime(2019, 1, 1)
            edate = datetime.datetime(2019, 12, 30)
        dates2use = pd.date_range(sdate, edate, freq='1D')
#        dates2use = None

        # Get data and use average value over time for now
        ds = AC.GetSpeciesConcDataset(wd=folder, dates2use=dates2use)
#        ds = ds.mean(dim='time')
        # Add a total nitrate variable
        ds = AC.AddChemicalFamily2Dataset(ds, fam=NO3Var)

    #
    vars2plot = ['SO4', 'SO4D1', 'SO4D2', 'SO4D3' ,'SO4D4',
                'NIT', 'NITs', 'NIT-all', 'NITD1', 'NITD2', 'NITD3', 'NITD4',]
    vars2plot = [SCprefix+i for i in vars2plot]
    # Setup a PDF
    savetitle = 'ARNA_spatial_striation_check_{}_II'.format(key)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # Plot up values at the surface
#    keys = ['Kas18', 'Shah22', 'Ye17', 'Andersen22a', 'Andersen22b']
#    key2plot = 'Kas18'
    for var2plot in vars2plot:

        for dt in ds.time:

            print(var2plot, dt)
            ReadableDate = AC.dt64_2_dt([dt.values])[0].strftime("%Y-%m-%d")

            # Select data to use in plotting
            ds2plot = dsD[key].copy().sel(lev=ds.lev[0])

            # Sub select the month
            ds2plot = ds2plot.sel(time=dt)

            kwargs = {}
            AC.quick_map_plot(ds2plot, var2plot=var2plot, verbose=verbose,
                              **kwargs)
            TitleStr = "Average surface value ('{}') on {}"
            plt.title(TitleStr.format(var2plot, ReadableDate))

            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
    #        if show_plot:
    #            plt.show()
            plt.close()

    # Do some memory management...
    gc.collect()

    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')



def check_nitrate_conc_by_month():
    """
    Check the nitrate concentration by month
    """

    # RunDict

    # Months to use
    dates2use = None
    sdate = datetime.datetime(2018, 6, 1)
    edate = datetime.datetime(2019, 12, 30)
    dates2use = pd.date_range(sdate, edate, freq='1D')

    # Other local variables
    SCprefix = 'SpeciesConc_'
    specs2use = ['O3', 'CO',  'NIT', 'NITs', 'NIT-all', 'NOx' ]
    vars2use = [SCprefix+i for i in specs2use]
    extra_burden_specs = vars2use
    extra_surface_specs = vars2use
    fam = 'NIT-all'
    use_time_in_trop = True
    rm_strat = True
    avg_over_time = False

    # For each month, get summary stats where available
    DataD = {}
    dfs = {}
    for key in RunDict.keys()

        # Get Species concentrations
        ds = GetSpeciesConcDataset(wd=RunDict[key])#, dates2use=dates2use)
        # Add nitrate family
        ds = AC.AddChemicalFamily2Dataset(ds, fam='NIT-all', prefix=SCprefix)
        ds = AC.AddChemicalFamily2Dataset(ds, fam='NOx', prefix=SCprefix)
        # Get the burdens
        StateMet = AC.get_StateMet_ds(wd=RunDict[key])#, dates2use=dates2use)
        S = AC.get_Gg_trop_burden(ds, vars2use=vars2use, StateMet=StateMet,
                               use_time_in_trop=use_time_in_trop,
                               avg_over_time=avg_over_time,
                               rm_strat=rm_strat,
                               sum_spatially=False,
                               debug=debug)
        #
        S = S.sum(['lat', 'lev', 'lon'])


        # Select a single run to attempt data extraction for
        folder = RunDict[key]

        # try to get data for each date
        for date in dates2use:
            try:
                __df = AC.get_stats4RunDict_as_df( RunDict={key:folder},
                          extra_burden_specs=extra_burden_specs,
                          extra_surface_specs=extra_surface_specs,
                          dates2use=[date])
                # Save the data
                dfs[date] = __df

            except:
                PrtStr = "WARNING: skipping date for '{}' ({})"
                print( PrtStr.format(key, date) )

        # Save the data where available for the dates
        df = pd.DataFrame()


        DataD[key]

    # - Summarise the numbers


    # -









def explore_print_statement_debug_of_aciduptake_photolysis():
    """
    Explore the "print statement" debugging...
    """
    # Data to use and its location?
    folder = ar.get_local_folder('RunRoot')
    folder += 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.'
    folder += 'DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags.TEMP_TEST/'
#    filename = 'geos.log.Andersen.v006'
#    filename = 'geos.log.Andersenv007'
    filename = 'geos.log.Andersenv008'

    def __convert2float(x):
        """ helper function to process print statement ouput """
        try:
            return float(x)
        except ValueError:
            return np.NaN

    # Get data from log file and store in data DataFrame
    df = pd.DataFrame()
    with open(folder+filename) as f:
        lines = [i for i in f]
        # Get all nitrate conc (molec/cm3)
        PrefixStr = ' @@@ NO3_CONC II (molec/cm3) - NIT(s)+NITD1-4'
        NIT_all_molec_cm3 = 'NIT-all (molec/cm3)'
        data = [i for i in lines if PrefixStr in i]
        data = [float(i.split(PrefixStr)[-1].strip()) for i in data]
        df[NIT_all_molec_cm3] = data
        # Get all NIT+NITS conc (molec/cm3)
        PrefixStr = '@@@ NO3_CONC I (molec/cm3) - NIT+NITs'
        NITs_molec_cm3 = 'NIT+NITs (molec/cm3)'
        data = [i for i in lines if PrefixStr in i]
        data = [float(i.split(PrefixStr)[-1].strip()) for i in data]
        df[NITs_molec_cm3] = data
        # Get all nites in nM/m3
        PrefixStr = '@@@ NO3_CONC IV (nM m^-3)'
        NITs_mN_m3 = 'NIT-all (nM/m3)'
        data = [i for i in lines if PrefixStr in i]
        data = [float(i.split(PrefixStr)[-1].strip()) for i in data]
        df[NITs_mN_m3] = data
        # Get JScale
        PrefixStr = '@@@ Jscale '
        JScale = 'JScale'
        data = [i for i in lines if PrefixStr in i]
        data = [float(i.split(PrefixStr)[-1].strip()) for i in data]
        df[JScale] = data
        # Get JNIT and JNITs
        PrefixStr = '@@@ JscaleNIT/s'
        JScaleNIT = 'JScale NIT'
        JScaleNITs = 'JScale NITs'
        data = [i for i in lines if PrefixStr in i]
        dataA = [i.split(PrefixStr)[-1].strip().split(' ')[-1] for i in data]
        dataA = [__convert2float(i) for i in dataA]
        dataB = [i.split(PrefixStr)[-1].strip().split(' ')[0] for i in data]
        dataB = [__convert2float(i) for i in dataB]
        df[JScaleNIT] = dataA
        df[JScaleNITs] = dataB
        # Get Channels
        PrefixStr = '@@@ JNITChanA/B'
        ChanA = 'Frac HONO'
        ChanB = 'Frac NO2'
        data = [i for i in lines if PrefixStr in i]
        dataA = [i.split(PrefixStr)[-1].strip().split(' ')[-1] for i in data]
        dataA = [__convert2float(i) for i in dataA]
        dataB = [i.split(PrefixStr)[-1].strip().split(' ')[0] for i in data]
        dataB = [__convert2float(i) for i in dataB]
        df[ChanA] = dataA
        df[ChanB] = dataB
    print(df)

    for col in df.columns:
        print(col)
        print(df[col].describe())
    df_BACKUP = df.copy()

    # Only show where nitrate >10 pptv
    pptv_in_molec_cm3 = 255440.55029432118  # '2.6E+05'
    __df = df.loc[df[NIT_all_molec_cm3] > pptv_in_molec_cm3, :]
    for col in __df.columns:
        print(col)
        print(__df[col].describe())

    # Only show where nitrate >100 pptv
    pptv_in_molec_cm3 = 2554405.502943212  # '2.6E+06'
    __df = df.loc[df[NIT_all_molec_cm3] > pptv_in_molec_cm3, :]
    for col in __df.columns:
        print(col)
        print(__df[col].describe())

    # Only show where nitrate >500 pptv
    pptv_in_molec_cm3 = 12772027.514716059  # '1.3E+07'
    __df = df.loc[df[NIT_all_molec_cm3] > pptv_in_molec_cm3, :]
    for col in __df.columns:
        print(col)
        print(__df[col].describe())

    # Only show where nitrate >1000 pptv
    pptv_in_molec_cm3 = 25543799.588881828  # '2.6E+07'
    __df = df.loc[df[NIT_all_molec_cm3] > pptv_in_molec_cm3, :]
    for col in __df.columns:
        print(col)
        print(__df[col].describe())


if __name__ == "__main__":
    main()

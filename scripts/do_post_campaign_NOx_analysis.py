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


def main():
    """
    Main driver function
    """
    # Get the core FAAM data
#    ar.get_FAAM_core4flightnum()
    # Get the ToF-CIMS data
    ds2 = ar.get_CIMS_data4flight()
    # Explore NOy/NIT after inclusion of acid uptake
    explore_NOy_with_acid_uptake()

    #
    explore_ARNA_period_with_acid_uptake()


def explore_NOy_with_acid_uptake():
    """
    Explore the global NOy budgets with acidic uptake to dust included
    """
    from AC_tools import species_mass
    # - Local variables
    # use a v12.9.1 compatible definition of NOx
#    NOySpecs = AC.GC_var('NOy')
    NOySpecs = [
        'NO', 'NO2', 'PAN', 'HNO3',
        #    'PMN', 'PPN', 'R4N2',
        'N2O5', 'HNO4',\
        'BrNO2', 'BrNO3',
        #    'MPN',
        #    'ISOPN', 'PROPNN',
        # 'MMN',\
        'NO3', 'HNO2', 'IONO', 'IONO2', 'INO', 'ClNO2', 'ClNO3'
    ]
    # include NITs
    NOySpecs += ['NIT', 'NITs']
    # Also consider dust uptake on NOx
    NOySpecsA = NOySpecs + ['NITD1', 'NITD2', 'NITD3', 'NITD4']
    # Set runs to use
    RunRoot = '/users/ts551/scratch/GC/rundirs/'
    RunStr = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020{}/'
    run_dict = {
        'BASE': RunRoot+RunStr.format('.BCs.repeat.III/spin_up/'),
        #    'ACID.III': RunRoot + RunStr.format('.DustUptake.III/spin_up/'),
        'ACID.IV': RunRoot+RunStr.format('.DustUptake.IV/OutputDir/'),
        'JNIT': RunRoot+RunStr.format('.DustUptake.IV.JNIT/OutputDir/'),
        'JNITx25': RunRoot+RunStr.format('.DustUptake.IV.JNIT.x25/OutputDir/'),
    }
#    run_dict = d
    # - Analysis

    # Get generic stats on runs
    dates2use = None
    extra_specs = NOySpecs + ['SO2', 'SO4']
    df = AC.get_general_stats4run_dict_as_df(run_dict=d,
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

    RunDir = '/users/ts551/scratch/GC/rundirs/'
    BASEprefix = 'geosfp_4x5_standard.v12.9.0.BASE.2019.2020.'
#     geosfp_4x5_standard.v12.9.0.BASE.2019.2020./OutputDir/

    ACIDprefix = 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.'
    ACIDprefix += 'DustUptake.'
    REF1 = 'BASE.BC'
    d = {
#    'JNITS': RunDir+BASEprefix+'ARNA.Nest.repeat.IV.JNIT.x25.Dust/',
#    'GFASx2': RunDir+BASEprefix+'ARNA.Nest.repeat.GFASx2/',
#    'JNITSx25':RunDir+BASEprefix+'ARNA.Nest.repeat.JNITs.x25/',
#    'ACID.BC': RunDir+ACIDprefix+'ARNA.DustUptake.BCs/', # failed run... re-start
    'BASE.BC': RunDir+BASEprefix+'ARNA.BCs.repeat/',
    'BASE.GFASx2.BC': RunDir+BASEprefix+'ARNA.GFASx2.BCs/',
    'ACID.BC': RunDir+ACIDprefix+'BCs/',
    #  'ACID.JNIT.BC' only has partial NetCDF output (for the campaign period)

#   'ACID.JNIT.BC': RunDir+ACIDprefix+'JNIT.BCs/',
    'ACID.JNITx25.BC': RunDir+ACIDprefix+'JNITx25.BCs/',
    'ACID.BC.Isotherm': RunDir+ACIDprefix+'JNIT.Isotherm.BCs.repeat.ON.II/',
    }
    #
    for key in d.keys():
        d[key] = d[key]+'OutputDir/'
#        d[key] = d[key]+'spin_up/'

    # Use spin up year - just for testing
    sdate = datetime.datetime(2018, 1, 1, )
    edate = datetime.datetime(2019, 1, 1, )
    # Use spun up year before the campaign (2019)
#    sdate = datetime.datetime(2019, 1, 1, )
#    edate = datetime.datetime(2020, 1, 1, )
    # Just use the dates of the ARNA campaign
#    sdate = datetime.datetime(2020, 1, 1, )
#    edate = datetime.datetime(2020, 3, 1, )
    # Just use the dates of the ARNA campaign - but in 2019
    sdate = datetime.datetime(2019, 1, 1, )
    edate = datetime.datetime(2019, 3, 1, )
    dates2use = pd.date_range(sdate, edate, freq='1H')

    # Get Key statistics for run
    ExSpecs =  [
    'HNO2', 'HNO3', 'NIT', 'NITs',
#    'NOy'
#    'NOx',
#    'NIT-all', 'NOy-gas','NOy',
    ]
    df = AC.get_general_stats4run_dict_as_df(run_dict=d,
                                             REF_wd=d[REF1],
                                             extra_surface_specs=ExSpecs,
                                             extra_burden_specs=ExSpecs,
                                             dates2use=dates2use,
                                             use_REF_wd4Met=True,
                                             debug=debug)

    # Calculate lightning source and add to pd.DataFrame
    var2use = 'EmisNO_Lightning'
    varName = 'Lightning (Tg N/yr)'
    dsD = {}
    for key in d.keys():
        dsD[key] = AC.get_HEMCO_diags_as_ds(wd=d[key],
                                            dates2use=dates2use)
    df2 = pd.DataFrame()
    for key in d.keys():
        ds = dsD[key]
        val = (ds[var2use].mean(dim='time').sum(dim='lev') * ds['AREA'] )
        val2 = val.values.sum() * 60 * 60 * 24 * 365 # => /yr
        df2.loc[key,varName] = val2*1E3/1E12
    # add into main DataFrame
    df = df.T
    df = pd.concat([df, df2], axis=1)
    df = df.T

    # Rename the keys to more readable names
    rename_dict = {
    'BASE.BC': 'BASE',
    'BASE.GFASx2.BC': 'BASE.BBx2',
    'ACID.BC': 'ACID',
    'ACID.JNITx25.BC': 'ACID.JNITx25'
    }
    for key in list(d.keys()):
        print(key, (key in rename_dict.keys()))
        if (key in rename_dict.keys()):
            NewKey = rename_dict[key]
            d[NewKey] = d[key]
            d.pop(key)

    # Get the model for all species by run
    prefix = 'SpeciesConc_'
    dsD = {}
    for key in d.keys():
        ds = AC.get_GEOSChem_files_as_ds(wd=d[key], dates2use=dates2use)
        ds = add_derived_GEOSChem_specs2ds(ds, prefix=prefix)
        dsD[key] = ds

#    StateMet = AC.get_StateMet_ds(wd=Folder)
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
            'NOy', 'HNO3', 'NIT-all','NOy-gas',
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
        elif spec == 'HNO2':
            plt.xlim(-0.001, 0.003*1E3)
        elif spec == 'HNO3':
            plt.xlim(-0.05, 0.5)
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

    # NO2 / NOx
#     spec = 'NO2:NOx'
#     for key in dsD.keys():
#         ds = dsD[key].copy()
#         dsI, units = __get_data4plotting(ds, spec='NO2')
#         dsII, units = __get_data4plotting(ds, spec='NOx')
#         ds = dsI / dsII
#         #  plot up...
#         ds.plot(y='lev', label=key)
#         __beautify_plot(spec=spec, units=units)
#     # Save to PDF
#     AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
#     if show_plot:
#         plt.show()
#     plt.close()

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
    specs2plot_AsLog = ['NO', 'NO2', 'NOx',]# 'HNO2']
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


    # plot up the difference between the runs
    REF = 'BASE'
    DIFF = 'ACID.BC.Isotherm'
    prefix = 'SpeciesConc_'
    vars2plot = ['O3', 'CO', 'NOx', 'NO2', 'NO', 'HNO2',
                'NOy', 'HNO3', 'NIT-all','NOy-gas',
                ]
    vars2plot = [prefix+i for i in vars2plot]
    plot_up_surface_diff_between_runs(dsD, REF=REF, DIFF=DIFF,
                                         vars2plot=vars2plot,
                                         prefix=prefix)




def tag_GC_simulations():
    """
    Setup tagging of NOx/HNO2 prod etc in KPP
    """
    # Diagnostics to use?
    diags = [
    'ProdHNO2fromHvNIT', 'ProdHNO2fromHvNITs', 'ProdHNO2fromHvNITD1',
    'ProdHNO2fromHvNITD2', 'ProdHNO2fromHvNITD3', 'ProdHNO2fromHvNITD4',
    'ProdNO2fromHvNIT', 'ProdNO2fromHvNITs', 'ProdNO2fromHvNITD1',
    'ProdNO2fromHvNITD2', 'ProdNO2fromHvNITD3', 'ProdNO2fromHvNITD4',
    'ProdNO2fromHONO', 'ProdHNO2fromOHandNO', 'ProdHNO2fromHET'
    ]
    prefix = 'TN{:0>3}'
    tags = [prefix.format(i+1) for i in range(len(diags))]
    # pair up numbering (so that runs with different diagnostics have same #s)?
    d = dict(zip(diags, tags))
    for key in d.keys():
        print('{} : {};'.format(key, d[key]) )
    # Also print out just using "P" as the prefix.
    for key in d.keys():
        print('P{} : {};'.format(d[key], d[key]) )



    # prepare other output for GEOS-Chem input files
    extr_str = 'ARNA_Standard'
    AC.print_out_lines_for_gckpp_file(tags=tags, extr_str=extr_str)
    AC.prt_lines4species_database_yml(tags=tags, extr_str=extr_str)

    ptr_str = '{:<11}= IGNORE; {}'
    d = dict(zip(diags, tags))
    for key in d.keys():
        print(ptr_str.format(d[key], '{'+key+'}' ) )






if __name__ == "__main__":
    main()

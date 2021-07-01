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
        'BASE': RunRoot + RunStr.format('.BCs.repeat.III/spin_up/'),
        #    'ACID.III': RunRoot + RunStr.format('.DustUptake.III/spin_up/'),
        'ACID.IV': RunRoot + RunStr.format('.DustUptake.IV/OutputDir/'),
        'JNIT': RunRoot + RunStr.format('.DustUptake.IV.JNIT/OutputDir/'),
        'JNITx25': RunRoot + RunStr.format('.DustUptake.IV.JNIT.x25/OutputDir/'),
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


if __name__ == "__main__":
    main()

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



    #
    explore_NOy_with_acid_uptake()


    #



def explore_NOy_with_acid_uptake():
    """
    Explore the NOy budgets with acidic uptake to dust included
    """
    from AC_tools import species_mass
    # - Local variables
    #
#    NOySpecs = AC.GC_var('NOy')
    NOySpecs = [
    'NO', 'NO2', 'PAN', 'HNO3',
#    'PMN', 'PPN', 'R4N2',
    'N2O5', 'HNO4',\
    'BrNO2', 'BrNO3',
#    'MPN',
#    'ISOPN', 'PROPNN',
    #'MMN',\
    'NO3', 'HNO2', 'IONO', 'IONO2', 'INO', 'ClNO2', 'ClNO3'
    ]
    # include NITs
    NOySpecs += ['NIT', 'NITs']

    #
    NOySpecsA = NOySpecs + ['NITD1', 'NITD2', 'NITD3', 'NITD4']

    # Set runs to use
    RunRoot = '/users/ts551/scratch/GC/rundirs/'
    RunStr = 'merra2_4x5_standard.v12.9.0.BASE.2019.2020{}/'
    d = {
    'BASE': RunRoot + RunStr.format('.BCs.repeat.III/spin_up/'),
#    'ACID.III': RunRoot + RunStr.format('.DustUptake.III/spin_up/'),
    'ACID.IV': RunRoot + RunStr.format('.DustUptake.IV/OutputDir/'),
    'JNIT': RunRoot + RunStr.format('.DustUptake.IV.JNIT/OutputDir/'),
    'JNITx25': RunRoot + RunStr.format('.DustUptake.IV.JNIT.x25/OutputDir/'),
    }
    run_dict = d
    # - Analysis

    # Get generic stats on runs
    dates2use = None
    extra_burden_specs = NOySpecs + ['SO2', 'SO4']
    df = AC.get_general_stats4run_dict_as_df(run_dict=d,
                                             dates2use=dates2use,
                                         extra_burden_specs=extra_burden_specs)

    # get NOy
    avg_over_time = True # Note: burdens area averaged overtime
    prefix = 'SpeciesConc_'
    dates2use = None
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
    keys2use = [ 'BASE', 'ACID.III', 'ACID.IV',][::-1]
    keys2use = ['BASE', 'ACID.IV', 'JNIT', 'JNITx25'][::-1]
    for key in keys2use:
#    for key in d.keys():
        print(key, d[key])
        # Use the appropriate list of species
        if ('ACID' in key) or ('NIT' in key):
            specs2use = NOySpecsA
        else:
            specs2use = NOySpecs
        vars2use = [prefix+i for i in specs2use]
             # Get StateMet object for 1st of the runs and use this for all runs
#            if use_REF_wd4Met:
                # Set working directory for shared variables
#                if isinstance(REF_wd, type(None)):
#                    REF_wd = run_dict[ list(run_dict.keys())[0] ]
#                StateMet = get_StateMet_ds(wd=REF_wd, dates2use=dates2use)
#            else:
        StateMet = AC.get_StateMet_ds(wd=run_dict[key], dates2use=dates2use)
        # Average burden over time
        ds = dsD[key]#.mean(dim='time', keep_attrs=True)

        # reduce datasets to region of analsysi


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
    df = df.T
    print(df[keys2use])
    # Add percent difference
    REF = 'BASE'
    for col in [i for i in df.columns if i != REF]:
        df[col] = (df[col] - df[REF]) / df[REF]*100
    print(df[keys2use])




if __name__ == "__main__":
    main()
#!/usr/bin/python
"""
Driver for analysis of "4 pptv HONO world" following the ARNA campaign
"""
import arna as ar
import sys


def main():
    """
    Main driver for 4 pptv HONO world analysis
    """
    # Do the core stats for the runs
    do_analysis_of_4pptv_HONO_world()

    # plot up key changes
    plt_spatial_changes_in_4pptv_HONO_world()


def do_analysis_of_4pptv_HONO_world():
    """
    Do analysis of a model world where HONO is a minimum of 4 pptv
    """
    # Model runs to use?
    RunSet = 'PostHONO'
    GC_version = 'v13.4'
    RunDict = get_RunDict_of_HONO_runs(RunSet=RunSet, GC_version=GC_version)
    # Temporally use one run for the references values (e.g. statemet)
    use_REF_wd4Met = True
    REF_wd = RunDict['min4pptvHONO']

    # Set dates to use (as model run progresses)
    dates2use = [
    datetime.datetime(2018, 1, 1 ),
    datetime.datetime(2018, 2, 1 ), # Initial month of spin up
    ]
    # Use spin up year?
#    dates2use = [datetime.datetime(2018, 1+i, 1) for i in range(12)]

    # -- Stats on model runs
    # Get generic stats
    extra_burden_specs = ['NOx', 'NIT-all', 'HNO2',  'NOy', 'HNO3',  ]
    extra_surface_specs = extra_burden_specs
#    df = AC.get_general_stats4run_dict_as_df(run_dict=RunDict,
    df = get_stats4RunDict_as_df(RunDict=RunDict,
                                             dates2use=dates2use,
                                             use_REF_wd4Met=use_REF_wd4Met,
                                             REF_wd=REF_wd,
                                         extra_burden_specs=extra_burden_specs,
                                     extra_surface_specs=extra_surface_specs,
                                             verbose=verbose,
                                             debug=debug,
                                             )
    # Get extra stats relate to HONO
    DiagVars = [
    'CO-lifetime', 'OH-lifetime', 'NOx-lifetime', 'NO2-lifetime',
    'Ox-lifetime', 'HNO2-lifetime',
    'HNO2:NOx', 'HNO2:HNO3',
    'HNO3:NOx', 'HNO3:NO2', 'HNO3:NIT',
    'NIT:NOx', 'NO2:NO',
    ]
    df2 = get_stats4RunDict_as_df(RunDict=RunDict,
                                             dates2use=dates2use,
                                             use_REF_wd4Met=use_REF_wd4Met,
                                             REF_wd=REF_wd,
                                             DiagVars=DiagVars,
                                             verbose=verbose,
                                             debug=debug,
                                             )
    # CO-lifetime

    # OH-lifetime

    # OH concentration

    # HNO3:NIT

    # HNO3:NO2

    # HNO3:NOx

    # HONO:NOx

    # HONO:HNO3

    # NIT:NOx

    # NO2:NO

    # Save
    SaveName = 'ARNA_Key_Stats_on_model_runs_post_HONO'
    df.to_csv('{}{}'.format(SaveName))

    # Save out core stats for manuscript
    corestats = []
    df = df[corestats].to_csv(SaveName)



def plt_spatial_changes_in_4pptv_HONO_world():
    """
    Plot up changes in key metrics spatially and zonally
    """
    # Model runs to use?
    RunDict = get_RunDict_of_HONO_runs()
    # Set dates to use (as model run progresses)
    dates2use = [
    datetime.datetime(2018, 1, 1 ),
    datetime.datetime(2018, 2, 1 ),
    ]
    # Get Species conc
    dsD = {}
    for key in RunDict.keys():
        dsD[key] = AC.GetSpeciesConcDataset(wd=RunDict[key],
                                            dates2use=dates2use)
    # Add families of interest
    families2use = ['NOx', 'NIT-all', 'NOy']
    for key in dsD.keys():
        ds = dsD[key]
        for fam in families2use:
            ds = AC.AddChemicalFamily2Dataset(ds, fam=fam)
        dsD[key] = ds

    # rename variables in dataset
    prefix = 'SpeciesConc_'
    OldNames = [i for i in ds.data_vars if (prefix in i )]
    NewNames = [ i.split(prefix)[-1] for i in OldNames ]
    for key in dsD.keys():
        ds = dsD[key]
        ds = ds.rename( name_dict=dict(zip(OldNames, NewNames)) )
        dsD[key] = ds

    # Get MetState object


    # -- Graphics on model runs
    pcent = True

    REF1 = 'Base'
    DIFF = '4pptHONO'
    savetitle = 'ARNA_spatial_HONO4pptv'
    specs2plot = [ 'O3', 'HNO2',  'HNO3',  ] + families2use

    # Updates setings for percentage plotting
    if pcent:
        kwargs = {'vmin': -100, 'vmax': 100}
        savetitle += '_pcent'
    else:
        kwargs = {}

    # Plot up HONO differences (surface)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    for var2plot in specs2plot:

        #
        ds1 = dsD[REF1][specs2plot].copy().mean(dim='time')
        ds2 = dsD[DIFF][specs2plot].copy().mean(dim='time')

        # Calculate difference
        if pcent:
            ds2plot = (ds2 - ds1) / ds1 *100
        else:
            ds2plot = ds2 - ds1

        # Select the surface
        ds2plot = ds2plot.sel( lev=ds2plot.lev.values[0] )

        # Units and scaling for plot?
        if pcent:
            scaleby = 1
            units = '%'
        else:
            units, scaleby = AC.tra_unit(var2plot, scale=True)

        TitleStr = "Diff of {} ('{}' vs. '{}') in '{}'"
        title = TitleStr.format(var2plot, DIFF, REF1, units)

        AC.quick_map_plot(ds2plot[[var2plot]]*scaleby, var2plot=var2plot,
                          verbose=verbose, title=title, **kwargs)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()
    #        del ds2plot


    # zonal

    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def get_RunDict_of_HONO_runs(RunSet='PostHONO', GC_version='v13.4'):
    """
    Retrieve the model runs
    """
    # Which runs to use
    RunRoot = ar.get_local_folder('RunRoot')
    if (RunSet == 'PostHONO') and (GC_version == 'v13.4'):
        RunStr = 'gc_4x5_47L_geosfp_fullchem.v13.4.0-rc.2'
        RunDict = {
        'Base': '{}{}{}/'.format(RunRoot, RunStr,'.orig.1monthTest'),
        'min4pptvHONO': '{}{}{}/'.format(RunRoot, RunStr,
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
        BASEStr = ''
        RunStr = ''

    else:
        PrtStr = "Unkonwn RunSet ('{}') and GC verions ('{}')"
        print(PrtStr.format(RunSet, GC_version))
        sys.exit(0)
    # Include the output directory folder in directory strings
    for key in RunDict.keys():
        RunDict[key] = '{}{}/'.format(RunDict[key], 'OutputDir')
    return RunDict


from AC_tools import GetSpeciesConcDataset, AddChemicalFamily2Dataset, species_mass, get_ref_spec, get_Gg_trop_burden, get_StateMet_ds, tra_unit, get_avg_2D_conc_of_X_weighted_by_Y, GC_var, get_ProdLoss_ds


def get_stats4RunDict_as_df(
                            RunDict=None,
                            extra_str='',
                            REF1=None,
                            REF2=None, REF_wd=None, res='4x5',
                            trop_limit=True,
                            save2csv=True,
                            extra_burden_specs=[],
                            extra_surface_specs=[],
                            GC_version='v12.6.0',
                            DiagVars=[],
                            use_time_in_trop=True, rm_strat=True,
                            dates2use=None, round=3,
                            use_REF_wd4Met=False,
                            IncConcAfterChemDiags=True,
                            IncProdLossDiags=True,
                            verbose=False, debug=False):
    """
    Get various stats on a set of runs in a dictionary ({name: location})

    Parameters
    ----------
    RunDict (dict): dicionary of run names and locations
    use_REF_wd4Met (bool): use a reference working directory for shared values?
    REF_wd (str): directory to use to extract shared variables
    REF1 (str): name of (1st) run in dictionary to to % change calculations from
    REF2 (str): name of (2nd) run in dictionary to to % change calculations from
    prefix (str):  string to include as a prefix in saved csv's filename
    extra_str (str):  string to include as a suffx in saved csv's filename
    save2csv (bool): save dataframe as a csv file
    trop_limit (bool): limit analysis to the troposphere?
    trop_mask (nd.array): 3D or 4D boolean array where stratosphere is False
    rm_strat (bool): remove the stratospheric values
    extra_burden_specs (list): list of extra species to give trop. burden stats on
    extra_surface_specs (list): list of extra species to give surface conc. stats on
    res (str): resolution of the modul output (e.g. 4x5, 2x2.5, 0.125x0.125)
    round (int): number of decimal places to round dataframe too

    Returns
    -------
    (pd.DataFrame)
    """
    # - Define local variables
    # Mass unit scaling
    mass_scale = 1E3
    mass_unit = 'Tg'
    # Mixing ratio (v/v) scaling?
    ppbv_unit = 'ppbv'
    ppbv_scale = 1E9
    pptv_unit = 'pptv'
    pptv_scale = 1E12
    # Setup lists and check for families (NOy, NIT-all, Iy, Cly, Bry, ... )
    core_specs = ['O3',  'NO', 'NO2']
    prefix = 'SpeciesConc_'
    ALLSp = core_specs + extra_burden_specs + extra_surface_specs
    # Also add all families and specs in 'DiagVars' to this list
    if len(DiagVars) >= 1:
        RatioVars = [i for i in DiagVars if ':' in i]
        RatioVars = [i.split(':') for i in RatioVars]
        RatioVars = [item for sublist in RatioVars for item in sublist]
        LifetimeVars = [i for i in DiagVars if '-' in i]
        LifetimeVars = [i.split(':') for i in LifetimeVars]
        LifetimeVars = [item for sublist in LifetimeVars for item in sublist]
        ALLSp = ALLSp + RatioVars + LifetimeVars
    ALLSp = list(set(ALLSp))
    FamilyNames = GC_var('FamilyNames')
    families2use = []
    for FamilyName in FamilyNames:
        if FamilyName in ALLSp:
            families2use += [FamilyName]
    if verbose or debug:
        print(families2use)
    # Invert RunDict
    RunDict_r = {v: k for k, v in list(RunDict.items())}

    # Core dataframe for storing calculated stats on runs
    df = pd.DataFrame()
    # - Get core data required
    # Get all of the 'SpeciesConcs' for runs as list of datasets
    dsD = {}
    for key in RunDict.keys():
        ds = GetSpeciesConcDataset(wd=RunDict[key], dates2use=dates2use)
        # Add families to Dataset
        if len(families2use) >= 1:
            for fam in families2use:
                ds = AddChemicalFamily2Dataset(ds, fam=fam, prefix=prefix)
        dsD[key] = ds

    # Get the StateMet object(s)
    MolecVar = 'Met_MOLCES'
    dsS = {}
    if use_REF_wd4Met:
        # Set working directory for shared variables
        if isinstance(REF_wd, type(None)):
            REF_wd = RunDict[list(RunDict.keys())[0]]
        StateMet = get_StateMet_ds(wd=REF_wd, dates2use=dates2use)
        StateMet = AC.add_molec_den2ds(StateMet)
        dsS[RunDict_r[REF_wd]] = StateMet
    else:
        for key in RunDict.keys():
            StateMet = get_StateMet_ds(wd=RunDict[key], dates2use=dates2use)
            StateMet = AC.add_molec_den2ds(StateMet)
            dsS[key] = StateMet

    # - Get burdens for core species
    avg_over_time = True  # Note: burdens area averaged overtime
    specs2use = list(set(core_specs+['CO']+extra_burden_specs))
    vars2use = [prefix+i for i in specs2use]
    for key in RunDict.keys():
        # Get StateMet object
        if use_REF_wd4Met:
            StateMet = dsS[RunDict_r[REF_wd]]
        else:
            StateMet = dsS[key]
        # Average burden over time
        ds = dsD[key]  # .mean(dim='time', keep_attrs=True)
        S = get_Gg_trop_burden(ds, vars2use=vars2use, StateMet=StateMet,
                               use_time_in_trop=use_time_in_trop,
                               avg_over_time=avg_over_time,
                               rm_strat=rm_strat,
                               debug=debug)
        # convert to ref spec equivalent (e.g. N for NO2, C for ACET)
        for spec in specs2use:
            ref_spec = get_ref_spec(spec)
            val = S[prefix+spec]
            S[prefix+spec] = val/species_mass(spec)*species_mass(ref_spec)
        # Upate varnames
        varnames = ['{} burden ({})'.format(i, mass_unit) for i in specs2use]
        S = S.rename(index=dict(zip(list(S.index.values), varnames)))
        # Save the values for run to central DataFrame
        df[key] = S

    # Transpose dataframe
    df = df.T

    # Scale units
    for col_ in df.columns:
        if 'Tg' in col_:
            df.loc[:, col_] = df.loc[:, col_].values/mass_scale
    # Transpose back to variables as index
    df = df.T

    # - Add Ozone production and loss...
    if IncProdLossDiags:
        try:
            pass
        except KeyError:
            pass

    # - Add lifetime calculations for species
    PtrStr1 = "Calculating lifetime for diag ('{}') for '{}' variable"
    ErrStr1 = "Loss diagnostic not found ({}), skipped lifetime calc ('{}')"
    ErrStr2 = "Prod/Loss files not found ('{}'), skipped lifetime calc ('{}')"
    lifetimes2calc = [i for i in DiagVars if ('lifetime' in i.lower())]
    PLprefix = 'Loss'
    prefix = 'SpeciesConc_'
    if (len(lifetimes2calc) >= 1):
        for key in RunDict.keys():
            for var2calc in lifetimes2calc:
                species2calc = var2calc.split('-lifetime')[0]
                print(var2calc, species2calc)
                if verbose:
                    print(PtrStr1.format(var2calc, species2calc))

                # Find loss value from prod/loss diagnostic
                PLvar = '{}_{}'.format(PLprefix, species2calc)
                try:
                    ds = get_ProdLoss_ds(wd=RunDict[key], dates2use=dates2use)
                    ds = ds[[PLvar]]

                except KeyError:
                    # If variable not found, skip lifetime calculation
                    print(ErrStr1.format(PLvar, var2calc))

                except AssertionError:
                    # If files not found, skip lifetime calculation
                    print(ErrStr2.format(key, var2calc))

            # find species/family concentration

            # calculate the tropospheric lifetime using burden and loss rate
            pass

    # - Add Ratio calculations
    PtrStr = "Calculating ratio for diag ('{}') for '{}' vs '{}'"
    long_nameStr = "Dry mixing ratio of species '{}' ('{}':'{}')"
    ratios2calc = [i for i in DiagVars if (':' in i)]
    SpecsAndFams = []
    prefix = 'SpeciesConc_'
    if (len(ratios2calc) >= 1):
        for key in RunDict.keys():
            for var2calc in ratios2calc:
                var1 = '{}{}'.format(prefix, var2calc.split(':')[0] )
                var2 = '{}{}'.format(prefix, var2calc.split(':')[-1] )
                if verbose:
                    print(PtrStr.format(var2calc, var1, var2))
                ds = dsD[key]
                ds[var2calc] = ds[var1].copy()
                ds[var2calc] = ds[var1] / ds[var2]
                attrs = ds[var1].attrs
                attrs['long_name'] = long_nameStr.format(var2calc, var1, var2)
                dsD[key] = ds
                # Get StateMet object
                if use_REF_wd4Met:
                    StateMet = dsS[RunDict_r[REF_wd]]
                else:
                    StateMet = dsS[key]
                # Calculate molecular weighted values
                avg = ds[[var2calc]] * StateMet[MolecVar]
                # Only consider troposphere?
                if rm_strat:
                    # Loop by spec
                    if use_time_in_trop:
                        avg = rm_fractional_troposphere(avg,
                                                        vars2use=[var2calc],
                                                        StateMet=StateMet)
                    else:
                        avg = avg[var2calc].where(trop_mask)

                # Weight over lon?
#                dstemp.sum(dim=['lon']) / StateMet[MolecVar].sum(dim=['lon'])
                # Weight values and save to dataframe
                avg = avg.sum() / StateMet[MolecVar].sum()
                if debug:
                    print(avg, avg[var2calc])
                    print( avg[var2calc].values )
                avg =  avg[var2calc].values
                if debug:
                    print(avg)
                df.loc[var2calc, key] = avg


    # - Add lightning source if in HEMCO output
#     var2use = 'EmisNO_Lightning'
#     varName = 'Lightning (Tg N/yr)'
#     try:
#     # TODO -  add check for HEMCO NetCDF files in the output folder
# #    if True:
#         dsH = {}
#         for key in run_dict.keys():
#             dsH[key] = get_HEMCO_diags_as_ds(wd=run_dict[key])
#             ds = ds[key]
#             val = (ds[var2use].mean(dim='time').sum(dim='lev') * ds['AREA'] )
#             val2 = val.values.sum() * 60 * 60 * 24 * 365 # => /yr
#             df.loc[key,varName] = val2*1E3/1E12
#     except KeyError:
#         pass
#         df.loc[key,varName] = np.nan

    # - Surface concentrations
    specs2use = list(set(core_specs+['N2O5']+extra_surface_specs))
    prefix = 'SpeciesConc_'
    # Loop by run and get stats
    for key in RunDict.keys():
        if debug:
            print(key)
        ds = dsD[key].copy()
        # Select surface and average over time
        ds = ds.mean(dim='time')
        ds = ds.isel(lev=ds.lev == ds.lev[0])
        for spec in specs2use:
            # Get units and scaling
            units, scale = tra_unit(spec, scale=True)
            # Surface ozone
            varname = '{} surface ({})'.format(spec, units)
            var = prefix+spec
            # Save values on a per species basis to series
            val = get_avg_2D_conc_of_X_weighted_by_Y(ds, Xvar=var, Yvar='AREA')
            # Save calculated values to dataframe
            df.loc[varname, key] = val

    # Transpose dataframe
    df = df.T
    # Scale units
    for col_ in df.columns:
        if 'ppb' in col_:
            df.loc[:, col_] = df.loc[:, col_].values*ppbv_scale
        if 'ppt' in col_:
            df.loc[:, col_] = df.loc[:, col_].values*pptv_scale

    # - OH concentrations if in NetCDF output
    if IncConcAfterChemDiags:
        try:
            pass
        except KeyError:
            pass


    # - Processing and save?
    # Calculate % change from base case for each variable
    if not isinstance(REF1, type(None)):
        for col_ in df.columns:
            pcent_var = col_+' (% vs. {})'.format(REF1)
            df[pcent_var] = (df[col_]-df[col_][REF1]) / df[col_][REF1] * 100
    if not isinstance(REF2, type(None)):
        for col_ in df.columns:
            pcent_var = col_+' (% vs. {})'.format(REF2)
            df[pcent_var] = (df[col_]-df[col_][REF2]) / df[col_][REF2] * 100

    # Transpose back to variables as index
    df = df.T
    # Re-order columns
    df = df.reindex(sorted(df.columns), axis=1)
    # Reorder index
    df = df.T.reindex(sorted(df.T.columns), axis=1).T
    # Now round the numbers
    df = df.round(round)
    # Save csv to disk
    if save2csv:
        csv_filename = '{}_summary_statistics{}.csv'.format(prefix, extra_str)
        df.to_csv(csv_filename)
    # Return the DataFrame too
    return df


if __name__ == "__main__":
    main()

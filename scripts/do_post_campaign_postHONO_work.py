#!/usr/bin/python
"""
Driver for analysis of "4 pptv HONO world" following the ARNA campaign
"""
import arna as ar
import sys

from AC_tools import GetSpeciesConcDataset, AddChemicalFamily2Dataset, species_mass, get_ref_spec, get_Gg_trop_burden, get_StateMet_ds, tra_unit, get_avg_2D_conc_of_X_weighted_by_Y, GC_var, get_ProdLoss_ds, add_molec_den2ds, rm_fractional_troposphere, constants, GetConcAfterChemDataset

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
#    GC_version = 'v13.4'
    GC_version = 'v12.9'
    RunDict = get_RunDict_of_HONO_runs(RunSet=RunSet, GC_version=GC_version)
    # Temporally use one run for the references values (e.g. statemet)
    use_REF_wd4Met = True
    REF_wd = RunDict['min4pptvHONO']

    # Set dates to use (as model run progresses)
    dates2use = [
    # First month of spin up
    datetime.datetime(2018, 1, 1 ), # Initial month of spin up
    #
    datetime.datetime(2018, 2, 1 ),
#    datetime.datetime(2018, 3, 1 ),
    # first month of spun up ouput
#    datetime.datetime(2018, 1, 1 ),
#    datetime.datetime(2018, 2, 1 ), # Initial month of spin up
    ]
    # Use spin up year?
#    dates2use = [datetime.datetime(2018, 1+i, 1) for i in range(12)]

    # -- Stats on model runs
    # Get generic stats
    extra_burden_specs = ['NOx', 'NIT-all', 'HNO2',  'NOy', 'HNO3',  ]
    extra_surface_specs = extra_burden_specs
    # Get extra stats relate to HONO
    DiagVars = [
    'CO-lifetime', 'OH-lifetime', 'NOx-lifetime', 'NO2-lifetime',
    'CH4-lifetime', 'O3-lifetime',
#    'Ox-lifetime',
    'HNO2-lifetime',
    'HNO2:NOx', 'HNO2:HNO3',
    'HNO3:NOx', 'HNO3:NO2', 'HNO3:NIT',
    'NIT:NOx', 'NO2:NO',
#    'OH-trop-avg',
    'Cl-trop-avg', 'HNO2-trop-avg', 'NOx-trop-avg'
    ]

    # Quick test for model output
    for key in RunDict.keys():
        wd = RunDict[key]
        print(key, RunDict[key])
        files2use = glob.glob( '{}*SpeciesConc.*2018*'.format(wd) )
        print( len(files2use) )
        print(files2use)
        print()

#    df = AC.get_general_stats4run_dict_as_df(run_dict=RunDict,
    use_time_in_trop = True
    df = get_stats4RunDict_as_df(RunDict=RunDict,
                                 dates2use=dates2use,
                                 use_REF_wd4Met=use_REF_wd4Met,
                                 REF_wd=REF_wd,
                                 extra_burden_specs=extra_burden_specs,
                                 extra_surface_specs=extra_surface_specs,
                                 DiagVars=DiagVars,
                                 use_time_in_trop=use_time_in_trop,
                                 verbose=verbose,
                                 debug=debug,
                                 )

    # Change Row ordering
    df = df.T
    FirstDiags = [
    'CH4-lifetime (years)', 'O3 burden (Tg)', 'O3 surface (ppbv)',
    'NOy burden (Tg)',
    'HNO2 surface (ppbv)', 'HNO2:NOx', 'NOy surface (ppbv)',
    'OH surface (molec cm-3)', 'NIT-all burden (Tg)'
    ]
    OtherDiags = [i for i in df.columns if (i not in FirstDiags)]
    df = df[FirstDiags+OtherDiags]
    df = df.T

    # Save
    column_ordering = [
    'Base',
    '4pptHONO','min4pptvHONO',
    'NIThv',
    'NIThvCap',
    'N2O5',
    'OH+NO2',
    'HalNitratesx10' ,
    ]
    SaveName = 'ARNA_Key_Stats_on_model_runs_post_HONO'
    SaveNameAll = '{}_{}_{}'.format(SaveName, RunSet, GC_version)
    SaveNameAll = AC.rm_spaces_and_chars_from_str(SaveNameAll)
    df[column_ordering].to_csv('{}{}'.format(SaveNameAll, '.csv'))

    # Save out ratio values
    REF = 'Base'
    cols2use = [i for i in df.columns if (i != REF)]
    df_diff = df[cols2use]
    for __col in cols2use:
        df_diff.loc[:, __col ] = df_diff.loc[:, __col ] / df[REF]
    df_diff.to_csv('{}{}'.format(SaveNameAll, '_ratio.csv'))

    # Save out core stats for manuscript table
    corestats = []
    df[corestats].to_csv(SaveName)



def plt_spatial_changes_in_4pptv_HONO_world(pcent=True,
                                            REF1='Base',
                                            DIFF='4pptHONO' ):
    """
    Plot up changes in key metrics spatially and zonally
    """
    # Settings
    pcent = True
    REF1 = 'Base'
    DIFF = '4pptHONO'
    # Model runs to use?
    RunSet = 'PostHONO'
    GC_version = 'v12.9'
    RunDict = get_RunDict_of_HONO_runs(RunSet=RunSet, GC_version=GC_version)
    RunDict = {REF1: RunDict[REF1], DIFF: RunDict[DIFF],}
    # Set dates to use (as model run progresses)
#    dates2use = [
#    datetime.datetime(2018, 1, 1 ),
#    datetime.datetime(2018, 2, 1 ),
#    ]
#    dates2use = [datetime.datetime(2019, 1+i, 1) for i in range(12)]
    dates2use = [datetime.datetime(2018, 1+i, 1) for i in range(4)]

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
    savetitle = 'ARNA_spatial_HONO4pptv'
    specs2plot = [ 'O3', 'HNO2',  'HNO3', 'CO' ] + families2use
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
        RunStr = 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.DustUptake'
        RunStr += '.JNIT.Isotherm.BCs.repeat.ON.II.diags'
        RunStr2 = RunStr + '.v2.J00.HourlyOutput'
        RunDict = {
        'Base': '{}{}{}/'.format(RunRoot, RunStr,'.v2.J00.HourlyOutput.2018'),
        'min4pptvHONO': '{}{}{}/'.format(RunRoot, RunStr2,
                                         '.HONO4pptv' ),
        '4pptHONO': '{}{}{}/'.format(RunRoot, RunStr2,
                                     '.HONO4pptv.all' ),
        'HalNitratesx10': '{}{}{}/'.format(RunRoot, RunStr2,
                                           '.2018.HalNitrate'),
        'NIThvCap': '{}{}{}/'.format(RunRoot, RunStr,
                                  '.v3.0.H2O.cap2J50.HONO100.2018'),
        # The below runs do not yet have enough output to analysis
        # But should do before 9am 17th March
        'NIThv': '{}{}{}/'.format(RunRoot, RunStr,
                                  '.v3.0.H2O.AcidII.100HONO.2018/'),
        'OH+NO2': '{}{}{}/'.format(RunRoot, RunStr2, '.2018.NO2andOH'),
        # The below probably will not have enough useful data by then
#        'Amedro2020': '{}{}{}/'.format(RunRoot, RunStr2,
#                                      '.2018.NO2andOH.Amedro'),
#        'NOxSink': '{}{}{}/'.format(RunRoot, RunStr, ''),
#         'HO2+NO':  '{}{}{}/'.format(RunRoot, RunStr2,
#                                      '.v2.J00.HourlyOutput.2018.HO2NO'),
         'N2O5': '{}{}{}/'.format(RunRoot, RunStr2, '.2018.N2O5'),
        }
    else:
        PrtStr = "Unkonwn RunSet ('{}') and GC verions ('{}')"
        print(PrtStr.format(RunSet, GC_version))
        sys.exit(0)
    # Include the output directory folder in directory strings
    for key in RunDict.keys():
        RunDict[key] = '{}{}/'.format(RunDict[key], 'OutputDir')
    return RunDict


def get_stats4RunDict_as_df(RunDict=None,
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
    core_specs = ['O3', 'CO', 'NO', 'NO2']
    HOx_vars = ['HO2', 'OH', 'HOx']
    SCprefix = 'SpeciesConc_'
    CACsuffix = 'concAfterChem'
    prefix = SCprefix
    ALLSp = core_specs + extra_burden_specs + extra_surface_specs
    # Also add all families and specs in 'DiagVars' to this list
    if len(DiagVars) >= 1:
        RatioVars = [i.split(':') for i in DiagVars if ':' in i]
        RatioVars = [item for sublist in RatioVars for item in sublist]
        LifetimeVars = [i for i in DiagVars if '-lifetime' in i]
        LifetimeVars = [i.split('-lifetime')[0] for i in LifetimeVars]
#        LifetimeVars = [item for sublist in LifetimeVars for item in sublist]
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
        StateMet = add_molec_den2ds(StateMet)
        dsS[RunDict_r[REF_wd]] = StateMet
    else:
        for key in RunDict.keys():
            StateMet = get_StateMet_ds(wd=RunDict[key], dates2use=dates2use)
            StateMet = add_molec_den2ds(StateMet)
            dsS[key] = StateMet
    # Create a tropospheric mask
#    trop_mask = create4Dmask4trop_level(StateMet=StateMet)

    # - Get burdens for core species
    avg_over_time = True  # Note: burdens area averaged overtime
    specs2use = list(set(core_specs+extra_burden_specs+ALLSp))
    specs2use = [i for i in specs2use if (i not in HOx_vars)]
    vars2use = [prefix+i for i in specs2use]
    BurdenStr = '{} burden ({})'
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
        varnames = [BurdenStr.format(i, mass_unit) for i in specs2use]
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
    PLprefix = 'Loss'
    vars2use = ['Loss_Ox', 'Prod_Ox']
    if IncProdLossDiags:
        for key in RunDict.keys():
            try:
                # Get StateMet object
                if use_REF_wd4Met:
                    StateMet = dsS[RunDict_r[REF_wd]]
                else:
                    StateMet = dsS[key]
                # Retrieve Prod-loss diagnostics
                ds = get_ProdLoss_ds(wd=RunDict[key], dates2use=dates2use)
                # Only consider troposphere?
                if rm_strat:
                    # Loop by spec
                    if use_time_in_trop:
                        ds = rm_fractional_troposphere(ds,
                                                        vars2use=vars2use,
                                                        StateMet=StateMet)
                    else:
                        ds = ds[PLvar].where(trop_mask)
                for var in vars2use:
                    loss = ds[[var]]
                    # Convert units to be species appropriate
                    # molecules/cm3/s-1 =>  mol/s-1
                    loss *= StateMet['Met_AIRVOL'] * 1E6 / constants('AVG')
                    # mol/s-1 => g/s => Tg/s
                    loss *= species_mass('O3') / 1E12
                    # Convert to Tg/year
                    loss = loss*60*60*24*365
                    units = 'Tg/year'
                    loss = loss.mean(dim='time').sum()[var].values
                    # Save to stats dataframe
                    SaveVar = '{} ({})'.format(var, units)
                    df.loc[SaveVar, key] = loss
                # Include net chemical production
                SaveVar = '{} ({})'.format('POx-LOx', units)
                LVar = '{} ({})'.format('Prod_Ox', units)
                PVar = '{} ({})'.format('Loss_Ox', units)
                df.loc[SaveVar, key] = df.loc[PVar, key] - df.loc[LVar, key]
            except KeyError:
                print('Key error whilst retrieving P/L diags for Ox')

    # - Add lifetime calculations for species
    PtrStr1 = "Calculating lifetime for diag ('{}') for '{}' variable"
    ErrStr1 = "Loss diagnostic not found ({}), skipped lifetime calc ('{}')"
    ErrStr2 = "Prod/Loss files not found ('{}'), skipped lifetime calc ('{}')"
    lifetimes2calc = [i for i in DiagVars if ('lifetime' in i.lower())]
    if (len(lifetimes2calc) >= 1):
        for key in RunDict.keys():
            # Get StateMet object
            if use_REF_wd4Met:
                StateMet = dsS[RunDict_r[REF_wd]]
            else:
                StateMet = dsS[key]
            # Loop and calculate lifetime one species at a time
            for var in lifetimes2calc:
                species2calc = var.split('-lifetime')[0]
                print(var, species2calc)
                if verbose:
                    print(PtrStr1.format(var, species2calc))

                # Find loss value from prod/loss diagnostic
                PLvar = '{}_{}'.format(PLprefix, species2calc)
                try:
                    ds = get_ProdLoss_ds(wd=RunDict[key], dates2use=dates2use)
                    ds = ds[[PLvar]]
                    # Only consider troposphere?
                    if rm_strat:
                        # Loop by spec
                        if use_time_in_trop:
                            ds = rm_fractional_troposphere(ds,
                                                            vars2use=[PLvar],
                                                            StateMet=StateMet)
                        else:
                            ds = ds[PLvar].where(trop_mask)
                    loss = ds[[PLvar]]
                    # Use burden calculated already
                    if debug:
                        print('WARNING: Check units for tropospheric burden')
                    try:
                        BurdenVar = BurdenStr.format(species2calc, 'Tg')
                        burden = df.loc[ BurdenVar, key]
                    except KeyError:
                        ErrStr = 'WARNING: variable not found in df ({})'
                        print(ErrStr.format(BurdenVar))

                    # Convert units to be species appropriate
                    # molecules/cm3/s-1 =>  mol/s-1
                    loss *= StateMet['Met_AIRVOL'] * 1E6 / constants('AVG')
                    # mol/s-1 => g/s => Tg/s
                    loss *= species_mass(species2calc) / 1E12

                    # (e.g. years for CH4, NO2 in minutes ....)
                    lifetime = burden / np.nansum( loss[PLvar].values )
                    LifeimeInDays = ['CO', 'Ox', 'NOx', 'NO', 'NO2']
                    if (species2calc in LifeimeInDays):
                        lifetime = lifetime /60/60/24
                        units = 'days'
                    elif species2calc == 'CH4':
                        lifetime = lifetime /60/60/24/365
                        units = 'years'
                    else:
                        units = 's'
                    # Save to stats dataframe
                    LifetimeVar = '{} ({})'.format(var, units)
                    df.loc[LifetimeVar, key] = lifetime

                except KeyError:
                    # If variable not found, skip lifetime calculation
                    print(ErrStr1.format(PLvar, var))

                except AssertionError:
                    # If files not found, skip lifetime calculation
                    print(ErrStr2.format(key, var))

    # - Add Ratio calculations
    PtrStr = "Calculating ratio for diag ('{}') for '{}' vs '{}'"
    long_nameStr = "Dry mixing ratio of species '{}' ('{}':'{}')"
    ratios2calc = [i for i in DiagVars if (':' in i)]
    prefix = SCprefix
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

                # Weight values and save to dataframe
                # Weight over lon?
#                dstemp.sum(dim=['lon']) / StateMet[MolecVar].sum(dim=['lon'])
                # Weight by molecules
                avg = avg.sum() / StateMet[MolecVar].sum()
                if debug:
                    print(avg)
                    print(avg[var2calc])
                    print( avg[var2calc].values )
                avg = avg[var2calc].values
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
    prefix = SCprefix
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

    # - OH concentrations if in NetCDF output
    if IncConcAfterChemDiags:
        # Hardcore stast on HOx
        vars2use = HOx_vars
        AttStr = '{} concentration immediately after chemistry'
        ErrStr = "Failed to include '{}' diagnostics in df output ('{}')"
        dsCC = {}
        for key in RunDict.keys():
            try:
                ds = GetConcAfterChemDataset(wd=RunDict[key],
                                             dates2use=dates2use)
                # Add family value of HOx into  dataset
                NewVar = '{}{}'.format('HOx', suffix)
                OHvar = '{}{}'.format('OH', suffix)
                HO2var = '{}{}'.format('HO2', suffix)
                ds[NewVar] = ds[OHvar].copy()
                ds[NewVar] = ds[OHvar] +  ds[HO2var]
                attrs = ds[OHvar].attrs
                attrs['long_name'] = AttStr.format('HOx')
                units = attrs['units']
                ds[NewVar].attrs = attrs
                # rename to drop suffix
                OldVars = [i for i in  ds.data_vars if suffix in i]
                NewVars = [i.split(suffix)[0] for i in OldVars]
                ds = ds.rename(name_dict=dict(zip(OldVars, NewVars)))
                dsCC[key] = ds

#            except KeyError:
            except:
                print(ErrStr.format(suffix, key))

        # Include surface weighted values in core df
        for key in RunDict.keys():
            ds = dsCC[key].copy()
            #
            ds = ds.isel(lev=(ds.lev==ds.lev[0])).mean(dim='time')
            for var in vars2use:
                varname = '{} surface ({})'.format(var, units)
                # Save values on a per species basis to series
                val = get_avg_2D_conc_of_X_weighted_by_Y(ds, Xvar=var,
                                                         Yvar='AREA')
                # Save calculated values to dataframe
                df.loc[varname, key] = val

    # - Tropospherically weighted averages for species (e.g. oxidants)
    PtrStr = "Calculating average trop conc ('{}') for '{}' "
    concs2calc = [i for i in DiagVars if ('-trop-avg' in i)]
    if (len(ratios2calc) >= 1):
        for key in RunDict.keys():
            ds = dsD[key]
            # Get StateMet object
            if use_REF_wd4Met:
                StateMet = dsS[RunDict_r[REF_wd]]
            else:
                StateMet = dsS[key]
            # Loop by variable to calculate
            for var2calc in concs2calc:
                species2calc = var2calc.split('-trop-avg')[0]
                # Special case for HOx (HOx, HO2, OH)
                if (species2calc in HOx_vars):
                    prefix = CACsuffix
                    Pstr = "WARNING: skipping calc for '{}' for model run({})"
                    print(Pstr.format(species2calc, key))
                else:
                    prefix = SCprefix
                    print(species2calc)
#                    try
                    # Calculate molecular weighted values
                    dsVar = '{}{}'.format(prefix, species2calc)
                    avg = ds[[dsVar]].copy() * StateMet[MolecVar]
                    # Only consider troposphere?
                    if rm_strat:
                        # Loop by spec
                        if use_time_in_trop:
                            avg = rm_fractional_troposphere(avg,
                                                           vars2use=[dsVar],
                                                            StateMet=StateMet)
                        else:
                            avg = avg[dsVar].where(trop_mask)
                # Weight values and save to dataframe
                # Weight over lon?
#                dstemp.sum(dim=['lon']) / StateMet[MolecVar].sum(dim=['lon'])
                # Weight by molecules
                avg = avg.sum() / StateMet[MolecVar].sum()
                # What scaling / units to use?
                units = tra_unit(spec)
                SaveVar = '{} ({})'.format(var2calc, units)
                if debug:
                    print(avg)
                    print(avg[dsVar])
                    print( avg[dsVar].values )
                avg = avg[dsVar].values
                if debug:
                    print(avg)
                df.loc[var2calc, key] = avg

    # Transpose dataframe
    df = df.T
    # Scale units
    for col_ in df.columns:
        if 'ppb' in col_:
            df.loc[:, col_] = df.loc[:, col_].values*ppbv_scale
        if 'ppt' in col_:
            df.loc[:, col_] = df.loc[:, col_].values*pptv_scale

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


def mk_figure2_HONO_surface_conc():
    """
    Make figure 2 for manuscript on post HONO world
    """
    # Which model output data to plot?
    RunRoot = ar.get_local_folder('RunRoot')
    RunStr = 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.DustUptake'
    RunStr += '.JNIT.Isotherm.BCs.repeat.ON.II.diags.v2.J00.HourlyOutput'
    folder = '{}{}/{}/'.format(RunRoot, RunStr, 'OutputDir')
    # Which dates to use?
    dates2use = [datetime.datetime(2019, 1+i, 1) for i in range(12)]
    # Get data and a surface average
    ds = AC.GetSpeciesConcDataset(wd=folder, dates2use=dates2use)
    ds = ds.mean(dim='time')
    ds = ds.isel(lev=ds.lev == ds.lev[0])
    # plot and save
    savename = 'ARNA_figure02_annual_average_HNO2'
    var2plot='SpeciesConc_HNO2'
    scale = 1E12
    units = 'pptv'
    kwargs={'log':True}
    ds = ds[[var2plot]].sum(dim='lev') * scale
    # TODO: make plot logarithmic, clear colour bar axis

    AC.quick_map_plot(ds, var2plot='SpeciesConc_HNO2', show_plot=False,
                      savename=savename,
#                      save_plot=False,
                      verbose=verbose, kwargs=kwargs)
#    AC.save_plot





if __name__ == "__main__":
    main()

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
        datetime.datetime(2018, 1, 1),  # Initial month of spin up
        #
        datetime.datetime(2018, 2, 1),
        #    datetime.datetime(2018, 3, 1 ),
        # first month of spun up ouput
        #    datetime.datetime(2018, 1, 1 ),
        #    datetime.datetime(2018, 2, 1 ), # Initial month of spin up
    ]
    # Use spin up year?
#    dates2use = [datetime.datetime(2018, 1+i, 1) for i in range(12)]

    # -- Stats on model runs
    # Get generic stats
    extra_burden_specs = ['NOx', 'NIT-all', 'HNO2',  'NOy', 'HNO3', ]
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
        files2use = glob.glob('{}*SpeciesConc.*2018*'.format(wd))
        print(len(files2use))
        if verbose:
            print(files2use)
        print()

#    df = AC.get_general_stats4run_dict_as_df(run_dict=RunDict,
    use_time_in_trop = True
    df = AC.get_stats4RunDict_as_df(RunDict=RunDict,
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
        '4pptHONO', 'min4pptvHONO',
        'NIThv',
        'NIThvCap',
        'N2O5',
        'OH+NO2',
        'HalNitratesx10',
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
        df_diff.loc[:, __col] = df_diff.loc[:, __col] / df[REF]
    df_diff.to_csv('{}{}'.format(SaveNameAll, '_ratio.csv'))

    # Save out core stats for manuscript table
    corestats = []
    df[corestats].to_csv(SaveName)


def plt_spatial_changes_in_4pptv_HONO_world(pcent=True,
                                            REF1='Base',
                                            DIFF='4pptHONO'):
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
    RunDict = {REF1: RunDict[REF1], DIFF: RunDict[DIFF], }
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
    OldNames = [i for i in ds.data_vars if (prefix in i)]
    NewNames = [i.split(prefix)[-1] for i in OldNames]
    for key in dsD.keys():
        ds = dsD[key]
        ds = ds.rename(name_dict=dict(zip(OldNames, NewNames)))
        dsD[key] = ds

    # Get MetState object

    # -- Graphics on model runs
    savetitle = 'ARNA_spatial_HONO4pptv'
    specs2plot = ['O3', 'HNO2',  'HNO3', 'CO'] + families2use
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
            ds2plot = (ds2 - ds1) / ds1 * 100
        else:
            ds2plot = ds2 - ds1

        # Select the surface
        ds2plot = ds2plot.sel(lev=ds2plot.lev.values[0])

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
            'Base': '{}{}{}/'.format(RunRoot, RunStr, '.orig.1monthTest'),
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
            'Base': '{}{}{}/'.format(RunRoot, RunStr, '.v2.J00.HourlyOutput.2018'),
            'min4pptvHONO': '{}{}{}/'.format(RunRoot, RunStr2,
                                             '.HONO4pptv'),
            '4pptHONO': '{}{}{}/'.format(RunRoot, RunStr2,
                                         '.HONO4pptv.all'),
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
    var2plot = 'SpeciesConc_HNO2'
    scale = 1E12
    units = 'pptv'
    kwargs = {'log': True}
    ds = ds[[var2plot]].sum(dim='lev') * scale
    # TODO: make plot logarithmic, clear colour bar axis

    AC.quick_map_plot(ds, var2plot='SpeciesConc_HNO2', show_plot=False,
                      savename=savename,
                      #                      save_plot=False,
                      verbose=verbose, kwargs=kwargs)
#    AC.save_plot


if __name__ == "__main__":
    main()

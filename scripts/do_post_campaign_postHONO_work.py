#!/usr/bin/python
"""
Driver for analysis of "4 pptv HONO world" following the ARNA campaign
"""
import arna as ar


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
    RunDict = get_RunDict_of_HONO_runs()
    # Temporally use one run for the references values (e.g. statemet)
    use_REF_wd4Met = True
    REF_wd = RunDict['4pptHONO']

    # Set dates to use (as model run progresses)
    dates2use = [
    datetime.datetime(2018, 1, 1 ),
    datetime.datetime(2018, 2, 1 ),
    ]

    # -- Stats on model runs
    # Get generic stats
    extra_burden_specs = ['NOx', 'NIT-all', 'HNO2',  'NOy', 'HNO3',  ]
    extra_surface_specs = extra_burden_specs
    df = AC.get_general_stats4run_dict_as_df(run_dict=RunDict,
                                             dates2use=dates2use,
                                             use_REF_wd4Met=use_REF_wd4Met,
                                             REF_wd=REF_wd,
                                         extra_burden_specs=extra_burden_specs,
                                     extra_surface_specs=extra_surface_specs,
                                             )

def plt_spatial_changes_in_4pptv_HONO_world()
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


def get_RunDict_of_HONO_runs():
    """
    Retrieve the model runs
    """
    # Which runs to use
    RunRoot = ar.get_local_folder('RunRoot')
    RunStr = 'gc_4x5_47L_geosfp_fullchem.v13.4.0-rc.2'
    RunDict = {
    'Base': '{}{}{}/'.format(RunRoot,RunStr,'.orig.1monthTest'),
    '4pptHONO': '{}{}{}/'.format(RunRoot,RunStr,'.orig.ARNA.HONO.4pptv'),
    }
    # Include the output directory folder in directory strings
    for key in RunDict.keys():
        RunDict[key] = '{}{}/'.format(RunDict[key], 'OutputDir')
    return RunDict


if __name__ == "__main__":
    main()

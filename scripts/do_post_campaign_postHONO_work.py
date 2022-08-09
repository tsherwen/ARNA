#!/usr/bin/python
"""
Driver for analysis of "4 pptv HONO world" following the ARNA campaign
"""
import arna as ar
import sys

from matplotlib.colors import LogNorm


from AC_tools import GetSpeciesConcDataset, AddChemicalFamily2Dataset, species_mass, get_ref_spec, get_Gg_trop_burden, get_StateMet_ds, tra_unit, get_avg_2D_conc_of_X_weighted_by_Y, GC_var, get_ProdLoss_ds, add_molec_den2ds, rm_fractional_troposphere, constants, GetConcAfterChemDataset, get_StateMet_ds, get_DryDep_ds, get_WetLossConv_ds, get_WetLossLS_ds

from arna import get_local_folder, get_tags_for_NOx_HONO, get_DryDepAndWetDep_ds


def main():
    """
    Main driver for 4 pptv HONO world analysis
    """
    # --- Local settings to pass to all analysis/plotting functions
    RunSet = 'PostHONO'
    res = '4x5'
    GC_version = 'v12.9'
#    GC_version = 'v13.4'
#    sdate = datetime.datetime(2018, 1, 1) # Beginning for spin-up year
    sdate = datetime.datetime(2019, 1, 1) # Beginning for analysis year
#    edate = datetime.datetime(2018, 3, 1) # 3 months into spin-up
#    edate = datetime.datetime(2018, 6, 1) # 6 months into spin-up
#    edate = datetime.datetime(2018, 10, 1) # 10 months into spin-up
    edate = datetime.datetime(2019, 12, 31) # End of analysis year
#    edate = datetime.datetime(2018, 12, 31) # End of spin-up year
    dates2use = pd.date_range(sdate, edate, freq='1D')
#    dates2use = None

    # Base for analysis and then perturbation
    REF1 = 'Base'
#    DIFF = '4pptHONO'
#    DIFF = 'min4pptHONO'
#    DIFF = 'Iso.Delq.NoJCap'
    DIFF = 'Iso.UnlimitedAll'

    # Get dictionary of model runs and their names (RunDict)
    RunDict = ar.get_dict_of_GEOSChem_model_output(RunSet=RunSet,
                                                   GC_version=GC_version,
                                                   res=res,
                                                   folder4netCDF=True)

    # TEMP: mannual set to use v6 diags and HO2+NOv2 here
    RunRoot = ar.get_local_folder('RunRoot')
    RunBaseStr1 = 'geosfp_4x5_aciduptake.v12.9.0'
    RunStr1 = RunBaseStr1 + '.ARNA.Isotherm.Diags.v6.HO2nNOv2.1'
    RunStr2 = '{}{}{}/OutputDir/'.format(RunRoot, RunStr1, '.IsopH4')
    RunStr1 = '{}{}{}/OutputDir/'.format(RunRoot, RunStr1, '')
    RunStr3 = RunBaseStr1 + '.BASE.2019.2020.ARNA.DustUptake.'
    RunStr3 += 'JNIT.Isotherm.BCs.repeat.ON.II.diags.v5.0.H2O.AcidII.100HONO'
    RunStr3 = '{}{}{}/OutputDir/'.format(RunRoot, RunStr3, '.2018.pH7')
    RunStr4 = RunBaseStr1 + '.ARNA.Isotherm.Diags.v6'
    BaseRunStr5 =  RunStr4 +  '.HO2nNOv2.5'
    RunStr5 = '{}{}{}/OutputDir/'.format(RunRoot, BaseRunStr5, '')
    RunStr6 = '{}{}{}/OutputDir/'.format(RunRoot, BaseRunStr5, '.IsopH4')
    RunStr7 = '{}{}{}/OutputDir/'.format(RunRoot, BaseRunStr5,
                                          '.Iso.UnlimitedpH')
    RunStr8 = '{}{}{}/OutputDir/'.format(RunRoot, RunStr4, '.Ye2017')
    RunStr9 = '{}{}{}/OutputDir/'.format(RunRoot, RunStr4, '.Ye2017online')
    RunStr10 = '{}{}{}/OutputDir/'.format(RunRoot, RunStr4, '.Ye2017.LL')
    RunStr11 = '{}{}{}/OutputDir/'.format(RunRoot, RunStr4,'.Iso.UnlimitedpH')
    RunStr12 = '{}{}{}/OutputDir/'.format(RunRoot, RunStr4,'.Iso.UnlimitedAll')
    RunStr13 = '{}{}{}/OutputDir/'.format(RunRoot, RunStr4,
                                          '.Iso.UnlimitedAll.repeat')

    RunDict2 = {
    'Base':RunDict['Base'], 'min4pptHONO': RunDict['min4pptHONO'],
#    'HO2+NOv2.1': RunStr1,
    'HO2+NO': RunStr5, # HO2+NO v5
    'Iso.Delq.CapJ50.pH4.HO2+NO': RunStr6, # HO2+NO v5
#    'Isov5.pH7': RunStr3, # Note, this is not the latest isotherm code.
#    'Iso.HO2+NO': RunStr2, # WARNING: This is using HO2+NO code with a bug
    'Iso.Delq.NoJCap.HO2+NO': RunStr7,
#    'Iso.Delq.NoJCap': RunStr11, # This has missing months in analysis year
    'Ye17': RunStr8,
    'Ye17.online': RunStr9,
    'Ye17.LL': RunStr10,
#    'Iso.UnlimitedAll.Base600': RunStr12,
    'Iso.UnlimitedAll': RunStr13, # This has the
    }

    RunDict = RunDict2

    # Kludge: temp
    keys2del = [
    'Ye17.online', 'Ye17', 'Ye17.LL',
    'HO2+NOv2.1',
#    'min4pptHONO',
    'HO2+NOv2.5',
    ]
    for key in keys2del:
        try:
            del RunDict[key]
        except KeyError:
            print("WARNING: key not dropped from RunDict ('{}')".format(key))

#
# 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags.v2.J00.HourlyOutput.2018.HO2andNOv2.BetaStandard.1day'
#
# 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags.v2.J00.HourlyOutput.2018.HO2andNOv2.BetaZero.1day'

    # Do checks on which runs have output?
#    df = GetModelOutputFileDates(wd=wd)

    # --- Calls to plotting and analysis functions
    # Do the core stats for the runs
    df = do_analysis_of_4pptv_HONO_world(dates2use=dates2use,
                                         REF1=REF1, DIFF=DIFF,
                                         RunDict=RunDict)

    # plot up key changes
    plt_spatial_changes_in_4pptv_HONO_world(pcent=True, REF1=REF1, DIFF=DIFF,
                                            dates2use=dates2use,
                                            RunDict=RunDict)
    plt_spatial_changes_in_4pptv_HONO_world(pcent=False, REF1=REF1, DIFF=DIFF,
                                            dates2use=dates2use,
                                            RunDict=RunDict)
    # Plot up all variables changes (in percent)
    plt_spatial_changes_in_4pptv_HONO_world(pcent=True, dates2use=dates2use,
                                            PltAllVars=True,
                                            RunDict=RunDict)

    # Plot up the JScale changes
    plt_JScale_analysis(RunDict=RunDict, dates2use=dates2use)

    # Plot up HONO production routes - N/A
    # See functions in do_post_campaign_NOx_analysis.py
    # (e.g.'plt_key_NOx_budget_terms' )
    plt_HO2_NO_branching_ratio_spatially()

    # Plot up the NO + HO2 + H2O => HNO3 channel
    plt_NO_HO2_channel_verues_H2O()
    # And the OH + NO2 channel...
    plt_OH_NO2_channel_verues_H2O()


    # Do a full analysis of variables for JNIT run
    DIFF = None
    for REF1 in RunDict.keys():
        variable_analyis4JNIT(REF1=REF1, DIFF=DIFF, NOxD=NOxD,
                              RunDict=RunDict, dates2use=dates2use)

    # plot the masks being used here spatially
    CheckMasksSpatially()


def AddMask2Dataset(ds, MaskName=None, res=None):
    """
    Apply a mask to a xr.dataset
    """
    # Retrieve Mask for MaskName.
    # NOTE: This should be applicable to any xr.dataset, but right now only
    #       GEOS-Chem standard arrays will fit.
    mask = RetrieveMaskFromName(MaskName=MaskName, res=res)

    # Ensure mask is in format to be able to apply to data
    # Ensure the mask can be used for masking of Datasets
    print('TODO: convert mask to xr.dataset? With lat and lon')
    print('WARNING: Kludge to consider first dimensions applied as a stopgap')
    mask = mask[...,0]

    # Manually construct NetCDF?
    # Use existing Dataset as the template
    var2copy = list(ds.data_vars)[0]
    ds[MaskName] = ds[var2copy].copy()
    attrs = ds[var2copy].attrs
    ds[MaskName] = (('lon', 'lat'), mask)
    attrs['long_name'] = "Mask of '{}'".format(MaskName)
    ds[MaskName].attrs = attrs

    return ds


def CheckMasksSpatially(ds=None, ):
    """
    Plot up masks to check the spatial extents these
    """
    if isinstance(ds, type(None)):
        print('WARNING')
        sys.exit(0)

    # which masks to test?
    MaskNames = [
            'global',
            'local_CVAO_area', 'Cape_Verde_Flying', 'inflow_CVAO_area',
            'CONUS',
            'Ocean Tropics', 'Tropics', 'Ocean', 'Land',
            ]
    MaskNames = [
    'tropics', 'Mid lats', 'south pole', 'north pole', 'Global', 'Ocean',
    'Ocn. Trop.', 'Ex. Tropics', 'NH', 'SH', 'Ice', 'Land', 'lat40_2_40',
    'Land Tropics', 'surface', 'Ocean Sur.', 'Land Sur.', 'Ice Sur.',
    '50S-50N', 'Ocn. 50S-50N',
#    'North Sea',  'Irish Sea', 'Black Sea', 'location',
    'Mediterranean Sea',
    'EU', 'Land Tropics Sur.',
    'Boreal Land',
#    'Alps',
    'France', 'CONUS',
     'Cape_Verde_Flying',
    'local_CVAO_area', 'inflow_CVAO_area',
    ]

    # Get some basic data to use for masking tests
#    get_surface_area()
    if isinstance(wd, type(None)):
        wd = '/users/ts551/scratch/GC/rundirs/P_ARNA/'
        wd +='geosfp_4x5_aciduptake.v12.9.0.ARNA.Isotherm.Diags.v6.HO2nNOv2.5/'
        wd += '/OutputDir/'

    # Use ozone surface data
    ds = AC.GetSpeciesConcDataset(wd=wd, dates2use=dates2use)
    ds = ds.isel(lev=(ds.lev==ds.lev[0]))
    ds = ds.mean(dim='time')
    ds = ds.squeeze()

#    get_surface_area
    # Loop and apply masks
    dsD = {}
    for MaskName in MaskNames:
        print(MaskName)
        #
        ds = AddMask2Dataset(ds.copy(), MaskName=MaskName, res=res)

        # Also mask for bottom 10 km? - N/A for now

        # Save masked array
#        dsD[MaskName] = ds.copy() * ds[MaskName]
        dsD[MaskName] = ds.copy().where(ds[MaskName] != 1.0)

    # Loop and plot and save to a PDF
    savetitle = 'Spatial_plots_of_masks_from_AC_Tools'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    var2plot = 'SpeciesConc_O3'
    kwargs = {}
    for MaskName in MaskNames:

        ds = dsD[MaskName].copy()
        #
        AC.quick_map_plot(ds, var2plot=var2plot, verbose=verbose, **kwargs)
        plt.title(MaskName)

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')
    #



def RetrieveMaskFromName(MaskName=None, res='4x5',
                         mask3D=True, mask2D=False, mask4D=False,
                         trop_limit=False):
    """
    Retrieve mask in format to use with xr.Datasets from provided name
    """

    mask = AC.mask_all_but(region=MaskName, res=res,
                           mask2D=mask2D, mask3D=mask3D, mask4D=mask4D,
                           use_multiply_method=False, trop_limit=False)

    return mask



def variable_analyis4JNIT(dates2use=None, RunDict=None, NOxD=None,
                          REF1=None, DIFF=None,
                          AvgOverTime=True):
    """
    Explore the key variables for nitrate photolysis
    """
    # Local variables
    SCprefix = 'SpeciesConc_'

#    masks = [None, 'Atlantic', 'ARNA', 'FIREX-AQ', ]
#    masks = ['Global', 'Atlantic', 'ARNA', 'FIREX-AQ', ]
#    masks = ['global', 'Atlantic', 'ARNA', 'FIREX-AQ', ]
#    masks = ['global', 'Ocean Tropics', 'Tropics', 'Ocean', 'Land' ]
    masks = ['global',
            'local_CVAO_area', 'Cape_Verde_Flying', 'CONUS',
            'Ocean Tropics', 'Tropics', 'Ocean', 'Land',
            ]
    res = '4x5' # Hardcode resolution for now. This should not be needed

    # Get entire NOx budget as a Dataset object
    # (Note: unit processing etc has should be done already - check this!)
    if isinstance(NOxD, type(None)):
        ExtraConcVars = [
            'NOx', 'NOy', 'NIT-all', 'Cly', 'Bry', 'Iy',
            'DST-all', 'DSTAL-all'
            ]
        ExtraConcVars += [SCprefix+i for i in ['O3'] ]
        NOxD = ar.get_NOx_budget_ds_dict_for_runs(dates2use=dates2use,
                                                  RunDict=RunDict,
                                                  ExtraConcVars=ExtraConcVars,
                                                  debug=debug)
    # Set xr.Dataset to use
    if not isinstance(DIFF, type(None)):
        # Calculate difference
        ds = NOxD[DIFF] - NOxD[REF1]
        AnalysisStr = '_{}_vs_{}_'.format(DIFF, REF1)
        # Ratio?

        # Percentage?

    else:
        # Use the reference run for analysis
        ds = NOxD[REF1].copy()
        AnalysisStr = REF1

    # Just use a time average
    if AvgOverTime:
        ds = ds.mean(dim='time')

    # Limit or average data vertically.
    # hardcode selecting surface for now.
    # Should have the option to molecule weight over alt
    lev2use = ds.lev.values[0]
    ds = ds.sel(lev=lev2use)

    # Apply masks to dataset
    dsD = {}
    for MaskName in masks:
        print(MaskName)
        #
        ds = AddMask2Dataset(ds.copy(), MaskName=MaskName, res=res)

        # Also mask for bottom 10 km? - N/A for now

        # Save masked array
#        dsD[MaskName] = ds.copy() * ds[MaskName]
        dsD[MaskName] = ds.copy().where(ds[MaskName] != 1.0)


    # What variables to get stats on?


    # Get Stats on the different
    df = pd.DataFrame()
    for key in dsD.keys():
        print(key)
        ds = dsD[key]
        SaveName = 'TEST_{}'.format(key)
        SaveName = AC.rm_spaces_and_chars_from_str(SaveName)
        ds.to_netcdf('{}.nc'.format(SaveName))
        #
#        vars2use = list(ds.data_vars)[:10]
        vars2use = list(ds.data_vars)


        # apply stats (e.g max, min, mean, area weighted?)
        for var in vars2use:
            avg = AC.get_avg_2D_conc_of_X_weighted_by_Y(ds, Xvar=var,
                                                        Yvar='AREA')
            # Update naming and units for variables
            if SCprefix in var:
                SpecName = var.split(SCprefix)[-1]
                units, scaleby = AC.tra_unit(VarName, scale=True)
                VarName = '{} ({})'.format(SpecName, units)
                avg = avg * scaleby
            else:
                VarName = var

            df.loc[key, VarName ] = avg

    # Save output the values with all the masking
    SaveName = 'ARNA_data_masked_by_area_{}.csv'.format(AnalysisStr)
    df.to_csv(SaveName)


    # - Plot up values
    # Setup PDF to hold surface plots

    #



def plt_JScale_analysis(RunDict=None, dates2use=None):
    """
    Plot up Jscale analysis vertically and at different horizontal levels
    """
    # Set data to use?
    if isinstance(RunDict, type(None)):
        keys2use = ['Base', 'min4pptHONO', 'HO2+NOv2.5',
                    'Isov5.pH4.HO2+NOv2.5',
#                    'Isov5'
                    ]

    # Include J100 as a reference
    RunStr = '/users/ts551/scratch/GC/rundirs/P_ARNA/'
    RunStr += 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.'
    J100 = 'DustUptake.JNITx100.BCs/OutputDir/'
#    J50 = 'DustUptake.JNITx50.BCs/OutputDir/'
    RunDict['J100'] = J100
#    RunDict['J50'] = J50

    # - Vertical
    plt_vertical_JScale_analysis()


    # -  Surface
    levels2use


    pass


def plt_vertical_JScale_analysis(RunDict=None, dates2use=None):
    """
    Plot vertical comparisons of Jscale
    """
    #
    regions = [
    None, 'local_CVAO_area', 'Cape_Verde_Flying', 'CONUS', 'Atlantic'
    ]
#    regions = [None, 'local_CVAO_area']
    # Local variables
    PressVar = 'Met_PMID'
    JScaleVar = 'Jscale'

    # Get StateMet collection - for unit conversions...
    StateD = {}
    for key in RunDict.keys():
        StateMet = AC.get_StateMet_ds(wd=RunDict[key],
                                         dates2use=dates2use)
        # Calculate number of molecules
        MolecVar = 'Met_MOLCES'
        try:
            StateMet[MolecVar]
        except KeyError:
            StateMet = add_molec_den2ds(StateMet, MolecVar=MolecVar)

        StateD[key] = StateMet

    # Get Photolysis rates
    JValD = {}
    for key in RunDict.keys():
        try:
            ds = AC.GetJValuesDataset(wd=RunDict[key],
                                              dates2use=dates2use)
            # Add JScale variable?

            ds[JScaleVar] = ds['Jval_NIT'].copy()
            ds[JScaleVar] = ds[JScaleVar] / ds['Jval_HNO3']
            JValD[key] = ds

        except AssertionError:
            print('WARNING: No JVal diagnostics found for {}'.format(key))

    # Select the average vertical values
    dfD = {}
    avg_over_time = True
    for key in RunDict.keys():

        # Setup dictionary to store masked data
        PlotD = {}
        for region in regions:

            # Retrieve relevent data
            ds = JValD[key][[JScaleVar]].copy()
            StateMet = StateD[key].copy()
            ds[PressVar] = StateMet[PressVar].copy()

            # Limit to region if requested
            if region != None:
                d = ar.get_analysis_region(region)
                x0, x1, y0, y1 = (d['x0'], d['x1'], d['y0'], d['y1'])
                # Reduce plotting arrays by the values
                bool1 = ((ds.lon >= y0) & (ds.lon <= y1)).values
                bool2 = ((ds.lat >= x0) & (ds.lat <= x1)).values
                # Cut by lon, then lat
                ds = ds.isel(lon=bool1)
                ds = ds.isel(lat=bool2)
#                StateMet =
#                ds =
            else:
                pass

            # Take averages over lat
            ds = ds.mean(dim='lat')
#            StateMet = StateMet.mean(dim='lat')
            # weighted by molecules?
            # Weight by molecules over lon
#            ds = ds * StateMet[MolecVar]
#            ds = ds.sum(dim=['lon']) / StateMet[MolecVar].sum(dim=['lon'])

            ds = ds.mean(dim='lon')
#            StateMet = StateMet.mean(dim='lon')

            # Select a specific month? (e.g. February) or average overtime
            if avg_over_time:
                ds = ds.mean(dim='time')
#                StateMet = StateMet.mean(dim='lon')
            else:
                pass

            # Store plotting data as a DataFrame
            df = ds[[JScaleVar, PressVar]].to_dataframe()
            PlotD[region] = df
        # Store dictionary of dataframes for plotting
        dfD[key] = PlotD

        # Do some memory management...
        gc.collect()

    # Plot up
#    sns.set_context(context)
    savetitle = 'ARNA_vertical_Jscale'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    for region in regions:
        fig, ax = plt.subplots()
        print(region)

        for key in RunDict.keys():

            df = dfD[key][region]
            #
#            StateMet = StateD[key]
            X = df[JScaleVar]
            Y = df[PressVar]

            # Limit plots to 10 km (~270 hPa)
            __bool = Y>260.0
            X = X[__bool]
            Y = Y[__bool]

            # Plot up
            ax.plot(X, Y, label=key)
#            ax = plt.gca()
            ax.set_xlabel('Jscale (xJHNO3)')
            ax.set_ylabel('Pressure (hPa)')
            ax.invert_yaxis()

        # Add a title and legend
        plt.title('{}'.format(region))
        plt.legend()

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()


    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plt_HO2_NO_branching_ratio_spatially(wd=None, RunDict=None,
                                         dates2use=None, StateMet=None):
    """
    Plot up the HO2+NO branch ratio spatially using model output.
    """
    # Which dates and model runs to use?
    if isinstance(wd, type(None)):
        RunRoot = ar.get_local_folder('RunRoot')
        Runstr = 'geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.DustUptake'
        Runstr += '.JNIT.Isotherm.BCs.repeat.ON.II.diags.v2.J00.HourlyOutput'
        Runstr += '.2018.HO2andNOv2.BetaStandard.1day/'
#        Runstr += '.2018.HO2andNOv2.BetaZero.1day/'
        wd = '{}{}/OutputDir/'.format(RunRoot, Runstr)
    if isinstance(RunDict, type(None)):
        RunDict = {'BetaStandard': wd}
    if isinstance(dates2use, type(None)):
        dates2use = None
    # Get StateMet object
    if isinstance(StateMet, type(None)):
        StateMet = get_StateMet_ds(wd=wd)
#        dates2use = [datetime.datetime(2018, 1+i, 1) for i in range(3)]

    # Get entire NOx budget as a Dataset object
    NOxD = ar.get_NOx_budget_ds_dict_for_runs(dates2use=dates2use,
                                              RunDict=RunDict)

    # Add ratio of channels
    ChanA = 'ProdNOnHO2ChannelA'
    ChanB = 'ProdNOnHO2ChannelB'
    ChanRatio = 'ProdNOnHO2ChanB_Pcent'
    vars2plot = [ChanRatio, ChanA, ChanB]
    for key in NOxD.keys():
        ds = NOxD[key]
        attrs = ds[ChanA].attrs
        ds[ChanRatio] = ds[ChanA].copy()
        ds[ChanRatio] = ds[ChanB] / ( ds[ChanA] + ds[ChanB] ) * 100
        attrs['long_name'] = 'Braching ratio NO+HO2 (B/A+B, %)'
        ds[ChanRatio].attrs = attrs
        NOxD[key] = ds

    # - Plot up zonal and surface values and ratio of branches
    # Loop by model output for a given run
    for key in NOxD.keys():

        # Setup plot
        savetitle = 'ARNA_spatial_HO2_NO_Channels_{}'.format(key)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        kwargs = {}

        # Get data
        ds = NOxD[key].copy()
        ds = ds.mean(dim='time')
        ds = ds.sel(lev=ds.lev.values[0])

        # plot surface
        for var2plot in vars2plot:

            #    for var2plot in [ 'Prod_HNO2']:
            AC.quick_map_plot(ds, var2plot=var2plot, verbose=verbose, **kwargs)
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
    #        del ds2plot

            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()

        # Get data
        ds = NOxD[key].copy()
        ds = ds.mean(dim='time')

        # plot zonally
        for var2plot in vars2plot:

            fig, ax = plt.subplots()
            im = AC.ds2zonal_plot(ds, var2plot=var2plot, StateMet=ds,
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

        # Save PDF
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')



def get_GEOSChem_H2O(units='molec/cm3', wd=None, rm_strat=True,
                     rtn_units=False, MolecVar='Met_MOLCES',
                     MetH2OVar='Met_AVGW', NewH2OVar='H2O',
                     StateMet=None):
    """
    Retrieve array of H2O from GEOS-Chem output
    """
    # Get StateMet object
    if isinstance(StateMet, type(None)):
        StateMet = get_StateMet_ds(wd=wd)
    # Get Water vapour in molecules/cm3
    # NOTE: Butkovskaya used: 4 × 10E-17 molec cm-3 (∼50% relative humidity)
    try:
        StateMet[MolecVar]
    except:
        StateMet = add_molec_den2ds(StateMet)
    # Select molecules variables and Water vapor volume mixing ratio
    ds = StateMet[[MolecVar, MetH2OVar]]
    if rm_strat:
        ds = rm_fractional_troposphere(ds, vars2use=[MolecVar],
                                          StateMet=StateMet)
    ds[NewH2OVar] = ds[MetH2OVar].copy()
    attrs = ds[MetH2OVar].attrs
    # Convert units
    if (units == 'molec/cm3'):
        ds[NewH2OVar] = ds[NewH2OVar] * ds[MolecVar]
        long_name = '[Water vapour] calculated from H2O mixing ratio'
        long_name += ' (w/r/t dry air) in {}'
        attrs['long_name'] = long_name.format(units)
        attrs['units'] = units
        ds[MetH2OVar].attrs = attrs
    else:
        units = 'vol H2O/vol dry air'
    # return the xr.dataset object
    if rtn_units:
        return {NewH2OVar:ds, 'units':units}
    else:
        return ds


def get_GEOSChem_StateMet_vars4rxns(wd=None, rm_strat=True,
                                    AvgOverTime=False,
                                    rtn1Darrays=True,
                                    LimitTemp=True,
                                    debug=False):
    """
    Return a dictionary of core variables used for
    """
    hPa2TORR = 0.7500616827
#    try:
    if True:
        if isinstance(wd, type(None)):
            # Get model runs to use for getin tropospheric H2O values
            RunDict = ar.get_dict_of_GEOSChem_model_output(RunSet='PostHONO',
                                                           GC_version='v12.9',
                                                           res='4x5',
                                                           folder4netCDF=True)
            wd = RunDict['min4pptHONO']
        # Get Water vapour in molecules/cm3
        StateMet = get_StateMet_ds(wd=wd)
        ds = AC.get_GEOSChem_H2O(wd=wd, StateMet=StateMet)
        NewVar = 'H2O'
        ds = ds[[NewVar]].copy()
        if AvgOverTime:
           ds = ds.mean(dim='time')
        H2O = ds[NewVar].values.ravel()
        # Get global water, temperature and pressure
        # Mask the troposphere(?)
        if rm_strat:
            StateMet = rm_fractional_troposphere(StateMet, vars2use=None,
                                                 StateMet=StateMet.copy())
        TEMP = StateMet['Met_T'] #+ 273.15
        PRESS = StateMet['Met_PMID']
        PRESS_TORR = StateMet['Met_PMID'] * hPa2TORR

        # Also get key reactants?

        # Average over time(?), then flatten the array
        ArrayList = [TEMP, PRESS, PRESS_TORR]
        for n, ds in enumerate(ArrayList):
            ds = ArrayList[n]
            if AvgOverTime:
                ds = ds.mean(dim='time')
            if rtn1Darrays:
                ds = ds.values.ravel()
            ArrayList[n] = ds
        TEMP, PRESS, PRESS_TORR = ArrayList

        # Ensure all arrays are 1D
#         ArrayList = [H2O, TEMP, PRESS, PRESS_TORR]
#         for n, ds in enumerate(ArrayList):
#             ds = ArrayList[n]
#             if (len(ds.shape) > 1):
#                 ds = ds.flatten()
#                 ArrayList[n] = ds

        # Only consider temperatures >200 Kelvin
        if LimitTemp:
#            if rtn1Darrays:
            TEMP[ TEMP<200 ] = np.NaN
#            else:
            # Only consider array locations with non NaN values
            idx = np.argwhere(~np.isnan(TEMP))
            TEMP = TEMP[idx].ravel()
            H2O = H2O[idx].ravel()
            PRESS_TORR = PRESS_TORR[idx].ravel()
            PRESS = PRESS[idx].ravel()

        # Sort the arrays by temperature
        if rtn1Darrays:
            if debug:
                print([i.shape for i in (H2O, TEMP, PRESS, PRESS_TORR)])
            idx = np.argsort(TEMP, axis=0)
            H2O = H2O[idx]
            TEMP = TEMP[idx]
            PRESS_TORR = PRESS_TORR[idx]
            PRESS = PRESS[idx]
#    except:
#         # Just use some reasonable manual numbers for the troposphere
#         TEMP = np.linspace(200, 298)
#         PRESS = np.linspace(150, 1013., len(TEMP))
#         PRESS_TORR = PRESS * 0.7500616827

    # Return as a dictionary of numpy arrays
    d = {
    'TEMP': TEMP,
    'PRESS': PRESS, 'PRESS_TORR': PRESS_TORR,
    'H2O':H2O,
#    'M':M,
    }

    return d


def plt_NO_HO2_channel_as_Arrhenius_form(ArrheniusPlt=False):
    """
    Plot up NO and HO2 reaction as an Arrhenius plot
    """
    # Get H2O, PRESS, T, etc from GEOS-Chem output
    d = get_GEOSChem_StateMet_vars4rxns(AvgOverTime=True, rtn1Darrays=True)
    PRESS_TORR = d['PRESS_TORR']
    PRESS = d['PRESS']
    H2O = d['H2O']
    TEMP = d['TEMP']
    # General rate for HO2+NO=>NO+OH
    k1 = AC.GCARR(3.30E-12, 0.0E+00, 270.0, TEMP=TEMP)
    # Calc new Alpha and Beta calculation
    RateDict = Calculate_Butkovskaya_channel_spit_P_T(TEMP=TEMP, H2O=H2O,
                                                      PRESS_TORR=PRESS_TORR )

    # Setup plot and colours colours
    savetitle = 'ARNA_NO_HO2_Arrhenius_channels_V'
    context = "paper"
    font_scale = 1
    import seaborn as sns
    sns.color_palette('colorblind')
    sns.set_context(context)
    plt.close('all')
    colors = AC.get_CB_color_cycle()

    # - Plot the Chemical space
    # Arrhenius plot?
    if ArrheniusPlt:
        xlabel = '1 / Temperature (K)'
        savetitle += '_1overK'
        xvalues = 1/TEMP
    else:
        xvalues = TEMP
        xlabel = 'Temperature (K)'
        xmin, xmax = 220, 320
    ylabel = 'Reaction rate molecules cm$^{-3}$ s$^{-1}$'
    # Plot Base rate
    plt.scatter(xvalues, k1, label='k1', color=colors[0])
    # Plot alpha and beta channels
#    vars2plot = RateDict.keys()
#    vars2plot = ['k1a', 'k1b']
#    for n, key in enumerate(vars2plot):
#        plt.scatter(xvalues, RateDict[key], label=key, color=colors[n+1])

    # Plot New rate
#     plt.scatter(xvalues, (RateDict['Beta']*k1), label='k1*beta',
#                 color=colors[1])

    # Plot New rate
#    plt.scatter(xvalues, (RateDict['Beta']/100*k1), label='k1*beta/100',
#                color=colors[1], alpha=0.3)

    # Plot New rate
    plt.scatter(xvalues, (RateDict['k1b.beta']), label='k1*beta/100',
                color=colors[1], alpha=0.3)

    plt.scatter(xvalues, (RateDict['k1b.H2O']), label='k1b.H2O',
                color=colors[1], alpha=0.3)


    # Beautify and then save
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    kwargs = {'bbox_inches':'tight'}
    AC.save_plot(title=savetitle, **kwargs)
    plt.close('all')

    # - And the branching ratio
    fig, ax = plt.subplots()
    savetitle = 'ARNA_NO_HO2_Arrhenius_branching_ratio'
    if ArrheniusPlt:
        savetitle += '_1overK'

    # Ratio
    plt.scatter(xvalues, k1/k1, label='k1/k1', color=colors[0])
    plt.scatter(xvalues, RateDict['k1b']/k1, label='k1b/k1', color=colors[1])
    plt.scatter(xvalues, (k1*RateDict['Beta'])/k1,
                label='k1*Beta', color=colors[2])

    # Beautify and then save
    plt.xlabel(xlabel)
    plt.ylabel('Branching ratio')
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    kwargs = {'bbox_inches':'tight'}
    AC.save_plot(title=savetitle, **kwargs)
    plt.close('all')

    # - Just plot k1b
    fig, ax = plt.subplots()
    savetitle = 'ARNA_NO_HO2_Arrhenius_k1b'
    if ArrheniusPlt:
        savetitle += '_1overK'
#    plt.scatter(xvalues, RateDict['k1b'], label='k1b', color=colors[1])
    plt.scatter(xvalues, (RateDict['Beta']*k1),
                label='k1*Beta', color=colors[2])

    # Beautify and then save
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    kwargs = {'bbox_inches':'tight'}
    AC.save_plot(title=savetitle, **kwargs)
    plt.close('all')

    # Check ratios of rates

    # - And the branching ratio
    fig, ax = plt.subplots()
    savetitle = 'ARNA_NO_HO2_Arrhenius_normalised_rate_fractions'
    if ArrheniusPlt:
        savetitle += '_1overK'
#    plt.scatter(xvalues, RateDict['k1a']/RateDict[var], label='k1b',
#                color=colors[1])

    # Plot all the rates
    vars2plot = ['k1b.beta', 'k1b.H2O', 'k1b.H2O.cap']
    for n, var in enumerate( vars2plot ):
        plt.scatter(xvalues, (RateDict[var]/RateDict['k1a']),
                    label=var, color=colors[n+1])

    # Beautify and then save
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    kwargs = {'bbox_inches':'tight'}
    AC.save_plot(title=savetitle, **kwargs)
    plt.close('all')



def Calculate_Butkovskaya_channel_spit_P_T(TEMP=None, H2O=None,
                                           PRESS_TORR=None):
    """
    Calculate the HO2+NO rates following Butkovskaya et al [2005,2007, 2009]
    """
    # Base rate
    k1a = AC.GCARR(3.30E-12, 0.0E+00, 270.0, TEMP=TEMP)
    # rate k1b ()
    # NOTE: In the 223- 300 K range the rate constant of the HNO3-forming
    #       channel can be expressed as [Butkovskaya et al 2005]
    k1b = 6.4E-17 * np.exp( (1644. / TEMP))
    # Beta
    Beta = 530. / TEMP + 6.4E-4 * PRESS_TORR - 1.73
    k1b_beta = k1a * (Beta/100.)
    # Calculate Enhancement
    F = 1.0 + ( 2E-17 * H2O )
    F_capped = F.copy()
    F_capped[F_capped>8.] = 8.0

    #
    k1b_H2O = k1a * (Beta/100.) * F
    k1b_H2O_cap = k1a * (Beta/100.) * F_capped
    # return the calculated rates
    d = {
    'k1a':k1, 'Beta': Beta, 'k1b': k1b,
    'k1b.beta': k1b_beta,
    'k1b.H2O': k1b_H2O,
    'k1b.H2O.cap': k1b_H2O_cap,
    }
    return d


def plt_NO_HO2_channel_verues_H2O(rm_strat=True):
    """
    Plot up the enhancement seen in the HNO3 channel [Butkovskaya et al 2009]
    """
    # Get H2O, PRESS, T, etc from GEOS-Chem output
    d = get_GEOSChem_StateMet_vars4rxns(AvgOverTime=True)
#    PRESS_TORR = d['PRESS_TORR']
#    PRESS = d['PRESS']
    H2O = d['H2O']
#    TEMP = d['TEMP']
    TEMP = 298

    # General rate for HO2+NO=>NO+OH
    k1 = AC.GCARR(3.30E-12, 0.0E+00, 270.0, TEMP=TEMP)
    # Base rate for HO2+NO+H2O=>HNO3+H2O
    k1bw = 6E-13  #
    #
    Butkovskaya2005k = 6.4E-17 * np.exp( (1644/TEMP) )

    # Now consider NO + HO2 + H2O channel
    F = 1 + (2E-17 * H2O)
    # Only consider enhancement values constrained by observations (up to ~8)
#    F[ F > 8 ] = 8.0

    # - Now plot up the rates and enhancement factors vs. [H2O]
    context = "paper"
    font_scale = 1
    import seaborn as sns
    sns.color_palette('colorblind')
    sns.set_context(context)
    plt.close('all')
    # Set colours
    colors = AC.get_CB_color_cycle()

    # Plot up scaling factor
    plt.plot(H2O, F, color=colors[1], label='Butkovskaya2009')
    title = 'Scaling factor of NO+HO2+H2O=>HNO3+H2O'
    title += '\n at GEOS-Chem tropospheric [H2O]'
    plt.title(title)
    plt.xlabel('Water concentration ([H2O], molecules cm$^{-3}$)' )
    plt.ylabel('Enhancement factor')
    plt.vlines(x=4E17, ymin=min(F), ymax=max(F), label='50% humidty', ls='--',
               color=colors[4])
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    kwargs = {'bbox_inches':'tight'}
    AC.save_plot(title='ARNA_NO_HO2_channel_HNO3_factor')
    plt.close('all')

    # Plot up rate constant
    plt.plot(H2O, k1bw*F, color=colors[2], label='Butkovskaya2009')
    title = 'Rate constant of NO+HO2+H2O=>HNO3+H2O'
    title += '\n at GEOS-Chem tropospheric [H2O]'
    plt.title(title)
    plt.xlabel('Water concentration ([H2O], molecules cm$^{-3}$)' )
    plt.ylabel('Reaction rate (cm3 molecules$^{-1}$ s$^{-1}$)')
    plt.hlines(y=k1bw, xmin=min(H2O), xmax=max(H2O),
               label='k(1b) (HO2+NO+H2O=>HNO3)',
               ls='--', color=colors[5],)

    plt.hlines(y=Butkovskaya2005k, xmin=min(H2O), xmax=max(H2O),
               label='k(1b) [Butkovskaya2009]',
               ls='--', color=colors[3],)

    plt.hlines(y=k1, xmin=min(H2O), xmax=max(H2O),
               label='k1 (NO+HO2=>OH+NO2)', ls='--',
               color=colors[6])
    ymin, ymax = plt.ylim()
    plt.vlines(x=4E17, ymin=ymin, ymax=ymax, label='50% humidty',
             ls='--', color=colors[4])
#    plt.legend()
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    kwargs = {'bbox_inches':'tight'}
    AC.save_plot(title='ARNA_NO_HO2_channel_HNO3_rate', **kwargs)
    plt.close('all')


def plt_OH_NO2_channel_verues_H2O():
    """
    Plot Explore the OH + NO2 channel
    """
    # Get model runs to use for getin tropospheric H2O values
    RunDict = ar.get_dict_of_GEOSChem_model_output(RunSet='PostHONO',
                                                   GC_version='v12.9',
                                                   res='4x5',
                                                   folder4netCDF=True)
    wd = RunDict['min4pptHONO']
    StateMet = AC.get_StateMet_ds(wd=wd)
    StateMet = AC.add_molec_den2ds(StateMet)
    MolecVar = 'Met_MOLCES'
    TEMP = StateMet['Met_T'].mean(dim='time').copy().values.flatten()
    idx = np.argsort(TEMP)
    NUMDEN = StateMet[MolecVar].mean(dim='time').copy().values.flatten()
    ConcAC = AC.GetConcAfterChemDataset(wd=wd)
    SpecConc = AC.GetSpeciesConcDataset(wd=wd)
    # Get Tropospheric H2O, OH and NO2 values
    # H2O
    ds = AC.get_GEOSChem_H2O(wd=wd)
    NewVar = 'H2O'
    ds = ds[[NewVar]].mean(dim='time').copy()
    H2O = ds[NewVar].values.flatten()
    # N2 & O2
    O2 = NUMDEN*0.2095
    N2 = NUMDEN*0.78084
    # OH
    OH = ConcAC['OHconcAfterChem'].mean(dim='time').copy()
    OH = OH.values.flatten()
    # NO2
    NO2 = (SpecConc['SpeciesConc_NO2']*StateMet[MolecVar])
    NO2 = NO2.mean(dim='time').copy()
    NO2 = NO2.values.flatten()
    print([i.shape for i in (N2, O2, H2O, NUMDEN)])

    # - Amedro et al 2022
    Amedrok = calc_Amedro_rate_constant(N2=N2, O2=O2, H2O=H2O, TEMP=TEMP)
    # Plot up
    plt.scatter(TEMP[idx], Amedrok[idx], )
    plt.legend()
    AC.save_plot(title='ARNA_OH_NO2_rate_Amedro')
    plt.close('all')

    # - Base values
    A0 = 1.80E-30
    B0 = 3.0E+00
    C0 = 0.0
    A1 = 2.8E-11
    FV = 0.6
    k = AC.GCJPLPR(A0=A0, B0=B0, A1=A1, C0=C0, FV=FV, TEMP=TEMP)

    # Plot Base
    plt.scatter(TEMP[idx], k[idx], )
    plt.legend()
    AC.save_plot(title='ARNA_OH_NO2_rate_BASE')
    plt.close('all')

    # Plot

    # -  Fred values

#        k*(1+0.25*H2Ot/5.75e17)


    AC.save_plot(title='ARNA_OH_NO2_rate_Fred')
    plt.close('all')

    #

    # plot for the different H2O values

#        vals = k*(1+0.25*H2Ot/5.75e17)

    # Plot K


    # Plot resultant rate


    #



#    , 0.0, 0.0, 0.6, 0.0, 0.0

    # -

#    k*(1+0.25*H2Ot/5.75e17)





def calc_Amedro_rate_constant(N2=None, O2=None, H2O=None,
                              NUMDEN=1E4, TEMP=298.0, PRESS=1000.0):
    """
    Reproduce the FORTRAN calculation here
    """
    # K0 values by species
    K0N2 = 2.6E-30
    K0O2 = 2.0E-30
    K0H2O = 15.9E-30
    # And K infinity
    KBINF = 6.3E-11
    # Temperature dependence
    Fc = 39E-2
    Tm = 3.6
    Tq = 3.6
    To = 3.4
    n = 0.0
    # Calculate N2 and O2 locally if nor provided
    if isinstance(N2, type(None)):
        N2 = NUMDEN*0.78084
    if isinstance(O2, type(None)):
        O2 = NUMDEN*0.2095

    # Note this approach considers the presence of H2O, N2, O2
    # First line of equation
    PARTI_I = ( ( N2*K0N2*(TEMP/300.)**-Tm   ) +
             ( O2*K0O2*(TEMP/300.)**-Tq   ) +
             ( H2O*K0H2O*(TEMP/300.)**-To )   )
    PARTI_II = ( ( N2*K0N2*(TEMP/300.)**-Tm   ) +
              ( O2*K0O2*(TEMP/300.)**-Tq   ) +
              ( H2O*K0H2O*(TEMP/300.)**-To )   )
    PARTI_III = ( ( N2*K0N2*(TEMP/300.)**-Tm   ) +
               ( O2*K0O2*(TEMP/300.)**-Tq   ) +
               ( H2O*K0H2O*(TEMP/300.)**-To )   )
    k = PARTI_I * ( NUMDEN*KBINF*(TEMP/300.)**-n )
    # Second line
    k = k / (PARTI_II*NUMDEN + (KBINF*(TEMP/300.)**-n))

    # Calculate broadening factor
    LOGF = np.log10( (PARTI_III*NUMDEN)/ (KBINF*(TEMP/300.)**-n) )
    LOGF = 1 + ( LOGF / (75E-2 - 1.27 *  np.log10(Fc) ) )**2
    LOGF = np.log10( Fc ) / LOGF

    # Apply calculated pressure boardening
    k = k * 10.**(LOGF)
    return k


def do_analysis_of_4pptv_HONO_world(RunDict=None,
                                    RunSet='PostHONO', GC_version='v12.9',
                                    res='4x5', dates2use=None,
                                    use_REF_wd4Met=True, REF_wd=None,
                                    REF1='Base',  DIFF='min4pptHONO'):
    """
    Do general analysis of a model world where HONO is a minimum of 4 pptv
    """
    # Model runs to use?
    if isinstance(RunDict, type(None)):
        RunDict = ar.get_dict_of_GEOSChem_model_output(RunSet=RunSet,
                                                       GC_version=GC_version,
                                                       res=res,
                                                       folder4netCDF=True)
    # Set dates to use (as model run progresses)
    if isinstance(dates2use, type(None)):
        dates2use = [
        # First month of spin up
        datetime.datetime(2018, 1, 1),  # Initial month of spin up
        ]

    # -- Stats on model runs
    # Get generic stats
    extra_burden_specs = ['NOx', 'NIT-all', 'HNO2',  'NOy', 'HNO3', ]
    extra_surface_specs = extra_burden_specs
    # Get extra stats relate to HONO
    DiagVars = [
        'CO-lifetime', 'OH-lifetime', 'NOx-lifetime', 'NO2-lifetime',
        'CH4-lifetime',
#        'O3-lifetime',
        #    'Ox-lifetime',
        'HNO2-lifetime',
        'HNO2:NOx', 'HNO2:HNO3',
        'HNO3:NOx', 'HNO3:NO2',
        'NIT:HNO3', 'NIT:NOx', 'NO:NO2',
        #    'OH-trop-avg',
        'Cl-trop-avg', 'HNO2-trop-avg', 'NOx-trop-avg'
    ]

    # Quick test for model output
    for key in RunDict.keys():
        wd = RunDict[key]
        print(key, RunDict[key])
        try:
            df = GetModelOutputFileDates(wd=wd, ReturnDates=False)
            print(df.count() )
        except:
#    except OutOfBoundsDatetime:
            pass
        print()


#    df = AC.get_general_stats4run_dict_as_df(run_dict=RunDict,
    use_time_in_trop = True
    use_REF_wd4Met = True
    REF_wd = RunDict[DIFF]
    df = AC.get_stats4RunDict_as_df(RunDict=RunDict,
                                    dates2use=dates2use,
                                    use_REF_wd4Met=use_REF_wd4Met,
                                    REF_wd=REF_wd,
                                    extra_burden_specs=extra_burden_specs,
                                    extra_surface_specs=extra_surface_specs,
                                    DiagVars=DiagVars,
                                    use_time_in_trop=use_time_in_trop,
                                    verbose=False,
                                    debug=False,
                                    )

    # Change Row ordering
    df = df.T
    FirstDiags = [
        'CH4-lifetime (years)', 'O3 burden (Tg)', 'O3 surface (ppbv)',
        'NOy burden (Tg)',
        'HNO2 surface (ppbv)', 'HNO2:NOx', 'NOy surface (ppbv)',
        'OH surface (molec/cm3)', 'NIT-all burden (Tg)',
        'HNO3:NOx', 'HNO2:NOx', 'NIT:NOx',

    ]
    OtherDiags = [i for i in df.columns if (i not in FirstDiags)]
    df = df[FirstDiags+OtherDiags]
    df = df.T

    # Save
    try:
        column_ordering = [
            'Base',
            'min4pptHONO',
            'NIThv',
            'NIThvCap',
            'N2O5',
            'OH+NO2',
            'HalNitratesx10',
    #        'min4pptHONOday'
    #        '4pptHONO',
        ]
        column_ordering +=[i for i in df.columns if (i not in column_ordering)]
        df = df[column_ordering]
    except KeyError:
        pass
    SaveName = 'ARNA_Key_Stats_on_model_runs_post_HONO'
    SaveNameAll = '{}_{}_{}'.format(SaveName, RunSet, GC_version)
    SaveNameAll = AC.rm_spaces_and_chars_from_str(SaveNameAll)
    df.to_csv('{}{}'.format(SaveNameAll, '.csv'))

    # Save out ratio values
    cols2use = [i for i in df.columns if (i != REF1)]
    df_diff = df[cols2use]
    for __col in cols2use:
        df_diff.loc[:, __col] = df_diff.loc[:, __col] / df[REF1]
    df_diff.to_csv('{}{}'.format(SaveNameAll, '_ratio.csv'))

    # Save out core stats for manuscript table
    corestats = []
    df[corestats].to_csv(SaveName)

    # Return the dataframe for ipython data analysis
    return df


def GetModelOutputFileDates(wd=None,
                                    FileStr=None,
                                    FileStrs=[],
                                    CheckForAllOutput=True,
                                    ReturnDates=True,
                                    IncBoundaryConditions=False,
                                    IncHEMCO=False,
#                                    NetCDFSuffix='.nc4'
                                    freq='24H',
                                    ):
    """
    Return list of dates with complete NetCDF output files
    """
    # Use requested categories as a list
    if len(FileStrs) >1:
        pass
    else:
        FileStrs = [FileStr]

    # Check number of complete output files for all file string types?
    if CheckForAllOutput:
        #
        FilesInDir = glob.glob('{}/{}'.format(wd, '*GEOSChem.*') )
        categories = list(set([ i.split('.')[-3] for i in FilesInDir]))
        FileStrs = categories
        #
    # Include (or remove) BoundaryConditions file string
    if IncBoundaryConditions:
        FileStrs += ['BoundaryConditions']
    else:
        try:
            FileStrs.pop( FileStrs.index('BoundaryConditions') )
        except ValueError:
            pass
    # As above, but for HEMCO
    if IncHEMCO:
        FileStrs += ['HEMCO_diagnostics']
    else:
        try:
            FileStrs.pop( FileStrs.index('HEMCO_diagnostics') )
        except ValueError:
            pass

    def convert_datetime2days(input):
        return datetime.datetime(*input.timetuple()[:3])

    # Now extract dates there are files for
    dtStr = 'datetime'
    MinDate = datetime.datetime(5000, 1, 1)
    MaxDate = datetime.datetime(1000, 1, 1)
    d = {}
    for FileStr in FileStrs:
        FilesNames = glob.glob( '{}/*.{}.*'.format(wd, FileStr) )
        FilesNames = list(sorted(FilesNames))
#        df['RawStrs'] = dates
        # Select format
        if ('HEMCO' in FileStr):
            format = '%Y%m%d%H%M'
            DateStrLength = 15
        else:
            format = '%Y%m%d_%H%Mz'
            DateStrLength = 14
        # Strip to get dates
        FilesNames = [i.split('/')[-1] for i in FilesNames]
        dates = [i.split('.')[-2][-DateStrLength:] for i in FilesNames]
        dates = pd.to_datetime(dates, format=format).values
        dates = [convert_datetime2days(i) for i in AC.dt64_2_dt(dates)]
        dates = list(sorted(set(dates)))
        d[FileStr] = pd.DataFrame(dates, index=dates)

        # Minimum date
        __MinDate = min(dates)
        if __MinDate < MinDate:
            MinDate = __MinDate
        # Maximum date
        __MaxDate = max(dates)
        if __MaxDate > MaxDate:
            MaxDate = __MaxDate

#        FileStr
#        df.loc[FileStr] = S
#        print(len(files2use))
#        if verbose:
#            print(files2use)
#        print()

    # Check all the same dates are available for all output
    index = pd.date_range(MinDate, MaxDate, freq=freq)
    df = pd.DataFrame(index=index)

    #
    for key in d.keys():
#        print(key)
        try:
#            df = pd.concat([df,  d[key] ], axis=1, join="inner")
            df[key] = d[key]
        except:
            print(key)


    # Only retain rows where there are values
    df = df.dropna(how='all')

    if ReturnDates:
        return df.index.values
    else:
        return df


def plt_spatial_changes_in_4pptv_HONO_world(pcent=True,
                                            REF1='Base',
                                            DIFF='4pptHONO',
                                            RunSet='PostHONO',
                                            res='4x5',
                                            dates2use=None,
                                            GC_version='v12.9',
                                            RunDict=None,
                                            PltAllVars=False):
    """
    Plot up changes in key metrics spatially and zonally
    """
    # Retrieve dictionary of model runs and the location of the output data
    if isinstance(RunDict, type(None)):
        RunDict = ar.get_dict_of_GEOSChem_model_output(RunSet=RunSet,
                                                       GC_version=GC_version,
                                                       res=res,
                                                       folder4netCDF=True)
        RunDict = {REF1: RunDict[REF1], DIFF: RunDict[DIFF], }
    # Set dates to use (as model run progresses)
    if isinstance(dates2use, type(None)):
         dates2use = [ datetime.datetime(2018, 1, 1 ), ]

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

    # Rename variables in dataset
    prefix = 'SpeciesConc_'
    OldNames = [i for i in ds.data_vars if (prefix in i)]
    NewNames = [i.split(prefix)[-1] for i in OldNames]
    for key in dsD.keys():
        ds = dsD[key]
        ds = ds.rename(name_dict=dict(zip(OldNames, NewNames)))
        dsD[key] = ds

    # Get MetState object
    StateMet = AC.get_StateMet_ds(wd=RunDict[REF1])

    # Get HOx values as xr.Dataset too
    CACsuffix = 'concAfterChem'
    HOx_specs = 'HO2', 'OH', 'HOx'
    dsCAC = {}
    for key in RunDict.keys():
        ds = GetConcAfterChemDataset(wd=RunDict[key],
                                        dates2use=dates2use)
        # Add
        ds = AC.add_HOx_to_CAC_ds(ds, UpdateHOxUnits=True, StateMet=StateMet)
        # Update CAC names
        NewNames = [i.split(CACsuffix)[0] for i in ds.data_vars]
        ds = ds.rename(name_dict=dict(zip(ds.data_vars, NewNames)) )
        dsCAC[key] = ds

    # -- Graphics on model runs
    SaveStr = 'ARNA_spatial_HONO4pptv_{}_vs_{}'
    savetitle = SaveStr.format(REF1, DIFF)
    specs2plot = ['O3', 'HNO2',  'HNO3', 'CO', 'OH'] + families2use
    if PltAllVars:
        specs2plot = list(set(list(ds.data_vars)+specs2plot)  )
        SaveStr += '_AllVars'
    # Updates setings for percentage plotting
    if pcent:
        kwargs = {'vmin': -100, 'vmax': 100}
        savetitle += '_pcent'
    else:
        kwargs = {}

    # Plot up HONO differences (surface)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    for var2plot in specs2plot:

        # Select the reference and difference datasets to use
        if (var2plot in HOx_specs):
            # NOTE: Averaging all arrays over time for now...
            ds1 = dsCAC[REF1][[var2plot]].copy().mean(dim='time')
            ds2 = dsCAC[DIFF][[var2plot]].copy().mean(dim='time')
        else:
            # NOTE: Averaging all arrays over time for now...
            ds1 = dsD[REF1][[var2plot]].copy().mean(dim='time')
            ds2 = dsD[DIFF][[var2plot]].copy().mean(dim='time')
        if 'time' in list(StateMet.dims):
            StateMet = StateMet.mean(dim='time')

        # Calculate difference
        if pcent:
            ds2plot = (ds2 - ds1) / ds1 * 100
        else:
            ds2plot = ds2 - ds1

        # Select the surface
        __ds2plot = ds2plot.sel(lev=ds2plot.lev.values[0])

        # Use the same colourmaps for all plots
        vmin = ds2plot[var2plot].values.min()
        vmax = ds2plot[var2plot].values.max()
        cmap = AC.get_colormap(np.array([vmin,vmax]))
        kwargs = {'cmap': cmap}
        # If no changes is seen (e.g. )
        if (vmin == np.NaN) and (vmax == np.NaN):
            if verbose or debug:
                print('WARNING: Skipping ploting of {}'.format(var2plot))
            continue  # continue past here

        # Units and scaling for plot?
        if pcent:
            scaleby = 1
            units = '%'
        else:
            units, scaleby = AC.tra_unit(var2plot, scale=True)

        TitleStr = "{} diff of {} ('{}' vs. '{}') in '{}'"
        title = TitleStr.format('Surface', var2plot, DIFF, REF1, units)

        AC.quick_map_plot(__ds2plot[[var2plot]]*scaleby, var2plot=var2plot,
                          verbose=verbose, title=title, **kwargs)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # zonal
        fig, ax = plt.subplots()
        im = AC.ds2zonal_plot(ds2plot.copy()*scaleby, var2plot=var2plot,
                              StateMet=StateMet,
                              fig=fig, ax=ax, **kwargs)
        # Add a colourbar
        fig.colorbar(im, orientation="horizontal", pad=0.2, label=units)
        # Add a title
        title = TitleStr.format('Zonal', var2plot, DIFF, REF1, units)
        plt.title(title)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()


    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


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
    vmin = 0.1
    vmax = 1000
    norm = LogNorm(vmin=vmin, vmax=vmax)
    kwargs = {'norm': norm}
    #    kwargs = {'log': True}
    ds = ds[[var2plot]].sum(dim='lev') * scale
    # TODO: make plot logarithmic, clear colour bar axis
    AC.quick_map_plot(ds, var2plot='SpeciesConc_HNO2', show_plot=False,
                      savename=savename,
                      #                      save_plot=False,
                      verbose=verbose, kwargs=kwargs)



def plot_SpecSubsetConc_comparisons()
    """
    """

    # --- Datasets to use?
    RunDict = {}
    RunDict['Iso.Unlimited']  = '/mnt/lustre/users/ts551/GC/rundirs/P_ARNA/geosfp_4x5_aciduptake.v12.9.0.ARNA.Isotherm.Diags.v9.Iso.UnlimitedAll/OutputDir/'

    RunDict['Base'] = '/mnt/lustre/users/ts551/GC/rundirs/P_ARNA/geosfp_4x5_aciduptake.v12.9.0.ARNA.Isotherm.Diags.v9.Base/OutputDir/'


    # --- ARNA
    dsD = {}
    for key in RunDict.keys():
        folder = RunDict[ key ]

        # Get the data of species to plot
        # HNO2, NOx, NO, NO2, CO, ozone,
        data = get_model_noon_vertical_conc(folder)

        dsD[key] = data


    # Setup a PDF to store values
    savetitle = 'ARNA_SpecConcSub_Noon_mean'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)


    # Units to use?
    units_d = {
        'CO': 'ppbv', 'O3': 'ppbv', 'NO2': 'pptv', 'NO': 'pptv', 'NOx': 'pptv',
        'HNO2': 'pptv', 'HONO': 'pptv',
        'NIT': 'pptv', 'NITs': 'pptv', 'SO4s': 'pptv', 'SO4': 'pptv',
        'NH4': 'pptv',
        'SO4-all': 'pptv', 'NIT-all': 'pptv',
    }

    # Plot up
    vars2plot = 'HNO2', 'NOx', 'NO', 'NO2', 'CO'
    for var2plot in vars2plot:

        # Loop model data
        for key in dsD.keys():
            # Get units and scaling
            units = units_d[var2plot]
            scaleby = AC.get_unit_scaling(units)

            # Select data to plot
            X = dsD[key][var2plot] * scaleby
            Y = dsD[key].lev

            # plot up
            plt.plot(X, Y, label=key)

        # Beautify plot
        plt.title( var2plot)

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def get_model_noon_vertical_conc( folder=None ):
    """
    Retrieve the vertical noon concentration
    """

    # Location or area? (could average a lat/lon extent or select location)

    # Select the dates to extract vertical profile for
    campaign = 'ARNA'
    if campaign == 'ARNA':
        dates2use = [datetime.datetime(2020, 2, 4)]
        dates2use += [datetime.datetime(2020, 2, 5)]
        dates2use += [datetime.datetime(2020, 2, 6)]
        dates2use += [datetime.datetime(2020, 2, 7)]
        dates2use += [datetime.datetime(2020, 2, 8)]
        dates2use += [datetime.datetime(2020, 2, 11)]
        dates2use += [datetime.datetime(2020, 2, 12)]
        dates2use += [datetime.datetime(2020, 2, 13)]

    # open just the dates of the campaign
    file_str='GEOSChem.SpeciesConcSubset.*.nc4'
    ds = AC.get_GEOSChem_files_as_ds(file_str=file_str, dates2use=dates2use,
                                     wd=folder)

    # Select the average noon profiles
    if campaign == 'ARNA':
        TZ = -1
    dt = AC.dt64_2_dt( ds.time.values  )
    dt = [AC.add_hrs(i, TZ ) for i in dt]
    ds = ds.assign( {'time': dt} )
    __bool = ds['time.hour'] == 12
    ds = ds.isel(time=__bool).mean(dim='time')

    # Select the location
    d = ar.get_analysis_region('Cape_Verde_Flying')
    x0, x1, y0, y1 = d['x0'], d['x1'], d['y0'], d['y1']
    # Set values region
    bool1 = ((ds.lon >= x0) & (ds.lon <= x1)).values
    bool2 = ((ds.lat >= y0) & (ds.lat <= y1)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    # Average over location
    ds = ds.mean(dim=('lat', 'lon'))

    # Add in families afterwards
    ds = AC.AddChemicalFamily2Dataset(ds, fam='NOx')
    ds = AC.AddChemicalFamily2Dataset(ds, fam='NIT-all')

    # remove prefix
    SCprefix = 'SpeciesConc_'
    variables = [i for i in ds.data_vars if SCprefix in i]
    species = [i.split(SCprefix)[-1] for i in variables]
    name_dict = dict(zip(variables, species))
    ds = ds.rename( name_dict=name_dict )
    ds = ds[species] # Only return species


    # return the dataset
    return ds


def calc_pressure_feild4ds():
    """
    Calculate offline the pressure using the NetCDF core variables

    Notes
    ---
     - The calculation to perform to get box top or bottom pressure

     Pbot = ( hyai(L  ) * PS(I,J) ) + hybi(L  )
     Ptop = ( hyai(L+1) * PS(I,J) ) + hybi(L+1)

     - More detail on the pressure calculation  http://wiki.seas.harvard.edu/geos-chem/index.php/Overview_of_History_diagnostics
    """

    pass







if __name__ == "__main__":
    main()

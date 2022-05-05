#!/usr/bin/python
"""
Driver for analysis of "4 pptv HONO world" following the ARNA campaign
"""
import arna as ar
import sys

from AC_tools import GetSpeciesConcDataset, AddChemicalFamily2Dataset, species_mass, get_ref_spec, get_Gg_trop_burden, get_StateMet_ds, tra_unit, get_avg_2D_conc_of_X_weighted_by_Y, GC_var, get_ProdLoss_ds, add_molec_den2ds, rm_fractional_troposphere, constants, GetConcAfterChemDataset, get_StateMet_ds, get_DryDep_ds, get_WetLossConv_ds, get_WetLossLS_ds

from arna import get_local_folder, get_tags_for_NOx_HONO


def main():
    """
    Main driver for 4 pptv HONO world analysis
    """
    # --- Local settings to pass to all analysis/plotting functions
    RunSet = 'PostHONO'
    res = '4x5'
    GC_version = 'v12.9'
#    GC_version = 'v13.4'
#    dates2use = None
#    dates2use = [datetime.datetime(2018, 1+i, 1) for i in range(1)]
    # first 3 months
#    dates2use = [datetime.datetime(2018, 1+i, 1) for i in range(3)]
    dates2use = [datetime.datetime(2018, 1+i, 1) for i in range(10)]
    # spin up year
    dates2use = [datetime.datetime(2018, 1+i, 1) for i in range(12)]
    # 6 months, after 6 months of spin up.
#    dates2use = [datetime.datetime(2018, 7+i, 1) for i in range(6)]
    # Spun up year
    dates2use = [datetime.datetime(2019, 1+i, 1) for i in range(12)]

    REF1 = 'Base'
#    DIFF = '4pptHONO'
    DIFF = 'min4pptHONO'
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
    RunStr4 = RunBaseStr1 + 'ARNA.Isotherm.Diags.v6'
    BaseRunStr5 =  RunStr4 +  '.HO2nNOv2.5'
    RunStr5 =  '{}{}{}/OutputDir/'.format(RunRoot, BaseRunStr5, '')
    RunStr6 =  '{}{}{}/OutputDir/'.format(RunRoot, BaseRunStr5, '.IsopH4')
    RunStr7 =  '{}{}{}/OutputDir/'.format(RunRoot, BaseRunStr5,
                                          '.Iso.UnlimitedpH')
    RunStr8 =  '{}{}{}/OutputDir/'.format(RunRoot, RunStr4, '.Ye2017')
    RunStr9 =  '{}{}{}/OutputDir/'.format(RunRoot, RunStr4, '.Ye2017online')

    RunDict2 = {
    'Base':RunDict['Base'], 'min4pptHONO': RunDict['min4pptHONO'],
    'HO2+NOv2.1': RunStr1,
    'HO2+NOv2.5': RunStr5,
    'Isov5.pH4.HO2+NOv2.5': RunStr6,
    'Isov5': RunStr3,
#    'Iso.HO2+NO': RunStr2, # WARNING: This is using HO2+NO code with a bug
    'Iso.Unlimited': RunStr7,
    'Ye17': RunStr8,
    'Ye17.online': RunStr9,
    }


    RunDict = RunDict2

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


    # Plot up HONO production routes - N/A
    # See functions in do_post_campaign_NOx_analysis.py
    # (e.g.'plt_key_NOx_budget_terms' )
    plt_HO2_NO_branching_ratio_spatially()

    # Plot up the NO + HO2 + H2O => HNO3 channel
    plt_NO_HO2_channel_verues_H2O()
    # And the OH + NO2 channel...
    plt_OH_NO2_channel_verues_H2O()


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
        'OH surface (molec/cm3)', 'NIT-all burden (Tg)'
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
        ds = AC.add_HOx_to_CAC_ds(ds, UpdateHO2units=True, StateMet=StateMet)
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
    kwargs = {'log': True}
    ds = ds[[var2plot]].sum(dim='lev') * scale
    # TODO: make plot logarithmic, clear colour bar axis
    AC.quick_map_plot(ds, var2plot='SpeciesConc_HNO2', show_plot=False,
                      savename=savename,
                      #                      save_plot=False,
                      verbose=verbose, kwargs=kwargs)


if __name__ == "__main__":
    main()

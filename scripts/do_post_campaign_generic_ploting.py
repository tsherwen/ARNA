#!/usr/bin/python
"""
Driver for generic plotting functions following the ARNA campaign
"""
import arna as ar


def main():
    """
    Main driver function
    """
#    from arna import get_FAAM_core4flightnum, get_filters_data4flight, get_GEOSCF4flightnum, add_derived_FAAM_flags2df4flight, get_GEOSChem4flightnum, plt_flightpath_spatially_over_CVAO, get_CIMS_data4flight

    # Seaborn context for plots?
#    context = 'talk'
    context = 'paper'

    # - Plot up all comparisons by altitude together
    ar.plt_comp_by_alt_4ARNA_together(context=context,
                                      res='0.25x0.3125',
                                      just_SLR=True,)
    ar.plt_comp_by_alt_4ARNA_together(context=context,
                                      res='0.25x0.3125',
                                      just_SLR=False)
    ar.plt_comp_by_alt_4ARNA_together(context=context,
                                      res='4x5', RunSet=None,
                                      just_SLR=False,)
    # temp for testing
    flight_nums = []
    RunSet = None
    res = '4x5'
    NOxAsLog = True
#    CoreRunsOnly = False
    CoreRunsOnly = True
    savetitle = 'ARNA_altitude_binned_combined_file_{}'.format(res)
    ar.plt_comp_by_alt_4ARNA_together(context=context,
                                      res=res, RunSet=RunSet,
                                      flight_nums=flight_nums,
                                      savetitle=savetitle,
                                      just_SLR=False,
                                      NOxAsLog=NOxAsLog,
                                      CoreRunsOnly=CoreRunsOnly,
                                      debug=True)
    # Output the same bulk plots for the ACID runs
    RunSet = 'ACID'
    res = '4x5'
    savetitle = 'ARNA_altitude_binned_combined_file_{}_{}'
    savetitle =  savetitle.format(res, RunSet)
    ar.plt_comp_by_alt_4ARNA_together(context=context,
                                      res=res, RunSet=RunSet,
                                      flight_nums=flight_nums,
                                      savetitle=savetitle,
                                      NOxAsLog=NOxAsLog,
                                      CoreRunsOnly=CoreRunsOnly,
                                      just_SLR=False, debug=True)

    # The same plots as above, but split by their own PDF file..
    # NOTE: below function fails with a ValueError
    ar.plt_comp_by_alt_4ARNA_all(just_SLR=True, context=context,
                                 RunSet='FP-Nest',inc_GEOSChem=True,
                                 just_plot_GEOS_Chem=True,
                                 res='0.25x0.3125', close_pdf=True)

    # Plot up data for SLRs with and without dust
#    ar.plt_comp_by_alt_4ARNA_all(just_SLR=False, context=context)
#    ar.plt_comp_by_alt_4ARNA_all(just_SLR=True, context=context)
    ar.plt_comp_by_alt_4ARNA_all(just_SLR=False, context=context,
                                 RunSet='FP-Nest',inc_GEOSChem=True,
                                 just_plot_GEOS_Chem=True,
                                 res='0.25x0.3125', close_pdf=True)

    ar.plt_comp_by_alt_4ARNA_flights_CIMS(context=context,
                                          RunSet='FP-Nest',
                                          res='0.25x0.3125',
                                          inc_GEOSChem=True,
                                          just_SLR=True)

    ar.plt_comp_by_alt_4ARNA_all_DUST(plt_model=False, context=context)
    ar.plt_comp_by_alt_4ARNA_all_DUST(plt_model=True, context=context)
    ar.plt_comp_by_alt_4ARNA_CIMS_all_DUST(context=context)

    # - Plot up core comparisons by flight as
    # As timeseries ...
    ar.plt_ts_comp4ARNA_flights(inc_GEOSChem=False, context=context)
    ar.plt_ts_comp4ARNA_flights(inc_GEOSChem=True, context=context,
                                just_plot_GEOS_Chem=True,
                                RunSet='FP-Nest', res='0.25x0.3125')
    # By altitude (and by flight)
    ar.plt_comp_by_alt_4ARNA_flights(context=context)
    ar.plt_comp_by_alt_4ARNA_flights(inc_GEOSChem=True, context=context,
                                     just_plot_GEOS_Chem=True,
                                     RunSet='FP-Nest', res='0.25x0.3125')

    # - Plot up ToF-CIMS data by flight as
    # As timeseries ...
    ar.plt_ts_comp4ARNA_flights_CIMS(context=context)
    ar.plt_ts_comp4ARNA_flights_CIMS(context=context,
                                     RunSet='FP-Nest',
                                     res='0.25x0.3125',
                                     inc_GEOSChem=True,
                                     flight_nums=flight_nums,
                                     LatVar='LAT',
                                     LonVar='LON',)
    # By altitude (and by flight)
    ar.plt_comp_by_alt_4ARNA_flights_CIMS(context=context)
    ar.plt_comp_by_alt_4ARNA_flights_CIMS(context=context,
                                          RunSet='FP-Nest',
                                          res='0.25x0.3125',
                                          inc_GEOSChem=True,
                                          )


    # - Plot up nitrate aerosol data by flight as
    # As timeseries ...
    ar.plt_ts_comp4ARNA_flights_filters(context=context)
    ar.plt_ts_comp4ARNA_flights_filters(context=context,
                                        RunSet='FP-Nest',
                                        res='0.25x0.3125',
                                        inc_GEOSChem=True,
                                        LatVar='LAT',
                                        LonVar='LON',)
    ar.plt_ts_comp4ARNA_flights_filters(context=context,
                                        res='4x5',
                                        inc_GEOSChem=True,
                                        LatVar='LAT',
                                        LonVar='LON',
                                        )


    # Plot up nitrate, JNIT, and their project
#    AC.mk_tri_NO3_JNIT_combination_plt()



    # - Plot up SWAS data by flight
    # Plot up SWAS data
    ar.plt_ts_comp4ARNA_flights_SWAS(context=context)


    # - Plot up PCASP/CDP data by flight as
    # NOTE: CAS data being ignored currently due to issue with mirror window
    # As timeseries ...
#    ar.plt_ts_comp4ARNA_flights_PCASP()


    # - Plot up velocity and Roll, amongst other core physical vars by flight
    # As timeseries ...
    ar.plt_ts_comp4ARNA_flights_PHYSICAL_VARS(context=context)
    ar.plt_ts_comp4ARNA_flights_PHYSICAL_VARS(context=context,
                                              just_plot_GEOS_Chem=True,
                                              inc_GEOSChem=True,
                                              res='0..25x0.3125',
                                              RunSet='FP-Nest')
    # Plot up the temperature data from Hannah Price
    # N/A? this is only for 2019. Data to be worked up for 2020.


    # - Plot up SWAS data by flight
    # Plot a comparison of NOy
    ar.plt_ts_comp4ARNA_flights_NOy_ALL(context=context)
    ar.plt_ts_comp4ARNA_flights_NOy_ALL(context=context,
                                        RunSet='FP-Nest',
                                        res='0.25x0.3125',
                                        inc_GEOSChem=True,
                                        LatVar='LAT',
                                        LonVar='LON',)

    # - Other misc. plotting tasks
#    explore_high_ozone_near_CVAO()
#    extract_GEOS54all_ARNA_flights()

    # Evaluate the high resolution modelling region
    ar.evaluate_regional_grid4GEOSChem()

    # Also plot up for related biomass-burning flights in MOYA campaign
    ar.plt_ts_comp4MOYA_flights()
    ar.plt_ts_comp4MOYA_flights_PHYSICAL_VARS()






def explore_high_ozone_near_CVAO():
    """
    High ozone seen in some locations just in nested GEOS output
    """
    import seaborn as sns
    import gc
    # - plot up locations of high ozone by flight.
    # Use the high res model
    RunSet = 'FP-Nest'
    res = '0.25x0.3125'
    # Which flights to plot? - Just use non-transit ARNA flights
    flights_nums = [218, 219, 220, 221, 222, 223, 224, 225,]
    flight_IDs = ['C{}'.format(i) for i in flights_nums]
    # Get model data
    dfs_mod_GC = {}
    for flight_ID in flight_IDs:
        dfs_mod_GC[flight_ID] = ar.get_GEOSChem4flightnum(flight_ID=flight_ID,
                                                          res=res,
                                                          RunSet=RunSet,)
    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        dfs_obs[flight_ID] = ar.get_FAAM_core4flightnum(flight_ID=flight_ID)

    # Now plot high ozone locations by flight
    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_flighttrack_high_ozone'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    for flight_ID in flight_IDs:
        #    for flight_ID in flight_IDs[-1:]:
        # Get dataframes for flight
        df_obs = dfs_obs[flight_ID]
        df_mod = dfs_mod_GC[flight_ID]
        # subselect locations where ozone > 70 ppbv
        df2plot = df_mod['FP-Nest']
        LatVar = 'LAT'
        LonVar = 'LON'
        cut_off = 70
        df2plot = df2plot.loc[df2plot['O3'] > cut_off*1E-09, :]
        title_str = 'Locations with modelled ozone >{} ppbv during {}'
        if (df2plot.shape[0] != 0):
            prt_str = 'Plotting {} for values >{} ppbv'
            print(prt_str.format(flight_ID, cut_off))

            title = title_str.format(cut_off, flight_ID)
            ar.plt_flightpath_spatially_over_CVAO(df=df2plot,
                                                  flight_ID=flight_ID,
                                                  title=title,
                                                  LatVar=LatVar, LonVar=LonVar)

            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        else:
            prt_str = 'WARNING: Not plotting for {} as no O3 ppbv values >{}'
            print(prt_str.format(flight_ID, cut_off))
    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')

    # - Check the average modelled ozone vertical profile at Cape verde

    #
    RunSet = 'FP-Nest'
    res = '0.25x0.3125'
    RunDict = ar.get_dict_of_GEOSChem_model_output(res=res, RunSet=RunSet)
    Folder = '{}/OutputDir/'.format(RunDict[RunSet])

    #
    ds = AC.get_GEOSChem_files_as_ds(wd=Folder)
#    StateMet = AC.get_StateMet_ds(wd=Folder)
    ModelAlt = AC.gchemgrid('c_km_geos5')
    prefix = 'SpeciesConc_'
    specs2plot = [i for i in ds.data_vars if prefix in i]
    specs2plot = [i.split(prefix)[-1] for i in specs2plot][::-1]
#    specs2plot = ['O3', 'CO', 'NO2'][::-1]

    # Setup PDF to save PDF plots to
    import seaborn as sns
    sns.set(color_codes=True)
    sns.color_palette('colorblind')
    sns.set_context(context)
    savetitle = 'ARNA_vertical_above_CVAO_GEOSChem_campaign_01'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    save2png = False
    show_plot = False
    font_scale = 1

    from funcs4obs import gaw_2_loc
    for spec in specs2plot[:100]:
#    for spec in specs2plot[-100:]:
#    for spec in specs2plot[-100:-90]:
        site = 'CVO'
        lat, lon, alt, TZ = gaw_2_loc(site)
        ds_tmp = ds[prefix+spec].sel(lat=lat, lon=lon, method='nearest')
        ds_tmp = ds_tmp.mean(dim='time')
        ds_tmp['lev'] = ModelAlt
        units, scalby = AC.tra_unit(spec, scale=True)
        ds_tmp *= scalby
        #  plot up...
        ds_tmp.plot(y='lev')
        plt.ylabel('{} ({})'.format('Altitude', 'km'))
        plt.xlabel('{} ({})'.format(spec, units))
        plt.title('Vertical Profile at CVAO during ARNA-2')
        plt.ylim(0, 15)
        if spec == 'O3':
            plt.xlim(-20, 200)
        elif spec == 'NO2':
            plt.xlim(-0.05, 0.2)
        elif spec == 'CO':
            plt.xlim(40, 160)
        if save2png:
            save_str = 'ARNA_vertical_above_CVAO_GEOSChem_{}'
            AC.save_plot(save_str.format(spec))
            AC.close_plot()
        else:
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            if show_plot:
                plt.show()
            plt.close()
        # Do some memory management...
        gc.collect()

    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def extract_GC_data4CVAO():
    """
    Extract model data for CVAO
    """
    # - Get photolysis surface data for Simone
    RunRoot = '/mnt/lustre/users/ts551/GC/rundirs/'
    folder = 'geosfp_4x5_standard.v12.9.0.BASE.2019.2020.ARNA.BCs.'
    folder += 'TEST.PF_Jrates.JVALS.GLOBAL/OutputDir/'
    folder = RunRoot + folder
    # Open the dataset
    ds = AC.GetJValuesDataset(wd=folder)
    # Extract for CVAO
    site = 'CVO'
    lat, lon, alt, TZ = gaw_2_loc(site)
    ds_tmp = ds.sel(lat=lat, lon=lon, method='nearest')
    ds_tmp = ds_tmp.sel(lev=ds_tmp.lev.values[0])
    # save out to csv.
    vars2use = [i for i in ds.data_vars if 'UV' in i]
    vars2use += [i for i in ds.data_vars if 'Jval_' in i]
    vars2use = list(set(vars2use))
    # Reduce to the variables of intereest and then save to disk
    ds_tmp = ds_tmp[vars2use].squeeze()
    del ds_tmp['lat']
    del ds_tmp['lon']
    del ds_tmp['lev']
    ds_tmp = AC.save_ds2disk_then_reload(ds_tmp)
    df = ds_tmp.to_dataframe()
    df.to_csv('GC_JValues_Collection_4x5_FP-GLOBAL_CVAO.csv')


def test_new_planeflight_Jrate_output():
    """
    Test online rxn. output via plane-flight diagnostic from GEOS-Chem
    """
    from funcs4obs import gaw_2_loc
    # - Setup sites to use
    df = pd.DataFrame()
    GAW_sites = [
        'ASK', 'BRW', 'CGO', 'CMN', 'CPT', 'CVO', 'JFJ', 'LAU', 'MHD', 'MLO',
        'MNM', 'NMY', 'SMO', 'SPO', 'THD'
    ]

    # - Get plane flight output
    RunRoot = '/mnt/lustre/users/ts551/GC/rundirs/'
    folder = RunRoot+'/geosfp_4x5_standard.v12.9.0.BASE.2019.2020.ARNA.BCs.TEST.PF_Jrates.REA.VI/'
    files2use = list(sorted(glob.glob(folder + '/TEST_1day/*plane*')))
    file2use = files2use[0]
    # Get Header infomation from first file
    vars, sites = AC.get_pf_headers(file2use, debug=debug)
    # Extract all points from file
    dfs = []
    for file2use in files2use:
        print(file2use)
        df, NIU = AC.pf_csv2pandas(file=file2use, vars=vars, epoch=True,
                                   r_vars=True)
        dfs += [df]
    # Append the dataframes together
    df = dfs[0].append(dfs[1:])
    df = AC.DF_YYYYMMDD_HHMM_2_dt(df, rmvars=None, epoch=False)
    df.index.name = None

    # Process and save csv files by date
    filename = 'GC_planeflight_data_FP-GLOBAL_JVALS_ARNA1_{}.csv'
    for file2use in files2use:
        date = file2use.split('plane.log.')[-1]
        print(file2use)
        vars, sites = AC.get_pf_headers(file2use, debug=debug)
        df, NIU = AC.pf_csv2pandas(file=file2use, vars=vars, epoch=True,
                                   r_vars=True)
        # Update the datetime index
        df = AC.DF_YYYYMMDD_HHMM_2_dt(df, rmvars=None, epoch=False)
        df.index.name = None
        df.to_csv(filename.format(date))

    # - Get output from NetCDF diagnostics
    # Get Get J-values
    folder = RunRoot+'/geosfp_4x5_standard.v12.9.0.BASE.2019.2020.ARNA.BCs.TEST.PF_Jrates.REA.III/OutputDir/'
    ds = AC.GetJValuesDataset(wd=folder)
    # And get StateMet
    FileStr = 'GEOSChem.StateMet.*'
    ds2 = xr.open_mfdataset(folder+FileStr)
    ds2 = ds2['Met_PMID']
    ds = xr.merge([ds2, ds])

    s_d = {}
    for site in GAW_sites:
        # Get: lat, lon, alt (press), timezone (UTC)
        lat, lon, alt, TZ = gaw_2_loc(site)
        # Use nearest grid box
#                lev2use = AC.find_nearest_value(HPa_r, alt)

        # Get the closest location
        ds_tmp = ds.sel(lat=lat, lon=lon, method='nearest')
        # then select height
#                lev=ds.lev.values[lev2use]
        lev2use = ds_tmp['Met_PMID'].mean(dim='time')
        first_level = int(lev2use.values[0])
        lev2use = AC.find_nearest(lev2use.values, alt)
#                int(lev2use.sel(lev=alt, method='nearest').values)
        ds_tmp = ds_tmp.isel(lev=lev2use)
        print(site, alt, lev2use, lev2use == 0)
        # extract data
        prefix = 'Jval_'
        vars2use = [i for i in ds_tmp.data_vars if prefix in i]
        ds_tmp = ds_tmp[vars2use]
        name_dict = zip(vars2use, [i.split(prefix)[-1] for i in vars2use])
        name_dict = dict(name_dict)
        ds_tmp = ds_tmp.rename(name_dict=name_dict)
        #
        for coord in [i for i in ds_tmp.coords if i != 'time']:
            del ds_tmp[coord]
        S = ds_tmp.to_dataframe()
        s_d[site] = S

    # - plot up
    specs2plot = ['NO2', 'HNO3', 'HNO2', 'BrO', 'IO', 'CH2I2', 'O3', 'PAN']
    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_test_Jvals'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    for site in GAW_sites:
        print(site)
        #
        nc_data = s_d[site]
        pf_data = df.loc[df['TYPE'] == site, :]
        for spec in specs2plot:
            print(spec)

            # Plot up PF
            pf_spec = None
            if spec == 'BrO':
                pf_spec = 'JVL_028'
            elif spec == 'O3':
                pf_spec = 'JVL_002'
            elif spec == 'NO2':
                pf_spec = 'JVL_011'
            elif spec == 'IO':
                pf_spec = 'JVL_116'
            elif spec == 'CH2I2':
                pf_spec = 'JVL_123'
            elif spec == 'PAN':
                pf_spec = 'JVL_059'
            elif spec == 'HNO3':
                pf_spec = 'JVL_016'
            elif spec == 'HNO2':
                pf_spec = 'JVL_015'
            elif spec == 'O1D':
                pf_spec = 'JVL_002' # Also ???
            else:
                print('case not setup for species: {}'.format(spec))
            if isinstance(pf_spec, str):
                print(spec, pf_spec)
                data = pf_data[pf_spec]
                plt.plot(data.index, data.values, label='PF')
            # Plot up NC
            data = nc_data[spec]
            plt.plot(data.index, data.values, label='NetCDF')
            plt.legend()
            plt.title('{} @ {}'.format(spec, site))
            # save
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()

    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')

    # - Plot up differences between J rates in nested and global run
    savetitle = 'ARNA_Jvals_Global_vs_nest_model'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    RunRoot = '/users/ts551/scratch/GC/rundirs/'
    RunStr = 'geosfp_4x5_standard.v12.9.0.BASE.2019.2020.'
    folderNest = RunRoot + RunStr + 'ARNA.Nest.repeat.JVALS/'
    folderGlobal = RunRoot + RunStr + 'ARNA.BCs.TEST.PF_Jrates.JVALS.GLOBAL/'
    files2useNest = list(sorted(glob.glob(folderNest + '/*plane*log*')))
    files2useGlobal = list(sorted(glob.glob(folderGlobal + '/*plane*log*')))
    for nfile2use, file2useNest in enumerate(files2useNest):
        date = file2use.split('plane.log.')[-1]
        print(file2use)
        file2useGlobal = files2useGlobal[nfile2use]
        file2use_dict = {'Nest': file2useNest, 'Global': file2useGlobal}
        dfs = {}
        for key in file2use_dict.keys():
            file2use = file2use_dict[key]
            vars, sites = AC.get_pf_headers(file2use, debug=debug)
            df, NIU = AC.pf_csv2pandas(file=file2use, vars=vars, epoch=True,
                                       r_vars=True)
            # Update the datetime index
            df = AC.DF_YYYYMMDD_HHMM_2_dt(df, rmvars=None, epoch=False)
            df.index.name = None
            dfs[key] = df
        # Now plot
        for nkey, key in enumerate(list(dfs.keys())):
            dfs[key]['JVL_134'].plot(label=key, color=['blue', 'red'][nkey])
            plt.title('JHNO3 for flighttrack on {}'.format(date))
        # save
        plt.ylabel('J-rate (s$^{-1}$)')
        plt.legend()
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()
    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def evaluate_regional_grid4GEOSChem(show_plot=False, dpi=320):
    """
    Evaluate the regional variable high(er) resolution (flex)grid for GEOS-chem
    """
    # - Load FAAM flighttrack data
    # Import all of the flights by campaign and plot on grid
    # All FAAM campaigns - included in James and Freya's analysis
    # NOTE: this will only use the downloaded files for now.
    campaign = 'ARNA-2'
    dfARNA = get_flighttracks4campaign(campaign)
    campaigns = [
        'ARNA-2', 'ACISIS-5', 'ACRUISE', 'ACISIS-4', 'MOYA-2', 'ACISIS-3',
        'ACISIS-2',
        #    'Clarify', # FAAM files not downloaded
        'MOYA-1',
        #    'ACISIS-1' # FAAM files not downloaded
    ]
    dfs = [get_flighttracks4campaign(i) for i in campaigns]

    # - Plot up high resolution modelling region around ARNA-2 flights
    savetitle = 'ARNA_high_resolution_model_grid'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    # Plot up a blank global background
    fig, ax = plt_highres_modelling_region(plot_blank_data=True)
    # Add the campaign flights to this
    campaign = 'ARNA-2'
    df = dfARNA
    LatVar = 'LAT_GIN'
    LonVar = 'LON_GIN'
    lats = df[LatVar].values
    lons = df[LonVar].values
    color_list = AC.get_CB_color_cycle()
    projection = ccrs.PlateCarree
    # Now scatter points on plot
    ax = add_scatter_points2cartopy_ax(ax=ax, lons=lons, lats=lats,
                                       color=color_list[0],
                                       label=campaign)

    fig.legend(loc=7)
    fig.suptitle('Flight-tracks during ARNA-2 campaign')
    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
    if show_plot:
        plt.show()
    plt.close()

    # - Plot up high resolution modelling region around BB FAAM flights
    fig, ax = plt_highres_modelling_region(plot_blank_data=True)
    for n_campaign, campaign in enumerate(campaigns):
        df = dfs[campaigns.index(campaign)]
        print(campaign)
        lats = df[LatVar].values
        lons = df[LonVar].values
        ax = add_scatter_points2cartopy_ax(ax=ax, lons=lons, lats=lats,
                                           color=color_list[n_campaign],
                                           label=campaign)
    title = "Flight-tracks during FAAM campaigns in biomass burning analysis"
    fig.suptitle(title)
    fig.legend(loc=7)
    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
    if show_plot:
        plt.show()
    plt.close()

    # - Now evaluate biomass burning emissions for grid
    earth0_folder = '/mnt/lustre/groups/chem-acm-2018/earth0_data/'
    HEMCO_folder = earth0_folder + 'GEOS//ExtData/HEMCO/'
    # Get GFAS emissions for the year

    # Plot by seasons - TODO

    # Plot for feb 2020
    folder = HEMCO_folder+'/GFAS/v2018-09/2020/'
    filename = 'GFAS_202002.nc'
    ds = xr.open_dataset(folder+filename)
    # Update lon to be in degrees West -
    var2use = 'cofire'
    ds = ds[[var2use]]
    ds = ds.assign_coords({'lon': ds.lon.values - 180})
    # Update name and scaling
    Uvar2use = '{} (1E-9 {})'.format(ds[var2use].long_name, ds[var2use].units)
    ds = ds.rename({var2use: Uvar2use})
    var2use = Uvar2use
    ds = ds[[var2use]].mean(dim='time') * 1E9
    # Remove zero data
    arr = ds[var2use].values
    arr[arr <= 0] = np.NaN
    ds[var2use].values = arr
    # And Roll the variables too
#    ds = ds.roll(lon=-int(len(ds.lon)/2))
    arr = ds[var2use].values
    arr = np.roll(arr, -int(len(ds.lon)/2), axis=1)
    ds[var2use].values = arr
    # Plot up the data
    fig, ax = plt_highres_modelling_region(ds=ds, var2use=var2use,
                                           plot_blank_data=False,
                                           rm_colourbar=False)

    fig.suptitle('Biomass burning emissions (GFAS) - Feb 2020 (ARNA-2)')
    del ds
    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
    if show_plot:
        plt.show()
    plt.close()

    # - Now evaluate dust emission for grid
    # Plot by seasons - TODO

    # Plot for February
    folder = HEMCO_folder+'/OFFLINE_DUST/v2019-01/0.5x0.625/2019/02/'
    files2use = glob.glob(folder+'*nc')
    ds = xr.open_mfdataset(files2use)
    # Combine all dust emissions
    var2use = 'Total dust emission (kg/m2/s)'
    ds[var2use] = ds['EMIS_DST1'].copy()
    ds[var2use] = ds[var2use].values + ds['EMIS_DST2']
    ds[var2use] = ds[var2use].values + ds['EMIS_DST3']
    ds[var2use] = ds[var2use].values + ds['EMIS_DST4']
    ds = ds[[var2use]].mean(dim='time')
    # Remove zero data
    arr = ds[var2use].values
    arr[arr <= 0] = np.NaN
    ds[var2use].values = arr
    # Plot up the data
    fig, ax = plt_highres_modelling_region(ds=ds, var2use=var2use,
                                           plot_blank_data=False,
                                           rm_colourbar=False)

    fig.suptitle('Dust emissions (online) - Feb *2019* (ARNA-2)')
    # Save to PDF
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
    if show_plot:
        plt.show()
    plt.close()

    # - Now evaluate NOx emission for grid?

    # - Others variables to plot / consider?
    # Night lights?

    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


if __name__ == "__main__":
    main()

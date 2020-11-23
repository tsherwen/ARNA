#!/usr/bin/python
"""
Driver for generic plotting functions following the ARNA campaign
"""
import arna as ar

def main():
    """
    Main driver function
    """
    # - run the extraction of ARNA flights from GEOS-CF data
#    context = 'talk'
    context = 'paper'
    # Plot up core FAAM data
#    extract_GEOS54all_ARNA_flights()
    ar.plt_timeseries_comp4ARNA_flights(inc_GEOSChem=False, context=context)
    ar.plt_comp_by_alt_4ARNA_flights(context=context)

    # plot up ToF-CIMS data
    ar.plt_comp_by_alt_4ARNA_flights_CIMS(context=context)
    ar.plt_timeseries_comp4ARNA_flights_CIMS(context=context)

    # Plot up nitrate aerosol data
    ar.plt_timeseries_comp4ARNA_flights_filters(context=context)

    #  Plot up PCASP/CDP date
    # NOTE: CAS data being ignored currently due to issue with mirror window
#    ar.plt_timeseries_comp4ARNA_flights_PCASP()

    # Plot up SWAS data
    ar.plt_timeseries_comp4ARNA_flights_SWAS(context=context)

    # Plot up vertical velocity and Roll, amongst other core physical vars
    ar.plt_timeseries_comp4ARNA_flights_PHYSICAL_VARS(context=context)

    # Plot up the temperature data from Hannah Price
    # N/A? this is only for 2019. Data to be worked up for 2020.

    # Plot a comparison
    ar.plt_timeseries_comp4ARNA_flights_NOy_ALL(context=context)

    # Plot up data for SLRs with and without dust
    ar.plt_comp_by_alt_4ARNA_all(just_SLR=False, context=context)
    ar.plt_comp_by_alt_4ARNA_all(just_SLR=True, context=context)
    ar.plt_comp_by_alt_4ARNA_all_DUST(plt_model=False, context=context)
    ar.plt_comp_by_alt_4ARNA_all_DUST(plt_model=True, context=context)
    ar.plt_comp_by_alt_4ARNA_CIMS_all_DUST(context=context)

    # Evaluate the high resolution modelling region
    ar.evaluate_regional_grid4GEOSChem()

    # Also plot up for related biomass-burning flights in MOYA campaign
    ar.plt_timeseries_comp4MOYA_flights()
    ar.plt_timeseries_comp4MOYA_flights_PHYSICAL_VARS()


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
    projection=ccrs.PlateCarree
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
    for n_campaign, campaign in enumerate( campaigns ):
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
    ds = ds.assign_coords({'lon':ds.lon.values -180})
    # Update name and scaling
    Uvar2use = '{} (1E-9 {})'.format(ds[var2use].long_name, ds[var2use].units)
    ds = ds.rename({var2use:Uvar2use})
    var2use = Uvar2use
    ds = ds[[var2use]].mean(dim='time') *1E9
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
    fig, ax = plt_highres_modelling_region(ds=ds,var2use=var2use,
                                           plot_blank_data=False,
                                           rm_colourbar=False )

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
    ds[var2use] = ds[var2use].values +ds['EMIS_DST2']
    ds[var2use] = ds[var2use].values +ds['EMIS_DST3']
    ds[var2use] = ds[var2use].values +ds['EMIS_DST4']
    ds = ds[[var2use]].mean(dim='time')
    # Remove zero data
    arr = ds[var2use].values
    arr[arr <= 0] = np.NaN
    ds[var2use].values = arr
    # Plot up the data
    fig, ax = plt_highres_modelling_region(ds=ds,var2use=var2use,
                                           plot_blank_data=False,
                                           rm_colourbar=False )

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
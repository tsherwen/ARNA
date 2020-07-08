"""
Plotting functions for ARNA campaign work
"""

import os
import sys
import glob
import gc
import numpy as np
import pandas as pd
import xarray as xr
import xesmf as xe
import AC_tools as AC
from netCDF4 import Dataset
from datetime import datetime as datetime_
import datetime as datetime
import time
from time import gmtime, strftime
import matplotlib.pyplot as plt
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
import seaborn as sns
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
#from cartopy.mpl.ticker import LatitudeLocator, LongitudeLocator
from shapely.geometry.polygon import LinearRing
import matplotlib

# import local ARNA module functions
from . GEOS import *
from . core import *
from . utils import *


def set_limits4ar_plotted_range(var, verbose=False):
    """
    Retrieve plotting limits for a given species
    """
    # Default is that there are no plotting limits
    vmin = None
    vmax = None
    # Dust (GEOS-5 product)
    if 'Dust' in var:
        vmin = 5
        vmax = 135
    # All particulate matter
    elif 'PM' in var:
        vmin = 5
        vmax = 50
    # All particulate matter
    elif 'CO' in var:
        vmin = 80
        vmax = 150
    # All particulate matter
    elif 'O3' in var:
        vmin = 0
        vmax = 50
    # NO2
    elif 'NO2' in var:
        vmin = 0
        vmax = 30/1E3
    # NOy
    elif 'NOy' in var:
        vmin = 0
        vmax = 1
    else:
        if verbose:
            print('WARNING: No range set for {}'.format(var))
    return vmin, vmax


def get_species_units(input):
    """
    Retrieve units for GEOS-CF species and derived values
    """
    if 'NOy/' in input:
        units = 'unitless'
    elif 'pm25' in input:
        units = '$\mu$g m$^{-3}$'
    elif 'Dust' in input:
        units = '$\mu$g m$^{-3}$'
    elif 'NOy' in input:
        units = 'ppbv'
    else:
        units = 'v/v'
    return units


def get_vmin_value4var(input):
    """
    Set a vmin value for a variable (below which values are set to NaN)
    """
    #
    if 'Dust' in input:
        vmin = 15
    elif 'NOy' in input:
        vmin = 0.5
    else:
        print('WARNING: vmin case not set for variable! ({})'.format(input))
    return vmin


def get_cmap4var(input):
    """
    Get a colour-map for a specific variable.
    """
    #
    if 'Dust' in input:
        if input == 'Dust':
            nticks=9-1
#            nticks=10-1
#             nticks=11-1
#             nticks=12-1
#             ticks = np.linspace(vmin2, vmax2, nticks+1)
#             ticks = [15*(i+1) for i in range(9)]
            ticks = [15*(i) for i in range(11)]
        # for PM2.5 dust
        else:
            # Get min and max values to plot
            vmin2, vmax2 = set_limits4ar_plotted_range(input)
            vmin2 = get_vmin_value4var(input)
            # Now set the number of ticks and range
            nticks=8-1
            ticks = np.linspace(vmin2, vmax2, nticks+1)
        cmap = get_reduced_cmap(cmap='Reds', maxval=0.75)
#        cmap = get_reduced_cmap(cmap='Greens', maxval=0.75)
#        cmap = get_reduced_cmap(cmap='Greys', maxval=0.75)
        cmap = AC.mk_discrete_cmap( nticks=nticks, cmap=cmap)
    elif input == 'NOy':
        ticks = None
        nticks = 6-1
        cmap = AC.mk_discrete_cmap( nticks=nticks, cmap='Blues')
    else:
        print('WARNING: vmin case not set for variable! ({})'.format(input))
    return cmap, ticks, nticks


def plt_smaller_area_around_CVAO(ds, LonVar='lon',
                                 LatVar='lat',
                                 use_local_CVAO_area=False):
    """
    Plot up the area around Cape Verde
    """
    # Make a new variable using an old one as a template
    var2copy = [i for i in ds.data_vars][0]
    NewVar = 'Cape_Verde_LOCAL'
    dsCV = ds.copy()
    dsCV[NewVar] = ds[var2copy].copy()
    # Set all values to zero
    arr = dsCV[NewVar].values
    arr[:] = 1
    dsCV[NewVar].values = arr
    # Set extents - from main dataset
    x0 = float(ds[LonVar].min())
    x1 = float(ds[LonVar].max())
    y0 = float(ds[LatVar].min())
    y1 = float(ds[LatVar].max())
    extents = (x0, x1, y0, y1)
    # Sub-select data
    dsCV = select_smaller_area_around_Cape_Verde(dsCV)
    # Select surface and average over time
    dsCV = dsCV.sel(lev=dsCV.lev[0])
    dsCV = dsCV.mean(dim='time')
    # Now plot
    quick_map_plt_CV_1layer(dsCV, var2plot=NewVar, extents=extents,
                            use_local_CVAO_area=use_local_CVAO_area,
                            extra_str='ARNA', save_plot=True)


def plot_individual_spec_alt_slices(ds, folder='./',
                                    extr_title_str=None,
                                    vars2use=None):
    """
    Plot up other core species individually
    """
    # Use all variables in dataset if
    if isinstance(vars2use, type(None)):
        vars2use = list(ds.data_vars)
    # Loop by variables  and plot for each time and level
    for var2plot in vars2use:
        # Get the LateX for of the species name
        try:
            LaTeX_spec = AC.latex_spec_name(var2plot)
        except KeyError:
            print('WARNING: not converted {} to LaTeX form'.format(var2plot))
            LaTeX_spec = var2plot
        # Then loop by times in output
        dts = AC.dt64_2_dt( ds.time.values )
#        for n_time, t in enumerate( dts[:3] ): # unhash if testing
        for n_time, t in enumerate( dts ):
            # Sub select for time
            ds_tmp = ds[[var2plot]].sel(time=ds.time.values[n_time])
            # Create a date string to use to save file
            date_str = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
            date_str = date_str.format(t.year, t.month, t.day, t.hour,
                                       t.minute)
            title_date = t.strftime('%Y/%m/%d %H:%M')
            for lev2use in ds.lev.values:
#            for lev2use in ds.lev.values[:2]: # unhash if testing
                print( var2plot, lev2use, title_date)
                ds_tmpII = ds_tmp[[var2plot]].sel(lev=lev2use)
                # Setup a string for the title
                title = '[{}] @ {:.0f}hPa on {}'
                title = title.format(LaTeX_spec, lev2use, title_date)
                if not isinstance(extr_title_str, type(None)):
                    title += '\n '+ extr_title_str
                # Force use of standard name as long name
                attrs = ds_tmpII[var2plot].attrs
                attrs['long_name'] = LaTeX_spec
                ds_tmpII[var2plot].attrs = attrs
                # Set extra string for filename
                extra_str = 'ARNA_lev_{:.0f}_hPa_{}'.format(lev2use, date_str)
                # now plot
                quick_map_plt_CV_1layer(ds_tmpII, var2plot=var2plot,
                                        use_local_CVAO_area=True,
                                        extra_str=extra_str, extend='both',
                                        title=title, folder=folder,
                                        save_plot=True )
                # Do some clean up
                plt.close('all')
                del ds_tmpII
            del ds_tmp
            gc.collect()


def plot_up_latitudinally_sampled_locs(ds, var2use='noy', extents=None,
                                       add_detailed_map=True,
                                       add_flyable_range_as_circle=True,
                                       add_ARNA_locs=True,
                                       ):
    """
    plot up locations that data is sliced by lon
    """
    # Use example data
    if isinstance(ds, type(None)):
        folder = '/Users/tomassherwen/Google_Drive/Data/ARNA/GEOS_CF/'
        folder += '/data_GEOSCF_2019_12_14/'
        filename = 'ARNA_GEOSCF_chm_inst_1hr_g1440x721_p23_Cape_Verde_2019_353_noy_'
        filename += 'lvls_1000_900_800_700_600_500.nc'
        ds = xr.open_dataset( folder + filename )
    # Local area analysed as Cape Verde
    x0 = -30
    x1 =-10
    y0 = 0
    y1 = 25
    # - Select the data
    # Just get an example dataset
    ds = ds[[var2use]]
    # Select a single level and time
    ds = ds.sel(time=ds.time[0])
    ds = ds.sel(lev=ds.lev[0])
    # Set values region
    bool1 = ((ds.lon >= x0) & (ds.lon <= x1)).values
    bool2 = ((ds.lat >= y0) & (ds.lat <= y1)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    # Set all values to 1
    arr = ds[var2use].values
    arr[:] = np.NaN
    ds[var2use].values = arr
    # Plot the data
    projection = ccrs.PlateCarree()
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection=projection, aspect='auto', alpha=0.5)
    # Plot up the dummy data
    LatVar = 'lat'
    LonVar = 'lon'
    ds[var2use].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                             transform=ccrs.PlateCarree())
    # Now plot as a linear ring
    lons = (-15, -30, -30, -15)
    lats2plot = [12, 13, 14, 15, 16, 17]
    for lat in lats2plot:
        # Set lats
        lats = (lat, lat+0.125, lat+0.125, lat)
        #
        ring = LinearRing(list(zip(lons, lats)))
        ax.add_geometries([ring], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='green',
                          zorder=10, linestyle='-',
                          )
    # Add some grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=.5, color='gray', alpha=0.25, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    # Mark known places to help geo-locate viewers
    if add_ARNA_locs:
        locs2plot  = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
        for loc2plot in locs2plot:
            lon, lat, alt = AC.get_loc(loc2plot)
            # Now plot up locations
            ax.plot(lon, lat, 'bo', markersize=5, markerfacecolor='none',
                    markeredgewidth=2,
                    zorder=10,
                    markeredgecolor='black',
                    transform=ccrs.PlateCarree())
            # Add a label for the location?
#            ax.text(lon, lat+0.25, loc2plot, transform=ccrs.PlateCarree())
    # Add the flyable range of the FAAM BAE146
    if add_flyable_range_as_circle:
        # Approximate from James' max distance
        # ( 16.8331-13 ) *110667.45
        locs4circles = 'Dakar', 'Sao Vicente Airport',
        for loc in locs4circles:
            # Get locations to centre circle on
            lon, lat, alt = AC.get_loc(loc)
            # Radius in degrees
#            radius = 16.8331-13
            radius = 21 - 16.8331
            # Plot up circle
            ax.add_patch(mpatches.Circle(xy=[lon, lat],
                                         radius=radius,
                                         transform=projection,
                                         facecolor='none',
                                         edgecolor='grey',
                                         linestyle=':',
                                         zorder=10
                                         ))
    # Get limits of plotting data
    if isinstance(extents, type(None)):
        x0 = float(ds[LonVar].min())
        x1 = float(ds[LonVar].max())
        y0 = float(ds[LatVar].min())
        y1 = float(ds[LatVar].max())
        extents = (x0, x1, y0, y1)
    ax.set_extent(extents, crs=ccrs.PlateCarree())
    # Beautify the figure/plot
    if add_detailed_map:
        # Add borders to map
        ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                       facecolor='none', zorder=50)
        # Also add minor islands (inc. Cape Verde)
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor=None,
                                                facecolor='none')
        ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)
    # Save plot
    savename = 'spatial_plot_Cape_Verde_plotted_data_lons'
    savename = AC.rm_spaces_and_chars_from_str(savename)
    plt.savefig(savename+'.png', dpi=dpi)


def plot_up_longitudinally_sampled_locs(ds, var2use='noy', extents=None,
                                       add_detailed_map=True,
                                       add_flyable_range_as_circle=True,
                                       add_ARNA_locs=True,
                                       ):
    """
    plot up locations that data is sliced by lon
    """
    # Use example data
    if isinstance(ds, type(None)):
        folder = '/Users/tomassherwen/Google_Drive/Data/ARNA/GEOS_CF/'
        folder += '/data_GEOSCF_2019_12_14/'
        filename = 'ARNA_GEOSCF_chm_inst_1hr_g1440x721_p23_Cape_Verde_2019_353_noy_'
        filename += 'lvls_1000_900_800_700_600_500.nc'
        ds = xr.open_dataset( folder + filename )
    # Local area analysed as Cape Verde
    x0 = -30
    x1 =-10
    y0 = 0
    y1 = 25
    # - Select the data
    # Just get an example dataset
    ds = ds[[var2use]]
    # Select a single level and time
    ds = ds.sel(time=ds.time[0])
    ds = ds.sel(lev=ds.lev[0])
    # Set values region
    bool1 = ((ds.lon >= x0) & (ds.lon <= x1)).values
    bool2 = ((ds.lat >= y0) & (ds.lat <= y1)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    # Set all values to 1
    arr = ds[var2use].values
    arr[:] = np.NaN
    ds[var2use].values = arr
    # Plot the data
    projection = ccrs.PlateCarree()
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection=projection, aspect='auto', alpha=0.5)
    # Plot up the dummy data
    LatVar = 'lat'
    LonVar = 'lon'
    ds[var2use].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                            transform=ccrs.PlateCarree())
    # Now plot as a linear ring
    lats = (5, 35, 35, 5)
    lons2plot = [-18 ,-19.5, -21, -22.5, -24, -25.5]
    for lon in lons2plot:
        # Set lats
        lons = (lon, lon+0.125, lon+0.125, lon)
        #
        ring = LinearRing(list(zip(lons, lats)))
        ax.add_geometries([ring], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='green',
                          zorder=10, linestyle='-',
                          )
    # Add some grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=.5, color='gray', alpha=0.25, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    # Mark known places to help geo-locate viewers
    if add_ARNA_locs:
        locs2plot  = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
        for loc2plot in locs2plot:
            lon, lat, alt = AC.get_loc(loc2plot)
            # Now plot up locations
            ax.plot(lon, lat, 'bo', markersize=5, markerfacecolor='none',
                    markeredgewidth=2,
                    zorder=10,
                    markeredgecolor='black',
                    transform=ccrs.PlateCarree())
            # Add a label for the location?
#            ax.text(lon, lat+0.25, loc2plot, transform=ccrs.PlateCarree())
    # Add the flyable range of the FAAM BAE146
    if add_flyable_range_as_circle:
#        n_points = 1000
        # Approximate from James' max distance
        # ( 16.8331-13 ) *110667.45
        locs4circles = 'Dakar', 'Sao Vicente Airport',
        for loc in locs4circles:
            # Get locations to centre circle on
            lon, lat, alt = AC.get_loc(loc)
            # Radius in degrees
#            radius = 16.8331-13
            radius = 21 - 16.8331
            # Plot up circle
            ax.add_patch(mpatches.Circle(xy=[lon, lat],
                                         radius=radius,
                                         transform=projection,
                                         facecolor='none',
                                         edgecolor='grey',
                                         linestyle=':',
                                         zorder=10
                                         ))
    # Get limits of plotting data
    if isinstance(extents, type(None)):
        x0 = float(ds[LonVar].min())
        x1 = float(ds[LonVar].max())
        y0 = float(ds[LatVar].min())
        y1 = float(ds[LatVar].max())
        extents = (x0, x1, y0, y1)
    ax.set_extent(extents, crs=ccrs.PlateCarree())
    # Beautify the figure/plot
    add_detailed_map = True
    if add_detailed_map:
        # Add borders to map
        ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                       facecolor='none', zorder=50)
        # Also add minor islands (inc. Cape Verde)
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor=None,
                                                facecolor='none')
        ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)
    # Save the plot to disk
    savename = 'spatial_plot_Cape_Verde_plotted_data_lons'
    savename = AC.rm_spaces_and_chars_from_str(savename)
    plt.savefig(savename+'.png', dpi=dpi)


def plt_timeseries4ds(ds, region='Cape_Verde', extr_str='',
                      vars2use=None, year=2018, verbose=False,
                      show_plot=False, dpi=320):
    """
    Plot timeseries of data at different heights
    """
    # Which variables to include in analysis?
    if not isinstance(vars2use, list):
        vars2use = [i for i in ds.data_vars]
    # - Now plot up species as PDf based on level
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context("talk", font_scale=0.75)
    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_timeseries_{}_{}'.format(region, year)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    color = 'red'
    alpha = 0.3
    xtickrotation = 45
    # Now loop and plot
    for lev in ds.lev.values:
        # Get vars for levels
        #        suffix = '-{:.0f}hPa'.format(lev)
        #        vars2plot = [i for i in df.columns if i.endswith(suffix) ]
        for var in vars2use:
            #            var = df_var.split(suffix)[0]
            print(lev, var)
            # Get the LateX for of the species name
            try:
                LaTeX_spec = AC.latex_spec_name(var)
            except KeyError:
                print('WARNING: not converted {} to LaTeX form'.format(var))
                LaTeX_spec = var
            # Retrive units
            units = ds[var].units
            # Set title
            title_str = 'Timeseries of [{}] @ {:.0f}hPa in Feb {}'
            title = title_str.format(LaTeX_spec, lev, year)
            # Setup the plot
            fig = plt.figure(figsize=(10, 6))
            ax = fig.add_subplot(111)
            # Setup a dictionary of colours
#            color_dict = dict(zip(lvls, AC.get_CB_color_cycle() ) )
            # Select data for species and level
            da = ds.sel(lev=lev)[var]
            # For sorting by time
            da = da.sortby('time')
            # Get mean, min and max
            mean_vals = da.mean(dim=['lat', 'lon'], skipna=True).values
            min_vals = da.min(dim=['lat', 'lon'], skipna=True).values
            max_vals = da.max(dim=['lat', 'lon'], skipna=True).values
            time = AC.dt64_2_dt(da.time.values)
            # Now plot up as a time series
            if mean_vals[np.isfinite(mean_vals)].shape[0] == 0:
                prtstr = 'WARNING: Not plotting {} @ {:.0f}hPa because all NaNs!'
                print(prtstr.format(var, lev))
                pass
            else:
                ax.plot(time, mean_vals, color=color, label=var,
                        #                linewidth=lw, ls=ls
                        )
                ax.fill_between(time, mean_vals, min_vals,
                                alpha=alpha, facecolor=color)
                ax.fill_between(time, mean_vals, max_vals,
                                alpha=alpha, facecolor=color)
                # Add a title
                plt.title(title)
                # Add a legend
                plt.legend()
                # Rotate dates by 45 degrees
#                xticks = ax.get_xticks()
#                labels = ax.get_xticklabels()
#                ax.set_xticks(xticks)
#                ax.set_xticklabels(labels, rotation=xtickrotation)
#                ax.set_xticklabels(labels)
#                print(labels)
#                ax.set_xticklabels(xticks, rotation=xtickrotation)
                # Make sure that all titles and labels are readable
                try:
                    plt.tight_layout()
                except:
                    pass
                # Save to PDF
                AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
                if show_plot:
                    plt.show()
            plt.close()
            del da
    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def PDF_on_species_in_ds4lvls(ds, region='Cape_Verde', extr_str='',
                              vars2use=None, year=2018, verbose=False,
                              show_plot=False, dpi=320):
    """
    plot the stats on species at different heights as a PDF
    """
    # Which variables to include in analysis?
    if not isinstance(vars2use, list):
        vars2use = [i for i in ds.data_vars]
    # Setup a dataframe to store data in
    df = pd.DataFrame()
    # Loop by variable
    vars2use = [i for i in ds.data_vars]
    for var in vars2use:
        # Loop by level
        lvls = np.array(ds.lev.values)
        for lvl in lvls:
            print(var, lvl)
            # Setup a name to save data to
            varname = '{}-{:.0f}hPa'.format(var, lvl)
            # Flatten the data
            vals = ds[var].sel(lev=lvl).values.flatten()
            S = pd.Series(vals)
            # Save to the main DataFrame
            df[varname] = S
            del vals, S

    # - Now plot up species as PDf based on level
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context("talk")
    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_PDF_of_concs_by_level_{}_{}'.format(region, year)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    # Now loop and plot
    for var in vars2use:
        #    for var in vars2use[:2]:
        # Get vars for levels
        vars2plot = [i for i in df.columns if i.startswith(var+'-')]
        # Get the LateX for of the species name
        try:
            LaTeX_spec = AC.latex_spec_name(var)
        except KeyError:
            print('WARNING: not converted {} to LaTeX form'.format(var))
            LaTeX_spec = var
        # Retrive units
        units = ds[var].units
        # Set title
        title = 'PDF of [{}] for hPa levels in Feb {}'.format(LaTeX_spec, year)
        # Setup the plot
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        # Setup a dictionary of colours
        lvls = np.array(ds.lev.values)
        color_dict = dict(zip(lvls, AC.get_CB_color_cycle()))
        # Plot up by level
        for lev2use in lvls:
            levstr = '{:.0f}hPa'.format(lev2use)
            var2plot = [i for i in vars2plot if levstr in i]
            ass_str = 'ERROR: There must only be one variable per level!'
            assert len(var2plot) == 1, ass_str
            if verbose:
                print(var, lev2use, var2plot)
            # Retrieve the data
            arr = df[var2plot].values
            # only consider values that are not NaN
            arr = arr[~np.isnan(arr)]
            # Plot up and add title
            axlabel = '{} ({})'.format(LaTeX_spec, units)
            ax = sns.distplot(arr, axlabel=axlabel, label=levstr,
                              color=color_dict[lev2use], ax=ax)
            # Make sure the scaling is correct
            ax.autoscale()
        # Add a title
        plt.title(title)
        # Add a legend
        plt.legend()
        # Make sure that all titles and labels are readable
        plt.tight_layout()
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()
    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plt_avg_spatial_by_lvl(ds, year=2018, vars2use=None,
                           verbose=False,
                           use_local_CVAO_area=False,
                           show_plot=False, dpi=320):
    """
    Make a PDF of average spatial concentrations by level
    """
    # Setup PDF to save plots to.
    savetitle = 'ARNA_avg_spatial_concs_{}'.format(year)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    # Which variables to plot
    if not isinstance(vars2use, list):
        vars2use = [i for i in ds.data_vars]
    vars2use = ['NOy', 'PM2.5(dust)', ]
    # Now loop and plot
    for var in vars2use:
        # Plot up by level
        for lev2use in ds.lev.values:
            if verbose:
                print(var, lev2use)
            # Get units for species
            units = ds[var].units
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lev=lev2use).mean(dim='time')
            # Get the LateX for of the species name
            try:
                LaTeX_spec = AC.latex_spec_name(var)
            except KeyError:
                print('WARNING: not converted {} to LaTeX form'.format(var))
                LaTeX_spec = var
            # Set title
            title = 'Average [{}] @ {:.0f}hPa in Feb {}'.format(
                LaTeX_spec, lev2use, year)
            # Plot up and add title
            quick_map_plt_CV_1layer(ds_tmp, var2plot=var, title=title,
                                    use_local_CVAO_area=use_local_CVAO_area,
                                    save_plot=False, units=units)
            del ds_tmp
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()
    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plot_average_spatial_concs4lon(ds, year=2018, vars2use=None,
                                   use_local_CVAO_area=False,
                                   show_plot=False, dpi=320,
                                   verbose=False,):
    """
    Make a PDF of average spatial concentrations by level
    """
    # Setup PDF to save plots to.
    savetitle = 'ARNA_avg_spatial_concs_{}'.format(year)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
    # Which variables to plot
    if not isinstance(vars2use, list):
        vars2use = [i for i in ds.data_vars]
    vars2use = ['NOy', 'PM2.5(dust)', ]
    # Now loop and plot
    for var in vars2use:
#    for var in vars2use[:2]:  # For testing
        # Plot up by level
        #        for lev2use in ds.lev.values:
        for lon2use in list(ds.lon.values)  :
            if verbose:
                print(var, lev2use)
            # Get units for species
            units = ds[var].units
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lon=lon2use).mean(dim='time')
            # Get the LateX for of the species name
            try:
                LaTeX_spec = AC.latex_spec_name(var)
            except KeyError:
                print('WARNING: not converted {} to LaTeX form'.format(var))
                LaTeX_spec = var
            # Set title
            title = 'Average [{}] @ {:.0f}hPa in Feb {}'.format(
                LaTeX_spec, lev2use, year)
            # Plot up and add title
#            quick_map_plt_CV_1layer(ds_tmp, var2plot=var, title=title,
#                                             use_local_CVAO_area=use_local_CVAO_area,
#                                             save_plot=False, units=units)
            # vertical plot
            del ds_tmp
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()
    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def quick_map_plt_CV_1layer(ds, var2plot=None, extra_str='',
                            projection=ccrs.PlateCarree(),
                            save_plot=True, show_plot=False,
                            savename=None, units=None, title=None,
                            LatVar='lat', LonVar='lon', fig=None,
                            ax=None, extents=None, region='Cape_Verde',
                            use_local_CVAO_area=True,
                            add_flyable_range_as_circle=True,
                            add_flyable_range=False,
                            add_detailed_map=True, add_ARNA_locs=True,
                            extend='neither', folder='./', dpi=320):
    """
    Plot up a quick spatial plot of data using cartopy

    Parameters
    -------
    ds (xr.Dataset): dataset object holding data to plot
    var2plot (str): variable to plot within the dataset
    LatVar, LonVar (str): variables to use for latitude and longitude
    save_plot (bool): save the plot as a .png ?
    show_plot (bool): show the plot on screen
    dpi (int): resolution to use for saved image (dots per square inch)
    savename (str): name to use for png of saved .png
    extra_str (str): extra string to append to save .png
    projection (cartopy.crs obj.):  projection to use
    fig (figure instance): matplotlib figure instance
    ax (axis instance): axis object to use

    Returns
    -------
    (None)
    """
    # Use the 1st data variable if not variable given
    if isinstance(var2plot, type(None)):
        pstr = 'WARNING: No variable to plot(var2plot), trying 1st data_var'
        print(pstr)
        var2plot = list(ds.data_vars)[0]

    # Setup figure and axis and plot
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=(10, 6))
#        fig = plt.figure()
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111, projection=projection, aspect='auto')
    # Setup plotted range
    vmin, vmax = set_limits4ar_plotted_range(var2plot)
    # Now plot up
    ds[var2plot].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                             transform=ccrs.PlateCarree(),
                             vmin=vmin, vmax=vmax, extend=extend)
    # Add some grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=.5, color='gray', alpha=0.25, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    # Limit plot to Cape Verde region
    if use_local_CVAO_area:
        x0 = -30
        x1 =-10
        y0 = 0
        y1 = 25
        extents = (x0, x1, y0, y1)
    # Mark known places to help geo-locate viewers
    if add_ARNA_locs:
#        colours = AC.get_CB_color_cycle()
        locs2plot  = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
        for loc2plot in locs2plot:
            lon, lat, alt = AC.get_loc(loc2plot)
            # Now plot up locations
            ax.plot(lon, lat, 'bo', markersize=5, markerfacecolor='none',
                    markeredgewidth=2,
                    zorder=10,
                    markeredgecolor='black',
                    transform=ccrs.PlateCarree())
            # Add a label for the location?
#            ax.text(lon, lat+0.25, loc2plot, transform=ccrs.PlateCarree())
    # Add a box to show the flyable range
    if add_flyable_range:
        # Get the minimum
        d = get_max_flying_range4BAE146()
        min_lon = d['min_lon']
        max_lon = d['max_lon']
        min_lat = d['min_lat']
        max_lat = d['max_lat']
        # Create new lists
        lons = [min_lon, min_lon, max_lon, max_lon]
        lats = [min_lat, max_lat, max_lat, min_lat]
        # Now plot as a linear ring
        ring = LinearRing(list(zip(lons, lats)))
        ax.add_geometries([ring], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='grey',
                          zorder=10, linestyle=':',
                          )
    if add_flyable_range_as_circle:
#        n_points = 1000
        # Approximate from James' max distance
        # ( 16.8331-13 ) *110667.45
        locs4circles = 'Dakar', 'Sao Vicente Airport',
        for loc in locs4circles:
            # Get locations to centre circle on
            lon, lat, alt = AC.get_loc(loc)
            # Radius in degrees
#            radius = 16.8331-13
            radius = 21 - 16.8331
            # Plot up circle
            ax.add_patch(mpatches.Circle(xy=[lon, lat],
                                         radius=radius,
                                         transform=projection,
                                         facecolor='none',
                                         edgecolor='grey',
                                         linestyle=':',
                                         zorder=10
                                         ))
    # Get limits of plotting data
    if isinstance(extents, type(None)):
        x0 = float(ds[LonVar].min())
        x1 = float(ds[LonVar].max())
        y0 = float(ds[LatVar].min())
        y1 = float(ds[LatVar].max())
        extents = (x0, x1, y0, y1)
    ax.set_extent(extents, crs=ccrs.PlateCarree())
    # Beautify the figure/plot
    if add_detailed_map:
        # Add borders for countries
        ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                       facecolor='none', zorder=50)
        # Also add minor islands (inc. Cape Verde)
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor=None,
                                                facecolor='none')
        ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)
    # Update the colour bar lanel
    if not isinstance(units, type(None)):
        im = ax.images
        cb = im[-1].colorbar
        cb.set_label('{}'.format(units))
    # Add a generic title if one is not provided
    if not isinstance(title, type(None)):
        plt.title(title)
    # Save the plot?
    if save_plot:
        if isinstance(savename, type(None)):
            savename = 'spatial_plot_{}_{}'.format(var2plot, extra_str)
        savename = AC.rm_spaces_and_chars_from_str(savename)
        plt.savefig(folder+savename+'.png', dpi=dpi)
    if show_plot:
        plt.show()



def plot_spatial_concs_2layer(ds, show_plot=False, folder=None,
                              var2plot1='NOy', var2plot2='PM2.5(dust)',
                              extr_title_str='', region='Cape_Verde',
                              add_max_vals_as_txt=False,
                              testing_mode=False,
                              verbose=False, testing=False ):
    """
    Plot up a two layer plot on a single map for given levels
    """
    # Local variables
    try:
        LaTeX_spec1 = AC.latex_spec_name(var2plot1)
    except KeyError:
        LaTeX_spec1 = var2plot1
    try:
        LaTeX_spec2 = AC.latex_spec_name(var2plot2)
    except KeyError:
        LaTeX_spec2 = var2plot2
    # Set data below threshold to zero based on variable name
    ds = set_values_below_range2NaNs4spec(var=var2plot1, ds=ds)
    ds = set_values_below_range2NaNs4spec(var=var2plot2, ds=ds)
    # Set lists of variables to loop and plot
    if testing_mode:
        levs2use = [700]
        times2use = ds.time[:4]
    else:
        levs2use = ds.lev.values
        times2use =  ds.time.values
    # Plot up by level and time
    for time2use in times2use:
        for lev2use in levs2use:
            # Get time as human readable string
            dstr = AC.dt64_2_dt([time2use])[0].strftime('%Y/%m/%d %H:%M')
            # Print out status
            pstr = "'plotting 2layer @ {:.0f}hPa on {}"
            print(pstr.format(lev2use, dstr) )
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lev=lev2use, time=time2use)
            # Set title
            title = '[{}] & [{}] @ {:.0f}hPa on {}'.format(
                LaTeX_spec1, LaTeX_spec2, lev2use, dstr)
            # Add extra string to existing title string
            title += '\n '+ extr_title_str
            # Save plots
            extra_str = 'lev_{}_dt_{}'.format( lev2use, dstr )
            quick_map_plt_2layer(ds_tmp, var2plot1=var2plot1,
                                  folder=folder, region=region,
                                  var2plot2=var2plot2, title=title,
                                  add_max_vals_as_txt=add_max_vals_as_txt,
                                  save_plot=True, extra_str=extra_str)
            # Tidy up...
            plt.close('all')


def set_values_below_range2NaNs4spec(var=None, ds=None):
    """
    To improve aesthetics of plots, values below a certain threshold are removed
    """
    # Limit plotted NOy values to those above 0.5 pptv
    if var == 'NOy':
        arr = ds[var].values
        arr[arr<0.5] = np.NaN
        ds[var].values = arr
    # Limit Dust values to those about
    elif 'Dust' in var:
        arr = ds[var].values
        arr[arr<15] = np.NaN
        ds[var].values = arr
    else:
        pstr = "WARNING: No case set for variable '{}', so not restricting array values"
        print(pstr.format(var))
    return ds


def plt_spatial_2layer_vertical_lon(ds, show_plot=False, folder=None,
                                    var2plot1='NOy',
                                    var2plot2='PM2.5(dust)',
                                    extr_title_str=None,
                                    testing_mode=False,
                                    ):
    """
    Plot up a two layer plot on a single map for given levels
    """
    # Local variables
    try:
        LaTeX_spec1 = AC.latex_spec_name(var2plot1)
    except KeyError:
        LaTeX_spec1 = var2plot1
    try:
        LaTeX_spec2 = AC.latex_spec_name(var2plot2)
    except KeyError:
        LaTeX_spec2 = var2plot2
    # Set data below threshold to zero based on variable name
    ds = set_values_below_range2NaNs4spec(var=var2plot1, ds=ds)
    ds = set_values_below_range2NaNs4spec(var=var2plot2, ds=ds)
    # Set lists of variables to loop and plot
    if testing_mode:
        lons2use = [-24.]
        times2use = ds.time[:4]
    else:
        lons2use = ds.lon.values
        times2use = ds.time.values
    # Plot by time and lon
    for time2use in times2use:
        for lon2use in lons2use:
            # Get time as human readable string
            dstr = AC.dt64_2_dt([time2use])[0].strftime('%Y/%m/%d %H:%M')
            # Print out status
            pstr = "'plotting 2layer @ {:.1f}$^{}$W on {}"
            print( pstr.format( lon2use*-1, '{\circ}', dstr) )
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lon=lon2use, time=time2use)
            # Set title
            title_str = '[{}] & [{}] @ {:.1f}$^{}$W on {}'
            title = title_str.format( LaTeX_spec1, LaTeX_spec2, lon2use*-1,
                                    '{\circ}', dstr)
            if not isinstance(extr_title_str, type(None)):
                title += extr_title_str
            # Save plots
            extra_str = 'lon_{}E_dt_{}'.format( lon2use, dstr )
            # Update the long_names - var1
            attrs = ds_tmp[var2plot1].attrs
            attrs['long_name'] = LaTeX_spec1
            ds_tmp[var2plot1].attrs = attrs
            # Update the long_names - var2
            attrs = ds_tmp[var2plot2].attrs
            attrs['long_name'] = LaTeX_spec2
            ds_tmp[var2plot2].attrs = attrs
            # Now call plotter
            quick_lon_plot_2layer(ds_tmp, var2plot1=var2plot1,
                                  folder=folder,
                                  var2plot2=var2plot2, title=title,
                                  save_plot=True, extra_str=extra_str
                                  )
            plt.close('all')


def plt_spatial_2layer_vertical_lat(ds, show_plot=False, folder=None,
                                    var2plot1='NOy',
                                    var2plot2='PM2.5(dust)',
                                    extr_title_str=None,
                                    testing_mode=False,
                                    ):
    """
    Plot up a two layer plot on a single map for given levels
    """
    # Local variables
    try:
        LaTeX_spec1 = AC.latex_spec_name(var2plot1)
    except KeyError:
        LaTeX_spec1 = var2plot1
    try:
        LaTeX_spec2 = AC.latex_spec_name(var2plot2)
    except KeyError:
        LaTeX_spec2 = var2plot2
    # Set data below threshold to zero based on variable name
    ds = set_values_below_range2NaNs4spec(var=var2plot1, ds=ds)
    ds = set_values_below_range2NaNs4spec(var=var2plot2, ds=ds)
    # Set lists of variables to loop and plot
    if testing_mode:
        lats2use = [16.]
        times2use = ds.time[:4]
    else:
        lats2use = ds.lat.values
        times2use = ds.time.values
    # Plot up by level and time
    for time2use in times2use:
        for lat2use in lats2use:
            # Get time as human readable string
            dstr = AC.dt64_2_dt([time2use])[0].strftime('%Y/%m/%d %H:%M')
            # Print out status
            pstr = "'plotting 2layer @ {:.1f}$^{}$N on {}"
            print( pstr.format( lat2use, '{\circ}', dstr) )
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lat=lat2use, time=time2use)
            # Set title
            title_str = '[{}] & [{}] @ {:.1f}$^{}$N on {}'
            title = title_str.format( LaTeX_spec1, LaTeX_spec2, lat2use,
                                    '{\circ}', dstr)
            if not isinstance(extr_title_str, type(None)):
                title += extr_title_str
            # Save plots
            extra_str = 'lat_{}N_dt_{}'.format( lat2use, dstr )
            # Update the long_names - var1
            attrs = ds_tmp[var2plot1].attrs
            attrs['long_name'] = LaTeX_spec1
            ds_tmp[var2plot1].attrs = attrs
            # Update the long_names - var2
            attrs = ds_tmp[var2plot2].attrs
            attrs['long_name'] = LaTeX_spec2
            ds_tmp[var2plot2].attrs = attrs
            # Now call plotter
            quick_lat_plt_2layer(ds_tmp, var2plot1=var2plot1,
                                  folder=folder,
                                  var2plot2=var2plot2, title=title,
                                  save_plot=True, extra_str=extra_str
                                  )
            plt.close('all')



def plot_CVAO_region_on_global_map(ds, var2use='NOy'):
    """
    Plot a global map to show ARNA campaign region
    """
    # - extents to use
    # extracted data from OPeNDAP
#     x0 = -30
#     x1 =-10
#     y0 = 0
#     y1 = 25
    # Local area analysed as Cape Verde
    x0 = -30
    x1 =-10
    y0 = 0
    y1 = 25
    # - Select the data
    # Just get an example dataset
    ds = ds[[var2use]]
    # Select a single level and time
    ds = ds.sel(time=ds.time[0])
    ds = ds.sel(lev=ds.lev[0])
    # Set values region
    bool1 = ((ds.lon >= x0) & (ds.lon <= x1)).values
    bool2 = ((ds.lat >= y0) & (ds.lat <= y1)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    # Set all values to 1
    arr = ds[var2use].values
    arr[:] = 1
    ds[var2use].values = arr
    # Plot the data
    projection = ccrs.Robinson()
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection=projection, aspect='auto', alpha=0.5)
    LatVar = 'lat'
    LonVar = 'lon'
    ds[var2use].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                             transform=ccrs.PlateCarree())
    # Beautify the figure/plot
    ax.coastlines()
    ax.set_global()
    # Force global perspective
    ax.set_global() # this will force a global perspective
    # Save
    savename = 'spatial_plot_Cape_Verde_flying_area'
    savename = AC.rm_spaces_and_chars_from_str(savename)
    plt.savefig(savename+'.png', dpi=dpi)


def quick_map_plt_2layer(ds, var2plot1=None, var2plot2=None, extra_str='',
                         projection=ccrs.PlateCarree(), folder=None,
                         save_plot=True, show_plot=False, savename=None,
                         units=None, title=None,
                         LatVar='lat', LonVar='lon', fig=None, ax=None,
                         extents=None,
                         add_ARNA_locs=True,
                         use_local_CVAO_area=True, region='Cape_Verde',
                         extend='both',
                         add_flyable_range_as_box=False,
                         add_flyable_range_as_circle=True,
                         add_detailed_map=True,
                         add_max_vals_as_txt=False,
                         dpi=320):
    """
    Plot up a quick spatial plot of data using cartopy

    Parameters
    -------
    ds (xr.Dataset): dataset object holding data to plot
    var2plot (str): variable to plot within the dataset
    LatVar, LonVar (str): variables to use for latitude and longitude
    save_plot (bool): save the plot as a .png ?
    show_plot (bool): show the plot on screen
    dpi (int): resolution to use for saved image (dots per square inch)
    savename (str): name to use for png of saved .png
    extra_str (str): extra string to append to save .png
    projection (cartopy.crs obj.):  projection to use
    fig (figure instance): matplotlib figure instance
    ax (axis instance): axis object to use
    add_ARNA_locs (bool):
    use_local_CVAO_area (bool):

    Returns
    -------
    (None)
    """
    # Use the 1st data variable if not variable given
    if isinstance(var2plot1, type(None)):
        pstr = 'WARNING: No variable to plot (var2plot), trying 1st data_var'
        print(pstr)
        var2plot1 = list(ds.data_vars)[0]
#        var2plot1 = 'NOy'
    if isinstance(var2plot2, type(None)):
        pstr = 'WARNING: No variable to plot (var2plot), trying 1st data_var'
        print(pstr)
        var2plot2 = list(ds.data_vars)[0]
#        var2plot2 = 'PM2.5(dust)'
    # Setup figure and axis and plot
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=(10, 6))
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111, projection=projection, aspect='auto')

    # - Plot first var
    alpha = 0.5
    # Setup plotted range
    vmin1, vmax1 = set_limits4ar_plotted_range(var2plot1)
    vmin1 = get_vmin_value4var(var2plot1)
    units1 = get_species_units(var2plot1)
    cmap, ticks, nticks = get_cmap4var(var2plot1)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap':cmap, 'extend':extend, }
    else:
        cbar_kwargs  = {'ticks':ticks, 'cmap':cmap, 'extend':extend,  }
    # Now plot up var1
    ds[var2plot1].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                              transform=ccrs.PlateCarree(),
                              vmin=vmin1, vmax=vmax1,
                              zorder=1, alpha=alpha,
                              cmap=cmap,
                              cbar_kwargs=cbar_kwargs,
#                             extend=extend,
                              )
    # Update the units on the colour bar panel
    im = ax.images
    cb = im[-1].colorbar
    try:
        LaTeX_spec1 = AC.latex_spec_name(var2plot1)
    except KeyError:
        LaTeX_spec1 = var2plot1
    cb.set_label('{} ({})'.format(LaTeX_spec1, units1))
    # - Plot second var
    alpha = 0.4
    # Now plot up var 2
    vmin2, vmax2 = set_limits4ar_plotted_range(var2plot2)
    vmin2 = get_vmin_value4var(var2plot2)
    units2 = get_species_units(var2plot2)
    cmap, ticks, nticks = get_cmap4var(var2plot2)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap': cmap, 'extend' : extend, }
    else:
        cbar_kwargs  = {'ticks': ticks, 'cmap': cmap, 'extend' : extend, }
    # Now plot up var2
    ds[var2plot2].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                              transform=ccrs.PlateCarree(),
                              vmin=vmin2, vmax=vmax2,
                              zorder=1, alpha=alpha, cmap=cmap,
                              cbar_kwargs=cbar_kwargs,
#                             extend=extend,
                              )
    # Update the units on the colour bar panel
    im = ax.images
    cb = im[-1].colorbar
    try:
        LaTeX_spec2 = AC.latex_spec_name(var2plot2)
    except KeyError:
        LaTeX_spec2 = var2plot2
    cb.set_label('{} ({})'.format(LaTeX_spec2, units2))
    # Update the xaxis ticks (lon) to deg W
#     update_lon_units = True
#     if update_lon_units:
#         xticks = ax.get_xticks()
#         xticklabels = ax.get_xticklabels()
#         ax.set_xticklabels([str(i)*-1 for i in xticks])

    # - Update plot aesthetics
    # Add some grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=.5, color='gray', alpha=0.25, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    # Just plot over the CVAO region?
    if use_local_CVAO_area and (region != 'Cape_Verde_Flying'):
        x0 = -30
        x1 =-10
        y0 = 0
        y1 = 25
        extents = (x0, x1, y0, y1)
    elif (region == 'Cape_Verde_Flying'):
        x0 = -29.1
        x1 =-15.9
        y0 = 11.9
        y1 = 21.1
        extents = (x0, x1, y0, y1)
    # Add extra lat and lon grid libnes
    if (region == 'Cape_Verde_Flying'):
        # Which X tickst to use?
        # x axis
        xticks = np.arange(ds.lon.values.min(), ds.lon.values.max(), 0.2   )
        # y axis
        yticks = np.arange(ds.lat.values.min(), ds.lat.values.max(), 0.2   )
        # --- New approach
#        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#        LON_FORMATTER = LONGITUDE_FORMATTER(dms=True)#, auto_hide=True)
#        LON_FORMATTER.set_locs(xticks)
#        LAT_FORMATTER = LATITUDE_FORMATTER(dms=True)#, auto_hide=True)
#        LAT_FORMATTER.set_locs(yticks)
#        gl.xformatter = LON_FORMATTER
#        gl.yformatter = LAT_FORMATTER
        # --- Old approach
#        xmajor = np.arange(ds.lon.values.min(), ds.lon.values.max(), 1   )
#        xminor = [i for i in xticks if i not in xmajor]
#        xminor_locator = mticker.FixedLocator(xminor)
#        xmajor_locator = mticker.FixedLocator(xmajor)
#        ax.xaxis.set_major_locator(xmajor_locator)
#        ax.xaxis.set_minor_locator(xminor_locator)
        gl.xlocator = mticker.FixedLocator(xticks) # last working setting...
#        ymajor = np.arange(ds.lat.values.min(), ds.lat.values.max(), 1 )
#        yminor = [i for i in yticks if i not in ymajor]
        gl.ylocator = mticker.FixedLocator(yticks) # last working setting...
#        get_labels
#        ymajor_locator = mticker.FixedLocator(ymajor)
#        yminor_locator = mticker.FixedLocator(yminor)
#        ax.yaxis.set_major_locator(ymajor_locator)
#        ax.yaxis.set_minor_locator(yminor_locator)
        # tight off the main labels.
#        gl.xlabels_bottom = False
#        gl.ylabels_left = False
        gl.xlabels_bottom = True # last working setting...
        gl.ylabels_left = True # last working setting...
        # Add main labels
#        for lon in xmajor:
#            buffer = 0.25
#            ax.text( y0, lon, '{}'.format(lon), fontsize=10, alpha=0.5,
#                    horizontalalignment='center' )
#        for lat in ymajor:
#            buffer = 0.25
#            ax.text( lat, x0, '{}'.format(lon), fontsize=10, alpha=0.5,
#                    horizontalalignment='center' )
        # Make axis label text smaller
        gl.xlabel_style = {'size': 6, 'rotation' :90 }
        gl.ylabel_style = {'size': 6, }

    # Mark a known place to help us geo-locate ourselves
    if add_ARNA_locs:
        colours = AC.get_CB_color_cycle()
        locs2plot  = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
        for loc2plot in locs2plot:
            lon, lat, alt = AC.get_loc(loc2plot)
            ax.plot(lon, lat, 'bo', markersize=5, markerfacecolor='none',
                    markeredgewidth=2,
#                    markeredgecolor=colours[0],
                    markeredgecolor='black',
                    transform=ccrs.PlateCarree())
#            ax.text(lon, lat+0.25, loc2plot, transform=ccrs.PlateCarree())

    # Add a box to show the flyable range
    if add_flyable_range_as_box:
        # Get the minimum
        d = get_max_flying_range4BAE146()
        min_lon = d['min_lon']
        max_lon = d['max_lon']
        min_lat = d['min_lat']
        max_lat = d['max_lat']
        # Create new lists
        lons = [min_lon, min_lon, max_lon, max_lon]
        lats = [min_lat, max_lat, max_lat, min_lat]
        # Now plot as a linear ring
        ring = LinearRing(list(zip(lons, lats)))
        ax.add_geometries([ring], ccrs.PlateCarree(),
                          facecolor='none', edgecolor='grey',
                          zorder=10, linestyle=':',
                          )
    if add_flyable_range_as_circle:
#        n_points = 1000
        # Approximate from James' max distance
        # ( 16.8331-13 ) *110667.45
        locs4circles = 'Dakar', 'Sao Vicente Airport',
        for loc in locs4circles:
            # Get locations to centre circle on
            lon, lat, alt = AC.get_loc(loc)
            # Radius in degrees
#            radius = 16.8331-13
            radius = 21 - 16.8331
            # Plot up circle
            ax.add_patch(mpatches.Circle(xy=[lon, lat],
                                         radius=radius,
#                                         color='red',
#                                         alpha=0.3,
                                         transform=projection,
                                         facecolor='none',
#                                         edgecolor='grey',
                                         edgecolor='black',
                                         linestyle=':',
                                         linewidth=3.0,
                                         zorder=100
                                         ))
    # Get limits of plotting data
    if isinstance(extents, type(None)):
        x0 = float(ds[LonVar].min())
        x1 = float(ds[LonVar].max())
        y0 = float(ds[LatVar].min())
        y1 = float(ds[LatVar].max())
        extents = (x0, x1, y0, y1)
    ax.set_extent(extents, crs=ccrs.PlateCarree())
    # Beautify the figure/plot
    if add_detailed_map:
        # Add borders for countries
        ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                       facecolor='none', zorder=50)
        # Also add minor islands (inc. Cape Verde)
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                edgecolor=None,
                                                facecolor='none')
        ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)

    # Update the xaxis ticks (lon) to deg W
#    update_lon_units = True
#    if update_lon_units:
#         xticks = ax.get_xticks()
#         xticklabels = ax.get_xticklabels()
#         ax.set_xticklabels([str(i)*-1 for i in xticks])

        # Trun of cartopy axis
#        gl.xlabels_bottom = False
        #
        # Add main labels
#         for lon in np.arange(-27, -12, 3):
#             buffer = 0.25
#             ax.text(x0-buffer, lon, '{}'.format(lon), fontsize=10, alpha=0.5,
#                     horizontalalignment='center', transform=ccrs.PlateCarree())


    # Add the grid box with maximum NOy*Dust (print value)
    if (region == 'Cape_Verde_Flying') and add_max_vals_as_txt:
        # Sub select data of interest
        try:
            ds_tmp = ds[['NOy', 'Dust', 'NOy*Dust']]
            # Get extents
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            # Reduce the dataset to this size.
            print(xlim, ylim)
            bool1 = ((ds_tmp.lon >= xlim[0]) & (ds_tmp.lon <= xlim[1])).values
            bool2 = ((ds_tmp.lat >= ylim[0]) & (ds_tmp.lat <= ylim[1])).values
            # Cut by lon, then lat
            ds_tmp = ds_tmp.isel(lon=bool1)
            ds_tmp = ds_tmp.isel(lat=bool2)
            # Find the
            da = ds_tmp[ 'NOy*Dust']
            print( da.where(da==da.max(), drop=True).squeeze() )
            lon = da.where(da==da.max(), drop=True).squeeze().lon.values
            lat = da.where(da==da.max(), drop=True).squeeze().lat.values
            # Add a cross on the map.
            radius = ( 21 - 16.8331 ) / 4
            # Plot up circle
            projection = ccrs.PlateCarree()
            ax.add_patch(mpatches.Circle(xy=[lon, lat],
                                         radius=radius,
                                         transform=projection,
                                         facecolor='none',
                                         edgecolor='black',
                                         linestyle='-',
                                         linewidth=3.0,
                                         zorder=100
                                         ))
            # Add a label saying the location and values for NOy and Dust
            # Add label for the airport
#            xtext = (xlim[1]-xlim[0])/2
#            ytext = (ylim[1]-ylim[0])/2
#            print(xtext, ytext)
#            xtext, ytext = -22.2, 20.
            xtext, ytext = -16.2, 20.
#            print(xtext, ytext)
            #
            ds_tmp = ds_tmp.sel(lon=lon)
            ds_tmp = ds_tmp.sel(lat=lat)
            NOy_at_loc = float(ds_tmp['NOy'].values)
            Dust_at_loc = float(ds_tmp['Dust'].values)
            lon_nav_units = convert_decimcal_degress2nav_format([lon*-1])[0]
            lat_nav_units = convert_decimcal_degress2nav_format([lat])[0]
#            plot_txt = 'NOy*Dust max. @ {}N {}E'.format(lat_nav_units, lon_nav_units)
            plot_txt = 'NOy*Dust max. @ {}N {}W'.format(lat_nav_units,
                                                        lon_nav_units)
            # Also add the values for NOy+Dust
            units1 = 'ppbv'
            units2 = '$\mu$g m$^{-3}$'
            pstr = '\n (NOy={:.2f} {}, Dust={:.1f} {})'
            plot_txt += pstr.format(NOy_at_loc, units1, Dust_at_loc, units2)

            # Now add to plot
            ax.text(xtext, ytext, plot_txt,
#                     '{}'.format(loc_),
                    fontsize=10,
#                    alpha=0.5,
#                    horizontalalignment='center',
                    horizontalalignment='right',
                    transform=projection,
                    zorder=100,
                    )
        except:
            print(ds.data_vars)
            pstr = 'WARNING: not adding waypoint as NOy*Dust var not present'
            print(pstr)
    # Add a generic title if one is not provided
    if not isinstance(title, type(None)):
        plt.title(title)
    # Save the plot?
    if save_plot:
        if isinstance(savename, type(None)):
            savename = 'spatial_plot_{}_{}_{}_{}'
            savename = savename.format(region, var2plot1, var2plot2, extra_str)
        savename = AC.rm_spaces_and_chars_from_str(savename)
        if isinstance(folder, type(None)):
            folder = './'
        plt.savefig(folder+savename+'.png', dpi=dpi)
    if show_plot:
        plt.show()



def quick_lon_plot_2layer(ds, var2plot1=None, var2plot2=None, extra_str='',
                          save_plot=True, show_plot=False, savename=None,
                          units=None, title=None, folder=None,
                          LatVar='lat', LonVar='lon', LevVar='lev',
                          fig=None, ax=None, extents=None,
                          add_ARNA_locs=True, use_local_CVAO_area=True,
                          extend='both', ylim=(0, 10), xlim=(5, 30),
                          dpi=320):
    """
    Plot up a quick longitude-altitude plot of data using cartopy

    Parameters
    -------
    ds (xr.Dataset): dataset object holding data to plot
    var2plot (str): variable to plot within the dataset
    LatVar, LonVar (str): variables to use for latitude and longitude
    save_plot (bool): save the plot as a .png ?
    show_plot (bool): show the plot on screen
    dpi (int): resolution to use for saved image (dots per square inch)
    savename (str): name to use for png of saved .png
    extra_str (str): extra string to append to save .png
    projection (cartopy.crs obj.):  projection to use
    fig (figure instance): matplotlib figure instance
    ax (axis instance): axis object to use
    add_ARNA_locs (bool):
    use_local_CVAO_area (bool):

    Returns
    -------
    (None)
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # Use the 1st data variable if not variable given
    if isinstance(var2plot1, type(None)):
        print('WARNING: No variable to plot (var2plot), trying 1st data_var')
#        var2plot1 = list(ds.data_vars)[0]
        var2plot1 = 'NOy'
    if isinstance(var2plot2, type(None)):
        print('WARNING: No variable to plot (var2plot), trying 1st data_var')
#        var2plot2 = list(ds.data_vars)[0]
        var2plot2 = 'PM2.5(dust)'
    # Setup figure and axis and plot
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=(12, 6))
    if isinstance(ax, type(None)):
        ax = fig.add_subplot(111)
    # Local variables (metres to kft)
    m2kft = 3.281
    # Convert lev to metres using hydrostatic aprox
    convert2km = True
    if convert2km:
        attrs = ds.lev.attrs
        print(attrs)
        if attrs['units'] == 'millibar':
            hPa_l = list(ds.lev.values).copy()
            ds.lev.values = AC.hPa_to_Km(hPa_l)
            attrs['units'] = 'km'
        else:
            pass
    # Then convert into kft?
    convert2kft = True
    if convert2kft:
        attrs = ds.lev.attrs
        print(attrs)
        if attrs['units'] == 'km':
#            hPa_l = list(ds.lev.values).copy()
            ds.lev.values = ds.lev.values * m2kft
            attrs['units'] = 'kft'
        else:
            pass
        #
        ylim = (ylim[0]*m2kft, ylim[1]*m2kft)

    # Setup devider axis for colourbar
    divider = make_axes_locatable(ax)
    # - Plot first var
    # Setup plotted range
    alpha = 0.5
    vmin1, vmax1 = set_limits4ar_plotted_range(var2plot1)
    vmin1 = get_vmin_value4var(var2plot1)
    units1 = get_species_units(var2plot1)
    cmap, ticks, nticks = get_cmap4var(var2plot1)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap': cmap, 'extend' : extend, 'pad':0.075, }
    else:
        cbar_kwargs  = {
        'ticks': ticks, 'cmap': cmap, 'extend' : extend, 'pad':0.075,
        }
    # Now plot up var1 - using pcolormesh
    cbar_ax = divider.append_axes("right", "2%", pad="1%")
    # Now plot
    ds[var2plot1].plot.pcolormesh(x=LatVar, y=LevVar, ax=ax,
                                 vmin=vmin1, vmax=vmax1,
                                 zorder=1, alpha=alpha,
#                                 origin='lower',
                                 yincrease=True,
                                 cmap=cmap,
#                                 extend=extend,
                                 cbar_ax=cbar_ax,
                                 cbar_kwargs=cbar_kwargs,
                                 )
    # Remove the title
    ax.set_title("")
    # - Plot second var
    # Setup a new axis for the colour bar
    cbar_ax = divider.append_axes("right", "2%", pad="7%")
    # Now plot up var 2
    alpha = 0.4
    vmin2, vmax2 = set_limits4ar_plotted_range(var2plot2)
    vmin2 = get_vmin_value4var(var2plot2)
    units2 = get_species_units(var2plot2)
    cmap, ticks, nticks = get_cmap4var(var2plot2)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap': cmap, 'extend' : 'both', }
    else:
        cbar_kwargs  = {'ticks': ticks, 'cmap': cmap, 'extend' : 'both', }
    # using pcolormesh
    ds[var2plot2].plot.pcolormesh(x=LatVar, y=LevVar, ax=ax,
                                 vmin=vmin2, vmax=vmax2,
                                 zorder=1, alpha=alpha,
                                 cmap=cmap,
#                                 origin='lower',
                                 yincrease=True,
                                 cbar_kwargs=cbar_kwargs,
                                 cbar_ax=cbar_ax,
                                 )
    # Remove the title
    ax.set_title("")
    # Add a box to show flyable area
    # 8km alt (and lon range)
    d = get_max_flying_range4BAE146()
    max_lat = d['max_lat']
    min_lat = d['min_lat']
    if convert2km and (not convert2kft):
        max_alt = d['max_alt'] /1E3
    elif convert2kft:
        max_alt = (d['max_alt'] / 1E3 * m2kft)
    else:
        pass
    max_alt_axis_coord = float(max_alt)/float(ylim[-1])
    print(max_alt_axis_coord)
    xrange = float(xlim[-1] - xlim[0])
    # Plot this up
    ax.axvline(x=max_lat, linestyle=':', alpha=0.5, color='grey', zorder=100,
               ymin=0, ymax=max_alt_axis_coord, linewidth=3.0)
    ax.axvline(x=min_lat, linestyle=':', alpha=0.5, color='grey', zorder=100,
               ymin=0, ymax=max_alt_axis_coord, linewidth=3.0)
    ax.axhline(y=max_alt, linestyle=':', alpha=0.5, color='grey', zorder=100,
               xmin=(min_lat - xlim[0] ) / xrange,
               xmax=(max_lat - xlim[0] ) / xrange,
               linewidth=3.0)

    # Add locations for airports
    locs2plot = [
    'Praia Airport', 'Dakar',  'Gran Canaria Airport', 'Sao Vicente Airport',
    'Lisbon Airport',  'Paris (Charles de Gaulle) Airport'
    ]
    for n, loc_ in enumerate( locs2plot ):
        lon_, lat_, NIU = AC.get_loc(loc_)
        # If lat in plotted range, then plot
        if (lat_ > xlim[0]) and (lat_ < xlim[-1]):
            ax.axvline(x=lat_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100, linewidth=3.0)
            # Add label for airports
            if n % 2 == 0:
                buffer = 0
            else:
                buffer = -0.5
            # Make a special case for Sao Vicente
            if loc_ == 'Sao Vicente Airport':
                buffer = -0.25
            # Set a base height for
            base = 9
            if convert2kft:
                base = base * m2kft
                buffer = buffer *3*1.5
            # Add label for the airport
            ax.text(lat_, base+buffer, '{}'.format(loc_), fontsize=10,
                    alpha=0.5,
                    horizontalalignment='center' )

    # Add lines for kft heights
    if convert2kft:
        hPa_heights = [1000, 900, 800, 700, 600, 500]
        km_heights = AC.hPa_to_Km(hPa_heights)
        kft_heights = [i*m2kft for i in km_heights]
        for n, height_ in enumerate( kft_heights ):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-2.5, height_, '{:.0f} hPa'.format(hPa_heights[n]),
                    fontsize=10,
                    alpha=0.5 )
    else:
        # Add lines for kft heights
        kft_heights = [20000, 15000, 10000, 5000]
        m_heights = [i/m2kft/1E3 for i in kft_heights]
        for n, height_ in enumerate( m_heights ):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-5, height_,
                    '{:.0f} kft'.format(kft_heights[n]/1E3),
                    fontsize=10, alpha=0.5 )

    # Beautify the figure/plot
    # Limit the yaxis
    ax.set_ylim(ylim)
    # Limit the xaxis
    ax.set_xlim(xlim)

    # Add hPA labels. - temporally remove .
#     press2plot = [1000, 900, 800, 700, 600, 500, 400]
#     ind = [ hPa_l.index(i) for i in press2plot ]
#     ax2.set_yticks(ds.lev[ind])
#     ax2.set_yticklabels([str(int(i)) for i in press2plot])
#     ax2.set_ylabel('Pressure [hPa]')

    # Add minor ticks for the x axis
#    ax.xaxis.grid(True, which='minor')

    # Add a generic title if one is not provided
    if not isinstance(title, type(None)):
#        ax.set_title(title)
        fig.suptitle(title)
    # Force tight layout
    plt.tight_layout(rect=(0, 0, 1, 0.95),)
#    plt.tight_layout()
    # Save the plot?
    if save_plot:
        if isinstance(savename, type(None)):
            sstr = 'spatial_plot_{}_{}_{}'
            savename = sstr.format(var2plot1, var2plot2, extra_str)
        savename = AC.rm_spaces_and_chars_from_str(savename)
        if isinstance(folder, type(None)):
            folder = './'
        plt.savefig(folder+savename+'.png', dpi=dpi)
    if show_plot:
        plt.show()


def plt_comp_by_alt_4ARNA_all(dpi=320, just_SLR=True, show_plot=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    import seaborn as sns
    # Which flights to plot?
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
    # Just use non-transit ARNA flights
    flights_nums = [
#    217,
    218, 219, 220, 221, 222, 223, 224, 225,
    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_FAAM_core4flightnum(flight_ID=flight_ID)
#        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_obs[flight_ID] = df
    # Observations
    dfs_mod = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_mod[flight_ID] = df
    # Combine to a single dataframe
    df_mod = pd.concat([dfs_mod[i] for i in dfs_mod.keys()], axis=0)
    df_obs = pd.concat([dfs_obs[i] for i in dfs_obs.keys()], axis=0)
    # Only consider data during SLRs?
    if just_SLR:
        df_obs = df_obs.loc[df_obs['IS_SLR'] == True, :]
        df_mod = df_mod.loc[df_mod['IS_SLR'] == True, :]
        extr_str = '_JUST_SLR'
    else:
        extr_str = ''
    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_altitude_binned_{}{}'.format('ALL', extr_str)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # - Plot up location of flights
    if just_SLR:
        title = "Flight tracks for 'Straight and Level' Runs during ARNA"
    else:
        title = 'Flight tracks for all flights during ARNA'
    plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID,
                                       title=title)
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
    plt.close()

    # - Put observations and vars to plot into a dictionary
    sns.set(color_codes=True)
    sns.set_context("paper", font_scale=0.75)
    # Force alt to be in units of km
    ALT_var = 'Altitude (km)'
    Y_unit = ALT_var
    df_mod[ALT_var] = AC.hPa_to_Km( df_mod['model-lev'].values )
    df_obs[ALT_var] = df_obs['ALT_GIN'].values / 1E3
    #
    data_d = {'GEOS-CF': df_mod, 'Obs.':df_obs}

    # - Now plot up flight time series plots by variable
    if just_SLR:
        title_str = "Altitude binned '{}' ({}) for all 'Straight and Level Runs'"
    else:
        title_str =  "Altitude binned '{}' ({}) for all flights"
    # Setup color dictinoary
    color_dict = {'GEOS-CF': 'red', 'Obs.':'k'}
    unit_d = {}
    mod2obs_varnames = {
    'CO':'CO_AERO', 'O3':'O3_TECO', 'NO2':'no2_mr', 'NO':'no_mr',
    'HNO2':'hono_mr',
    'NOx':'NOx'
    }
    units_d = {
    'CO':'ppbv', 'O3':'ppbv', 'NO2':'pptv', 'NO':'pptv', 'NOx':'pptv',
    'HNO2':'pptv', 'HONO':'pptv',
    }
    range_d = {
    'CO':(50, 400), 'O3':(-10, 100), 'NO2':(-50, 500), 'NO':(-50, 500),
    'NOx':(-50, 500),
    'HNO2':(-60, 60), 'HONO':(-60, 60),
    }
    NOx_specs = ['HNO2', 'NOx', 'NO', 'NO2', 'HONO']
    # - by variable
    runs = list(sorted(data_d.keys()))
    # Which variables to use?
    vars2plot = list(sorted(mod2obs_varnames.keys()))[::-1]
    vars2plot = ['CO', 'O3', 'NOx', 'NO2', 'NO', 'HNO2']
    print(vars2plot)
    print(df_obs.columns)
    vars2plot = [
    i for i in vars2plot if mod2obs_varnames[i] in df_obs.columns
    ]
    # What bins should be used?
    bins = [0.5*i for i in np.arange(15)]
    for var2plot in vars2plot:
        fig = plt.figure()
        ax = plt.gca()
        # Now loop data
        for n_key, key_ in enumerate(runs):
            print(n_key, key_, var2plot )
            #
            if key_ == 'Obs.':
                varname = mod2obs_varnames[var2plot]
            else:
                varname = var2plot
            # Setup an axis label
            units = units_d[var2plot]
            xlabel = '{} ({})'.format( var2plot, units )
            # Add alt to DataFrame
            df = pd.DataFrame({
            var2plot: data_d[key_][varname], ALT_var: data_d[key_][ALT_var]
            })
            #
            if key_ != 'Obs.':
                scaleby = AC.get_unit_scaling(units)
                df[var2plot] = df[var2plot].values * scaleby

            # drop any NaNs from the DataFrame
            s_shape = df.shape
            df.dropna(axis=0, how='any', inplace=True)
            if s_shape != df.shape:
                pcent = (float(df.shape[0]) - s_shape[0])/s_shape[0] * 100.
                pstr_dtr = 'WANRING dropped values - shape {}=>{} ({:.2f})'
                print(pstr_dtr.format(s_shape, df.shape, pcent))
            # Plot up as binned boxplots using existing function
            try:
                AC.binned_boxplots_by_altitude(df=df, fig=fig, ax=ax,
                                               var2bin_by=ALT_var,
                                               label=key_, xlabel=xlabel,
                                               binned_var=var2plot,
                                               num_of_datasets=len(runs),
                                               bins=bins,
                                               widths = 0.15,
                                               dataset_num=n_key,
                                               color=color_dict[key_])
            except:
                pass
            # Make NOx species be on a log scale
            if (var2plot in NOx_specs):
                ax.set_xscale('log')
            #                     ax.set_xlim( (1E-5, 1E3) )
                ax.set_xlim( (0.3, 400) )
            else:
                ax.set_xscale('linear')
            # Beautify plot
            plt.legend()
            plt.title(title_str.format(var2plot, units, flight_ID ))
            plt.xlim(range_d[var2plot])

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()

    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plt_comp_by_alt_4ARNA_all_DUST(dpi=320, just_SLR=True,
                                   plt_model=False, show_plot=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    import seaborn as sns
    # Which flights to plot?
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
    # Just use non-transit ARNA flights
    flights_nums = [
#    217,
    218, 219, 220, 221, 222, 223, 224, 225,
    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_FAAM_core4flightnum(flight_ID=flight_ID)
#        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_obs[flight_ID] = df
    # Observations
    dfs_mod = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_mod[flight_ID] = df
    # Combine to a single dataframe
    df_mod = pd.concat([dfs_mod[i] for i in dfs_mod.keys()], axis=0)
    df_obs = pd.concat([dfs_obs[i] for i in dfs_obs.keys()], axis=0)
    # Only consider data during SLRs?
    if just_SLR:
        df_obs = df_obs.loc[df_obs['IS_SLR'] == True, :]
        df_mod = df_mod.loc[df_mod['IS_SLR'] == True, :]
        extr_str = '_JUST_SLR_DUST'
    else:
        extr_str = '_DUST'
    # - Setup data objects or plotting
    # Force alt to be in units of km
    ALT_var = 'Altitude (km)'
    Y_unit = ALT_var
    df_mod[ALT_var] = AC.hPa_to_Km( df_mod['model-lev'].values )
    df_obs[ALT_var] = df_obs['ALT_GIN'].values / 1E3
    # Plot up just observations? Or model too?
    data_d = {
    'Obs.':  df_obs.loc[df_obs['IS_DUST'] == False, :],
    'Obs. (Dust)': df_obs.loc[df_obs['IS_DUST'] == True, :],
    }
    if plt_model:
        data_d['GEOS-CF'] = df_mod.loc[df_mod['IS_DUST'] == False, :]
        data_d['GEOS-CF (dust)'] = df_mod.loc[df_mod['IS_DUST'] == True, :]
        extr_str += '_inc_MODEL'

    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_altitude_binned_{}{}'.format('ALL', extr_str)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # - Plot up location of flights
    if just_SLR:
        title = "Flight tracks for 'Straight and Level Runs' during ARNA"
    else:
        title = 'Flight tracks for all flights during ARNA'
    plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID,
                                       title=title)
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
    plt.close()

    # - Put observations and vars to plot into a dictionary
    sns.set(color_codes=True)
    sns.set_context("paper", font_scale=0.75)

    # - Now plot up flight time series plots by variable
    if just_SLR:
        title_str = "Altitude binned '{}' ({}) for all 'Straight and Level Runs'"
    else:
        title_str =  "Altitude binned '{}' ({}) for all flights"
    # Setup color dictinoary
    color_dict = {'GEOS-CF': 'red', 'Obs.':'k'}
    colors2use = AC.get_CB_color_cycle()
    runs2color = [i for i in data_d.keys() if i not in color_dict.keys()]
    for n_run, run in enumerate( runs2color ):
        color_dict[run] = colors2use[n_run]
    unit_d = {}
    mod2obs_varnames = {
    'CO':'CO_AERO', 'O3':'O3_TECO', 'NO2':'no2_mr', 'NO':'no_mr',
    'HNO2':'hono_mr',
    'NOx':'NOx'
    }
    units_d = {
    'CO':'ppbv', 'O3':'ppbv', 'NO2':'pptv', 'NO':'pptv', 'NOx':'pptv',
    'HNO2':'pptv', 'HONO':'pptv',
    }
    range_d = {
    'CO':(50, 400), 'O3':(-10, 100), 'NO2':(-50, 500), 'NO':(-50, 500),
    'NOx':(-50, 500),
    'HNO2':(-60, 60), 'HONO':(-60, 60),
    }
    NOx_specs = ['HNO2', 'NOx', 'NO', 'NO2', 'HONO']
    # - by variable
    runs = list(sorted(data_d.keys()))
    # Which variables to use?
    vars2plot = list(sorted(mod2obs_varnames.keys()))[::-1]
    vars2plot = ['CO', 'O3', 'NOx', 'NO2', 'NO', 'HNO2']
    print(vars2plot)
    print(df_obs.columns)
    vars2plot = [
    i for i in vars2plot if mod2obs_varnames[i] in df_obs.columns
    ]
    # What bins should be used?
    bins = [0.5*i for i in np.arange(15)]
    for var2plot in vars2plot:
        fig = plt.figure()
        ax = plt.gca()
        # Now loop data
        for n_key, key_ in enumerate(runs):
            print(n_key, key_, var2plot )
            #
            if ('Obs.' in key_):
                varname = mod2obs_varnames[var2plot]
            else:
                varname = var2plot
            # Setup an axis label
            units = units_d[var2plot]
            xlabel = '{} ({})'.format( var2plot, units )
            # Add alt to DataFrame
            df = pd.DataFrame({
            var2plot: data_d[key_][varname], ALT_var: data_d[key_][ALT_var]
            })
            print(df.describe())
            # Scale the model values to pptv/ppbv
            if ('Obs.' not in key_):
                scaleby = AC.get_unit_scaling(units)
                df[var2plot] = df[var2plot].values * scaleby

            # drop any NaNs from the DataFrame
            s_shape = df.shape
            df.dropna(axis=0, how='any', inplace=True)
            if s_shape != df.shape:
                pcent = (float(df.shape[0]) - s_shape[0])/s_shape[0] * 100.
                pstr_dtr = 'WANRING dropped values - shape {}=>{} ({:.2f} %)'
                print(pstr_dtr.format(s_shape, df.shape, pcent))
            # Plot up as binned boxplots using existing function
            try:
                AC.binned_boxplots_by_altitude(df=df, fig=fig, ax=ax,
                                               var2bin_by=ALT_var,
                                               label=key_, xlabel=xlabel,
                                               binned_var=var2plot,
                                               num_of_datasets=len(runs),
                                               bins=bins,
                                               widths=0.15,
                                               dataset_num=n_key,
                                               color=color_dict[key_])
            except:
                pass
            # Make NOx species be on a log scale
            if (var2plot in NOx_specs):
                ax.set_xscale('log')
            #                     ax.set_xlim( (1E-5, 1E3) )
                ax.set_xlim( (0.3, 400) )
            else:
                ax.set_xscale('linear')
            # Beautify plot
            plt.legend()
            plt.title(title_str.format(var2plot, units, flight_ID ))
            plt.xlim(range_d[var2plot])

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()

    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plt_comp_by_alt_4ARNA_CIMS_all_DUST(dpi=320, just_SLR=True,
                                        plt_model=False, show_plot=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    import seaborn as sns
    # Which flights to plot?
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
    # Just use non-transit ARNA flights
    flights_nums = [
#    217, # Missing data for C217 (NOy)
    218, 219, 220,
#    221, # Missing data for C221 (NOy)
    222, 223,
#    224,  # Missing data for C221 (BrO... )
#    225,  # Missing data for C221 (BrO... )
    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_CIMS_data4flight(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_obs[flight_ID] = df
    # Observations
    dfs_mod = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_mod[flight_ID] = df
    # Combine to a single dataframe
    df_mod = pd.concat([dfs_mod[i] for i in dfs_mod.keys()], axis=0)
    df_obs = pd.concat([dfs_obs[i] for i in dfs_obs.keys()], axis=0)
    # Only consider data during SLRs?
    if just_SLR:
        df_obs = df_obs.loc[df_obs['IS_SLR'] == True, :]
        df_mod = df_mod.loc[df_mod['IS_SLR'] == True, :]
        extr_str = '_JUST_SLR_DUST'
    else:
        extr_str = '_DUST'
    # - Setup data objects or plotting
    # Force alt to be in units of km
    ALT_var = 'Altitude (km)'
    Y_unit = ALT_var
    df_mod[ALT_var] = AC.hPa_to_Km( df_mod['model-lev'].values )
    df_obs[ALT_var] = df_obs['ALT_GIN'].values / 1E3
    # Plot up just observations? Or model too?
    data_d = {
    'Obs.':  df_obs.loc[df_obs['IS_DUST'] == False, :],
    'Obs. (Dust)': df_obs.loc[df_obs['IS_DUST'] == True, :],
    }
    if plt_model:
        data_d['GEOS-CF'] = df_mod.loc[df_mod['IS_DUST'] == False, :]
        data_d['GEOS-CF (dust)'] = df_mod.loc[df_mod['IS_DUST'] == True, :]
        extr_str += '_inc_MODEL'

    # Setup PDF to save PDF plots to
    savetitle = 'ARNA_altitude_binned_{}_CIMS_{}'.format('ALL', extr_str)
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # - Plot up location of flights
    if just_SLR:
        title = "Flight tracks for 'Straight and Level Runs' during ARNA"
    else:
        title = 'Flight tracks for all flights during ARNA'
    plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID,
                                       title=title)
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
    plt.close()

    # - Put observations and vars to plot into a dictionary
    sns.set(color_codes=True)
    sns.set_context("paper", font_scale=0.75)

    # - Now plot up flight time series plots by variable
    if just_SLR:
        title_str = "Altitude binned '{}' ({}) for all 'Straight and Level Runs'"
    else:
        title_str =  "Altitude binned '{}' ({}) for all flights"
    # Setup color dictinoary
    color_dict = {'GEOS-CF': 'red', 'Obs.':'k'}
    colors2use = AC.get_CB_color_cycle()
    runs2color = [i for i in data_d.keys() if i not in color_dict.keys()]
    for n_run, run in enumerate( runs2color ):
        color_dict[run] = colors2use[n_run]
    NOx_specs = ['HNO2', 'NOx', 'NO', 'NO2', 'HONO']
    mod2obs_varnames = {
    'BrO':'BrO',
    'HNO3':'HNO3',
    'HNO2':'HONO',
#        'CO':'CO_AERO', 'O3':'O3_TECO', 'NO2':'no2_mr', 'NO':'no_mr',
#        'HNO2':'hono_mr',
#        'NOx':'NOx'
    }
    units_d = {
#        'CO':'ppbv', 'O3':'ppbv', 'NO2':'pptv', 'NO':'pptv', 'NOx':'pptv',
    'BrO':'pptv', 'HNO3':'pptv', 'HNO2':'pptv', 'HONO':'pptv',
    }
    range_d = {
#        'CO':(50, 400), 'O3':(-10, 100), 'NO2':(-50, 500), 'NO':(-50, 500),
#        'NOx':(-50, 500),
    'HNO2':(-10, 60),
    'HNO3':(-30, 1500),
    'BrO':(-0.2, 1.0),
    'HONO':(-10, 60),
    }
    # - by variable
    runs = list(sorted(data_d.keys()))
    # Which variables to use?
    vars2plot = list(sorted(mod2obs_varnames.keys()))[::-1]
#    vars2plot = ['CO', 'O3', 'NOx', 'NO2', 'NO', 'HNO2']
    print(vars2plot)
    print(df_obs.columns)
    vars2plot = [
    i for i in vars2plot if mod2obs_varnames[i] in df_obs.columns
    ]
    # What bins should be used?
    bins = [0.5*i for i in np.arange(15)]
    for var2plot in vars2plot:
        fig = plt.figure()
        ax = plt.gca()
        # Now loop data
        for n_key, key_ in enumerate(runs):
            print(n_key, key_, var2plot )
            #
            if ('Obs.' in key_):
                varname = mod2obs_varnames[var2plot]
            else:
                varname = var2plot
            # Setup an axis label
            units = units_d[var2plot]
            xlabel = '{} ({})'.format( var2plot, units )
            # Add alt to DataFrame
            df = pd.DataFrame({
            var2plot: data_d[key_][varname], ALT_var: data_d[key_][ALT_var]
            })
            print(df.describe())
            # Scale the model values to pptv/ppbv
            if ('Obs.' not in key_):
                scaleby = AC.get_unit_scaling(units)
                df[var2plot] = df[var2plot].values * scaleby

            # drop any NaNs from the DataFrame
            s_shape = df.shape
            df.dropna(axis=0, how='any', inplace=True)
            if s_shape != df.shape:
                pcent = (float(df.shape[0]) - s_shape[0])/s_shape[0] * 100.
                pstr_dtr = 'WANRING dropped values - shape {}=>{} ({:.2f} %)'
                print(pstr_dtr.format(s_shape, df.shape, pcent))
            # Plot up as binned boxplots using existing function
            try:
                AC.binned_boxplots_by_altitude(df=df, fig=fig, ax=ax,
                                               var2bin_by=ALT_var,
                                               label=key_, xlabel=xlabel,
                                               binned_var=var2plot,
                                               num_of_datasets=len(runs),
                                               bins=bins,
                                               widths=0.15,
                                               dataset_num=n_key,
                                               color=color_dict[key_])
            except:
                pass
            # Make NOx species be on a log scale
            if (var2plot in NOx_specs):
                ax.set_xscale('log')
            #                     ax.set_xlim( (1E-5, 1E3) )
                ax.set_xlim( (0.3, 400) )
            else:
                ax.set_xscale('linear')
            # Beautify plot
            plt.legend()
            plt.title(title_str.format(var2plot, units, flight_ID ))
            plt.xlim(range_d[var2plot])

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()

    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plt_comp_by_alt_4ARNA_flights(dpi=320, just_SLR=True, show_plot=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    import seaborn as sns
    # Which flights to plot?
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
    # Just use non-transit ARNA flights
    flights_nums = [
    217, 218, 219, 220, 221, 222, 223, 224, 225,
    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_FAAM_core4flightnum(flight_ID=flight_ID)
#        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        dfs_obs[flight_ID] = df
    # Observations
    dfs_mod = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        dfs_mod[flight_ID] = df

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        df_mod = dfs_mod[flight_ID]
        # Only consider data during SLRs?
        if just_SLR:
            df_obs = df_obs.loc[df_obs['IS_SLR'] == True, :]
            df_mod = df_mod.loc[df_mod['IS_SLR'] == True, :]
            extr_str = '_JUST_SLR'
        else:
            extr_str = ''
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_altitude_binned_{}'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID)
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Put observations and vars to plot into a dictionary
        sns.set(color_codes=True)
        sns.set_context("paper", font_scale=0.75)
        # Force alt to be in units of km
        ALT_var = 'Altitude (km)'
        Y_unit = ALT_var
        df_mod[ALT_var] = AC.hPa_to_Km( df_mod['model-lev'].values )
        df_obs[ALT_var] = df_obs['ALT_GIN'].values / 1E3
        #
        data_d = {'GEOS-CF': df_mod, 'Obs.':df_obs}

        # - Now plot up flight time series plots by variable
        title_str =  "Altitude binned '{}' ({}) during flight '{}'"
        # Setup color dictinoary
        color_dict = {'GEOS-CF': 'red', 'Obs.':'k'}
        unit_d = {}
        mod2obs_varnames = {
        'CO':'CO_AERO', 'O3':'O3_TECO', 'NO2':'no2_mr', 'NO':'no_mr',
        'HNO2':'hono_mr',
        'NOx':'NOx'
        }
        units_d = {
        'CO':'ppbv', 'O3':'ppbv', 'NO2':'pptv', 'NO':'pptv', 'NOx':'pptv',
        'HNO2':'pptv', 'HONO':'pptv',
        }
        range_d = {
        'CO':(50, 400), 'O3':(-10, 100), 'NO2':(-50, 500), 'NO':(-50, 500),
        'NOx':(-50, 500),
        'HNO2':(-60, 60), 'HONO':(-60, 60),
        }
        NOx_specs = ['HNO2', 'NOx', 'NO', 'NO2', 'HONO']
        # - by variable
        runs = list(sorted(data_d.keys()))
		# Which variables to use?
#        vars2plot = list(sorted(mod2obs_varnames.keys()))[::-1]
        vars2plot = ['CO', 'O3', 'NOx', 'NO2', 'NO', 'HNO2']
        print(vars2plot)
        print(df_obs.columns)
        vars2plot = [
        i for i in vars2plot if mod2obs_varnames[i] in df_obs.columns
        ]
        # What bins should be used?
        bins = [0.5*i for i in np.arange(15)]
        for var2plot in vars2plot:
            fig = plt.figure()
            ax = plt.gca()
            # Now loop data
            for n_key, key_ in enumerate(runs):
                print(n_key, key_, var2plot )
                #
                if key_ == 'Obs.':
                    varname = mod2obs_varnames[var2plot]
                else:
                    varname = var2plot
                # Setup an axis label
                units = units_d[var2plot]
                xlabel = '{} ({})'.format( var2plot, units )
                # Add alt to DataFrame
                df = pd.DataFrame({
                var2plot: data_d[key_][varname], ALT_var: data_d[key_][ALT_var]
                })
                #
                if key_ != 'Obs.':
                    scaleby = AC.get_unit_scaling(units)
                    df[var2plot] = df[var2plot].values * scaleby

                # drop any NaNs from the DataFrame
                s_shape = df.shape
                df.dropna(axis=0, how='any', inplace=True)
                if s_shape != df.shape:
                    pcent = (float(df.shape[0]) - s_shape[0])/s_shape[0] * 100.
                    pstr_dtr = 'WANRING dropped values - shape {}=>{} ({:.2f})'
                    print(pstr_dtr.format(s_shape, df.shape, pcent))
                # Plot up as binned boxplots using existing function
                try:
                    AC.binned_boxplots_by_altitude(df=df, fig=fig, ax=ax,
                                                   var2bin_by=ALT_var,
                                                   label=key_, xlabel=xlabel,
                                                   binned_var=var2plot,
                                                   num_of_datasets=len(runs),
                                                   bins=bins,
                                                   widths = 0.15,
                                                   dataset_num=n_key,
                                                   color=color_dict[key_])
                except:
                    pass
                # Make NOx species be on a log scale
                if (var2plot in NOx_specs):
                    ax.set_xscale('log')
                #                     ax.set_xlim( (1E-5, 1E3) )
                    ax.set_xlim( (0.3, 400) )
                else:
                    ax.set_xscale('linear')
                # Beautify plot
                plt.legend()
                plt.title(title_str.format(var2plot, units, flight_ID ))
                plt.xlim(range_d[var2plot])

            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()

        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_comp_by_alt_4ARNA_flights_CIMS(dpi=320, just_SLR=False,
                                       show_plot=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context("paper", font_scale=0.75)
    # Which flights to plot?
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
    # Just use non-transit ARNA flights
    flights_nums = [
#    217, # Missing data for C217 (NOy)
    218, 219, 220,
#    221, # Missing data for C221 (NOy)
    222, 223,
#    224,  # Missing data for C221 (BrO... )
#    225,  # Missing data for C221 (BrO... )
    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]

    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
#        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        dfs_mod[flight_ID] = df

    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_CIMS_data4flight(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        dfs_obs[flight_ID] = df

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        df_mod = dfs_mod[flight_ID]
        # Only consider data during SLRs?
        if just_SLR:
            df_obs = df_obs.loc[df_obs['IS_SLR'] == True, :]
            df_mod = df_mod.loc[df_mod['IS_SLR'] == True, :]
            extr_str = '_JUST_SLR'
        else:
            extr_str = ''
        # Setup PDF to save PDF plots to
        savetitle_str = 'ARNA_altitude_binned_{}_CIMS{}'
        savetitle = savetitle_str.format(flight_ID, extr_str)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID)
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Put observations and vars to plot into a dictionary
        sns.set(color_codes=True)
        sns.set_context("paper", font_scale=0.75)
        # Force alt to be in units of km
        ALT_var = 'Altitude (km)'
        Y_unit = ALT_var
        df_mod[ALT_var] = AC.hPa_to_Km( df_mod['model-lev'].values )
        df_obs[ALT_var] = df_obs['ALT_GIN'].values / 1E3
        #
        data_d = {'GEOS-CF': df_mod, 'Obs.':df_obs}

        # - Now plot up flight time series plots by variable
        title_str =  "Altitude binned '{}' ({}) during flight '{}'"
        # Setup color dictinoary
        color_dict = {'GEOS-CF': 'red', 'Obs.':'k'}
        unit_d = {}
        mod2obs_varnames = {
        'BrO':'BrO',
        'HNO3':'HNO3',
        'HNO2':'HONO',
#        'CO':'CO_AERO', 'O3':'O3_TECO', 'NO2':'no2_mr', 'NO':'no_mr',
#        'HNO2':'hono_mr',
#        'NOx':'NOx'
        }
        units_d = {
#        'CO':'ppbv', 'O3':'ppbv', 'NO2':'pptv', 'NO':'pptv', 'NOx':'pptv',
        'BrO':'pptv', 'HNO3':'pptv', 'HNO2':'pptv', 'HONO':'pptv',
        }
        range_d = {
#        'CO':(50, 400), 'O3':(-10, 100), 'NO2':(-50, 500), 'NO':(-50, 500),
#        'NOx':(-50, 500),
        'HNO2':(-10, 60),
        'HNO3':(-30, 1500),
        'BrO':(-0.2, 1.0),
        'HONO':(-10, 60),
        }
        # - by variable
        runs = list(sorted(data_d.keys()))
		# Which variables to use?
        vars2plot = mod2obs_varnames.keys()
        print(vars2plot)
        print(df_obs.columns)
        vars2plot = [
        i for i in vars2plot if mod2obs_varnames[i] in df_obs.columns
        ]
        # What bins should be used?
        bins = [0.5*i for i in np.arange(15)]
        for var2plot in vars2plot:
            fig = plt.figure()
            ax = plt.gca()
            # Now loop data
            for n_key, key_ in enumerate(runs):
                print(n_key, key_, var2plot )
                #
                if key_ == 'Obs.':
                    varname = mod2obs_varnames[var2plot]
                else:
                    varname = var2plot
                # Setup an axis label
                units = units_d[var2plot]
                xlabel = '{} ({})'.format( var2plot, units )
                # Add alt to DataFrame
                df = pd.DataFrame({
                var2plot: data_d[key_][varname], ALT_var: data_d[key_][ALT_var]
                })
                #
                if key_ != 'Obs.':
                    scaleby = AC.get_unit_scaling(units)
                    df[var2plot] = df[var2plot].values * scaleby

                # drop any NaNs from the DataFrame
                s_shape = df.shape
                df.dropna(axis=0, how='any', inplace=True)
                if s_shape != df.shape:
                    pcent = (float(df.shape[0]) - s_shape[0])/s_shape[0] * 100.
                    pstr_dtr = 'WANRING dropped values - shape {}=>{} ({:.2f})'
                    print(pstr_dtr.format(s_shape, df.shape, pcent))
                # Plot up as binned boxplots using existing function
                try:
                    AC.binned_boxplots_by_altitude(df=df, fig=fig, ax=ax,
                                                   var2bin_by=ALT_var,
                                                   label=key_, xlabel=xlabel,
                                                   binned_var=var2plot,
                                                   num_of_datasets=len(runs),
                                                   bins=bins,
                                                   widths = 0.15,
                                                   dataset_num=n_key,
                                                   color=color_dict[key_])
                # Make NOx species be on a log scale
#                 if spec in NOx_specs:
#                     ax.set_xscale('log')
#                     ax.set_xlim( (1E-5, 1E3) )
#                 else:
                    ax.set_xscale('linear')
                except:
                    pass

                # Beautify plot
                plt.legend()
                plt.title(title_str.format(var2plot, units, flight_ID ))
                plt.xlim(range_d[var2plot])

            # Save to PDF
        #        fig.legend(loc='best', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
        #        plt.legend()
        #        plt.tight_layout()
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            if show_plot:
                plt.show()
            plt.close()

        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_flightpath_spatially_over_CVAO(df, LatVar='LAT_GIN', LonVar='LON_GIN',
                                       title=None, flight_ID='', dpi=320):
    """
    Plot flightpath spatially in CVAO region
    """
    # Reset sns for spatial plot
    sns.reset_orig()
    # New figure
    fig = plt.figure()
    # Get lat and lons
    lons = df[LonVar].values
    lats = df[LatVar].values
    # Get dates/datetimes of flight
    sdate_str = df.index.min().strftime('%x').strip()
    edate_str = df.index.max().strftime('%x').strip()
    stime_str = df.index.min().strftime('%H:%M').strip()
    etime_str = df.index.max().strftime('%H:%M').strip()
    # Make title string
    if isinstance(title, type(None)):
        if sdate_str != edate_str:
            title_str = 'Flight track for ARNA flight {} ({} {}-{} {})'
            title = title_str.format(flight_ID, sdate_str, stime_str,
                                     edate_str, etime_str)
        else:
            title_str = 'Flight track for ARNA flight {} ({}, {}-{})'
            title = title_str.format(flight_ID, sdate_str, stime_str,
                                     etime_str)
    # Setup projection for plotting
    projection = ccrs.PlateCarree
    central_longitude = 0
    fig = plt.figure(dpi=dpi, facecolor='w', edgecolor='k')
    # Setup a cartopy projection for plotting
    ax = fig.add_subplot(111,
                         projection=projection(
                         central_longitude=central_longitude)
                         )
    # Beautify
    ax.gridlines()
    # Add borders for countries
    ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                   facecolor='none', zorder=50)
    # Also add minor islands (inc. Cape Verde)
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                            edgecolor=None,
                                            facecolor='none')
    ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)
    # Plot settings
    marker='o'
    s=2
    alpha = 0.75
    cmap_list = AC.get_CB_color_cycle()
    # Now plot locations as scatter points on plot
    ax.scatter(lons, lats, color=cmap_list[0], s=s, marker=marker,
               alpha=alpha,
               label=flight_ID,
               transform=projection(), zorder=999
               )
    # Local area analysed as Cape Verde
    x0 = -30
    x1 =-10
    y0 = 0
    y1 = 25
    extents = (x0, x1, y0, y1)
    # Get limits of plotting data
    if isinstance(extents, type(None)):
        x0 = float(ds[LonVar].min())
        x1 = float(ds[LonVar].max())
        y0 = float(ds[LatVar].min())
        y1 = float(ds[LatVar].max())
        extents = (x0, x1, y0, y1)
    ax.set_extent(extents, crs=ccrs.PlateCarree())
    # Add a title to the plot
    plt.title(title)
    plt.tight_layout()


def plt_timeseries4ARNA_flight(df_obs=None, df_mod=None,
                               obs_label='Obs.',
                               mod_scale=1, obs_scale=1,
                               obs_adjustby=0, mod_adjustby=0,
                               ylim=(None, None),
                               units='ppbv', var2plot='CO',
                               obs_var2plot='CO_AERO',
                               mod_var2plot='CO',
                               mod_label='GEOS-CF',
                               plt_alt_as_shadow=True,
                               plt_dust_as_backfill=True,
                               plt_BB_as_backfill=False,
                               aspect_ratio=0.25,
                               flight_ID='C216',
                               yscale='linear',
                               invert_yaxis=False,
                               title=None,
                               plt_legend=True,
                               ):
    """
    Plot up a timeseries of observations and model for a given flight
    """
    # Get the beginning and end of the flight from the extracted model times
    xylim_min = AC.add_minutes( df_mod.index.min(), -15)
    xylim_max = AC.add_minutes( df_mod.index.max(), 15 )
    xticks = df_mod.resample('15T' ).mean().index.values
    xticks = AC.dt64_2_dt( xticks )
    xticks_labels = [ i.strftime('%H:%M') for i in xticks]
    # Get dates/datetimes of flight
    sdate_str = df_obs.index.min().strftime('%x').strip()
    edate_str = df_obs.index.max().strftime('%x').strip()
    stime_str = df_obs.index.min().strftime('%H:%M').strip()
    etime_str = df_obs.index.max().strftime('%H:%M').strip()
    # Now use Seaborn settings
    sns.set(color_codes=True)
    sns.set_context("paper", font_scale=0.75)
    # Set a shared title string
    # Setup the figure
    w, h = matplotlib.figure.figaspect(aspect_ratio)
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_subplot(111)
    # Plot up where
    if plt_dust_as_backfill:
        colour_plot_background_by_bool(df=df_obs, ax=ax,
                                       bool2use='IS_DUST',
                                       color='sandybrown',
                                       alpha=0.25,
                                       label='Dust (non-MBL)')
    if plt_BB_as_backfill:
        # TODO add plotting of biomass as backfill too
        pass

    # Plot up the observations and model...
    if not isinstance(obs_var2plot, type(None)):
        df_obs = df_obs[[obs_var2plot]].dropna()
        plt.plot(df_obs.index,
                 (df_obs[obs_var2plot].values+obs_adjustby)*obs_scale,
                 label=obs_label, color='k' )
    if not isinstance(mod_var2plot, type(None)):
        plt.plot(df_mod.index,
                 (df_mod[ mod_var2plot ].values+mod_adjustby)*mod_scale,
                 label=mod_label, color='red' )
    # Beautify plot
    if isinstance(title, str):
        plt.title(title)
#    elif isinstance(title, type(None)):
    else:
        title_str =  "Timeseries of '{}' ({}) during flight '{}' on {}"
        plt.title(title_str.format(var2plot, units, flight_ID, sdate_str ))

    plt.yscale(yscale)
    plt.ylim(ylim)
    plt.ylabel( '{} ({})'.format( var2plot, units) )
    plt.xlim(xylim_min, xylim_max)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks_labels, rotation=45)
#    print(xticks_labels)
    # Invert the second y-axis
    if plt_alt_as_shadow:
        # Add a shadow of the altitude
        ax2 = ax.twinx()
        mod_var2plot = 'model-lev'
        ax2.plot(df_mod.index, df_mod[ mod_var2plot ].values,
                 label='Altitude',
                 color='grey', zorder=100, alpha=0.25  )
        ax2.set_ylabel('Altitude (hPa)')
        ax2.grid(None)
        ax2.invert_yaxis()
        # Force use of the same ticks
        ax2.set_xticks(xticks)
        ax2.set_xticklabels(xticks_labels, rotation=45)
    # Invert the y-xais?
    if invert_yaxis:
        plt.gca().invert_yaxis()
    # Save to PDF
    if plt_legend:
        fig.legend(loc='best', bbox_to_anchor=(1,1),
                   bbox_transform=ax.transAxes)
    plt.tight_layout()


def plt_timeseries4ARNA_flight_point_obs(df_obs=None, df_mod=None,
                                         obs_label='Obs.',
                                         mod_scale=1, obs_adjustby=0,
                                         ylim=(None, None),
                                         units='ppbv', var2plot='CO',
                                         obs_var2plot='CO_AERO',
                                         mod_var2plot='CO',
                                         mod_label='GEOS-CF',
                                         plt_alt_as_shadow=True,
                                         aspect_ratio=0.25,
                                         flight_ID='C216',
                                         yscale='linear',
                                         invert_yaxis=False,
                                         plt_dust_as_backfill=True,
                                         plt_errorbar=False,
                                         obs_var2plot_err='',
                                         ):
    """
    Plot up a timeseries of observations and model for a given flight
    """
    # Get the beginning and end of the flight from the extracted model times
    xylim_min = AC.add_minutes( df_mod.index.min(), -15)
    xylim_max = AC.add_minutes( df_mod.index.max(), 15 )
    xticks = df_mod.resample('15T' ).mean().index.values
    xticks = AC.dt64_2_dt( xticks )
    xticks_labels = [ i.strftime('%H:%M') for i in xticks]
    # Get dates/datetimes of flight
    sdate_str = df_mod.index.min().strftime('%x').strip()
    edate_str = df_mod.index.max().strftime('%x').strip()
    stime_str = df_mod.index.min().strftime('%H:%M').strip()
    etime_str = df_mod.index.max().strftime('%H:%M').strip()
    # Now use Seaborn settings
    sns.set(color_codes=True)
    sns.set_context("paper", font_scale=0.75)
    # Set a shared title string
    title_str =  "Timeseries of '{}' ({}) during flight '{}' on {}"
    # Setup the figure
    w, h = matplotlib.figure.figaspect(aspect_ratio)
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_subplot(111)
    # Plot up the model data...
    if not isinstance(mod_var2plot, type(None)):
        plt.plot(df_mod.index, df_mod[ mod_var2plot ].values*mod_scale,
                 label=mod_label, color='red' )

    # Plot up where
    if plt_dust_as_backfill:
        colour_plot_background_by_bool(df=df_obs, ax=ax,
                                       bool2use='IS_DUST',
                                       color='sandybrown',
                                       alpha=0.25,
                                       label='Dust (non-MBL)')
    # Now plot up the observations
    if plt_errorbar:
        ax = plt.gca()
        ax.errorbar(df_obs.index, df_obs[obs_var2plot].values+obs_adjustby,
                    yerr=df_obs[obs_var2plot_err].values, fmt='-o',
                    color='k', capsize=3.0, label=obs_label)
    else:
#        df_obs = df_obs[[obs_var2plot]].dropna()
        plt.scatter(df_obs.index, df_obs[obs_var2plot].values+obs_adjustby,
                    label=obs_label, color='k' )

    # Beautify plot
    plt.title(title_str.format(var2plot, units, flight_ID, sdate_str ))
    plt.yscale(yscale)
    plt.ylim(ylim)
    plt.ylabel( '{} ({})'.format( var2plot, units) )
    plt.xlim(xylim_min, xylim_max)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks_labels, rotation=45)
#    print(xticks_labels)
    # Invert the second y-axis
    if plt_alt_as_shadow:
        # Add a shadow of the altitude
        ax2 = ax.twinx()
        mod_var2plot = 'model-lev'
        ax2.plot(df_mod.index, df_mod[ mod_var2plot ].values,
                 label='Altitude',
                 color='grey', zorder=100, alpha=0.25  )
        ax2.set_ylabel('Altitude (hPa)')
        ax2.grid(None)
        ax2.invert_yaxis()
        # Force use of the same ticks
        ax2.set_xticks(xticks)
        ax2.set_xticklabels(xticks_labels, rotation=45)
    # Invert the y-xais?
    if invert_yaxis:
        plt.gca().invert_yaxis()
    # Save to PDF
    fig.legend(loc='best', bbox_to_anchor=(1,1),
               bbox_transform=ax.transAxes)
    plt.tight_layout()


def add_deriv_vars2df(df):
    """
    Add derived/combined model variables to dataframe
    """
    # Add HOBr+Br2 (measured by the ToF-CIMS)
    try:
        df['Br2+HOBr'] = df['Br2']+df['HOBr']
    except KeyError:
        print("Derived variable not added to dataframe ('Br2+HOBr')")
    # total SO4
    try:
        df['SO4.total'] = df['SO4']+df['SO4s']
    except KeyError:
        print("Derived variable not added to dataframe ('SO4.total')")
    # total NIT
    try:
        df['NIT.total'] = df['NIT']+df['NITs']
    except KeyError:
        print("Derived variable not added to dataframe ('NIT.total')")
    return df


def colour_plot_background_by_bool(df, ax=None, bool2use=None, color='grey',
                                   alpha=0.25, label=None):
    """
    Colour background of plot for boolean locations on axis
    """
    # fill between plot max and min Y coordinate, for locations bool==True
    ymin, ymax = 0, 1
    x = df.index
    ax.fill_between(x, ymin, ymax, where=df[bool2use].values,
                    color=color, alpha=alpha, label=label,
                    transform=ax.get_xaxis_transform())
    return ax


def add_derived_FAAM_flags2df4flight(df=None, df_FAAM=None, flight_ID='C217'):
    """
    Add the derived FAAM variables to dataframe (e.g. IS_DUST)
    """
    # Get the relevant FAAM dataframe is not provided
    if isinstance(df_FAAM, type(None)):
#        df = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
        df_FAAM = get_FAAM_core4flightnum(flight_ID=flight_ID,
                                          resample_data=True )
    # Just add is the "is dust" and "is SLR"
    vars2use = ['IS_DUST', 'IS_SLR']
    df = pd.concat([df, df_FAAM[vars2use]],  axis=1)
    return df


def add_secs2duplicate_index_values(df):
    """
    duplicate values for index are not limit panadas operations
    """
    # Add a number to each ro
    RowVar = 'Row'
    df[RowVar] = np.arange(df.shape[0])
    # Loop index and add a second to any duplicates
    NewIndex = []
    last_index = None
    penultimate_index = None
    for Row in df[RowVar].values:
        index_val = df.loc[df[RowVar]== Row, :].index.values[0]
        index_val = AC.dt64_2_dt([index_val])[0]
        duplicates = df.loc[ df.index == index_val,:].shape[0]
        if duplicates >1:
            print(Row, df.loc[df[RowVar]==Row, : ])
            index_val = AC.add_secs(index_val, 1)
            # Update the index
            index_values =  df.index.values
            index_values[Row] = index_val
            df.index = index_values
    return df


def plt_timeseries_comp4ARNA_flights_CIMS(dpi=320, show_plot=False):
    """
    Plot up timeseries comparisons between core observations and model data
    """
    import seaborn as sns
    # Which flights to plot?
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
	# Just use non-transit ARNA flights
    flights_nums = [
#    217,
    218, 219, 220, 221, 222, 223, 224, 225,
    ]
    # Use flightnumbers with both NOy and halogens data
#    flights_nums = [
#    217, # Missing data for C217 (NOy)
#    218, 219, 220,
#    221, # Missing data for C221 (NOy)
#    222, 223,
#    224,  # Missing data for C221 (BrO... )
#    225,  # Missing data for C221 (BrO... )
#    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
        # Add the derived variables to the dataframe
        df = add_deriv_vars2df(df=df)
        dfs_mod[flight_ID] = df

    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_CIMS_data4flight(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        dfs_obs[flight_ID] = df

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        df_mod = dfs_mod[flight_ID]
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_timeseries_flighttrack_{}_CIMS'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        try:
            plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID)
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format('spatial', flight_ID))
        # - Now plot up flight time series plots by variable
        # - Plot up Bromine monoxide
        try:
            units = 'ppt'
            var2plot = 'BrO'
            obs_var2plot = 'BrO'
            mod_var2plot = 'BrO'
            mod_label = 'GEOS-CF'
            mod_scale = 1E12
            ylim = (-0.2, 1)
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))

        # - Plot up Nitric acid
        try:
            units = 'ppt'
            var2plot = 'HNO3'
            obs_var2plot = 'HNO3'
            mod_var2plot = 'HNO3'
            mod_label = 'GEOS-CF'
            mod_scale = 1E12
            ylim = (-30, 1500)
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))

        # - Plot up Nitrous acid
        try:
            units = 'ppt'
            var2plot = 'HONO'
            obs_var2plot = 'HONO'
            mod_var2plot = 'HNO2'
            mod_label = 'GEOS-CF'
            mod_scale = 1E12
            ylim = (-10, 60)
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))

        # - Plot up HCN
        try:
            units = 'ppt'
            var2plot = 'HCN'
            obs_var2plot = 'HCN'
            mod_var2plot = None
            mod_label = 'GEOS-CF'
            mod_scale = None
            ylim = (-10, 120)
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))

        # - Add a plot with masking for biomass burning.
        try:
            #
            units = 'ppt'
            var2plot = 'HCN'
            obs_var2plot = 'HCN'
            mod_var2plot = None
            mod_label = 'GEOS-CF'
            mod_scale = None
            ylim = (-10, 120)
            title_str =  "Timeseries of '{}' ({}) during flight '{}' - {}"
            ext_str = 'Biomass flagged as HCN @ background+3$\sigma$'
            title = title_str.format(var2plot, units, flight_ID, ext_str )
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       plt_legend=False,
                                       title=title,
                                       )
            # Add a flat for biomass burning
            # hardware the stats for now
            sigma = 20.082764920179258
            background = 7.5
            threshold = background+(3*sigma)
            df_obs = add_biomass_flag2df(df_obs, CIMSdf=df_obs,
                                         threshold=threshold,
                                         flight_ID=flight_ID)
            # Colour background based on this
            bool2use = 'IS_BB'
            ax = plt.gca()
            colour_plot_background_by_bool(df=df_obs, ax=ax, bool2use=bool2use,
                                           color='firebrick', alpha=0.25,
                                           label='Biomass Burning')
            fig = plt.gcf()
            fig.legend(loc='best', bbox_to_anchor=(1,1),
                       bbox_transform=ax.transAxes)
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))
        # - Add a plot with masking for biomass burning.
        try:
            #
            units = 'ppt'
            var2plot = 'HCN'
            obs_var2plot = 'HCN'
            mod_var2plot = None
            mod_label = 'GEOS-CF'
            mod_scale = None
            ylim = (-10, 120)
            title_str =  "Timeseries of '{}' ({}) during flight '{}' - {}"
            ext_str = 'Biomass flagged as HCN @ background+2$\sigma$'
            title = title_str.format(var2plot, units, flight_ID, ext_str )
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       plt_legend=False,
                                       title=title,
                                       )
            # Add a flat for biomass burning
            # hardware the stats for now
            sigma = 20.082764920179258
            background = 7.5
            threshold = background+(2*sigma)
            df_obs = add_biomass_flag2df(df_obs, CIMSdf=df_obs,
                                         threshold=threshold,
                                         flight_ID=flight_ID)
            # Colour background based on this
            bool2use = 'IS_BB'
            ax = plt.gca()
            colour_plot_background_by_bool(df=df_obs, ax=ax, bool2use=bool2use,
                                           color='firebrick', alpha=0.25,
                                           label='Biomass Burning')
            fig = plt.gcf()
            fig.legend(loc='best', bbox_to_anchor=(1,1),
                       bbox_transform=ax.transAxes)
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))
        # - Add a plot with masking for biomass burning.
        try:
            #
            units = 'ppt'
            var2plot = 'HCN'
            obs_var2plot = 'HCN'
            mod_var2plot = None
            mod_label = 'GEOS-CF'
            mod_scale = None
            ylim = (-10, 120)
            title_str =  "Timeseries of '{}' ({}) during flight '{}' - {}"
            ext_str = 'Biomass flagged as HCN @ background+1$\sigma$'
            title = title_str.format(var2plot, units, flight_ID, ext_str )
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       plt_legend=False,
                                       title=title,
                                       )
            # Add a flat for biomass burning
            # hardware the stats for now
            sigma = 20.082764920179258
            background = 7.5
            threshold = background+(1*sigma)
            df_obs = add_biomass_flag2df(df_obs, CIMSdf=df_obs,
                                         threshold=threshold,
                                         flight_ID=flight_ID)
            # Colour background based on this
            bool2use = 'IS_BB'
            ax = plt.gca()
            colour_plot_background_by_bool(df=df_obs, ax=ax, bool2use=bool2use,
                                           color='firebrick', alpha=0.25,
                                           label='Biomass Burning')
            fig = plt.gcf()
            fig.legend(loc='best', bbox_to_anchor=(1,1),
                       bbox_transform=ax.transAxes)
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))
        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_timeseries_comp4ARNA_flights_filters(dpi=320, show_plot=False,
                                             debug=False):
    """
    Plot up timeseries comparisons between filter samples and model data
    """
    import seaborn as sns
    # Which flights to plot?
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
	# Just use non-transit ARNA flights
    flights_nums = [
#    217,
    218, 219, 220, 221, 222, 223, 224, 225,
    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod = {}
    for flight_ID in flight_IDs:
        if debug:
            print(flight_ID)
        df = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
        # Add the derived variables to the dataframe
        df = add_deriv_vars2df(df=df)
        dfs_mod[flight_ID] = df
    # Observations
    dfs_obs = get_filters_data4flight()
    dfs_obs = add_secs2duplicate_index_values(dfs_obs)

    # Convert observation units into model units
    # unit on recipt were 'nanomoles/m3', which were updated to ug/m3
    # model units? 'pptv'
    NIT_obs_var = 'NO3.total'
    data =  dfs_obs[NIT_obs_var].values
    dfs_obs[NIT_obs_var] = AC.convert_ug_per_m3_2_ppbv(data, spec='NIT') *1E3
    SO4_obs_var = 'SO4.total'
    data =  dfs_obs[SO4_obs_var].values
    dfs_obs[SO4_obs_var] = AC.convert_ug_per_m3_2_ppbv(data, spec='SO4') *1E3

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs.loc[ dfs_obs['Flight']==flight_ID, :]
        df_mod = dfs_mod[flight_ID]
        # add the dust and SLR flags to the core dataframe
        df_obs = add_derived_FAAM_flags2df4flight(df=df_obs,
                                                  flight_ID=flight_ID)
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_timeseries_flighttrack_{}_FILTERS'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_mod, flight_ID=flight_ID,
                                           LatVar='model-lat',
                                           LonVar='model-lon',
                                           )
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()
        # - Now plot up flight time series plots by variable
        vars2plot =  'SO4', 'NIT'

        # - Plot up nitrate
        units = 'pptv'
        var2plot = 'NO3'
        mod_var2plot = 'NIT.total'
        obs_var2plot = NIT_obs_var
        mod_label = 'GEOS-CF'
        mod_scale = 1E12
        obs_var2plot_err = None
        plt_errorbar = False
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_timeseries4ARNA_flight_point_obs(var2plot=var2plot, units=units,
                                             obs_var2plot=obs_var2plot,
                                             mod_scale=mod_scale,
                                             mod_label=mod_label,
                                             mod_var2plot=mod_var2plot,
                                             ylim=ylim,
                                             df_mod=df_mod, df_obs=df_obs,
                                             flight_ID=flight_ID,
                                             obs_var2plot_err=obs_var2plot_err,
                                             plt_errorbar=plt_errorbar,
                                             )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Plot up sulfate
        units = 'pptv'
        var2plot = 'SO4'
        mod_var2plot = 'SO4.total'
        obs_var2plot = SO4_obs_var
        mod_label = 'GEOS-CF'
        mod_scale = 1E12
        obs_var2plot_err = None
        plt_errorbar = False
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_timeseries4ARNA_flight_point_obs(var2plot=var2plot, units=units,
                                             obs_var2plot=obs_var2plot,
                                             mod_scale=mod_scale,
                                             mod_label=mod_label,
                                             mod_var2plot=mod_var2plot,
                                             ylim=ylim,
                                             df_mod=df_mod, df_obs=df_obs,
                                             flight_ID=flight_ID,
                                             obs_var2plot_err=obs_var2plot_err,
                                             plt_errorbar=plt_errorbar,
                                             )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_timeseries_comp4ARNA_flights_SWAS(dpi=320, show_plot=False,
                                          debug=False):
    """
    Plot up timeseries comparisons between SWAS observations and model data
    """
    import seaborn as sns
    # Which flights to plot?
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
	# Just use non-transit ARNA flights
    flights_nums = [
    217, 218, 219, 220, 221, 222, 223, 224, 225,
    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod = {}
    for flight_ID in flight_IDs:
        if debug:
            print(flight_ID)
        df = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
        # Add the derived variables to the dataframe
        df = add_deriv_vars2df(df=df)
        dfs_mod[flight_ID] = df

    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_SWAS_data4flight(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        dfs_obs[flight_ID] = df

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        df_mod = dfs_mod[flight_ID]
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_timeseries_flighttrack_{}_SWAS'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_mod, flight_ID=flight_ID,
                                           LatVar='model-lat',
                                           LonVar='model-lon',
                                           )
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()
        # - Now plot up flight time series plots by variable
        vars2plot =  'ALD2', 'ACET', 'C2H6', 'C3H8'

        # - Plot up acetone
        units = 'ppt'
        var2plot = 'ACET'
        mod_var2plot = 'ACET'
        obs_var2plot = map_SWAS_var2GEOS_var(mod_var2plot, invert=True)
        mod_label = 'GEOS-CF'
        mod_scale = 1E9
        obs_var2plot_err = obs_var2plot+'_uncertainty'
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_timeseries4ARNA_flight_point_obs(var2plot=var2plot, units=units,
                                             obs_var2plot=obs_var2plot,
                                             mod_scale=mod_scale,
                                             mod_label=mod_label,
                                             mod_var2plot=mod_var2plot,
                                             ylim=ylim,
                                             df_mod=df_mod, df_obs=df_obs,
                                             flight_ID=flight_ID,
                                             obs_var2plot_err=obs_var2plot_err,
                                             plt_errorbar=True
                                             )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Plot up acetaldehyde
        units = 'ppt'
        var2plot = 'ALD2'
        mod_var2plot = 'ALD2'
        obs_var2plot = map_SWAS_var2GEOS_var(mod_var2plot, invert=True)
        mod_label = 'GEOS-CF'
        mod_scale = 1E9
        obs_var2plot_err = obs_var2plot+'_uncertainty'
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_timeseries4ARNA_flight_point_obs(var2plot=var2plot, units=units,
                                             obs_var2plot=obs_var2plot,
                                             mod_scale=mod_scale,
                                             mod_label=mod_label,
                                             mod_var2plot=mod_var2plot,
                                             ylim=ylim,
                                             df_mod=df_mod, df_obs=df_obs,
                                             flight_ID=flight_ID,
                                             obs_var2plot_err=obs_var2plot_err,
                                             plt_errorbar=True
                                             )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Plot up ethane
        units = 'ppt'
        var2plot = 'C2H6'
        mod_var2plot = 'C2H6'
        obs_var2plot = map_SWAS_var2GEOS_var(mod_var2plot, invert=True)
        mod_label = 'GEOS-CF'
        mod_scale = 1E9
        obs_var2plot_err = obs_var2plot+'_uncertainty'
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_timeseries4ARNA_flight_point_obs(var2plot=var2plot, units=units,
                                             obs_var2plot=obs_var2plot,
                                             mod_scale=mod_scale,
                                             mod_label=mod_label,
                                             mod_var2plot=mod_var2plot,
                                             ylim=ylim,
                                             df_mod=df_mod, df_obs=df_obs,
                                             flight_ID=flight_ID,
                                             obs_var2plot_err=obs_var2plot_err,
                                             plt_errorbar=True
                                             )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Plot up propane
        units = 'ppt'
        var2plot = 'C3H8'
        mod_var2plot = 'C3H8'
        obs_var2plot = map_SWAS_var2GEOS_var(mod_var2plot, invert=True)
        mod_label = 'GEOS-CF'
        mod_scale = 1E9
        obs_var2plot_err = obs_var2plot+'_uncertainty'
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_timeseries4ARNA_flight_point_obs(var2plot=var2plot, units=units,
                                             obs_var2plot=obs_var2plot,
                                             mod_scale=mod_scale,
                                             mod_label=mod_label,
                                             mod_var2plot=mod_var2plot,
                                             ylim=ylim,
                                             df_mod=df_mod, df_obs=df_obs,
                                             flight_ID=flight_ID,
                                             obs_var2plot_err=obs_var2plot_err,
                                             plt_errorbar=True
                                             )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_timeseries_comp4ARNA_flights_PHYSICAL_VARS(dpi=320, show_plot=False):
    """
    Plot up timeseries comparisons between physical variables and model data
    """
    import seaborn as sns
    # Which flights to plot?
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
	# Just use non-transit ARNA flights
    flights_nums = [
    217, 218, 219, 220, 221, 222, 223, 224, 225,
    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod = {}
    for flight_ID in flight_IDs:
        dfs_mod[flight_ID] = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        dfs_obs[flight_ID] = get_FAAM_core4flightnum(flight_ID=flight_ID )

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        df_mod = dfs_mod[flight_ID]
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_timeseries_flighttrack_{}_PHYSICAL_VARS'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Plot up ROLL_GIN
        try:
            units = 'ADD THIS'
    #        VarName = 'PCAS2CON'
    #        FlagName = 'PCAS2_FLAG'
            var2plot = 'ROLL_GIN'
            obs_var2plot = 'ROLL_GIN'
            mod_var2plot = None
            mod_label = 'GEOS-CF'
            mod_scale = None
    #        ylim = (10, 100)
    #        ylim = None
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
    #                                   ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       plt_legend=False,
                                       )
            # Colour in SLRs
            ax = plt.gca()
            colour_plot_background_by_bool(df=df_obs, ax=ax, bool2use='IS_SLR',
                                           color='blue', alpha=0.25,
                                           label='SLR')
            fig = plt.gcf()
            fig.legend(loc='best', bbox_to_anchor=(1,1),
                       bbox_transform=ax.transAxes)
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot ROLL_GIN')

        # - Plot up VELD_GIN
        try:
            units = 'ADD THIS'
    #        VarName = 'PCAS2CON'
    #        FlagName = 'PCAS2_FLAG'
            var2plot = 'VELD_GIN'
            obs_var2plot = 'VELD_GIN'
            mod_var2plot = None
            mod_label = 'GEOS-CF'
            mod_scale = None
    #        ylim = (10, 100)
            ylim = None
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       plt_legend=False,
                                       )
            # Colour in SLRs
            ax = plt.gca()
            colour_plot_background_by_bool(df=df_obs, ax=ax, bool2use='IS_SLR',
                                           color='blue', alpha=0.25,
                                           label='SLR')
            fig = plt.gcf()
            fig.legend(loc='best', bbox_to_anchor=(1,1),
                       bbox_transform=ax.transAxes)

            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot VELD_GIN')

        # - Plot up PCASP
        try:
            units = '# cm$^{-1}$'
    #        VarName = 'PCAS2CON'
    #        FlagName = 'PCAS2_FLAG'
            var2plot = 'PCAS'
            obs_var2plot = 'PCAS2CON'
            mod_var2plot = None
            mod_label = 'GEOS-CF'
            mod_scale = None
    #        ylim = (10, 100)
            ylim = None
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up PCASP surface area
        try:
            units = 'x10$^{-12}$ m${^3}$/cm$^{-3}$'
    #        VarName = 'PCAS2CON'
    #        FlagName = 'PCAS2_FLAG'
            var2plot = 'PCASP (total surface area)'
            obs_var2plot = 'PCASP-total-surface'
            mod_var2plot = None
            mod_label = 'GEOS-CF'
            obs_scale = 1E12
            mod_scale = None
    #        ylim = (10, 100)
            ylim = None
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       obs_scale=obs_scale,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up PCASP surface area
        try:
            units = 'x10$^{-9}$ m${^3}$/cm$^{-3}$'
            var2plot = 'CDP (total surface area)'
            obs_var2plot = 'CDP-total-surface'
            mod_var2plot = None
            mod_label = 'GEOS-CF'
            obs_scale = 1E9
            mod_scale = None
    #        ylim = (10, 100)
            ylim = None
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       obs_scale=obs_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up temperature
        units = '$^{\circ}$C'
        var2plot = 'Temperature'
        obs_var2plot = 'TAT_DI_R'
        mod_var2plot = 'T'
        mod_label = 'GEOS-CF'
        mod_scale = 1
        ylim = (-25, 35)
        obs_adjustby = -273.15
        # Call timeseries plotter function
        plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                   obs_var2plot=obs_var2plot,
                                   mod_scale=mod_scale,
                                   mod_label=mod_label,
                                   mod_var2plot=mod_var2plot,
                                   ylim=ylim,
                                   df_mod=df_mod, df_obs=df_obs,
                                   flight_ID=flight_ID,
                                   obs_adjustby=obs_adjustby,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()


        # - Plot up Eastward wind
        units = 'm s$^{-1}$'
        var2plot = 'Eastward wind'
        obs_var2plot = 'U_C'
        mod_var2plot = 'U'
        mod_label = 'GEOS-CF'
        mod_scale = 1
        ylim = (-25, 25)
        # Call timeseries plotter function
        plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                   obs_var2plot=obs_var2plot,
                                   mod_scale=mod_scale,
                                   mod_label=mod_label,
                                   mod_var2plot=mod_var2plot,
                                   ylim=ylim,
                                   df_mod=df_mod, df_obs=df_obs,
                                   flight_ID=flight_ID,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Plot up Northward wind
        units = 'm s$^{-1}$'
        var2plot = 'Northward wind'
        obs_var2plot = 'V_C'
        mod_var2plot = 'V'
        mod_label = 'GEOS-CF'
        mod_scale = 1
        ylim = (-25, 25)
        # Call timeseries plotter function
        plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                   obs_var2plot=obs_var2plot,
                                   mod_scale=mod_scale,
                                   mod_label=mod_label,
                                   mod_var2plot=mod_var2plot,
                                   ylim=ylim,
                                   df_mod=df_mod, df_obs=df_obs,
                                   flight_ID=flight_ID,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Plot up Latitude
        try:
            units = '$^{\circ}$N'
            var2plot = 'Latitude'
            obs_var2plot = 'LAT_GIN'
            mod_var2plot = 'model-lat'
            mod_label = 'GEOS-CF'
            mod_scale = 1
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up Longitude
        try:
            units = '$^{\circ}$E'
            var2plot = 'Longitude'
            obs_var2plot = 'LON_GIN'
            mod_var2plot = 'model-lon'
            mod_label = 'GEOS-CF'
            mod_scale = 1
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up altitude
        try:
            units = 'hPa'
            var2plot = 'Altitude'
            mod_var2plot = 'model-lev'
            obs_var2plot = 'PS_RVSM' # Use pressure measurement
#            obs_var2plot = 'ALT_GIN' # Use GPS altitude?
#            vals = AC.hPa_to_Km( df_obs['ALT_GIN'].values/1E3, reverse=True )
            mod_label = 'GEOS-CF'
            mod_scale = 1
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       invert_yaxis=True,
                                       plt_alt_as_shadow=False,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
            # Colour in SLRs
            ax = plt.gca()
            colour_plot_background_by_bool(df=df_obs, ax=ax, bool2use='IS_SLR',
                                           color='blue', alpha=0.25,
                                           label='SLR')
            plt.legend()
        except:
            print('Failed to plot {}'.format(var2plot))


        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_timeseries_comp4ARNA_flights(dpi=320, show_plot=False):
    """
    Plot up timeseries comparisons between core observations and model data
    """
    import seaborn as sns
    # Which flights to plot?
#    flights_nums = [ 216, 217, 218, 219, 220, 221, 222, 223, 224, 225 ]
	# Just use non-transit ARNA flights
    flights_nums = [
    217, 218, 219, 220, 221, 222, 223, 224, 225,
    ]
    flight_IDs = [ 'C{}'.format(i) for i in flights_nums ]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod = {}
    for flight_ID in flight_IDs:
        dfs_mod[flight_ID] = get_GEOSCF_output4flightnum(flight_ID=flight_ID )
    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        dfs_obs[flight_ID] = get_FAAM_core4flightnum(flight_ID=flight_ID )

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        df_mod = dfs_mod[flight_ID]
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_timeseries_flighttrack_{}'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Now plot up flight time series plots by variable
        # - Plot up carbon monoxide
        units = 'ppbv'
        var2plot = 'CO'
        obs_var2plot = 'CO_AERO'
        mod_var2plot = 'CO'
        mod_label = 'GEOS-CF'
        mod_scale = 1E9
        ylim = (50, 400)
        # Call timeseries plotter function
        plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                   obs_var2plot=obs_var2plot,
                                   mod_scale=mod_scale, mod_label=mod_label,
                                   mod_var2plot=mod_var2plot,
                                   ylim=ylim,
                                   df_mod=df_mod, df_obs=df_obs,
                                   flight_ID=flight_ID,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Plot up ozone
        units = 'ppbv'
        var2plot = 'Ozone'
        obs_var2plot = 'O3_TECO'
        mod_var2plot = 'O3'
        mod_label = 'GEOS-CF'
        mod_scale = 1E9
        ylim = (10, 100)
        # Call timeseries plotter function
        plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                   obs_var2plot=obs_var2plot,
                                   mod_scale=mod_scale, mod_label=mod_label,
                                   mod_var2plot=mod_var2plot,
                                   ylim=ylim,
                                   df_mod=df_mod, df_obs=df_obs,
                                   flight_ID=flight_ID,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
        plt.close()

        # - Plot up NO2
        try:
            units = 'pptv'
            var2plot = 'NO2'
            obs_var2plot = 'no2_mr'
            mod_var2plot = 'NO2'
            mod_label = 'GEOS-CF'
            mod_scale = 1E12
            ylim = (0.3, 400)
            yscale = 'log'
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       yscale=yscale,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot NO2')

        # - Plot up NO
        try:
            units = 'pptv'
            var2plot = 'NO'
            obs_var2plot = 'no_mr'
            mod_var2plot = 'NO'
            mod_label = 'GEOS-CF'
            mod_scale = 1E12
            ylim = (0.3, 400)
            yscale = 'log'
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       yscale=yscale,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot NO')

        # - Plot up NOx
        try:
            units = 'pptv'
            var2plot = 'NOx'
            mod_var2plot = 'NOx'
            obs_var2plot = mod_var2plot
            mod_label = 'GEOS-CF'
            mod_scale = 1E12
            ylim = (0.3, 400)
            yscale = 'log'
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       yscale=yscale,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot NOx')

        # - Plot up HNO2
        try:
            units = 'pptv'
            var2plot = 'HONO'
            obs_var2plot = 'hono_mr'
            mod_var2plot = 'HNO2'
            mod_label = 'GEOS-CF'
            mod_scale = 1E12
            ylim = (0.3, 400)
            yscale = 'log'
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       yscale=yscale,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot HONO')

        # - Plot up NOx / 'PCASP-total-surface'
        try:
            #
            obs_var2plot = 'NOx/PCASP'
            ratio = df_obs['NOx'] / df_obs['PCASP-total-surface']
            df_obs[obs_var2plot] = ratio
            # Then plot
            units = 'ratio'
            var2plot = obs_var2plot
            mod_var2plot = None
            obs_var2plot = obs_var2plot
            mod_label = 'GEOS-CF'
            mod_scale = None
#            ylim = (0.3, 400)
            yscale = 'log'
            # Call timeseries plotter function
            plt_timeseries4ARNA_flight(var2plot=var2plot, units=units,
                                       obs_var2plot=obs_var2plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       mod_var2plot=mod_var2plot,
                                       ylim=ylim,
                                       df_mod=df_mod, df_obs=df_obs,
                                       flight_ID=flight_ID,
                                       yscale=yscale,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi)
            plt.close()
        except:
            print('Failed to plot NOx')


        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plot_up_flight_locations_from_FAAM_website(d=None,
                                               LatVar='latitude',
                                               LonVar='longitude',
                                               folder='./'):
    """
    Plot up flight locations for all ARNA flights on a map
    """
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    # Get flight locations as a dataframe
    if isinstance(d, type(None)):
        d = get_ARNA_flights_as_dfs()
    # Setup plot
    projection=ccrs.PlateCarree
    central_longitude = 0
    fig = plt.figure(dpi=dpi, facecolor='w', edgecolor='k')
    # Setup a cartopy projection for plotting
    ax = fig.add_subplot(111,
                         projection=projection(
                             central_longitude=central_longitude)
                         )
    # Beautify
    ax.gridlines()
    # Add borders for countries
    ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                   facecolor='none', zorder=50)
    # Also add minor islands (inc. Cape Verde)
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                            edgecolor=None,
                                            facecolor='none')
    ax.add_feature(land_10m, edgecolor='grey', facecolor='none', zorder=50)
    # Plot settings
    marker='o'
    s=2
    alpha = 0.75
    cmap_list = AC.get_CB_color_cycle()
    #
    flight_IDs = list(sorted(d.keys()))
    for flight_ID_n, flight_ID in enumerate( flight_IDs ):
        df = d[flight_ID]
        # Save the data as a csv file?
        save2csv = True
        if save2csv:
            df.to_csv('TEMP_{}.csv'.format(flight_ID) )
        # Get lats and lons for flight
        lats = df[LatVar]
        lons = df[LonVar]
        # Now plot locations as scatter points on plot
        ax.scatter(lons, lats, color=cmap_list[flight_ID_n], s=s,
                   marker=marker,
                   alpha=alpha,
                   label=flight_ID,
                   transform=projection(), zorder=999
                   )
    # Local area analysed as Cape Verde
    x0 = -30
    x1 =-10
    y0 = 0
    y1 = 25
    extents = (x0, x1, y0, y1)
    # Get limits of plotting data
    if isinstance(extents, type(None)):
        x0 = float(ds[LonVar].min())
        x1 = float(ds[LonVar].max())
        y0 = float(ds[LatVar].min())
        y1 = float(ds[LatVar].max())
        extents = (x0, x1, y0, y1)
    ax.set_extent(extents, crs=ccrs.PlateCarree())
    # Include a legend
    plt.legend(ncol=3, loc="lower left", fontsize='small', handletextpad=0.05)
    # Add a tiutle
    plt.suptitle( u'ARNA campaign flights (excluding transits)')
    # Save the plot
    filename = 'ARNA_campaign_flight_locations_ALL'
    plt.savefig(folder+filename, pad_inches=0.25,
                bbox_inches='tight'
                )
    plt.close()


def quick_lat_plt_2layer(ds, var2plot1=None, var2plot2=None, extra_str='',
                         save_plot=True, show_plot=False, savename=None,
                         units=None, title=None, folder=None,
                         LatVar='lat', LonVar='lon', LevVar='lev',
                         fig=None, ax=None, extents=None,
                         add_ARNA_locs=True, use_local_CVAO_area=True,
                         extend='both', ylim=(0, 10), xlim=(-30, -16.5),
                         dpi=320):
    """
    Plot up a quick latitude-altitude plot of data using cartopy

    Parameters
    -------
    ds (xr.Dataset): dataset object holding data to plot
    var2plot (str): variable to plot within the dataset
    LatVar, LonVar (str): variables to use for latitude and longitude
    save_plot (bool): save the plot as a .png ?
    show_plot (bool): show the plot on screen
    dpi (int): resolution to use for saved image (dots per square inch)
    savename (str): name to use for png of saved .png
    extra_str (str): extra string to append to save .png
    projection (cartopy.crs obj.):  projection to use
    fig (figure instance): matplotlib figure instance
    ax (axis instance): axis object to use
    add_ARNA_locs (bool):
    use_local_CVAO_area (bool):

    Returns
    -------
    (None)
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # Use the 1st data variable if not variable given
    if isinstance(var2plot1, type(None)):
        print('WARNING: No variable to plot (var2plot), trying 1st data_var')
        var2plot1 = list(ds.data_vars)[0]
        var2plot1 = 'NOy'
    if isinstance(var2plot2, type(None)):
        print('WARNING: No variable to plot (var2plot), trying 1st data_var')
#        var2plot2 = list(ds.data_vars)[0]
        var2plot2 = 'PM2.5(dust)'
    # Setup figure and axis and plot
    if isinstance(fig, type(None)):
        fig = plt.figure(figsize=(12, 6))
#        fig = plt.figure()
    if isinstance(ax, type(None)):
#        ax = fig.add_subplot(111, projection=projection, aspect='auto')
        ax = fig.add_subplot(111)
#    if isinstance(xlim, type(None)):
#        xlim = (ds.lat.values.min(), ds.lat.values.max())
    # Local variables (metres to kft)
    m2kft = 3.281
    # Add a second y axis - tempoaroiliy hash out
#    ax2 = ax.twinx()
#        ax.set_yscale('log')
#        ax.set_xscale('linear')
    # Convert lev to metres using hydrostatic aprox
    convert2km = True
    if convert2km:
        attrs = ds.lev.attrs
        print(attrs)
        if attrs['units'] == 'millibar':
            hPa_l = list(ds.lev.values).copy()
            ds.lev.values = AC.hPa_to_Km(hPa_l)
            attrs['units'] = 'km'
        else:
            pass
    # Then convert into kft?
    convert2kft = True
    if convert2kft:
        attrs = ds.lev.attrs
        print(attrs)
        if attrs['units'] == 'km':
#            hPa_l = list(ds.lev.values).copy()
            ds.lev.values = ds.lev.values * m2kft
            attrs['units'] = 'kft'
        else:
            pass
        #
        ylim = (ylim[0]*m2kft, ylim[1]*m2kft)

    # Divider axis for colourbar
    divider = make_axes_locatable(ax)
    # - Plot first var
    #
#    cbar_kwargs = {  'extend': extend }
    # Setup plotted range
    alpha = 0.5
    vmin1, vmax1 = set_limits4ar_plotted_range(var2plot1)
    vmin1 = get_vmin_value4var(var2plot1)
    units1 = get_species_units(var2plot1)
    cmap, ticks, nticks = get_cmap4var(var2plot1)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap': cmap, 'extend' : 'both', 'pad':0.075,}
    else:
        cbar_kwargs  = {
        'ticks': ticks, 'cmap': cmap, 'extend' : 'both','pad':0.075,
        }
    # Now plot up var1 - using pcolormesh
    cbar_ax = divider.append_axes("right", "2%", pad="1%")
    # Now plot
    ds[var2plot1].plot.pcolormesh(x=LonVar, y=LevVar, ax=ax,
                                 vmin=vmin1, vmax=vmax1,
                                 zorder=1, alpha=alpha,
#                                 origin='lower',
                                 yincrease=True,
                                 cmap=cmap,
#                                 extend=extend,
                                 cbar_ax=cbar_ax,
                                 cbar_kwargs=cbar_kwargs,
                                 )
    # Remove the title
    ax.set_title("")

    # - Plot second var
    # Now plot up var 2
    alpha = 0.4
    vmin2, vmax2 = set_limits4ar_plotted_range(var2plot2)
    vmin2 = get_vmin_value4var(var2plot2)
    units2 = get_species_units(var2plot2)
    cmap, ticks, nticks = get_cmap4var(var2plot2)
    if isinstance(ticks, type(None)):
        cbar_kwargs  = { 'cmap': cmap, 'extend' : 'both', }
    else:
        cbar_kwargs  = {'ticks': ticks, 'cmap': cmap, 'extend' : 'both', }
    # Setup a colour bar axis
    cbar_ax = divider.append_axes("right", "2%", pad="7%")
    # using pcolormesh
    im = ds[var2plot2].plot.pcolormesh(x=LonVar, y=LevVar, ax=ax,
                                 vmin=vmin2, vmax=vmax2,
                                 zorder=1, alpha=alpha, cmap=cmap,
#                                 origin='lower',
                                 yincrease=True,
                                 cbar_kwargs=cbar_kwargs,
                                 cbar_ax=cbar_ax,
                                 )
    # Remove the title
    ax.set_title("")

    # Add a box to show flyable area
    # 8km alt (and lon range)
    d = get_max_flying_range4BAE146()
#    max_lat = d['max_lat']
#    min_lat = d['min_lat']
    max_lon = d['max_lon']
    min_lon = d['min_lon']
    if convert2km and (not convert2kft):
        max_alt = d['max_alt'] /1E3
    elif convert2kft:
        max_alt = (d['max_alt'] / 1E3 * m2kft)
    else:
        pass
#    print(max_lat, min_lat, max_alt  )
    max_alt_axis_coord = float(max_alt)/float(ylim[-1])
    print(max_alt_axis_coord)
    xrange = float(xlim[-1] - xlim[0])
    # Plot this up
    ax.axvline(x=max_lon, linestyle=':', alpha=0.5, color='grey', zorder=100,
               ymin=0, ymax=max_alt_axis_coord, linewidth=3.0)
    ax.axvline(x=min_lon, linestyle=':', alpha=0.5, color='grey', zorder=100,
               ymin=0, ymax=max_alt_axis_coord, linewidth=3.0)
    ax.axhline(y=max_alt, linestyle=':', alpha=0.5, color='grey', zorder=100,
#               xmin=min_lon/xrange, xmax=max_lon/xrange,
               xmin=(min_lon - xlim[0] ) / xrange,
               xmax=(max_lon - xlim[0] ) / xrange,
               linewidth=3.0)

    # Add locations for airports
#    locs2plot = ['DSS', 'RAI', 'VXE', 'LPA', 'LIS', 'CDG']
    locs2plot = [
    'Praia Airport',  'Sao Vicente Airport', 'Dakar',
#    'Gran Canaria Airport', 'Lisbon Airport',  'Paris (Charles de Gaulle) Airport'
    ]
    for n, loc_ in enumerate( locs2plot ):
        lon_, lat_, NIU = AC.get_loc(loc_)
        # If lat in plotted range, then plot
        if (lon_ > xlim[0]) and (lon_ < xlim[-1]):
    #        print(lat_)
            ax.axvline(x=lon_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=3.0)
            # Add label for airports
    #        if loc_ in ['RAI','DSS' ]:
            if n % 2 == 0:
                buffer = 0
            else:
                buffer = -0.5
            # Make a special case for Sao Vicente
            if loc_ == 'Sao Vicente Airport':
                buffer = -0.25

            # Set a base height for
            base = 9
            if convert2kft:
                base = base * m2kft
                buffer = buffer *3*1.5

            # Add label for the airport
            ax.text(lon_, base+buffer, '{}'.format(loc_), fontsize=10,
                    alpha=0.5,
                    horizontalalignment='center' )

    #        ax.annotate(loc_, xy=(lat_, 5), xytext=(lat_+buffer, 5+2),
    #            arrowprops=dict(facecolor='black', shrink=0.05))

    # Add lines for kft heights
    if convert2kft:
        hPa_heights = [1000, 900, 800, 700, 600, 500]
        km_heights = AC.hPa_to_Km(hPa_heights)
        kft_heights = [i*m2kft for i in km_heights]
        for n, height_ in enumerate( kft_heights ):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-3.5, height_, '{:.0f} hPa'.format(hPa_heights[n]),
                    fontsize=10,
                    alpha=0.5 )
    else:
        # Add lines for kft heights
        kft_heights = [20000, 15000, 10000, 5000]
        m_heights = [i/m2kft/1E3 for i in kft_heights]
        for n, height_ in enumerate( m_heights ):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-5, height_,
                    '{:.0f} kft'.format(kft_heights[n]/1E3),
                    fontsize=10, alpha=0.5 )

    # Beautify the figure/plot
    # Limit the yaxis
    ax.set_ylim(ylim)
    # Limit the xaxis
    ax.set_xlim(xlim)

    # Add hPA labels. - temporally remove .
#     press2plot = [1000, 900, 800, 700, 600, 500, 400]
#     ind = [ hPa_l.index(i) for i in press2plot ]
#     ax2.set_yticks(ds.lev[ind])
#     ax2.set_yticklabels([str(int(i)) for i in press2plot])
#     ax2.set_ylabel('Pressure [hPa]')

    # Add minor ticks for the x axis
#    ax.xaxis.grid(True, which='minor')

    # Update the lon axis labels
#    update_lon_units = True
#    if update_lon_units:
#        ax.tick_params(axis='x', which='both', labelbottom='off')
#        im.tick_params(axis='x', which='both', labelbottom='off')
#        pass

    # Add a generic title if one is not provided
    if not isinstance(title, type(None)):
#        ax.set_title(title)
        fig.suptitle(title)
    # Force tight layout
    plt.tight_layout(rect=(0, 0, 1, 0.95),)
#    plt.tight_layout()
    # Save the plot?
    if save_plot:
        if isinstance(savename, type(None)):
            sstr = 'spatial_plot_{}_{}_{}'
            savename = sstr.format(var2plot1, var2plot2, extra_str)
        savename = AC.rm_spaces_and_chars_from_str(savename)
        if isinstance(folder, type(None)):
            folder = './'
        plt.savefig(folder+savename+'.png', dpi=dpi)
    if show_plot:
        plt.show()

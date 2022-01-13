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
from shapely.geometry.polygon import LinearRing
import matplotlib
# Import local ARNA module functions
from . GEOS import *
from . GEOS_Chem import *
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
            nticks = 9-1
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
            nticks = 8-1
            ticks = np.linspace(vmin2, vmax2, nticks+1)
        cmap = get_reduced_cmap(cmap='Reds', maxval=0.75)
#        cmap = get_reduced_cmap(cmap='Greens', maxval=0.75)
#        cmap = get_reduced_cmap(cmap='Greys', maxval=0.75)
        cmap = AC.mk_discrete_cmap(nticks=nticks, cmap=cmap)
    elif input == 'NOy':
        ticks = None
        nticks = 6-1
        cmap = AC.mk_discrete_cmap(nticks=nticks, cmap='Blues')
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
                                    testing_mode=False,
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
        dts = AC.dt64_2_dt(ds.time.values)
        if testing_mode:
            dts = dts[:3]
        for n_time, t in enumerate(dts):
            # Sub select for time
            ds_tmp = ds[[var2plot]].sel(time=ds.time.values[n_time])
            # Create a date string to use to save file
            date_str = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
            date_str = date_str.format(t.year, t.month, t.day, t.hour,
                                       t.minute)
            title_date = t.strftime('%Y/%m/%d %H:%M')
            levs2use = ds.lev.values
            if testing_mode:
                levs2use = levs2use[:3]
            for lev2use in levs2use:
                print(var2plot, lev2use, title_date)
                ds_tmpII = ds_tmp[[var2plot]].sel(lev=lev2use)
                # Setup a string for the title
                title = '[{}] @ {:.0f}hPa on {}'
                title = title.format(LaTeX_spec, lev2use, title_date)
                if not isinstance(extr_title_str, type(None)):
                    title += '\n ' + extr_title_str
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
                                        save_plot=True)
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
        filename = 'ARNA_GEOSCF_chm_inst_1hr_g1440x721_p23_Cape_Verde_2019_'
        filename += '353_noy_lvls_1000_900_800_700_600_500.nc'
        ds = xr.open_dataset(folder + filename)
    # Local area analysed as Cape Verde
    d = get_analysis_region('local_CVAO_area')
    x0, x1, y0, y1 = d['x0'], d['x1'], d['y0'], d['y1']
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
        locs2plot = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
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
        filename = 'ARNA_GEOSCF_chm_inst_1hr_g1440x721_p23_Cape_Verde'
        filename += '_2019_353_noy_lvls_1000_900_800_700_600_500.nc'
        ds = xr.open_dataset(folder + filename)
    # Local area analysed as Cape Verde
    d = get_analysis_region('local_CVAO_area')
    x0, x1, y0, y1 = d['x0'], d['x1'], d['y0'], d['y1']
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
    lons2plot = [-18, -19.5, -21, -22.5, -24, -25.5]
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
        locs2plot = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
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


def plt_ts4ds(ds, region='Cape_Verde', extr_str='',
              vars2use=None, year=2018, verbose=False,
              show_plot=False, dpi=320, context="talk",
              font_scale=0.75):
    """
    Plot timeseries of data at different heights
    """
    # Which variables to include in analysis?
    if not isinstance(vars2use, list):
        vars2use = [i for i in ds.data_vars]
    # - Now plot up species as PDf based on level
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context(context, font_scale=font_scale)
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
            prtstr = 'WARNING: Not plotting {} @ {:.0f}hPa because all NaNs!'
            if mean_vals[np.isfinite(mean_vals)].shape[0] == 0:
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
                AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
                if show_plot:
                    plt.show()
            plt.close()
            del da
    # Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def PDF_on_species_in_ds4lvls(ds, region='Cape_Verde', extr_str='',
                              vars2use=None, year=2018, verbose=False,
                              show_plot=False, dpi=320,
                              context="talk"):
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
    sns.set_context(context)
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
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
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
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
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
        for lon2use in list(ds.lon.values):
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
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
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
        d = get_analysis_region('local_CVAO_area')
        extents = (d['x0'], d['x1'], d['y0'], d['y1'])
    # Mark known places to help geo-locate viewers
    if add_ARNA_locs:
        #        colours = AC.get_CB_color_cycle()
        locs2plot = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
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
                              verbose=False, testing=False):
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
        times2use = ds.time.values
    # Plot up by level and time
    for time2use in times2use:
        for lev2use in levs2use:
            # Get time as human readable string
            dstr = AC.dt64_2_dt([time2use])[0].strftime('%Y/%m/%d %H:%M')
            # Print out status
            pstr = "'plotting 2layer @ {:.0f}hPa on {}"
            print(pstr.format(lev2use, dstr))
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lev=lev2use, time=time2use)
            # Set title
            title = '[{}] & [{}] @ {:.0f}hPa on {}'.format(
                LaTeX_spec1, LaTeX_spec2, lev2use, dstr)
            # Add extra string to existing title string
            title += '\n ' + extr_title_str
            # Save plots
            extra_str = 'lev_{}_dt_{}'.format(lev2use, dstr)
            quick_map_plt_2layer(ds_tmp, var2plot1=var2plot1,
                                 folder=folder, region=region,
                                 var2plot2=var2plot2, title=title,
                                 add_max_vals_as_txt=add_max_vals_as_txt,
                                 save_plot=True, extra_str=extra_str)
            # Tidy up...
            plt.close('all')


def set_values_below_range2NaNs4spec(var=None, ds=None):
    """
    To improve aesthetics of plots, values below a threshold are removed
    """
    # Limit plotted NOy values to those above 0.5 pptv
    if var == 'NOy':
        arr = ds[var].values
        arr[arr < 0.5] = np.NaN
        ds[var].values = arr
    # Limit Dust values to those about
    elif 'Dust' in var:
        arr = ds[var].values
        arr[arr < 15] = np.NaN
        ds[var].values = arr
    else:
        pstr = "WARNING: No case set for '{}', so not restricting array values"
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
            print(pstr.format(lon2use*-1, '{\circ}', dstr))
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lon=lon2use, time=time2use)
            # Set title
            title_str = '[{}] & [{}] @ {:.1f}$^{}$W on {}'
            title = title_str.format(LaTeX_spec1, LaTeX_spec2, lon2use*-1,
                                     '{\circ}', dstr)
            if not isinstance(extr_title_str, type(None)):
                title += extr_title_str
            # Save plots
            extra_str = 'lon_{}E_dt_{}'.format(lon2use, dstr)
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
            print(pstr.format(lat2use, '{\circ}', dstr))
            # Select for level and variable, and average over time
            ds_tmp = ds.sel(lat=lat2use, time=time2use)
            # Set title
            title_str = '[{}] & [{}] @ {:.1f}$^{}$N on {}'
            title = title_str.format(LaTeX_spec1, LaTeX_spec2, lat2use,
                                     '{\circ}', dstr)
            if not isinstance(extr_title_str, type(None)):
                title += extr_title_str
            # Save plots
            extra_str = 'lat_{}N_dt_{}'.format(lat2use, dstr)
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
#    d = get_analysis_region('OPeNDAP_download_area')
    # Local area analysed as Cape Verde
    d = get_analysis_region('local_CVAO_area')
    x0, x1, y0, y1 = d['x0'], d['x1'], d['y0'], d['y1']
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
    # Force global perspective
    ax.set_global()  # this will force a global perspective
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
        cbar_kwargs = {'cmap': cmap, 'extend': extend, }
    else:
        cbar_kwargs = {'ticks': ticks, 'cmap': cmap, 'extend': extend, }
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
        cbar_kwargs = {'cmap': cmap, 'extend': extend, }
    else:
        cbar_kwargs = {'ticks': ticks, 'cmap': cmap, 'extend': extend, }
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
        d = get_analysis_region(region)  # 'local_CVAO_area'
        extents = (d['x0'], d['x1'], d['y0'], d['y1'])
    elif (region == 'Cape_Verde_Flying'):
        d = get_analysis_region(region)
        extents = (d['x0'], d['x1'], d['y0'], d['y1'])
    # Add extra lat and lon grid libnes
    if (region == 'Cape_Verde_Flying'):
        # Which X tickst to use?
        # x axis
        xticks = np.arange(ds.lon.values.min(), ds.lon.values.max(), 0.2)
        # y axis
        yticks = np.arange(ds.lat.values.min(), ds.lat.values.max(), 0.2)
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
        gl.xlocator = mticker.FixedLocator(xticks)  # last working setting...
#        ymajor = np.arange(ds.lat.values.min(), ds.lat.values.max(), 1 )
#        yminor = [i for i in yticks if i not in ymajor]
        gl.ylocator = mticker.FixedLocator(yticks)  # last working setting...
#        get_labels
#        ymajor_locator = mticker.FixedLocator(ymajor)
#        yminor_locator = mticker.FixedLocator(yminor)
#        ax.yaxis.set_major_locator(ymajor_locator)
#        ax.yaxis.set_minor_locator(yminor_locator)
        # tight off the main labels.
#        gl.xlabels_bottom = False
#        gl.ylabels_left = False
        gl.xlabels_bottom = True  # last working setting...
        gl.ylabels_left = True  # last working setting...
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
        gl.xlabel_style = {'size': 6, 'rotation': 90}
        gl.ylabel_style = {'size': 6, }

    # Mark a known place to help us geo-locate ourselves
    if add_ARNA_locs:
        colours = AC.get_CB_color_cycle()
        locs2plot = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
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
                                         transform=projection,
                                         facecolor='none',
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
            da = ds_tmp['NOy*Dust']
            print(da.where(da == da.max(), drop=True).squeeze())
            lon = da.where(da == da.max(), drop=True).squeeze().lon.values
            lat = da.where(da == da.max(), drop=True).squeeze().lat.values
            # Add a cross on the map.
            radius = (21 - 16.8331) / 4
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
        cbar_kwargs = {'cmap': cmap, 'extend': extend, 'pad': 0.075, }
    else:
        cbar_kwargs = {
            'ticks': ticks, 'cmap': cmap, 'extend': extend, 'pad': 0.075,
        }
    # Now plot up var1 - using pcolormesh
    cbar_ax = divider.append_axes("right", "2%", pad="1%")
    # Now plot
    ds[var2plot1].plot.pcolormesh(x=LatVar, y=LevVar, ax=ax,
                                  vmin=vmin1, vmax=vmax1,
                                  zorder=1, alpha=alpha,
                                  yincrease=True,
                                  cmap=cmap,
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
        cbar_kwargs = {'cmap': cmap, 'extend': 'both', }
    else:
        cbar_kwargs = {'ticks': ticks, 'cmap': cmap, 'extend': 'both', }
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
        max_alt = d['max_alt'] / 1E3
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
               xmin=(min_lat - xlim[0]) / xrange,
               xmax=(max_lat - xlim[0]) / xrange,
               linewidth=3.0)

    # Add locations for airports
    locs2plot = [
        'Praia Airport', 'Dakar',  'Gran Canaria Airport',
        'Sao Vicente Airport',
        'Lisbon Airport',  'Paris (Charles de Gaulle) Airport'
    ]
    for n, loc_ in enumerate(locs2plot):
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
                buffer = buffer * 3*1.5
            # Add label for the airport
            ax.text(lat_, base+buffer, '{}'.format(loc_), fontsize=10,
                    alpha=0.5,
                    horizontalalignment='center')

    # Add lines for kft heights
    if convert2kft:
        hPa_heights = [1000, 900, 800, 700, 600, 500]
        km_heights = AC.hPa_to_Km(hPa_heights)
        kft_heights = [i*m2kft for i in km_heights]
        for n, height_ in enumerate(kft_heights):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-2.5, height_, '{:.0f} hPa'.format(hPa_heights[n]),
                    fontsize=10,
                    alpha=0.5)
    else:
        # Add lines for kft heights
        kft_heights = [20000, 15000, 10000, 5000]
        m_heights = [i/m2kft/1E3 for i in kft_heights]
        for n, height_ in enumerate(m_heights):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-5, height_,
                    '{:.0f} kft'.format(kft_heights[n]/1E3),
                    fontsize=10, alpha=0.5)

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


def plt_comp_by_alt_4ARNA_all(dpi=320, just_SLR=True, show_plot=False,
                              context="paper", font_scale=0.75,
                              just_plot_GEOS_Chem=False,
                              inc_GEOSChem=False,
                              res='4x5', RunSet=None,
                              flight_nums=[],
                              savetitle=None,
                              pdff=None,
                              close_pdf=False,
                              PltPointObs=False,
                              JustPlotModel=False,
                              CoreRunsOnly=False,
                              PltSpatialLocsOfdata=True,
                              vars2plot=None,
                              NOxAsLog=False,
                              debug=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    import seaborn as sns
    # Which flights to plot? - Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [218, 219, 220, 221, 222, 223, 224, 225, ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_FAAM_core4flightnum(flight_ID=flight_ID)
#        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_obs[flight_ID] = df
    # Get the point observations by flight if these are being plotted
    if PltPointObs:
        dfP = get_filters_data4flight(all_flights=True)
        dfP = dfP.loc[dfP['Flight'].isin(flight_IDs), :]
    # Model - GEOS-CF (online)
    dfs_mod_CF = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF4flightnum(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_mod_CF[flight_ID] = df
    gc.collect()
    # Model - GEOS-Chem (offline)
    if inc_GEOSChem:
        dfs_mod_GC = {}
        for flight_ID in flight_IDs:
            dfs = get_GEOSChem4flightnum(flight_ID=flight_ID,
                                         res=res,
                                         CoreRunsOnly=CoreRunsOnly,
                                         RunSet=RunSet,)
            for key in dfs.keys():
                df = dfs[key]
                df = add_derived_FAAM_flags2df4flight(df=df,
                                                      flight_ID=flight_ID)
                df['flight_ID'] = flight_ID

                dfs[key] = df
            dfs_mod_GC[flight_ID] = dfs
    gc.collect()
#            del dfs
    # TODO: Combine to a single DataFrame to plot GEOS-CF and GEOS-Chem
    # Kludge - for now just plot GEOS-Chem
    if just_plot_GEOS_Chem:
        dfs_mod = dfs_mod_GC
        print('TEMP DEBUG:', dfs_mod_GC.keys(), dfs_mod.keys(), dfs.keys())
    else:
        print('TODO: setup plotting of GEOS-CF and GEOS-Chem')
        dfs_mod = dfs_mod_CF
    # Which 'RunSet' to use
    if RunSet == 'FP-Nest':
        run2use = RunSet
    else:
        run2use = list(dfs_mod[list(dfs_mod.keys())[0]].keys())[0]
    if debug:
        print(dfs_mod.keys())
        print(list(dfs_mod.keys())[0])
        print(dfs_mod[list(dfs_mod.keys())[0]])
        print(dfs_mod[list(dfs_mod.keys())[0]][run2use])
        print(dfs_mod[list(dfs_mod.keys())[0]][run2use].head())
    # Combine to a single dataframe (dictionary lists are by flight )
    df_mod = pd.concat([dfs_mod[i][run2use] for i in dfs_mod.keys()], axis=0)

    dfs_mod_ALL = {}
    if len(dfs.keys()) > 1:
        for key in dfs.keys():
            ModByFlight = [dfs_mod[i][key] for i in dfs_mod.keys()]
            dfs_mod_ALL[key] = pd.concat(ModByFlight, axis=0)
    if debug:
        print('TEMP DEBUG:', df_mod)
    df_obs = pd.concat([dfs_obs[i] for i in dfs_obs.keys()], axis=0)
    # Only consider data during SLRs?
    if just_SLR:
        df_obs = df_obs.loc[df_obs['IS_SLR'] == True, :]
        df_mod = df_mod.loc[df_mod['IS_SLR'] == True, :]
        extr_str = '_JUST_SLR'
    else:
        extr_str = ''
    # Setup PDF to save PDF plots to
    if isinstance(savetitle, type(None)):
        savetitle = 'ARNA_altitude_binned_{}{}'.format('ALL', extr_str)
    if isinstance(pdff, type(None)):
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # - Plot up location of flights
    if just_SLR:
        title = "Flight tracks for 'Straight and Level' Runs during ARNA"
    else:
        title = 'Flight tracks for all flights during ARNA'
    if PltSpatialLocsOfdata:
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID,
                                           title=title)
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

    # - Put observations and vars to plot into a dictionary
    sns.set(color_codes=True)
    sns.set_context(context, font_scale=font_scale)
    # Force alt to be in units of km
    ALT_var = 'Altitude (km)'
    Y_unit = ALT_var
    key4GEOSCF = 'GEOS-CF'
    print('WARNING: Kludged below to only plot GEOS-Chem alt var')
    if just_plot_GEOS_Chem:
        df_mod[ALT_var] = AC.hPa_to_Km(df_mod['PRESS'].values)
        if len(dfs.keys()) > 1:
            for key in dfs_mod_ALL.keys():
                df_tmp = dfs_mod_ALL[key]
                df_tmp[ALT_var] = AC.hPa_to_Km(df_tmp['PRESS'].values)
                dfs_mod_ALL[key] = df_tmp
    else:
        df_mod[ALT_var] = AC.hPa_to_Km(df_mod['model-lev'].values)
    df_obs[ALT_var] = df_obs['ALT_GIN'].values / 1E3
    #
    if just_plot_GEOS_Chem:
        if (len(list(dfs_mod.keys())) > 1):
            # If more than one set of model output provided?
            data_d = AC.merge_two_dicts({'Obs.': df_obs}, dfs_mod_ALL)
        else:
            data_d = {RunSet: df_mod, 'Obs.': df_obs}
    else:
        data_d = {'GEOS-CF': df_mod, 'Obs.': df_obs}

    # - Now plot up flight time series plots by variable
    if just_SLR:
        title_str = "Altitude binned FAAM Core '{}' ({})"
        title_str += "for all 'Straight+Level Runs'"
    else:
        title_str = "Altitude binned '{}' ({}) for all flights"
    # Setup color dictionary for plotting...
    colors2use = AC.get_CB_color_cycle()
    color_dict = {
        'GEOS-CF': 'red',
        'Obs.': 'k',
        RunSet: colors2use[0],
        'FP-Nest-JNITx25': colors2use[1],
        'FP-Nest-BBx2': colors2use[2],
    }
    runs2color = [i for i in data_d.keys() if i not in color_dict.keys()]
    for n_run, run in enumerate(runs2color):
        color_dict[run] = colors2use[n_run]
    # And conversion scales and units for variables
    unit_d = {}
    mod2obs_varnames = {
        'CO': 'CO_AERO', 'O3': 'O3_TECO', 'NO2': 'NO2_pptV', 'NO': 'NO_pptV',
        'HNO2': 'HONO_pptV',
        'NOx': 'NOx'
    }
    units_d = {
        'CO': 'ppbv', 'O3': 'ppbv', 'NO2': 'pptv', 'NO': 'pptv', 'NOx': 'pptv',
        'HNO2': 'pptv', 'HONO': 'pptv',
        'NIT': 'pptv', 'NITs': 'pptv', 'SO4s': 'pptv', 'SO4': 'pptv',
        'NH4': 'pptv',
        'SO4-all': 'pptv', 'NIT-all': 'pptv',
    }
    range_d = {
        'CO': (50, 400), 'O3': (-10, 100), 'NO2': (-50, 500), 'NO': (-50, 500),
        'NOx': (-50, 500),
        'HNO2': (-60, 60), 'HONO': (-60, 60),
    }
    NOx_specs = ['HNO2', 'NOx', 'NO', 'NO2', 'HONO']
    # - by variable
    runs = list(sorted(data_d.keys()))
    if debug:
        print(runs)
    if JustPlotModel:
        runs.pop(runs.index('Obs.'))
    # Which variables to use?
    if isinstance(vars2plot, type(None)):
        vars2plot = list(sorted(mod2obs_varnames.keys()))[::-1]
        vars2plot = ['CO', 'O3', 'NOx', 'NO2', 'NO', 'HNO2']
        vars2plot = [
            i for i in vars2plot if mod2obs_varnames[i] in df_obs.columns
        ]
    if debug:
        print(vars2plot)
        print(df_obs.columns)

    # What bins should be used?
    bins = [0.5*i for i in np.arange(15)]
    for var2plot in vars2plot:
        fig = plt.figure()
        ax = plt.gca()
        # Now loop data
        for n_key, key_ in enumerate(runs):
            print(n_key, key_, var2plot)
            #
            if key_ == 'Obs.':
                varname = mod2obs_varnames[var2plot]
            else:
                varname = var2plot
            # Setup an axis label
            units = units_d[var2plot]
            xlabel = '{} ({})'.format(var2plot, units)
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
            print(df.head())
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
            # Beautify plot
            plt.legend()
            plt.title(title_str.format(var2plot, units, flight_ID))
            # Make NOx species be on a log scale
            xscale = 'linear'
            if (var2plot in NOx_specs) and NOxAsLog:
                xscale = 'log'
            ax.set_xscale(xscale)
            if xscale == 'log':
                if var2plot in ['HONO', 'HNO2']:
                    xlim = (0.03, 400)
                else:
                    xlim = (0.3, 400)
                ax.set_xlim(xlim)
            else:
                try:
                    plt.xlim(range_d[var2plot])
                except KeyError:
                    pass

        # Add point / time limited observations?
        if PltPointObs:
            VarStr = 'Total_{}_ppt'
            if 'so4' in var2plot.lower():
                PointVar = VarStr.format('SO4')
            elif 'nit' in var2plot.lower():
                PointVar = VarStr.format('NO3')
            elif 'nh4' in var2plot.lower():
                PointVar = VarStr.format('NH4')
            else:
                PrtStr = "WARING: Plotting not setup for '{}' filter data"
                print(PrtStr.format(var2plot))
            AltVar = 'Average_altitude_m'
            plt.scatter(dfP[PointVar].values, dfP[AltVar].values/1E3,
                        s=25, color='K', label='Filters',
                        alpha=0.75, zorder=20)
            # TODO - add plotting of error bars for the observation here.
            plt.legend()

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        if show_plot:
            plt.show()
        plt.close()

    # - Save entire pdf
    if close_pdf:
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plt_comp_by_alt_4ARNA_all_DUST(dpi=320, just_SLR=True, flight_nums=[],
                                   plt_model=False, show_plot=False,
                                   context="paper", font_scale=0.75,
                                   NOxAsLog=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    import seaborn as sns
    # Which flights to plot? - Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [218, 219, 220, 221, 222, 223, 224, 225, ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_FAAM_core4flightnum(flight_ID=flight_ID)
#        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_obs[flight_ID] = df
    gc.collect()
    # Model
    dfs_mod = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF4flightnum(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_mod[flight_ID] = df
    gc.collect()
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
    df_mod[ALT_var] = AC.hPa_to_Km(df_mod['model-lev'].values)
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
    AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
    plt.close()

    # - Put observations and vars to plot into a dictionary
    sns.set(color_codes=True)
    sns.set_context(context, font_scale=font_scale)

    # - Now plot up flight time series plots by variable
    if just_SLR:
        title_str = "Altitude binned '{}' ({}) for all 'Straight+Level Runs'"
    else:
        title_str = "Altitude binned '{}' ({}) for all flights"
    # Setup color dictinoary
    color_dict = {'GEOS-CF': 'red', 'Obs.': 'k'}
    colors2use = AC.get_CB_color_cycle()
    runs2color = [i for i in data_d.keys() if i not in color_dict.keys()]
    for n_run, run in enumerate(runs2color):
        color_dict[run] = colors2use[n_run]
    unit_d = {}
    mod2obs_varnames = {
        'CO': 'CO_AERO', 'O3': 'O3_TECO', 'NO2': 'NO2_pptV', 'NO': 'NO_pptV',
        'HNO2': 'HONO_pptV',
        'NOx': 'NOx'
    }
    units_d = {
        'CO': 'ppbv', 'O3': 'ppbv', 'NO2': 'pptv', 'NO': 'pptv', 'NOx': 'pptv',
        'HNO2': 'pptv', 'HONO': 'pptv',
    }
    range_d = {
        'CO': (50, 400), 'O3': (-10, 100), 'NO2': (-50, 500), 'NO': (-50, 500),
        'NOx': (-50, 500),
        'HNO2': (-60, 60), 'HONO': (-60, 60),
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
            print(n_key, key_, var2plot)
            #
            if ('Obs.' in key_):
                varname = mod2obs_varnames[var2plot]
            else:
                varname = var2plot
            # Setup an axis label
            units = units_d[var2plot]
            xlabel = '{} ({})'.format(var2plot, units)
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

            # Beautify plot
            plt.legend()
            plt.title(title_str.format(var2plot, units, flight_ID))
            # Make NOx species be on a log scale
            xscale = 'linear'
            if (var2plot in NOx_specs) and NOxAsLog:
                xscale = 'log'
            ax.set_xscale(xscale)
            if xscale == 'log':
                if var2plot in ['HONO', 'HNO2']:
                    xlim = (0.03, 400)
                else:
                    xlim = (0.3, 400)
                ax.set_xlim(xlim)
            else:
                try:
                    plt.xlim(range_d[var2plot])
                except KeyError:
                    pass

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        if show_plot:
            plt.show()
        plt.close()

    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')
    gc.collect()


def plt_comp_by_alt_4ARNA_CIMS_all_DUST(dpi=320, just_SLR=True,
                                        plt_model=False, show_plot=False,
                                        flight_nums=[],
                                        res='4x5', RunSet=None,
                                        context="paper", font_scale=0.75,
                                        savetitle=None, pdff=None,
                                        PltSpatialLocsOfdata=True,
                                        filter_by_dust=False,
                                        just_plot_GEOS_Chem=False,
                                        inc_GEOSChem=False,
                                        CoreRunsOnly=False,
                                        close_pdf=True,
                                        NOxAsLog=False,
                                        debug=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    import seaborn as sns
    # Which flights to plot? - Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [217, 218, 219, 220, 221, 222, 223, 224, 225, ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Setup dictionary of Observation (CIMS) dataframes
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_CIMS_data4flight(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_obs[flight_ID] = df
    # Setup dictionary of Model dataframes
    # Model - GEOS-CF (online)
    dfs_mod_CF = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF4flightnum(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_mod_CF[flight_ID] = df
    # Model - GEOS-Chem (offline)
    if inc_GEOSChem:
        dfs_mod_GC = {}
        for flight_ID in flight_IDs:
            dfs = get_GEOSChem4flightnum(flight_ID=flight_ID,
                                         res=res,
                                         CoreRunsOnly=CoreRunsOnly,
                                         RunSet=RunSet,)
            for key in dfs.keys():
                df = dfs[key]
                df = add_derived_FAAM_flags2df4flight(df=df,
                                                      flight_ID=flight_ID)
                df['flight_ID'] = flight_ID

                dfs[key] = df
            dfs_mod_GC[flight_ID] = dfs
#            del dfs
    gc.collect()

    # Kludge - for now just plot GEOS-Chem
    if just_plot_GEOS_Chem:
        dfs_mod = dfs_mod_GC
    else:
        print('TODO: setup plotting of GEOS-CF and GEOS-Chem')
        dfs_mod = dfs_mod_CF
    # Which 'RunSet' to use
    if RunSet == 'FP-Nest':
        run2use = RunSet
    else:
        run2use = list(dfs_mod[list(dfs_mod.keys())[0]].keys())[0]
    # Combine to a single dataframe
    df_mod = pd.concat([dfs_mod[i][run2use] for i in dfs_mod.keys()], axis=0)
    if debug:
        print('TEMP DEBUG:', df_mod)
    # Now setup a dictionary for all model runs
    dfs_mod_ALL = {}
    if len(dfs.keys()) > 1:
        for key in dfs.keys():
            ModByFlight = [dfs_mod[i][key] for i in dfs_mod.keys()]
            dfs_mod_ALL[key] = pd.concat(ModByFlight, axis=0)
    if debug:
        print('TEMP DEBUG:', df_mod)
    df_obs = pd.concat([dfs_obs[i] for i in dfs_obs.keys()], axis=0)

    # Only consider data during straight and level runs (SLRs)?
    if just_SLR:
        df_obs = df_obs.loc[df_obs['IS_SLR'] == True, :]
        df_mod = df_mod.loc[df_mod['IS_SLR'] == True, :]
        extr_str = '_JUST_SLR_DUST'
    else:
        extr_str = ''
    if filter_by_dust:
        extr_str += '_DUST'

    # Setup PDF to save PDF plots to
    if isinstance(savetitle, type(None)):
        savetitle = 'ARNA_altitude_binned_{}_CIMS_{}'.format('ALL', extr_str)
    if isinstance(pdff, type(None)):
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # - Setup data objects or plotting
    # Force alt to be in units of km
    ALT_var = 'Altitude (km)'
    Y_unit = ALT_var
    print('WARNING: Kludged below to only plot GEOS-Chem alt var')
    if just_plot_GEOS_Chem:
        #        df_mod[ALT_var] = AC.hPa_to_Km(df_mod['PRESS'].values)
        df_mod[ALT_var] = AC.hPa_to_Km(df_mod['PRESS'].values)
        if len(dfs.keys()) > 1:
            for key in dfs_mod_ALL.keys():
                df_tmp = dfs_mod_ALL[key]
                df_tmp[ALT_var] = AC.hPa_to_Km(df_tmp['PRESS'].values)
                dfs_mod_ALL[key] = df_tmp
    else:
        df_mod[ALT_var] = AC.hPa_to_Km(df_mod['model-lev'].values)
    df_obs[ALT_var] = df_obs['ALT_GIN'].values / 1E3
#    df_mod[ALT_var] = AC.hPa_to_Km(df_mod['model-lev'].values)
#    df_obs[ALT_var] = df_obs['ALT_GIN'].values / 1E3
    # Plot up just observations? Or model too?
    # Add filtering for Dust (NOTE: only setup for GEOS-CF plots)
    if filter_by_dust:
        data_d = {
            'Obs.':  df_obs.loc[df_obs['IS_DUST'] == False, :],
            'Obs. (Dust)': df_obs.loc[df_obs['IS_DUST'] == True, :],
        }
        if plt_model:
            data_d['GEOS-CF'] = df_mod.loc[df_mod['IS_DUST'] == False, :]
            data_d['GEOS-CF (dust)'] = df_mod.loc[df_mod['IS_DUST'] == True, :]
            extr_str += '_inc_MODEL'

    else:
        data_d = {'Obs.':  df_obs}
        if just_plot_GEOS_Chem:
            if (len(list(dfs_mod.keys())) > 1):
                # If more than one set of model output provided?
                data_d = AC.merge_two_dicts({'Obs.': df_obs}, dfs_mod_ALL)
            else:
                data_d = {RunSet: df_mod, 'Obs.': df_obs}
#            data_d = {RunSet: df_mod, 'Obs.': df_obs}
            # Add extra runs via dfs_mod dictionary
#            AC.merge_two_dicts(data_d, dfs_mod)
        else:
            data_d = {'GEOS-CF': df_mod, 'Obs.': df_obs}

    # - Plot up location of flights
    if just_SLR:
        title = "Flight tracks for 'Straight+Level Runs (SLRs)' during ARNA"
    else:
        title = 'Flight tracks for all flights during ARNA'
    if PltSpatialLocsOfdata:
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID,
                                           title=title)
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

    # - Put observations and vars to plot into a dictionary
    sns.set(color_codes=True)
    sns.set_context(context, font_scale=font_scale)

    # - Now plot up flight time series plots by variable
    if just_SLR:
        title_str = "Altitude binned CIMS '{}' ({})"
        title_str += "for all 'Straight+Level Runs'"
    else:
        title_str = "Altitude binned '{}' ({}) for all flights"
    # Setup color dictinoary
    colors2use = AC.get_CB_color_cycle()
    color_dict = {
        'GEOS-CF': 'red',
        'Obs.': 'k',
        RunSet: colors2use[0],
        'FP-Nest-JNITx25': colors2use[1],
        'FP-Nest-BBx2': colors2use[2],
    }
    runs2color = [i for i in data_d.keys() if i not in color_dict.keys()]
    for n_run, run in enumerate(runs2color):
        color_dict[run] = colors2use[n_run]
    NOx_specs = ['HNO2', 'NOx', 'NO', 'NO2', 'HONO']
    mod2obs_varnames = {
        'BrO': 'BrO',
        'HNO3': 'HNO3',
        'HNO2': 'HONO',
    }
    units_d = {
        'BrO': 'pptv', 'HNO3': 'pptv', 'HNO2': 'pptv', 'HONO': 'pptv',
    }
    range_d = {
        'HNO2': (-10, 60),
        'HNO3': (-30, 1500),
        'BrO': (-0.2, 1.0),
        'HONO': (-10, 60),
    }
    # - by variable
    print(data_d.keys())
    runs = list(data_d.keys())
    # Which variables to use?
    vars2plot = list(sorted(mod2obs_varnames.keys()))[::-1]
#    vars2plot = ['CO', 'O3', 'NOx', 'NO2', 'NO', 'HNO2']
    if debug:
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
            print(n_key, key_, var2plot)
            if ('Obs.' in key_):
                varname = mod2obs_varnames[var2plot]
            else:
                varname = var2plot
            # Setup an axis label
            units = units_d[var2plot]
            xlabel = '{} ({})'.format(var2plot, units)
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
            xscale = 'linear'
            if (var2plot in NOx_specs) and NOxAsLog:
                xscale = 'log'
            ax.set_xscale(xscale)
            # Beautify plot
            plt.legend()
            plt.title(title_str.format(var2plot, units, flight_ID))
            if xscale == 'log':
                if var2plot in ['HONO', 'HNO2']:
                    xlim = (0.03, 400)
                else:
                    xlim = (0.3, 400)
                ax.set_xlim(xlim)
            else:
                try:
                    plt.xlim(range_d[var2plot])
                except KeyError:
                    pass

        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        if show_plot:
            plt.show()
        plt.close()

    # - Save entire pdf
    if close_pdf:
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plt_comp_by_alt_4ARNA_flights(dpi=320, just_SLR=True, show_plot=False,
                                  RunSet=None, res='4x5', flight_nums=[],
                                  just_plot_GEOS_Chem=False,
                                  inc_GEOSChem=False,
                                  context="paper", font_scale=0.75,
                                  NOxAsLog=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    import seaborn as sns
    # Which flights to plot?
    # Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [
            # Just use non-transit ARNA flights
            #        216,
            #        217,
            218, 219, 220, 221, 222, 223, 224, 225,
        ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_FAAM_core4flightnum(flight_ID=flight_ID)
#        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        dfs_obs[flight_ID] = df
    # Observations
#    dfs_mod = {}
#    for flight_ID in flight_IDs:
#        df = get_GEOSCF4flightnum(flight_ID=flight_ID)
#        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
#        dfs_mod[flight_ID] = df

    # Model - GEOS-CF (online)
    dfs_mod_CF = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF4flightnum(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        dfs_mod_CF[flight_ID] = df
    # Model - GEOS-Chem (offline)
    if inc_GEOSChem:
        #        RunSet='MERRA2-0.5-initial'
        #        res='0.5x0.625'
        #        RunSet='MERRA2-BC'
        #        res='4x5'
        #        RunSet='FP-Nest'
        #        res='0.25x0.3125'
        dfs_mod_GC = {}
        for flight_ID in flight_IDs:
            dfs = get_GEOSChem4flightnum(flight_ID=flight_ID,
                                         res=res,
                                         RunSet=RunSet,)
            for key in dfs.keys():
                df = dfs[key]
                df = add_derived_FAAM_flags2df4flight(df=df,
                                                      flight_ID=flight_ID)
                dfs[key] = df
            dfs_mod_GC[flight_ID] = dfs
            del dfs

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        df_obs = dfs_obs[flight_ID]
        df_mod_CF = dfs_mod_CF[flight_ID]
        if inc_GEOSChem:
            if just_plot_GEOS_Chem:
                dfs_mod = dfs_mod_GC[flight_ID]
                mod_label_master = RunSet

            else:
                dfs_mod_GC4flight = dfs_mod_GC[flight_ID]
                dfs_mod = {'GEOS-CF': df_mod_CF}
                for key in list(dfs_mod_GC4flight.keys()):
                    dfs_mod[key] = dfs_mod_GC4flight[key]
                mod_label_master = 'GEOS-CF'

        else:
            mod_label_master = 'GEOS-CF'
            dfs_mod = {'GEOS-CF': df_mod_CF}
        print(('!'*30, mod_label_master, dfs_mod.keys()))

        # Get observations and model timeseries data as a DataFrame
#        df_obs = dfs_obs[flight_ID]
#        df_mod = dfs_mod[flight_ID]
        # Only consider data during SLRs?
        if just_SLR:
            df_obs = df_obs.loc[df_obs['IS_SLR'] == True, :]
            for key in dfs_mod.keys():
                df_mod = dfs_mod[key]
                df_mod = df_mod.loc[df_mod['IS_SLR'] == True, :]
                dfs_mod[key] = df_mod
            extr_str = '_JUST_SLR'
        else:
            extr_str = ''
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_altitude_binned_{}'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID)
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Put observations and vars to plot into a dictionary
        sns.set(color_codes=True)
        sns.set_context(context, font_scale=font_scale)
        # Force alt to be in units of km
        ALT_var = 'Altitude (km)'
        Y_unit = ALT_var
        key4GEOSCF = 'GEOS-CF'
        for key in dfs_mod.keys():
            #        if key in dfs_mod.keys():
            df_mod = dfs_mod[key]
            if key4GEOSCF == key:
                df_mod[ALT_var] = AC.hPa_to_Km(df_mod['model-lev'].values)
            else:
                df_mod[ALT_var] = AC.hPa_to_Km(df_mod['PRESS'].values)
            dfs_mod[key] = df_mod
        df_obs[ALT_var] = df_obs['ALT_GIN'].values / 1E3
        # Online consider offline model output?
        if just_plot_GEOS_Chem:
            #            GEOSChem_varname = sorted(list(dfs_mod.keys())
            GCvarname = 'FP-Nest'
            data_d = {GCvarname: dfs_mod[GCvarname], 'Obs.': df_obs}
        else:
            data_d = AC.merge_two_dicts(dfs_mod, {'Obs.': df_obs})
        print('Plotting model runs: ', data_d.keys())
        # - Now plot up flight time series plots by variable
        title_str = "Altitude binned '{}' ({}) during flight '{}'"
        # Setup color dictionary
        color_dict = {'GEOS-CF': 'red', 'Obs.': 'k'}
        CB_color_cycle = AC.get_CB_color_cycle()
        for n_key, key in enumerate(list(data_d.keys())):
            if key not in color_dict.keys():
                color_dict[key] = CB_color_cycle[n_key]
        unit_d = {}
        mod2obs_varnames = {
            'CO': 'CO_AERO', 'O3': 'O3_TECO', 'NO2': 'NO2_pptV',
            'NO': 'NO_pptV', 'HNO2': 'HONO_pptV', 'NOx': 'NOx'
        }
        units_d = {
            'CO': 'ppbv', 'O3': 'ppbv', 'NO2': 'pptv', 'NO': 'pptv', 'NOx': 'pptv',
            'HNO2': 'pptv', 'HONO': 'pptv',
        }
        range_d = {
            'CO': (50, 400), 'O3': (-10, 100), 'NO2': (-50, 500), 'NO': (-50, 500),
            'NOx': (-50, 500),
            'HNO2': (-60, 60), 'HONO': (-60, 60),
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
        print('Plotting:', runs)
        bins = [0.5*i for i in np.arange(15)]
        for var2plot in vars2plot:
            fig = plt.figure()
            ax = plt.gca()
            # Now loop data
            for n_key, key_ in enumerate(runs):
                print(n_key, key_, var2plot)
                #
                if key_ == 'Obs.':
                    varname = mod2obs_varnames[var2plot]
                else:
                    varname = var2plot
                # Setup an axis label
                units = units_d[var2plot]
                xlabel = '{} ({})'.format(var2plot, units)
                # Add alt to DataFrame
                df = pd.DataFrame({
                    var2plot: data_d[key_][varname],
                    ALT_var: data_d[key_][ALT_var]
                })
                # Scale the modelled values to the same units
                if key_ != 'Obs.':
                    scaleby = AC.get_unit_scaling(units)
                    df[var2plot] = df[var2plot].values * scaleby
                print(df.head())
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
                                                   widths=0.15,
                                                   dataset_num=n_key,
                                                   color=color_dict[key_])
                except:
                    pass
                # Beautify plot
                plt.legend()
                plt.title(title_str.format(var2plot, units, flight_ID))
                # Make NOx species be on a log scale
                xscale = 'linear'
                if (var2plot in NOx_specs) and NOxAsLog:
                    xscale = 'log'
                ax.set_xscale(xscale)
                if xscale == 'log':
                    if var2plot in ['HONO', 'HNO2']:
                        xlim = (0.03, 400)
                    else:
                        xlim = (0.3, 400)
                    ax.set_xlim(xlim)
                else:
                    try:
                        plt.xlim(range_d[var2plot])
                    except KeyError:
                        pass

            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            if show_plot:
                plt.show()
            plt.close()

        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_comp_by_alt_4ARNA_together(dpi=320, just_SLR=True, show_plot=False,
                                   RunSet='FP-Nest', res='0.25x3125',
                                   flight_nums=[], savetitle=None,
                                   just_plot_GEOS_Chem=False,
                                   inc_GEOSChem=False,
                                   context="paper", font_scale=0.75,
                                   NOxAsLog=False,
                                   CoreRunsOnly=False,
                                   verbose=True, debug=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    PrtStr = "NOTE: Plotting RunSet: '{}' @ res. of: '{}' (CoreRunsOnly?: {})"
    print( PrtStr.format(RunSet, res, CoreRunsOnly) )
    # Setup the pdf file to use
    if isinstance(savetitle, type(None)):
        savetitle = 'ARNA_altitude_binned_combined_file_{}'.format(res)
        if just_SLR:
            savetitle += '_JUST_SLR'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

#    if isinstance(RunSet, type(None)):
#        RunSet='FP-Nest'
    # Call the standard species plotter
    plt_comp_by_alt_4ARNA_all(just_SLR=just_SLR, context=context,
                              RunSet=RunSet, res=res,
                              inc_GEOSChem=True,
                              just_plot_GEOS_Chem=True,
                              CoreRunsOnly=CoreRunsOnly,
                              savetitle=savetitle,
                              flight_nums=flight_nums,
                              pdff=pdff,
                              close_pdf=False,
                              PltPointObs=False,
                              NOxAsLog=NOxAsLog,
                              debug=debug,
                              )

    # Now plot the aerosol species
    vars2plot = ['NIT-all', 'SO4-all', 'NH4']
    if 'acid' in str(RunSet).lower():
        vars2plot = ['NIT-all', 'NH4']
    plt_comp_by_alt_4ARNA_all(just_SLR=just_SLR, context=context,
                              RunSet=RunSet, res=res,
                              inc_GEOSChem=True,
                              just_plot_GEOS_Chem=True,
                              CoreRunsOnly=CoreRunsOnly,
                              savetitle=savetitle,
                              flight_nums=flight_nums,
                              pdff=pdff,
                              close_pdf=False,
                              PltPointObs=True,
                              JustPlotModel=True,
                              PltSpatialLocsOfdata=False,
                              vars2plot=vars2plot,
                              debug=debug,
                              )

    # Call the CIMS plotter
    # NOTE the HNO3 data is not final.
    plt_comp_by_alt_4ARNA_CIMS_all_DUST(context=context,
                                        inc_GEOSChem=True,
                                        just_plot_GEOS_Chem=True,
                                        CoreRunsOnly=CoreRunsOnly,
                                        plt_model=True,
                                        RunSet=RunSet, res=res,
                                        savetitle=savetitle,
                                        flight_nums=flight_nums,
                                        pdff=pdff,
                                        PltSpatialLocsOfdata=False,
                                        close_pdf=False,
                                        NOxAsLog=NOxAsLog,
                                        debug=debug,
                                        )
    # Add total NOy (inc. filters)
    vars2plot = ['NOy', 'NOy-gas']
    plt_comp_by_alt_4ARNA_NOy(vars2plot=vars2plot,
                              flight_nums=flight_nums,
                              context=context,
                              RunSet=RunSet, res=res,
                              inc_GEOSChem=True,
                              inc_GEOSCF=False,
                              CoreRunsOnly=CoreRunsOnly,
                              savetitle=savetitle,
                              pdff=pdff,
                              close_pdf=False,
                              show_plot=False,
                              debug=debug,
                              )

    # And SWAS ...

    # Add any other plotter calls here  ...

    # - Save entire pdf
    AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plt_comp_by_alt_4ARNA_flights_CIMS(dpi=320, just_SLR=False,
                                       show_plot=False,
                                       RunSet=None, res='4x5', flight_nums=[],
                                       just_plot_GEOS_Chem=False,
                                       inc_GEOSChem=False,
                                       context="paper", font_scale=0.75,
                                       NOxAsLog=False):
    """
    Plot up altitude binned comparisons between core obs. and model data
    """
    # Which flights to plot?
    # Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [
            #    217, # Missing data for C217 (NOy)
            218, 219, 220,
            #    221, # Missing data for C221 (NOy)
            222, 223,
            #    224,  # Missing data for C221 (BrO... )
            #    225,  # Missing data for C221 (BrO... )
        ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]

    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        df = get_CIMS_data4flight(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        dfs_obs[flight_ID] = df

    # Model
    dfs_mod_CF = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF4flightnum(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        dfs_mod_CF[flight_ID] = df

    # Model - GEOS-Chem (offline)
    if inc_GEOSChem:
        #        RunSet='MERRA2-0.5-initial'
        #        res='0.5x0.625'
        #        RunSet='MERRA2-BC'
        #        res='4x5'
        #        RunSet='FP-Nest'
        #        res='0.25x0.3125'
        dfs_mod_GC = {}
        for flight_ID in flight_IDs:
            dfs = get_GEOSChem4flightnum(flight_ID=flight_ID, res=res,
                                         RunSet=RunSet,)
            for key in dfs.keys():
                df = dfs[key]
                df = add_derived_FAAM_flags2df4flight(df=df,
                                                      flight_ID=flight_ID)
                dfs[key] = df
            dfs_mod_GC[flight_ID] = dfs
            del dfs

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
#        df_mod = dfs_mod[flight_ID]
        df_mod_CF = dfs_mod_CF[flight_ID]
        if inc_GEOSChem:
            if just_plot_GEOS_Chem:
                dfs_mod = dfs_mod_GC[flight_ID]
                mod_label_master = RunSet

            else:
                dfs_mod_GC4flight = dfs_mod_GC[flight_ID]
                dfs_mod = {'GEOS-CF': df_mod_CF}
                for key in list(dfs_mod_GC4flight.keys()):
                    dfs_mod[key] = dfs_mod_GC4flight[key]
                mod_label_master = 'GEOS-CF'

        else:
            mod_label_master = 'GEOS-CF'
            dfs_mod = {'GEOS-CF': df_mod_CF}
        print(('!'*30, mod_label_master, dfs_mod.keys()))

        # Only consider data during SLRs?
        if just_SLR:
            df_obs = df_obs.loc[df_obs['IS_SLR'] == True, :]
#            df_mod = df_mod.loc[df_mod['IS_SLR'] == True, :]
            for key in dfs_mod.keys():
                df_mod = dfs_mod[key]
                df_mod = df_mod.loc[df_mod['IS_SLR'] == True, :]
                dfs_mod[key] = df_mod
            extr_str = '_JUST_SLR'
        else:
            extr_str = ''
        # Setup PDF to save PDF plots to
        savetitle_str = 'ARNA_altitude_binned_{}_CIMS{}'
        savetitle = savetitle_str.format(flight_ID, extr_str)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID)
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Put observations and vars to plot into a dictionary
        sns.set(color_codes=True)
        sns.set_context(context, font_scale=font_scale)
        # Force alt to be in units of km
        ALT_var = 'Altitude (km)'
        Y_unit = ALT_var
#        df_mod[ALT_var] = AC.hPa_to_Km( df_mod['model-lev'].values )
        key4GEOSCF = 'GEOS-CF'
        for key in dfs_mod.keys():
            #        if key in dfs_mod.keys():
            df_mod = dfs_mod[key]
            if key4GEOSCF == key:
                df_mod[ALT_var] = AC.hPa_to_Km(df_mod['model-lev'].values)
            else:
                df_mod[ALT_var] = AC.hPa_to_Km(df_mod['PRESS'].values)
            dfs_mod[key] = df_mod
        df_obs[ALT_var] = df_obs['ALT_GIN'].values / 1E3
        # Online consider offline model output?
        if just_plot_GEOS_Chem:
            #            GEOSChem_varname = sorted(list(dfs_mod.keys())
            GCvarname = 'FP-Nest'
            data_d = {GCvarname: dfs_mod[GCvarname], 'Obs.': df_obs}
        else:
            data_d = AC.merge_two_dicts(dfs_mod, {'Obs.': df_obs})
        print('Plotting model runs: ', data_d.keys())
#        data_d = {'GEOS-CF': df_mod, 'Obs.':df_obs}

        # - Now plot up flight time series plots by variable
        title_str = "Altitude binned '{}' ({}) during flight '{}'"
        # Setup color dictinoary
        color_dict = {'GEOS-CF': 'red', 'Obs.': 'k'}
        CB_color_cycle = AC.get_CB_color_cycle()
        for n_key, key in enumerate(list(data_d.keys())):
            if key not in color_dict.keys():
                color_dict[key] = CB_color_cycle[n_key]
        unit_d = {}
        mod2obs_varnames = {
            'BrO': 'BrO',
            'HNO3': 'HNO3',
            'HNO2': 'HONO',
        }
        units_d = {
            'BrO': 'pptv', 'HNO3': 'pptv', 'HNO2': 'pptv', 'HONO': 'pptv',
        }
        range_d = {
            'HNO2': (-10, 60),
            'HNO3': (-30, 1500),
            'BrO': (-0.2, 1.0),
            'HONO': (-10, 60),
        }
        NOx_specs = ['HNO2', 'NOx', 'NO', 'NO2', 'HONO']
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
                print(n_key, key_, var2plot)
                #
                if key_ == 'Obs.':
                    varname = mod2obs_varnames[var2plot]
                else:
                    varname = var2plot
                # Setup an axis label
                units = units_d[var2plot]
                xlabel = '{} ({})'.format(var2plot, units)
                # Add alt to DataFrame
                df = pd.DataFrame({
                    var2plot: data_d[key_][varname],
                    ALT_var: data_d[key_][ALT_var]
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
                                                   widths=0.15,
                                                   dataset_num=n_key,
                                                   color=color_dict[key_])
                except:
                    pass
                # Beautify plot
                plt.legend()
                plt.title(title_str.format(var2plot, units, flight_ID))
                # Make NOx species be on a log scale
                xscale = 'linear'
                if (var2plot in NOx_specs) and NOxAsLog:
                    xscale = 'log'
                ax.set_xscale(xscale)
                if xscale == 'log':
                    if var2plot in ['HONO', 'HNO2']:
                        xlim = (0.03, 400)
                    else:
                        xlim = (0.3, 400)
                    ax.set_xlim(xlim)
                else:
                    try:
                        plt.xlim(range_d[var2plot])
                    except KeyError:
                        pass

            # Save to PDF
        #        fig.legend(loc='best', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
        #        plt.legend()
        #        plt.tight_layout()
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
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
    marker = 'o'
    s = 2
    alpha = 0.75
    cmap_list = AC.get_CB_color_cycle()
    # Now plot locations as scatter points on plot
    ax.scatter(lons, lats, color=cmap_list[0], s=s, marker=marker,
               alpha=alpha,
               label=flight_ID,
               transform=projection(), zorder=999
               )
    # Local area analysed as Cape Verde
    d = get_analysis_region('local_CVAO_area')
    extents = (d['x0'], d['x1'], d['y0'], d['y1'])
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


def plt_ts4ARNA_flight(df_obs=None, dfs_mod=None,
                       obs_label='Obs.',
                       mod_scale=1, obs_scale=1,
                       obs_adjustby=0, mod_adjustby=0,
                       ylim=(None, None),
                       units='ppbv', var2plot='CO',
                       ObsVar2Plot='CO_AERO',
                       ModVar2Plot='CO',
                       mod_label='GEOS-CF',
                       plt_alt_shadow=True,
                       plt_dust_as_backfill=True,
                       plt_BB_as_backfill=False,
                       aspect_ratio=0.25,
                       flight_ID='C216',
                       yscale='linear',
                       invert_yaxis=False,
                       title=None,
                       plt_legend=True,
                       context="paper", font_scale=0.75,
                       debug=False):
    """
    Plot up a timeseries of observations and model for a given flight
    """
    # Now use Seaborn settings
    sns.set(color_codes=True)
    sns.set_context(context, font_scale=font_scale)
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
    if not isinstance(ObsVar2Plot, type(None)):
        df_obs = df_obs[[ObsVar2Plot]].dropna()
        plt.plot(df_obs.index,
                 (df_obs[ObsVar2Plot].values+obs_adjustby)*obs_scale,
                 label=obs_label, color='k')
    # Setup model(s) variables lists regardless of whether these are plotted
    CB_cycle = AC.get_CB_color_cycle()
    mods2plot = list(dfs_mod.keys())
    mod_colours = dict(zip(CB_cycle, dfs_mod.keys()))
    ModelVar_CF = 'GEOS-CF'
    # If GEOS-CF is in the model list, then plot this last
    # This is because it used for the shadow altitude
    if debug:
        print(mods2plot)
    if len(mods2plot) >= 1:
        if ModelVar_CF in mods2plot:
            mods2plot.pop(mods2plot.index(ModelVar_CF))
            mods2plot += [ModelVar_CF]
        for n_key, key in enumerate(sorted(list(dfs_mod.keys()))):
            mod_colours[key] = CB_cycle[n_key]
    # Plot the model if requested...
    if not isinstance(ModVar2Plot, type(None)):
        mod_colours[ModelVar_CF] = 'red'
        # Now just loop and plot
        if debug:
            print(mods2plot)
        for mod2plot in mods2plot:
            df_mod = dfs_mod[mod2plot]
            # Exc. points in the model dataframe where there is no model output
            df_mod = df_mod.loc[~df_mod['model-lev'].isnull(), :]
            # Use the key of the dataframe dictionary as label, if >1 dataframe
            if len(mods2plot) > 1:
                mod_label = mod2plot
            # Plot model data
            plt.plot(df_mod.index,
                     (df_mod[ModVar2Plot].values+mod_adjustby)*mod_scale,
                     label=mod_label,
                     color=mod_colours[mod2plot])
    else:
        # Use the first model input
        df_mod = dfs_mod[mods2plot[0]]
        # Exc. points in the model dataframe where there is no model output
        df_mod = df_mod.loc[~df_mod['model-lev'].isnull(), :]

    # Get the beginning and end of the flight from the extracted model times
    xylim_min = AC.add_minutes(df_mod.index.min(), -15)
    xylim_max = AC.add_minutes(df_mod.index.max(), 15)
    xticks = df_mod.resample('15T').mean().index.values
    xticks = AC.dt64_2_dt(xticks)
    xticks_labels = [i.strftime('%H:%M') for i in xticks]
    # Get dates/datetimes of flight (using observed time in preference)
    try:
        sdate_str = df_obs.index.dropna().min().strftime('%x').strip()
        edate_str = df_obs.index.dropna().max().strftime('%x').strip()
        stime_str = df_obs.index.dropna().min().strftime('%H:%M').strip()
        etime_str = df_obs.index.dropna().max().strftime('%H:%M').strip()
    except ValueError:
        sdate_str = df_mod.index.dropna().min().strftime('%x').strip()
        edate_str = df_mod.index.dropna().max().strftime('%x').strip()
        stime_str = df_mod.index.dropna().min().strftime('%H:%M').strip()
        etime_str = df_mod.index.dropna().max().strftime('%H:%M').strip()
    # Set a shared title string and fill with variables specific to flight
    if isinstance(title, str):
        plt.title(title)
    else:
        title_str = "Timeseries of '{}' ({}) during flight '{}' on {}"
        plt.title(title_str.format(var2plot, units, flight_ID, sdate_str))
    # Beautify plot
    plt.yscale(yscale)
    plt.ylim(ylim)
    plt.ylabel('{} ({})'.format(var2plot, units))
    plt.xlim(xylim_min, xylim_max)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks_labels, rotation=45)
    # Invert the second y-axis
    if plt_alt_shadow:
        # Add a shadow of the altitude
        ax2 = ax.twinx()
        ModVar2Plot = 'model-lev'
        ax2.plot(df_mod.index, df_mod[ModVar2Plot].values,
                 label='Altitude',
                 color='grey', zorder=100, alpha=0.25)
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
        fig.legend(loc='best', bbox_to_anchor=(1, 1),
                   bbox_transform=ax.transAxes)
    plt.tight_layout()


def plt_ts4ARNA_flt_period_obs(df_obs=None,
                               dfs_mod=None,
                               dfs_mod_period=None,
                               obs_label='Obs.',
                               mod_scale=1, obs_adjustby=0,
                               ylim=(None, None),
                               units='ppbv', var2plot='CO',
                               ObsVar2Plot='CO_AERO',
                               ModVar2Plot='CO',
                               mod_label='GEOS-CF',
                               plt_alt_shadow=True,
                               aspect_ratio=0.25,
                               flight_ID='C216',
                               yscale='linear',
                               invert_yaxis=False,
                               plt_dust_as_backfill=True,
                               plt_errorbar=False,
                               ObsVar2PlotErr='',
                               StartVar='Start_time',
                               EndVar='End_time',
                               context="paper", font_scale=0.75,
                               title=None,
                               debug=False,
                               ):
    """
    Plot up a timeseries of observations and model for a given flight
    """
    # Now use Seaborn settings
    sns.set(color_codes=True)
    sns.set_context(context, font_scale=font_scale)
    # Setup the figure
    w, h = matplotlib.figure.figaspect(aspect_ratio)
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_subplot(111)
    # Setup model(s) variables lists regardless of whether these are plotted
    CB_cycle = AC.get_CB_color_cycle()
    mods2plot = list(dfs_mod.keys())
    mod_colours = dict(zip(CB_cycle, dfs_mod.keys()))
    ModelVar_CF = 'GEOS-CF'
    # If GEOS-CF is in the model list, then plot this last
    # This is because it used for the shadow altitude
    if debug:
        print(mods2plot)
    if len(mods2plot) >= 1:
        if ModelVar_CF in mods2plot:
            mods2plot.pop(mods2plot.index(ModelVar_CF))
            mods2plot += [ModelVar_CF]
            df_mod = dfs_mod[ModelVar_CF]
        else:
            df_mod = dfs_mod[list(dfs_mod.keys())[0]]
        for n_key, key in enumerate(sorted(list(dfs_mod.keys()))):
            mod_colours[key] = CB_cycle[n_key]

    # Plot up the model data if requested...
    if not isinstance(ModVar2Plot, type(None)):
        mod_colours[ModelVar_CF] = 'red'
        # Now just loop and plot
        if debug:
            print(mods2plot)
        for mod2plot in mods2plot:
            df_mod_period = dfs_mod_period[mod2plot]
            xmin = df_mod_period[StartVar]
            xmax = df_mod_period[EndVar]
            plt.hlines(df_mod_period[ModVar2Plot].values*mod_scale,
                       xmin, xmax,
                       label=mod2plot, color=mod_colours[mod2plot])

    else:
        # Use the first model input
        # This needs to be GEOS-CF (for model-lev) variable
        df_mod = dfs_mod[mods2plot[0]]
        # Exc. points in the model dataframe where there is no model output
        df_mod = df_mod.loc[~df_mod['model-lev'].isnull(), :]

    # Plot up where
    if plt_dust_as_backfill:
        colour_plot_background_by_bool(df=df_mod, ax=ax,
                                       bool2use='IS_DUST',
                                       color='sandybrown',
                                       alpha=0.25,
                                       label='Dust (non-MBL)')
    # Now plot up the observations
    if plt_errorbar:
        ax = plt.gca()
        ax.errorbar(df_obs.index, df_obs[ObsVar2Plot].values+obs_adjustby,
                    yerr=df_obs[ObsVar2PlotErr].values, fmt='-o',
                    color='k', capsize=3.0, label=obs_label)
    else:
        xmin = df_obs[StartVar]
        xmax = df_obs[EndVar]
        plt.hlines(df_obs[ObsVar2Plot].values+obs_adjustby,
                   xmin, xmax,
                   label=obs_label, color='k')

    # Exclude points in the model dataframe where there is no model output
    df_mod = df_mod.loc[~df_mod['model-lev'].isnull(), :]
    # Get the beginning and end of the flight from the extracted model times
    xylim_min = AC.add_minutes(df_mod.index.min(), -15)
    xylim_max = AC.add_minutes(df_mod.index.max(), 15)
    xticks = df_mod.resample('15T').mean().index.values
    xticks = AC.dt64_2_dt(xticks)
    xticks_labels = [i.strftime('%H:%M') for i in xticks]
    # Get dates/datetimes of flight
    sdate_str = df_mod.index.min().strftime('%x').strip()
    edate_str = df_mod.index.max().strftime('%x').strip()
    stime_str = df_mod.index.min().strftime('%H:%M').strip()
    etime_str = df_mod.index.max().strftime('%H:%M').strip()
    # Beautify plot
    # Set a shared title string and fill with variables specific to flight
    if isinstance(title, str):
        plt.title(title)
    else:
        title_str = "Timeseries of '{}' ({}) during flight '{}' on {}"
        plt.title(title_str.format(var2plot, units, flight_ID, sdate_str))
    plt.yscale(yscale)
    plt.ylim(ylim)
    plt.ylabel('{} ({})'.format(var2plot, units))
    plt.xlim(xylim_min, xylim_max)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks_labels, rotation=45)
    if plt_alt_shadow:
        # Add a shadow of the altitude
        ax2 = ax.twinx()
        ModVar2Plot = 'model-lev'
        ax2.plot(df_mod.index, df_mod[ModVar2Plot].values,
                 label='Altitude',
                 color='grey', zorder=100, alpha=0.25)
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
    fig.legend(loc='best', bbox_to_anchor=(1, 1),
               bbox_transform=ax.transAxes)
    plt.tight_layout()


def plt_ts4ARNA_flt_point_obs(df_obs=None, df_mod=None,
                              obs_label='Obs.',
                              mod_scale=1, obs_adjustby=0,
                              ylim=(None, None),
                              units='ppbv', var2plot='CO',
                              ObsVar2Plot='CO_AERO',
                              ModVar2Plot='CO',
                              mod_label='GEOS-CF',
                              plt_alt_shadow=True,
                              aspect_ratio=0.25,
                              flight_ID='C216',
                              yscale='linear',
                              invert_yaxis=False,
                              plt_dust_as_backfill=True,
                              plt_errorbar=False,
                              ObsVar2PlotErr='',
                              context="paper", font_scale=0.75,
                              title=None,
                              ):
    """
    Plot up a timeseries of observations and model for a given flight
    """
    # Now use Seaborn settings
    sns.set(color_codes=True)
    sns.set_context(context, font_scale=font_scale)
    # Exclude points in the model dataframe where there is no model output
    df_mod = df_mod.loc[~df_mod['model-lev'].isnull(), :]
    # Get the beginning and end of the flight from the extracted model times
    xylim_min = AC.add_minutes(df_mod.index.min(), -15)
    xylim_max = AC.add_minutes(df_mod.index.max(), 15)
    xticks = df_mod.resample('15T').mean().index.values
    xticks = AC.dt64_2_dt(xticks)
    xticks_labels = [i.strftime('%H:%M') for i in xticks]
    # Get dates/datetimes of flight
    sdate_str = df_mod.index.min().strftime('%x').strip()
    edate_str = df_mod.index.max().strftime('%x').strip()
    stime_str = df_mod.index.min().strftime('%H:%M').strip()
    etime_str = df_mod.index.max().strftime('%H:%M').strip()
    # Setup the figure
    w, h = matplotlib.figure.figaspect(aspect_ratio)
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_subplot(111)
    # Plot up the model data...
    if not isinstance(ModVar2Plot, type(None)):
        plt.plot(df_mod.index, df_mod[ModVar2Plot].values*mod_scale,
                 label=mod_label, color='red')

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
        ax.errorbar(df_obs.index, df_obs[ObsVar2Plot].values+obs_adjustby,
                    yerr=df_obs[ObsVar2PlotErr].values, fmt='-o',
                    color='k', capsize=3.0, label=obs_label)
    else:
        #        df_obs = df_obs[[ObsVar2Plot]].dropna()
        plt.scatter(df_obs.index, df_obs[ObsVar2Plot].values+obs_adjustby,
                    label=obs_label, color='k')

    # Beautify plot
    if isinstance(title, str):
        plt.title(title)
    else:
        title_str = "Timeseries of '{}' ({}) during flight '{}' on {}"
        plt.title(title_str.format(var2plot, units, flight_ID, sdate_str))
    plt.yscale(yscale)
    plt.ylim(ylim)
    plt.ylabel('{} ({})'.format(var2plot, units))
    plt.xlim(xylim_min, xylim_max)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks_labels, rotation=45)
    # Invert the second y-axis
    if plt_alt_shadow:
        # Add a shadow of the altitude
        ax2 = ax.twinx()
        ModVar2Plot = 'model-lev'
        ax2.plot(df_mod.index, df_mod[ModVar2Plot].values,
                 label='Altitude',
                 color='grey', zorder=100, alpha=0.25)
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
    fig.legend(loc='best', bbox_to_anchor=(1, 1),
               bbox_transform=ax.transAxes)
    plt.tight_layout()


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
        #        df = get_GEOSCF4flightnum(flight_ID=flight_ID)
        df_FAAM = get_FAAM_core4flightnum(flight_ID=flight_ID,
                                          resample_data=True)
    # Just add is the "is dust" and "is SLR"
    vars2use = ['IS_DUST', 'IS_SLR']
    df = pd.concat([df, df_FAAM[vars2use]],  axis=1)
    return df


def add_secs2duplicate_index_values(df):
    """
    Duplicate values for index are not limit panadas operations
    """
    # Add a number to each ro
    RowVar = 'Row'
    df[RowVar] = np.arange(df.shape[0])
    # Loop index and add a second to any duplicates
    NewIndex = []
    last_index = None
    penultimate_index = None
    for Row in df[RowVar].values:
        index_val = df.loc[df[RowVar] == Row, :].index.values[0]
        index_val = AC.dt64_2_dt([index_val])[0]
        duplicates = df.loc[df.index == index_val, :].shape[0]
        if duplicates > 1:
            print(Row, df.loc[df[RowVar] == Row, :])
            index_val = AC.add_secs(index_val, 1)
            # Update the index
            index_values = df.index.values
            index_values[Row] = index_val
            df.index = index_values
    return df


def plt_ts_comp4ARNA_flights_CIMS(dpi=320, context='paper',
                                  flight_nums=[],
                                  RunSet=None, res='4x5',
                                  just_plot_GEOS_Chem=False,
                                  inc_GEOSChem=False,
                                  LatVar='model-lat',
                                  LonVar='model-lon',
                                  show_plot=False):
    """
    Plot up timeseries comparisons between core observations and model data
    """
    import seaborn as sns
    # Which flights to plot? - Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [218, 219, 220, 221, 222, 223, 224, 225, ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod_CF = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF4flightnum(flight_ID=flight_ID)
        # Add the derived variables to the dataframe
        df = add_deriv_vars2df(df=df)
        dfs_mod_CF[flight_ID] = df
    # Model - GEOS-Chem (offline)
    if inc_GEOSChem:
        dfs_mod_GC = {}
        for flight_ID in flight_IDs:
            df = get_GEOSChem4flightnum(flight_ID=flight_ID, res=res,
                                        RunSet=RunSet,)
            # Add the derived variables to the dataframe
            df = add_deriv_vars2df(df=df)
            dfs_mod_GC[flight_ID] = df
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
#        df_obs = dfs_obs[flight_ID]
#        df_mod = dfs_mod[flight_ID]
        df_obs = dfs_obs[flight_ID]
        df_mod_CF = dfs_mod_CF[flight_ID]
        if inc_GEOSChem:
            if just_plot_GEOS_Chem:
                dfs_mod = dfs_mod_GC[flight_ID]
                mod_label_master = RunSet

            else:
                dfs_mod_GC4flight = dfs_mod_GC[flight_ID]
                dfs_mod = {'GEOS-CF': df_mod_CF}
                for key in list(dfs_mod_GC4flight.keys()):
                    dfs_mod[key] = dfs_mod_GC4flight[key]
                mod_label_master = 'GEOS-CF'
        else:
            mod_label_master = 'GEOS-CF'
            dfs_mod = {'GEOS-CF': df_mod_CF}
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_timeseries_flighttrack_{}_CIMS'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        try:
            plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID,
                                               LatVar=LatVar, LonVar=LonVar,)
            # Save to PDF
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format('spatial', flight_ID))
        # - Now plot up flight time series plots by variable
        try:
            units = 'ppt'
            var2plot = 'BrO'
            ObsVar2Plot = 'BrO'
            ModVar2Plot = 'BrO'
            mod_label = mod_label_master
            mod_scale = 1E12
            ylim = (-0.2, 1)
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))

        # - Plot up Nitric acid
        try:
            units = 'ppt'
            var2plot = 'HNO3'
            ObsVar2Plot = 'HNO3'
            ModVar2Plot = 'HNO3'
            mod_label = mod_label_master
            mod_scale = 1E12
            ylim = (-30, 1500)
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))

        # - Plot up Nitrous acid
        try:
            units = 'ppt'
            var2plot = 'HONO'
            ObsVar2Plot = 'HONO'
            ModVar2Plot = 'HNO2'
            mod_label = mod_label_master
            mod_scale = 1E12
            ylim = (-10, 25)
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))

        # - Plot up HCN
        try:
            units = 'ppt'
            var2plot = 'HCN'
            ObsVar2Plot = 'HCN'
            ModVar2Plot = None
            mod_label = mod_label_master
            mod_scale = None
            ylim = (-10, 120)
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))

        # - Add a plot with masking for biomass burning.
        try:
            #
            units = 'ppt'
            var2plot = 'HCN'
            ObsVar2Plot = 'HCN'
            ModVar2Plot = None
            mod_label = mod_label_master
            mod_scale = None
            ylim = (-10, 120)
            title_str = "Timeseries of '{}' ({}) during flight '{}' - {}"
            ext_str = 'Biomass flagged as HCN @ background+3$\sigma$'
            title = title_str.format(var2plot, units, flight_ID, ext_str)
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               plt_legend=False,
                               title=title,
                               context=context,
                               )
            # Add a flat threshold for biomass burning
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
            fig.legend(loc='best', bbox_to_anchor=(1, 1),
                       bbox_transform=ax.transAxes)
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))
        # - Add a plot with masking for biomass burning.
        try:
            #
            units = 'ppt'
            var2plot = 'HCN'
            ObsVar2Plot = 'HCN'
            ModVar2Plot = None
            mod_label = mod_label_master
            mod_scale = None
            ylim = (-10, 120)
            title_str = "Timeseries of '{}' ({}) during flight '{}' - {}"
            ext_str = 'Biomass flagged as HCN @ background+2$\sigma$'
            title = title_str.format(var2plot, units, flight_ID, ext_str)
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               plt_legend=False,
                               title=title,
                               context=context,
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
            fig.legend(loc='best', bbox_to_anchor=(1, 1),
                       bbox_transform=ax.transAxes)
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))
        # - Add a plot with masking for biomass burning.
        try:
            #
            units = 'ppt'
            var2plot = 'HCN'
            ObsVar2Plot = 'HCN'
            ModVar2Plot = None
            mod_label = mod_label_master
            mod_scale = None
            ylim = (-10, 120)
            title_str = "Timeseries of '{}' ({}) during flight '{}' - {}"
            ext_str = 'Biomass flagged as HCN @ background+1$\sigma$'
            title = title_str.format(var2plot, units, flight_ID, ext_str)
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               plt_legend=False,
                               title=title,
                               context=context,
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
            fig.legend(loc='best', bbox_to_anchor=(1, 1),
                       bbox_transform=ax.transAxes)
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))
        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def only_use_filter_times(df, FILTERdf=None, flight_ID='C225',
                          average4period=True,
                          StartVar='Start_time', EndVar='End_time',
                          ):
    """
    Chop out only the values for the filter sample periods
    """
    # Loop by filter
    dfs = []
    for filter in FILTERdf['Filter'].values:
        #        print( filter )
        FILTERdf_tmp = FILTERdf.loc[FILTERdf['Filter'] == filter]
        dt = FILTERdf_tmp.index.values
        sdate = FILTERdf_tmp[StartVar].values[0]
        edate = FILTERdf_tmp[EndVar].values[0]
        # Select just the sampling period
        df4filter = df.loc[df.index >= sdate, :]
        df4filter = df4filter.loc[df4filter.index <= edate, :]
        if average4period:
            df4filter = df4filter.mean(axis=0)
        # Add start time and end time to the dataframe
        df4filter[StartVar] = sdate
        df4filter[EndVar] = edate
        # Save to list
        dfs += [df4filter.copy()]
        del df4filter
    # Return averaged over filter sample period
    if average4period:
        return pd.DataFrame(dfs, index=FILTERdf.index)
    else:
        return pd.DataFrame(pd.concat(dfs, axis=0))


def mk_combined_NOy_obs_variable(FAAMdf=None, CIMSdf=None, Filtersdf=None,
                                 bin_dfs4filter_periods=False, NOyVar='NOy',
                                 RTN_Filtersdf=True, flight_ID='C225'):
    """
    Make a combined NOy to compare against modelled NOy
    ---- Notes
     - GEOS-CF NOy is
       no_no2_hno3_hno4_hono_2xn2o5_pan_organicnitrates_aerosolnitrates
    """
    # Bin data to filter observation periods
    if bin_dfs4filter_periods:
        print('WARNING: binning of dataframes to filter perios not setup!')
    # Use the filter dataframe as the basis to add otther species too too
    # NOTE: conversion to pptv hs already been done!
    SpeciesVar = 'Total_{}_{}'.format('NO3', 'ppt')
    Filtersdf[NOyVar] = Filtersdf[SpeciesVar].copy()
    # add NO, NO2 - from Faam obs
    FAAMdf = FAAMdf.copy().replace(np.NaN, 0)
    var2use = 'NO_pptV'
    try:
        Filtersdf.loc[:, NOyVar] += FAAMdf[var2use].values
    except KeyError:
        pstr = "WARNING: '{}' not added to '{}' variable for flight '{}'"
        print(pstr.format(var2use, NOyVar, flight_ID))
    var2use = 'NO2_pptV'
    try:
        Filtersdf.loc[:, NOyVar] += FAAMdf[var2use].values
    except KeyError:
        pstr = "WARNING: '{}' not added to '{}' variable for flight '{}'"
        print(pstr.format(var2use, NOyVar, flight_ID))
    # Add Nitrate rack HONO too?
    # Add HNO3 from CIMS
    var2use = 'HNO3'
    CIMSdf = CIMSdf.copy().replace(np.NaN, 0)
    try:
        Filtersdf.loc[:, NOyVar] += CIMSdf[var2use].values
    except KeyError:
        pstr = "WARNING: '{}' not added to '{}' variable for flight '{}'"
        print(pstr.format(var2use, NOyVar, flight_ID))
    # Add HONO from CIMS for now
    var2use = 'HONO'
    try:
        Filtersdf.loc[:, NOyVar] += CIMSdf[var2use].values
    except KeyError:
        pstr = "WARNING: '{}' not added to '{}' variable for flight '{}'"
        print(pstr.format(var2use, NOyVar, flight_ID))
    if RTN_Filtersdf:
        return Filtersdf
    else:
        return Filtersdf[[NOyVar]]


def plt_ts_comp4MOYA_flights(dpi=320, inc_GEOSChem=False,
                             show_plot=False, context='talk'):
    """
    Plot up timeseries comparisons between core observations and model data
    """
    import seaborn as sns
    # Which flights to plot?
    flight_IDs = ['C006', 'C007']
    sdate_d = {
        'C006': datetime.datetime(2017, 3, 1),
        'C007': datetime.datetime(2017, 3, 2),
    }
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod = {}
#    dfs_mod_GC = {}
    for flight_ID in flight_IDs:
        # No GEOS-CF availiblie for MOYA flights (pre 2018)
        #        dfs_mod[flight_ID] = get_GEOSCF4flightnum(flight_ID=flight_ID)
        #
        RunSet = 'FP-MOYA-Nest'
        res = '0.25x0.3125'
        sdate = sdate_d[flight_ID]
        dfs_mod_GC = get_GEOSChem4flightnum(flight_ID=flight_ID,
                                            res=res,
                                            RunSet=RunSet,
                                            sdate=sdate,
                                            )
        dfs_mod[flight_ID] = dfs_mod_GC[list(dfs_mod_GC.keys())[0]]

#    dfs_mod = dfs_mod_GC
    # Observations
    folder = '/users/ts551/scratch/data/FAAM/core_faam_NetCDFs/'
    dfs_obs = {}
    for flight_ID in flight_IDs:
        dfs_obs[flight_ID] = get_FAAM_core4flightnum(flight_ID=flight_ID,
                                                     folder=folder)
    # - Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        # Setup PDF to save PDF plots to
        savetitle = 'MOYA_timeseries_flighttrack_{}'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Now plot up flight time series plots by variable
        # - Plot up carbon monoxide
        units = 'ppbv'
        var2plot = 'CO'
        ObsVar2Plot = 'CO_AERO'
        ModVar2Plot = 'CO'
        mod_label = RunSet
        mod_scale = 1E9
        ylim = (50, 400)
        # Call timeseries plotter function
        plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                           ObsVar2Plot=ObsVar2Plot,
                           mod_scale=mod_scale,
                           mod_label=mod_label,
                           ModVar2Plot=ModVar2Plot,
                           ylim=ylim,
                           dfs_mod={RunSet: dfs_mod[flight_ID]},
                           df_obs=df_obs,
                           flight_ID=flight_ID,
                           context=context,
                           )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up ozone
        units = 'ppbv'
        var2plot = 'Ozone'
        ObsVar2Plot = 'O3_TECO'
        ModVar2Plot = 'O3'
        mod_label = RunSet
        mod_scale = 1E9
        ylim = (10, 100)
        # Call timeseries plotter function
        plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                           ObsVar2Plot=ObsVar2Plot,
                           mod_scale=mod_scale,
                           mod_label=mod_label,
                           ModVar2Plot=ModVar2Plot,
                           ylim=ylim,
                           dfs_mod={RunSet: dfs_mod[flight_ID]},
                           df_obs=df_obs,
                           flight_ID=flight_ID,
                           context=context,
                           )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_ts_comp4MOYA_flights_PHYSICAL_VARS(dpi=320, show_plot=False,
                                           debug=False,
                                           context='paper'):
    """
    Plot up overview for MOYA campaign flights
    """
    import seaborn as sns
    # Which flights to plot?
    flight_IDs = ['C006', 'C007']
    sdate_d = {
        'C006': datetime.datetime(2017, 3, 1),
        'C007': datetime.datetime(2017, 3, 2),
    }
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod = {}
#    dfs_mod_GC = {}
    for flight_ID in flight_IDs:
        # No GEOS-CF availiblie for MOYA flights (pre 2018)
        #        dfs_mod[flight_ID] = get_GEOSCF4flightnum(flight_ID=flight_ID)
        #
        RunSet = 'FP-MOYA-Nest'
        res = '0.25x0.3125'
        sdate = sdate_d[flight_ID]
        dfs_mod_GC = get_GEOSChem4flightnum(flight_ID=flight_ID,
                                            res=res,
                                            RunSet=RunSet,
                                            sdate=sdate,
                                            )
        dfs_mod[flight_ID] = dfs_mod_GC[list(dfs_mod_GC.keys())[0]]

#    dfs_mod = dfs_mod_GC
    # Observations
    folder = '/users/ts551/scratch/data/FAAM/core_faam_NetCDFs/'
    dfs_obs = {}
    for flight_ID in flight_IDs:
        dfs_obs[flight_ID] = get_FAAM_core4flightnum(flight_ID=flight_ID,
                                                     folder=folder)

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        df_mod = dfs_mod[flight_ID]
        # Setup PDF to save PDF plots to
        savetitle = 'MOYA_timeseries_flighttrack_{}_PHYSICAL_VARS'.format(
            flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up ROLL_GIN
        try:
            units = 'ADD THIS'
    #        VarName = 'PCAS2CON'
    #        FlagName = 'PCAS2_FLAG'
            var2plot = 'ROLL_GIN'
            ObsVar2Plot = 'ROLL_GIN'
            ModVar2Plot = None
            mod_label = RunSet
            mod_scale = None
    #        ylim = (10, 100)
    #        ylim = None
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               #                                   ylim=ylim,
                               dfs_mod={mod_label: df_mod},
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               plt_legend=False,
                               context=context,
                               )
            # Colour in SLRs
            ax = plt.gca()
            colour_plot_background_by_bool(df=df_obs, ax=ax, bool2use='IS_SLR',
                                           color='blue', alpha=0.25,
                                           label='SLR')
            fig = plt.gcf()
            fig.legend(loc='best', bbox_to_anchor=(1, 1),
                       bbox_transform=ax.transAxes)
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot ROLL_GIN')

        # - Plot up VELD_GIN
        try:
            units = 'ADD THIS'
    #        VarName = 'PCAS2CON'
    #        FlagName = 'PCAS2_FLAG'
            var2plot = 'VELD_GIN'
            ObsVar2Plot = 'VELD_GIN'
            ModVar2Plot = None
            mod_label = RunSet
            mod_scale = None
    #        ylim = (10, 100)
            ylim = None
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod={mod_label: df_mod},
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               plt_legend=False,
                               context=context,
                               )
            # Colour in SLRs
            ax = plt.gca()
            colour_plot_background_by_bool(df=df_obs, ax=ax, bool2use='IS_SLR',
                                           color='blue', alpha=0.25,
                                           label='SLR')
            fig = plt.gcf()
            fig.legend(loc='best', bbox_to_anchor=(1, 1),
                       bbox_transform=ax.transAxes)

            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot VELD_GIN')

        # - Plot up temperature
        units = '$^{\circ}$C'
        var2plot = 'Temperature'
        ObsVar2Plot = 'TAT_DI_R'
        ModVar2Plot = 'T'
        mod_label = RunSet
        mod_scale = 1
        ylim = (-25, 35)
        obs_adjustby = -273.15
        # Call timeseries plotter function
        plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                           ObsVar2Plot=ObsVar2Plot,
                           mod_scale=mod_scale,
                           mod_label=mod_label,
                           ModVar2Plot=ModVar2Plot,
                           ylim=ylim,
                           dfs_mod={mod_label: df_mod},
                           df_obs=df_obs,
                           flight_ID=flight_ID,
                           obs_adjustby=obs_adjustby,
                           context=context,
                           )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up Eastward wind
        units = 'm s$^{-1}$'
        var2plot = 'Eastward wind'
        ObsVar2Plot = 'U_C'
        ModVar2Plot = 'U'
        mod_label = RunSet
        mod_scale = 1
        ylim = (-25, 25)
        # Call timeseries plotter function
        plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                           ObsVar2Plot=ObsVar2Plot,
                           mod_scale=mod_scale,
                           mod_label=mod_label,
                           ModVar2Plot=ModVar2Plot,
                           ylim=ylim,
                           dfs_mod={mod_label: df_mod},
                           df_obs=df_obs,
                           flight_ID=flight_ID,
                           context=context,
                           )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up Northward wind
        units = 'm s$^{-1}$'
        var2plot = 'Northward wind'
        ObsVar2Plot = 'V_C'
        ModVar2Plot = 'V'
        mod_label = RunSet
        mod_scale = 1
        ylim = (-25, 25)
        # Call timeseries plotter function
        plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                           ObsVar2Plot=ObsVar2Plot,
                           mod_scale=mod_scale,
                           mod_label=mod_label,
                           ModVar2Plot=ModVar2Plot,
                           ylim=ylim,
                           dfs_mod={mod_label: df_mod},
                           df_obs=df_obs,
                           flight_ID=flight_ID,
                           context=context,
                           )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up Latitude
        try:
            units = '$^{\circ}$N'
            var2plot = 'Latitude'
            ObsVar2Plot = 'LAT_GIN'
            ModVar2Plot = 'LAT'
            mod_label = RunSet
            mod_scale = 1
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               dfs_mod={mod_label: df_mod},
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up Longitude
        try:
            units = '$^{\circ}$E'
            var2plot = 'Longitude'
            ObsVar2Plot = 'LON_GIN'
            ModVar2Plot = 'LON'
            mod_label = RunSet
            mod_scale = 1
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               dfs_mod={mod_label: df_mod},
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up altitude
        try:
            units = 'hPa'
            var2plot = 'Altitude'
            ModVar2Plot = 'PRESS'
            ObsVar2Plot = 'PS_RVSM'  # Use pressure measurement
#            ObsVar2Plot = 'ALT_GIN' # Use GPS altitude?
#            vals = AC.hPa_to_Km( df_obs['ALT_GIN'].values/1E3, reverse=True )
            mod_label = RunSet
            mod_scale = 1
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               dfs_mod={mod_label: df_mod},
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               invert_yaxis=True,
                               plt_alt_shadow=False,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
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


def get_derived_total_NOy4flights(flight_nums=[], res='4x5', RunSet=None,
                                  inc_GEOSChem=True, CoreRunsOnly=False,
                                  inc_GEOSCF=False,
                                  debug=False):
    """
    Get derived total NOy from observations (and equivalent model values)
    """
    # Which flights to plot? - Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [219, 220, 221, 222, 223, 224, 225, ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Setup Observation (filters) dataframe
    dfs_obs = get_filters_data4flight()
    dfs_obs = add_secs2duplicate_index_values(dfs_obs)
    # Model
    if inc_GEOSCF:
        mod_label_master = 'GEOS-CF'
        dfs_mod_CF = {}
        for flight_ID in flight_IDs:
            df = get_GEOSCF4flightnum(flight_ID=flight_ID)
            # Add the derived variables to the dataframe
            df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
            df = add_deriv_vars2df(df=df)
            # Extra actions for this specific function
            df['flight_ID'] = flight_ID
            dfs_mod_CF[flight_ID] = df.copy()
            dfs_mod_CF[flight_ID] = df
    # Model - GEOS-Chem (offline)
    if inc_GEOSChem:
        mod_label_master = RunSet
        dfs_mod_GC = {}
        for flight_ID in flight_IDs:
            dfs = get_GEOSChem4flightnum(flight_ID=flight_ID, res=res,
                                         CoreRunsOnly=CoreRunsOnly,
                                         RunSet=RunSet,)
            for key in dfs.keys():
                df = dfs[key]
                df = add_derived_FAAM_flags2df4flight(df=df,
                                                      flight_ID=flight_ID)
                df = add_deriv_vars2df(df=df)
                dfs[key] = df
            dfs_mod_GC[flight_ID] = dfs

    # KLUDGE -  just plot GEOS-Chem for now
    # (previously just plotted just GEOS-CF)
    print('WARNING: TODO setup to plot multiple model runs')
    dfs_mod = dfs_mod_GC
#    dfs_mod_period = dfs_mod_GC_period

    # Setup dictionary of observation dataframes (CIMS)
    dfs_CIMS = {}
    dfs_CIMS_period = {}
    for flight_ID in flight_IDs:
        df = get_CIMS_data4flight(flight_ID=flight_ID)
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_CIMS[flight_ID] = df.copy()
        # Just consider values during filter observation period
        FILTERdf = dfs_obs.loc[dfs_obs['Flight'] == flight_ID, :]
        df = only_use_filter_times(df=df, FILTERdf=FILTERdf,
                                   flight_ID=flight_ID)
        dfs_CIMS_period[flight_ID] = df
    # Setup dictionary of observation dataframes (FAAM)
    dfs_FAAM = {}
    dfs_FAAM_period = {}
    for flight_ID in flight_IDs:
        df = get_FAAM_core4flightnum(flight_ID=flight_ID)
        df['flight_ID'] = flight_ID
        dfs_FAAM[flight_ID] = df.copy()
        # Just consider values during filer observation period
        FILTERdf = dfs_obs.loc[dfs_obs['Flight'] == flight_ID, :]
        df = only_use_filter_times(df=df, FILTERdf=FILTERdf,
                                   flight_ID=flight_ID)
        dfs_FAAM_period[flight_ID] = df
    # Observations - combine to make a single NOy variable
    dfs_obs_NOy = {}
    for flight_ID in flight_IDs:
        # Get observations and model timeseries data as a DataFrame
        Filtersdf = dfs_obs.loc[dfs_obs['Flight'] == flight_ID, :]
        CIMSdf = dfs_CIMS_period[flight_ID]
        FAAMdf = dfs_FAAM_period[flight_ID]
        # Now make a combined NOy variable.
        df = mk_combined_NOy_obs_variable(Filtersdf=Filtersdf,
                                          FAAMdf=FAAMdf,
                                          CIMSdf=CIMSdf,
                                          flight_ID=flight_ID,
                                          )
        dfs_obs_NOy[flight_ID] = df

    # KLUDGE -  just plot GEOS-Chem for now
    # (previously just plotted just GEOS-CF)
    print('WARNING: TODO setup to plot multiple model runs')
    dfs_mod = dfs_mod_GC

    # Only consider the model for the observation (filter) period
    dfs_mod_period = {}
    for flight_ID in flight_IDs:
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs.loc[dfs_obs['Flight'] == flight_ID, :]
#        df_mod = dfs_mod[flight_ID]
        # Kludge! force use of single set of model output for now...
        print('WARNING: TODO setup to plot multiple model runs')
        ModelVarName = list(dfs_mod[flight_ID].keys())[0]
        dfs_mod4flight = dfs_mod[flight_ID]  # [ModelVarName]
        for key in dfs_mod4flight.keys():
            df_mod = dfs_mod4flight[key]
            df_mod = only_use_filter_times(df=df_mod, FILTERdf=df_obs,
                                           flight_ID=flight_ID)
            dfs_mod4flight[key] = df_mod
        dfs_mod_period[flight_ID] = dfs_mod4flight

    # Return a list of dictionaries
    gc.collect()
    RtnList = dfs_obs, dfs_mod, dfs_mod_period, dfs_obs_NOy, dfs_CIMS_period, dfs_FAAM_period
    return RtnList


def plt_comp_by_alt_4ARNA_NOy(vars2plot=['NOy', 'NOy-gas'],
                              flight_nums=[],
                              res='4x5', RunSet=None,
                              inc_GEOSChem=True,
                              inc_GEOSCF=False,
                              CoreRunsOnly=False,
                              pdff=None, savetitle=None,
                              show_plot=False, close_pdf=True,
                              debug=False, context='paper', dpi=320,
                              ):
    """
    Plot up a comparison of derived NOy for all flights
    """
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context(context)
    # Which flights to plot? - Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [219, 220, 221, 222, 223, 224, 225, ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # Get model and observation data for filter periods
    DataList = get_derived_total_NOy4flights(flight_nums=flight_nums,
                                             res=res, RunSet=RunSet,
                                             CoreRunsOnly=CoreRunsOnly,
                                             inc_GEOSChem=inc_GEOSChem,
                                             inc_GEOSCF=inc_GEOSCF,
                                             debug=debug)
#    dfs_obs, dfs_mod, dfs_mod_period, dfs_obs_NOy, dfs_CIMS_period, dfs_FAAM_period = DataList
    NIU, NIU, dfs_mod_period, dfs_obs_NOy, NIU, NIU = DataList

    # Make a single DataFrame
    Mods2Plot = list(dfs_mod_period[flight_IDs[0]].keys())
    dfsMod = {}
    for mod in Mods2Plot:
        ModByFlight = [dfs_mod_period[i][mod] for i in dfs_mod_period.keys()]
        dfsMod[mod] = pd.concat(ModByFlight, axis=0)
    # Convert observations to points?
    dfObs = pd.concat([dfs_obs_NOy[i] for i in dfs_obs_NOy.keys()], axis=0)

    # Setup PDF to save PDF plots to
    if isinstance(savetitle, type(None)):
        savetitle = 'ARNA_altitude_binned_{}{}'.format('NOy', '')
    if isinstance(pdff, type(None)):
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # Loop by species to plot
    colors2use = AC.get_CB_color_cycle()
    markers2use = [
        'o', 'v', '^', '<', '>', 'p', '*', 'h', 'H', 'D', 'd', 'P',
        'X'
        '8', 's',
    ]
    for n_var2plot, var2plot in enumerate(vars2plot):
        # Plot observations
        plt.scatter(dfObs['NOy'].values,
                    dfObs['Average_altitude_m'].values/1E3,
                    label='Obs.',
                    color='k',
                    alpha=0.75)

        # loop by model and plot
        for n_Mod2Plot, Mod2Plot in enumerate(Mods2Plot):
            # Get pressure variable and plot
            df2plot = dfsMod[Mod2Plot]
            plt.scatter(df2plot[var2plot].values*1E12,
                        AC.hPa_to_Km(df2plot['PRESS'].values),
                        label=Mod2Plot,
                        color=colors2use[n_Mod2Plot],
                        marker=markers2use[n_Mod2Plot],
                        alpha=0.6
                        )
            # TODO: Add option to optionally plot model as a line

        # Beautify plot
        TitleStr = "Vertical distribution of '{}' during all ARNA-2 flights"
        plt.title(TitleStr.format(var2plot))
        ALT_var = 'Altitude (km)'
        plt.ylabel(ALT_var)
        plt.xlabel('{} (pptv)'.format('NOy'))
        plt.legend()
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        if show_plot:
            plt.show()
        plt.close()

    # Save entire pdf
    if close_pdf:
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
    plt.close('all')


def plt_ts_comp4ARNA_flights_NOy_ALL(dpi=320, show_plot=False,
                                     flight_nums=[],
                                     res='4x5', RunSet=None,
                                     inc_GEOSChem=True,
                                     inc_GEOSCF=False,
                                     just_plot_GEOS_Chem=False,
                                     LatVar='model-lat',
                                     LonVar='model-lon',
                                     plt_alt_shadow=False,
                                     PltSpatialLocsOfdata=True,
                                     debug=False, context='paper'):
    """
    Plot up timeseries comparisons between filter samples and model data
    """
    import seaborn as sns
    # Which flights to plot? - Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [219, 220, 221, 222, 223, 224, 225, ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # Get model and observation data for filter periods
    DataList = get_derived_total_NOy4flights(flight_nums=flight_nums,
                                             res=res, RunSet=RunSet,
                                             inc_GEOSChem=inc_GEOSChem,
                                             inc_GEOSCF=inc_GEOSCF,
                                             debug=debug)
    dfs_obs, dfs_mod, dfs_mod_period, dfs_obs_NOy, dfs_CIMS_period, dfs_FAAM_period = DataList

    # TODO - Update code to make temp variables below redundant
    mod_label_master = None

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs.loc[dfs_obs['Flight'] == flight_ID, :]
        print('WARNING: TODO setup to plot multiple model runs')
#        df_mod = dfs_mod[flight_ID]
        df_mod = dfs_mod[flight_ID][list(dfs_mod[flight_ID].keys())[0]]
        df_mod_period = dfs_mod_period[flight_ID]
        df_obs_NOy = dfs_obs_NOy[flight_ID]
        df_CIMS_period = dfs_CIMS_period[flight_ID]
        df_FAAM_period = dfs_FAAM_period[flight_ID]
        # add the dust and SLR flags to the core dataframe
        df_obs = add_derived_FAAM_flags2df4flight(df=df_obs,
                                                  flight_ID=flight_ID)
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_timeseries_flighttrack_{}_NOy'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        if PltSpatialLocsOfdata:
            plt_flightpath_spatially_over_CVAO(df=df_mod, flight_ID=flight_ID,
                                               LatVar=LatVar,
                                               LonVar=LonVar,
                                               )
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()
        # - Now plot up flight time series plots by variable
        vars2plot = 'NOy'

        # - Plot up NOy
        units = 'pptv'
        var2plot = 'NOy'
        ModVar2Plot = 'NOy'
        ObsVar2Plot = 'NOy'
        mod_label = mod_label_master
        mod_scale = 1E12
        ObsVar2PlotErr = None
        plt_errorbar = False
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_ts4ARNA_flt_period_obs(var2plot=var2plot, units=units,
                                   ObsVar2Plot=ObsVar2Plot,
                                   mod_scale=mod_scale,
                                   mod_label=mod_label,
                                   ModVar2Plot=ModVar2Plot,
                                   ylim=ylim,
                                   dfs_mod=dfs_mod[flight_ID],
                                   df_obs=df_obs_NOy,
                                   dfs_mod_period=dfs_mod_period[flight_ID],
                                   plt_alt_shadow=plt_alt_shadow,
                                   flight_ID=flight_ID,
                                   ObsVar2PlotErr=ObsVar2PlotErr,
                                   plt_errorbar=plt_errorbar,
                                   context=context,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up NOy (just gas phase species)
#         units = 'pptv'
#         var2plot = 'NOy-all'
#         ModVar2Plot = 'NOy-all'
#         ObsVar2Plot = 'NOy'
#         mod_label = mod_label_master
#         mod_scale = 1E12
#         ObsVar2PlotErr = None
#         plt_errorbar = False
# #        ylim = (-0.2, 1)
#         ylim = None
#         # Call timeseries plotter function
#         plt_ts4ARNA_flt_period_obs(var2plot=var2plot, units=units,
#                                    ObsVar2Plot=ObsVar2Plot,
#                                    mod_scale=mod_scale,
#                                    mod_label=mod_label,
#                                    ModVar2Plot=ModVar2Plot,
#                                    ylim=ylim,
#                                    dfs_mod=dfs_mod[flight_ID],
#                                    df_obs=df_obs_NOy,
#                                    dfs_mod_period=dfs_mod_period[flight_ID],
#                                    plt_alt_shadow=plt_alt_shadow,
#                                    flight_ID=flight_ID,
#                                    ObsVar2PlotErr=ObsVar2PlotErr,
#                                    plt_errorbar=plt_errorbar,
#                                    context=context,
#                                    )
#         # Save to PDF and close the plot
#         AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
#         plt.close()

        # - Plot up NOy-HNO3
        units = 'pptv'
        var2plot = 'NOy (model-HNO3)'
        ModVar2Plot = 'NOy-HNO3'
        ObsVar2Plot = 'NOy'
        mod_label = mod_label_master
        mod_scale = 1E12
        ObsVar2PlotErr = None
        plt_errorbar = False
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_ts4ARNA_flt_period_obs(var2plot=var2plot, units=units,
                                   ObsVar2Plot=ObsVar2Plot,
                                   mod_scale=mod_scale,
                                   mod_label=mod_label,
                                   ModVar2Plot=ModVar2Plot,
                                   ylim=ylim,
                                   dfs_mod=dfs_mod[flight_ID],
                                   df_obs=df_obs_NOy,
                                   dfs_mod_period=dfs_mod_period[flight_ID],
                                   plt_alt_shadow=plt_alt_shadow,
                                   flight_ID=flight_ID,
                                   ObsVar2PlotErr=ObsVar2PlotErr,
                                   plt_errorbar=plt_errorbar,
                                   context=context,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up 'NOy-HNO3-PAN'
        units = 'pptv'
        var2plot = 'NOy (model-HNO3-PAN)'
        ModVar2Plot = 'NOy-HNO3'
        ObsVar2Plot = 'NOy'
        mod_label = mod_label_master
        mod_scale = 1E12
        ObsVar2PlotErr = None
        plt_errorbar = False
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_ts4ARNA_flt_period_obs(var2plot=var2plot, units=units,
                                   ObsVar2Plot=ObsVar2Plot,
                                   mod_scale=mod_scale,
                                   mod_label=mod_label,
                                   ModVar2Plot=ModVar2Plot,
                                   ylim=ylim,
                                   dfs_mod=dfs_mod[flight_ID],
                                   df_obs=df_obs_NOy,
                                   dfs_mod_period=dfs_mod_period[flight_ID],
                                   plt_alt_shadow=plt_alt_shadow,
                                   flight_ID=flight_ID,
                                   ObsVar2PlotErr=ObsVar2PlotErr,
                                   plt_errorbar=plt_errorbar,
                                   context=context,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up 'NOy-HNO3-PAN'
        units = 'pptv'
        var2plot = 'NOy (NOx+HONO+NIT(s))'
        ModVar2Plot = 'NOy-Limited'
        ObsVar2Plot = 'NOy'
        mod_label = mod_label_master
        mod_scale = 1E12
        ObsVar2PlotErr = None
        plt_errorbar = False
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_ts4ARNA_flt_period_obs(var2plot=var2plot, units=units,
                                   ObsVar2Plot=ObsVar2Plot,
                                   mod_scale=mod_scale,
                                   mod_label=mod_label,
                                   ModVar2Plot=ModVar2Plot,
                                   ylim=ylim,
                                   dfs_mod=dfs_mod[flight_ID],
                                   df_obs=df_obs_NOy,
                                   dfs_mod_period=dfs_mod_period[flight_ID],
                                   plt_alt_shadow=plt_alt_shadow,
                                   flight_ID=flight_ID,
                                   ObsVar2PlotErr=ObsVar2PlotErr,
                                   plt_errorbar=plt_errorbar,
                                   context=context,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up HONO
        try:
            units = 'pptv'
            var2plot = 'HONO'
            ModVar2Plot = 'HNO2'
            ObsVar2Plot = 'HONO'
            mod_label = mod_label_master
            mod_scale = 1E12
            ObsVar2PlotErr = None
            plt_errorbar = False
    #        ylim = (-0.2, 1)
            ylim = None
            title_str = "Timeseries of '{}' ({}) during flight '{}' - {}"
            ext_str = 'ToF CIMS'
            title = title_str.format(var2plot, units, flight_ID, ext_str)
            # Call timeseries plotter function
            plt_ts4ARNA_flt_period_obs(var2plot=var2plot,
                                       units=units,
                                       ObsVar2Plot=ObsVar2Plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       ModVar2Plot=ModVar2Plot,
                                       ylim=ylim,
                                       dfs_mod=dfs_mod[flight_ID],
                                       df_obs=df_CIMS_period,
                                       dfs_mod_period=dfs_mod_period[flight_ID],
                                       plt_alt_shadow=plt_alt_shadow,
                                       flight_ID=flight_ID,
                                       ObsVar2PlotErr=ObsVar2PlotErr,
                                       plt_errorbar=plt_errorbar,
                                       context=context,
                                       title=title,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))
        # - Plot up HNO3
        try:
            units = 'pptv'
            var2plot = 'HNO3'
            ModVar2Plot = 'HNO3'
            ObsVar2Plot = 'HNO3'
            mod_label = mod_label_master
            mod_scale = 1E12
            ObsVar2PlotErr = None
            plt_errorbar = False
    #        ylim = (-0.2, 1)
            ylim = None
            title_str = "Timeseries of '{}' ({}) during flight '{}' - {}"
            ext_str = 'ToF CIMS'
            title = title_str.format(var2plot, units, flight_ID, ext_str)
            # Call timeseries plotter function
            plt_ts4ARNA_flt_period_obs(var2plot=var2plot,
                                       units=units,
                                       ObsVar2Plot=ObsVar2Plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       ModVar2Plot=ModVar2Plot,
                                       ylim=ylim,
                                       dfs_mod=dfs_mod[flight_ID],
                                       df_obs=df_CIMS_period,
                                       dfs_mod_period=dfs_mod_period[flight_ID],
                                       plt_alt_shadow=plt_alt_shadow,
                                       flight_ID=flight_ID,
                                       ObsVar2PlotErr=ObsVar2PlotErr,
                                       plt_errorbar=plt_errorbar,
                                       context=context,
                                       title=title,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))
        # - Plot up HONO (from nitrate rack)
        try:
            units = 'pptv'
            var2plot = 'HONO'
            ObsVar2Plot = 'HONO_pptV'
            ModVar2Plot = 'HNO2'
            mod_label = mod_label_master
            mod_scale = 1E12
#            yscale = 'log'
#             yscale = 'linear'
#             if yscale == 'log':
#                 ylim = (0.3, 400)
#             if yscale == 'linear':
#                 ylim = (-20, 70)
            ylim = None
            ObsVar2PlotErr = None
            plt_errorbar = False
            title_str = "Timeseries of '{}' ({}) during flight '{}' - {}"
            ext_str = 'FAAM Nitrate rack data'
            title = title_str.format(var2plot, units, flight_ID, ext_str)
            # Call timeseries plotter function
            plt_ts4ARNA_flt_period_obs(var2plot=var2plot,
                                       units=units,
                                       ObsVar2Plot=ObsVar2Plot,
                                       mod_scale=mod_scale,
                                       mod_label=mod_label,
                                       ModVar2Plot=ModVar2Plot,
                                       ylim=ylim,
                                       dfs_mod=dfs_mod[flight_ID],
                                       df_obs=df_FAAM_period,
                                       dfs_mod_period=dfs_mod_period[flight_ID],
                                       plt_alt_shadow=plt_alt_shadow,
                                       flight_ID=flight_ID,
                                       ObsVar2PlotErr=ObsVar2PlotErr,
                                       plt_errorbar=plt_errorbar,
                                       context=context,
                                       title=title,
                                       )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except KeyError:
            pstr = "WARNING: '{}' not plotted for '{}' - not in DataDrame"
            print(pstr.format(var2plot, flight_ID))

        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_ts_comp4ARNA_flights_filters(dpi=320, show_plot=False,
                                     flight_nums=[],
                                     res='4x5', RunSet=None,
                                     inc_GEOSChem=False,
                                     just_plot_GEOS_Chem=False,
                                     LatVar='model-lat',
                                     LonVar='model-lon',
                                     plt_alt_shadow=False,
                                     debug=False, context='paper'):
    """
    Plot up timeseries comparisons between filter samples and model data
    """
    import seaborn as sns
    # Which flights to plot? - Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [218, 219, 220, 221, 222, 223, 224, 225, ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod_CF = {}
    for flight_ID in flight_IDs:
        df = get_GEOSCF4flightnum(flight_ID=flight_ID)
        # Add the derived variables to the dataframe
        df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
        df = add_deriv_vars2df(df=df)
        dfs_mod_CF[flight_ID] = df
        mod_label_master = 'GEOS-CF'

    # Model - GEOS-Chem (offline)
    if inc_GEOSChem:
        mod_label_master = RunSet
        dfs_mod_GC = {}
        for flight_ID in flight_IDs:
            dfs = get_GEOSChem4flightnum(flight_ID=flight_ID, res=res,
                                         RunSet=RunSet,)
            for key in dfs.keys():
                df = dfs[key]
                df = add_derived_FAAM_flags2df4flight(df=df,
                                                      flight_ID=flight_ID)
                df = add_deriv_vars2df(df=df)
                dfs[key] = df
            dfs_mod_GC[flight_ID] = dfs

#     dfs_mod = {}
#     for flight_ID in flight_IDs:
#         if debug:
#             print(flight_ID)
#         df = get_GEOSCF4flightnum(flight_ID=flight_ID)
#         # Add the derived variables to the dataframe
#         df = add_derived_FAAM_flags2df4flight(df=df, flight_ID=flight_ID)
#         df = add_deriv_vars2df(df=df)
#         dfs_mod[flight_ID] = df

    # KLUDGE -  just plot GEOS-Chem for now
    # (previously just plotted just GEOS-CF)
    print('WARNING: TODO setup to plot multiple model runs')
    dfs_mod = dfs_mod_GC

    # Observations
    dfs_obs = get_filters_data4flight()
    dfs_obs = add_secs2duplicate_index_values(dfs_obs)
    # Only consider the model for the observation (filter) period
    dfs_mod_period = {}
    for flight_ID in flight_IDs:
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs.loc[dfs_obs['Flight'] == flight_ID, :]
#        df_mod = dfs_mod[flight_ID]
        # Kludge! force use of single set of model output for now...
        print('WARNING: TODO setup to plot multiple model runs')
        ModelVarName = list(dfs_mod[flight_ID].keys())[0]
        dfs_mod4flight = dfs_mod[flight_ID]  # [ModelVarName]
        for key in dfs_mod4flight.keys():
            df_mod = dfs_mod4flight[key]
            df_mod = only_use_filter_times(df=df_mod,
                                           FILTERdf=df_obs,
                                           flight_ID=flight_ID)
            dfs_mod4flight[key] = df_mod
        dfs_mod_period[flight_ID] = dfs_mod4flight

    # Convert observation units into model units
    # unit on recipt were 'nanomoles/m3', which were updated to ug/m3
    # model units? 'pptv'
    NewUnits = 'ug_m-3'
    UncertaintyStr = 'Total_{}_uncertainty_{}'
    SpeciesStr = 'Total_{}_{}'
    # NOTE: TODO, include uncertainty for plots
    ['Cl', 'NO3', 'NO2', 'SO4', 'C2O4', 'Na', 'K', 'NH4', 'Ca', 'Mg']

    NIT_obs_var = SpeciesStr.format('NO3', 'ppt')
    SO4_obs_var = SpeciesStr.format('SO4', 'ppt')
    Cl_obs_var = SpeciesStr.format('Cl', 'ppt')
    NO2_obs_var = SpeciesStr.format('NO2', 'ppt')
    C2O4_obs_var = SpeciesStr.format('C2O4', 'ppt')
    NH4_obs_var = SpeciesStr.format('NH4', 'ppt')

    # No need to convert as using ppt values (provided in new dataset)
    obs2ModName = {
        NIT_obs_var: 'NIT', SO4_obs_var: 'SO4', Cl_obs_var: 'Cl',
        NO2_obs_var: 'NO2'
    }
#     for var2use in [NIT_obs_var, SO4_obs_var, NH4_var2use]:
#         spec = obs2ModName[var2use]
#         data = dfs_obs[var2use].values
#         dfs_obs[var2use] = AC.convert_ug_per_m3_2_ppbv(data, spec=spec)*1E3

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs.loc[dfs_obs['Flight'] == flight_ID, :]
#        df_mod = dfs_mod[flight_ID]
        # Kludge! force use of single set of model output for now...
        print('WARNING: TODO setup to plot multiple model runs')
        df_mod = dfs_mod[flight_ID][ModelVarName]
        print(df_mod)
#        df_mod_period = dfs_mod_period[flight_ID]
        # add the dust and SLR flags to the core dataframe
        df_obs = add_derived_FAAM_flags2df4flight(df=df_obs,
                                                  flight_ID=flight_ID)
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_timeseries_flighttrack_{}_FILTERS'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        print(df_mod)
#        print(df_mod.columns)
#        print(df_mod.head())
        plt_flightpath_spatially_over_CVAO(df=df_mod, flight_ID=flight_ID,
                                           LatVar=LatVar, LonVar=LonVar,)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()
        # - Now plot up flight time series plots by variable
        vars2plot = 'SO4', 'NIT'

        # - Plot up nitrate
        units = 'pptv'
        var2plot = 'NO3'
        ModVar2Plot = 'NIT.total'
        ObsVar2Plot = NIT_obs_var
        mod_label = mod_label_master
        mod_scale = 1E12
        ObsVar2PlotErr = None
        plt_errorbar = False
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_ts4ARNA_flt_period_obs(var2plot=var2plot, units=units,
                                   ObsVar2Plot=ObsVar2Plot,
                                   mod_scale=mod_scale,
                                   mod_label=mod_label,
                                   ModVar2Plot=ModVar2Plot,
                                   ylim=ylim,
                                   dfs_mod=dfs_mod[flight_ID],
                                   df_obs=df_obs,
                                   dfs_mod_period=dfs_mod_period[flight_ID],
                                   plt_alt_shadow=plt_alt_shadow,
                                   flight_ID=flight_ID,
                                   ObsVar2PlotErr=ObsVar2PlotErr,
                                   plt_errorbar=plt_errorbar,
                                   context=context,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up sulfate
        units = 'pptv'
        var2plot = 'SO4'
        ModVar2Plot = 'SO4.total'
        ObsVar2Plot = SO4_obs_var
        mod_label = mod_label_master
        mod_scale = 1E12
        ObsVar2PlotErr = None
        plt_errorbar = False
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_ts4ARNA_flt_period_obs(var2plot=var2plot, units=units,
                                   ObsVar2Plot=ObsVar2Plot,
                                   mod_scale=mod_scale,
                                   mod_label=mod_label,
                                   ModVar2Plot=ModVar2Plot,
                                   ylim=ylim,
                                   dfs_mod=dfs_mod[flight_ID],
                                   df_obs=df_obs,
                                   dfs_mod_period=dfs_mod_period[flight_ID],
                                   plt_alt_shadow=plt_alt_shadow,
                                   flight_ID=flight_ID,
                                   ObsVar2PlotErr=ObsVar2PlotErr,
                                   plt_errorbar=plt_errorbar,
                                   context=context,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up Ammonium
        units = 'pptv'
        var2plot = 'NH4'
        ModVar2Plot = 'NH4'
        ObsVar2Plot = NH4_obs_var
        mod_label = mod_label_master
        mod_scale = 1E12
        ObsVar2PlotErr = None
        plt_errorbar = False
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_ts4ARNA_flt_period_obs(var2plot=var2plot, units=units,
                                   ObsVar2Plot=ObsVar2Plot,
                                   mod_scale=mod_scale,
                                   mod_label=mod_label,
                                   ModVar2Plot=ModVar2Plot,
                                   ylim=ylim,
                                   dfs_mod=dfs_mod[flight_ID],
                                   df_obs=df_obs,
                                   dfs_mod_period=dfs_mod_period[flight_ID],
                                   plt_alt_shadow=plt_alt_shadow,
                                   flight_ID=flight_ID,
                                   ObsVar2PlotErr=ObsVar2PlotErr,
                                   plt_errorbar=plt_errorbar,
                                   context=context,
                                   )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_ts_comp4ARNA_flights_SWAS(dpi=320, show_plot=False,
                                  flight_nums=[],
                                  context='paper', debug=False):
    """
    Plot up timeseries comparisons between SWAS observations and model data
    """
    import seaborn as sns
    # Which flights to plot?
    # Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [
            217, 218, 219, 220, 221, 222, 223, 224, 225,
        ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model
    dfs_mod = {}
    for flight_ID in flight_IDs:
        if debug:
            print(flight_ID)
        df = get_GEOSCF4flightnum(flight_ID=flight_ID)
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
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()
        # - Now plot up flight time series plots by variable
        vars2plot = 'ALD2', 'ACET', 'C2H6', 'C3H8'

        # - Plot up acetone
        units = 'ppt'
        var2plot = 'ACET'
        ModVar2Plot = 'ACET'
        ObsVar2Plot = map_SWAS_var2GEOS_var(ModVar2Plot, invert=True)
        mod_label = 'GEOS-CF'
        mod_scale = 1E9
        ObsVar2PlotErr = ObsVar2Plot+'_uncertainty'
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_ts4ARNA_flt_point_obs(var2plot=var2plot, units=units,
                                  ObsVar2Plot=ObsVar2Plot,
                                  mod_scale=mod_scale,
                                  mod_label=mod_label,
                                  ModVar2Plot=ModVar2Plot,
                                  ylim=ylim,
                                  df_mod=df_mod, df_obs=df_obs,
                                  flight_ID=flight_ID,
                                  ObsVar2PlotErr=ObsVar2PlotErr,
                                  plt_errorbar=True,
                                  context=context,
                                  )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up acetaldehyde
        units = 'ppt'
        var2plot = 'ALD2'
        ModVar2Plot = 'ALD2'
        ObsVar2Plot = map_SWAS_var2GEOS_var(ModVar2Plot, invert=True)
        mod_label = 'GEOS-CF'
        mod_scale = 1E9
        ObsVar2PlotErr = ObsVar2Plot+'_uncertainty'
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_ts4ARNA_flt_point_obs(var2plot=var2plot, units=units,
                                  ObsVar2Plot=ObsVar2Plot,
                                  mod_scale=mod_scale,
                                  mod_label=mod_label,
                                  ModVar2Plot=ModVar2Plot,
                                  ylim=ylim,
                                  df_mod=df_mod, df_obs=df_obs,
                                  flight_ID=flight_ID,
                                  ObsVar2PlotErr=ObsVar2PlotErr,
                                  plt_errorbar=True,
                                  context=context,
                                  )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up ethane
        units = 'ppt'
        var2plot = 'C2H6'
        ModVar2Plot = 'C2H6'
        ObsVar2Plot = map_SWAS_var2GEOS_var(ModVar2Plot, invert=True)
        mod_label = 'GEOS-CF'
        mod_scale = 1E9
        ObsVar2PlotErr = ObsVar2Plot+'_uncertainty'
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_ts4ARNA_flt_point_obs(var2plot=var2plot, units=units,
                                  ObsVar2Plot=ObsVar2Plot,
                                  mod_scale=mod_scale,
                                  mod_label=mod_label,
                                  ModVar2Plot=ModVar2Plot,
                                  ylim=ylim,
                                  df_mod=df_mod, df_obs=df_obs,
                                  flight_ID=flight_ID,
                                  ObsVar2PlotErr=ObsVar2PlotErr,
                                  plt_errorbar=True,
                                  context=context,
                                  )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up propane
        units = 'ppt'
        var2plot = 'C3H8'
        ModVar2Plot = 'C3H8'
        ObsVar2Plot = map_SWAS_var2GEOS_var(ModVar2Plot, invert=True)
        mod_label = 'GEOS-CF'
        mod_scale = 1E9
        ObsVar2PlotErr = ObsVar2Plot+'_uncertainty'
#        ylim = (-0.2, 1)
        ylim = None
        # Call timeseries plotter function
        plt_ts4ARNA_flt_point_obs(var2plot=var2plot, units=units,
                                  ObsVar2Plot=ObsVar2Plot,
                                  mod_scale=mod_scale,
                                  mod_label=mod_label,
                                  ModVar2Plot=ModVar2Plot,
                                  ylim=ylim,
                                  df_mod=df_mod, df_obs=df_obs,
                                  flight_ID=flight_ID,
                                  ObsVar2PlotErr=ObsVar2PlotErr,
                                  plt_errorbar=True,
                                  context=context,
                                  )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def plt_ts_comp4ARNA_flights_PHYSICAL_VARS(dpi=320, show_plot=False,
                                           inc_GEOSChem=False,
                                           RunSet=None, res='4x5',
                                           just_plot_GEOS_Chem=False,
                                           flight_nums=[],
                                           context='paper'):
    """
    Plot up timeseries comparisons between physical variables and model data
    """
    import seaborn as sns
    # Which flights to plot?
    # Just use non-transit ARNA flights
    if len(flight_nums) == 0:
        flight_nums = [
            217, 218, 219, 220, 221, 222, 223, 224, 225,
        ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        dfs_obs[flight_ID] = get_FAAM_core4flightnum(flight_ID=flight_ID)

    # Model
    dfs_mod_CF = {}
    for flight_ID in flight_IDs:
        dfs_mod_CF[flight_ID] = get_GEOSCF4flightnum(flight_ID=flight_ID)
    if inc_GEOSChem:
        #        RunSet='MERRA2-0.5-initial'
        #        res='0.5x0.625'
        #        RunSet='MERRA2-BC'
        #        res='4x5'
        #        RunSet='FP-Nest'
        #        res='0.25x0.3125'
        dfs_mod_GC = {}
        for flight_ID in flight_IDs:
            dfs_mod_GC[flight_ID] = get_GEOSChem4flightnum(flight_ID=flight_ID,
                                                           res=res,
                                                           RunSet=RunSet,)

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        df_mod_CF = dfs_mod_CF[flight_ID]
        if inc_GEOSChem:
            if just_plot_GEOS_Chem:
                dfs_mod = dfs_mod_GC[flight_ID]
                mod_label_master = RunSet

            else:
                dfs_mod_GC4flight = dfs_mod_GC[flight_ID]
                dfs_mod = {'GEOS-CF': df_mod_CF}
                for key in list(dfs_mod_GC4flight.keys()):
                    dfs_mod[key] = dfs_mod_GC4flight[key]
                mod_label_master = 'GEOS-CF'

        else:
            mod_label_master = 'GEOS-CF'
            dfs_mod = {'GEOS-CF': df_mod_CF}

        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_timeseries_flighttrack_{}_PHYSICAL_VARS'
        savetitle = savetitle.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up ROLL_GIN
        try:
            units = 'ADD THIS'
    #        VarName = 'PCAS2CON'
    #        FlagName = 'PCAS2_FLAG'
            var2plot = 'ROLL_GIN'
            ObsVar2Plot = 'ROLL_GIN'
            ModVar2Plot = None
            mod_label = mod_label_master
            mod_scale = None
    #        ylim = (10, 100)
    #        ylim = None
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               #                                   ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               plt_legend=False,
                               context=context,
                               )
            # Colour in SLRs
            ax = plt.gca()
            colour_plot_background_by_bool(df=df_obs, ax=ax, bool2use='IS_SLR',
                                           color='blue', alpha=0.25,
                                           label='SLR')
            fig = plt.gcf()
            fig.legend(loc='best', bbox_to_anchor=(1, 1),
                       bbox_transform=ax.transAxes)
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot ROLL_GIN')

        # - Plot up VELD_GIN
        try:
            units = 'ADD THIS'
    #        VarName = 'PCAS2CON'
    #        FlagName = 'PCAS2_FLAG'
            var2plot = 'VELD_GIN'
            ObsVar2Plot = 'VELD_GIN'
            ModVar2Plot = None
            mod_label = mod_label_master
            mod_scale = None
    #        ylim = (10, 100)
            ylim = None
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               plt_legend=False,
                               context=context,
                               )
            # Colour in SLRs
            ax = plt.gca()
            colour_plot_background_by_bool(df=df_obs, ax=ax, bool2use='IS_SLR',
                                           color='blue', alpha=0.25,
                                           label='SLR')
            fig = plt.gcf()
            fig.legend(loc='best', bbox_to_anchor=(1, 1),
                       bbox_transform=ax.transAxes)

            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot VELD_GIN')

        # - Plot up PCASP
        try:
            units = '# cm$^{-1}$'
    #        VarName = 'PCAS2CON'
    #        FlagName = 'PCAS2_FLAG'
            var2plot = 'PCAS'
            ObsVar2Plot = 'PCAS2CON'
            ModVar2Plot = None
            mod_label = mod_label_master
            mod_scale = None
    #        ylim = (10, 100)
            ylim = None
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up PCASP surface area
        try:
            units = 'x10$^{-12}$ m${^3}$/cm$^{-3}$'
    #        VarName = 'PCAS2CON'
    #        FlagName = 'PCAS2_FLAG'
            var2plot = 'PCASP (total surface area)'
            ObsVar2Plot = 'PCASP-total-surface'
            ModVar2Plot = None
            mod_label = mod_label_master
            obs_scale = 1E12
            mod_scale = None
    #        ylim = (10, 100)
            ylim = None
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               obs_scale=obs_scale,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up PCASP surface area
        try:
            units = 'x10$^{-9}$ m${^3}$/cm$^{-3}$'
            var2plot = 'CDP (total surface area)'
            ObsVar2Plot = 'CDP-total-surface'
            ModVar2Plot = None
            mod_label = mod_label_master
            obs_scale = 1E9
            mod_scale = None
    #        ylim = (10, 100)
            ylim = None
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               obs_scale=obs_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up temperature
        units = '$^{\circ}$C'
        var2plot = 'Temperature'
        ObsVar2Plot = 'TAT_DI_R'
        ModVar2Plot = 'T'
        mod_label = mod_label_master
        mod_scale = 1
        ylim = (-25, 35)
        obs_adjustby = -273.15
        # Call timeseries plotter function
        plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                           ObsVar2Plot=ObsVar2Plot,
                           mod_scale=mod_scale,
                           mod_label=mod_label,
                           ModVar2Plot=ModVar2Plot,
                           ylim=ylim,
                           dfs_mod=dfs_mod,
                           df_obs=df_obs,
                           flight_ID=flight_ID,
                           obs_adjustby=obs_adjustby,
                           context=context,
                           )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up Eastward wind
        units = 'm s$^{-1}$'
        var2plot = 'Eastward wind'
        ObsVar2Plot = 'U_C'
        ModVar2Plot = 'U'
        mod_label = mod_label_master
        mod_scale = 1
        ylim = (-25, 25)
        # Call timeseries plotter function
        plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                           ObsVar2Plot=ObsVar2Plot,
                           mod_scale=mod_scale,
                           mod_label=mod_label,
                           ModVar2Plot=ModVar2Plot,
                           ylim=ylim,
                           dfs_mod=dfs_mod,
                           df_obs=df_obs,
                           flight_ID=flight_ID,
                           context=context,
                           )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up Northward wind
        units = 'm s$^{-1}$'
        var2plot = 'Northward wind'
        ObsVar2Plot = 'V_C'
        ModVar2Plot = 'V'
        mod_label = mod_label_master
        mod_scale = 1
        ylim = (-25, 25)
        # Call timeseries plotter function
        plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                           ObsVar2Plot=ObsVar2Plot,
                           mod_scale=mod_scale,
                           mod_label=mod_label,
                           ModVar2Plot=ModVar2Plot,
                           ylim=ylim,
                           dfs_mod=dfs_mod,
                           df_obs=df_obs,
                           flight_ID=flight_ID,
                           context=context,
                           )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up Latitude
        try:
            units = '$^{\circ}$N'
            var2plot = 'Latitude'
            ObsVar2Plot = 'LAT_GIN'
            ModVar2Plot = 'model-lat'
            mod_label = mod_label_master
            mod_scale = 1
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up Longitude
        try:
            units = '$^{\circ}$E'
            var2plot = 'Longitude'
            ObsVar2Plot = 'LON_GIN'
            ModVar2Plot = 'model-lon'
            mod_label = mod_label_master
            mod_scale = 1
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot {}'.format(var2plot))

        # - Plot up altitude
        try:
            units = 'hPa'
            var2plot = 'Altitude'
            if mod_label_master == 'GEOS-CF':
                ModVar2Plot = 'model-lev'
            else:
                ModVar2Plot = 'PRESS'
            ObsVar2Plot = 'PS_RVSM'  # Use pressure measurement
#            ObsVar2Plot = 'ALT_GIN' # Use GPS altitude?
#            vals = AC.hPa_to_Km( df_obs['ALT_GIN'].values/1E3, reverse=True )
            mod_label = mod_label_master
            mod_scale = 1
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               invert_yaxis=True,
                               plt_alt_shadow=False,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
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


def plt_ts_comp4ARNA_flights(dpi=320, inc_GEOSChem=False,
                             show_plot=False, context='talk',
                             RunSet=None, res='4x5',
                             just_plot_GEOS_Chem=False,
                             flight_nums=[],
                             ):
    """
    Plot up timeseries comparisons between core observations and model data
    """
    import seaborn as sns
    # Which flights to plot?
    if len(flight_nums) == 0:
        flight_nums = [
            # Just use non-transit ARNA flights
            #        216,
            #        217,
            218, 219, 220, 221, 222, 223, 224, 225,
        ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # - Loop by flight and retrieve the files as dataframes (mod + obs)
    # Model - GEOS-CF (online)
    dfs_mod_CF = {}
    for flight_ID in flight_IDs:
        dfs_mod_CF[flight_ID] = get_GEOSCF4flightnum(flight_ID=flight_ID)
    # Model - GEOS-Chem (offline)
    if inc_GEOSChem:
        #        RunSet='MERRA2-0.5-initial'
        #        res='0.5x0.625'
        #        RunSet='MERRA2-BC'
        #        res='4x5'
        #        RunSet='FP-Nest'
        #        res='0.25x0.3125'
        dfs_mod_GC = {}
        for flight_ID in flight_IDs:
            dfs_mod_GC[flight_ID] = get_GEOSChem4flightnum(flight_ID=flight_ID,
                                                           res=res,
                                                           RunSet=RunSet,)
    # Observations
    dfs_obs = {}
    for flight_ID in flight_IDs:
        dfs_obs[flight_ID] = get_FAAM_core4flightnum(flight_ID=flight_ID)

    # -  Now plot up
    for flight_ID in flight_IDs:
        print(flight_ID)
        # Get observations and model timeseries data as a DataFrame
        df_obs = dfs_obs[flight_ID]
        df_mod_CF = dfs_mod_CF[flight_ID]
        if inc_GEOSChem:
            if just_plot_GEOS_Chem:
                dfs_mod = dfs_mod_GC[flight_ID]
                mod_label_master = RunSet

            else:
                dfs_mod_GC4flight = dfs_mod_GC[flight_ID]
                dfs_mod = {'GEOS-CF': df_mod_CF}
                for key in list(dfs_mod_GC4flight.keys()):
                    dfs_mod[key] = dfs_mod_GC4flight[key]
                mod_label_master = 'GEOS-CF'

        else:
            mod_label_master = 'GEOS-CF'
            dfs_mod = {'GEOS-CF': df_mod_CF}
        # Setup PDF to save PDF plots to
        savetitle = 'ARNA_timeseries_flighttrack_{}'.format(flight_ID)
        pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)
        # - Plot up location of flights
        plt_flightpath_spatially_over_CVAO(df=df_obs, flight_ID=flight_ID)
        # Save to PDF
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Now plot up flight time series plots by variable
        # - Plot up carbon monoxide
        units = 'ppbv'
        var2plot = 'CO'
        ObsVar2Plot = 'CO_AERO'
        ModVar2Plot = 'CO'
        mod_label = mod_label_master
        mod_scale = 1E9
        ylim = (50, 400)
        # Call timeseries plotter function
        plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                           ObsVar2Plot=ObsVar2Plot,
                           mod_scale=mod_scale,
                           mod_label=mod_label,
                           ModVar2Plot=ModVar2Plot,
                           ylim=ylim,
                           dfs_mod=dfs_mod,
                           df_obs=df_obs,
                           flight_ID=flight_ID,
                           context=context,
                           )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up ozone
        units = 'ppbv'
        var2plot = 'Ozone'
        ObsVar2Plot = 'O3_TECO'
        ModVar2Plot = 'O3'
        mod_label = mod_label_master
        mod_scale = 1E9
        ylim = (10, 100)
        # Call timeseries plotter function
        plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                           ObsVar2Plot=ObsVar2Plot,
                           mod_scale=mod_scale,
                           mod_label=mod_label,
                           ModVar2Plot=ModVar2Plot,
                           ylim=ylim,
                           dfs_mod=dfs_mod,
                           df_obs=df_obs,
                           flight_ID=flight_ID,
                           context=context,
                           )
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

        # - Plot up NO2
        try:
            units = 'pptv'
            var2plot = 'NO2'
            ObsVar2Plot = 'NO2_pptV'
            ModVar2Plot = 'NO2'
            mod_label = mod_label_master
            mod_scale = 1E12
#            yscale = 'log'
            yscale = 'linear'
            if yscale == 'log':
                ylim = (0.3, 400)
            if yscale == 'linear':
                ylim = (-20, 400)
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               yscale=yscale,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot NO2')

        # - Plot up NO
        try:
            units = 'pptv'
            var2plot = 'NO'
            ObsVar2Plot = 'NO_pptV'
            ModVar2Plot = 'NO'
            mod_label = mod_label_master
            mod_scale = 1E12
#            yscale = 'log'
            yscale = 'linear'
            if yscale == 'log':
                ylim = (0.3, 400)
            if yscale == 'linear':
                ylim = (-10, 200)
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               yscale=yscale,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot NO')

        # - Plot up NOx
        try:
            units = 'pptv'
            var2plot = 'NOx'
            ModVar2Plot = 'NOx'
            ObsVar2Plot = ModVar2Plot
            mod_label = mod_label_master
            mod_scale = 1E12
#            yscale = 'log'
            yscale = 'linear'
            if yscale == 'log':
                ylim = (0.3, 400)
            if yscale == 'linear':
                ylim = (-20, 400)
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               yscale=yscale,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot NOx')

        # - Plot up HNO2
        try:
            units = 'pptv'
            var2plot = 'HONO'
            ObsVar2Plot = 'HONO_pptV'
            ModVar2Plot = 'HNO2'
            mod_label = mod_label_master
            mod_scale = 1E12
#            yscale = 'log'
            yscale = 'linear'
            if yscale == 'log':
                ylim = (0.3, 400)
            if yscale == 'linear':
                ylim = (-20, 70)
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod=dfs_mod,
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               yscale=yscale,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
            plt.close()
        except:
            print('Failed to plot HONO')

        # - Plot up NOx / 'PCASP-total-surface'
        try:
            #
            ObsVar2Plot = 'NOx/PCASP'
            ratio = df_obs['NOx'] / df_obs['PCASP-total-surface']
            df_obs[ObsVar2Plot] = ratio
            # Then plot
            units = 'ratio'
            var2plot = ObsVar2Plot
            ModVar2Plot = None
            ObsVar2Plot = ObsVar2Plot
            mod_label = mod_label_master
            mod_scale = None
            yscale = 'linear'
            # Call timeseries plotter function
            plt_ts4ARNA_flight(var2plot=var2plot, units=units,
                               ObsVar2Plot=ObsVar2Plot,
                               mod_scale=mod_scale,
                               mod_label=mod_label,
                               ModVar2Plot=ModVar2Plot,
                               ylim=ylim,
                               dfs_mod={mod_label: df_mod},
                               df_obs=df_obs,
                               flight_ID=flight_ID,
                               yscale=yscale,
                               context=context,
                               )
            # Save to PDF and close the plot
            AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
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
    marker = 'o'
    s = 2
    alpha = 0.75
    cmap_list = AC.get_CB_color_cycle()
    #
    flight_IDs = list(sorted(d.keys()))
    for flight_ID_n, flight_ID in enumerate(flight_IDs):
        df = d[flight_ID]
        # Save the data as a csv file?
        save2csv = True
        if save2csv:
            df.to_csv('TEMP_{}.csv'.format(flight_ID))
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
    d = get_analysis_region('local_CVAO_area')
    extents = (d['x0'], d['x1'], d['y0'], d['y1'])
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
    plt.suptitle(u'ARNA campaign flights (excluding transits)')
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
    if isinstance(ax, type(None)):
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
        cbar_kwargs = {'cmap': cmap, 'extend': 'both', 'pad': 0.075, }
    else:
        cbar_kwargs = {
            'ticks': ticks, 'cmap': cmap, 'extend': 'both', 'pad': 0.075,
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
        cbar_kwargs = {'cmap': cmap, 'extend': 'both', }
    else:
        cbar_kwargs = {'ticks': ticks, 'cmap': cmap, 'extend': 'both', }
    # Setup a colour bar axis
    cbar_ax = divider.append_axes("right", "2%", pad="7%")
    # using pcolormesh
    im = ds[var2plot2].plot.pcolormesh(x=LonVar, y=LevVar, ax=ax,
                                       vmin=vmin2, vmax=vmax2,
                                       zorder=1, alpha=alpha, cmap=cmap,
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
        max_alt = d['max_alt'] / 1E3
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
               xmin=(min_lon - xlim[0]) / xrange,
               xmax=(max_lon - xlim[0]) / xrange,
               linewidth=3.0)

    # Add locations for airports
#    locs2plot = ['DSS', 'RAI', 'VXE', 'LPA', 'LIS', 'CDG']
    locs2plot = [
        'Praia Airport',  'Sao Vicente Airport', 'Dakar',
        #    'Gran Canaria Airport', 'Lisbon Airport',  'Paris (Charles de Gaulle) Airport'
    ]
    for n, loc_ in enumerate(locs2plot):
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
                buffer = buffer * 3*1.5

            # Add label for the airport
            ax.text(lon_, base+buffer, '{}'.format(loc_), fontsize=10,
                    alpha=0.5,
                    horizontalalignment='center')

    #        ax.annotate(loc_, xy=(lat_, 5), xytext=(lat_+buffer, 5+2),
    #            arrowprops=dict(facecolor='black', shrink=0.05))

    # Add lines for kft heights
    if convert2kft:
        hPa_heights = [1000, 900, 800, 700, 600, 500]
        km_heights = AC.hPa_to_Km(hPa_heights)
        kft_heights = [i*m2kft for i in km_heights]
        for n, height_ in enumerate(kft_heights):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-3.5, height_, '{:.0f} hPa'.format(hPa_heights[n]),
                    fontsize=10,
                    alpha=0.5)
    else:
        # Add lines for kft heights
        kft_heights = [20000, 15000, 10000, 5000]
        m_heights = [i/m2kft/1E3 for i in kft_heights]
        for n, height_ in enumerate(m_heights):
            ax.axhline(y=height_, linestyle='--', alpha=0.5, color='grey',
                       zorder=100,
                       linewidth=1.0)
            # Add label for heights
            ax.text(xlim[1]-5, height_,
                    '{:.0f} kft'.format(kft_heights[n]/1E3),
                    fontsize=10, alpha=0.5)

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


def plt_highres_modelling_region(ds=None, var2use='NOy', LatVar='lat',
                                 LonVar='lon',
                                 plot_blank_data=False, rm_colourbar=True,
                                 add_box4highres_region=True,
                                 projection=ccrs.PlateCarree):
    """
    Plot a map to show ARNA campaign region
    """
    # Reset sns for spatial plot
    sns.reset_orig()
    # Plot some blank data?
    if plot_blank_data:
        if isinstance(ds, type(None)):
            #            NASA_data = get_local_folder('NASA_data')
            #            folder = NASA_data + 'nature_run/O3/'
            #            filename = 'nature_run_lev_72_res_0.125_spec_O3_2014_034_ctm.nc'
            folder = '/mnt/lustre/groups/chem-acm-2018/earth0_data/'
            folder += 'Data4Papers/Data_Full/'
            folder += '201609_ACP_Sherwen_GC_Present_day_halogens_Externally_shared_2x25_files/HAL/'
            filename = 'ctm.nc'
            ds = xr.open_dataset(folder+filename)
            var2use = 'IJ_AVG_S__O3'
            ds = ds[[var2use]].mean(dim='time')
            ds[var2use][:] = np.NaN
            # Select a single level
            ds = ds.sel(model_level_number=ds.model_level_number[0])
            del ds['model_level_number']
            ds[[var2use]].attrs = {}
            # Rename lat and long
            ds = ds.rename({'latitude': LatVar, 'longitude': LonVar})
    # - Extents to use
    # area of extracted data from OPeNDAP
#    d = get_analysis_region('OPeNDAP_download_area')
    # Context area - local area
    d = get_analysis_region('context_area')
    # Local area analysed as Cape Verde
#     d = get_analysis_region('local_CVAO_area')
    # extract dictionary to local variables
    x0, x1, y0, y1 = d['x0'], d['x1'], d['y0'], d['y1']
    # - Select the data
    # Set values region
    bool1 = ((ds.lon >= x0) & (ds.lon <= x1)).values
    bool2 = ((ds.lat >= y0) & (ds.lat <= y1)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    # Plot the data
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection=projection(), aspect='auto',
                         alpha=0.5)
    ds[var2use].plot.imshow(x=LonVar, y=LatVar, ax=ax,
                            transform=ccrs.PlateCarree())
    # Beautify the figure/plot
    ax.coastlines()
    # Add borders for countries
    ax.add_feature(cfeature.BORDERS, edgecolor='grey',
                   facecolor='none', zorder=50)
    # Also add minor islands (inc. Cape Verde)
#    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
#                                            edgecolor=None,
#                                            facecolor='none')
    # Force map perspective to be on West/Northern Africa
    ax.set_extent((x0, x1, y0, y1), projection())
    #
    if add_box4highres_region:
        d = get_analysis_region('model_nest_area')
        x0, x1, y0, y1 = d['x0'], d['x1'], d['y0'], d['y1']
        #
        ax = add_LinearRing2cartopy_map(ax=ax, x0=x0, x1=x1, y0=y0, y1=y1)

#    ax.set_global() # this will force a global perspective
    # Remove the colout mpa from the plot
    if rm_colourbar:
        im = ax.images
        cb = im[-1].colorbar
        cb.remove()
        del cb
    return fig, ax


def add_LinearRing2cartopy_map(ax=None, x0=15, x1=-35, y0=0, y1=34,
                               projection=ccrs.PlateCarree,
                               linewidth=5, zorder=100, linestyle='--',
                               edgecolor='green', facecolor='none',
                               ):
    """
    Add a LinearRing to cartopy plot
    """
    # Add rectangles for location of the high resolution grid
    # Now plot as a linear ring
    # Initial choice
#    x0=15, x1=-32, y0=0, y1=32.5,
    # 4x5 grid choice
#    x0=15, x1=-35, y0=0, y1=34,
    lons = (x0, x1, x1, x0)
    lats = (y0, y0, y1, y1)
    ring = LinearRing(list(zip(lons, lats)))
    ax.add_geometries([ring], projection(), facecolor=facecolor,
                      edgecolor=edgecolor, zorder=zorder,
                      linestyle=linestyle, linewidth=linewidth)
    return ax


def add_scatter_points2cartopy_ax(ax=None, projection=ccrs.PlateCarree,
                                  marker='o', s=2, lats=None, lons=None,
                                  alpha=0.5, color='red', zorder=999,
                                  label=None):
    """
    Add scatter points to cartopy axis
    """
    # Now scatter points on plot
    ax.scatter(lons, lats, color=color, s=s,
               marker=marker, alpha=alpha,
               transform=projection(),
               zorder=zorder, label=label)
    return ax


def plot_up_pNO3_photolysis_params():
    """
    Plot up the various pNO3 photolysis parameterisation
    """
    # - Simone's one copied from google colab notes
    no3 = np.arange(0., 10000, 0.01)
    f = 0.7*4000./(1.+0.7*no3)/2800

    photo = 2800*f+1.*(1-f)
    plt.loglog(no3*1e-9*(14+16+16+16)*1e6, photo, label='fbulk=1')

    photo = 2800*f+0.*(1-f)
    plt.loglog(no3*1e-9*(14+16+16+16)*1e6, photo, label='fbulk=0')

    x = np.array([1, 70])*1e-9*(14+16+16+16)*1e6
    plt.fill_between(x, 56, 2800, color='r', alpha=0.5)

    plt.legend()

    plt.xlabel('bulk [NO$_{3}^{-}$] (ug m$^{-3}$)')
    plt.ylabel('enhancement factor (f)')

    plt.show()

    # Check numbers
    no3 = np.arange(0, 70)
    f = 0.7*4000./(1.+0.7*no3)
    print(min(f), max(f))
    print(0.7*4000./(1+0.7*16.))

    # - And ones from Ye et al

    # - And...

    # What about other nitrate on other species? We have nitrate on sulfate, seasalt, and dust in the GEOS-Chem model.


def plt_quick_ts4df(df, vars2plot=None, savetitle=None, save2pdf=True,
                    dpi=320, context='paper', debug=False):
    """
    Plot up a quick timeseries plot for variables in DataFrame
    """
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context(context)
    # which variables to plot?
    if isinstance(vars2plot, type(None)):
        vars2plot = ['LAT', 'LON', 'PRESS', 'V', 'U', 'T', 'OH', 'O3', 'CO']
    if isinstance(savetitle, type(None)):
        savetitle = 'timeseries_plot_of_df'
    pdff = AC.plot2pdfmulti(title=savetitle, open=True, dpi=dpi)

    # - Plot entire data
    for var in vars2plot:
        try:
            units, scaleby, adjustby = AC.tra_unit(var, scale=True,
                                                   adjustment=True)
        except:
            units = ''
            scaleby = 1
            adjustby = 0
        if debug:
            print(var, units, scaleby, adjustby)
            df[var].dropna().describe()

        Y = (df[var].values * scaleby) + adjustby
        X = df.index.values
        plt.plot(X, Y)
        plt.ylabel("'{}' ({})".format(var, units))
        ax = plt.gca()
        labels = ax.get_xticklabels()
        ax.set_xticklabels(labels, rotation=45)
        # Beatify
        plt.title("'Timeseries plot of '{}' ({})".format(var, units))
        # Save to PDF and close the plot
        AC.plot2pdfmulti(pdff, savetitle, dpi=dpi, tight=True)
        plt.close()

    # Save to PDF
    if save2pdf:
        # Save entire pdf
        AC.plot2pdfmulti(pdff, savetitle, close=True, dpi=dpi)
        plt.close('all')


def mk_trifigure_NO3_JNIT_combination_plt(context='paper'):
    """
    Make a tri plot of nitrate, JNIT, and their product
    """
    #
    import seaborn as sns
    sns.set(color_codes=True)
    sns.set_context(context)
    # Plottting settings (shouldn't these be arguemnets ... )
    RunSet = 'ACID'
    res = '4x5'
    trop_limit = False
    dates2use = [datetime.datetime(2019, i+1, 1) for i in range(12)]
    ColorList = AC.get_CB_color_cycle()
    # Get data into a single dictionary
    RunDict = ar.get_dict_of_GEOSChem_model_output(res=res, RunSet=RunSet)
    NOxD = get_NOx_budget_ds_dict_for_runs(RunDict=RunDict,
                                           #                                           RunSet=RunSet, res=res,
                                           trop_limit=trop_limit,
                                           dates2use=dates2use)
    runs2use = [
        'Acid-4x5-Isotherm.v2',
        #        'Acid-4x5-J25',
        'Acid-4x5-J50',
        #        'Acid-4x5-Isotherm-BBx3'
        'Acid-4x5-J50-AfBBx3-NH3x3',
    ]
    for key in list(sorted(NOxD.keys())):
        if (key not in runs2use):
            del NOxD[key]

    # Add data to plot into a single DataFrame
    dfs = {}
    for key in NOxD.keys():
        ds = NOxD[key]
        ds = NOxD[key].sel(lev=ds.lev.values[0]).mean(dim='time')
        df = pd.DataFrame()

        # Get Nitrate concentration
        NITvar = 'SpeciesConc_NIT'
        data = ds[NITvar].values.flatten()
        df[NITvar] = data * 1E12

        # JNIT
        JNITvar = 'Jval_NIT'
        data1 = ds[JNITvar].values.flatten()
        JHNO3var = 'Jval_HNO3'
        data2 = ds[JHNO3var].values.flatten()
        JScaleVar = 'Jscale'
        df[JScaleVar] = data1/data2

        # JNIT * NIT
        ProdJNIT = 'ProdHNO2fromHvNIT-all'
        data = ds[ProdJNIT].values.flatten()
        df[ProdJNIT] = data

        dfs[key] = df

    # Plot up seperately
    if debug:
        for var in [NITvar, JScaleVar, ProdJNIT]:
            fig, ax = plt.subplots()
            for nKey, key in enumerate(dfs.keys()):
                sns.histplot(data=dfs[key], x=var, ax=ax, kde=True, label=key,
                             color=ColorList[nKey])
                plt.title(" '{}' in '{}'".format(var, key))
#                FileStr = 'ARNA_PDF_of_NIT_JNIT_and_product_{}_{}'
#                AC.save_plot(FileStr.format(key, var), tight=True)
                FileStr = 'ARNA_PDF_of_NIT_JNIT_and_product_{}'
                AC.save_plot(FileStr.format(var), tight=True)
                plt.close('all')

    # Plot up as a single plot
    fig, [ax1, ax2, ax3] = plt.subplots(nrows=3, ncols=1)
    for nKey, key in enumerate(dfs.keys()):
        sns.histplot(data=dfs[key], x=NITvar, ax=ax1, kde=True, label=key,
                     color=ColorList[nKey])
        ax1.legend()
        sns.histplot(data=dfs[key], x=JScaleVar, ax=ax2, kde=True, label=key,
                     color=ColorList[nKey])
        sns.histplot(data=dfs[key], x=ProdJNIT, ax=ax3, kde=True, label=key,
                     color=ColorList[nKey])
    AC.save_plot('ARNA_PDF_of_NIT_JNIT_and_product', tight=True)
    plt.close('all')


def mk_Andersen2021_figure_02(dpi=720, figsize=(7, 3), aspect=None, ms=None,
                              fontsize=None, labelsize=None, tight=True):
    """
    Make figure 2 for Andersen et al. (2021)

    Parameters
    ----------

    Returns
    -------
    (None)

    Notes
    -------
     - Original version author: Simone T. Andersen (Fri Nov  5th)
    """
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    # Boxes showing previous measurements
    e1 = Rectangle((0.04, 150), 1.96, 300, color="cyan", fill=True,
                   label="Ye et al. (2016)")
    e2 = Rectangle(((1*0.0409*10**3), 1), (5.5*0.0409*10**3), 29,
                   color="green", fill=True, label="Romer et al. (2018)")
    e3 = Rectangle(((40*0.0409*10**3), 0), (210*0.0409*10**3), 10,
                   color="coral", fill=True, label="Shi et al. (2021)")

    # Dataset of aircraft SLR measurements
    folder = get_local_folder('ARNA_data') + '/Filters/'
    filename = 'FAAM_data_for_missing_HONO_source_plot.csv'
    df1 = pd.read_csv(folder+filename, index_col=0)

    # Dataset of ground HONO and pNO3 data combined with modelled J-rates
    filename = 'CV_data_for_missing_HONO_source_plot.csv'
    df2 = pd.read_csv(folder+filename, index_col=0,
                             parse_dates=True, dayfirst=True)
    df2 = df2.loc[(df2.index.hour > 10) & (df2.index.hour < 16)]

    # Ye et al. (2017) individual measurements
    filename = 'Ye2017_data.csv'
    df3 = pd.read_csv(folder+filename)
    df3["Enhancement"] = df3["JpNO3"]/(7*10**-7)

    # Plotting settings
    if isinstance(fontsize, type(None)):
#        fontsize = 20
        fontsize = 8
    if isinstance(labelsize, type(None)):
#        labelsize = 16
        labelsize = fontsize/5*4
    if isinstance(ms, type(None)):
#        ms = 10
        ms = 6

    # Creating legend elements
    from matplotlib.patches import Patch
    legend_elements = [Line2D([0], [0], marker='o', color='white',
                              label='Sea-salt',
                              markerfacecolor="blue", ms=ms),
                       Line2D([0], [0], marker='o', color='white',
                              label='Dust',
                              markerfacecolor='red', ms=ms),
                       Line2D([0], [0], marker='o', color='white',
                              label=r'Biomass burning',
                              markerfacecolor='black', ms=ms),
                       Line2D([0], [0], marker='o', color='white',
                              label=r'Free troposphere',
                              markerfacecolor='lime', ms=ms),
                       Line2D([0], [0], marker='o', color='white',
                              label='Mixed dust/biomass burning',
                              markerfacecolor='orange', ms=ms),
                       Line2D([0], [0], marker='o', color='white',
                              label='CVAO 2019',
                              markerfacecolor='purple', ms=ms),
                       Line2D([0], [0], marker='o', color='white',
                              label='Ye et al (2017)',
                              markerfacecolor='lightgrey', ms=ms),
                       Patch(facecolor='cyan', edgecolor='cyan',
                             label="Ye et al (2016)"),
                       Patch(facecolor='green', edgecolor='green',
                             label="Romer et al (2018)"),
                       Patch(facecolor='coral', edgecolor='coral',
                             label="Shi et al (2021)"),
                       Line2D([0], [0], color='black', linestyle="dashed",
                              lw=4, label='Ye et al (2017) fit'),
                       Line2D([0], [0], color='black', lw=4,
                              label='Langmuir fit (this study)')]

    # Creating colour scheme for scatterplots - should probably consider changing the colours...
    c1_new = np.where(df1["Categories"] == 1, "blue", "orange")
    c2_new = np.where(df1["Categories"] == 2, "lime", c1_new)
    c3_new = np.where(df1["Categories"] == 4, "black", c2_new)
    c4_new = np.where(df1["Categories"] == 3, "red", c3_new)

    # Creating the Langmuir fit and the Ye et al. fit.
    x = np.linspace(0.1, 10000, 100000)
    y = (100*18.8*0.9)/(1+0.9*x)
    y2 = ((6.1*10**-4*np.log(1+4.4*10**-1*x))/x)-3.5*10**-5
    y3 = y2/(7*10**-7)

    # Plotting all the data
    f = plt.figure(figsize=figsize, dpi=dpi)
    if not isinstance(aspect, type(None)):
        AC.adjustFigAspect(f, aspect=aspect)
    ax1 = plt.subplot2grid((1, 2), (0, 0), colspan=1)
    ax2 = plt.subplot2grid((1, 2), (0, 1), colspan=1)

    ax1.errorbar(x=3600*df1["J_HNO3"]*df1["NO3_ppt"],
                 y=3600*df1["Missing_HONO_source"],
                 xerr=(df1["Photolysis_nitrate_uncertainty"]
                       * 3600*df1["J_HNO3"]*df1["NO3_ppt"]),
                 yerr=df1["Missing_HONO_source_uncertainty"]*3600,
                 ecolor=c4_new, fmt="none", alpha=0.2)
    ax1.scatter(x=3600*df1["J_HNO3"]*df1["NO3_ppt"],
                y=3600*df1["Missing_HONO_source"],
                c=c4_new, zorder=100)

    __x = df2["JVL_016"]*3600*((df2["Nitrate.ug_m3"])/(62*4.09*10**-5))
    ax1.scatter(x=__x, y=df2["Missing_HONO_source_per_hour"], color="purple")
    ax1.tick_params(axis="both", labelsize=labelsize)
    ax1.set_xlim(0, 5.5)
    ax1.set_ylim(0, 180)
    ax1.text(0.05, 0.95, "(A)", fontsize=fontsize, transform=ax1.transAxes)
    ax1.set_ylabel("Missing HONO source (ppt h$^{-1}$)", fontsize=fontsize)
    label = r"$\it{j}_\mathrm{HNO_3}$ $\times$ "
    label += r"[pNO$_3^-$]$_\mathrm{bulk}$ (ppt h$^{-1}$)"
    ax1.set_xlabel(label, fontsize=fontsize)

    ax2.add_artist(e1)
    ax2.add_artist(e2)
    ax2.add_artist(e3)
    ax2.errorbar(df1["Total_NO3_nmoles_m-3"],
                 df1["Enhancement"],
                 xerr=df1["Total_NO3_nmoles_m-3_uncertainty"],
                 yerr=(df1["Enhancement_uncertainty_%"]
                       * df1["Enhancement"]),
                 ecolor=c4_new, fmt="none", alpha=0.2)
    ax2.scatter(df1["Total_NO3_nmoles_m-3"],
                df1["Enhancement"],
                c=c4_new, zorder=100)
    ax2.scatter(x=(df2["Nitrate.ug_m3"]/(0.062)),
                y=df2["Missing_HONO_source"] /
                (df2["JVL_016"] *
                 (df2["Nitrate.ug_m3"]/(62*4.09*10**-5))),
                color="purple", zorder=99)
    ax2.scatter(x=df3["nmoles_m-3"],
                y=df3["Enhancement"],
                color="lightgrey", zorder=99)
    ax2.plot(x, y, color="black")
    ax2.plot(x, y3, color="black", linestyle="dashed")
    ax2.set_xlim(0.1, 10000)
#    ax2.set_ylim(0, 20000) # Manually set the y-axis extents?
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.tick_params(axis="both", labelsize=labelsize)
    label = "[pNO$_3^-$]$_\mathrm{bulk}$ (10$^{-9}$ moles m$^{-3}$)"
    ax2.set_xlabel(label, fontsize=fontsize)
    ax2.set_ylabel("Enhancement Factor ($\it{f}$)", fontsize=fontsize)
    ax2.text(0.05, 0.95, "(B)", fontsize=fontsize, transform=ax2.transAxes)

    label = r"$\it{f}$ = $\frac{\it{a Q^0 K_{ads}}}"
    label += r"{1 + \it{K_{ads}} [\mathrm{pNO_3^-}]_\mathrm{bulk}}$"
    ax2.text(20, 2000, label, fontsize=fontsize)

    ax2.legend(handles=legend_elements,
             fontsize=fontsize/5*4,
             loc=(1.05, 0.05),
             frameon=False)
    # Save figure
    AC.save_plot('ARNA_Andersen_figure_02', dpi=dpi, tight=tight)


def plt_seasonal_comparoisons_of_nitrate():
    """
    Make a plot of seasonal nitrate at CVAO
    """
    #
    RunSet = 'ACID'
    res = '4x5'
    CoreRunsOnly = False
    RunDict = get_dict_of_GEOSChem_model_output(res=res, RunSet=RunSet,
                                                CoreRunsOnly=CoreRunsOnly,
                                                folder4netCDF=True)
    # Choose runs to use
    d = {}
    d['BC-BASE'] =  '/users/ts551/scratch/GC/rundirs/P_ARNA//geosfp_4x5_standard.v12.9.0.BASE.2019.2020.ARNA.BCs.repeat//OutputDir/'
    d['Acid-4x5-J00'] = '/users/ts551/scratch/GC/rundirs/P_ARNA//geosfp_4x5_aciduptake.v12.9.0.BASE.2019.2020.ARNA.DustUptake.JNIT.Isotherm.BCs.repeat.ON.II.diags.v2.J00//OutputDir/'
    dates2use = [datetime.datetime(2019, 1+i, 1) for i in range(12)]
    dsD = {}
    for key in d.keys():
        dsD[key] = AC.GetSpeciesConcDataset(wd=d[key],
                                              dates2use=dates2use)

    # Get model
    from Prj_NITs_analysis import get_CVAO_NITs_data
#    df_obs = get_CVAO_NITs_data()
    # Update dates
    dates = AC.dt64_2_dt( df_obs.index.values )
    df_obs.index = [ update_year(i, year=2019) for i in df_obs.index ]

    # Create a 'NIT-all'
    prefix = 'SpeciesConc_'
    for key in d.keys():
        ds = dsD[key]
        ds = AC.AddChemicalFamily2Dataset(ds, fam='NIT-all', prefix=prefix)
        dsD[key] = ds

    # Select for CVAO
    from funcs4obs import gaw_2_loc
    site = 'CVO'
    lat, lon, alt, TZ = gaw_2_loc(site)
    for key in d.keys():
        ds = dsD[key]
        ds = ds.sel(lat=lat, lon=lon, method='nearest')
        ds = ds.sel(lev=ds.lev.values[0] )
        dsD[key] = ds
    # remove the redundent coordinates
    for key in d.keys():
        ds = dsD[key]
        del ds['lat']
        del ds['lon']
        del ds['lev']
        dsD[key] = ds

    # plot up
    import seaborn as sns
    sns.set(color_codes=True)
    var2plot = 'SpeciesConc_NIT-all'
    fig, ax = plt.subplots()
    colors = AC.get_CB_color_cycle()
    for nKey, key in enumerate( d.keys() ):
        ds = dsD[key]
        AC.BASIC_seasonal_plot(dates=ds.time.values,
                               color=colors[nKey],
                               data=ds[var2plot].values*1E12,
                               ax=ax,
                               label=key)

    # Add observations
    AC.BASIC_seasonal_plot(dates=df_obs.index, color='k',
                           data=df_obs['Nitrate, pptv'].values, ax=ax,
                           label='Obs.')

    plt.ylabel('All nitrate (inc. on dust)')
    plt.legend()
    plt.title('Seasonal cycle of ntirate at CVAO (pptv)')
    AC.save_plot(title='NITs_comparison_CVAO_model_obs', dpi=720)
    plt.close('all')

    # rename obs
    species = [i.split(prefix)[-1] for i in ds.data_vars]
    for key in d.keys():
        ds = dsD[key]
        ds = ds.rename(name_dict=dict(zip(ds.data_vars, species)))
        dsD[key] = ds

    # plot up a stacked plot of nitrate by season.
    key = 'Acid-4x5-J00'
    ds = dsD[key]
    del ds['ilev']
#    vars2plot = [i for i in ds.data_vars if 'NIT' in i ]
    vars2plot = ['NITD4', 'NITD3', 'NITD2', 'NITD1', 'NITs', 'NIT',]
    df = ds[vars2plot+ ['NIT-all']].to_dataframe()
    df = df * 1E12 # Conevrt to pptv
#    for var in

    plt.close('all')
    fig, ax = plt.subplots()
    df[vars2plot].plot.area()

    # Add monthly mean for obs.
    df_obs = df_obs.resample('1M').mean()
    df_obs['Nitrate, pptv'].plot(label='obs', color='k', ls='--', zorder=99)
    plt.autoscale(enable=True, axis='y')
    plt.legend()

    AC.save_plot(title='NITs_comparison_CVAO_by_species', dpi=720)
    plt.close('all')


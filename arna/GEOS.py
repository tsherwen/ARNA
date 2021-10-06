"""
Functions to interact with GEOS-CF/GEOS-5/GMAO output during the ARNA campaign
"""
import os
import sys
import glob
import gc
import requests
import re
import wget
import numpy as np
import pandas as pd
import xarray as xr
import xesmf as xe
import AC_tools as AC
import datetime as datetime
import time
from time import gmtime, strftime
from bs4 import BeautifulSoup
from PIL import Image, ImageDraw
import PIL

# Import from elsewhere in ARNA module
from . core import *
from . utils import *
from . observations import *


def get_data_for_fcasts4datetimes(dts=None):
    """
    Download all forecast data in bulk for given datetimes (dts)
    """
    # Use the last 18 days, unless specific dates are provided
    if isinstance(dts, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
        # Use the last 18 days
        dts = [AC.add_days(dt, i*-1) for i in range(1, 18)]
    # Loop by date and download data if available and not downloaded.
    for dt in dts:
        dstr = dt.strftime('%Y/%m/%d %H:%M')
        pstr = 'Checking historic GEOS5-CF forecasts are downloaded for {}'
        print(pstr.format(dstr))
        # Check whether files are available
        # TODO - Add check for GEOS5/CF for a specific date
        # Check if the files have been downloaded
        n_failed_doys = check4failed_downloads(dt=dt, retry=False,
                                               inc_GEOS5_Nv_collection=True)
        print(n_failed_doys)
        # If there files are not present, then download in series
        # TODO: Confirm why this call repeated - does this matter?
        n_failed_doys = check4failed_downloads(dt=dt, retry=True,
                                               inc_GEOS5_Nv_collection=True)
        print(n_failed_doys)


def get_expanded_variable_GMAO_files(verbose=True, debug=False):
    """
    Get the expanded variable NetCDF files from GMAO
    """
    # What is the root URL for the data?
    URL = 'https://gmao.gsfc.nasa.gov/gmaoftp/geoscf/ARNA/'
    # Access the html via requests
    r = requests.get(URL)
    html = r.text
    soup = BeautifulSoup(html, "html.parser")
    links = soup.findAll(href=re.compile("\.nc4$"))
    # Where to save the files?
    folder = get_local_folder('NASA_data')
    folder += 'GEOS_CF/ARNA/assim/expanded_variable_files/'
    AC.mk_folder(folder, verbose=verbose)
    # Check that the URL was correct
    assert len(links) > 0, 'WARNING: No links found @ URL ({})'.format(URL)
    # Get the one image
    for link in links:
        filename = link['href']
        # Check it file present...
        if os.path.isfile(folder+filename):
            pstr = 'WARNING: Not downloading GMAO file as it exists! ({})'
            print(pstr.format(folder+filename))
        # If not, download
        else:
            print(URL+filename)
            wget.download(URL+filename, folder+filename)


def mk_test_plots_from_GMAO_output(dt=None, load_all_dates=True):
    """
    Make a set of test plots for the expanded NASA GMAO plots
    """
    # Which day(s) to use for testing output?
    if isinstance(dt, type(None)) and (not load_all_dates):
        # Use the last 18 days, unless specific dates are provided
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
        # Setup the glob string for the files
        dstr = dt.strftime('%Y/%m/%d')
        gfstr = 'GEOS-CF.v01.rpl.ACSIS-ARNA.{}*.nc4'.format(dstr)

    elif load_all_dates:
        gfstr = 'GEOS-CF.v01.rpl.ACSIS-ARNA.*.nc4'
    else:
        # Setup the glob string for the files
        dstr = dt.strftime('%Y/%m/%d')
        gfstr = 'GEOS-CF.v01.rpl.ACSIS-ARNA.{}*.nc4'.format(dstr)
    # Where are the files?
    folder = get_local_folder('NASA_data')
    folder += 'GEOS_CF/ARNA/assim/expanded_variable_files/'
    # Load the NetCDF files
    files2use = glob.glob(folder+gfstr)
    ds = xr.open_mfdataset(files2use)
    # Update the units for lev
    HPa_l = AC.get_GEOSCF_vertical_levels(native_levels=True)
    attrs = ds.lev.attrs.copy()
    attrs['units'] = 'hPa'
    attrs['standard_name'] = 'Pressure'
    attrs['long_name'] = 'Pressure'
    ds.lev.attrs = attrs
    ds.lev.values = [HPa_l[int(i)] for i in ds.lev.values - 1]

    # Update units
    ds = convert_GEOSCF_units(ds=ds, debug=True)
    # Now do a handful of quick plots
#    extra_str = 'GMAO_EXPANDED_TEST'
#    vars2plot = [i for i in ds.data_vars]
    vars2plot = ['O3', 'CO', 'Cl', 'U', 'V', 'T', 'IO', 'BrO']
    title_date = 'avg 24/25th jan'
    folder = './'
#    for var2plot in vars2plot[:3]:
    for var2plot in vars2plot:
        #        for lev2use in ds.lev.values:
        #        for lev2use in [72, 51]:
        for lev2use in [525., 985.]:
            #            for lev2use in ds.lev.values[:2]: # unhash if testing
            print(var2plot, lev2use, title_date)
#            ds_tmp = ds[[var2plot]].sel(lev=lev2use)
            bool1 = ds.lev.values == lev2use
            ds_tmp = ds[[var2plot]].isel(lev=bool1)
            ds_tmp.squeeze()
            #
            ds_tmp = ds_tmp.mean(dim='time')
#            if not isinstance(extr_title_str, type(None)):
#                title += '\n '+ extr_title_str
            # Get the LateX for of the species name
            try:
                LaTeX_spec = AC.latex_spec_name(var2plot)
            except KeyError:
                pstr = 'WARNING: not converted {} to LaTeX form'
                print(pstr.format(var2plot))
                LaTeX_spec = var2plot
            # Force use of standard name as long name
            attrs = ds[var2plot].attrs
            attrs['long_name'] = LaTeX_spec
            ds_tmp[var2plot].attrs = attrs
            # Setup a string for the title
            title = '[{}] @ level {:0>2} on {}'
            title = title.format(LaTeX_spec, lev2use, title_date)
            # Set extra string for filename
            extra_str = 'ARNA_lev_{:0>2}_{}'.format(lev2use,
                                                    'GMAO_EXPANDED_TEST')
            # Now plot
            quick_map_plt_CV_1layer(ds_tmp, var2plot=var2plot,
                                    use_local_CVAO_area=True,
                                    extra_str=extra_str,
                                    extend='both',
                                    title=title,
                                    folder=folder,
                                    save_plot=True)

            # Do some garbage collecting
            gc.collect()


def regrid_GEOS5_files_in_folder(folder=None, dt=None, doys2use=None,
                                 year=None,
                                 collection='inst3_3d_aer_Np',
                                 remake_files=False, debug=False):
    """
    Regrid all GEOS-5 NetCDF files in folder to GEOS-CF format
    """
    #
    # Use the latest downloaded GEOS-5 folder if one is not set.
    if isinstance(folder, type(None)):
        # Use yesterday as the date if this is not provided
        if isinstance(dt, type(None)):
            # Just use yesterday for Now
            Tnow = AC.time2datetime([gmtime()])[0]
            # Get the 5-day forecast at noon...
            dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
            # Use yesterday
            dt = AC.add_days(dt, -1)
        # Get the GEOS5 folder for a given date
        folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                         collection=collection)

    # Get current year if a year not specified
    if isinstance(year, type(None)):
        year = AC.time2datetime([gmtime()])[0].year
    # Which files are available
    files_in_folder = glob.glob(folder+'ARNA*_{}_*.nc'.format(year))
    # Make sure already regridded files are not used
    files_in_folder = [i for i in files_in_folder if 'REGRIDDED' not in i]
    # Which doys to use
    if isinstance(doys2use, type(None)):
        doys2use = [i.split('/')[-1] for i in files_in_folder]
        doys2use = [i.split(str(year))[-1][1:4] for i in doys2use]
        doys2use = list(sorted((set(doys2use))))
    # Save name for the regridded file?
    writestr = 'ARNA_GEO5_Combined_Cape_Verde_{}_{:0>3}_REGRIDDED_ALL_lvls.nc'
    # Loop by doy and process
    for doy in doys2use:
        # Select files for a given doy
        doy_str = '_{:0>3}_'.format(doy)
        files2use = [i for i in files_in_folder if doy_str in i]
        if debug:
            print(len(files2use), [i.split('/')[-1] for i in files2use], '')
        # Set the name of the file to save
        filename2save = writestr.format(year, doy)
        if os.path.isfile(folder + filename2save) and not remake_files:
            if debug:
                pstr = 'WARNING: not regridding file as file already present'
                print(pstr)
        else:
            # Open all files as a single xarray dataset
            ds = xr.open_mfdataset(files2use)
            if debug:
                print([i for i in ds.time.values])
            # Convert files to a single file with variable for dust
            ds = convert_GEOS5_dust2_single_var_in_ug_m3(ds=ds,
                                                         rtn_just_new_var=True)
            # Regridd file to GEOS-CF template file
            # NOTE: this is hardcoded for now
            tplate_folder = get_local_folder('NASA_data')
            tplate_folder += 'GEOS_CF/ARNA/fcast/2019_11_06/data/'
            tplate_fname = 'ARNA_GEOSCF_chm_inst_1hr_g1440x721_p23_Cape_Verde_'
            tplate_fname += '2019_315_pm25du_rh35_gcc_lvls_1000_850_500.nc'
            # Regrid GEOS5 to GEOS-CF
            # TODO: add option to not vertically regrid
            ds = regrid_GEOS52GEOSCF_coordinates(ds, collection=collection,
                                                 tplate_folder=tplate_folder,
                                                 tplate_fname=tplate_fname)
            # Save as a new file
            ds.to_netcdf(folder+filename2save)


def regrid_GEOS52GEOSCF_coordinates(ds=None, collection='inst3_3d_aer_Np',
                                    tplate_folder=None, tplate_fname=None,
                                    #                                    regrid_vertically=False,
                                    lat_var='lat', lon_var='lon',
                                    vars2regrid=None):
    """
    Regrid GEOS-5 output to be on the same grid as GEOS-CF
    """
    # Which variables to regrid?
    if isinstance(vars2regrid, type(None)):
        vars2regrid = list(ds.data_vars)

    # - Regrid horizontally (lat, lon)
    # Get lats and lons to regrid to
    dsT = xr.open_dataset(tplate_folder+tplate_fname)
    lat = dsT[lat_var]
    lon = dsT[lon_var]
    # Create a new dataset to template to from this
    ds_out = xr.Dataset({
        # 'time': ( ['time'], ds['time'] ),
        'lat': (['lat'], lat),
        'lon': (['lon'], lon),
    })
    # Create a regridder (to be reused )
    regridder = xe.Regridder(ds, ds_out, 'bilinear', reuse_weights=True)
    ds_l = []
    for var2use in vars2regrid:
        # Create a dataset to re-grid into
        ds_out = xr.Dataset({
            #            'time': ( ['time'], ds['time'] ),
            'lat': (['lat'], lat),
            'lon': (['lon'], lon),
        })
        # Get a DataArray
        dr = ds[var2use]
        # Build regridder
        dr_out = regridder(dr)
        # Important note: Extra dimensions must be on the left, i.e. (time, lev, lat, lon) is correct but (lat, lon, time, lev) would not work. Most data sets should have (lat, lon) on the right (being the fastest changing dimension in the memory). If not, use DataArray.transpose or numpy.transpose to preprocess the data.
        # Exactly the same as input?
        xr.testing.assert_identical(dr_out['time'], ds['time'])
        # Save variable
        ds_l += [dr_out]

    # Combine variables
    dsN = xr.Dataset()
    for n, var2use in enumerate(vars2regrid):
        dsN[var2use] = ds_l[n]
        # use attributes of input dataset
        dsN[var2use].attrs = ds[var2use].attrs
    # Clean up
#    regridder.clean_weight_file()

    # - Regrid vertically (alt)
#    if regrid_vertically:
    if collection == 'inst3_3d_aer_Nv':
        # use bulk values for altitude
        lon, lat, alt = AC.get_latlonalt4res(res='4x5',
                                             full_vertical_grid=True)
        # NOTE: GEOS-5 output has its z axis flipped.
        dsN.lev.values = alt[::-1]
        # Which
        HPa_l = AC.get_GEOSCF_vertical_levels(native_levels=True)
        hPa_as_km = [i for i in AC.hPa_to_Km(HPa_l)]
        # Now interpolate
        dsN = dsN.interp(lev=hPa_as_km)
        # Put the units back as millibar
        dsN.lev.values = HPa_l
    # Return the dataset
    return dsN


def convert_GEOS5_dust2_single_var_in_ug_m3(ds=None, rtn_just_new_var=True):
    """
    Convert GEOS5 dust (in v/v) into ug/m3
    """
    # Name for combined data
    NewVar = 'Dust'
    # If more than 3 dust tracers, then sum up the dust tracers
    if len([i for i in ds.data_vars if 'du' in i]) > 3:
        dust_tracers = ['du{:0>3}'.format(i) for i in range(1, 6)]
        ds[NewVar] = ds[dust_tracers[0]].copy()
        for tra in dust_tracers[1:]:
            ds[NewVar].values = ds[NewVar].values + ds[tra].values
    else:
        ds = ds.rename({'du': NewVar})
    # Now convert to ug/3m
    print(ds.time)
    data = ds[NewVar].values
    ds[NewVar].values = AC.convert_spec_v_v_2_ugm3(data=data, spec='DST1')
    # Update the units
    attrs = ds[NewVar].attrs
    attrs['units'] = '$\mu$g m$^{-3}$'
    attrs['long_name'] = 'Dust ($\Sigma$ bins 1-5)'
    ds[NewVar].attrs = attrs
    if rtn_just_new_var:
        return ds[[NewVar]]
    else:
        return ds


def get_most_recent_GEOSCF_data(dt=None,
                                limit_lvls=True, limit_lons=False,
                                filestr='ARNA*{}*lvls_10*.nc',
                                collection='chm_inst_1hr_g1440x721_p23',
                                ):
    """
    Get most recent GEOS-CF output (downloaded and processed)
    """
    if isinstance(dt, type(None)):
        # Just use yesterday for Now
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Get the date/location string for data
    filename = filestr.format(dt.year)
    folder2use = get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                         collection=collection)
    glob_str = '{}/{}'.format(folder2use, filename)
    files2use = glob.glob(glob_str)
    files2use = list(sorted(files2use))
    asstr = 'WARNING no files found matching {}'
    assert len(files2use), asstr.format(glob_str)
    # Open all the files as a single dataset
    ds = xr.open_mfdataset(files2use)
    print(ds.time)
    # Update the variable names to be in GEOS-Chem format
    ds = update_GEOSCF_names2GEOSChem_names(ds=ds)
    # Update the units to GEOS-Chem/interpretable ones
    ds = convert_GEOSCF_units(ds=ds)
    # Add extra analysis species
    ds = add_extra_derived_species(ds=ds)
    return ds


def get_GEOS_data_folder4dt(dt=None, product='GEOS_CF', host=None,
                            mode='fcast', collection='inst3_3d_aer_Np',
                            inc_collection=True):
    """
    Get the data folder location for a given GEOS product and datetime
    """
    # Where is the data for the product?
    NASA_data = get_local_folder('NASA_data', host=host)
    folder = '{}/{}/ARNA/'.format(NASA_data, product)
    # And for the datetime
    date_str = '{}_{:0>2}_{:0>2}_{:0>2}z/'
    date_str = date_str.format(dt.year, dt.month, dt.day, dt.hour)
    # Return the string including collection and production mode
    if inc_collection:
        return '{}/{}/{}/data.{}/'.format(folder, mode, date_str, collection)
    else:
        return '{}/{}/{}/'.format(folder, mode, date_str)


def get_latest_GEOSCF_fcast_data(dt=None, just_check_yesterday=True,
                                 limit_lvls=False,
                                 vars2use=None, doys2use=None,
                                 retry_download_if_data_not_present=True,
                                 limit_lons=False, rm_existing_file=False,
                                 limit_lats=False, debug=True):
    """
    Get the latest GEOS-CF forecast data (if new data available)
    """
    if debug:
        pstr = 'limit_lons={}, limit_lvls={}, limit_lats={} - in {}'
        func = 'get_latest_GEOSCF_fcast_data'
        print(pstr.format(limit_lons, limit_lvls, limit_lats, func))
    # Set date to use if not provided
    if isinstance(dt, type(None)):
        if just_check_yesterday:
            # Just use yesterday for Now
            TNow = AC.time2datetime([gmtime()])[0]
            dt = AC.add_days(TNow, -1)
            dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
        else:
            # What is the last data there is data locally for?
            last_local_data_date = what_is_latest_data_locally(debug=debug)
            assert_err = 'Stopping as no data for last 5 days'
            assert not isinstance(last_local_data_date, type(None)), assert_err
            # When was the last forecast?
            last_fcast_start = when_should_GEOSCF_have_last_run_from()
            # Use this as the datetime to download
            dt = last_fcast_start
            ass_str = 'The most recent data is already locally available'
            assert last_fcast_start > last_local_data_date, ass_str
            # Download data if available from nasa
#            lastest_GEOSCF_avail = check_if_latest_GEOSCF_is_available(last_fcast_start)
    else:
        # Using the provided datetime (dt)
        pass
    # Setup a text string with for dt
    dstr = dt.strftime('%Y/%m/%d %H:%M')

    # -
    # Is the data available
    dt_is_lastest_run = is_dt_latest_GEOSCF(dt=dt)
    # If it is, just download it
    if dt_is_lastest_run:
        # - Retrieve the latest GEOS-CF data - sliced by alt, lat, or lon
        download_GEOSCF_fcast_data4date(dt=dt,
                                        rm_existing_file=rm_existing_file,
                                        vars2use=vars2use,
                                        doys2use=doys2use,
                                        limit_lvls=limit_lvls,
                                        limit_lons=limit_lons,
                                        limit_lats=limit_lats)
    elif retry_download_if_data_not_present:
        # Minutes to wait
        times2wait = 30
        min2wait = 4
        pstr = 'WARNING: data not availible ({}), so trying every {} min.({}x)'
        print(pstr.format(dstr, min2wait, times2wait))
        for time2wait in range(times2wait):
            pstr = 'WARNING: Wait for {:>2} min @ {}, then retrying download'
            TNow = AC.time2datetime([gmtime()])[0]
            print(pstr.format(min2wait, TNow.strftime('%Y/%m/%d %H:%M')))
            time.sleep(min2wait*60)
            pstr = 'Finished waiting @ {}'
            print(pstr.format(TNow.strftime('%Y/%m/%d %H:%M')))
            # Check if the data is availible
            dt_is_lastest_run = is_dt_latest_GEOSCF(dt)
            if dt_is_lastest_run:
                print('Data now available! Attempting download.')
                # - Retrieve the latest GEOS-CF data - sliced by alt, lat, or lon
                download_GEOSCF_fcast_data4date(dt=dt,
                                                rm_existing_file=rm_existing_file,
                                                vars2use=vars2use,
                                                doys2use=doys2use,
                                                limit_lvls=limit_lvls,
                                                limit_lons=limit_lons,
                                                limit_lats=limit_lats)
                print('WARNING: Stopping. Completed attempt at data download.')
                sys.exit()
            else:
                print('Data still not available...')
    else:
        print('ERROR: Data not available for date so cannot attempt download')


def get_latest_GEOSCF_fcast_data_alt_slice(dt=None, debug=True):
    """
    Wrapper to get alt slices from the latest GEOS-CF forecast data
    """
    # Hardcode options here for now.
    limit_lvls = True
    limit_lons = False
    # Call the main fcast data retrieval function
    get_latest_GEOSCF_fcast_data(dt=dt,
                                 limit_lvls=limit_lvls,
                                 limit_lons=limit_lons,
                                 #                                 limit_lats=limit_lats
                                 )


def check4GEOSCF_failed_downloads(dt=None, folder=None, n_doys=6,
                                  var2use='noy',
                                  limit_lons=False, limit_lvls=False,
                                  limit_lats=False,
                                  debug=False):
    """
    Check that all the GEOSCF files have downloaded correctly
    """
    import os
    import glob
    # Print
    if debug:
        pstr = 'limit_lons={}, limit_lvls={}, limit_lats={} - in {}'
        func = 'check4GEOSCF_failed_downloads'
        print(pstr.format(limit_lons, limit_lvls, limit_lats, func))
    # Sliced by longitude
    if (limit_lvls == True) and (limit_lons == False) and (limit_lats == False):
        ext_str = '*lvls_1000_900_800_700_600_500*'
    # Sliced by altitude
    if (limit_lvls == False) and (limit_lons == True) and (limit_lats == False):
        ext_str = '*_lons_*'
    # Sliced by latitude
    if (limit_lvls == False) and (limit_lons == False) and (limit_lats == True):
        ext_str = '*_lats_*'
    # Also only consider one variable at a fime
    file_str = '*{}*{}*'.format(var2use, ext_str)
    # - Work out what the expected files for a date are
    expected_dates = [AC.add_days(dt, i) for i in range(n_doys)]
    dfT = pd.DataFrame(index=expected_dates)
    expected_doys = list(dfT.index.dayofyear.values)
    # -  Now test for files
    # Get file(s)
    files = glob.glob(folder+file_str)
    files = list(sorted(files))
    # Check if any files were found?
    if len(files) > 0:
        # Get size(s) of files  and the doys they are for
        doys = [i.split(str(dt.year))[-1][1:4] for i in files]
        sizes = [os.path.getsize(i) for i in files]
        # Compile into a dataframe
        df = pd.DataFrame({'doy': doys, 'sizes': sizes})
        if debug:
            print(df)
        # Get details on dimensions
        dims = []
        for file in files:
            try:
                dims2add = [dict(xr.open_dataset(file).dims)]
            except OSError:
                dims2add = [[]]
            # Add the dimensions to the save list
            if (len(dims2add[0]) != 4):
                dims2add = [{'lat': 0, 'lev': 0, 'lon': 0, 'time': 0}]
            dims += dims2add
        # Add the dimensions
        for n_file, file in enumerate(files):
            dims4file = dims[n_file]
            for key in dims4file.keys():
                df.loc[n_file, key] = dims4file[key]
        # - Add checks that all the lat, lon and, lev dims are the same length?
        if debug:
            print(df)
        coords2test = ['lon', 'lat', 'lev']
        for coord in coords2test:
            mode_ = df[coord].mode().values[0]
            # Add a column with a test for the coordinates
            var2use = 'eq_{}_mode'.format(coord)
            df[var2use] = df[coord] == mode_
        # - Do tests on doys and save those that fail to a list
        failed_doys = []
        # - Check if all coords are modal
        # Check that all coords are the same
        cols2use = [i for i in df.columns if 'mode' in i]
        if debug:
            print(df[cols2use].values.all())
            print('')
            print(df)
        if df[cols2use].values.all():
            pass
        else:
            for doy in df['doy'].values:
                tmp = df.loc[df['doy'] == doy, :]
                if debug:
                    print(tmp)
                if tmp[cols2use].values.all():
                    pass
                else:
                    failed_doys += [int(tmp['doy'].values.astype(int))]

        # - Check if the expected number of timestamps are in each time
        expected_times = [12.0, 24.0, 24.0, 24.0, 24.0, 13.0][:df.shape[0]]
        var2use = 'expected_time'
        df[var2use] = df['time'].values == expected_times
        # Save a list of the doys for which the retrieval failed for
        tmp_df = df.loc[df['expected_time'] == False, :]['doy']
        failed_doys += list(tmp_df.values.astype(int))
        # - Check if all the expected doys are present
        retrieved_doys = [int(i) for i in df['doy'].values]
        failed_doys += [i for i in expected_doys if (i not in retrieved_doys)]
        # Print a warning if there are doys included that are not expected
        unexpected_doys = [i for i in retrieved_doys if i not in expected_doys]
        if len(unexpected_doys) > 0:
            pstr = 'WARNING: Unexpected doys have been downloaded -'
            print(pstr, unexpected_doys)
        # - Return a list of the failed doys
        if debug:
            print(df, failed_doys)
        return list(sorted(set(failed_doys)))
    else:
        pstr = "WARNING: No files for '{}'  - using expected doys (folder:{})"
        print(pstr.format(var2use, folder))
        return expected_doys


def check4GEOS5_failed_downloads(dt=None, folder=None, n_doys=6, var2use='du',
                                 debug=False):
    """
    Check that all the GEOS5 files have downloaded correctly
    """
    import os
    import glob
    # Get file(s), size(s), and the doys they are for
    # Also only consider one variable at a fime
    files = glob.glob(folder+'*_{}_*.nc*'.format(var2use))
    # Make sure only none regridded files are checked
    files = [i for i in files if 'REGRIDDED' not in i]
    files = list(sorted(files))
    # - Work out what the expected files for a date are
    expected_dates = [AC.add_days(dt, i) for i in range(n_doys)]
    dfT = pd.DataFrame(index=expected_dates)
    expected_doys = list(dfT.index.dayofyear.values)
    # Check if any files were found?
    if len(files) > 0:
        # Get size(s) of files  and the doys they are for
        doys = [i.split(str(dt.year))[-1][1:4] for i in files]
        sizes = [os.path.getsize(i) for i in files]
        # Compile into a dataframe
        df = pd.DataFrame({'doy': doys, 'sizes': sizes})
        # Get details on dimensions
        dims = []
        for file in files:
            try:
                dims2add = [dict(xr.open_dataset(file).dims)]
            except OSError:
                dims2add = [[]]
            # Add the dimensions to the save list
            if debug:
                print(len(dims2add), dims2add)
            if (len(dims2add[0]) != 4):
                dims2add = [{'lat': 0, 'lev': 0, 'lon': 0, 'time': 0}]
            dims += dims2add
        # Add the dimensions
        for n_file, file in enumerate(files):
            dims4file = dims[n_file]
            for key in dims4file.keys():
                df.loc[n_file, key] = dims4file[key]
        # Do tests on other variables
    #    assert len(set(df['lat'])) == 1, 'WARNING: files do not have the same # of lats!'
    #    assert len(set(df['lon'])) == 1, 'WARNING: files do not have the same # of lats!'
    #    assert len(set(df['lev'])) == 1, 'WARNING: files do not have the same # of lats!'
        # - add checks that all the lat, lon and, lev dims are the same length?
        coords2test = ['lon', 'lat', 'lev']
        for coord in coords2test:
            mode_ = df[coord].mode().values[0]
            # Add a column with a test for the coordinates
            var2use = 'eq_{}_mode'.format(coord)
            df[var2use] = df[coord] == mode_
        # - Do tests on doys and save those that fail to a list
        failed_doys = []
        # - Check if all coords are modal
        # Check that all coords are the same
        cols2use = [i for i in df.columns if 'mode' in i]
        if df[cols2use].values.all():
            pass
        else:
            for doy in df['doy'].values:
                tmp = df.loc[df['doy'] == doy, :]
                if tmp[cols2use].values.all():
                    pass
                else:
                    failed_doys += [int(tmp['doy'].values)]

        # - Check if the expected number of timestamps are in each time
        expected_times = [4., 8., 8., 8., 8., 5.][:df.shape[0]]
        var2use = 'expected_time'
        df[var2use] = df['time'].values == expected_times
        # Save a list of the doys for which the retrieval failed for
        df_tmp = df.loc[df['expected_time'] == False, :]['doy']
        failed_doys += list(df_tmp.values.astype(int))
        # - Check if all the expected doys are present
        expected_doys = list(dfT.index.dayofyear.values)
        retrieved_doys = [int(i) for i in df['doy'].values]
        failed_doys += [i for i in expected_doys if (i not in retrieved_doys)]
        # Print a warning if there are days
        unexpected_doys = [i for i in retrieved_doys if i not in expected_doys]
        if len(unexpected_doys) > 0:
            pstr = 'WARNING: Unexpected doys have been downloaded -'
            print(pstr, unexpected_doys)
        # - Return a list of the failed doys
        if debug:
            print(failed_doys)
        return list(sorted(set(failed_doys)))
    else:
        pstr = "WARNING: No files found for GEOS5 so using expected doys (folder:{})"
        print(pstr.format(folder))
        return expected_doys


def get_latest_GEOSCF_fcast_data_lon_slice(dt=None, debug=True):
    """
    Wrapper to get lon slices from the latest GEOS-CF forecast data
    """
    # Hardcode options here for now.
    limit_lvls = False
    limit_lons = True
    # Call retrieval of forecast data
    get_latest_GEOSCF_fcast_data(dt=dt,
                                 limit_lvls=limit_lvls,
                                 limit_lons=limit_lons,
                                 #                                 limit_lats=limit_lats
                                 )


def get_latest_GEOSCF_fcast_data_lat_slice(dt=None, debug=True):
    """
    Wrapper to get lat slices from the latest GEOS-CF forecast data
    """
    # Hardcode options here for now.
    limit_lvls = False
    limit_lons = False
    limit_lats = True
    # Call retrieval of forecast data
    get_latest_GEOSCF_fcast_data(dt=dt,
                                 limit_lvls=limit_lvls,
                                 limit_lons=limit_lons,
                                 limit_lats=limit_lats
                                 )


def get_latest_GEOS5_fcast_data(dt=None, mode='fcast',
                                collection='inst3_3d_aer_Np', doys2use=None,
                                vars2use=None,
                                #                                collection='inst3_3d_aer_Nv',
                                #                                dt=None,
                                retry_download_if_data_not_present=True,
                                rm_existing_file=False,
                                debug=True):
    """
    Get the latest GEOS-5 forecast data (5-day started @ noon)
    """
    # Use yesterday at noon if no date provided
    if isinstance(dt, type(None)):
        TNow = AC.time2datetime([gmtime()])[0]
        dt = AC.add_days(TNow, -1)
        dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
    # Setup a text string with for dt
    dstr = dt.strftime('%Y/%m/%d %H:%M')
    # NOTE: this is current dealt with in get_GEOS5_as_ds_via_OPeNDAP (in fcast mode)
    try:
        ds = AC.get_GEOS5_as_ds_via_OPeNDAP(dt=dt, mode=mode,
                                            collection=collection)
        # Temporarily do not allow passing of date to OPeNDAP fetcher
#        ds = AC.get_GEOS5_as_ds_via_OPeNDAP(dt=None, mode=mode, collection=collection)
    # Check it has files...
    except OSError:
        print('GEOS5 - Failure loading data for ', dt)
        try:
            dt = AC.add_days(dt, -1)
            print('GEOS5 - Looking on the day before', dt)
            ds = AC.get_GEOS5_as_ds_via_OPeNDAP(dt=dt, mode=mode,
                                                collection=collection)
            # Temporarily do not allow passing of date to OPeNDAP fetcher
#            ds = AC.get_GEOS5_as_ds_via_OPeNDAP(dt=None, mode=mode, collection=collection)
        except:
            print('WARNING: FAILED to find GEOS5 latest dt... stopping now.')
    # Get initial date in forecast
    dt = AC.dt64_2_dt([ds['time'].values[0]])[0]
    # What would the folder be called for this noon forecast?
    NASA_data = get_local_folder('NASA_data')
    folder = NASA_data + 'GEOS_5/ARNA/fcast/'
    folderstr = '{}/{}_{:0>2}_{:0>2}_{:0>2}z/data.{}/'
    folder = folderstr.format(folder, dt.year, dt.month, dt.day, dt.hour,
                              collection)
    # Check if the folder is present, if not create it
    if os.path.isdir(folder):
        print('WARNING: folder already exists')
    else:
        os.makedirs(folder)
        print('Created folder for data ({})'.format(folder))

    # Set the variables to extract
    if isinstance(vars2use, type(None)):
        # Make the default variable du for inst3_3d_aer_Np (only one variabel regardless)
        if collection == 'inst3_3d_aer_Np':
            vars2use = ['du']
    # Check if the requested date is present.
    dt_is_lastest_run = is_dt_latest_GEOS5(dt=dt, collection=collection,
                                           mode=mode)
    override_GEOS5_check = True
    if dt_is_lastest_run or override_GEOS5_check:
        # Download this data if not already present
        get_GEOS5_data_cubes4collection(ds=ds, mode=mode,
                                        collection=collection,
                                        folder=folder,
                                        vars2use=vars2use, doys2use=doys2use,
                                        rm_existing_file=rm_existing_file,
                                        )

    elif retry_download_if_data_not_present:
        # Minutes to wait
        times2wait = 5
        min2wait = 4
        pstr = 'WARNING: data not availible ({}), so trying again every {} minutes ({}x)'
        print(pstr.format(dstr, min2wait, times2wait))
        for time2wait in range(times2wait):
            pstr = 'WARNING: Waiting for {:>2} min @ {}, then re-attempting download'
            TNow = AC.time2datetime([gmtime()])[0]
            print(pstr.format(min2wait, TNow.strftime('%Y/%m/%d %H:%M')))
            time.sleep(min2wait*60)
            print('Finished waiting @ {}'.format(TNow.strftime('%Y/%m/%d %H:%M')))
            # Check if the data is availible
            dt_is_lastest_run = is_dt_latest_GEOS5(dt)
            if dt_is_lastest_run:
                print('Data now available! Attempting download.')
                # Download this data if not already present
                get_GEOS5_data_cubes4collection(ds=ds, mode=mode,
                                                collection=collection,
                                                folder=folder,
                                                vars2use=vars2use,
                                                doys2use=doys2use,
                                                rm_existing_file=True,
                                                )

                print('WARNING: Stopping. Completed attempt at data download')
                sys.exit()
            else:
                print('Data still not available...')
    else:
        print('ERROR: Data not available for date so cannot attempt download')

    # Force a garbage clean up
    gc.collect()
    # - Regrid the files if all are present
    print('Checking if files sucessfully downloaded. If so, regridding them!')
    #
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                     collection=collection)
    failed_doys = check4GEOS5_failed_downloads(dt=dt, folder=folder)
    # If there are no failed days, then regrid
    if len(failed_doys) == 0:
        print('Attempting to regrid GEOS5 files')
        regrid_GEOS5_files_in_folder(folder=folder)
        print('Regridded GEOS5 files')
    else:
        print('Not all files present so, ')


def check4failed_downloads(dt=None, retry=True, n_doys=6,
                           inc_GEOS5_Nv_collection=False,
                           debug=False):
    """
    Check for failed downloads of GEOS-5/GEOS-CF files, then retry getting these
    """
    # Just use noon yesterday if no date tis provided
    if isinstance(dt, type(None)):
        TNow = AC.time2datetime([gmtime()])[0]
        dt = AC.add_days(TNow, -1)
        dt = datetime.datetime(dt.year, dt.month, dt.day, 12)

    # Local variables
    dt_fmt = dt.strftime('%Y/%m/%d %H:%M')
    sum_n_failed_doys = 0
    # What are the expected doys?
    expected_dates = [AC.add_days(dt, i) for i in range(n_doys)]
    dfT = pd.DataFrame(index=expected_dates)
    expected_doys = list(dfT.index.dayofyear.values)
    # - GEOS5 downloads
    collection = 'inst3_3d_aer_Np'
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                     collection=collection)
    failed_doys = check4GEOS5_failed_downloads(dt=dt, folder=folder)
    n_failed_doys = len(failed_doys)
    sum_n_failed_doys += n_failed_doys
    pstr = 'GEOS5 - there were {} failed dowloads on {}'
    print(pstr.format(n_failed_doys, dt_fmt))
    if (n_failed_doys > 0):
        print('The failed days were: ', failed_doys)
    if (n_failed_doys > 0) and retry:
        print('')
        get_latest_GEOS5_fcast_data(dt=dt, doys2use=failed_doys,
                                    rm_existing_file=True)

    # Also download the GEOS5 Nv collection
    if inc_GEOS5_Nv_collection:
        collection = 'inst3_3d_aer_Nv'
        folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                         collection=collection)
        vars2use = ['du{:0>3}'.format(i) for i in range(1, 6)]
        for var2use in vars2use:
            failed_doys = check4GEOS5_failed_downloads(dt=dt, folder=folder,
                                                       var2use=var2use)
            n_failed_doys = len(failed_doys)
    #        sum_n_failed_doys += n_failed_doys
            pstr = 'GEOS5 - there were {} failed dowloads for {} on {}'
            print(pstr.format(n_failed_doys, var2use, dt_fmt))
            if (n_failed_doys > 0):
                print('The failed days were: ', failed_doys)
            if (n_failed_doys > 0) and retry:
                print('')
                get_latest_GEOS5_fcast_data(dt=dt, doys2use=failed_doys,
                                            collection=collection,
                                            vars2use=[var2use],
                                            rm_existing_file=True)

    # - GEOSCF downloads
    collection = 'chm_inst_1hr_g1440x721_p23'
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                     collection=collection)
    vars2use = list(convert_GEOSCF_var2GEOSChem_name(rtn_dict=True).keys())
    # First check for lon_slice
    limit_lats = False
    limit_lons = True
    limit_lvls = False
    # Loop by variable
    for var2use in vars2use:
        failed_doys = check4GEOSCF_failed_downloads(dt=dt, folder=folder,
                                                    var2use=var2use,
                                                    limit_lons=limit_lons,
                                                    limit_lats=limit_lats,
                                                    limit_lvls=limit_lvls,
                                                    debug=debug)
        n_failed_doys = len(failed_doys)
        sum_n_failed_doys += n_failed_doys
        pstr = "GEOSCF (lon slice): {} failed downloads for '{}' on {}"
        print(pstr.format(n_failed_doys, var2use, dt_fmt))
        if (n_failed_doys > 0):
            print('The failed days were: ', failed_doys)
        if (n_failed_doys > 0) and retry:
            print('')
            get_latest_GEOSCF_fcast_data(dt=dt,
                                         vars2use=[var2use],
                                         doys2use=failed_doys,
                                         limit_lvls=limit_lvls,
                                         limit_lats=limit_lats,
                                         limit_lons=limit_lons,
                                         rm_existing_file=True,
                                         debug=debug)

    # Then lat_slice
    limit_lats = True
    limit_lons = False
    limit_lvls = False
    # Loop by variable
    for var2use in vars2use:
        failed_doys = check4GEOSCF_failed_downloads(dt=dt, folder=folder,
                                                    var2use=var2use,
                                                    limit_lons=limit_lons,
                                                    limit_lats=limit_lats,
                                                    limit_lvls=limit_lvls,
                                                    debug=debug)
        n_failed_doys = len(failed_doys)
        sum_n_failed_doys += n_failed_doys
        pstr = "GEOSCF (lat slice) - {} failed downloads for '{}' on {}"
        print(pstr.format(n_failed_doys, var2use, dt_fmt))
        if (n_failed_doys > 0):
            print('The failed days were: ', failed_doys)
        # If the doys do not line up... delete the files and start again

        # Download missing files
        if (n_failed_doys > 0) and retry:
            print('')
            get_latest_GEOSCF_fcast_data(dt=dt,
                                         vars2use=[var2use],
                                         doys2use=failed_doys,
                                         limit_lvls=limit_lvls,
                                         limit_lats=limit_lats,
                                         limit_lons=limit_lons,
                                         rm_existing_file=True,
                                         debug=debug)

    # Then alt_slice
    limit_lats = False
    limit_lons = False
    limit_lvls = True
    for var2use in vars2use:
        failed_doys = check4GEOSCF_failed_downloads(dt=dt, folder=folder,
                                                    var2use=var2use,
                                                    limit_lons=limit_lons,
                                                    limit_lats=limit_lats,
                                                    limit_lvls=limit_lvls,
                                                    debug=debug)
        n_failed_doys = len(failed_doys)
        sum_n_failed_doys += n_failed_doys
        pstr = "GEOSCF (alt slice) - {} failed dowloads for '{}' on {}"
        print(pstr.format(n_failed_doys, var2use, dt_fmt))
        if (n_failed_doys > 0):
            print('The failed days were: ', failed_doys)
        if (n_failed_doys > 0) and retry:
            print('')
            get_latest_GEOSCF_fcast_data(dt=dt,
                                         vars2use=[var2use],
                                         doys2use=failed_doys,
                                         limit_lvls=limit_lvls,
                                         limit_lats=limit_lats,
                                         limit_lons=limit_lons,
                                         rm_existing_file=True,
                                         debug=debug)
    # Return the total number of failed download doys
    return sum_n_failed_doys


def is_dt_latest_GEOSCF(dt=None, mode='fcast',
                        collection='chm_inst_1hr_g1440x721_p23'):
    """
    Check if the latest GEOS-CF run is available on OPeNDAP
    """
    # Just use noon yesterday if no date tis provided
    if isinstance(dt, type(None)):
        TNow = AC.time2datetime([gmtime()])[0]
        dt = AC.add_days(TNow, -1)
        dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
    # Check  containers
#    ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode, date=dt)
    # Temporarily do not allow passing of date to OPeNDAP fetcher
    ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode,
                                         date=None)
    last_start_date = AC.dt64_2_dt([ds.time.values[0]])[0]
    # Return check of files
    if last_start_date == dt:
        return True
    else:
        return False


def is_dt_latest_GEOS5(dt=None, mode='fcast',
                       collection='chm_inst_1hr_g1440x721_p23'):
    """
    Check if the latest GEOS-5 run is available on OPeNDAP
    """
    # Just use noon yesterday if no date tis provided
    if isinstance(dt, type(None)):
        TNow = AC.time2datetime([gmtime()])[0]
        dt = AC.add_days(TNow, -1)
        dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
    # Check  containers
    # Temporarily do not allow passing of date to OPeNDAP fetcher
    ds = AC.get_GEOS5_as_ds_via_OPeNDAP(collection=collection, mode=mode,
                                        dt=dt)
    last_start_date = AC.dt64_2_dt([ds.time.values[0]])[0]
    # Return check of files
    if last_start_date == dt:
        return True
    else:
        return False


def download_GEOSCF_fcast_data4date(dt=None, ds=None, mode='fcast',
                                    collection='chm_inst_1hr_g1440x721_p23',
                                    vars2use=None, doys2use=None,
                                    rm_existing_file=False,
                                    limit_lvls=True, limit_lons=False,
                                    limit_lats=False):
    """
    Retrieve the latest GEOS-CF data
    """
    # Get a pointer at the latest forecast data
    if isinstance(ds, type(None)):
        #        ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode, date=dt)
        ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode,
                                             date=None)
    # use the latest forecast data in the dataset if one not provided
    if isinstance(dt, type(None)):
        dt = AC.dt64_2_dt([ds.time.values[0]])[0]
    pstr = 'Attempting to download GEOS-CF data ({}) starting on {}'
    print(pstr.format(collection, dt.strftime('%Y/%m/%d %H:%M')))
    # Where to save the data
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                     collection=collection)
    # Check if the folder is present, if not create it
    if os.path.isdir(folder):
        print('WARNING: folder already exists ({})'.format(folder))
        print('(IF DATA IS PRESENT IT WILL BE OVERWRITEN - CHECK THIS!)')
        pstr = 'SETTINGS: limit_lvls={}; limit_lons={}; limit_lats={}'
        print(pstr.format(limit_lvls, limit_lons, limit_lats))
    else:
        os.makedirs(folder)
        print('Created folder for data ({})'.format(folder))
    # - Get the 3D data
    # Now retrieve the data -
    get_GEOSCF_data_cubes4collection(ds=ds, mode=mode, collection=collection,
                                     folder=folder, vars2use=vars2use,
                                     doys2use=doys2use, dt=dt,
                                     limit_lvls=limit_lvls,
                                     limit_lons=limit_lons,
                                     rm_existing_file=rm_existing_file,
                                     limit_lats=limit_lats)

    # - Get the 2D data
#    collection = 'chm_inst_1hr_g1440x721_p23'
#    ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode)
#    get_GEOSCF_data_cubes4collection(ds=ds, mode=mode, collection=collection)


def when_should_GEOSCF_have_last_run_from():
    """
    GEOS-CF is run every midnight
    """
    # Time now (GMT)
    TNow = AC.time2datetime([gmtime()])[0]
    # Therefore most recent forecast started when?
    dt = datetime.datetime(TNow.year, TNow.month, TNow.day, 12, )
    # Return as datetime
    return dt


def what_is_latest_data_locally(only_check_last5dates=True, debug=True):
    """
    Check what the latest data is locally
    """
    # Where to look for data
    folder = get_local_folder('NASA_data')
    folder += '/GEOS_CF/ARNA/fcast/'
    # What is the latest data that has been?
    subfolders = list(sorted(glob.glob(folder+'/*')))
    # Only look at the data in the last 5 folders for expendiency
    if only_check_last5dates:
        subfolders = subfolders[-5:]
    # Convert the subfolder names into datetimes
    dates4folders = [i.split(folder)[-1].split('_') for i in subfolders]
    print(dates4folders)
#    dates4folders = [i for i in dates4folders if ('z' not in i) ]
    dates4folders = [[int(ii) for ii in i[:3]] for i in dates4folders]
    dates4folders = [datetime.datetime(*i) for i in dates4folders]
    # Does the folder actually contain the data?
    CF_data_present_list = []
    for n_folder, folders in enumerate(subfolders):
        date = dates4folders[n_folder]
        if debug:
            pstr = "Checking data in folder for date '{}': {}"
            print(pstr.format(date, folder))
        CF_data_present_list += [is_GEOSCF_data_in_folder(folder)]
    #
    set_of_bools = list(set(CF_data_present_list))
    if len(set_of_bools) == 1:
        if set_of_bools[0] == True:
            if debug:
                print('Data present locally for all dates checked!')
            return dates4folders[-1]
        else:
            pstr = 'WARNING: No local data present 4 the last 5 dates checked!'
            print(pstr)
            return None


def is_GEOSCF_data_in_folder(folder):
    """
    check if GEOS-CF data is in the folder
    """
    # TODO - write a checker for GEOS-CF
    return False


def get_latest_GEOS5_diagnostics(dt=None,
                                 add_circles2GMAO_spatial_maps=True,
                                 verbose=False, debug=False):
    """
    Get the latest GEOS5 plots (datagrams, weather maps etc)
    """
    # Local variables
    prefix = 'ARNA'
    # If no day given, use previous day
    if isinstance(dt, type(None)):
        # Just use yesterday for Now
        TNow = AC.time2datetime([gmtime()])[0]
        dt = AC.add_days(TNow, -1)
        dt = datetime.datetime(dt.year, dt.month, dt.day, 12)
    # Folder to dave plots
    G5_folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                        inc_collection=False)
    folder2save = G5_folder + '/plots/plots.GMAO/'
    # Check the folder exists and create if not
    # Check if the folder is present, if not create it
    if os.path.isdir(folder2save):
        if verbose:
            print('WARNING: folder already exists ({})'.format(folder2save))
    else:
        os.makedirs(folder2save)
        if verbose:
            print('Created folder for data ({})'.format(folder2save))

    #  - Get the latest datagram plots
    plts2get = [
        # Dust at Cape Verde
        'du_16.7_-23.0', 'dumass_16.7_-23.0',
        # Nitrate, CO, BC at Cape Verde
        'nimass_16.7_-23.0', 'bc_16.7_-23.0', 'bcmass_16.7_-23.0',
        # Total aerosol at Cape Verde
        'total_16.7_-23.0', 'co_16.7_-23.0',
        # Dust and total at Dakar
        'du_14.4_-17.0', 'dumass_14.4_-17.0', 'total_14.4_-17.0',
    ]
    AC.get_GEOS5_datagram_plots(folder=folder2save, dt=dt, plts2get=plts2get,
                                prefix=prefix)
    # - Get Weather maps
    ptype = 'wxmaps'
    taus = [i*24 for i in range(6)]
    # Get 3 hourly values for the 1st 72 hours
    taus += [0+(3*i) for i in range(3*9)]
    taus = list(sorted(set(taus)))
    field = 'precip'
    AC.get_GEOS5_online_diagnostic_plots(folder=folder2save, ptype=ptype,
                                         field=field, dt=dt, taus=taus,
                                         prefix=prefix)

    # - Get composition maps
    ptype = 'chem2d'
    taus = [i*24 for i in range(6)]
    # Get 3 hourly values for the 1st 72 hours
    taus += [0+(3*i) for i in range(3*9)]
    taus = list(sorted(set(taus)))
    fields = ['duaot', 'cobbaf', 'bcsmass', 'niaot', 'nismass', ]
    for field in fields:
        AC.get_GEOS5_online_diagnostic_plots(folder=folder2save, ptype=ptype,
                                             dt=dt, field=field, taus=taus,
                                             prefix=prefix)

    # - Add a circle over Dakar and Sao Vicente
    if add_circles2GMAO_spatial_maps:
        files = glob.glob(folder2save+'*.png')
        files = [i for i in files if ('chem2d' in i) or ('wxmaps' in i)]
        for file in files:
            if debug:
                print('Adding Circles to file: {}'.format(file))
            # Open the image
            image = PIL.Image.open(file)
            draw = ImageDraw.Draw(image)
            # Add circles
            r = 39
            coords_list = [(774, 410), (840, 434)]
            for coords in coords_list:
                x, y = coords
                # Calculate the locations
                leftUpPoint = (x-r, y-r)
                rightDownPoint = (x+r, y+r)
                twoPointList = [leftUpPoint, rightDownPoint]
                draw.ellipse(twoPointList, fill=None,
                             outline='Grey',
                             width=3)
            # Save resulting image
            image.save(file)
    dt_str = dt.strftime('%Y/%m/%d %H:%M')
    print('Finished get_latest_GEOS5_diagnostics for {}'.format(dt_str))


def get_GEOSCF4flightnum(flight_ID='C225', resample_data=True):
    """
    Get the extracted GEOS-CF flight data for a specific flight
    """
    # Where is the extract GEOS-CF data?
    ARNA_data = get_local_folder('ARNA_data')
    folder = '{}/GEOS-CF/extracted_planeflight/'.format(ARNA_data)
    # Extract the data for a specific flight
    file2use = glob.glob(folder + '*_{}.csv'.format(flight_ID))[0]
    df = pd.read_csv(file2use)
    # NOTE: this is a kludge for an early version of the file
#    df = df.T
#    new_header = df.iloc[0] #grab the first row for the header
#    df = df[1:] #take the data less the header row
#    df.columns = new_header #set the header row as the df header
    # Make the datatime the index and remove and unneeded columns
    df.index = df['Datetime'].values
    df.index = pd.to_datetime(df.index.values)
    # Add temperature in deg C
    df['TempK'] = df['T'].copy()
    df['T'] = df['T'].values - 273.15
    # Add NOx as combined NO and NO2
    df['NOx'] = df['NO'].values + df['NO2'].values
    # Include a variable of NOy where HNO3 is removed
    # NOy = no_no2_hno3_hno4_hono_2xn2o5_pan_organicnitrates_aerosolnitrates
    df['NOy-HNO3'] = df['NOy'].values - df['HNO3'].values
    # Include a variable of NOy where HNO3 is removed
    df['NOy-HNO3-PAN'] = df['NOy'].values - \
        df['HNO3'].values - df['PAN'].values
    # gas-phase (exc. PAN, HNO3, HNO4, Org-NIT, N2O5)
    df['NOy-Limited'] = df['NO'].values + df['NO2'].values + \
        df['HNO2'].values + df['NIT'].values + df['NITs'].values
    # Resample the data?
    if resample_data:
        df = df.resample('1T').mean()
    return df


def extract_GEOS54all_ARNA_flights(flight_nums=[], debug=True):
    """
    Extract GEOS-CF model data for all ARNA flights
    """
    # Which flights to plot
    if len(flights_nums) == 0:
        flights_nums = [
            217, 218, 219, 220, 221, 222, 223, 224, 225,
        ]
    flight_IDs = ['C{}'.format(i) for i in flight_nums]
    # Loop by flight and extract the files
    for flight_ID in flight_IDs:
        if debug:
            print(flight_ID)
        # Extract data for flight
        df = extract_GEOS54ARNA_flight(flight_ID=flight_ID)
        # Save to csv.
        filename = 'ARNA_flightpath_extracted_from_GEOSCF_{}.csv'
        filename = filename.format(flight_ID)
        folder = './'
        df.to_csv(folder+filename)
        del df
        gc.collect()


def extract_GEOS54all_ARNA_surface_dates(testing_mode=False, debug=True):
    """
    Extract GEOS-CF model data for all ARNA surface dates
    """
    testing_mode = False
    ARNA_data = get_local_folder('ARNA_data')
    folder = '{}/GEOS-CF/extracted_surface4CVAO/'.format(ARNA_data)
    # Which dates to use?
    d = {
        'ARNA-1': (datetime.datetime(2019, 8, 6),
                   datetime.datetime(2019, 9, 3)),
        'ARNA-Winter-1': (datetime.datetime(2019, 11, 26),
                          datetime.datetime(2019, 12, 14)),
        'ARNA-2': (datetime.datetime(2020, 2, 5),
                   datetime.datetime(2020, 2, 27)),
        #    'CVAO-ALL' : (datetime.datetime(2018, 1, 1), datetime.datetime.now()),
    }
    # Variables for CVAO location
    location = 'CVO'
    LON, LAT, ALT = AC.get_loc(location)
    hPa_ = 985  # Note: this is the lowest level saved!
    # Loop by campaign and store data
    for campaign in list(sorted(d.keys()))[::-1]:
        # Start and end of campaign
        sdate, edate = d[campaign]
        # Make a dummy dataframe of
        date_range = pd.date_range(start=sdate, end=edate, freq='T')
        df = pd.DataFrame(index=date_range)
        LonVar = 'lon'
        LatVar = 'lat'
        TimeVar = 'time'
        PressVar = 'hPa'
        AltVar = PressVar
        df[LonVar] = LON
        df[LatVar] = LAT
        df[TimeVar] = df.index.values
        df[PressVar] = hPa_
        # Extract the standard output for dates
        # (only hourly as this is the natively save resolution)
#            dfE = extract_ds4df_locs(ds=dsGCF, df=df,)
        #
        # Add date column to dataframe

        def datetime2date(x):
            return datetime.datetime(x.year, x.month, x.day)
        df['date'] = df.index.map(datetime2date)
        dates2use = list(set(df.resample('D').sum().index.values))
        # Loop by collections to extract
        FileStr = 'ARNA_CVAO_extracted_GEOSCF_surface_{}hPa_{}_{}_{}_{}'
        folder = './'
        mode = 'assim'
        collections2use = [
            # Below are the core chemistry related collections
            'aqc_tavg_1hr_g1440x721_v1',
            'chm_tavg_1hr_g1440x721_v1',
            # Also a higher-resolution chemistry collection
            'htf_inst_15mn_g1440x721_x1',
            # And more model related output
            #        'met_tavg_1hr_g1440x721_x1',
            #        'xgc_tavg_1hr_g1440x721_x1',
        ]
#        collections2use += ['chm_inst_1hr_g1440x721_p23', 'met_inst_1hr_g1440x721_p23']
        for collection in collections2use[::-1]:
            if testing_mode:
                dates2use = dates2use[:2]
            # Open the OPenDAP dataset, get the data_vars, then close
            ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection,
                                                 mode=mode)
            vars2extract = list(ds.data_vars)
            del ds
            # remove the gocart variable as this has different dimensions
            if 'pm25_rh35_gocar' in vars2extract:
                vars2extract.pop(vars2extract.index('pm25_rh35_gocar'))
            # Loop by date
            for date in dates2use:
                date_dt = AC.dt64_2_dt([date])[0]
                date_str = '{}{:0>2}{:0>2}'.format(date_dt.year, date_dt.month,
                                                   date_dt.day)
                # Loop by variable to extract
                for var2extract in vars2extract:
                    print(campaign, collection, date_str, var2extract)
                    # Check if file exists already, if not then download
                    filename = FileStr.format(hPa_, collection, date_str,
                                              var2extract, campaign,)
                    filename = AC.rm_spaces_and_chars_from_str(filename)
                    file_exists = os.path.isfile(folder+filename)
#                     try:
#                         file_size = os.path.getsize(folder+filename)
#                     except OSError:
#                         file_size = 0
                    if file_exists:
                        FileStr = 'File already present, so skipping - {}'
                        print(FileStr.format(filename))
                    else:
                        try:
                            df_tmp = df.loc[df['date'] == date, :]
                            print(df_tmp.shape)
                            # Get variables from collection
                            _df = AC.extract_GEOSCF_assim4df(df_tmp,
                                                             PressVar=PressVar,
                                                             LonVar=LonVar,
                                                             LatVar=LatVar,
                                                             TimeVar=TimeVar,
                                                             collection=collection,
                                                             resample_df2ds_freq=True,
                                                             vars2extract=[
                                                                 var2extract],
                                                             spatial_buffer=1.5,
                                                             debug=True)
                            # Save to csv.
                            _df.to_csv(folder+filename+'.csv')
                            del _df
                            gc.collect()
                        except RuntimeError:
                            FileStr = 'RuntimeError, so skipping file - {}'
                            print(FileStr.format(filename))

        # Also extract the expanded output for ARNA-2
        extract_GEOSCF = False
        if (campaign == 'ARNA-2') and extract_GEOSCF:
            # Get the model data for the days of surface campaign
            dsGCF = get_GEOS_assim_expanded_dataset4ARNA(dts=date_range)
            # Extract nearest point using xarray function
            dsGCF = dsGCF.sel(lon=LON, lat=LAT, lev=hPa_, method='nearest')
            dsGCF = dsGCF.squeeze()
            del dsGCF['lon']
            del dsGCF['lat']
            del dsGCF['lev']
            dsGCF = AC.save_ds2disk_then_reload(dsGCF)
            # Save to csv.
            filename = 'ARNA_CVAO_extracted_GEOSCF_surface_{}hPa_{}_{}.csv'
            filename = filename.format(hPa_, '_expanded', campaign)
            filename = AC.rm_spaces_and_chars_from_str(filename)
            folder = './'
            df.to_csv(folder+filename)
            del df
            gc.collect()
        del df
        gc.collect()


def extract_GEOS54ARNA_flight(flight_ID='C225'):
    """
    Extract ARNA flightpath from GEOS-CF assimilation data
    """
    # -  Get the measurement flight tracks
    # Manually set FAAM flight file to use for now...
    filename = 'core_faam_*_*_r*_{}_1hz.nc'.format(flight_ID.lower())
    folder = '{}/CEDA/v2020_05/'.format(get_local_folder('ARNA_data'))
    file2use = glob.glob(folder+filename)
    assert len(file2use) == 1, 'WARNING: more that one file found!'
    ds = xr.open_dataset(file2use[0])
    # Get a dataframe of the coordinates to extract
    df = get_coordinates_from_NetCDF_file(ds=ds, falt_var='PS_RVSM',
                                          convert_m2hPa=False)
    dt = AC.dt64_2_dt([df.index.values[0]])[0]
    # - Extract the flight tracks from the model
    #  Get the model data for the days near the flight
    dsGCF = get_GEOS_assim_expanded_dataset4ARNA(dts=[dt])
    # Extract for known locations
    dfE = extract_ds4df_locs(ds=dsGCF, df=df,)
    return dfE


def get_GEOS_assim_expanded_dataset4ARNA(dts=None, update_lvl_units=True):
    """
    Get GEOS-CF assimilation data (expanded) for the ARNA campaign
    """
    # Where is the expanded GEOS-CF data
    NASA_data = get_local_folder('NASA_data')
    folder = '{}/GEOS_CF/ARNA/assim/expanded_variable_files/'.format(NASA_data)
    file_str = '*.nc4'
    # Make a dataframe of available files
    files = list(sorted(glob.glob(folder+file_str)))
    dates = [i.split('.nc4')[0][-14:] for i in files]
    dates2dt = [datetime.datetime.strptime(i, '%Y%m%d_%H%Mz') for i in dates]
    df = pd.DataFrame({'files': files}, index=dates2dt)
    # If a date is given, only open dates a day before and after
    if not isinstance(dts, type(None)):
        # Use one dat before and after
        if len(dts) == 1:
            sdate = AC.add_days(dts[0], -1)
            edate = AC.add_days(dts[0], 1)
        else:
            sdate = dts[0]
            edate = dts[-1]
        # Limit data to one day before or after
        bool1 = (df.index >= sdate) & (df.index <= edate)
        df = df.loc[bool1, :]
    # Get a list of those to extract
    files = df['files'].values.astype(list)
    # Open the NetCDF files as a xr.Dataset
    ds = xr.open_mfdataset(files)
    # Update the units for lev
    update_lvl_units = True
    if update_lvl_units:
        HPa_l = AC.get_GEOSCF_vertical_levels(native_levels=True)
        attrs = ds.lev.attrs.copy()
        attrs['units'] = 'hPa'
        attrs['standard_name'] = 'Pressure'
        attrs['long_name'] = 'Pressure'
        ds.lev.attrs = attrs
        # NOTE: Only 39 levels were saved offline
        #       so we're extracting these using there pythonic (-1) index
#        ds.lev.values = [ HPa_l[int(i)] for i in ds.lev.values -1]
        ds = ds.assign(lev=[HPa_l[int(i)] for i in ds.lev.values - 1])
    # Return the updated ds
    return ds


def get_GEOSCF_data_cubes4collection(ds=None, region='Cape_Verde',
                                     doys2use=None,
                                     limit_lvls=True, limit_lons=False,
                                     limit_lats=True,
                                     vars2use=None, dt=None,
                                     collection='chm_inst_1hr_g1440x721_p23',
                                     rm_existing_file=False,
                                     dates2use=None, mode='fcast',
                                     folder=None):
    """
    Extract cubes of data for region of interest during campaign period
    """
    # Which variables to use
    if isinstance(vars2use, type(None)):
        #    vars2use = ("pm25du_rh35_gcc", "pm25_rh35_gcc", "noy", "co", "no2", "o3",)
        vars2use = list(convert_GEOSCF_var2GEOSChem_name(rtn_dict=True).keys())
    #    vars2use = ('pm25su_rh35_gcc', 'pm25ni_rh35_gcc', 'pm25soa_rh35_gc', 'pm25ss_rh35_gcc', 'so2')

    # - Get the data as dataset via OPeNDAP
    if not isinstance(ds, type(None)):
        #        ds = C2BR2F4CF_as_ds_via_OPeNDAP(collection=collection, mode=mode, date=dt)
        # Temporarily do not allow passing of date to OPeNDAP fetcher
        ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(collection=collection, mode=mode,
                                             date=None)

    # - Only consider dates that are in dates2use list (if provided)
    # Use the datetime components via pandas as not yet fully in xarray.
    if not isinstance(dates2use, type(None)):
        df = pd.DataFrame(index=ds['time'].values)
        df['date'] = df.index.date
        date_bool = df['date'].isin(dates2use).values
        ds = ds.isel(time=date_bool)
    # Check that the initial retrieved data is the same as the requested on
    dt_str = dt.strftime('%Y/%m/%d %H:%M')
    t0_CF = ds.time.values[0]
    t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
    t0_CF_str = AC.dt64_2_dt([t0_CF])[0].strftime('%Y/%m/%d %H:%M')
    if dt_str != t0_CF_str:
        pstr = 'WARNING: provided date ({}) is not 1st date ({}) in dataset!'
        print(pstr.format(dt_str, t0_CF_str))
        print('WARNING: this data is being saved here: {}'.format(folder))

    # - Subset the dataset by region
    if region == 'Cape_Verde':
        # only consider the region the the plane will fly through +2 days transport
        # 20W-55E, 40S-40N,
        bool1 = ((ds.lon >= -35) & (ds.lon <= 10)).values
        bool2 = ((ds.lat >= 0) & (ds.lat <= 60)).values
        # Cut by lon, then lat
        ds = ds.isel(lon=bool1)
        ds = ds.isel(lat=bool2)
    elif region == 'global':
        pass
    else:
        print('WARNING: exiting as no region set')
        sys.exit()

    # - Limit the vertical levels
    if limit_lvls:
        #        lvls = [1000, 850, 500]
        lvls = [1000, 900, 800, 700, 600, 500]
        bool_lvls = [i in lvls for i in ds['lev']]
        ds = ds.isel(lev=bool_lvls)
        extr_str = 'lvls_'+'_'.join([str(i) for i in lvls])
    else:
        extr_str = 'ALL_lvls'

    # - Limit the number of longitudes?
    if limit_lons:
        lons2check = [-18, -19.5, -21, -22.5, -24, -25.5]
        idx = [AC.find_nearest(ds.lon.values, i) for i in lons2check]
        lons2use = ds.lon.values[idx]
        bool_lons = [i in lons2use for i in ds['lon'].values]
        ds = ds.isel(lon=bool_lons)
        extr_str += '_lons_'+'_'.join([str(i) for i in lons2use])
    else:
        pass
    # - Limit the number of latitudes?
    if limit_lats:
        lats2check = [12, 13, 14, 15, 16, 17]
        idx = [AC.find_nearest(ds.lat.values, i) for i in lats2check]
        lats2use = ds.lat.values[idx]
        bool_lats = [i in lats2use for i in ds['lat'].values]
        ds = ds.isel(lat=bool_lats)
        extr_str += '_lats_'+'_'.join([str(i) for i in lats2use])

    # Setup a print strin
    if limit_lons and limit_lvls:
        pstr = 'sliced by lvls+lons'
    elif limit_lvls:
        pstr = 'sliced by lvls'
    elif limit_lons:
        pstr = 'sliced by lons'
    elif limit_lats:
        pstr = 'sliced by lats'
    else:
        pstr = ''

    # -  Save the data by variable and day of year reduce packet size of transfers
    # Where to save?
    if isinstance(folder, type(None)):
        folder = './'
    # Save format for data (region, year, doy, variable)
    filestr = 'ARNA_GEOSCF_{}_{}_{}_{:0>3}_{}_{}.nc'
    if isinstance(doys2use, type(None)):
        doys2use = list(sorted(set(ds['time.dayofyear'].values)))
    # Loop by var and save
    for var in vars2use:
        # Now loop and save files in single day units
        print("Getting data via OPeNDAP for var '{}' - {}".format(var, pstr))
        for doy in doys2use:
            doy_int = int(float(doy))
            pstr = "Getting data via OPeNDAP for doy '{}' ({})"
            print(pstr.format(doy_int, var))
            # Subset for a single day
            ds_tmp = ds[var].sel(time=ds['time.dayofyear'] == doy)
            # If there is data in the dataset, then download then download this
            if len(ds_tmp['time'].values) > 0:
                # Get the year for the doy
                year = list(set(ds_tmp['time.year'].values))[0]
                # Set the name of the file to save
                fname = filestr.format(collection, region, year, doy_int, var,
                                       extr_str)
                # now save this single file
                if rm_existing_file:
                    AC.rm_file(folder=folder, filename=fname)
                ds_tmp.to_netcdf(folder+fname)
                # Do some memory management
                del ds_tmp
                gc.collect()
            else:
                pstr = "WARNING: no data for '{}' on doy ({}) - CALLING STOP!"
                print(pstr.format(var, doy_int))
                pstr = 'Please check the correct doys are downloaded for {}'
                print(pstr.format(dt_str))
                print('doys downloading: ', doys2use)
                sys.exit()


def get_GEOS5_data_cubes4collection(ds=None, region='Cape_Verde',
                                    doys2use=None, mode='fcast',
                                    vars2use=None,
                                    collection='inst3_3d_aer_Nv',
                                    folder=None, dates2use=None,
                                    #                                    limit_lvls=False, limit_lons=False
                                    rm_existing_file=False,
                                    ):
    """
    Extract cubes of data for region of interest during campaign period
    """
    # Which variables to use
    if isinstance(vars2use, type(None)):
        # Below are the dust variables for inst3_3d_aer_Nv
        vars2use = ['du{:0>3}'.format(i) for i in range(1, 6)]
    # Get the data as dataset via OPeNDAP
    if isinstance(ds, type(None)):
        ds = AC.get_GEOS5_as_ds_via_OPeNDAP(collection=collection, mode=mode,
                                            date=dt)
        # Temporarily do not allow passing of date to OPeNDAP fetcher
#        ds = AC.get_GEOS5_as_ds_via_OPeNDAP(collection=collection, mode=mode, date=None)

    # - Only consider dates that are in dates2use list (if provided)
    # Use the datetime components via pandas as not yet fully in xarray.
    if not isinstance(dates2use, type(None)):
        df = pd.DataFrame(index=ds['time'].values)
        df['date'] = df.index.date
        date_bool = df['date'].isin(dates2use).values
        ds = ds.isel(time=date_bool)

    # - Subset the dataset by region
    if region == 'Cape_Verde':
        # only consider the region the the plane will fly through +2 days transport
        # 20W-55E, 40S-40N,
        # -35E to 10E
        # 0N to 60N
        bool1 = ((ds.lon >= -35) & (ds.lon <= 10)).values
        bool2 = ((ds.lat >= 0) & (ds.lat <= 60)).values
        # Cut by lon, then lat
        ds = ds.isel(lon=bool1)
        ds = ds.isel(lat=bool2)
    elif region == 'global':
        pass
    else:
        print('WARNING: exiting as no region set')
        sys.exit()

    # - Limit the vertical levels
    extr_str = 'ALL_lvls'
    # -  Save the data by variable and day of year reduce packet size of transfers
    # Where to save?
    if isinstance(folder, type(None)):
        folder = './'
    # Save format for data (region, year, doy, variable)
    filestr = 'ARNA_GEOS5_{}_{}_{}_{:0>3}_{}_{}.nc'
    if isinstance(doys2use, type(None)):
        doys2use = list(set(ds['time.dayofyear'].values))
    # Make an extra string to print what slicing the data has had
    pstr = '- getting GEOS5 data as cuboids'
    # Loop by var and save
    for var in vars2use:
        # Now loop and save files in single day units
        print("Getting data via OPeNDAP for var '{}' {}".format(var, pstr))
        for doy in doys2use:
            doy_int = int(float(doy))
            pstr = "Getting data via OPeNDAP for doy '{}' ({})"
            print(pstr.format(doy_int, var))
            # Subset for a single day
            ds_tmp = ds[var].sel(time=ds['time.dayofyear'] == doy)
            # Get the year for the day in the model output
            year = list(set(ds_tmp['time.year'].values))[0]
            # now save this single file
            savename = filestr.format(collection, region, year, doy_int, var,
                                      extr_str)
            print(folder+savename)
            # Remove the file with this name, if it already there
            if rm_existing_file:
                AC.rm_file(folder=folder, filename=savename)
            # Now save the file
            ds_tmp.to_netcdf(folder+savename)
            # Do some memory management
            del ds_tmp
            gc.collect()


def get_GEOSCF_data_cubes4campaign(year=2018, region='Cape_Verde',
                                   doys2use=None,
                                   vars2use=None, folder=None,
                                   limit_lvls=True, limit_lons=False):
    """
    Extract cubes of data for region of interest during campaign period
    """
    # Which variables to use
    if isinstance(vars2use, type(None)):
        #    vars2use = ("pm25du_rh35_gcc", "pm25_rh35_gcc", "noy", "co", "no2", "o3",)
        vars2use = list(convert_GEOSCF_var2GEOSChem_name(rtn_dict=True).keys())
    #    vars2use = ('pm25su_rh35_gcc', 'pm25ni_rh35_gcc', 'pm25soa_rh35_gc', 'pm25ss_rh35_gcc', 'so2')

    # Setup lists of days to use
    dates2use = get_dates4campaigns(year=year)
    print('Using dates:', dates2use)

    # - Get the data as dataset via OPeNDAP
    collection = 'chm_inst_1hr_g1440x721_p23'
    ds = get_GEOSCF_assimlation_ds(collection=collection)

    # - Subset the dataset by date - Only consider dates that are in dates2use list
    # Use the datetime components via pandas as not yet fully in xarray.
    df = pd.DataFrame(index=ds['time'].values)
    df['date'] = df.index.date
    date_bool = df['date'].isin(dates2use).values
    ds = ds.isel(time=date_bool)

    # - Subset the dataset by region
    if region == 'Cape_Verde':
        # only consider the region the the plane will fly through +2 days transport
        # 20W-55E, 40S-40N,
        bool1 = ((ds.lon >= -35) & (ds.lon <= 10)).values
        bool2 = ((ds.lat >= 0) & (ds.lat <= 60)).values
        # Cut by lon, then lat
        ds = ds.isel(lon=bool1)
        ds = ds.isel(lat=bool2)
    elif region == 'global':
        pass
    else:
        print('WARNING: exiting as no region set')
        sys.exit()

    # - Limit the vertical levels
    extr_str = ''
    if limit_lvls:
        #        lvls = [1000, 850, 500]
        lvls = [1000, 900, 800, 700, 600, 500]
        bool_lvls = [i in lvls for i in ds['lev']]
        ds = ds.isel(lev=bool_lvls)
        extr_str += 'lvls_'+'_'.join([str(i) for i in lvls])
    else:
        extr_str += 'ALL_lvls'
    # - Limit the longitudes levels
    if limit_lons:
        lons2check = [-18, -19.5, -21, -22.5, -24, -25.5]
        idx = [AC.find_nearest(ds.lon.values, i) for i in lons2check]
        lons2use = ds.lon.values[idx]
        bool_lons = [i in lons2use for i in ds['lon'].values]
        ds = ds.isel(lon=bool_lons)
        extr_str += '_lons_'+'_'.join([str(i) for i in lons2use])
    else:
        pass

    # -  Save the data by variable and day of year reduce packet size of transfers
    # Where to save?
    if isinstance(folder, type(None)):
        folder = './'
    # Save format for data (region, year, doy, variable)
    filestr = 'ARNA_GEOSCF_chm_inst_1hr_g1440x721_p23_{}_{}_{:0>3}_{}_{}.nc'
    if isinstance(doys2use, type(None)):
        doys2use = list(set(ds['time.dayofyear'].values))
    # Make an extra string to print what slicing the data has had
    if limit_lons and limit_lvls:
        pstr = 'sliced by lvls+lons'
    elif limit_lvls:
        pstr = 'sliced by lvls'
    elif limit_lons:
        pstr = 'sliced by lons'
    else:
        pstr = ''

    # Loop by var and save
    for var in vars2use:
        # Now loop and save files in single day units
        print("Getting data via OPeNDAP for var '{}' {}".format(var, pstr))
        for doy in doys2use:
            doy_int = int(float(doy))
            pstr = "Getting data via OPeNDAP for doy '{}' ({})"
            print(pstr.format(doy_int, var))
            # Subset for a single day
            ds_tmp = ds[var].sel(time=ds['time.dayofyear'] == doy)
            # now copy locally the data
#            .copy()
            # now save this single file
            savename = filestr.format(region, year, doy_int, var, extr_str)
            print(folder+savename)
            ds_tmp.to_netcdf(folder+savename)
            # Do some memory management
            del ds_tmp
            gc.collect()


def get_GEOS5_data_cubes4campaign(ds=None, year=2018, region='Cape_Verde',
                                  doys2use=None,
                                  vars2use=None, collection='inst3_3d_aer_Nv',
                                  limit_lvls=False, limit_lons=False):
    """
    Extract cubes of data for region of interest during campaign period
    """
    # Which variables to use
    if isinstance(vars2use, type(None)):
        vars2use = ['du{:0>3}'.format(i) for i in range(1, 6)]
    # Setup lists of days to use
    dates2use = get_dates4campaigns(year=year)
    print('Using dates:', dates2use)
    # - Get the data as dataset via OPeNDAP
    if isinstance(ds, type(None)):
        ds = AC.get_GEOS5_as_ds_via_OPeNDAP(collection=collection)
    # - Subset the dataset by date - Only consider dates that are in dates2use list
    # Use the datetime components via pandas as not yet fully in xarray.
    df = pd.DataFrame(index=ds['time'].values)
    df['date'] = df.index.date
    date_bool = df['date'].isin(dates2use).values
    ds = ds.isel(time=date_bool)

    # - Subset the dataset by region
    if region == 'Cape_Verde':
        # only consider the region the the plane will fly through +2 days transport
        # 20W-55E, 40S-40N,
        bool1 = ((ds.lon >= -35) & (ds.lon <= 10)).values
        bool2 = ((ds.lat >= 0) & (ds.lat <= 60)).values
        # Cut by lon, then lat
        ds = ds.isel(lon=bool1)
        ds = ds.isel(lat=bool2)
    elif region == 'global':
        pass
    else:
        print('WARNING: exiting as no region set')
        sys.exit()
    # State that all data (cuboids are being extracted)
    extr_str = 'ALL_lvls'

    # -  Save the data by variable and day of year reduce packet size of transfers
    # Where to save?
    if isinstance(folder, type(None)):
        folder = './'
    # Save format for data (region, year, doy, variable)
    filestr = 'ARNA_GEOS5_{}_{}_{}_{:0>3}_{}_{}.nc'
    if isinstance(doys2use, type(None)):
        doys2use = list(set(ds['time.dayofyear'].values))
    # Make an extra string to print what slicing the data has had
    if limit_lons and limit_lvls:
        pstr = 'sliced by lvls+lons'
    elif limit_lvls:
        pstr = 'sliced by lvls'
    elif limit_lons:
        pstr = 'sliced by lons'
    else:
        pstr = ''

    # Loop by var and save
    for var in vars2use:
        # Now loop and save files in single day units
        print("Getting data via OPeNDAP for var '{}' {}".format(var, pstr))
        for doy in doys2use:
            doy_int = int(float(doy))
            pstr = "Getting data via OPeNDAP for doy '{}' ({})"
            print(pstr.format(doy_int, var))
            # Subset for a single day
            ds_tmp = ds[var].sel(time=ds['time.dayofyear'] == doy)
            # Get the year for the day in the model output
            year = list(set(ds_tmp['time.year'].values))[0]
            # now save this single file
            savename = filestr.format(collection, region, year, doy_int, var,
                                      extr_str)
            print(folder+savename)
            ds_tmp.to_netcdf(folder+savename)
            # Do some memory management
            del ds_tmp
            gc.collect()


def get_data_surface4campaign(year=2018, region='Cape_Verde', doys=None):
    """
    Extract surface data for region of interest during campaign period
    """
    # - Local variables
    # Get useful variables from the  3D cubes
#    vars2use = "pm25du_rh35_gcc", "pm25_rh35_gcc", "noy", "co",  "no2", "o3"
    vars2use = list(convert_GEOSCF_var2GEOSChem_name(rtn_dict=True).keys())
    # Add extra surface species
#    vars2use += [ 'hno3', 'nit', 'no', 'n2o5' 'dst4', 'dst3', 'dst2', 'dst1',]
    extra_vars = [
        'xyle', 'dst2', 'hno4', 'pm25su_rh35_gcc', 'pm25ni_rh35_gcc', 'ocpi', 'eoh',
        'benz', 'so2', 'rcho', 'h2o2', 'ald2', 'pm25soa_rh35_gc', 'n2o5', 'tolu',
        'pm25_rh35_gocar', 'dst4', 'no', 'ocpo', 'alk4', 'nh3', 'nh4', 'prpe', 'c2h6',
        'pm25bc_rh35_gcc', 'salc', 'hno3', 'ch4', 'c3h8', 'dst3', 'sala', 'bcpo',
        'dst1', 'macr', 'hcho', 'soap', 'acet', 'isop', 'soas', 'pm25oc_rh35_gcc',
        'mvk', 'nit', 'bcpi', 'pan', 'mek', 'pm25ss_rh35_gcc'
    ]
    vars2use += extra_vars

    # Also get other detail from the 2D fields
    vars2use = vars2use[::-1]
    dates2use = get_dates4campaigns(year=year)
    # Add to filename that the data is just the model surface
    extr_str = 'surface'
    # - Get the data as dataset via OPeNDAP
    collection = 'chm_tavg_1hr_g1440x721_v1'
    ds = get_GEOSCF_assimlation_ds(collection=collection)
    # - Subset the dataset by date - Only consider dates that are in dates2use list
    # Use the datetime components via pandas as not yet full yin xarray.
    df = pd.DataFrame(index=ds['time'].values)
    df['date'] = df.index.date
    date_bool = df['date'].isin(dates2use).values
    ds = ds.isel(time=date_bool)

    # - Subset the dataset by region
    if region == 'Cape_Verde':
        # only consider the region the the plane will fly through +2 days transport
        # 20W-55E, 40S-40N,
        bool1 = ((ds.lon >= -35) & (ds.lon <= 10)).values
        bool2 = ((ds.lat >= 0) & (ds.lat <= 60)).values
        # Cut by lon, then lat
        ds = ds.isel(lon=bool1)
        ds = ds.isel(lat=bool2)
    elif region == 'global':
        pass
    else:
        print('WARNING: exiting as no region set')
        sys.exit()

    # -  Save the date in the earth0
    # Where to save the extracted data?
    if isinstance(folder, type(None)):
        folder = './'
    # Save format for data (region, year, doy, variable)
    filestr = 'ARNA_GEOSCF_{}_{}_{}_{:0>3}_{}_{}.nc'
    # Loop by var and save
#    for var in vars2use:
    # Now loop and save files in single day units
#    print("Getting data via OPeNDAP for var '{}'".format(var))
    if isinstance(doys, type(None)):
        doys = list(set(ds['time.dayofyear'].values))
    for doy in doys:
        doy_int = int(float(doy))
        print("Getting data via OPeNDAP for doy '{}'".format(doy_int))
        # Subset for a single day
        ds_tmp = ds[vars2use].sel(time=ds['time.dayofyear'] == doy)
        # Now copy locally the data
        ds_tmp = ds_tmp.copy()
        # Now save this single file
        savename = filestr.format(
            collection, region, year, doy_int, 'ALL', extr_str)
        ds_tmp.to_netcdf(folder+savename)
        # Do some memory management
        del ds_tmp
        gc.collect()


def get_GEOSCF_assimlation_ds(collection='chm_inst_1hr_g1440x721_p23',
                              mode='assim'):
    """
    Wrapper to get the GEOS-CF assimilation data
    """
    # Get entire archive as a xr.Dataset
    ds = AC.get_GEOSCF_as_ds_via_OPeNDAP(mode=mode, collection=collection)
    return ds


def add_extra_derived_species(ds, debug=False):
    """
    Add extra derived variables of interest
    """
    # Add NOy divided by PM2.5 concentration
#     NewVar = 'NOy/PM2.5(dust)'
#     ds[NewVar] = ds['NOy'] / ds[ 'PM2.5(dust)']
#     attrs = ds[NewVar].attrs
#     attrs['units'] = 'unitless'
#     ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration
    NewVar = 'NOy*PM2.5(dust)'
    ds[NewVar] = ds['NOy'] * ds['PM2.5(dust)']
    attrs = ds[NewVar].attrs
    attrs['units'] = 'unitless'
    ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration (only for where dust >50 ug/m3)
#     CopyVar =  'NOy'
#     NewVar = 'NOy/PM2.5(dust>50ug/m3)'
#     ds[NewVar] = ds[CopyVar].copy()
#     vals2use = ds[ 'PM2.5(dust)'].copy().values
#     vals2use[vals2use<50] = np.NaN
#     ds[NewVar].values = ds['NOy'].values / vals2use
#     attrs = ds[NewVar].attrs
#     attrs['units'] = 'unitless'
#     ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration (only for where dust >25 ug/m3)
    CopyVar = 'NOy'
    NewVar = 'NOy/PM2.5(dust>10ug/m3)'
    ds[NewVar] = ds[CopyVar].copy()
    vals2use = ds['PM2.5(dust)'].copy().values
    vals2use[vals2use < 10] = np.NaN
    ds[NewVar].values = ds['NOy'].values / vals2use
    attrs = ds[NewVar].attrs
    attrs['units'] = 'unitless'
    ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration (only for where dust >25 ug/m3)
    CopyVar = 'NOy'
    NewVar = 'NOy/PM2.5(dust>25ug/m3)'
    ds[NewVar] = ds[CopyVar].copy()
    vals2use = ds['PM2.5(dust)'].copy().values
    vals2use[vals2use < 25] = np.NaN
    ds[NewVar].values = ds['NOy'].values / vals2use
    attrs = ds[NewVar].attrs
    attrs['units'] = 'unitless'
    ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration (only for where dust >25 ug/m3)
    CopyVar = 'NOy'
    NewVar = 'NOy(where>10ug/m3 dust)'
    ds[NewVar] = ds[CopyVar].copy()
    NOy = ds[NewVar].values
    dust = ds['PM2.5(dust)'].values
    NOy[dust < 10] = np.NaN
    ds[NewVar].values = NOy
    attrs = ds[NewVar].attrs
    attrs['units'] = 'ppbv'
    ds[NewVar].attrs = attrs

    # Add NOy divided by PM2.5 concentration (only for where dust >0.5 ug/m3)
    CopyVar = 'NOy'
    NewVar = 'NOy/PM2.5(dust>0.5ug/m3)'
    ds[NewVar] = ds[CopyVar].copy()
    vals2use = ds['PM2.5(dust)'].copy().values
    vals2use[vals2use < 0.5] = np.NaN
    ds[NewVar].values = ds['NOy'].values / vals2use
    attrs = ds[NewVar].attrs
    attrs['units'] = 'unitless'
    ds[NewVar].attrs = attrs
    # print data variables if debugging...
    if debug:
        print('-'*10, [i for i in ds.data_vars], '-'*10)

    return ds

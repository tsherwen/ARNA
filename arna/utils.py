"""
Utility functions for ARNA campaign/project work
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
import matplotlib


def mk_core_plot_folders_then_mv2webfiles(dt=None, mv2webfiles=True,
                                          debug=True):
    """
    Make core folders (+?? hours), then move these to webfiles
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # - mv the +24/+48 files into the core folders
    dstr = dt.strftime('%Y/%m/%d %H:%M')
    copy_files2core_plot_folders(dt=dt)
    # - Now move the folder to webfiles
    if mv2webfiles:
        TNow = AC.time2datetime([gmtime()])[0]
        pstr = "Started moving files for {} to webfiles @ {}"
        print(pstr.format(dstr, TNow.strftime('%Y/%m/%d %H:%M')))
        mv_plots2webfiles(dt=dt)
    # - Now move the files to google drive
    TNow = AC.time2datetime([gmtime()])[0]
    pstr = "Started moving files for {} to google drive @ {}"
    print(pstr.format(dstr, TNow.strftime('%Y/%m/%d %H:%M')))
    # Move the files
    mv_plots2google_drive(dt=dt, debug=debug)
    # print that the job is finished.
    TNow = AC.time2datetime([gmtime()])[0]
    pstr = "Finished moving files for {} to google drive @ {}"
    print(pstr.format(dstr, TNow.strftime('%Y/%m/%d %H:%M')))


def which_plot_folders_are_not_complete4dt(dt):
    """
    Check which plots have been made for a specific datetime (dt)
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Get the root plot folder
    folder = get_GEOS_data_folder4dt(
        dt=dt, product='GEOS_5', inc_collection=False)
    folder += '/plots/'
    # Hardwire folder names to check
    subfolders2check = [
        'plots.GMAO',
        'alt_slice',
        'lon_slice',
        'lat_slice',
        'alt_slice.zoomed',
        'alt_slice.individual'
    ]
    if not os.path.isdir(folder):
        return subfolders2check
    else:
        # Setup a dataframe to store values
        df = pd.DataFrame()
        # Loop subfolders and get info on folders
        for subfolder in subfolders2check:
            # Evaluate full folder name
            subfolder_full = '{}/{}'.format(folder, subfolder)
            if os.path.isdir(subfolder_full):
                d = AC.get_stats_on_files_in_folder_as_dict(
                    folder=subfolder_full)
            else:
                d = {'#': 0}
            df[subfolder] = pd.Series(d)
        # - Now decide if all the file are present
    #    df = df.T
        # Hardwire number of expected files
        nfiles4ARNA = {
            'plots.GMAO': 169,
            'alt_slice': 246,
            'lon_slice': 246,
            'lat_slice': 246,
            'alt_slice.zoomed': 246,
            'alt_slice.individual': 1230,
            #    'core.plus_024H.2020_01_22_12_00': 61,
            #    'core.plus_024H.2020_01_22_12_00': 61,
            #    'core.plus_030H.2020_01_22_18_00': 61,
            #    'core.plus_048H.2020_01_23_12_00': 61,
            #    'core.plus_054H.2020_01_23_18_00': 61,
            #    'core.plus_072H.2020_01_24_12_00': 61,
            #    'core.plus_078H.2020_01_24_18_00': 55, # Note this should be 61!
        }
        # Check the right number of files are present
        correct_num_col = 'Correct #'
        for col in df.columns:
            correct_num = df[col]['#'] == nfiles4ARNA[col]
            df.loc[correct_num_col, col] = correct_num
        # just select the name of the
        names = df.loc[:, df.T[correct_num_col].values == False].columns
        return list(names)


def copy_files2core_plot_folders(dt=None, verbose=True, debug=False):
    """
    Copy the files to core plot folders (+24, +48) for a given date
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Setup a string for the date of the forecast
    dt0_Str = dt.strftime('%Y/%m/%d %H:%M')
    # Also setup datetime objects for the +24 fcast and +48 fcasts
    dt24 = AC.add_days(dt, 1)
    dt30 = AC.add_hrs(dt, 24+6)
    dt48 = AC.add_days(dt, 2)
    dt54 = AC.add_hrs(dt, 48+6)
    dt72 = AC.add_days(dt, 3)
    dt78 = AC.add_hrs(dt, 72+6)
    # Get the root plot folder
    folder = get_GEOS_data_folder4dt(
        dt=dt, product='GEOS_5', inc_collection=False)
    folder += '/plots/'

    # - Copy the files for the +24 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt.year, dt.month, dt24.day, dt24.hour, dt24.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(24, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt24.year, dt24.month, dt24.day, dt24.hour, dt24.minute)
    # Get a list of files to copy
    files = glob.glob('{}/*/*{}.png'.format(folder, dstr))
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_024*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # print some information about the files found
    pstr = 'Moving {:>3} files to {} from forecast started on {}'
    print(pstr.format(len(files), fstr, dt0_Str))
    # Loop by files and move
    for file in files:
        filename = file.split('/')[-1]
        if debug:
            print(filename, file)
        os.popen('cp {} {}'.format(file, dfolder+filename))

    # - Copy the files for the +30 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt30.year, dt30.month, dt30.day, dt30.hour, dt30.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(30, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    files = glob.glob('{}/*/*{}.png'.format(folder, dstr))
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_030*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # Loop by files and move
    for file in files:
        filename = file.split('/')[-1]
        if debug:
            print(filename, file)
        os.popen('cp {} {}'.format(file, dfolder+filename))

    # - Copy the files for the +48 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt48.year, dt48.month, dt48.day, dt48.hour, dt48.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(48, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    files = glob.glob('{}/*/*{}.png'.format(folder, dstr))
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_048*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # Loop by files and move
    for file in files:
        filename = file.split('/')[-1]
        if debug:
            print(filename, file)
        os.popen('cp {} {}'.format(file, dfolder+filename))

    # - Copy the files for the +54 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt54.year, dt54.month, dt54.day, dt54.hour, dt54.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(54, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    files = glob.glob('{}/*/*{}.png'.format(folder, dstr))
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_054*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # Loop by files and move
    for file in files:
        filename = file.split('/')[-1]
        if debug:
            print(filename, file)
        os.popen('cp {} {}'.format(file, dfolder+filename))

    # - Copy the files for the +72 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt72.year, dt72.month, dt72.day, dt72.hour, dt72.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(72, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    files = glob.glob('{}/*/*{}.png'.format(folder, dstr))
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_072*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # Loop by files and move
    for file in files:
        filename = file.split('/')[-1]
        if debug:
            print(filename, file)
        os.popen('cp {} {}'.format(file, dfolder+filename))

    # - Copy the files for the +78 folder
    # Setup strings for date, filename, and folder
    dstr = '{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr = dstr.format(dt78.year, dt78.month, dt78.day, dt78.hour, dt78.minute)
    fstr = 'core.plus_{:0>3}H.{}'.format(78, dstr)
    dfolder = '{}/{}/'.format(folder, fstr)
    # Get a list of files to copy
    files = glob.glob('{}/*/*{}.png'.format(folder, dstr))
    # Also add the met plots into the core plots
    gstr = '*atlantic*00_078*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # Copy across the core z+0 datagram for total AOD at CVAO too.
    gstr = '*datagram_total*'
    files += glob.glob('{}/*/*{}*'.format(folder, gstr))
    # Loop by files and move
    for file in files:
        filename = file.split('/')[-1]
        if debug:
            print(filename, file)
        os.popen('cp {} {}'.format(file, dfolder+filename))


def get_Google_drive_ID4ARNA_forecast_plots():
    """
    Retrieve the google drive folder ID for the ARNA fcast folder
    """
    # Hardcoded for now...
    # TODO: retrieve this based on a walk through directories
    GD_ID = '1u4TAwsefzfFJ5OSa0PwiMgF_sanNR1QI'
    return GD_ID


def mv_plots2google_drive(dt=None, verbose=True, debug=False):
    """
    Move plots to a google drive account in an automated manner
    """
    from pydrive.auth import GoogleAuth
    from pydrive.drive import GoogleDrive
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Make the datetime string to be used for folders
    dt_str = '{}_{:0>2}_{:0>2}_{:0>2}z'
    dt_str = dt_str.format(dt.year, dt.month, dt.day, dt.hour)
    # What are the root folders for the data/plots
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                     inc_collection=False)
    folder += '/plots/'
    # Set location to store data
    gauth = GoogleAuth()
    drive = GoogleDrive(gauth)
    # Authenticate the connections
    gauth.CommandLineAuth()
    # Get the ID for the core
    fcast_GD_folder_ID = get_Google_drive_ID4ARNA_forecast_plots()
    # Create a folder for the date (dt)
    folder_metadata = {
        #        'title' : dt_str+'_earth0_PyDrive_sync',
        'title': dt_str,
        #
        "parents":  [{"id": fcast_GD_folder_ID}],
        # The mimetype defines this new file as a folder
        'mimeType': 'application/vnd.google-apps.folder'
    }
    GD_folder_root = drive.CreateFile(folder_metadata)
    GD_folder_root.Upload()
    # Get info on the folder
    root_folder_title = GD_folder_root['title']
    root_folder_id = GD_folder_root['id']
    # Loop through files and directories to upload
    for path, directories, files in os.walk(folder):
        # What is the subroot for the folder?
        subfolder = path.split(dt_str)[-1].split('/plots/')[-1]
        if verbose:
            print(path, directories, files)
        # Create the subfolder for plots
        metadata = {
            'title': subfolder,
            "parents":  [{"id": root_folder_id}],
            "mimeType": "application/vnd.google-apps.folder"
        }
        GD_subfolder = drive.CreateFile(metadata)
        GD_subfolder.Upload()
        GD_subfolder_title = GD_subfolder['title']
        GD_subfolder_id = GD_subfolder['id']
        if verbose:
            print(subfolder, GD_subfolder_title, GD_subfolder_id)
        # Loop and upload by file
        for file in files:
            if debug:
                print(file)
            # Metadata for file
            metadata = {
                'title': file,
                "parents":  [{"id": GD_subfolder_id}],
            }
            # Set a file space
            GD_file = drive.CreateFile(metadata)
            # point at the required file
            GD_file.SetContentFile('{}/{}'.format(path, file))
            GD_file.Upload()  # Upload the file.


def mv_GEOS_data_from_earth2viking(dt=None, user=None):
    """
    Move the GEOS-5/CF data from earth0 to viking
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Which user to use?
    if isinstance(user, type(None)):
        import getpass
        user = getpass.getuser()
    # - Get the locations of the GEOS-5 data
    # Get the locations on
    host = 'viking'
    collection = 'inst3_3d_aer_Np'
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                     collection=collection,
                                     host=host)
    # Make the local folder if it is not present
    AC.mk_folder(folder=folder)
    # Get the remote folder
    host = 'earth0'
    remote_folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                            collection=collection,
                                            host=host)
    # Check that the remote folder is present?  and the files are correct?
    # TODO...
    # Call scp via os
    host = 'earth0.york.ac.uk'
    call_str = "scp -r {}@{}:{}/* {}".format(user, host, remote_folder, folder)
    os.system(call_str)

    # - Get the locations of the GEOS-5 data
    # Get the locations on
    host = 'viking'
    collection = 'chm_inst_1hr_g1440x721_p23'
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                     collection=collection,
                                     host=host)
    # Make the local folder if it is not present
    AC.mk_folder(folder=folder)
    # Get the remote folder
    host = 'earth0'
    remote_folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_CF',
                                            collection=collection,
                                            host=host)
    # Check that the remote folder is present?  and the files are correct?
    # TODO...
    # Call scp via os
    host = 'earth0.york.ac.uk'
    call_str = "scp {}@{}:{}/* {}".format(user, host, remote_folder, folder)
    os.system(call_str)


def mv_plots2webfiles(dt=None, debug=True):
    """
    Move plots to York's 'webfiles' server in an automated manner
    """
    import paramiko
    import sys
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Make the datetime string to be used for folders
    dt_str = '{}_{:0>2}_{:0>2}_{:0>2}z'
    dt_str = dt_str.format(dt.year, dt.month, dt.day, dt.hour)
    # What are the root folders for the data/plots
    folder = get_GEOS_data_folder4dt(dt=dt, product='GEOS_5',
                                     inc_collection=False)
    folder += '/plots/'
    # Define data location and ports etc
    host = "webfiles.york.ac.uk"  # hard-coded
    port = 22
    transport = paramiko.Transport((host, port))
    password = "niangw?ndi"  # hard-coded (NOTE: no longer valid now on GitHub)
    username = "chem631"  # hard-coded
    # Connect to webfiles
    transport.connect(username=username, password=password)
    # Open a sftp link
    sftp = paramiko.SFTPClient.from_transport(transport)
    # Set location to store data
    dt_root_path = './SHERWEN_TOMAS/ARNA/fcast/{}/'.format(dt_str)
    # Create a folder for the date
    mkdir_if_not_present(sftp=sftp, path=dt_root_path)
    # Loop through files and directories to upload
    for path, directories, files in os.walk(folder):
        # What is the subroot for the folder?
        subfolder = path.split(dt_str)[-1].split('/plots/')[-1]
        dt_sub_path = '{}/{}/'.format(dt_root_path, subfolder)
        if debug:
            print(path, directories, files)
            print('')
            print(dt_sub_path)
        # open connection
        sftp = paramiko.SFTPClient.from_transport(transport)
        # Create the path if it doesn't exist
        mkdir_if_not_present(sftp=sftp, path=dt_sub_path)
        # Close connection
        sftp.close()
        # Loop and upload by file
        for file in files:
            # open connection
            sftp = paramiko.SFTPClient.from_transport(transport)
            # tranfer file
            sftp.put('{}/{}'.format(path, file), dt_sub_path+file)
            # Close connection
            sftp.close()


def mkdir_if_not_present(sftp=None, path=None, verbose=True):
    """
    Make a remote path via sftp if not present
    """
    try:
        sftp.chdir(path)  # Test if path exists
        if verbose:
            pstr = 'WARNING: remote path exists, so folder not created - {}'
            print(pstr.format(path))
    except IOError:
        sftp.mkdir(path)  # Create path
        sftp.chdir(path)
        if verbose:
            print('CREATED remote path that did not exist - {}'.format(path))


def mk_folder_structure4fcast(dt=None, mode='fcast', verbose=False):
    """
    Make folders required for a date forecast data & plots 4 ARNA
    """
    # Use yesterday's forecast at noon if others not available.
    if isinstance(dt, type(None)):
        Tnow = AC.time2datetime([gmtime()])[0]
        # Get the 5-day forecast at noon...
        dt = datetime.datetime(Tnow.year, Tnow.month, Tnow.day, 12, )
        # Use yesterday
        dt = AC.add_days(dt, -1)
    # Check that the forecast date is for the noon forecast
    assert dt.hour == 12, 'WARNING forecast date not starting at noon!'
    # Make the datetime string to be used for folders
    dt_str = '{}_{:0>2}_{:0>2}_{:0>2}z'
    dt_str = dt_str.format(dt.year, dt.month, dt.day, dt.hour)
    # Also setup datetime objects for the +24 fcast and +48 fcasts
    dt24 = AC.add_days(dt, 1)
    dt30 = AC.add_hrs(dt, 24+6)
    dt48 = AC.add_days(dt, 2)
    dt54 = AC.add_hrs(dt, 48+6)
    dt72 = AC.add_days(dt, 3)
    dt78 = AC.add_hrs(dt, 72+6)
    # What are the root folders for the data/plots
    NASA_data = get_local_folder('NASA_data')
    G5_folder = NASA_data + 'GEOS_5/ARNA/'
    GCF_folder = NASA_data + 'GEOS_CF/ARNA/'
    # Setup a list to populate with folders to make
    folders = []
    # - Setup folder structure for GEOS-CF structure
    # Setup folder for the data
    folders += ['{}/{}/{}'.format(GCF_folder, mode, dt_str)]
    # - Setup folder structure for GEOS-5 structure
    # Setup folder for the data
    folders += ['{}/{}/{}'.format(G5_folder, mode, dt_str)]
    # Setup folder for the plots
    ex_str = 'plots'
    plots_root = '{}/{}/{}/{}'.format(G5_folder, mode, dt_str, ex_str)
    folders += [plots_root]
    # those directly snapped from GMAO
    ex_str = 'plots.GMAO'
    folders += ['{}/{}'.format(plots_root, ex_str)]
    # Make strings for folders for +24H and +48H
    dstr = 'core.plus_{:0>3}H.{}_{:0>2}_{:0>2}_{:0>2}_{:0>2}'
    dstr24 = dstr.format(24, dt24.year, dt24.month,
                         dt24.day, dt24.hour, dt24.minute)
    dstr30 = dstr.format(30, dt30.year, dt30.month,
                         dt30.day, dt30.hour, dt30.minute)
    dstr48 = dstr.format(48, dt48.year, dt48.month,
                         dt48.day, dt48.hour, dt48.minute)
    dstr54 = dstr.format(54, dt54.year, dt54.month,
                         dt54.day, dt54.hour, dt54.minute)
    dstr72 = dstr.format(72, dt72.year, dt72.month,
                         dt72.day, dt72.hour, dt72.minute)
    dstr78 = dstr.format(78, dt78.year, dt78.month,
                         dt78.day, dt78.hour, dt78.minute)
    # List of subfolders to make
    subfolders = [
        'alt_slice', 'lon_slice', 'lat_slice', 'alt_slice.zoomed',
        'alt_slice.individual', dstr24, dstr30, dstr48, dstr54, dstr72, dstr78,
    ]
    folders += ['{}/{}'.format(plots_root, i) for i in subfolders]
    # - Now loop and make folders
    for folder2mk in folders:
        AC.mk_folder(folder=folder2mk, verbose=verbose)


def select_smaller_area_around_Cape_Verde(ds):
    """
    Sub-select a smaller region just around Cape Verde
    """
    # What is the smaller area around Cape Verde defined as?
#    lat_min = 15 # used for initial work
    lat_min = 14
    lat_max = 18
    lon_max = -21
    lon_min = -26
    # Create booleans to sub select
    bool1 = ((ds.lon >= lon_min) & (ds.lon <= lon_max)).values
    bool2 = ((ds.lat >= lat_min) & (ds.lat <= lat_max)).values
    # Cut by lon, then lat
    ds = ds.isel(lon=bool1)
    ds = ds.isel(lat=bool2)
    return ds


def update_GEOSCF_names2GEOSChem_names(ds):
    """
    Update the GEOS-CF names to be in GEOS-Chem form
    """
    # Get a dictionary of the GEOS-CF and GEOS-Chem names
    d = convert_GEOSCF_var2GEOSChem_name(rtn_dict=True)
    return ds.rename(name_dict=d)


def convert_GEOSCF_var2GEOSChem_name(input=None, rtn_dict=False):
    """
    Convert the GEOS-CG species name to GEOS-Chem name
    """
    d = {
        'noy': 'NOy', 'co': 'CO', 'no2': 'NO2', 'o3': 'O3',
        'pm25_rh35_gcc': 'PM2.5',
        'pm25du_rh35_gcc': 'PM2.5(dust)',
        #    'pm25su_rh35_gcc' : 'PM2.5(SO4)',
        #    'pm25ni_rh35_gcc' : 'PM2.5(NIT)',
        #    'pm25soa_rh35_gc' : 'PM2.5(SOA)',
        #    'pm25ss_rh35_gcc' : 'PM2.5(SSA)',
        #    'pm25bc_rh35_gcc' : 'PM2.5(BC)',
        #    'pm25oc_rh35_gcc' : 'PM2.5(OC)',
        'so2': 'SO2',
    }
    if rtn_dict:
        return d
    else:
        return d[input]


def convert_GEOSCF_units(ds, debug=False):
    """
    Convert from GEOS-CF units to GEOS-Chem units
    """
    vars2use = [i for i in ds.data_vars]
    for var in vars2use:
        # Get standard conversion units
        try:
            units, scaleby = AC.tra_unit(var, scale=True)
            if debug:
                print(var, units, scaleby)
        except KeyError:
            print('WARNING: setting {} to unitless!'.format(var))
            scaleby = 1
            units = 'unitless'
        # update the data by the scaling
        ds[var].values = ds[var].values.copy() * scaleby
        # updates units
        attrs = ds[var].attrs.copy()
        attrs['units'] = units
        ds[var].attrs = attrs
        # Do some garbage cleaning
        gc.collect()
    return ds


def mk_video_from_plts(FileStr=None):
    """
    Convert images frames into videos using ffmpeg
    """
    FileStr = 'spatial_plot_NOy_PM2_5_dust__lev_1000_0_dt_*.png'
    #
    print('WARNING: ffmpeg calls here have been switched off')
#    ffmpegstr = "ffmpeg -framerate 3 -pattern_type glob -i '{FileStr}' -c:v libx264 -pix_fmt yuv420p out.mp4"

# ffmpeg -framerate 3 -pattern_type glob -i 'spatial_plot_NOy_PM2_5_dust__lev_500_0_dt_2019_11_*.png' -c:v libx264 -pix_fmt yuv420p out_2019_11_06_5day_fcast_500hPa.mp4
    os.system(" ")


def get_reduced_cmap(cmap='Reds', minval=0.0, maxval=0.75, npoints=100):
    """
    Get a reduced colormap object (cmap)
    """
    # Load color map
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    # Make a truncated colour map
    Nstr = 'trunc({n},{a:.2f},{b:.2f})'
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        Nstr.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval,  npoints))
    )
    return cmap


def convert_decimcal_degress2nav_format(degrees_dec, just_degrees_minutes=True,
                                        LaTeX_format=True, debug=False):
    """
    convert decimal degree to degrees, minutes and then decimal seconds
    """
    # Setup a list to store values
    vals2rtn = []
    # Loop by degrees provided
    for deg in degrees_dec:
        degrees = int(deg)
        degrees_frac = abs(AC.myround(deg, base=0.01, integer=False) - degrees)
        minutes = AC.myround(degrees_frac * 60,  base=1, integer=False)
        if just_degrees_minutes:
            pass
        else:
            print('TODO - calculation not setup!')
        if debug:
            print('Conversions for: {:.7f}'.format(deg))
            # Into decimal location (the same)
            pstr = 'Decimal format for location: {:>3}{:0>5}'
            print(pstr.format(degrees, degrees_frac))
            # Into degrees, minutes and then decimal seconds
            if just_degrees_minutes:
                pstr = "Pilot format for location: {:>3}o {:0>2}'"
                print(pstr.format(degrees, minutes))
                print('')
            else:
                pstr = "Pilot format for location: {:>3}o {:0>2}.{:0>2}'"
                print(pstr.format(degrees, minutes, seconds_dec))
                print('')
        # Return in LaTeX form?
        if LaTeX_format:
            if just_degrees_minutes:
                degrees_str = str(degrees)
                minutes = '{:0>2}'.format(minutes)
                vals2rtn += [degrees_str+'$^{\circ}$'+"{}'".format(minutes)]
            else:
                print('TODO - calculation not setup!')
        else:
            if just_degrees_minutes:
                vals2rtn += ['{}{}'.forma(degrees, minutes)]
            else:
                print('TODO - calculation not setup!')
    return vals2rtn


def adjustFigAspect(fig, aspect=1):
    '''
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    '''
    xsize, ysize = fig.get_size_inches()
    minsize = min(xsize, ysize)
    xlim = .4*minsize/xsize
    ylim = .4*minsize/ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(left=.5-xlim,
                        right=.5+xlim,
                        bottom=.5-ylim,
                        top=.5+ylim)


def get_max_flying_range4BAE146(calculate_online=False):
    """
    Get the maximum flying extent of the BAE146 aircraft
    """
    if calculate_online:
        min_lat, max_lat, min_lon, max_lon = 90, -90, 180, -180
        locs2plot = 'Praia Airport', 'Dakar', 'Sao Vicente Airport',
        for loc2plot in locs2plot:
            lon, lat, alt = AC.get_loc(loc2plot)
            min_lat = min(min_lat, lat)
            max_lat = max(max_lat, lat)
            min_lon = min(min_lon, lon)
            max_lon = max(max_lon, lon)
        # What buffers to use for lat and lon
        # ~69 miles for latitude - so add buffer of 500/69 (7 degrees lat)
        lat_buffer = 7.0
        # ~55 miles for longitude - so add buffer of 9 degrees (500/55 = 9 degrees Lon)
        lon_buffer = 9.0
        # Update min and max lons
        min_lat -= lat_buffer
        max_lat += lat_buffer
        min_lon -= lon_buffer
        max_lon += lon_buffer
    else:
        # Use lats and lons.
        min_lat = 13
        max_lat = 21
        min_lon = -29
        max_lon = -21
    # Return data as a dictionary
    d = {
        'min_lat': min_lat,
        'max_lat': max_lat,
        'min_lon': min_lon,
        'max_lon': max_lon,
        'max_alt': 8000
    }
    return d


def get_dates4campaigns(year=2018, all_dates=False,
                        ground_CVAO_campaign=False):
    """
    Get a list of dates for campaigns
    """
    # Setup lists of days to use
    days2use = np.arange(1, 15)
    if isinstance(year, int):
        dates2use = [datetime.date(year, 2, i) for i in days2use]
    elif all_dates:
        dates2use = [datetime.date(2018, 2, i) for i in days2use]
        dates2use += [datetime.date(2019, 2, i) for i in days2use]
    else:
        print('WARNING: A integer year must be provided')
        sys.exit()
    # The 2019 August CVAO ground campaign (+2 days flights)
    if ground_CVAO_campaign:
        days2use = np.arange(1, 28)
        dates2use = [datetime.date(2019, 8, i) for i in days2use]
    return dates2use


def extract_ds4df_locs(ds=None, df=None, LonVar='lon', LatVar='lat',
                       TimeVar='time',
                       AltVar='hPa', dsAltVar='hPa',
                       dsLonVar='lon', dsLatVar='lat', dsTimeVar='time',
                       vars2extract=None, debug=False, testing_mode=False):
    """
    Extract a xr.Dataset as for given locations (aka offline planeflight)
    """
    # Check the DataFrame index is datetime.datetime
    ass_str = 'WARNING: DataFrame index is not of datetime.datetime dtype'
#    assert df.index.dtype == datetime.datetime, ass_str
    # Extract all of the data variables unless a specific list is provided
    if isinstance(vars2extract, type(None)):
        vars2extract = list(ds.data_vars)
    # - Create a data frame for values
    dfN = pd.DataFrame()
    # - Loop via time
    loop_via_time = False
    if loop_via_time:
        if testing_mode:
            times2use = df.index.values[:10]
        else:
            times2use = df.index.values
        for dftime in times2use:
            # get the times for a specific data
            df_tmp = df.loc[df.index == dftime, :]
            if debug:
                print(df_tmp.shape)
            lat_ = df_tmp[LatVar].values
            hPa_ = df_tmp[AltVar].values
            lon_ = df_tmp[LonVar].values
            # Extract nearest point using xarray function
            ds_tmp = ds.sel(lon=lon_, lat=lat_, lev=hPa_, time=dftime,
                            method='nearest')
            # loop and extract by data variable
            for data_var in vars2extract:
                dfN.loc[dftime, data_var] = ds_tmp[data_var].values
            del ds_tmp, df_tmp
    else:
        # get indexes =en masse then extract with these
        d = AC.calc_4D_idx_in_ds(ds=ds, df=df)
        # Run in testing mode?
        if testing_mode:
            times2use = df.index.values[:10]
        else:
            times2use = df.index.values
        # Loop by timestamp
        for n, time in enumerate(times2use):
            # get the times for a specific data
            lat_idx = d['lat'][n]
            lon_idx = d['lon'][n]
            lev_idx = d['hPa'][n]
            time_idx = d['time'][n]
            # Extract nearest point using isel function
            ds_tmp = ds.isel(lat=lat_idx, lon=lon_idx, time=time_idx,
                             lev=lev_idx)
#            d_tmp = ds_tmp.to_dict()
#            vals = ds_tmp.to_array().values
            vals = [ds_tmp[i].data for i in vars2extract]
            vals = np.array(vals)
            for nval, val in enumerate(vals):
                dfN.loc[vars2extract[nval], time] = vals[nval]
            # Add the model position coordinates...
            dfN.loc['model-lat', time] = float(ds_tmp['lat'].values)
            dfN.loc['model-lon', time] = float(ds_tmp['lon'].values)
            dfN.loc['model-lev', time] = float(ds_tmp['lev'].values)
            dfN.loc['model-time', time] = float(ds_tmp['time'].values)
            del ds_tmp, vals
        # Add the tracer names to the index
#        dfN.index = vars2extract
        # Make datetime the index
        dfN = dfN.transpose()
        # Save the datetime as a column too
        dfN['Datetime'] = dfN.index.values
        # Update the model datetime to be in datetime units
        dfN['model-time'] = pd.to_datetime(dfN['model-time'].values)
    return dfN


def summarise_stats_on_species_in_ds4lvls(ds, region='Cape_Verde', extr_str='',
                                          vars2use=None, year=2018):
    """
    Summarise the stats on species at different heights
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
            S = pd.Series(vals).describe()
            # Save to the main DataFrame
            df[varname] = S
            del vals, S
    # Save summary statistics
    savename = 'ARNA_stats_species_{}_{}_{}.csv'.format(region, year, extr_str)
    df.to_csv(savename)


def get_analysis_region(RegionName):
    """
    Function to store analysis regions extents (in degrees) for ARNA work
    """
    d = {
        # "context area"
        "context_area": {'x0': -40, 'x1': 30, 'y0': -5, 'y1': 45},
        # Download region for GEOS-CF/5
        'OPeNDAP_download_area': {'x0': -35, 'x1': 10, 'y0': 0, 'y1': 60},
        # Cape Verde Regional perspective
        # ... ?

        # Cape Verde Nested offline model region
        # (for both 0.25x0.315 and 0.5*x.625 resolutions)
        # ... ?
        # TODO: Needs checking on plots (-32 or -35 in plots , rest the same )
        'model_nest_area': {'x0': -32, 'x1': 15, 'y0': 0, 'y1': 34},
        # Cape Verde "local" perspective
        'local_CVAO_area': {'x0': -30, 'x1': -10, 'y0': 0, 'y1': 25},

        # Cape Verde flying area for ARNA
        'Cape_Verde_Flying': {'x0': -29.1, 'x1': -15.9, 'y0': 11.9,
                              'y1': 21.1},
    }
    return d[RegionName]


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


def get_local_folder(key, host=None, rtn_dict=False):
    """
    Hold folders in a dictionary and return specific variables or a dictionary
    """
    import platform
    # Get the host
    host = platform.node()
    # - Set locations of York HPC
    if ('viking' in host):
        earth0_data_shelf = '/mnt/lustre/groups/chem-acm-2018/'
        NASA_data = earth0_data_shelf + 'earth0_data/NASA/'
        HEMCO_data = earth0_data_shelf + 'earth0_data/GEOS/ExtData/HEMCO/'
        ARNA_data = '/users/ts551/scratch/data/ARNA/'
        DataRoot = '/users/ts551/scratch/data/'
        RunRoot = '/users/ts551/scratch/GC/rundirs/P_ARNA/'
    elif ('earth0' in host):
        NASA_data = '/work/data/NASA/'
        ARNA_data = ''
        DataRoot = ''
        RunRoot = ''
    elif ('Tomas' in host):
        NASA_data = '/work/data/NASA/'
        ARNA_data = '/Users/tomassherwen/tmp/ARNA/'
        DataRoot = '/Users/tomassherwen/Google_Drive/Data/'
        RunRoot = '/Users/tomassherwen/tmp/ARNA_TMP_RUNS/'
    else:
        print('NASA folder loction not known')
    # - Setup a dictionary for the variables
    d = {
        'NASA_data': NASA_data,
        #    'folder4plots': folder4plots,
        'HEMCO_data': HEMCO_data,
        'ARNA_data': ARNA_data,
        'RunRoot': RunRoot,
        'DataRoot': DataRoot,
    }
    return d[key]

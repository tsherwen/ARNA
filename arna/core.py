"""
Module of analysis and processing functions for the ARNA project
"""
import xarray as xr
import numpy as np
import pandas as pd


def main():
    """
    Main driver for analysis for the ARNA campaign
    """
    pass # This file is now treated as a module and not called directly


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


if __name__ == "__main__":
    main()

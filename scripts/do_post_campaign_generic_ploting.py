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
    # Plot up core FAAM data
#    extract_GEOS54all_ARNA_flights()
    ar.plt_timeseries_comp4ARNA_flights()
    ar.plt_comp_by_alt_4ARNA_flights()

    # plot up ToF-CIMS data
    ar.plt_comp_by_alt_4ARNA_flights_CIMS()
    ar.plt_timeseries_comp4ARNA_flights_CIMS()

    # Plot up nitrate aerosol data
    ar.plt_timeseries_comp4ARNA_flights_filters()

    #  Plot up PCASP/CDP date
    # NOTE: CAS data being ignored currently due to issue with mirror window
#    ar.plt_timeseries_comp4ARNA_flights_PCASP()

    # Plot up SWAS data
    ar.plt_timeseries_comp4ARNA_flights_SWAS()

    # Plot up vertical velocity and Roll, amongst other core physical vars
    ar.plt_timeseries_comp4ARNA_flights_PHYSICAL_VARS()

    # Plot up the temperature data from Hannah Price
    # N/A? this is only for 2019. Data to be worked up for 2020.

    # Plot a comparison
    ar.plt_timeseries_comp4ARNA_flights_NOy_ALL()

    # Plot up data for SLRs with and without dust
    ar.plt_comp_by_alt_4ARNA_all(just_SLR=False)
    ar.plt_comp_by_alt_4ARNA_all(just_SLR=True)
    ar.plt_comp_by_alt_4ARNA_all_DUST(plt_model=False)
    ar.plt_comp_by_alt_4ARNA_all_DUST(plt_model=True)
    ar.plt_comp_by_alt_4ARNA_CIMS_all_DUST()


if __name__ == "__main__":
    main()
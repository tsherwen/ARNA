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
#    extract_GEOS54all_ARNA_flights()
    ar.plt_timeseries_comparisons4ARNA_flights()
    ar.plt_alt_binned_comparisons4ARNA_flights()


if __name__ == "__main__":
    main()
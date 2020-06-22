"""
Driver script to do miscellaneous pieces of post campaign analysis
"""

#!/usr/bin/python
# import ARNA analysis/campaign code as a module
import arna as ar


def main():
    """
    Driver to do miscellaneous pieces of post campaign analysis
    """
    # Extract the extended GEOS-CF output for the FAAM flight-tracks
    ar.extract_GEOS54all_ARNA_flights()


if __name__ == "__main__":
    main()

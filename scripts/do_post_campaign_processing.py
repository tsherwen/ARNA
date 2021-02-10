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
    # Extract the expanded GEOS-CF output for the FAAM flight-tracks
#    ar.extract_GEOS54all_ARNA_flights()
    # Extract the expanded GEOS-CF output for the FAAM flight-tracks
    ar.extract_GEOS54all_ARNA_surface_dates()

    # Save model output to csv file
    ar.save_model_output2csv(RunSet='FP-MOYA-Nest', res='0.25x0.3125')
    ar.save_model_output2csv(RunSet='FP-Nest', res='0.25x0.3125')



if __name__ == "__main__":
    main()

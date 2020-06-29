"""
Driver script to do NOx/HONO (post campaign) analysis
"""
#!/usr/bin/python
# - Packages
import numpy as np
from time import gmtime, strftime
import time
import glob
import AC_tools as AC
import sys
import pandas as pd
# import ARNA analysis/campaign code as a module
import arna as ar




def main():
    """
    Main driver function
    """
    # Get the core FAAM data
#    ar.get_FAAM_core4flightnum()
    # Get the ToF-CIMS data
    ds2 = ar.get_CIMS_data4flight()




if __name__ == "__main__":
    main()
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

    # - Save model output to csv file for
    ar.save_model_output2csv(RunSet='FP-MOYA-Nest', res='0.25x0.3125')
    ar.save_model_output2csv(RunSet='FP-Nest', res='0.25x0.3125')

    # - Save out model output (j-rate, OH) for surface and airborne campaign
    campaigns = ['ARNA-1', 'ARNA-1-surface', 'CVAO-2015-surface']
    dfs = {}
    for campaign in campaigns:
        dfs[campaign] = ar.get_whole_related_campaign_data(campaign=campaign,
                                                           save2csv=True)

    # plot up test plots for campaign data (e.g. lat & lon vs time)
    for campaign in campaigns:
        df = dfs[campaign].copy()
        try:
            TYPE = 'CVO1'
            df = df.loc[df['TYPE'] == TYPE, :]
        except:
            pass

        savetitle = 'timeseries_plot_{}'.format(campaign)
        plt_quick_ts4df(df, savetitle=savetitle)

    # Generic checking on boundary runs for ARNA


if __name__ == "__main__":
    main()

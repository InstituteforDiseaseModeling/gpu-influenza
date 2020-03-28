""" gpuflu

"""
# Standard packages
import numpy
import pandas
import os
import sys
import argparse

# IDM packages
import gpuflu




def getParametersFromConfig( config ):
    """Read configuration file and return separate dataframes for 
       simulation parameters and for report parameters.
    Args:
     config : Name of configuration file or "demo" for executing
              the simulator demo. The configuration file must have
              ".conf" extension.
    Return:
     simulation and report parameter dataframes
    """ 

    # Read configuration file name
    configpath = os.path.splitext( config )
    configFilename  = configpath[0]
    configExtension = configpath[1]

    # Initialize simulation and report dataframes
    simulation_params = pandas.DataFrame()
    report_params = pandas.DataFrame()

    # Set parameters for demo
    if configFilename == "demo":
        nSims = 50
        simulation_paramsA = pandas.DataFrame( { "population"  : [ 10000 ],
                                                "connections" : [    50 ],
                                                "exposures"   : [     0 ],
                                                "infections"  : [     1 ],
                                                "time"        : [   100 ],
                                                "durationE"   : [     4 ],
                                                "stddevE"     : [     1 ],
                                                "durationI"   : [     5 ],
                                                "stddevI"     : [     2 ],
                                                "r0"          : [   2.3 ],
                                                "stddevr0"    : [   1.0 ]
                                               }
                                             )
        simulation_paramsB = pandas.DataFrame( { "population"  : [ 10000 ],
                                                "connections" : [    50 ],
                                                "exposures"   : [     0 ],
                                                "infections"  : [     1 ],
                                                "time"        : [   100 ],
                                                "durationE"   : [     4 ],
                                                "stddevE"     : [     1 ],
                                                "durationI"   : [     5 ],
                                                "stddevI"     : [     2 ],
                                                "r0"          : [   1.5 ],
                                                "stddevr0"    : [   1.0 ]
                                               }
                                             )
        simulation_params \
                   = simulation_params.append( [simulation_paramsA]*(nSims), 
                                               ignore_index=True 
                                              )
        simulation_params \
                   = simulation_params.append( [simulation_paramsB]*(nSims), 
                                               ignore_index=True 
                                              )
        numpy.random.seed(101)
        simulation_params["seed"] = numpy.random.randint( 0, 
                                                          numpy.iinfo(numpy.int32).max,
                                                          len(simulation_params),
                                                          dtype="int"
                                                         )

        report_params = pandas.DataFrame()
        report_params.loc[0, "type"      ] = "timeseries"
        report_params.loc[0, "subtype"   ] = "simple"
        report_params.loc[0, "title"     ] = "Infectious"
        report_params.loc[0, "label"     ] = "Infectious"
        report_params.loc[0, "xlabel"    ] = "day"
        report_params.loc[0, "xtype"     ] = "range"
        report_params.loc[0, "data"      ] = "I"
        report_params.loc[0, "sim_start" ] = 0
        report_params.loc[0, "sim_end"   ] = nSims
        report_params.loc[0, "report_id" ] = 0
        report_params.loc[0, "filename"  ] = "simA_infectious.png"
                                            
        report_params.loc[1, "type"      ] = "timeseries"
        report_params.loc[1, "subtype"   ] = "median"
        report_params.loc[1, "qlow"      ] = 0.1
        report_params.loc[1, "qhigh"     ] = 0.9
        report_params.loc[1, "title"     ] = "Infectious"
        report_params.loc[1, "label"     ] = "Infectious"
        report_params.loc[1, "xlabel"    ] = "day"
        report_params.loc[1, "xtype"     ] = "range"
        report_params.loc[1, "data"      ] = "I"
        report_params.loc[1, "sim_start" ] = 0
        report_params.loc[1, "sim_end"   ] = nSims
        report_params.loc[1, "report_id" ] = 1
        report_params.loc[1, "filename"  ] = "simA_infectious_median.png"

        report_params.loc[2, "type"      ] = "timeseries"
        report_params.loc[2, "subtype"   ] = "median"
        report_params.loc[2, "title"     ] = "Counts"
        report_params.loc[2, "label"     ] = "Exposed"
        report_params.loc[2, "xlabel"    ] = "day"
        report_params.loc[2, "xtype"     ] = "range"
        report_params.loc[2, "data"      ] = "E"
        report_params.loc[2, "sim_start" ] = 0
        report_params.loc[2, "sim_end"   ] = nSims
        report_params.loc[2, "report_id" ] = 2
        report_params.loc[2, "filename"  ] = "simA_EIR.png"

        report_params.loc[3, "type"      ] = "timeseries"
        report_params.loc[3, "subtype"   ] = "median"
        report_params.loc[3, "label"     ] = "Infectious"
        report_params.loc[3, "xlabel"    ] = "day"
        report_params.loc[3, "xtype"     ] = "range"
        report_params.loc[3, "data"      ] = "I"
        report_params.loc[3, "sim_start" ] = 0
        report_params.loc[3, "sim_end"   ] = nSims
        report_params.loc[3, "report_id" ] = 2

        report_params.loc[4, "type"      ] = "timeseries"
        report_params.loc[4, "subtype"   ] = "median"
        report_params.loc[4, "label"     ] = "Removed"
        report_params.loc[4, "xlabel"    ] = "day"
        report_params.loc[4, "xtype"     ] = "range"
        report_params.loc[4, "data"      ] = "R"
        report_params.loc[4, "sim_start" ] = 0
        report_params.loc[4, "sim_end"   ] = nSims
        report_params.loc[4, "report_id" ] = 2

        report_params.loc[5, "type"      ] = "timeseries"
        report_params.loc[5, "subtype"   ] = "cumsum_median"
        report_params.loc[5, "title"     ] = "Infectious (cumulative)"
        report_params.loc[5, "label"     ] = "Infectious"
        report_params.loc[5, "xlabel"    ] = "day"
        report_params.loc[5, "xtype"     ] = "range"
        report_params.loc[5, "data"      ] = "I"
        report_params.loc[5, "sim_start" ] = 0
        report_params.loc[5, "sim_end"   ] = nSims
        report_params.loc[5, "report_id" ] = 5
        report_params.loc[5, "filename"  ] = "simA_infectious_cumulative.png"

        report_params.loc[6, "type"      ] = "timeseries"
        report_params.loc[6, "subtype"   ] = "median"
        report_params.loc[6, "qlow"      ] = 0.1
        report_params.loc[6, "qhigh"     ] = 0.9
        report_params.loc[6, "title"     ] = "Infectious"
        report_params.loc[6, "label"     ] = "Ro_A"
        report_params.loc[6, "xlabel"    ] = "day"
        report_params.loc[6, "xtype"     ] = "range"
        report_params.loc[6, "data"      ] = "I"
        report_params.loc[6, "sim_start" ] = 0
        report_params.loc[6, "sim_end"   ] = nSims
        report_params.loc[6, "report_id" ] = 6
        report_params.loc[6, "filename"  ] = "simB_infectious_median.png"

        report_params.loc[7, "type"      ] = "timeseries"
        report_params.loc[7, "subtype"   ] = "median"
        report_params.loc[7, "qlow"      ] = 0.1
        report_params.loc[7, "qhigh"     ] = 0.9
        report_params.loc[7, "title"     ] = "Infectious"
        report_params.loc[7, "label"     ] = "Ro_B"
        report_params.loc[7, "xlabel"    ] = "day"
        report_params.loc[7, "xtype"     ] = "range"
        report_params.loc[7, "data"      ] = "I"
        report_params.loc[7, "sim_start" ] = nSims
        report_params.loc[7, "sim_end"   ] = 2*nSims-1
        report_params.loc[7, "report_id" ] = 6

        report_params.loc[8, "type"      ] = "csv"
        report_params.loc[8, "data"      ] = "I"
        report_params.loc[8, "sim_start" ] = 0
        report_params.loc[8, "sim_end"   ] = nSims
        report_params.loc[8, "report_id" ] = 8
        report_params.loc[8, "filename"  ] = "sim_infectious.csv"

        report_params.loc[9, "type"      ] = "csv"
        report_params.loc[9, "data"      ] = "I"
        report_params.loc[9, "sim_start" ] = nSims
        report_params.loc[9, "sim_end"   ] = 2*nSims-1
        report_params.loc[9, "report_id" ] = 8

        report_params.loc[10, "type"     ] = "histogram"
        report_params.loc[10, "subtype"  ] = "simple"
        report_params.loc[10, "bins"     ] = 10
        report_params.loc[10, "title"    ] = "Infectious (frequency)"
        report_params.loc[10, "ylabel"   ] = "day 10"
        report_params.loc[10, "xlabel"   ] = "Number of infectious individuals"
        report_params.loc[10, "col"      ] = "I_10"
        report_params.loc[10, "sim_start"] = 0
        report_params.loc[10, "sim_end"  ] = nSims
        report_params.loc[10, "report_id"] = 10
        report_params.loc[10, "filename" ] = "simA_infectious_hist.png"

        report_params.loc[11, "type"     ] = "histogram"
        report_params.loc[11, "subtype"  ] = "simple"
        report_params.loc[11, "bins"     ] = 20
        report_params.loc[11, "ylabel"   ] = "day 20"
        report_params.loc[11, "col"      ] = "I_20"
        report_params.loc[11, "sim_start"] = 0
        report_params.loc[11, "sim_end"  ] = nSims
        report_params.loc[11, "report_id"] = 10

        report_params.loc[12, "type"     ] = "histogram"
        report_params.loc[12, "subtype"  ] = "simple"
        report_params.loc[12, "bins"     ] = 50
        report_params.loc[12, "ylabel"   ] = "day 30"
        report_params.loc[12, "col"      ] = "I_30"
        report_params.loc[12, "sim_start"] = 0
        report_params.loc[12, "sim_end"  ] = nSims
        report_params.loc[12, "report_id"] = 10

        report_params.loc[13, "type"     ] = "timeseries"
        report_params.loc[13, "subtype"  ] = "simple"
        report_params.loc[13, "title"    ] = "Incidence"
        report_params.loc[13, "label"    ] = "Incidence"
        report_params.loc[13, "xlabel"   ] = "day"
        report_params.loc[13, "xtype"    ] = "range"
        report_params.loc[13, "data"     ] = "incidence"
        report_params.loc[13, "sim_start"] = 0
        report_params.loc[13, "sim_end"  ] = nSims
        report_params.loc[13, "report_id"] = 13
        report_params.loc[13, "filename" ] = "simA_incidence.png"

        report_params.loc[14, "type"      ] = "timeseries"
        report_params.loc[14, "subtype"   ] = "cumsum_median"
        report_params.loc[14, "title"     ] = "Incidence (cumulative)"
        report_params.loc[14, "label"     ] = "Incidence (cumulative)"
        report_params.loc[14, "xlabel"    ] = "day"
        report_params.loc[14, "xtype"     ] = "range"
        report_params.loc[14, "data"      ] = "incidence"
        report_params.loc[14, "sim_start" ] = 0
        report_params.loc[14, "sim_end"   ] = nSims
        report_params.loc[14, "report_id" ] = 14
        report_params.loc[14, "filename"  ] = "simA_incidence_cumulative.png"

        report_params.loc[15, "type"      ] = "timeseries"
        report_params.loc[15, "subtype"   ] = "median"
        report_params.loc[15, "qlow"      ] = 0.1
        report_params.loc[15, "qhigh"     ] = 0.9
        report_params.loc[15, "title"     ] = "Incidence"
        report_params.loc[15, "label"     ] = "Ro_A"
        report_params.loc[15, "xlabel"    ] = "day"
        report_params.loc[15, "xtype"     ] = "range"
        report_params.loc[15, "data"      ] = "incidence"
        report_params.loc[15, "sim_start" ] = 0
        report_params.loc[15, "sim_end"   ] = nSims
        report_params.loc[15, "report_id" ] = 15
        report_params.loc[15, "filename"  ] = "simB_incidence_median.png"

        report_params.loc[16, "type"      ] = "timeseries"
        report_params.loc[16, "subtype"   ] = "median"
        report_params.loc[16, "qlow"      ] = 0.1
        report_params.loc[16, "qhigh"     ] = 0.9
        report_params.loc[16, "title"     ] = "Incidence"
        report_params.loc[16, "label"     ] = "Ro_B"
        report_params.loc[16, "xlabel"    ] = "day"
        report_params.loc[16, "xtype"     ] = "range"
        report_params.loc[16, "data"      ] = "incidence"
        report_params.loc[16, "sim_start" ] = nSims
        report_params.loc[16, "sim_end"   ] = 2*nSims-1
        report_params.loc[16, "report_id" ] = 15


    # Set parameters based on configuration file
    elif configExtension == ".json":

        # Column names for simulations and report (they are used for cleaning
        # simulations and reports dataframes)
        reportParamNames = [ "def",
                             "sim_id",
                             "nsims",
                             "data",
                             "filename",
                             "label",
                             "qhigh",
                             "qlow",
                             "report_id",
                             "sim_end",
                             "sim_start",
                             "subtype",
                             "title",
                             "type",
                             "xlabel",
                             "xtype",
                             "ylabel"
                            ]
        simulationParamNames = [ "def", 
                                 "sim_id", 
                                 "nsims", 
                                 "population",
                                 "connections",
                                 "exposures",
                                 "infections",
                                 "time",
                                 "durationE",
                                 "stddevE",
                                 "durationI",
                                 "stddevI",
                                 "r0",
                                 "stddevr0",
                                 "seed"
                                ]


        # Read configuration file
        configParams = pandas.read_json( config, orient="records" )


        # Get simulation parameters
        sims = configParams[configParams["def"]=="simulation"].groupby("sim_id")
        sim_start = int(0)
        simulationConfigInfo = pandas.DataFrame() 
        for simGroupId in sims.groups:
            # Read and clean simulation parameters
            thisSim = sims.get_group( simGroupId )
            simId = thisSim["sim_id"].values[0]
            nSims = int( thisSim["nsims"].values[0] )
            #thisSimClean = thisSim.drop(columns=["def", "sim_id", "nsims"]).copy()
            thisSimClean = thisSim.drop( columns = reportParamNames, 
                                         errors  = "ignore"
                                        ).copy()

            # Add this simulation params to overall simulation params dataframe
            simulation_params \
                   = simulation_params.append( [thisSimClean]*nSims, 
                                               ignore_index=True 
                                              )

            # Add seed if necessary
            if "seed" in thisSim.columns:
                numpy.random.seed( int( thisSim["seed"].values[0] ) )
            simulation_params.loc[sim_start:sim_start+nSims,"seed"] \
                                 = numpy.random.randint( 0, 
                                                         numpy.iinfo(numpy.int32).max,
                                                         nSims,
                                                         dtype="int"
                                                        )

            # Prepare for next set of simulation parameters
            simulationConfigInfo.loc[simId, "start"] = int( sim_start )
            simulationConfigInfo.loc[simId, "end"  ] = int( sim_start + nSims )
            sim_start = sim_start + nSims

        # Get reports parameters
        reports = configParams[configParams["def"]=="report"]
        for reportIndex, report in reports.iterrows():
            # Read report parameters
            simId = report["sim_id"]
            simStart = int( simulationConfigInfo.loc[simId, "start"] )
            simEnd   = int( simulationConfigInfo.loc[simId, "end"  ] )
            report["sim_start"] = simStart
            report["sim_end"  ] = simEnd

            # Clean and add this report params to overall set of params
            thisReportClean = report.drop( labels = simulationParamNames,
                                           errors = "ignore"
                                          )
            report_params = report_params.append( thisReportClean )


    # Configuration file not supported
    else:
        sys.exit( "configuration file format not supported" )

    #print("****simulation_params = \n", simulation_params)
    #print("****report_params = \n", report_params)
    #quit()
    # Finalize and return
    return simulation_params, report_params




#-----------------------------------------------------------------------------
def buildReports( report_params, EIR ):
    """Process simulation results and generate output reports.
    Args:
     report_params : Pandas dataframe containing report generation parameters.
     EIR           : Simulation results (pandas dataframe)
    """
    reports = report_params.groupby( ["report_id"] ) 

    for reportId in reports.groups:
        report = reports.get_group( reportId )
        reportType = report["type"].values[0]

        if reportType == "timeseries":
            gpuflu.reports.timeseries( report, EIR )
        elif reportType == "histogram":
            gpuflu.reports.histogram( report, EIR )
        elif reportType == "csv":
            gpuflu.reports.csv( report, EIR )
        else:
            print( "report type ", reportType, " not supported" )

    return




#----------------------------------------------------------------------------
#    MAIN
#-----------------------------------------------------------------------------
if __name__ == "__main__":

    # Read arguments
    parser = argparse.ArgumentParser( prog        ="gpuflu",
                                      usage       ="python gpuflu.py [options]",
                                      description ="Run gpu-influenza" 
                                     )
    parser.add_argument( "-c", 
                         "--config", 
                         dest = "configfile",
                         type = str,
                         help = "configuration file",
                         required=True 
                        )
    args = parser.parse_args()

    simulationParameters, \
    reportParameters      = getParametersFromConfig( args.configfile )


    # Run simulation
    EIR = gpuflu.interface.gpu_influenza( simulationParameters )


    # Prepare output and exit
    buildReports( reportParameters, EIR )
    sys.exit()


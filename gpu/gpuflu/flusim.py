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


    # Set parameters for demo
    if configFilename == "demo":
        nSims = 20
        simulation_params = pandas.DataFrame( { "population"  : [ 10000 ],
                                                "connections" : [   500 ],
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
                                                "connections" : [   500 ],
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
                   = simulation_params.append( [simulation_params]*(nSims-1), 
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
        report_params.loc[0, "reportid"  ] = 0
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
        report_params.loc[1, "reportid"  ] = 1
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
        report_params.loc[2, "reportid"  ] = 2
        report_params.loc[2, "filename"  ] = "simA_EIR.png"

        report_params.loc[3, "type"      ] = "timeseries"
        report_params.loc[3, "subtype"   ] = "median"
        report_params.loc[3, "label"     ] = "Infectious"
        report_params.loc[3, "xlabel"    ] = "day"
        report_params.loc[3, "xtype"     ] = "range"
        report_params.loc[3, "data"      ] = "I"
        report_params.loc[3, "sim_start" ] = 0
        report_params.loc[3, "sim_end"   ] = nSims
        report_params.loc[3, "reportid"  ] = 2

        report_params.loc[4, "type"      ] = "timeseries"
        report_params.loc[4, "subtype"   ] = "median"
        report_params.loc[4, "label"     ] = "Removed"
        report_params.loc[4, "xlabel"    ] = "day"
        report_params.loc[4, "xtype"     ] = "range"
        report_params.loc[4, "data"      ] = "R"
        report_params.loc[4, "sim_start" ] = 0
        report_params.loc[4, "sim_end"   ] = nSims
        report_params.loc[4, "reportid"  ] = 2

        report_params.loc[5, "type"      ] = "timeseries"
        report_params.loc[5, "subtype"   ] = "cumsum_median"
        report_params.loc[5, "title"     ] = "Infectious (cumulative)"
        report_params.loc[5, "label"     ] = "Infectious"
        report_params.loc[5, "xlabel"    ] = "day"
        report_params.loc[5, "xtype"     ] = "range"
        report_params.loc[5, "data"      ] = "I"
        report_params.loc[5, "sim_start" ] = 0
        report_params.loc[5, "sim_end"   ] = nSims
        report_params.loc[5, "reportid"  ] = 5
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
        report_params.loc[6, "reportid"  ] = 6
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
        report_params.loc[7, "reportid"  ] = 6

        report_params.loc[8, "type"      ] = "csv"
        report_params.loc[8, "data"      ] = "I"
        report_params.loc[8, "sim_start" ] = 0
        report_params.loc[8, "sim_end"   ] = nSims
        report_params.loc[8, "reportid"  ] = 8
        report_params.loc[8, "filename"  ] = "sim_infectious.csv"

        report_params.loc[9, "type"      ] = "csv"
        report_params.loc[9, "data"      ] = "I"
        report_params.loc[9, "sim_start" ] = nSims
        report_params.loc[9, "sim_end"   ] = 2*nSims-1
        report_params.loc[9, "reportid"  ] = 8

    # Set parameters based on configuration file
    elif configExtension == ".config":
        sys.exit( "configuration file execution not supported" )


    # Configuration file not supported
    else:
        sys.exit( "configuration file format not supported" )


    # Finalize and return
    return simulation_params, report_params




#-----------------------------------------------------------------------------
def buildReports( report_params, EIR ):
    """Process simulation results and generate output reports.
    Args:
     report_params : Pandas dataframe containing report generation parameters.
     EIR           : Simulation results (pandas dataframe)
    """
    reports = report_params.groupby( ["reportid"] ) 

    for reportId in reports.groups:
        report = reports.get_group( reportId )
        reportType = report["type"].values[0]

        if reportType == "timeseries":
            gpuflu.reports.timeseries( report, EIR )
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


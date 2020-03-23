import numpy
import pandas
import datetime
import re
import matplotlib
from matplotlib import pyplot as plt
matplotlib.pyplot.switch_backend("Agg")




def timeseries( params, data ):
    """Generate timeseries figures for gpu-influenza simulations.
     Args:
      params : Pandas dataframe with parameters for the reports. Each row
               defines a report and each columns defines a parameter for
               reports. The following columns are supported:
               - subtype   : Defines the type of timeseries figure to 
                             generate. It can be:
                             "cumsum_median" : Cumulative sum of the median
                                               of selected data.
                             "median"        : Median of selected data.
                             "simple"        : (default) Plot selected data.
               - data      : Regular expression identifying the set of 
                             columns to be included in the data selection.
               - title     : Figure title.
               - label     : Data series label.
               - xlabel    : x-axis label.
               - xtype     : Type of ticks to be used in the x axis. If not
                             define use the name of the columns of the 
                             selected data. Otherwise, it can be:
                             "range" : Ordinal names are used in the x-axis.
               - ylabel    : y-axis label.
               - sim_start : Index of first simulation in include in the
                             data selection.
               - sim_end   : Index +1 of last simulation to include in the
                             data selection.
               - qlow      : Lower quantile. Only supported for "median"
                             subtype.
               - qhigh     : Upper quantile. Only supported for "median" 
                             subtype.
               - filename  : Name of output file.  
      data   : Simulation results (pandas dataframe).
    """

    # Read figure formatting parameters
    subtype  = params["subtype" ].values[0] if "subtype" in params else False
    title    = params["title"   ].values[0] if "title"   in params else False
    xlabel   = params["xlabel"  ].values[0] if "xlabel"  in params else False
    xtype    = params["xtype"   ].values[0] if "xtype"   in params else False
    ylabel   = params["ylabel"  ].values[0] if "ylabel"  in params else False
    filename = params["filename"].values[0]


    # Read data processing parameters
    qlow  = params["qlow" ].values[0] if "qlow"  in params else False
    qhigh = params["qhigh"].values[0] if "qhigh" in params else False


    # Figure properties
    propCycle = plt.rcParams["axes.prop_cycle"]
    colors = propCycle.by_key()["color"]
    iColor = 0
    nColors = len(colors)
    nx = 0


    # Draw figure
    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    for index, row in params.iterrows():

        xPrefix = row["data"     ] 
        start   = int( row["sim_start"] )
        end     = int( row["sim_end"  ] )
        label   = row["label"    ]

        x = data.filter( regex = "^" + xPrefix, axis = 1 ).iloc[start:end]

        lineColor = colors[iColor]
        iColor = iColor + 1


        if subtype == "cumsum_median":  #----------------------------- cumsum

            if xtype == "range":
                x.cumsum(axis=1)    \
                 .median(axis=0).transpose().plot( use_index = False,
                                                   label     = label,
                                                   color     = lineColor,
                                                   legend    = True,
                                                   ax        = ax
                                                  )
            else:
                x.cumsum(axis=1)    \
                 .median(axis=0).transpose().plot( label  = label,
                                                   color  = lineColor,
                                                   legend = False, 
                                                   ax     = ax
                                                  )

        elif subtype == "median": #----------------------------------- median

            if xtype == "range":  
                x.median(axis=0).transpose().plot( use_index = False,
                                                   label     = label,
                                                   color     = lineColor,
                                                   legend    = True,
                                                   ax        = ax
                                                  )
            else:
                x.median(axis=0).transpose().plot( label  = label,
                                                   color  = lineColor,
                                                   legend = False, 
                                                   ax     = ax
                                                  )

            if (not numpy.isnan(qlow)) and (not numpy.isnan(qhigh)):
                xlow  = x.quantile(q=qlow )
                xhigh = x.quantile(q=qhigh)

                xlow .plot( use_index = False,
                            label     = label + "(q=" + str(qlow) + ")",
                            legend    = True,
                            color     = lineColor,
                            ls        = "--",
                            lw        = 1,
                            ax        = ax
                           )
                xhigh.plot( use_index = False,
                            label     = label + "(q=" + str(qhigh) + ")",
                            legend    = True,
                            color     = lineColor,
                            ls        = ":",
                            lw        = 1,
                            ax        = ax
                           )
                ax.fill_between( range(0, len(xlow)), 
                                 xlow, 
                                 xhigh, 
                                 color = lineColor,
                                 alpha = 0.15 )

        else:    #--------------------------------------------------- simple
            if xtype == "range":
                x.transpose().plot( use_index = False,
                                    label     = label,
                                    legend    = False,
                                    ax        = ax
                                   )
            else:
                x.transpose().plot( label  = label,
                                    legend = False, 
                                    ax     = ax
                                   )


    # Figure formatting
    ax.set_title( title )
    if xlabel:
        ax.set_xlabel( xlabel )
    if ylabel:
        ax.set_ylabel( ylabel )


    # Finalize and return
    plt.savefig( filename )
    plt.close(fig)
    return


import numpy
import pandas
import datetime
import re
import matplotlib
from matplotlib import pyplot as plt
matplotlib.pyplot.switch_backend("Agg")

import gpuflu
from gpuflu.interface import gpu_influenza_sampleI as gpu_influenza_sampleI




#---- Parameters ---------------------------------------------------------------

# Simulations
nDays    = 200
nSims    = 11
baseDate = "Jan 15 2020"
seed     = 101


# Scenarios
scenario_1__r0  = 2
scenario_1__tag = "Ro = " + str(scenario_1__r0) + "\n(low)"

scenario_2__r0  = 3
scenario_2__tag = "Ro = " + str(scenario_2__r0) + "\n(medium)"

scenario_3__r0  = 4
scenario_3__tag = "Ro = " + str(scenario_3__r0) + "\n(high)"

scenario_4__r0  = 5
scenario_4__tag = "Ro = " + str(scenario_4__r0) + "\n(very high)"


scenario_1 = pandas.DataFrame( { "population"  : [        2200000 ],
                                 "connections" : [              4 ],
                                 "exposures"   : [              0 ],
                                 "infections"  : [              1 ],
                                 "time"        : [          nDays ],
                                 "durationE"   : [              5 ],
                                 "stddevE"     : [            1.5 ],
                                 "durationI"   : [              5 ],
                                 "stddevI"     : [            1.5 ],
                                 "r0"          : [ scenario_1__r0 ],
                                 "stddevr0"    : [            1.5 ],
                                 }
                                )

scenario_2 = pandas.DataFrame( { "population"  : [        2200000 ],
                                 "connections" : [              4 ],
                                 "exposures"   : [              0 ],
                                 "infections"  : [              1 ],
                                 "time"        : [          nDays ],
                                 "durationE"   : [              5 ],
                                 "stddevE"     : [            1.5 ],
                                 "durationI"   : [              5 ],
                                 "stddevI"     : [            1.5 ],
                                 "r0"          : [ scenario_2__r0 ],
                                 "stddevr0"    : [            1.5 ],
                                 }
                                )

scenario_3 = pandas.DataFrame( { "population"  : [        2200000 ],
                                 "connections" : [              4 ],
                                 "exposures"   : [              0 ],
                                 "infections"  : [              1 ],
                                 "time"        : [          nDays ],
                                 "durationE"   : [              5 ],
                                 "stddevE"     : [            1.5 ],
                                 "durationI"   : [              5 ],
                                 "stddevI"     : [            1.5 ],
                                 "r0"          : [ scenario_3__r0 ],
                                 "stddevr0"    : [            1.5 ],
                                 }
                                )

scenario_4 = pandas.DataFrame( { "population"  : [        2200000 ],
                                 "connections" : [              4 ],
                                 "exposures"   : [              0 ],
                                 "infections"  : [              1 ],
                                 "time"        : [          nDays ],
                                 "durationE"   : [              5 ],
                                 "stddevE"     : [            1.5 ],
                                 "durationI"   : [              5 ],
                                 "stddevI"     : [            1.5 ],
                                 "r0"          : [ scenario_4__r0 ],
                                 "stddevr0"    : [            1.5 ],
                                 }
                                )




#---- Run simulations ----------------------------------------------------------
numpy.random.seed( seed )

scenario_1_rep = scenario_1.append( [scenario_1]*(nSims-1), ignore_index=True )
scenario_2_rep = scenario_2.append( [scenario_2]*(nSims-1), ignore_index=True )
scenario_3_rep = scenario_3.append( [scenario_3]*(nSims-1), ignore_index=True )
scenario_4_rep = scenario_4.append( [scenario_4]*(nSims-1), ignore_index=True )

y_1 = gpu_influenza_sampleI( scenario_1_rep )
y_2 = gpu_influenza_sampleI( scenario_2_rep )
y_3 = gpu_influenza_sampleI( scenario_3_rep )
y_4 = gpu_influenza_sampleI( scenario_4_rep )

# Format columns
base = datetime.datetime.strptime( baseDate, '%b %d %Y')
date_list = [base + datetime.timedelta(days=x) for x in range(nDays)]

y_1.columns = date_list
y_2.columns = date_list
y_3.columns = date_list
y_4.columns = date_list

# Compute statistics
y_1_median = y_1.median()
y_2_median = y_2.median()
y_3_median = y_3.median()
y_4_median = y_4.median()




#---- Write output -------------------------------------------------------------

# Draw plots
plt.figure()
ax1 = plt.subplot(4,1,1)
y_1.transpose().plot( title  = "cases per day",
                      color  = "gray",
                      legend = False,
                      ax     = ax1
                     )
y_1_median.transpose().plot( color  = "brown",
                             legend = False,
                             ax     = ax1
                            )
ax1.set_ylabel( scenario_1__tag )


ax2 = plt.subplot(4,1,2, sharex=ax1)
y_2.transpose().plot( color  = "gray",
                      legend = False,
                      ax     = ax2
                     )
y_2_median.transpose().plot( color  = "brown",
                             legend = False,
                             ax     = ax2
                            )
ax2.set_ylabel( scenario_2__tag )


ax3 = plt.subplot(4,1,3, sharex=ax1)
y_3.transpose().plot( color  = "gray",
                      legend = False,
                      ax     = ax3
                     )
y_3_median.transpose().plot( color  = "brown",
                             legend = False,
                             ax     = ax3
                            )
ax3.set_ylabel( scenario_3__tag )


ax4 = plt.subplot(4,1,4, sharex=ax1)
y_4.transpose().plot( color  = "gray",
                      legend = False,
                      ax     = ax4
                     )
y_4_median.transpose().plot( color  = "brown",
                             legend = False,
                             ax     = ax4
                            )
ax4.set_ylabel( scenario_4__tag )

plt.savefig( "gpu_influenza_y_allSims.png" )
plt.close()


plt.figure()
ax = plt.subplot(1,1,1)
y_1_median.transpose().plot( title  = "cases per day",
                             label  = scenario_1__tag,
                             legend = True,
                             ax     = ax
                             )
y_2_median.transpose().plot( label  = scenario_2__tag,
                             legend = True, 
                             ax     = ax 
                            )
y_3_median.transpose().plot( label  = scenario_3__tag,
                             legend = True, 
                             ax     = ax 
                            )
y_4_median.transpose().plot( label  = scenario_4__tag,
                             legend = True, 
                             ax     = ax 
                            )

plt.savefig( "gpu_influenza_y_medianSims.png" )
plt.close()


# Write output files
scenario_1__y_outfile       = "gpu_influenza_sims_"                          \
                              + re.sub( "[^a-zA-Z\d]", "_", scenario_1__tag) \
                              + ".csv"
scenario_1__yMedian_outfile = "gpu_influenza_sims_median_"                   \
                              + re.sub( "[^a-zA-Z\d]", "_", scenario_1__tag) \
                              + ".csv"

scenario_2__y_outfile       = "gpu_influenza_sims_"                          \
                              + re.sub( "[^a-zA-Z\d]", "_", scenario_2__tag) \
                              + ".csv"
scenario_2__yMedian_outfile = "gpu_influenza_sims_median_"                   \
                              + re.sub( "[^a-zA-Z\d]", "_", scenario_2__tag) \
                              + ".csv"

scenario_3__y_outfile       = "gpu_influenza_sims_"                          \
                              + re.sub( "[^a-zA-Z\d]", "_", scenario_3__tag) \
                              + ".csv"
scenario_3__yMedian_outfile = "gpu_influenza_sims_median_"                   \
                              + re.sub( "[^a-zA-Z\d]", "_", scenario_3__tag) \
                              + ".csv"

scenario_4__y_outfile       = "gpu_influenza_sims_"                          \
                              + re.sub( "[^a-zA-Z\d]", "_", scenario_4__tag) \
                              + ".csv"
scenario_4__yMedian_outfile = "gpu_influenza_sims_median_"                   \
                              + re.sub( "[^a-zA-Z\d]", "_", scenario_4__tag) \
                              + ".csv"

y_1.to_csv( scenario_1__y_outfile )
y_2.to_csv( scenario_2__y_outfile )
y_3.to_csv( scenario_3__y_outfile )
y_4.to_csv( scenario_4__y_outfile )

y_1_median.to_csv( scenario_1__yMedian_outfile )
y_2_median.to_csv( scenario_2__yMedian_outfile )
y_3_median.to_csv( scenario_3__yMedian_outfile )
y_4_median.to_csv( scenario_4__yMedian_outfile )

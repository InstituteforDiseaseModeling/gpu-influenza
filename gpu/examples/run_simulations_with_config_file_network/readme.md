This example uses synthetic population data from the following source:

2010 RTI U.S. Synthetic Population Ver. 1.0
RTI International. May, 2014. Downloaded from internet
URL: https://www.epimodels.org/midas/pubsyntdata1.do  

Instructions:

- Download data for the county of interest from 
  https://fred.publichealth.pitt.edu/syn_pops

- Update `generate_network.py` with input file name

- Run `generate_network.py`

- Run `python3 ../../gpuflu/flusim.py -c multiple_r0_sim.json`

# gpu-influenza

Network-based compartmental model. In this model, nodes represent individuals of 
a given population, and edges represent physical contacts between individuals. 
The use of a network to model these contacts allows us to incorporate and 
analyze the social structure of the transmission. `gpu-influenza` is 
implemented in CUDA to support computationally efficient models for populations 
with hundreds of millions of nodes (e.g., >200,000 million nodes on NVIDIA Tesla 
P100 GPU). It also provides a Python interface that integrates relevant features
such as network generation, analysis, and testing.

In `gpu-influenza`, disease transmission is governed by a stochastic SEIR 
model. Between-node transmission is modeled as follows:

- **Exposure**: Nodes become exposed either by acquired transmission from a 
  neighboring node or as the result of an importation event. Exposed nodes 
  become infectious after a (lognormal) exposure duration. 
  
- **Infection**: Nodes become infected after they leave their exposure state. 
  Infectious nodes recover after a lognormal infectious duration. At each 
  timestep, each infectious node k transmits the infection to neighboring nodes 
  with probability determined by Ro and intervention-based attenuation factors.
  
- **Recovery/removed**: Nodes enter the recovered/removed compartment after they 
  complete their infection duration.
  
- **Susceptibility**: All nodes that are not exposed, infected, or recovered, are 
  susceptible. They are removed from the Susceptible compartment either by 
  importation or infection events.
  
- **Importations**: At each time step, new importations are added to the exposed 
  compartment following a Poisson distribution. Specific nodes that become 
  exposed due to importation events are chosen randomly from the pool of 
  susceptible nodes.
  
- **Initial conditions**: A set of initial exposed and infectious nodes are randomly 
  selected based on parameters Eo and Io, respectively.




## Installation


1. Download the library
 
     git clone https://github.com/InstituteforDiseaseModeling/gpu-influenza


2. Install dependencies

   Packages and libraries not on PyPI will not be installed by setuptools in 
   step 3 below. Hence, they should be installed separately. These packages 
   include CUDA and its dependencies. 

   A Conda environment definition file is provided with this package. This
   file, called `gpu_influenza_env.yml`. To use it for installing 
   dependencies, edit the  first line to set the name of the environment, and 
   then execute the following commands for installing and activating a 
   conda environment that satisfies all required dependencies:

     `conda env create -f gpu_influenza_env.yml`

     `conda activate [environment name]`

   NOTE: Environment activation command is `source activate [environment name]`
         instead of `conda activate [environment name]` in some systems. 


3. Install package

    `phyton3 setup.py develop`

   NOTES:
 
   a. Alternatively, you could run `python3 setup.py install` if you plan
      on using the package as is, without making any changes to the code.
             
   b. A more basic setup for running the model through a command line 
      interface can be created by running `make`. The resulting `flu` 
      binary is installed in the `/bin` directory.

4. Register environment as jupyter kernel (optional)

   You may need to manually add the new environment as a jupyter kernel. With the
   conda environment activated, this can be done as follows:

     `python -m ipykernel install --user --name [environment name] --display-name "Python ([environment name])"`




## Getting started

Model parameters are defined in a pandas dataframe (in the python interface). 
Future versions may also support parameter definitions via configuration
files (see a description of parameters in the Configuration section below). An
example of this configuration is in the `gpuflu/flusim.py` file, which can be 
executed as:
```
python3 ./flusim.py -c demo
```

Additional examples are in the `/examples` directory.

Command line execution of the core model is done by executing `./bin/flu`. For
a list of parameters and configuration of the command line interface run:

```
./bin/flu --help
```




## Configuration

**Disease:**
- `r0`                  : mean of Ro. 
- `r0_stddev`           : standard deviation of Ro.
- `exposure_duration`   : mean exposed duration (in time-step units).
- `exposure_stddev`     : variance of exposed duration.
- `infectious_duration` : mean infectious duration (in time-step units).
- `infectious_stddev`   : variance of infectious duration.

NOTE: All disease parameters are non-negative floating point numbers, except for 
      `r0` that is a floating point number larger than 0.

**Population and network**
- `population`  : Population size N (positive integer).
- `connections` : Maximum number of neighbors C that a node can have (positive 
                  integer).
- `attributes`  : Number of attributes M that each node has.
- `network`     : Name of file containing the network specification. This 
                  parameter is optional. If not provided, or if the file name 
                  does not exist, a random network is generated. If the file 
                  with the name specified by this parameter does not exist, then 
                  network is used to set the name of an output file containing 
                  the randomly generated network.
                  
  The network file consists of three parts (each part follows after the last 
  one). Given a population of size N with a maximum degree of C (maximum number
  of outgoing edges), the 3 sequential parts of this file are:
                  
  - Edge list: N rows of C 32-bit unsigned integers (`uint32_t`). The C entries
    of the ith row represent the potential infectees of individual i.

  - Edge weights: N rows of C 16-bit floating point values (`__half float`). 
    The C entries of the ith row represent the strength of connection 
    (transmission) from individual i to the corresponding individual whose id 
    is in the same slot in part 1.

  - Node attributes: N rows of M 16-bit unsigned integers. The M entries of the 
    ith row represent the value of attributes for node i. Examples of 
    attributes can be the age, zip-code, and household id associated to each 
    node.

**Initial and final conditions**
- `infections` : Initial infections.
- `exposures`  : Initial exposures.
- `time`       : Number of timesteps.
- `seed`       : Random number generator seed (integer).


**Execution mode/control**
- `quiet`: Disable verbose mode.
- `debug`: Enable debugging messages.


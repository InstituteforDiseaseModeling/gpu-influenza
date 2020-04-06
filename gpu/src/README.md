# gpu-influenza (SEIR)

### Command Line

```
bin/flu:
        --population:#       - population size [32]
        --connections:#      - per agent network degree [4]
        --exposures:#        - number of initial exposures [0]
        --infections:#       - number of initial infections [1]
        --time:#             - number of time steps [10]
        --durationE:#        - mean exposed duration (days) [5]
        --stddevE:#          - exposed duration std dev (days) [1.5]
        --durationI:#        - mean infectious duration (days) [5]
        --stddevI:#          - infectious duration std dev (days) [1.5]
        --r0:#               - mean r0 [2.5]
        --stddevr0:#         - r0 std dev [1.5]
        --seed:#             - random number seed [20190101]
        --input:filename     - name of file with network data [ none ]
        --output:filename    - name of file for transmission data [ 'transmissions.bin' ]
        --network:filename   - name of file for saving network data [ none ]
        --quiet              - do not display state during execution
        --debug              - display debugging information
        --vaccinations:#     - number of vaccination events (must be passed before other vaccine parameters) [0]
        --vacTimes:#,#       - times of vaccinaton events (must be a single number or comma separated list) [none]
        --vacEffect:#,#      - effectiveness of vaccine (must be a single number or comma separated list) [none]
        --vacCoverage:#,#    - coverage of vaccination (must be a single number or comma separated list) [none]
        --importMag:#        - magitude of the importations (model assumes exponential growth of the source) [0]
        --importDoubleTime:# - doubling time for the importations (model assumes exponential growth of the source) [7]
        --durationImport:#   - mean exposure time (left) for imported casses (days) [2.5]
        --stddevImport:#     - exposed time std dev of the imported cases (days) [1.5]
        --reseedTime:#       - maximum number of days without infection before ending sim (days) [0]
```
Options may be truncated, e.g. ```--pop:1024```, except ```--i``` and ```--in``` will set initial infections rather than set input filename and ```--d``` will set durationE rather than durationI or the debug flag. The truncated form corrisponds to the first match (based on order in the above list), because of this I may rename some of these options to make the truncation more useful.
Options may use ':' or '=', e.g. ```--pop:1024``` and ```--pop=1024``` are equivalent.

The latent period (exposed period), infectious period, and R0 (really beta, but the input is R0) are all drawn from lognormal distributions. The values to set here are the geometric parameters (median and scatter) not the mean and standard deviation of the logarithm.

## Data Structures

```cpp
struct SimulationState {
    uint32_t        population_size;
    uint32_t        max_connections;
    IND_TYPE*       network_m;
    CON_TYPE*       connect_m;
    curandState*    prngStates;
    IND_TYPE*       exposed_individuals_odd_time_steps;
    IND_TYPE*       exposed_individuals_even_time_steps;
    IND_TYPE*       infected_individuals_odd_time_steps;
    IND_TYPE*       infected_individuals_even_time_steps;
    TIMER_TYPE*     exposure_timers_odd_time_steps;
    TIMER_TYPE*     exposure_timers_even_time_steps;
    TIMER_TYPE*     infected_timers_odd_time_steps;
    TIMER_TYPE*     infected_timers_even_time_steps;
    INF_TYPE*       status_v;
    TxEvent*        tx_events;
    Counters*       counters;
};
```
|field|description|
|-----|-----------|
|population_size|scalar - number of agents (nodes) in the simulation|
|max_connections|scalar - maximum node degree|
|network_m|pointer to the "matrix" of edges - see below|
|connect_m|pointer to the "matrix" of connection weights - see below|
|prngStates|pointer to ```population_size``` ```curandState``` structs, one per agent|
|exposed_individuals_odd_time_steps|pointer to up to ```population_size``` entries of exposed agent indices, double buffered with ```exposed_individuals_even_time_steps```|
|exposed_individuals_even_time_steps|pointer to up to ```population_size``` entries of exposed agent indices, double buffered with ```exposed_individuals_odd_time_steps```|
|infected_individuals_odd_time_steps|pointer to up to ```population_size``` entries of infectious agent indices, double buffered with ```infected_individuals_even_time_steps```|
|infected_individuals_even_time_steps|pointer to up to ```population_size``` entries of infectious agent indices, double buffered with ```infected_individuals_odd_time_steps```|
|exposure_timers_odd_time_steps|pointer to up to ```population_size``` exposed (latent) timers, one for each entry of ```exposed_individuals_odd_time_steps```, double buffered with ```exposure_timers_even_time_steps```|
|exposure_timers_even_time_steps|pointer to up to ```population_size``` exposed (latent) timers, one for each entry of ```exposed_individuals_even_time_steps```, double buffered with ```exposure_timers_odd_time_steps```|
|infected_timers_odd_time_steps|pointer to up to ```population_size``` infection timers, one for each entry of ```infected_individuals_odd_time_steps```, double buffered with ```infected_timers_even_time_steps```|
|infected_timers_even_time_steps|pointer to up to ```population_size``` infection timers, one for each entry of ```infected_individuals_even_time_steps```, double buffered with ```infected_timers_odd_time_steps```|
|status_v|pointer to ```population_size``` entries of "status" indicating the state (susceptible: 0, exposed: 1, infectious: 2, recovered: 3) of i<sup>th</sup> agent|
|tx_events|pointer to population_size entries of ```TxEvent``` records|
|counters|pointer to ```Counters``` struct|

### Network "Matrix"
Consists of ```population_size``` rows of ```max_connections``` entries of indices to connected agents.
The j<sup>th</sup> entry of the i<sup>th</sup> row gives the index/id of a potential recipient of infection from agent i.

### Connection Weight "Matrix"
Consists of ```population_size``` rows of ```max_connections``` entries of weights to connected agents.
The j<sup>th</sup> entry of the i<sup>th</sup> row gives the connection weight from agent i to agent j.

### TxEvent struct
|field|notes|
|-----|-----|
|source|infection source index/id|
|target|infection recipient index/id|
|time_step|time step of infection transmission event|

## Network Specification File

A saved network file consists of two parts with part two following directly after part one:
Given a population of size N with a maximum degree of C (maximum number of outgoing edges)
1. edge specification
  - N rows of C 32-bit unsigned integers (uint32_t). The C entries of the ith row represent the potential infectees of individual i.
2. edge weights
  - N rows of C 16-bit floating point values (__half float). The C entries of the ith row represent the strength of connection (transmission) from individual i to the corresponding individual whose id is in the same slot in part 1

Example (N=6, C=2, and all weights 0.1):
```
0 -> 1 -> 3
 \-> 2 -> 4
      \-> 5
```
| C0| C1|
|---|---|
| edges |
| 1 | 2 |
| 3 | 0 |
| 4 | 5 |
| 0 | 0 |
| 0 | 0 |
| 0 | 0 |
| weights |
| 0.1 | 0.1 |
| 0.1 | 0.0 |
| 0.1 | 0.1 |
| 0.0 | 0.0 |
| 0.0 | 0.0 |
| 0.0 | 0.0 |


Currently this is written directly to disk with no metadata so information about N and C must be recorded separately.

## Transmission Recording on GPU

Some options:

* Given a total population, P, assume no more than P transmissions (might not be true since every individual may become infected at some point and there is the possibility of >1 transmission to a particular individual on a given time step).
  * Proposal: allocate P transmission event slots and record as they occur.
* Given a total population, P, and maximum network degree, C, allocate space for PxC transmissions (worst case).
  * Proposal #1: allocate PxC transmission event slots and record as they occur. Requires little post-processing.
  * Proposal #2: allocate memory for PxC transmission times, record time of transmission from individual A to individual B (must be connected in the network) in the corresponding slot. Requires some post-processing - scan for transmission times (>= 0), source is implied in the index to the array, destination has to be looked up from the connection matrix.

## Attributes and Interventions

Attributes are ATTR_TYPE (currently uint16_t) values assigned to each node. Attributes for a simulation are specified by a combination of an entry in a text file (`attributes.txt` for example) which consists of key:value pairs where the key is the attribute name used in the code and referenced by interventions and the value is the filename for the file containing the attribute data. Use the `--attributes:` command line option with a filename to specify the attributes file (can be shortened to `--attr:`). See format/schema below.

There are currently 4 interventions (vaccination is configured separately with command line options):

1. targeting nodes by a particular attribute/attribute value (`ATTRIBUTE` intervention).
2. targeting nodes by attribute/attribute value matching the attribute/attribute value of a currently infected node (`NODEMATCHING` intervention). E.g., if there are two infected nodes, one with "age" value 0 and another with "age" value 2, target all nodes with "age" value 0 or 2.
3. targeting all neighbors of currently infected nodes (`NEIGHBOR` intervention).
4. targeting neighbors of currently infected nodes with attribute/attribute value matching the infected node's attribute/attribute value (`NEIGHBROMATCHING` intervention). E.g., if an infected node has "daytime" attribute value 42, target the neighbors of that node which also have "daytime" value 42.

Interventions are specified with the `--campaign:` (also `--camp:`) command line argument with a filename. The filename specifies a text CSV file with the desired interventions. See format/schema below.

### Sample attribute.txt file

Attribute files may contain comment lines prefixed with `#`.

```text
# Age bin buckets, 0 = <20 (25%), 1 = 20-60 (50%), 2 = >60 (25%)
age:ages.bin
# Household identifiers, nodes are assigned to households with household sizes selected from the
# following distribution:
# 1 - 28.37%, 2 - 34.51%, 3 - 15.07%, 4 - 12.76%, 5 - 05.78%, 6 - 02.26%, 7 - 01.25%
household:households.bin
# Locales are determined by the following formula: schools ~1K students, workplaces ~100 employees, "hangouts" ~50 participants
daytime:locales.bin
```

### Attribute data

Attribute data is a plain binary file with `population_size` entries of ATTR_TYPE (see `flu-types.h`, currently `uint16_t`).

For example, if `population_size` is 10,000 the attribute data file, e.g., `ages.bin`, should be 20,000 bytes in size with 10,000 `uint16_t` entries.

### Interventions CSV file

A campaign (interventions) file consists of a header (required), optional single line comments, and one or more intervention specifications:

```text
tag,type,trigger,trigger_value,duration,factor,attribute,attribute_value
# When infected count >=2500, reduce transmission from all nodes with "ages" attribute == 1 for 60 days (time steps).
threshold,ATTRIBUTE,COUNT,2500,60,0.0625,ages,1
```

# Attribute Intervention

Select all nodes with a particular value for a given attribute, e.g., 42 for "household."
Reduce transmission probability for those nodes by a fixed factor.

## Scenario

- Population 64K (Bothell-ish), network degree 32, initial infections 4, duration 180
- **Baseline:** Assign individuals to <20, 20-60, and >60 age groups (values 0, 1, and 2 respectively).
- **Option 1:** Reduce transmission for working age individuals, 20-60 or attribute value 1, by 15/16 (0.0625) when infectious **count** passes 2500.
- **Option 2:** Reduce transmission for all individuals by 7/8 (0.125) at **time** 40.

`../../../bin/flu --pop:65536 --con:21 --inf:4 --time:180 --output:tx-baseline.bin --quiet`  
`../../../bin/flu --pop:65536 --con:21 --inf:4 --time:180 --attr:ages.txt --campaign:option-one.csv --output:tx-one.bin --quiet`  
`../../../bin/flu --pop:65536 --con:21 --inf:4 --time:180 --attr:ages.txt --campaign:option-two.csv --output:tx-two.bin --quiet`

### attributes.txt

Build with `age-distribution.py`

```
ages:ages.bin
```

### option-one.csv

```csv
tag,type,trigger,trigger_value,duration,factor,attribute,attribute_value
threshold,1,COUNT,2500,60,0.0625,ages,1
```

### option-two.csv

```csv
tag,type,trigger,trigger_value,duration,factor,attribute,attribute_value
scheduled,1,TIME,40,60,0.125,ages,1
```
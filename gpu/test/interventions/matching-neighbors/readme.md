# Neighbor Matching Intervention

Select all neighbroing nodes matching infectious nodes with a given attribute, e.g., "primary" (day-time location).
Reduce transmission probability for those nodes by a fixed factor.

## Scenario

- Population 64K (Bothell-ish), network degree 32, initial infections 4, duration 180
- **Baseline:** Assign individuals to <20, 20-60, and >60 age groups (values 0, 1, and 2 respectively). Subsequently assign to primary groups of ~128 individuals (~512 groups).
- **Option 1:** Reduce transmission for matching neighboring nodes by 15/16 (0.0625) when infectious **count** passes 250.
- **Option 2:** Reduce transmission for matching neighboring nodes by 7/8 (0.125) at **time** 25.

`../../../bin/flu --pop:65536 --con:32 --inf:4 --time:180 --input:network.bin --output:tx-baseline.bin --quiet`  
`../../../bin/flu --pop:65536 --con:32 --inf:4 --time:180 --input:network.bin --attr:attributes.txt --campaign:option-one.csv --output:tx-one.bin --quiet`  
`../../../bin/flu --pop:65536 --con:32 --inf:4 --time:180 --input:network.bin --attr:attributes.txt --campaign:option-two.csv --output:tx-two.bin --quiet`

### attributes.txt

Build with `ages-and-groups.py`

```
ages:ages.bin
primary:primary.bin
```

### network.bin

Build with `network.py`

### option-one.csv

```csv
tag,type,trigger,trigger_value,duration,factor,attribute,attribute_value
threshold,4,COUNT,250,60,0.0625,primary,0
```

### option-two.csv

```csv
tag,type,trigger,trigger_value,duration,factor,attribute,attribute_value
scheduled,4,TIME,25,60,0.125,primary,0
```

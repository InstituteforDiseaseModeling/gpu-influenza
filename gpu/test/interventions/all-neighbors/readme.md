# [All] Neighbors Intervention

Select all neighboring nodes for infectious nodes.
Reduce transmission probability for those nodes by a fixed factor.

## Scenario

- Population 64K (Bothell-ish), network degree 32, initial infections 4, duration 180
- **Baseline:** Random network, no interventions.
- **Option 1:** Reduce transmission for neighboring nodes, by 15/16 (0.0625) when infectious **count** passes 250.
- **Option 2:** Reduce transmission for all individuals by 7/8 (0.125) at **time** 25.

`../../../bin/flu --pop:65536 --con:32 --inf:4 --time:180 --output:tx-baseline.bin --quiet`  
`../../../bin/flu --pop:65536 --con:32 --inf:4 --time:180 --campaign:option-one.csv --output:tx-one.bin --quiet`  
`../../../bin/flu --pop:65536 --con:32 --inf:4 --time:180 --campaign:option-two.csv --output:tx-two.bin --quiet`

### option-one.csv

```csv
tag,type,trigger,trigger_value,duration,factor,attribute,attribute_value
threshold,2,COUNT,250,60,0.0625,n/a,0
```

### option-two.csv

```csv
tag,type,trigger,trigger_value,duration,factor,attribute,attribute_value
scheduled,2,TIME,25,60,0.125,n/a,0
```
# Intervention Expiration

Test "work from home" like intervention where everyone in the 20-60 age group has reduced transmission.
Option #1, expires after two weeks.
Option #2, expires after 60 days.


## Scenario

- Population 64K (Bothell-ish), network degree 32, initial infections 4, duration 180
- **Baseline:** Assign individuals to <20, 20-60, and >60 age groups (values 0, 1, and 2 respectively).
- **Option 1:** Reduce transmission for working age individuals, 20-60 or attribute value 1, by 15/16 (0.0625) when infectious **count** passes 2500, for **14 days**.
- **Option 2:** Reduce transmission for working age individuals, 20-60 or attribute value 1, by 15/16 (0.0625) when infectious **count** passes 2500, for **60 days**.

`../../../bin/flu --pop:65536 --con:32 --inf:4 --time:180 --output:tx-baseline.bin --quiet`  
`../../../bin/flu --pop:65536 --con:32 --inf:4 --time:180 --attr:attributes.txt --campaign:option-one.csv --output:tx-one.bin --quiet`  
`../../../bin/flu --pop:65536 --con:32 --inf:4 --time:180 --attr:attributes.txt --campaign:option-two.csv --output:tx-two.bin --quiet`

### attributes.txt

Build with `age-distribution.py`

```
# age attribute, distributed 1/4, 1/2, 1/4 to <20, 20-60, >60
ages:ages.bin
```

### option-one.csv

```csv
tag,type,trigger,trigger_value,duration,factor,attribute,attribute_value
# reduce transmission from working age individuals (age attribute == 1) for 14 days when infectious count >= 2500
threshold,attribute,COUNT,2500,14,0.0625,ages,1
```

### option-two.csv

```csv
tag,type,trigger,trigger_value,duration,factor,attribute,attribute_value
# reduce transmission from working age individuals (age attribute == 1) for 60 days when infectious count >= 2500
threshold,attribute,COUNT,2500,60,0.0625,ages,1
```
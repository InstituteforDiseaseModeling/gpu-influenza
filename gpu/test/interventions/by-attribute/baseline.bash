#! /usr/bin/env bash

../../../bin/flu --pop:65536 --con:21 --inf:4 --time:180 --output:tx-baseline.bin --quiet
../../../gpuflu/utils/tx-to-seir.py -i tx-baseline.bin -o baseline-seir.csv

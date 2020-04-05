#! /usr/bin/env bash

../../../bin/flu --pop:65536 --con:32 --inf:4 --time:180 --campaign:option-one.csv --output:tx-one.bin --quiet
../../../gpuflu/utils/tx-to-seir.py -i tx-one.bin -o seir-one.csv

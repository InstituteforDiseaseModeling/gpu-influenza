#! /usr/bin/env bash

../../../bin/flu --pop:65536 --con:32 --inf:4 --time:180 --campaign:option-two.csv --output:tx-two.bin --quiet
../../../gpuflu/utils/tx-to-seir.py -i tx-two.bin -o seir-two.csv

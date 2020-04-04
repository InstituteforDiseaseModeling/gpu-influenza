#! /usr/bin/env bash

../../../bin/flu --pop:65536 --con:21 --inf:4 --time:180 --attr:attributes.txt --campaign:option-one.csv --output:tx-one.bin --quiet
../../../gpuflu/utils/tx-to-seir.py -i tx-one.bin -o one-seir.csv

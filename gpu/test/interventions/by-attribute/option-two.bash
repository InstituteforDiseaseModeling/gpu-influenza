#! /usr/bin/env bash

../../../bin/flu --pop:65536 --con:21 --inf:4 --time:180 --attr:attributes.txt --campaign:option-two.csv --output:tx-two.bin --quiet
../../../gpuflu/utils/tx-to-seir.py -i tx-two.bin -o two-seir.csv

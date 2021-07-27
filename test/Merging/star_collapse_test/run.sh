#!/bin/bash

rm *.png
rm -r Dir*


echo "Running the Star Collapse test"

export CHARM_HOME=~/codes/charm/
export ENZO_E_HOME=~/codes/enzo-e/
export CELLO_ARCH="linux_gnu"
export CELLO_PREC="double"

$CHARM_HOME/bin/charmrun +p16 ++local $ENZO_E_HOME/bin/enzo-e star_collapse_test.in >star_collapse_test.out 2>&1 

echo "Standard output in star_collapse_test.out"

mpirun -np 16 python make_images.py > make_images.out 2>&1 

python data.py > data.out 2>&1
python multiple.py > multiple.out 2>&1

#!/bin/bash
mg5=/home/E/chaikangyu/packages/MG5_aMC_v2_6_7
dest=/home/E/chaikangyu/work/Hpair/events/lhe

a=3
for i in 1 2 3 4 5 6 7
do
    if [ $i -lt $a ]
    then
        echo "i<3"
    else
        echo "i>=3"
    fi
done

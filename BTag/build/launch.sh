#!/bin/sh

rootPath=/home/E/chaikangyu/work/Hpair/events/root

for i in 10 15 20 25
do
./BTagEff $rootPath/eeza_jj_40_${i}.root &
done

./BTagEff $rootPath/eeza_jj_50_10.root &
./BTagEff $rootPath/eeza_jj_50_15.root &
wait
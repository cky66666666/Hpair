#!/bin/sh

rootPath=/home/E/chaikangyu/work/Hpair/events/root

./BTagEff $rootPath/sig_aa_1_10_28.root &
./BTagEff $rootPath/sig_aa_1_10_28_1.root &
wait
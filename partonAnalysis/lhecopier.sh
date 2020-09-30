#!/bin/bash
lhepath=/mnt/d/packages/MG5_aMC_v2_6_7/hhj4bSignal/Events
destpath=/mnt/d/work/Hpair/events/lhe
for ((i=1;i<=9;i++))
do
    cd $lhepath/run_0$i
    gzip -d unweighted_events_decayed.lhe.gz
    mv unweighted_events_decayed.lhe hhj$i.lhe
    cp hhj$i.lhe $destpath
done
for ((i=10;i<=11;i++))
do
    cd $lhepath/run_$i
    gzip -d unweighted_events_decayed.lhe.gz
    mv unweighted_events_decayed.lhe hhj$i.lhe
    cp hhj$i.lhe $destpath
done

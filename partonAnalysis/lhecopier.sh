#!/bin/bash
lhepath=/mnt/d/packages/MG5_aMC_v2_6_7/hhj4bSignal/Events
#lhepath=/mnt/c/Users/CHAIKANGYU/Desktop/MG5_aMC_v2_6_6/hhj5/Events
destpath=/mnt/d/work/Hpair/events/lhe
for ((i=2;i<=4;i++))
do
    cd $lhepath/run_0$i
    gzip -d unweighted_events_decayed.lhe.gz
    mv unweighted_events_decayed.lhe sig_tt_$[i]_10.lhe
    #mv hhj$i.lhe hhj$[i-6]_1.lhe
    cp sig_tt_$[i]_10.lhe $destpath
done

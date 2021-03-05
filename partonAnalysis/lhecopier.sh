#!/bin/bash
mg5=/home/E/chaikangyu/packages/MG5_aMC_v2_6_7
dest=/home/E/chaikangyu/work/Hpair/events/lhe

for i in 2 8 9 10 11 12 13 14
do
    #cd $mg5/Events/run_0${i}
    if (($i==2))
    then
        cd $mg5/hhjSignal0/Events/run_02
        gzip -d unweighted_events_decayed.lhe.gz
        mv unweighted_events_decayed.lhe sig_aa_0_10_100.lhe
        cp sig* $dest
    elif (($i>2 && $i<10))
    then
        cd $mg5/hhjSignal/Events/run_0${i}
        gzip -d unweighted_events_decayed.lhe.gz
        mv unweighted_events_decayed.lhe sig_aa_$[i-10]_10_100.lhe
        cp sig* $dest
    else
        cd $mg5/hhjSignal/Events/run_${i}
        gzip -d unweighted_events_decayed.lhe.gz
        mv unweighted_events_decayed.lhe sig_aa_$[i-9]_10_100.lhe
        cp sig* $dest
    fi
done

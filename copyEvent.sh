#!/bin/bash
eventPath=/home/E/chaikangyu/packages/MG5_aMC_v2_6_7/bkg_aa_tthj/Events
rootDest=/home/E/chaikangyu/work/Hpair/events/root
hepmcDest=/home/H/chaikangyu/Hpair/hepmc

for i in 4 5
do
cd $eventPath/run_0${i}
mv tag_1_delphes_events.root bkg_aa_tthj_100_100_$[i-3].root
mv tag_1_pythia8_events.hepmc.gz bkg_aa_tthj_100_100_$[i-3].hepmc.gz
cp bkg_aa_tthj_100_100_$[i-3].root $rootDest/
cp bkg_aa_tthj_100_100_$[i-3].hepmc.gz $hepmcDest/
rm bkg_aa_*
done


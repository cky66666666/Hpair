#!/bin/bash
eventPath=/home/E/chaikangyu/packages/MG5_aMC_v2_6_7/hhjSignal/Events
rootDest=/home/E/chaikangyu/work/Hpair/events/root
hepmcDest=/home/H/chaikangyu/Hpair/hepmc
a=(8)
j=0
for i in 22
do
cd $eventPath/run_${i}
gzip -d unweighted_events_decayed.lhe.gz
mv unweighted_events_decayed.lhe sig_aa_0_${a[j]}_10_100.lhe
cp sig_aa_0_${a[j]}_10_100.lhe /home/E/chaikangyu/work/Hpair/events/lhe/sig_aa_0_${a[j]}_10_100.lhe
#mv tag_1_delphes_events.root bkg_aa_jjaa_100_100_$[i-3].root
#mv tag_1_pythia8_events.hepmc.gz bkg_aa_jjaa_100_100_$[i-3].hepmc.gz
#cp bkg_aa_jjaa_100_100_$[i-3].root $rootDest/
#p bkg_aa_jjaa_100_100_$[i-3].hepmc.gz $hepmcDest/
#m bkg_aa_*
((j+=1))
done


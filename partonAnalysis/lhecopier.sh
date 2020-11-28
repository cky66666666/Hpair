#!/bin/bash
mg5=/home/E/chaikangyu/packages/MG5_aMC_v2_6_7
dest=/home/E/chaikangyu/work/Hpair/events/lhe

cd $mg5/bkg_aa_bbaa/Events/run_01
gzip -d unweighted_events.lhe.gz
mv unweighted_events.lhe bkg_aa_bbaa_10_28.lhe
cp bkg_aa_bbaa_10_28.lhe $dest

cd $mg5/bkg_aa_bbaj/Events/run_01
gzip -d unweighted_events.lhe.gz
mv unweighted_events.lhe bkg_aa_bbaj_10_28.lhe
cp bkg_aa_bbaj_10_28.lhe $dest
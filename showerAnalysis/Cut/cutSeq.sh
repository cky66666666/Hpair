a=1
cutList=("deltaR_bb_max 2.0" "deltaR_bb_min 0.4" "deltaR_aa_max 2.0" "deltaR_aa_min 0.4" \
            "deltaR_ab_max 2.0" "deltaR_ab_min 0.4" "delta_maa 3" "delta_mbb 30" "ptb 40.0" "ptj 150" \
            "topness 1" "topness2 4")
for ((i=1;i<=12;i++))
do
    touch cut${a}_00.txt
    touch cut${a}_01.txt
    touch cut${a}_02.txt
    touch cut${a}_10.txt
    touch cut${a}_20.txt
    for ((j=0;j<i;j++))
    do
        echo ${cutList[j]} >> cut${a}_00.txt
        echo ${cutList[j]} >> cut${a}_01.txt
        echo ${cutList[j]} >> cut${a}_02.txt
        echo ${cutList[j]} >> cut${a}_10.txt
        echo ${cutList[j]} >> cut${a}_20.txt
    done
    echo "fakePhoton 0" >> cut${a}_00.txt
    echo "fakebJet 0" >> cut${a}_00.txt
    echo "fakePhoton 0" >> cut${a}_01.txt
    echo "fakebJet 1" >> cut${a}_01.txt
    echo "fakePhoton 0" >> cut${a}_02.txt
    echo "fakebJet 2" >> cut${a}_02.txt
    echo "fakePhoton 1" >> cut${a}_10.txt
    echo "fakebJet 0" >> cut${a}_10.txt
    echo "fakePhoton 2" >> cut${a}_20.txt
    echo "fakebJet 0" >> cut${a}_20.txt
    ((a+=1))
done


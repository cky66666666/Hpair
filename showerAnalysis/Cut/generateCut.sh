a=1
pt=(30 50 75 100 125)
for ((i=1;i<=4;i++))
do
    for ((j=0;j<i;j++))
    do
        touch cut${a}.txt
        for cut in "deltaR_bb_max 2.0" "deltaR_bb_min 0.4" "deltaR_aa_max 2.0" "deltaR_aa_min 0.4" \
            "deltaR_ab_max 2.0" "deltaR_ab_min 0.4" "delta_maa 3" "delta_mbb 35" "ptb 40.0" "ptj 150" \
            "nlpt ${pt[i]}" "nnlpt ${pt[j]}"
        do
            echo $cut >> cut${a}.txt
        done
        ((a+=1))
    done
done
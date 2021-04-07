i=37
for deltaR_bb_min in 0.3 0.4 0.5
do
    for deltaR_aa_min in 0.3 0.4 0.5
    do
        touch cut${i}.txt
        for cut in "deltaR_bb_max 2.0" "deltaR_bb_min ${deltaR_bb_min}" "deltaR_aa_max 2.0" "deltaR_aa_min ${deltaR_aa_min}" "deltaR_ab_max 2.0" "deltaR_ab_min 0.4" "delta_maa 3" "delta_mbb 25" "ptb 40.0" "ptj 100"
        do
            echo $cut >> cut${i}.txt
        done
        ((i+=1))
    done
done
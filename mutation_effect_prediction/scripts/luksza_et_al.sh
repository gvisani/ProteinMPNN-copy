

pdb_dir='/gscratch/spe/gvisan01/mutation_effect_predictions/luksza_et_al/pdbs/'

base_dir='/gscratch/spe/gvisan01/mutation_effect_predictions/luksza_et_al/'

# systems='luksza_et_al_tcr1_ec50_sat_mut luksza_et_al_tcr1_ec50_sat_mut_FILTERED_HIGH_EC50s'
# systems='luksza_et_al_tcr2_ec50_sat_mut luksza_et_al_tcr3_ec50_sat_mut luksza_et_al_tcr4_ec50_sat_mut luksza_et_al_tcr5_ec50_sat_mut luksza_et_al_tcr6_ec50_sat_mut'
systems='luksza_et_al_tcr7_ec50_sat_mut'

for system in $systems
    do

    python -u ../zero_shot_mutation_effect_prediction.py \
                    --csv_file $base_dir$system'.csv' \
                    --folder_with_pdbs $pdb_dir \
                    --output_dir $base_dir'proteinmpnn' \
                    --use_mt_structure 0 \
                    --num_seq_per_target 10 \
                    --batch_size 10 \
                    --mutant_column mutant \
                    --mutant_chain_column mutant_chain

done


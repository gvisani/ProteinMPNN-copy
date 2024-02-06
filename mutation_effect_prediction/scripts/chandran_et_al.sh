

pdb_dir='/gscratch/spe/gvisan01/mutation_effect_predictions/chandran_et_al/pdbs/'

base_dir='/gscratch/spe/gvisan01/mutation_effect_predictions/chandran_et_al/'

systems='chandran_et_al_peptide_A_and_G_scans'

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


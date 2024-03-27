


num_seq_per_target='10'
use_mt_structure_list='1'  # '0 1'

pdb_dir='/gscratch/spe/gvisan01/tcr_pmhc/pdbs/ATLAS/'

base_dir='/gscratch/spe/gvisan01/mutation_effect_predictions/ATLAS/'

for use_mt_structure in $use_mt_structure_list
    do

    python -u ../zero_shot_mutation_effect_prediction.py \
                    --csv_file $base_dir'ATLAS_cleaned.csv' \
                    --folder_with_pdbs $pdb_dir \
                    --output_dir $base_dir'proteinmpnn' \
                    --use_mt_structure $use_mt_structure \
                    --num_seq_per_target $num_seq_per_target \
                    --batch_size 10 \
                    --dms_column '[ log10(Kd_wt/Kd_mt) ]' \
                    --mutant_column mutant \
                    --mutant_chain_column chain \
                    --dont_run_inference 1
    
done



# num_seq_per_target_list='100 10 1'
# use_mt_structure_list='0 1'

num_seq_per_target=10
use_mt_structure=0

base_dir='/gscratch/spe/gvisan01/mutation_effect_predictions/T4_ddG/'

python -u $base_dir'standardize_T4_ddG.py' --csv_file $base_dir'T4_mutant_ddG.csv'

python -u ../zero_shot_mutation_effect_prediction.py \
                --csv_file $base_dir'T4_mutant_ddG_standardized.csv' \
                --folder_with_pdbs $base_dir'pdbs' \
                --output_dir $base_dir'proteinmpnn' \
                --use_mt_structure $use_mt_structure \
                --num_seq_per_target $num_seq_per_target \
                --batch_size 10 \
                --dms_column ddG \
                --mutant_column mutant \
                --mutant_chain_column mutant_chain \
                --dont_run_inference 1

# python -u ../zero_shot_mutation_effect_prediction.py \
#                 --csv_file $base_dir'T4_mutant_ddG_no_G_and_A_standardized.csv' \
#                 --folder_with_pdbs $base_dir'pdbs' \
#                 --output_dir $base_dir'proteinmpnn' \
#                 --use_mt_structure $use_mt_structure \
#                 --num_seq_per_target $num_seq_per_target \
#                 --batch_size 10 \
#                 --dms_column ddG \
#                 --mutant_column mutant \
#                 --mutant_chain_column mutant_chain \

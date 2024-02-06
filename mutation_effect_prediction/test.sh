


base_dir='/gscratch/spe/gvisan01/mutation_effect_predictions/SKEMPI/'

python -u zero_shot_mutation_effect_prediction.py \
                --csv_file ./test.csv \
                --folder_with_pdbs /gscratch/spe/kborisia/learning/DSMbind_HCNN_comparison/PDBs \
                --output_dir ./output \
                --use_mt_structure 0 \
                --num_seq_per_target 10 \
                --batch_size 10 \
                --wt_pdb_column PDB_filename \
                --mutant_column mutant \
                --mutant_chain_column mutant_chain










folder_with_pdbs="./pdbs/tax/"

path_to_fasta="./fasta_files/tax/1ao7_copy.fa"

output_dir="./outputs/tax/"
if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi


path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
chains_to_design="C"
#The first amino acid in the chain corresponds to 1 and not PDB residues index for now.
design_only_positions="1 2 3 4 5 6 7 8 9" #design only these residues; use flag --specify_non_fixed

python ../helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python ../helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

python ../helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$design_only_positions" --specify_non_fixed

python ../protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 100 \
        --sampling_temp "0.5" \
        --seed 37 \
        --batch_size 25 \
        --score_only 1


# python ../protein_mpnn_run.py \
#         --pdb_path "./pdbs/tax/1ao7.pdb" \
#         --pdb_path_chains "C" \
#         --out_folder $output_dir \
#         --num_seq_per_target 1 \
#         --sampling_temp "0.1" \
#         --seed 37 \
#         --batch_size 1 \
#         --path_to_fasta $path_to_fasta \
#         --score_only 1 


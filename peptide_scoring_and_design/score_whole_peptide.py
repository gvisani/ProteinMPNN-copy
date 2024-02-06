
import os, sys
import json

import argparse
from typing import *

def load_jsonl(jsonl_path: str) -> List[Dict]:
    with open(jsonl_path, 'r') as f:
        json_list = list(f)
    dict_list = []
    for json_str in json_list:
        dict_list.append(json.loads(json_str))
    return dict_list

def dump_jsonl(jsonl_data: List[Dict], jsonl_path: str) -> None:
    with open(jsonl_path, 'w+') as f:
        for entry in jsonl_data:
            f.write(json.dumps(entry) + '\n')

def set_peptide_in_parsed_jsonl(jsonl_input_path: str, peptide_seq: str, peptide_chain: str) -> str:
    '''
    Saves updated jsonl in a file in the same directory and with the same name, except for adding the peptide sequence to the end of the filename.
    Returns the path of said file.
    '''

    parsed_pdb_list = load_jsonl(jsonl_input_path)

    for parsed_pdb in parsed_pdb_list:

        if len(parsed_pdb[f'seq_chain_{peptide_chain}']) != len(peptide_seq):
            raise ValueError(f'Length of provided peptide ({len(peptide_seq)}) mismatches length of peptide in PDB ({len(parsed_pdb[f"seq_chain_{peptide_chain}"])}). Perhaps the chain (chain {peptide_chain}) is wrong?')

        # substitute in the desired peptide sequence
        parsed_pdb[f'seq_chain_{peptide_chain}'] = peptide_seq

        # change the peptide sequence within the overall sequence
        seq_chain_keys = sorted([key for key in parsed_pdb if key.startswith('seq_chain_')])
        parsed_pdb['seq'] = ''.join([parsed_pdb[key] for key in seq_chain_keys])

        # change the name to include the peptide as well --> this will imact the name of the output file so it includes the peptide sequence
        parsed_pdb['name'] = parsed_pdb['name'] + '_' + peptide_seq

    dump_jsonl(parsed_pdb_list, jsonl_input_path)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--folder_with_pdbs', type=str, required=True) # can score the peptide in multiple PDBs!
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--peptide_seq', type=str, required=True)
    parser.add_argument('--peptide_chain', type=str, default='C')

    # ProteinMPNN arguments
    parser.add_argument('--num_seq_per_target', type=int, default=10)
    parser.add_argument('--batch_size', type=int, default=10)

    args = parser.parse_args()

    print(os.listdir(args.folder_with_pdbs))

    os.makedirs(args.output_dir, exist_ok=True)

    path_for_parsed_chains = os.path.join(args.output_dir, 'parsed_pdbs.jsonl')
    path_for_assigned_chains = os.path.join(args.output_dir, 'assigned_pdbs.jsonl')
    path_for_fixed_positions = os.path.join(args.output_dir, 'fixed_pdbs.jsonl')

    # script call to parse pdbs
    os.system(f'python ../helper_scripts/parse_multiple_chains.py --input_path={args.folder_with_pdbs} --output_path={path_for_parsed_chains}')

    # change the amino-acid sequence of the peptidfe within the parsed pdb
    set_peptide_in_parsed_jsonl(path_for_parsed_chains, args.peptide_seq, args.peptide_chain)

    # other proteinmpnn helper scripts to make it focus only on the peptide
    os.system(f'python ../helper_scripts/assign_fixed_chains.py --input_path={path_for_parsed_chains} --output_path={path_for_assigned_chains} --chain_list "{args.peptide_chain}"')

    os.system(f'python ../helper_scripts/make_fixed_positions_dict.py --input_path={path_for_parsed_chains} --output_path={path_for_fixed_positions} --chain_list "{args.peptide_chain}" --position_list "{" ".join(str(i) for i in range(1, len(args.peptide_seq)+1))}" --specify_non_fixed')

    # score!
    os.system(f'python ../protein_mpnn_run.py \
                                --jsonl_path {path_for_parsed_chains} \
                                --chain_id_jsonl {path_for_assigned_chains} \
                                --fixed_positions_jsonl {path_for_fixed_positions} \
                                --out_folder {args.output_dir} \
                                --num_seq_per_target {args.num_seq_per_target} \
                                --batch_size {args.batch_size} \
                                --score_only 1')


# chains_to_design="C"
# #The first amino acid in the chain corresponds to 1 and not PDB residues index for now.
# design_only_positions="1 2 3 4 5 6 7 8 9" #design only these residues; use flag --specify_non_fixed

# python ../helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

# python ../helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$design_only_positions" --specify_non_fixed

# python ../protein_mpnn_run.py \
#         --jsonl_path $path_for_parsed_chains \
#         --chain_id_jsonl $path_for_assigned_chains \
#         --fixed_positions_jsonl $path_for_fixed_positions \
#         --out_folder $output_dir \
#         --num_seq_per_target 100 \
#         --sampling_temp "0.5" \
#         --seed 37 \
#         --batch_size 25 \
#         --score_only 1
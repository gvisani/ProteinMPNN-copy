

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO

import logomaker

import argparse

curr_dir = os.path.dirname(os.path.realpath(__file__))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--folder_with_pdbs', type=str, required=True)
    parser.add_argument('--output_dir_template', type=str, required=True)
    parser.add_argument('--sampling_temp', type=str, default='0.1')
    parser.add_argument('--num_seq_per_target', type=int, default=100)
    parser.add_argument('--seed', type=int, default=37)
    parser.add_argument('--batch_size', type=int, default=25)

    parser.add_argument('--peptide_chain', type=str, default='C')
    parser.add_argument('--peptide_resnum_start', type=int, default=1)
    parser.add_argument('--peptide_length', type=int, default=9)
    parser.add_argument("--model_name", nargs='+', type=str, default="v_48_020", help="ProteinMPNN model names, can be more than one, which will be saved in different outpur_dirs: v_48_002, v_48_010, v_48_020, v_48_030; v_48_010=version with 48 edges 0.10A noise")

    args = parser.parse_args()

    for model_name in args.model_name:
        output_dir = args.output_dir_template + f'_{model_name}' + f'_{args.sampling_temp}'
        os.makedirs(output_dir, exist_ok=True)

        path_for_parsed_chains = os.path.join(output_dir, 'parsed_pdbs.jsonl')
        path_for_assigned_chains = os.path.join(output_dir, 'assigned_pdbs.jsonl')
        path_for_fixed_positions = os.path.join(output_dir, 'fixed_pdbs.json')
        design_only_positions = ' '.join([str(i) for i in range(args.peptide_resnum_start, args.peptide_length + args.peptide_resnum_start)])

        os.system(f'python {curr_dir}/../helper_scripts/parse_multiple_chains.py --input_path={args.folder_with_pdbs} --output_path={path_for_parsed_chains}')

        os.system(f'python {curr_dir}/../helper_scripts/assign_fixed_chains.py --input_path={path_for_parsed_chains} --output_path={path_for_assigned_chains} --chain_list "{args.peptide_chain}"')

        os.system(f'python {curr_dir}/../helper_scripts/make_fixed_positions_dict.py --input_path={path_for_parsed_chains} --output_path={path_for_fixed_positions} --chain_list "{args.peptide_chain}" --position_list "{design_only_positions}" --specify_non_fixed')

        os.system(f'python {curr_dir}/../protein_mpnn_run.py \
                                --jsonl_path {path_for_parsed_chains} \
                                --chain_id_jsonl {path_for_assigned_chains} \
                                --fixed_positions_jsonl {path_for_fixed_positions} \
                                --out_folder {output_dir} \
                                --model_name {model_name} \
                                --num_seq_per_target {args.num_seq_per_target} \
                                --sampling_temp {args.sampling_temp} \
                                --seed {args.seed} \
                                --batch_size {args.batch_size}')

        # print(model_name)
        # pdbids = [f[:-4] for f in os.listdir(args.folder_with_pdbs) if f.endswith('.pdb')]
        # sequences_per_pdbid = {}
        # for pdbid in pdbids:
        #     print(pdbid)
        #     sequences = []
        #     with open(os.path.join(output_dir, 'seqs', f'{pdbid}.fa'), 'r') as f:
        #         for i, record in enumerate(SeqIO.parse(f, 'fasta')):
        #             if i > 0:
        #                 sequences.append(str(record.seq))
            
        #     sequences = list(set(sequences))
        #     print(f'Total unique sequences: {len(sequences)}')
        #     sequences_per_pdbid[pdbid] = sequences
                
        #     # make pwm, save it as image
        #     pwm = logomaker.alignment_to_matrix(sequences)
        #     pwm_df = pd.DataFrame(pwm, columns=list('ACDEFGHIKLMNPQRSTVWY'))
        #     for i in range(pwm_df.shape[0]):
        #         pwm_df.iloc[i] = pwm_df.iloc[i] / pwm_df.iloc[i].sum()

        #     # set nan values to zero
        #     pwm_df = pwm_df.fillna(0)
            
        #     fig, ax = plt.subplots(figsize=(10, 2))
        #     logomaker.Logo(pwm_df, ax=ax)
        #     plt.tight_layout()
        #     plt.savefig(os.path.join(output_dir, 'seqs', f'{pdbid}_logo.png'))
        #     plt.close()

        # print()


        # ## now put these sequences in a tsv file for tcrdock to use
        # header = 'organism	mhc_class	mhc	peptide	va	ja	cdr3a	vb	jb	cdr3b	pdbid	hcnn_model'
        # template_row = 'human	1	A*02:01	{seq}	TRAV21*01	TRAJ6*01	CAVRPTSGGSYIPTF	TRBV6-5*01	TRBJ2-2*01	CASSYVGNTGELFF	{pdbid}	'
        
        # with open(os.path.join(output_dir, f'proteinmpnn_samples_{model_name}_{args.sampling_temp}.tsv'), 'w+') as f:
        #     f.write(header + '\n')
        #     written_sequences = set()
        #     for pdbid, sequences in sequences_per_pdbid.items():
        #         for seq in sequences:
        #             if seq in written_sequences:
        #                 continue
        #             else:
        #                 f.write(template_row.format(seq=seq, pdbid=pdbid) + '\n')
        #                 written_sequences.add(seq)

        # ## test the tsv file is read successfully
        # df = pd.read_csv(os.path.join(output_dir, f'proteinmpnn_samples_{model_name}_{args.sampling_temp}.tsv'), sep='\t')
        # print(df.head())

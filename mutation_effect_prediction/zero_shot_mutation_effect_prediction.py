
'''

NOTE: ignoring insertions codes for the time being

'''

import os
import json
import gzip, pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
from typing import *

def optional_str(astr: Union[None, str]) -> Union[None, str]:
    if astr == 'None':
        return None
    else:
        return astr

def identify_single_mutation_sites(df, args):
    sites = []
    for i, row in df.iterrows():
        chain = row[args.mutant_chain_column]

        # skip multiple mutations
        if args.mutant_split_symbol in row[args.mutant_column]:
            continue
            
        resnum = int(row[args.mutant_column][1:-1])
        sites.append((chain, resnum))
    
    return sites

def format_chain_and_resnums_for_proteinmpnn(chains, resnums):
    assert len(chains) == len(resnums)
    chain_to_resnums_dict = {}
    for chain, resnum in zip(chains, resnums):
        if chain not in chain_to_resnums_dict:
            chain_to_resnums_dict[chain] = []
        chain_to_resnums_dict[chain].append(resnum)
    chains_str = ' '.join(chain_to_resnums_dict.keys())
    resnums_str = ','.join([' '.join(map(str, resnums)) for resnums in chain_to_resnums_dict.values()])
    return chains_str, resnums_str


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

def order_mutants_with_proteinmpnn_output(path_for_parsed_chains: str, path_for_assigned_chains: str, mutants: List[str], chains_of_mutants: List[str]):
    assert len(mutants) == len(chains_of_mutants), 'The number of mutants and chains should be the same'
    resnums = [int(mutant[1:-1]) for mutant in mutants]
    chain_to_mutants_resnums_tuple = {}
    for mutant, resnum, chain in zip(mutants, resnums, chains_of_mutants):
        if chain not in chain_to_mutants_resnums_tuple:
            chain_to_mutants_resnums_tuple[chain] = []
        chain_to_mutants_resnums_tuple[chain].append((mutant, resnum))
    # sort the mutants by resnums
    for chain in chain_to_mutants_resnums_tuple:
        chain_to_mutants_resnums_tuple[chain] = sorted(chain_to_mutants_resnums_tuple[chain], key=lambda x: x[1])

    parsed_chains = load_jsonl(path_for_parsed_chains)[0]
    assigned_chains = load_jsonl(path_for_assigned_chains)[0]

    assert len(assigned_chains) == 1, 'Only one PDB file should be processed at a time'
    pdb = list(assigned_chains.keys())[0]
    chains_list_of_lists = assigned_chains[pdb]
    
    # chains_in_order = []
    # for elem in chains_list_of_lists:
    #     assert len(elem) in {0, 1}, f"These should be singleton lists or empty lists? Not sure why exactly. This might come back to bite me in the ass in the future... just kidding, it's now! {elem}"
    #     if len(elem) == 1:
    #         chains_in_order.append(elem[0])

    def flatten(xss):
        return [x for xs in xss for x in xs]
    chains_in_order = flatten(chains_list_of_lists)

    # now sort the mutants by chain as well
    chain_mutant_in_order = []
    for chain in chains_in_order:
        if chain in chain_to_mutants_resnums_tuple:
            for mutant, resnum in chain_to_mutants_resnums_tuple[chain]:
                chain_mutant_in_order.append((chain, mutant))

    return chain_mutant_in_order
            



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--csv_file', type=str, required=True)
    parser.add_argument('--folder_with_pdbs', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)

    parser.add_argument('--dont_run_inference', type=int, default=0)

    # ProteinMPNN arguments
    parser.add_argument('--num_seq_per_target', type=int, default=100)
    parser.add_argument('--batch_size', type=int, default=25)

    parser.add_argument('--use_mt_structure', type=int, default=0,
                        help='0 for false, 1 for true. If toggled, compute logits for mutations on the corresponding mutated structure.')

    # arguments specitying how to read the csv file with the mutations
    parser.add_argument('--wt_pdb_column', type=str, default='wt_pdb')
    parser.add_argument('--mt_pdb_column', type=str, default='mt_pdb')
    parser.add_argument('--mutant_column', type=str, default='mutant', help='Column name with the mutation')
    parser.add_argument('--mutant_chain_column', type=optional_str, default='mutant_chain', help='Column name with the chain the mutation occurs on. If None, expects that all mutations occur on the same chain.')
    parser.add_argument('--mutant_split_symbol', type=str, default='|', help='Symbol used to split multiple mutations.')

    parser.add_argument('--dms_column', type=optional_str, default=None) # ddG
    parser.add_argument('--dms_label', type=optional_str, default=None) #r'stability effect, $\Delta\Delta G$')

    ## when predicting a binary outcome (e.g. Destabilizing vs. Neutral, Bind vs. Not-Bind, ...)
    parser.add_argument('--binary_dms_column', type=optional_str, default=None) #'effect')
    parser.add_argument('--binary_dms_pos_value', type=optional_str, default=None) #'Destabilizing')
    parser.add_argument('--binary_dms_label', type=optional_str, default=None) #'Stability Effect\n(Destabilizing = 1, Neutral = 0)')

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, 'zero_shot_predictions'), exist_ok=True)

    df = pd.read_csv(args.csv_file)

    identifier = f'num_seq_per_target={args.num_seq_per_target}-use_mt_structure={args.use_mt_structure}'

    ## single-point mutations: first get 20-dim conditional probas for all sites of interest, then get requested `log(p_mt) - log(p_wt)`

    if args.dont_run_inference:
        with gzip.open(os.path.join(args.output_dir, 'zero_shot_predictions', f'table_of_conditional_probas-{identifier}.gz'), 'rb') as f:
            table_of_conditional_probas = pickle.load(f)
    else:
        table_of_conditional_probas = {
            'WT': {},
            'MT': {}
        }

        path_for_parsed_chains = os.path.join(args.output_dir, 'parsed_pdbs.jsonl')
        path_for_assigned_chains = os.path.join(args.output_dir, 'assigned_pdbs.jsonl')
        path_for_fixed_positions = os.path.join(args.output_dir, 'fixed_pdbs.jsonl')

        seen_sites = set()
        for i, row in tqdm(df.iterrows(), total=len(df)):

            # # skip multiple mutations
            # if args.mutant_split_symbol in row[args.mutant_column]:
            #     continue

            if not isinstance(row[args.mutant_chain_column], str) or not isinstance(row[args.mutant_column], str):
                # sometimes the chain and mutant columns are NaN, in which case we skip the row
                continue
            
            chains = row[args.mutant_chain_column].split(args.mutant_split_symbol)
            mutants = row[args.mutant_column].split(args.mutant_split_symbol)
            resnums = [int(mutant[1:-1]) for mutant in mutants]
            chains_str, resnums_str = format_chain_and_resnums_for_proteinmpnn(chains, resnums)

            wt_pdb = row[args.wt_pdb_column]

            # make temporary directory with single pdb
            temp_pdbs_dir = os.path.join(args.output_dir, "temp_pdbs")
            os.makedirs(temp_pdbs_dir, exist_ok=True)
            os.system(f'cp {os.path.join(args.folder_with_pdbs, wt_pdb + ".pdb")} {temp_pdbs_dir}')

            if (wt_pdb, chains_str, resnums_str) in seen_sites:
                # still need to add mutants to the table_of_conditional_probas, but let us avoid running the whole proteinmpnn
                os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/parse_multiple_chains.py")} --input_path={temp_pdbs_dir} --output_path={path_for_parsed_chains}')
                os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/assign_fixed_chains.py")} --input_path={path_for_parsed_chains} --output_path={path_for_assigned_chains} --chain_list "{chains_str}"')
                chain_mutant_in_order = order_mutants_with_proteinmpnn_output(path_for_parsed_chains, path_for_assigned_chains, mutants, chains)
                if (wt_pdb, chains_str, resnums_str) in table_of_conditional_probas['WT']:
                    table_of_conditional_probas['WT'][(wt_pdb, chains_str, resnums_str)][1].append(chain_mutant_in_order)
                os.system(f'rm -r {temp_pdbs_dir}') # delete temporary directory with single pdb
                os.system(f'rm {path_for_parsed_chains}')
                os.system(f'rm {path_for_assigned_chains}')
                continue
                
            seen_sites.add((wt_pdb, chains_str, resnums_str))
            
            os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/parse_multiple_chains.py")} --input_path={temp_pdbs_dir} --output_path={path_for_parsed_chains}')
            os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/assign_fixed_chains.py")} --input_path={path_for_parsed_chains} --output_path={path_for_assigned_chains} --chain_list "{chains_str}"')
            os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/make_fixed_positions_dict.py")} --input_path={path_for_parsed_chains} --output_path={path_for_fixed_positions} --chain_list "{chains_str}" --position_list "{resnums_str}" --specify_non_fixed')
            os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../protein_mpnn_run.py")} \
                                                                                                --suppress_print 1 \
                                                                                                --jsonl_path {path_for_parsed_chains} \
                                                                                                --chain_id_jsonl {path_for_assigned_chains} \
                                                                                                --fixed_positions_jsonl {path_for_fixed_positions} \
                                                                                                --out_folder {args.output_dir} \
                                                                                                --num_seq_per_target {args.num_seq_per_target} \
                                                                                                --batch_size {args.batch_size} \
                                                                                                --conditional_probs_only 1')
            
            os.system(f'rm -r {temp_pdbs_dir}') # delete temporary directory with single pdb

            # parse output file to get the conditional probas
            output = np.load(os.path.join(args.output_dir, 'conditional_probs_only', wt_pdb + '.npz'))
            chain_mutant_in_order = order_mutants_with_proteinmpnn_output(path_for_parsed_chains, path_for_assigned_chains, mutants, chains)
            os.system(f'rm {path_for_parsed_chains}')
            os.system(f'rm {path_for_assigned_chains}')
            os.system(f'rm {path_for_fixed_positions}')

            if len(output['log_p'].shape) == 3: # why this again? is it the number of sequence per target thing?
                log_p = np.mean(output['log_p'], axis=0)[output['design_mask'].astype(np.bool_)]
            elif len(output['log_p'].shape) == 2:
                log_p = output['log_p'][output['design_mask'].astype(np.bool_)]
            else:
                raise ValueError()
            
            # if log_p.shape[0] > 1:
            #     print(f"Warning: there is a redundancy of {log_p.shape[0]} in logits, taking the first one assuming it's an inserion code thing and that the empty icode is desired..")
            # elif log_p.shape[0] == 0:
            #     print(f"Warning: site {(chains_str, resnums_str)} not found!")
            #     continue
            if log_p.shape[0] == 0:
                print(f"Warning: site {(chains_str, resnums_str)} not found!")
                continue

            assert log_p.shape[0] == len(chain_mutant_in_order), f"Length of log_p and chain_mutant_in_order should be the same, but are {log_p.shape[0]} and {len(chain_mutant_in_order)} respectively."

            # log_p = log_p[0]
            table_of_conditional_probas['WT'][(wt_pdb, chains_str, resnums_str)] = (log_p, [chain_mutant_in_order])


            if args.use_mt_structure:
                mt_pdb = row[args.mt_pdb_column]

                # make temporary directory with single pdb
                temp_pdbs_dir = os.path.join(args.output_dir, "temp_pdbs")
                os.makedirs(temp_pdbs_dir, exist_ok=True)
                os.system(f'cp {os.path.join(args.folder_with_pdbs, str(mt_pdb) + ".pdb")} {temp_pdbs_dir}')

                if (mt_pdb, chains_str, resnums_str) in seen_sites:
                    # still need to add mutants to the table_of_conditional_probas, but let us avoid running the whole proteinmpnn
                    os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/parse_multiple_chains.py")} --input_path={temp_pdbs_dir} --output_path={path_for_parsed_chains}')
                    os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/assign_fixed_chains.py")} --input_path={path_for_parsed_chains} --output_path={path_for_assigned_chains} --chain_list "{chains_str}"')
                    chain_mutant_in_order = order_mutants_with_proteinmpnn_output(path_for_parsed_chains, path_for_assigned_chains, mutants, chains)
                    if (wt_pdb, chains_str, resnums_str) in table_of_conditional_probas['MT']:
                        table_of_conditional_probas['MT'][(wt_pdb, chains_str, resnums_str)][1].append(chain_mutant_in_order)
                    os.system(f'rm -r {temp_pdbs_dir}') # delete temporary directory with single pdb
                    os.system(f'rm {path_for_parsed_chains}')
                    os.system(f'rm {path_for_assigned_chains}')
                    continue

                seen_sites.add((mt_pdb, chains_str, resnums_str))

                os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/parse_multiple_chains.py")} --input_path={temp_pdbs_dir} --output_path={path_for_parsed_chains}')
                os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/assign_fixed_chains.py")} --input_path={path_for_parsed_chains} --output_path={path_for_assigned_chains} --chain_list "{chains_str}"')
                os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/make_fixed_positions_dict.py")} --input_path={path_for_parsed_chains} --output_path={path_for_fixed_positions} --chain_list "{chains_str}" --position_list "{resnums_str}" --specify_non_fixed')
                os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../protein_mpnn_run.py")} \
                                                                                                    --suppress_print 1 \
                                                                                                    --jsonl_path {path_for_parsed_chains} \
                                                                                                    --chain_id_jsonl {path_for_assigned_chains} \
                                                                                                    --fixed_positions_jsonl {path_for_fixed_positions} \
                                                                                                    --out_folder {args.output_dir} \
                                                                                                    --num_seq_per_target {args.num_seq_per_target} \
                                                                                                    --batch_size {args.batch_size} \
                                                                                                    --conditional_probs_only 1')
                
                os.system(f'rm -r {temp_pdbs_dir}') # delete temporary directory with single pdb
                os.system(f'rm {path_for_parsed_chains}')
                os.system(f'rm {path_for_assigned_chains}')
                os.system(f'rm {path_for_fixed_positions}')

                # parse output file to get the conditional probas
                try:
                    output = np.load(os.path.join(args.output_dir, 'conditional_probs_only', str(mt_pdb) + '.npz'))
                except FileNotFoundError:
                    print(f"Warning: file {str(mt_pdb) + '.npz'} not found!")
                    continue
                
                if len(output['log_p'].shape) == 3:
                    log_p = np.mean(output['log_p'], axis=0)[output['design_mask'].astype(np.bool_)]
                elif len(output['log_p'].shape) == 2:
                    log_p = output['log_p'][output['design_mask'].astype(np.bool_)]
                else:
                    raise ValueError()
                
                # assert log_p.shape == (1, 21)
                # log_p = np.squeeze(log_p, axis=0)

                assert log_p.shape[0] == len(chain_mutant_in_order), f"Length of log_p and chain_mutant_in_order should be the same, but are {log_p.shape[0]} and {len(chain_mutant_in_order)} respectively."

                table_of_conditional_probas['MT'][(chains_str, resnums_str)] = (log_p, [chain_mutant_in_order])
        
        with gzip.open(os.path.join(args.output_dir, 'zero_shot_predictions', f'table_of_conditional_probas-{identifier}.gz'), 'wb') as f:
            pickle.dump(table_of_conditional_probas, f)

    alphabet = 'ACDEFGHIKLMNPQRSTVWYX'
    alphabet_dict = dict(zip(alphabet, range(21)))

    ## 2) now iterate over the rows again, and get `log(p_mt) - log(p_wt)` for every row
    log_p_mt = []
    log_p_wt = []
    log_p_mt__minus__log_p_wt = []
    for i, row in df.iterrows():
        # # skip multiple mutations
        # if args.mutant_split_symbol in row[args.mutant_column]:
        #     log_p_mt__minus__log_p_wt.append(np.nan)
        #     continue

        if not isinstance(row[args.mutant_chain_column], str) or not isinstance(row[args.mutant_column], str):
            # sometimes the chain and mutant columns are NaN, in which case we skip the row
            continue

        wt_pdb = row[args.wt_pdb_column]
        
        chains = row[args.mutant_chain_column].split(args.mutant_split_symbol)
        mutants = row[args.mutant_column].split(args.mutant_split_symbol)
        resnums = [int(mutant[1:-1]) for mutant in mutants]
        chains_str, resnums_str = format_chain_and_resnums_for_proteinmpnn(chains, resnums)

        # chain = row[args.mutant_chain_column]
        # resnum = int(row[args.mutant_column][1:-1])

        if (wt_pdb, chains_str, resnums_str) not in table_of_conditional_probas['WT']:
            log_p_mt.append(np.nan)
            log_p_wt.append(np.nan)
            log_p_mt__minus__log_p_wt.append(np.nan)
        elif args.use_mt_structure and (row[args.mt_pdb_column], chains_str, resnums_str) not in table_of_conditional_probas['MT']:
            log_p_mt.append(np.nan)
            log_p_wt.append(np.nan)
            log_p_mt__minus__log_p_wt.append(np.nan)
        else:
            log_p_wt_vec, chain_mutant_in_order_list_wt = table_of_conditional_probas['WT'][(wt_pdb, chains_str, resnums_str)]
            if args.use_mt_structure:
                mt_pdb = row[args.mt_pdb_column]
                log_p_mt_vec, chain_mutant_in_order_list_mt = table_of_conditional_probas['MT'][(mt_pdb, chains_str, resnums_str)]
                assert chain_mutant_in_order_list_wt == chain_mutant_in_order_list_mt, f"chain_mutant_in_order_wt and chain_mutant_in_order_mt should be the same, but are {chain_mutant_in_order_wt} and {chain_mutant_in_order_mt} respectively."
            else:
                log_p_mt_vec = log_p_wt_vec
            
            temp_log_p_wt = []
            temp_log_p_mt = []
            for i, (mutant, chain) in enumerate(zip(mutants, chains)):
                for chain_mutant_in_order_wt in chain_mutant_in_order_list_wt:
                    if (chain, mutant) in chain_mutant_in_order_wt:
                        index = chain_mutant_in_order_wt.index((chain, mutant))
                        break
                aa_wt = mutant[0]
                aa_mt = mutant[-1]
                temp_log_p_wt.append(log_p_wt_vec[index, alphabet_dict[aa_wt]])
                temp_log_p_mt.append(log_p_mt_vec[index, alphabet_dict[aa_mt]])
        
            log_p_wt.append(np.mean(temp_log_p_wt))
            log_p_mt.append(np.mean(temp_log_p_mt))
            log_p_mt__minus__log_p_wt.append(log_p_mt[-1] - log_p_wt[-1])

    df_out = df.copy()
    df_out['log_p_mt'] = log_p_mt
    df_out['log_p_wt'] = log_p_wt
    df_out['log_p_mt__minus__log_p_wt'] = log_p_mt__minus__log_p_wt


    ## 3) TODO handle multiple mutations

    df_out.to_csv(os.path.join(args.output_dir, 'zero_shot_predictions', args.csv_file.split('/')[-1].replace('.csv', f'-{identifier}.csv')), index=False)


    if args.dms_column is not None:
        args.dms_column = args.dms_column.strip('[ ').strip(' ]')
        from scipy.stats import pearsonr, spearmanr
        non_nan_mask = np.logical_and(np.isfinite(df_out[args.dms_column]), np.isfinite(df_out['log_p_mt__minus__log_p_wt']))
        pearson_r, pearson_pval = pearsonr(df_out[args.dms_column][non_nan_mask], df_out['log_p_mt__minus__log_p_wt'][non_nan_mask])
        spearman_r, spearman_pval = spearmanr(df_out[args.dms_column][non_nan_mask], df_out['log_p_mt__minus__log_p_wt'][non_nan_mask])

        dms_comparison_json = {
            'Pearson R': pearson_r,
            'Pearson R - pval': pearson_pval,
            'Spearman R': spearman_r,
            'Spearman R - pval': spearman_pval
        }

        with open(os.path.join(args.output_dir, 'zero_shot_predictions', args.csv_file.split('/')[-1].replace('.csv', f'-{identifier}.json')), 'w+') as f:
            json.dump(dms_comparison_json, f, indent=2)




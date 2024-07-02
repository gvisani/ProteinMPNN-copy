
'''

NOTE: ignoring insertions codes for the time being

'''

import os, sys
import json
import gzip, pickle
import numpy as np
import pandas as pd
import glob
from tqdm import tqdm
import argparse
from typing import *

from copy import deepcopy

from utils import get_offset_between_proteinmpnn_residue_index_and_PDB_residue_index

try:
    from protein_holography_web.utils.dms_plots import dms_scatter_plot
    CAN_PLOT = True
except:
    CAN_PLOT = False
    pass

alphabet = 'ACDEFGHIKLMNPQRSTVWYX'
alphabet_dict = dict(zip(alphabet, range(21)))


def make_filename(model_version, pdb, chain, resnums):
    return f"{model_version}__{pdb}__{chain}__{','.join([str(x) for x in resnums])}"

def parse_filename(name):
    model_version, pdb, chain, resnums = name[:-4].split('__')
    resnums = [int(x) for x in resnums.split(',')]
    return model_version, pdb, chain, resnums

def get_file_that_matches_specs(inference_dir, model_version, pdb, chain, resnum):
    candidate_files = glob.glob(os.path.join(inference_dir, f"{model_version}__{pdb}__{chain}__*.npz"))
    for file in candidate_files:
        curr_resnums = parse_filename(os.path.basename(file))[3]
        if resnum in curr_resnums:
            return file
    print(f"{model_version}__{pdb}__{chain}__*.npz")
    print(candidate_files)

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

def format_chain_and_resnums_for_proteinmpnn(pdbfile, chains, resnums):

    chain_to_offset_dict = get_offset_between_proteinmpnn_residue_index_and_PDB_residue_index(pdbfile)

    assert len(chains) == len(resnums), f'Length of chains and resnums should be the same, but are {len(chains)} and {len(resnums)} respectively. The chains are {chains} and the resnums are {resnums}.'
    chain_to_resnums_dict = {}
    for chain, resnum in zip(chains, resnums):
        if chain not in chain_to_resnums_dict:
            chain_to_resnums_dict[chain] = []
        chain_to_resnums_dict[chain].append(resnum + chain_to_offset_dict[chain][resnum])
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

    parsed_chains = load_jsonl(path_for_parsed_chains)
    # assigned_chains = load_jsonl(path_for_assigned_chains)[0]

    assert len(parsed_chains) == 1, "Only one PDB file should be processed at a time. It's likely that the temp_pdbs folder has an extra PDB file that should not be there."
    # pdb = list(assigned_chains.keys())[0]
    # chains_list_of_lists = assigned_chains[pdb]
    
    ## the commented-out code below was wrong
    # chains_in_order = []
    # for elem in chains_list_of_lists:
    #     assert len(elem) in {0, 1}, f"These should be singleton lists or empty lists? Not sure why exactly. This might come back to bite me in the ass in the future... just kidding, it's now! {elem}"
    #     if len(elem) == 1:
    #         chains_in_order.append(elem[0])

    # def flatten(xss):
    #     return [x for xs in xss for x in xs]

    # chains_in_order = flatten(chains_list_of_lists)

    chains_in_order = []
    for key in parsed_chains[0]:
        if key.startswith('seq_chain_'):
            chains_in_order.append(key[10:])

    # now sort the mutants by chain as well
    chain_mutant_in_order = []
    for chain in chains_in_order:
        if chain in chain_to_mutants_resnums_tuple:
            for mutant, resnum in chain_to_mutants_resnums_tuple[chain]:
                chain_mutant_in_order.append((chain, mutant))

    return chain_mutant_in_order


def make_prediction(args, inference_dir, pdb, chain, resnum):

    path_for_parsed_chains = os.path.join(args.output_dir, 'parsed_pdbs.jsonl')
    path_for_assigned_chains = os.path.join(args.output_dir, 'assigned_pdbs.jsonl')
    path_for_fixed_positions = os.path.join(args.output_dir, 'fixed_pdbs.jsonl')

    if not os.path.exists(os.path.join(args.folder_with_pdbs, pdb + ".pdb")):
        print(f"Warning: {pdb} not found in {args.folder_with_pdbs}. Skipping...", file=sys.stderr)
        return

    # make temporary directory with single pdb
    temp_pdbs_dir = os.path.join(args.output_dir, "temp_pdbs")
    if os.path.exists(temp_pdbs_dir):
        print(f"Warning: {temp_pdbs_dir} already exists, deleting it...")
        os.system(f'rm -r {temp_pdbs_dir}')
    
    os.makedirs(temp_pdbs_dir)
    os.system(f'cp {os.path.join(args.folder_with_pdbs, pdb + ".pdb")} {temp_pdbs_dir}')

    try:

        # chains = chains_str.split(args.mutant_split_symbol)
        # resnums = [int(resnum) for resnum in resnums_str.split(args.mutant_split_symbol)]
        chains = [chain]
        resnums = [resnum]
        orig_chain = deepcopy(chain)
        orig_resnums = deepcopy(resnums)
        chains_str_for_proteinmpnn, resnums_str_for_proteinmpnn = format_chain_and_resnums_for_proteinmpnn(os.path.join(args.folder_with_pdbs, pdb + ".pdb"), chains, resnums)

        os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/parse_multiple_chains.py")} --input_path={temp_pdbs_dir} --output_path={path_for_parsed_chains}')
        os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/assign_fixed_chains.py")} --input_path={path_for_parsed_chains} --output_path={path_for_assigned_chains} --chain_list "{chains_str_for_proteinmpnn}"')
        os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../helper_scripts/make_fixed_positions_dict.py")} --input_path={path_for_parsed_chains} --output_path={path_for_fixed_positions} --chain_list "{chains_str_for_proteinmpnn}" --position_list "{resnums_str_for_proteinmpnn}" --specify_non_fixed')
        os.system(f'python {os.path.join(os.path.dirname(os.path.realpath(__file__)), f"../protein_mpnn_run.py")} \
                                                                                            --suppress_print 1 \
                                                                                            --model_name {args.model_name} \
                                                                                            --jsonl_path {path_for_parsed_chains} \
                                                                                            --chain_id_jsonl {path_for_assigned_chains} \
                                                                                            --fixed_positions_jsonl {path_for_fixed_positions} \
                                                                                            --out_folder {args.output_dir} \
                                                                                            --num_seq_per_target {args.num_seq_per_target} \
                                                                                            --batch_size {args.batch_size} \
                                                                                            --conditional_probs_only 1')
        

        # parse output file to get the conditional probas
        proteinmpnn_output = os.path.join(args.output_dir, 'conditional_probs_only', pdb + '.npz')
        output = np.load(proteinmpnn_output)

        # chain_mutant_in_order = order_mutants_with_proteinmpnn_output(path_for_parsed_chains, path_for_assigned_chains, mutants, chains)
    
    except Exception as e:
        os.system(f'rm -r {temp_pdbs_dir}') # delete temporary directory with single pdb
        os.system(f'rm {path_for_parsed_chains}')
        os.system(f'rm {path_for_assigned_chains}')
        os.system(f'rm {path_for_fixed_positions}')
        print()
        print(f'Caught the following error, skipping {pdb} {chain} {resnum}:', file=sys.stderr)
        print(e, file=sys.stderr)
        print()
        return


    os.system(f'rm -r {temp_pdbs_dir}') # delete temporary directory with single pdb
    os.system(f'rm {proteinmpnn_output}')
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
        print(f"Warning: site {(chains_str_for_proteinmpnn, resnums_str_for_proteinmpnn)} not found!")
        return

    # wt_aas_being_designed = ''.join([alphabet[index] for index in output['S'][output['design_mask'].astype(np.bool_)]])
    # wt_aas_true = ''.join([m[0] for c, m in chain_mutant_in_order])
    # assert wt_aas_being_designed == wt_aas_true, f"wt_aas_being_designed and wt_aas_true should be the same, but are {wt_aas_being_designed} and {wt_aas_true} respectively."
    # assert log_p.shape[0] == len(chain_mutant_in_order), f"Length of log_p and chain_mutant_in_order should be the same, but are {log_p.shape[0]} and {len(chain_mutant_in_order)} respectively. The mutants are {mutants} and the chains are {chains}. The PDB is {wt_pdb}."

    wt_aas_being_designed = ''.join([alphabet[index] for index in output['S'][output['design_mask'].astype(np.bool_)]])

    os.makedirs(inference_dir, exist_ok=True)
    np.savez(os.path.join(inference_dir, f"{make_filename(args.model_name, pdb, orig_chain, orig_resnums)}.npz"),
                logps=log_p,
                resnums=np.array(orig_resnums), # using the resnum as requested, not the ones for ProteinMPNN! Because this has to match the resnum in the csv file
                wt_aas=np.array(list(wt_aas_being_designed)),
                chain=orig_chain)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--csv_file', type=str, required=True)
    parser.add_argument('--folder_with_pdbs', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)

    parser.add_argument('--dont_run_inference', type=int, default=0)

    # ProteinMPNN arguments
    parser.add_argument('--model_name', type=str, default='v_48_020')
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

    parser.add_argument('--dms_column', type=optional_str, nargs='+', default=None) # ddG
    parser.add_argument('--dms_label', type=optional_str, default=None) #r'stability effect, $\Delta\Delta G$')

    ## when predicting a binary outcome (e.g. Destabilizing vs. Neutral, Bind vs. Not-Bind, ...)
    parser.add_argument('--binary_dms_column', type=optional_str, default=None) #'effect')
    parser.add_argument('--binary_dms_pos_value', type=optional_str, default=None) #'Destabilizing')
    parser.add_argument('--binary_dms_label', type=optional_str, default=None) #'Stability Effect\n(Destabilizing = 1, Neutral = 0)')

    args = parser.parse_args()

    print()
    print('Use MT structure:', type(args.use_mt_structure), args.use_mt_structure)
    print()

    # args.output_dir = args.output_dir + '_v2'

    # prepare output directory
    os.makedirs(args.output_dir, exist_ok=True)
    inference_dir = os.path.join(args.output_dir, 'inference')
    os.makedirs(inference_dir, exist_ok=True)
    zero_shot_predictions_dir = os.path.join(args.output_dir, 'zero_shot_predictions')
    os.makedirs(zero_shot_predictions_dir, exist_ok=True)

    df = pd.read_csv(args.csv_file)

    identifier = f'num_seq_per_target={args.num_seq_per_target}-use_mt_structure={args.use_mt_structure}'

    csv_filename_out = os.path.basename(args.csv_file).split('/')[-1].replace('.csv', f'-{identifier}.csv')
    if not csv_filename_out.endswith('.csv'):
        csv_filename_out += '.csv'
    
    if os.path.exists(os.path.join(zero_shot_predictions_dir, csv_filename_out.replace('.csv', '_correlations.json'))):
        # exit with OK exit code
        exit(0)


    # get pdbs
    pdbs = set()
    if args.use_mt_structure:
        for pdb in df[args.mt_pdb_column]:
            pdbs.add(pdb)
    for pdb in df[args.wt_pdb_column]:
        pdbs.add(pdb)
    pdbs = list(pdbs)

    # download necessary if they are not found in folder
    os.makedirs(args.folder_with_pdbs, exist_ok=True) # make folder if does not exist, so it doesn't have to be made beforehand if it's empty
    pdbs_in_folder = [file[:-4] for file in os.listdir(args.folder_with_pdbs)]
    # print(pdbs_in_folder)
    # pdbs_to_download = set(pdbs) - set(pdbs_in_folder)
    # if len(pdbs_to_download) > 0:
    #     print(f'Downloading the following PDBs: {pdbs_to_download}')
    #     for pdb in pdbs_to_download:
    #         try:
    #             urllib.request.urlretrieve(f'http://files.rcsb.org/download/{pdb}.pdb', os.path.join(args.folder_with_pdbs, f'{pdb}.pdb'))
    #         except Exception as e:
    #             print(f'Error downloading {pdb}: {e}')
    #             continue
    
    # get chains and resnums for each PDB (single call per individual chain)
    pdb_to_chain_to_resnums = {}
    for i, row in df.iterrows():
        if args.use_mt_structure:
            row_pdbs = [row[args.wt_pdb_column], row[args.mt_pdb_column]]
        else:
            row_pdbs = [row[args.wt_pdb_column]]
        
        if not isinstance(row[args.mutant_chain_column], str) or not isinstance(row[args.mutant_column], str):
            # sometimes the chain and mutant columns are NaN, in which case we skip the row
            continue
        
        chains = row[args.mutant_chain_column].split(args.mutant_split_symbol)
        mutants = row[args.mutant_column].split(args.mutant_split_symbol)
        assert len(chains) == len(mutants)

        for chain, mutant in zip(chains, mutants):
            resnum = int(mutant[1:-1])
            for pdb in row_pdbs:
                if not isinstance(pdb, str): # excluded nan pdbs
                    continue
                if pdb not in pdb_to_chain_to_resnums:
                    pdb_to_chain_to_resnums[pdb] = {}
                if chain not in pdb_to_chain_to_resnums[pdb]:
                    pdb_to_chain_to_resnums[pdb][chain] = set()
                
                pdb_to_chain_to_resnums[pdb][chain].add(resnum)
        
        # chains_str = row[args.mutant_chain_column]

        # chains = chains_str.split(args.mutant_split_symbol)
        # mutants = row[args.mutant_column].split(args.mutant_split_symbol)
        # assert len(chains) == len(mutants)

        # resnums = [int(mutant[1:-1]) for mutant in mutants]
        # resnums_str = args.mutant_split_symbol.join([str(resnum) for resnum in resnums])

        # for pdb in row_pdbs:
        #     if pdb not in pdb_to_chain_to_resnums:
        #         pdb_to_chain_to_resnums[pdb] = {}
        #     if chains_str not in pdb_to_chain_to_resnums[pdb]:
        #         pdb_to_chain_to_resnums[pdb][chains_str] = set()
            
        #     pdb_to_chain_to_resnums[pdb][chains_str].add(resnums_str)

    
    print(pdb_to_chain_to_resnums)

    tuples = []
    for pdb in pdb_to_chain_to_resnums:
        for chain, resnums_set in pdb_to_chain_to_resnums[pdb].items():
            resnums_list = list(sorted(list(resnums_set)))
            for resnum in resnums_list:
                tuples.append((pdb, chain, resnum))

    ## run inference!!
    from time import time
    start = time()
    if not args.dont_run_inference:
        for pdb, chain, resnum in tqdm(tuples):
            # print(f'Running inference for {pdb} {chain} {resnum}')
            make_prediction(args, inference_dir, pdb, chain, resnum)
    end = time()
    print(f'Inference took {end - start} seconds')
    

    # ## run inference!!
    # from time import time
    # start = time()
    # if not args.dont_run_inference:
    #     for pdb in tqdm(pdb_to_chain_to_resnums):
    #         for chains_str, resnums_set in pdb_to_chain_to_resnums[pdb].items():
    #             ## split resnums into chunks of at most 20 otherwise I might get "File name too long" error
    #             resnums_list = list(sorted(list(resnums_set)))
    #             for resnums_str in resnums_list:
    #                 print(f'Running inference for {pdb} {chains_str} {resnums_list}')
    #                 make_prediction(args, inference_dir, pdb, chains_str, resnums_str)
    # end = time()
    # print(f'Inference took {end - start} seconds')

    print()
    print('Use MT structure:', type(args.use_mt_structure), args.use_mt_structure)
    print()

    ## parse the .npz files and save the results as a new column
    log_proba_wt_all = []
    log_proba_mt_all = []
    print('Parsing .npz files...')
    for i, row in tqdm(df.iterrows(), total=len(df)):

        wt_pdb = row[args.wt_pdb_column]

        if not isinstance(row[args.mutant_chain_column], str) or not isinstance(row[args.mutant_column], str):
            # sometimes the chain and mutant columns are NaN, in which case we skip the row
            print(f'WARNING: chain ({row[args.mutant_chain_column]}) and/or mutant ({row[args.mutant_column]}) is not a string at row {i}.', file=sys.stderr)
            log_proba_wt_all.append(np.nan)
            log_proba_mt_all.append(np.nan)
            continue

        chains = row[args.mutant_chain_column].split(args.mutant_split_symbol)
        mutants = row[args.mutant_column].split(args.mutant_split_symbol)
        assert len(chains) == len(mutants)

        ## average results across multiple mutations
        temp_log_proba_wt = []
        temp_log_proba_mt = []
        for chain, mutant in zip(chains, mutants):
            aa_wt = mutant[0]
            aa_mt = mutant[-1]
            resnum = int(mutant[1:-1])

            if aa_mt == 'X':
                print(f'WARNING: Mutant amino-acid is X at {wt_pdb} {chain} {resnum}.', file=sys.stderr)
                temp_log_proba_wt.append(np.nan)
                temp_log_proba_mt.append(np.nan)
                continue

            wt_file = get_file_that_matches_specs(inference_dir, args.model_name, wt_pdb, chain, resnum)
            if wt_file is None:
                print(f'WARNING: No file found for {wt_pdb} {chain} {resnum}.', file=sys.stderr)
                temp_log_proba_wt.append(np.nan)
                temp_log_proba_mt.append(np.nan)
                continue

            wt_data = np.load(wt_file)

            if args.use_mt_structure:
                mt_pdb = row[args.mt_pdb_column]
                mt_file = get_file_that_matches_specs(inference_dir, args.model_name, mt_pdb, chain, resnum)
                if mt_file is None:
                    print(f'WARNING: No file found for {mt_pdb} {chain} {resnum}.', file=sys.stderr)
                    temp_log_proba_wt.append(np.nan)
                    temp_log_proba_mt.append(np.nan)
                    continue
                mt_data = np.load(mt_file)
            else:
                mt_pdb = wt_pdb
                mt_data = wt_data
            
            if wt_data['logps'].shape[0] != wt_data['resnums'].shape[0]:
                print(f"WARNING: Shape mismatch in wt_data. pes is {wt_data['logps'].shape[0]}, resnums is {wt_data['resnums'].shape[0]}. Skipping {wt_pdb} {chain} {resnum}.", file=sys.stderr)
                temp_log_proba_wt.append(np.nan)
                temp_log_proba_mt.append(np.nan)
                continue
            
            try:
                np.where(wt_data['resnums'] == resnum)[0][0]
            except IndexError:
                print(f"WARNING: Resnum {resnum} not found in wt_data. Skipping {wt_pdb} {chain} {resnum}.", file=sys.stderr)
                temp_log_proba_wt.append(np.nan)
                temp_log_proba_mt.append(np.nan)
                continue
                
            # check that the wildtype amino-acids as they are in the csv file match the amino-acids in the structure
            aa_wt_in_structure = wt_data['wt_aas'][np.where(wt_data['resnums'] == resnum)[0][0]]
            # assert aa_wt_in_structure == aa_wt, f'Wildtype residue mismatch! {aa_wt_in_structure} != {aa_wt} at {wt_pdb} {chain} {resnum}'
            if aa_wt_in_structure != aa_wt:
                print(f'\nWARNING: Wildtype residue mismatch! {aa_wt_in_structure} != {aa_wt} at {wt_pdb} {chain} {resnum}\n', file=sys.stderr)
                temp_log_proba_wt.append(np.nan)
                temp_log_proba_mt.append(np.nan)
                continue

            if args.use_mt_structure:
                aa_mt_in_structure = mt_data['wt_aas'][np.where(mt_data['resnums'] == resnum)[0][0]]
                # assert aa_mt_in_structure == aa_mt, f'Mutant residue mismatch! {aa_mt_in_structure} != {aa_mt} at {wt_pdb} {chain} {resnum}'
                if aa_mt_in_structure != aa_mt:
                    print(f'\nWARNING: Mutant residue mismatch! {aa_mt_in_structure} != {aa_mt} at {mt_pdb} {chain} {resnum}\n', file=sys.stderr)
                    temp_log_proba_wt.append(np.nan)
                    temp_log_proba_mt.append(np.nan)
                    continue
            

            wt_logp = wt_data['logps'][np.where(wt_data['resnums'] == resnum)[0][0], alphabet_dict[aa_wt]]
            mt_logp = mt_data['logps'][np.where(mt_data['resnums'] == resnum)[0][0], alphabet_dict[aa_mt]]

            temp_log_proba_wt.append(wt_logp)
            temp_log_proba_mt.append(mt_logp)
        
        log_proba_wt_all.append(np.mean(temp_log_proba_wt))
        log_proba_mt_all.append(np.mean(temp_log_proba_mt))

    df['log_p_wt'] = log_proba_wt_all
    df['log_p_mt'] = log_proba_mt_all
    df['log_p_mt__minus__log_p_wt'] = df['log_p_mt'] - df['log_p_wt']


    df.to_csv(os.path.join(zero_shot_predictions_dir, csv_filename_out), index=False)


    if args.dms_column is not None and CAN_PLOT:
        for dms_column in args.dms_column:
            dms_column = dms_column.strip('[ ').strip(' ]')
            dms_filepath = os.path.join(zero_shot_predictions_dir, f'correlation_with_{dms_column}-{csv_filename_out.replace(".csv", ".png")}')
            (pearson_r, pearson_pval), (spearman_r, spearman_pval) = dms_scatter_plot(df,
                                                                                  dms_column, 'log_p_mt__minus__log_p_wt',
                                                                                  dms_label=None, pred_label=r'ProteinMPNN Prediction',
                                                                                  filename = dms_filepath)
    
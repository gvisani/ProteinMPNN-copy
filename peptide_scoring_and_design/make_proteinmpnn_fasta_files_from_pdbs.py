
'''

Leverage the script `parse_multiple_chains.py` to make fasta files containing the amino-acid sequence within PDBs.
We are using that script because that is the script that ProteinMPNN uses to parse PDB files, so we can be sure that the sequence variants we put in are going to be compatible.

'''

import os
import json
import argparse

alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder_with_pdbs', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    args = parser.parse_args()

    filepath_with_parsed_chains = './parsed_pdbs.jsonl'

    os.system(f"python ../helper_scripts/parse_multiple_chains.py --input_path {args.folder_with_pdbs} --output_path {filepath_with_parsed_chains}") # will generate a file called "parsed_pdbs.jsonl"

    
    # now open the "parsed_pdbs.jsonl" file and make a fasta file for each pdb that has been parsed

    with open(filepath_with_parsed_chains, 'r') as f:
        json_list = list(f)
    
    for json_str in json_list:
        parsed = json.loads(json_str)
        sequence = []
        for chain_i in range(parsed['num_of_chains']):
            chain = alphabet[chain_i]
            chain_seq = parsed[f'seq_chain_{chain}']
            sequence.append(chain_seq.replace('-', 'X'))
        
        pdbid = parsed['name']
        sequence = '/'.join(sequence)

        with open(os.path.join(args.output_dir, pdbid + '.fa'), 'w+') as f:
            f.write(f'>{pdbid}\n')
            f.write(sequence)
    
    os.system(f'rm {filepath_with_parsed_chains}')







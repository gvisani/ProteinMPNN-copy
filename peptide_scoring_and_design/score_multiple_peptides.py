
import os
import numpy as np
import pandas as pd
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--csv_file', type=str, required=True)
    parser.add_argument('--folder_with_pdbs', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)

    parser.add_argument('--pdb_column', type=str, default='pdbs')
    parser.add_argument('--chain_column', type=str, default='chains')
    parser.add_argument('--peptides_column', type=str, default='sequence')
    parser.add_argument('--affinity_column', type=str, default='-log10(Kd)', help='Column in the csv file that we want to correlate with. Usually an affinity score.')

    parser.add_argument('--collect_scores_only', type=int, default=0, choices=[0, 1], help='If True, then only collect scores from the output directory and do not score the peptides again.')

    # ProteinMPNN arguments
    parser.add_argument('--num_seq_per_target', type=int, default=10)
    parser.add_argument('--batch_size', type=int, default=10)

    args = parser.parse_args()

    df = pd.read_csv(args.csv_file)

    identifier = f'num_seq_per_target={args.num_seq_per_target}'

    pdbs_list = df[args.pdb_column].tolist()
    chains_list = df[args.chain_column].tolist()
    peptides = df[args.peptides_column].tolist()

    if not args.collect_scores_only:

        # score!
        for pdbs_str, chains_str, peptide in zip(pdbs_list, chains_list, peptides):

            pdbs = pdbs_str.split('|')
            chains = chains_str.split('|')
            assert len(pdbs) == len(chains)
            assert len(set(chains)) == 1, "All chains should be the same because it's most convenient with how ProteinMPNN works."
            peptide_chain = chains[0]

            # put the requested pdbs for the peptide of interest in a temporary folder
            temp_folder = os.path.join(args.output_dir, 'temp_pdbs')
            os.makedirs(temp_folder, exist_ok=True)
            for pdb in pdbs:
                os.system(f"cp {os.path.join(args.folder_with_pdbs, pdb+'.pdb')} {temp_folder}")
            
            print(os.listdir(temp_folder))

            os.system(f'python score_whole_peptide.py \
                                        --folder_with_pdbs {temp_folder} \
                                        --output_dir {args.output_dir} \
                                        --peptide_seq {peptide} \
                                        --peptide_chain {peptide_chain} \
                                        --num_seq_per_target {args.num_seq_per_target} \
                                        --batch_size {args.batch_size}')
            
            # remove the temporary pdbs
            os.system(f'rm -rf {temp_folder}')
    
    # collect scores
    all_score_files = os.listdir(os.path.join(args.output_dir, 'score_only'))
    scores = []
    global_scores = []
    all_scores_as_str = []
    all_global_scores_as_str = []
    for pdbs_str, chains_str, peptide in zip(pdbs_list, chains_list, peptides):
        pdbs = pdbs_str.split('|')
        peptide_scores_list = []
        global_scores_list = []
        score_files = list(filter(lambda filename: filename.split('_')[1] == peptide, all_score_files))
        for score_file in score_files:
            # keep only scores with the right pdbs
            curr_pdb = score_file.split('_')[0]
            if curr_pdb not in pdbs:
                continue
            score_data = np.load(os.path.join(args.output_dir, 'score_only', score_file))
            peptide_scores_list.append(score_data['score'])
            global_scores_list.append(score_data['global_score'])

        scores.append(-np.mean(np.hstack(peptide_scores_list))) # score is mean -log(p)
        global_scores.append(-np.mean(np.hstack(global_scores_list))) # score is mean -log(p)
        all_scores_as_str.append('|'.join([str(np.mean(score)) for score in peptide_scores_list]))
        all_global_scores_as_str.append('|'.join([str(np.mean(score)) for score in global_scores_list]))
    
    # save scores in new columns in a new csv file
    df_out = df.copy()
    df_out['-log(p_pep)'] = scores
    df_out['-log(p_global)'] = global_scores
    df_out['all -log(p_pep)'] = all_scores_as_str
    df_out['all -log(p_global)'] = all_global_scores_as_str
    os.makedirs(os.path.join(args.output_dir, 'zero_shot_predictions'), exist_ok=True)
    df_out.to_csv(os.path.join(args.output_dir, 'zero_shot_predictions', args.csv_file.split('/')[-1].replace('.csv', f'-{identifier}.csv')), index=False)


    # produce scatterplot and compute correlations! - p_pep
    from scipy.stats import pearsonr, spearmanr
    pearson = pearsonr(scores, df_out[args.affinity_column])
    spearman = spearmanr(scores, df_out[args.affinity_column])
    print('----- p_pep -----')
    print('Pearson r:\t%.3f\t(%.5f)' % (pearson[0], pearson[1]))
    print('Spearmanr r:\t%.3f\t(%.5f)' % (spearman[0], spearman[1]))

    import matplotlib.pyplot as plt
    plt.figure(figsize=(4, 3))
    plt.scatter(df_out[args.affinity_column], scores, label='Pr: %.3f\nSr: %.3f' % (pearson[0], spearman[0]))
    plt.xlabel(args.affinity_column)
    plt.ylabel('ProteinMPNN -log(p_pep)')
    plt.legend()
    plt.tight_layout()
    os.makedirs(os.path.join(args.output_dir, 'zero_shot_predictions'), exist_ok=True)
    plt.savefig(os.path.join(args.output_dir, 'zero_shot_predictions', 'scatterplot_' + args.csv_file.split('/')[-1].replace('.csv', f'-{identifier}-with_-log(p_pep).png')))


    # produce scatterplot and compute correlations! - p_global
    from scipy.stats import pearsonr, spearmanr
    pearson = pearsonr(global_scores, df_out[args.affinity_column])
    spearman = spearmanr(global_scores, df_out[args.affinity_column])
    print('----- p_global -----')
    print('Pearson r:\t%.3f\t(%.5f)' % (pearson[0], pearson[1]))
    print('Spearmanr r:\t%.3f\t(%.5f)' % (spearman[0], spearman[1]))

    import matplotlib.pyplot as plt
    plt.figure(figsize=(4, 3))
    plt.scatter(df_out[args.affinity_column], global_scores, label='Pr: %.3f\nSr: %.3f' % (pearson[0], spearman[0]))
    plt.xlabel(args.affinity_column)
    plt.ylabel('ProteinMPNN -log(p_global)')
    plt.legend()
    plt.tight_layout()
    os.makedirs(os.path.join(args.output_dir, 'zero_shot_predictions'), exist_ok=True)
    plt.savefig(os.path.join(args.output_dir, 'zero_shot_predictions', 'scatterplot_' + args.csv_file.split('/')[-1].replace('.csv', f'-{identifier}-with_-log(p_global).png')))









import os
import numpy as np
import pandas as pd
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--csv_file', type=str, required=True)
    parser.add_argument('--folder_with_pdbs', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--peptide_chain', type=str, default='C')

    parser.add_argument('--peptides_column', type=str, default='sequence', help='Column in the csv file to correlate with')
    parser.add_argument('--affinity_column', type=str, default='-log10(Kd)', help='Column in the csv file that we want to correlate with. Usually an affinity score.')

    # ProteinMPNN arguments
    parser.add_argument('--num_seq_per_target', type=int, default=100)
    parser.add_argument('--batch_size', type=int, default=25)

    args = parser.parse_args()

    df = pd.read_csv(args.csv_file)

    # get peptides
    peptides = df[args.peptides_column].tolist()

    # score!
    for peptide in peptides:
        os.system(f'python score_whole_peptide.py \
                                    --folder_with_pdbs {args.folder_with_pdbs} \
                                    --output_dir {args.output_dir} \
                                    --peptide_seq {peptide} \
                                    --peptide_chain {args.peptide_chain} \
                                    --num_seq_per_target {args.num_seq_per_target} \
                                    --batch_size {args.batch_size}')
    
    # collect scores
    all_score_files = os.listdir(os.path.join(args.output_dir, 'score_only'))
    scores = []
    global_scores = []
    for peptide in peptides:
        peptide_scores_list = []
        global_scores_list = []
        score_files = list(filter(lambda filename: filename.split('_')[1] == peptide, all_score_files))
        for score_file in score_files:
            score_data = np.load(os.path.join(args.output_dir, 'score_only', score_file))
            peptide_scores_list.append(score_data['score'])
            global_scores_list.append(score_data['global_score'])
        scores.append(-np.mean(np.hstack(peptide_scores_list))) # score is mean -log(p)
        global_scores.append(-np.mean(np.hstack(global_scores_list))) # score is mean -log(p)
    
    # save scores in new columns in the same csv file
    df['ProteinMPNN_-log(p_pep)'] = scores
    df['ProteinMPNN_-log(p_global)'] = global_scores
    df.to_csv(args.csv_file, index=False)

    # produce scatterplot and compute correlations! - p_pep
    from scipy.stats import pearsonr, spearmanr
    pearson = pearsonr(scores, df[args.affinity_column])
    spearman = spearmanr(scores, df[args.affinity_column])
    print('----- p_pep -----')
    print('Pearson r:\t%.3f\t(%.5f)' % (pearson[0], pearson[1]))
    print('Spearmanr r:\t%.3f\t(%.5f)' % (spearman[0], spearman[1]))

    import matplotlib.pyplot as plt
    plt.figure(figsize=(4, 3))
    plt.scatter(df[args.affinity_column], scores, label='Pr: %.3f\nSr: %.3f' % (pearson[0], spearman[0]))
    plt.xlabel(args.affinity_column)
    plt.ylabel('ProteinMPNN -log(p_pep)')
    plt.legend()
    plt.tight_layout()
    os.makedirs(os.path.join(args.output_dir, 'binding_affinity_plots'), exist_ok=True)
    plt.savefig(os.path.join(args.output_dir, 'binding_affinity_plots', 'scatterplot_with_-log(p_pep).png'))


    # produce scatterplot and compute correlations! - p_global
    from scipy.stats import pearsonr, spearmanr
    pearson = pearsonr(global_scores, df[args.affinity_column])
    spearman = spearmanr(global_scores, df[args.affinity_column])
    print('----- p_global -----')
    print('Pearson r:\t%.3f\t(%.5f)' % (pearson[0], pearson[1]))
    print('Spearmanr r:\t%.3f\t(%.5f)' % (spearman[0], spearman[1]))

    import matplotlib.pyplot as plt
    plt.figure(figsize=(4, 3))
    plt.scatter(df[args.affinity_column], global_scores, label='Pr: %.3f\nSr: %.3f' % (pearson[0], spearman[0]))
    plt.xlabel(args.affinity_column)
    plt.ylabel('ProteinMPNN -log(p_global)')
    plt.legend()
    plt.tight_layout()
    os.makedirs(os.path.join(args.output_dir, 'binding_affinity_plots'), exist_ok=True)
    plt.savefig(os.path.join(args.output_dir, 'binding_affinity_plots', 'scatterplot_with_-log(p_global).png'))








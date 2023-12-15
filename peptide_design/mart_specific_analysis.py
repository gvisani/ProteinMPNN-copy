


import os
import numpy as np
import pandas as pd
import argparse

def get_scores(args, peptides, pdbs):
    all_score_files = os.listdir(os.path.join(args.output_dir, 'score_only'))
    scores = []
    global_scores = []
    for peptide in peptides:
        peptide_scores_list = []
        global_scores_list = []
        score_files = list(filter(lambda filename: filename.split('_')[1] == peptide and filename.split('_')[0] in pdbs, all_score_files))
        for score_file in score_files:
            score_data = np.load(os.path.join(args.output_dir, 'score_only', score_file))
            peptide_scores_list.append(score_data['score'])
            global_scores_list.append(score_data['global_score'])
        scores.append(-np.mean(np.hstack(peptide_scores_list))) # score is mean -log(p)
        global_scores.append(-np.mean(np.hstack(global_scores_list))) # score is mean -log(p)
    return scores, global_scores

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--csv_file', type=str, default='./binding_affinity/mart/data.csv')
    parser.add_argument('--output_dir', type=str, default='./outputs/mart')

    parser.add_argument('--peptides_column', type=str, default='sequence', help='Column in the csv file to correlate with')
    parser.add_argument('--affinity_column', type=str, default='-log10(Kd)', help='Column in the csv file that we want to correlate with. Usually an affinity score.')

    args = parser.parse_args()

    df = pd.read_csv(args.csv_file)

    # get peptides
    peptides = df[args.peptides_column].tolist()

    ## gather which pdbs and peptides are associated with which motif
    GIG_pdbs = None
    DRG_pdbs = None
    GIG_peptides = []
    DRG_peptides = []
    for i, row in df.iterrows():
        if row['motif'] == 'GIG':
            GIG_pdbs = row['structures'].split('_')
            GIG_peptides.append(row['sequence'])
        elif row['motif'] == 'DRG':
            DRG_pdbs = row['structures'].split('_')
            DRG_peptides.append(row['sequence'])
    
    ## scores with GIG structures only
    GIG_pep_scores, GIG_global_scores = get_scores(args, peptides, GIG_pdbs)
    df['ProteinMPNN_-log(p_pep)_GIG'] = GIG_pep_scores
    df['ProteinMPNN_-log(p_global)_GIG'] = GIG_global_scores

    ## scores with DRG structures only
    DRG_pep_scores, DRG_global_scores = get_scores(args, peptides, DRG_pdbs)
    df['ProteinMPNN_-log(p_pep)_DRG'] = DRG_pep_scores
    df['ProteinMPNN_-log(p_global)_DRG'] = DRG_global_scores

    ## scores with respective structures - just need to slice the lists
    pep_scores_respective_pdbs = []
    global_scores_respective_pdbs = []
    for i, row in df.iterrows():
        if row['motif'] == 'GIG':
            pep_scores_respective_pdbs.append(GIG_pep_scores[i])
            global_scores_respective_pdbs.append(GIG_global_scores[i])
        elif row['motif'] == 'DRG':
            pep_scores_respective_pdbs.append(DRG_pep_scores[i])
            global_scores_respective_pdbs.append(DRG_global_scores[i])
    df['ProteinMPNN_-log(p_pep)_respective'] = pep_scores_respective_pdbs
    df['ProteinMPNN_-log(p_global)_respective'] = global_scores_respective_pdbs

    df.to_csv(args.csv_file, index=False)

    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr, spearmanr

    for p_str in ['p_pep', 'p_global']:

        ncols = 4
        nrows = 1
        fig, axs = plt.subplots(figsize=(4*ncols, 4*nrows), ncols=ncols, nrows=nrows)

        for col, (ax, scores, title) in enumerate(zip(axs,
                                                      [df[f'ProteinMPNN_-log({p_str})_GIG'], df[f'ProteinMPNN_-log({p_str})_DRG'], df[f'ProteinMPNN_-log({p_str})_respective'], df[f'ProteinMPNN_-log({p_str})']],
                                                      ['GIG template', 'DRG template', 'Respective struc.', 'All struc.'])):

            pearson = pearsonr(scores, df[args.affinity_column])
            spearman = spearmanr(scores, df[args.affinity_column])

            GIG_mask = df['motif'] == 'GIG'
            DRG_mask = df['motif'] == 'DRG'
            ax.scatter(df[args.affinity_column][GIG_mask], scores[GIG_mask], color='tab:cyan', label='Pr: %.3f\nSr: %.3f' % (pearson[0], spearman[0]))
            ax.scatter(df[args.affinity_column][DRG_mask], scores[DRG_mask], color='tab:pink')

            ax.set_xlabel(args.affinity_column)
            if col == 0: ax.set_ylabel(f'ProteinMPNN -log({p_str})')
            ax.set_title(title)
            ax.legend()
        
        plt.tight_layout()
        os.makedirs(os.path.join(args.output_dir, 'binding_affinity_plots'), exist_ok=True)
        plt.savefig(os.path.join(args.output_dir, 'binding_affinity_plots', f'scatterplots_broken_down_by_motif_with_-log({p_str}).png'))


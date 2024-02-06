

base_dir='/gscratch/spe/gvisan01/mutation_effect_predictions/full_pep/'

# system_and_csv_files=(tax#tax_peptide_kd.csv nyeso#nyeso_peptide_kd.csv mart#mart_peptide_kd_ALL.csv mart#mart_peptide_kd_GIG.csv mart#mart_peptide_kd_DRG.csv mart#mart_peptide_kd_RESPECTIVE.csv)
system_and_csv_files=(mart#mart_peptide_kd_GIG.csv mart#mart_peptide_kd_DRG.csv mart#mart_peptide_kd_RESPECTIVE.csv)

folder_with_pdbs='/gscratch/spe/gvisan01/tcr_pmhc/pdbs/'

IFS=$IFS# # to use '#' as a delimiter

for system_and_csv_file in "${system_and_csv_files[@]}"
    do
    echo $system_and_csv_file

    # get system and csv file, the bash code is rather convoluted
    temp=()
    read -a temp <<< "$system_and_csv_file"
    system="${temp[0]}"
    csv_file="${temp[1]}"

    echo $system
    echo $csv_file

    python -u score_multiple_peptides.py \
                        --csv_file $base_dir$system'/'$csv_file \
                        --folder_with_pdbs $folder_with_pdbs \
                        --output_dir $base_dir$system'/proteinmpnn' \
                        --pdb_column 'pdbs' \
                        --chain_column 'chains' \
                        --peptides_column 'sequence' \
                        --affinity_column'=-log10(Kd)' \
                        --num_seq_per_target 10 \
                        --batch_size 10 \
                        --collect_scores_only 1

done

IFS="${IFS:0:3}" # reset to default values after loop


# num_seq_per_target_list='100 10 1'
# use_mt_structure_list='0 1'
# systems='single_point_mutation__A6__HLA-Astar02column01__LLFGYPVYV single_point_mutation__DMF5__HLA-Astar02column01__AAGIGILTV single_point_mutation__JM22__HLA-Astar02column01__GILGFVFTL single_point_mutation__LC13__HLA-Bstar08column01__FLRGRAYGL'


num_seq_per_target_list='10'
use_mt_structure_list='0 1'

systems='single_point_mutation__A6__HLA-Astar02column01__LLFGYPVYV'

pdb_dir='/gscratch/spe/gvisan01/tcr_pmhc/pdbs/ATLAS/'

for system in $systems
    do

    base_dir='/gscratch/spe/gvisan01/mutation_effect_predictions/ATLAS_TCR/'$system'/'

    for num_seq_per_target in $num_seq_per_target_list
        do
        for use_mt_structure in $use_mt_structure_list
            do

            python -u ../zero_shot_mutation_effect_prediction.py \
                            --csv_file $base_dir$system'.csv' \
                            --folder_with_pdbs $pdb_dir \
                            --output_dir $base_dir'proteinmpnn' \
                            --use_mt_structure $use_mt_structure \
                            --num_seq_per_target $num_seq_per_target \
                            --batch_size 10 \
                            --dms_column '[ log10(Kd_wt/Kd_mt) ]' \
                            --mutant_column TCR_mut \
                            --mutant_chain_column TCR_PDB_chain \

        done
    done
done

# ## tax
# python kd_correlation_analysis.py --csv_file ./binding_affinity/tax/data.csv --folder_with_pdbs ./pdbs/tax --output_dir ./outputs/tax

# ## ny-eso
# python kd_correlation_analysis.py --csv_file ./binding_affinity/nyeso/data.csv --folder_with_pdbs ./pdbs/nyeso --output_dir ./outputs/nyeso

## mart
python kd_correlation_analysis.py --csv_file ./binding_affinity/mart/data.csv --folder_with_pdbs ./pdbs/mart --output_dir ./outputs/mart
python mart_specific_analysis.py --csv_file ./binding_affinity/mart/data.csv --output_dir ./outputs/mart
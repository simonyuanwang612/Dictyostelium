import os
import sys
import argparse
import pandas as pd

# usage statement and input descriptions
parser = argparse.ArgumentParser(
    description='A post-processing script to merge files from Skyline\'s Group Comparison report format.\ '
                'This script will return three outputs: a merged file with all the combined data, a file containing\ '
                'just the differential testing results, and a file containing just the abundance values.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('annotation_file', type=str,
                    help='a file with each of the Group Comparison filenames and their annotation')
parser.add_argument('acquisition_file', type=str,
                    help='a file with each of the individual acquisition filenames and their annotation')
parser.add_argument('--output_path', default=os.getcwd(), type=str,
                    help='specify an output path for the merged result file')


# parse arguments from command line
args = parser.parse_args()
annotation_file = args.annotation_file
acquisition_file = args.acquisition_file
output_path = args.output_path

cwd = os.getcwd()

# iterate over the filenames, add a sample group column, and merge them together
annotation_df = pd.read_csv(annotation_file)
files2merge = annotation_df['filename']
merge_df = pd.DataFrame()
for file in files2merge:
    temp_df = pd.read_csv(os.path.join(cwd, file))
    annotation = annotation_df[annotation_df['filename'] == file]['annotation'].iloc[0]
    temp_df['Group Comparison'] = [annotation] * len(temp_df)
    merge_df = merge_df.append(temp_df)



# split into a differential testing dataframe and an abundance dataframe
allcols = list(merge_df.columns)
abundance_cols = [s for s in allcols if "Abundance" in s]
difftest_df = merge_df.drop(abundance_cols, axis=1)

acquisition_df = pd.read_csv(acquisition_file)
difftest_cols = ['Fold Change Result', 'Adjusted P-Value', 'Group Comparison']
abundance_df = pd.melt(merge_df.drop(difftest_cols, axis=1),
                       id_vars=['Protein', 'Peptide', 'Peptide Modified Sequence', 'MS Level', 'histone mark'],
                       value_vars=abundance_cols,
                       var_name='Acquisition',
                       value_name='Normalized Abundance')
abundance_df = pd.merge(abundance_df, acquisition_df, how='inner', on='Acquisition')

# drop the resulting NaN from melting the abundance dataframe
abundance_df = abundance_df.dropna(subset=['Normalized Abundance'])
abundance_df = abundance_df.drop_duplicates()

# write a "percentages" dataframe
percentage_df = abundance_df
temp = percentage_df.groupby(['Protein', 'Peptide', 'Acquisition'])['Normalized Abundance'].agg('sum').reset_index()
temp.rename(columns={'Normalized Abundance': 'Peptide Sum Abnd'}, inplace=True)

percentage_df = pd.merge(percentage_df, temp, how='left', on=['Protein', 'Peptide', 'Acquisition'])
percentage_df['Percentage Histone Mark'] = (percentage_df['Normalized Abundance']/percentage_df['Peptide Sum Abnd']) * 100


# write out the parsed results
merge_df.to_csv(path_or_buf=os.path.join(output_path, 'merged_skyline_groupcomparisons.csv'),
                index=False)
difftest_df.to_csv(path_or_buf=os.path.join(output_path, 'merged_skyline_groupcomparisons_difftest.csv'),
                index=False)
abundance_df.to_csv(path_or_buf=os.path.join(output_path, 'merged_skyline_groupcomparisons_abundances.csv'),
                index=False)
percentage_df.to_csv(path_or_buf=os.path.join(output_path, 'merged_skyline_groupcomparisons_percentages.csv'),
                index=False)

sys.stdout.write("Finished merging Skyline Group Comparison exports. Wrote out three parsed result dataframes.\n")
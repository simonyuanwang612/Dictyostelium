import os
import sys
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

sns.set_style("whitegrid")


# usage statement and input descriptions
parser = argparse.ArgumentParser(
    description='A post-processing script to make bar plots from the merge_groupcomparisons.py abundance output.\ '
                'This script will build, for each peptide, a ',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('abundances_file', type=str,
                    help='the merged_skyline_groupcomparisons_abundances.csv file that is output from running \
                         merge_groupcomparisons.py.')
parser.add_argument('--output_path', default=os.getcwd(), type=str,
                    help='specify an output path for the bar plots')

# parse arguments from command line
args = parser.parse_args()
abundances_file = args.abundances_file
output_dir = args.output_path

# read in group comparison results from skyline
quant_df = pd.read_csv(os.path.join(os.getcwd(), abundances_file))
#quant_df['Total Area Fragment'] = quant_df['Total Area Fragment'].replace(np.nan, 0)

# for each modified peptide
sys.stderr.write("Building plots for each modified peptide...")
for pepmodseq in tqdm(quant_df['Peptide Modified Sequence'].drop_duplicates()):

    temp_df = quant_df[quant_df['Peptide Modified Sequence'] == pepmodseq]
    temp_df = temp_df.drop_duplicates()
    # print(test_df)

    protein = temp_df['Protein'].iloc[0]

    # prettify the protein name for this peptide
    # TODO use SeqIO to parse protein names in a less stupid way
    new_protein = ''
    if "|" in protein:
        if "," in protein:
            protein_list = protein.split(',')
            new_protein = []
            for prot in protein_list:
                temp = str(prot.split("|")[-1])
                new_protein.append(temp)
            new_protein = ', '.join(new_protein)
        else:
            new_protein = protein.split("|")[-1]
    protein = new_protein

    # TODO remove this hardcoded DICDI cleanup
    protein = protein.replace('_DICDI', '')

    # get the name of this histone mark and add it to the protein name
    if temp_df['histone mark'].any():
        his_mark = temp_df['histone mark'].iloc[0]
        his_mark = his_mark.replace('[', '')
        his_mark = his_mark.replace(']', '')
        title = protein + ' ' + his_mark
    else:
        title = protein + ' ' + pepmodseq

    # make a barchart
    sns.barplot(x='Sample Group', y='Normalized Abundance', data=temp_df,
                capsize=.1, ci="sd")
    sns.swarmplot(x='Sample Group', y='Normalized Abundance', data=temp_df,
                  color="0", alpha=.5)
    plt.xticks(fontsize=15)
    plt.title(title + "\n", fontsize=20)
    plt.savefig(os.path.join(output_dir, title + '.svg'))
    plt.close()

# TODO group the bars by mark, with a bar for each sample group?

import sys
import os
import pandas as pd
import re
import numpy as np
from Bio import SeqIO
from pyteomics import parser
from tqdm import tqdm
import itertools

sys.stdout.write("Imported required packages successfully.\n")

# set the master table for modification mass shifts
MZSHIFT_DICT = {'[nPR]':'[+56]',
                '[AC]': '[+42]',
                '[PR]': '[+56]',
                '[ME1]': '[+70]',
                '[ME2]': '[+28]',
                '[ME3]': '[+42]',
                '[PH]': '[+80]',
                '[nME3PR]': '[+98.1]',
                '[nME2PR]': '[+84.1]',
                '[nME1PR]': '[+126.1]',
                '[nPR2]': '[+112.1]',
                '[nAC]': '[+98]',
                '[nACPR]': '[+84.1]'}

##
## TODO PARSER
##

# read in FASTA
fasta = os.path.join(os.getcwd(), "../../../resources/fasta/Histones.fasta")

# read in a user-defined list of modifications for each protein
mod_list = pd.read_csv(os.path.join(os.getcwd(), "../../histone_assay/data/histone_mod_mastersheet.csv"))

## OPTIONS?
# digestion enzyme/rule

##


## digest protein(s) into peptides
def digest_proteins(fasta_file):
    protein_df = pd.DataFrame()  # Initialize a dataframe to store results
    counter = 0
    for protein in tqdm(SeqIO.parse(fasta_file, "fasta")):

        protein_sequence = str(protein.seq).upper()

        # cleave initial "start" methionine if present
        if protein_sequence[0] == "M":
            protein_sequence = protein_sequence[1:]

        # digest proteins into peptides following Arg-C cleavage (r'R')
        # glu-c == r'[DE]'
        # peptides = list(parser.cleave(protein_sequence, parser.expasy_rules['arg-c']))
        peptides = list(parser.cleave(protein_sequence, '[DE]'))

        # get the aa residue positions for the peptides to later link residue PTM positions to the peptide
        initial_index = []
        terminal_index = []
        for pep in peptides:
            initial_aa_pos = 1 + int(protein_sequence.index(pep))  # +1 for zero-based array indexing
            terminal_aa_pos = initial_aa_pos + len(pep)
            initial_index.append(initial_aa_pos)
            terminal_index.append(terminal_aa_pos)

        # add this protein to the dataframe
        new_df = pd.DataFrame({'protein': [protein.id] * len(peptides),
                               'peptide_sequence': peptides,
                               'initial_aa_index': initial_index,
                               'terminal_aa_index': terminal_index})

        protein_df = protein_df.append(new_df)

    protein_df = protein_df[protein_df['peptide_sequence'].notna()]  # drop any nan peptides(?)
    protein_df = protein_df[protein_df['peptide_sequence'].str.len() >= 4]  # drop peptides <4 aa
    protein_df = protein_df[protein_df['peptide_sequence'].str.len() <= 40]  # drop peptides >40 aa

    sys.stdout.write("Finished in silico digestion and generating peptides.\n")

    return(protein_df)

peptide_df = digest_proteins(fasta)


## group mods by position to separate protein-level PTM list into peptide-level PTM list

# force all characters to upper just in case
mod_list['mod_id'] = [i.upper() for i in mod_list['mod_id']]

# parse amino acid and residue position out of list with regex
index_list = []
for index, row in mod_list.iterrows():
    mod_aa = row['mod_id']
    aa_index = int(re.findall(r'[0-9]+', (re.findall(r'[A-Z]+[0-9]+', mod_aa))[0])[0])
    index_list.append(aa_index)

# make dataframe with the PTM identity and the residue position for the PTM mark
mod_list['mod_index'] = index_list

# initialize an empty dataframe to store the results
mastermod_df = pd.DataFrame(columns=['protein', 'mod_index', 'mod_id'])

# iterate over each peptide, linking the sequence to its residue modifications
for index, row in protein_df.iterrows():

    this_protein = row['protein']

    # subset for just the mods on this protein
    ptm_df = mod_list[mod_list['protein'] == this_protein]

    # retrieve all mods supplied by user
    mods = pd.DataFrame(ptm_df[(ptm_df['mod_index'] >= row['initial_aa_index']) &
                               (ptm_df['mod_index'] <= row['terminal_aa_index'])])

    ##
    ## begin generating modification annotations
    ##

    ## propionylate N-terminal K with biological mods (Kme1, Kme2, Kac)
    pep_start_aa = row['peptide_sequence'][0] + str(row['initial_aa_index'])
    if mods.mod_id.str.contains(pep_start_aa).any() & (row['peptide_sequence'][0] == "K"):
        mods_to_fix = mods[mods['mod_id'].str.contains(pep_start_aa)]
        for temp_index, mod in mods_to_fix.iterrows():

            # replace n-term Kme1 with me + 2 prop
            if "[ME1]" in mod['mod_id']:
                new_modid = (pep_start_aa + "[nME1PR2]")

            # replace n-term Kme2 with me + 1 prop
            elif "[ME2]" in mod['mod_id']:
                new_modid = (pep_start_aa + "[nME2PR2]")

            # replace n-term Kme3 with me + 1 prop
            elif "[ME3]" in mod['mod_id']:
                new_modid = (pep_start_aa + "[nME3PR]")

            # replace n-term Kac with ac + 1 prop
            elif "[AC]" in mod['mod_id']:
                new_modid = (pep_start_aa + "[nACPR]")

            else:
                new_modid = mod['mod_id']

            prot = mods_to_fix[mods_to_fix['mod_id'] == mod['mod_id']]['protein']

            new_modindex = int(mods_to_fix[mods_to_fix['mod_id'] == mod['mod_id']]['mod_index'])
            replacement_values = [this_protein, new_modid, new_modindex]
            mods.loc[mods['mod_id'] == mod['mod_id']] = replacement_values

    elif mods.mod_id.str.contains(r'\[PH\]').any():
        mods_to_fix = mods[mods['mod_id'].str.contains(r'[STY]')]
        if mods_to_fix > 1:
            for temp_index, mod in mods_to_fix.iterrows():
                temp_mod_id = np.nan
                new_row = pd.Series([this_protein, '', mod['mod_index']], index=mods.columns)
                mods = mods.append(new_row, ignore_index=True)
        else:
            temp_mod_id = ''

    ## mark all N-terminal residues for propionylation
    if row['peptide_sequence'][0] == 'K':
        temp_mod_id = row['peptide_sequence'][0] + str(row['initial_aa_index']) + '[nPR2]'
    else:
        temp_mod_id = row['peptide_sequence'][0] + str(row['initial_aa_index']) + '[nPR]'
    new_row = pd.Series([this_protein, temp_mod_id, row['initial_aa_index']], index=mods.columns)
    mods = mods.append(new_row, ignore_index=True)

    ## propionylate all other lysine residues
    # first get the indices for each lysine in the sequence
    k_index = [m.start() for m in re.finditer('K', row['peptide_sequence'][1:])]  # skip N-term K
    k_index = [x + 1 for x in k_index]  # re-base +1 to account for skipping N-term residue

    # then make a K[PR] modification for each lysine
    for k in k_index:
        temp_mod_id = 'K' + str(k + row['initial_aa_index']) + '[PR]'
        new_row = pd.Series([this_protein, temp_mod_id, k + row['initial_aa_index']],
                            index=mods.columns)
        mods = mods.append(new_row, ignore_index=True)

    ##
    ## end generating modification annotations
    ##

    # group modifications by position into list-of-lists
    temp = mods.groupby('mod_index')['mod_id'].apply(list).reset_index()

    temp['PeptideSeq'] = row['peptide_sequence']
    temp['protein'] = this_protein
    temp['initial_aa_index'] = row['initial_aa_index']
    temp = temp.loc[temp.astype(str).drop_duplicates().index]

    mastermod_df = mastermod_df.append(temp)

sys.stdout.write("Finished generating modifications.\n")

##
## map shorthand histone marks to the peptide aa residue +[m/z] shift
##

# get the Cartesian product (every possible combination of values) from a list of lists
# from https://stackoverflow.com/questions/533905/get-the-cartesian-product-of-a-series-of-lists
mastermod_df['pep_prot_pair'] = mastermod_df['protein'] + mastermod_df['PeptideSeq']

# initialize an empty dataframe to store the results of loop below
out_df = pd.DataFrame(columns=['Protein Name', 'PeptideSeq', 'Peptide Note',
                               'PepModSeq', 'Peptide Modified Sequence'])

# iterate over each peptide (specific to a protein to avoid shared peptide problems)
for pep_prot_pair in mastermod_df['pep_prot_pair']:

    # parse some values for this peptide
    peptide = mastermod_df[mastermod_df['pep_prot_pair'] == pep_prot_pair]['PeptideSeq'][0]
    initial_aa = int(mastermod_df[mastermod_df['pep_prot_pair'] == pep_prot_pair]['initial_aa_index'][0])
    somelist = list(mastermod_df[mastermod_df['pep_prot_pair'] == pep_prot_pair]['mod_id'])

    # initialize empty lists to store the peptide notes (histone marks) and the modified sequences
    temp_pepnote = []
    temp_pepmod = []
    temp_pepmz = []

    # now generate all the modification combinations for this peptide
    for element in list(itertools.product(*somelist)):

        pepmod = ''.join(element)

        # remove any empty modifications resulting from phosphorylation status
        element = [t for t in element if t]

        # make a new copy of the peptide sequence for adding modified residues
        new_peptide = peptide
        pepmz = peptide

        for i in reversed(element):
            # get the mod aa index relative to the peptide string by subtracting initial pep aa index
            aa_index = int(re.findall('\d+', i)[0]) - initial_aa

            # get the mod without its residue index
            mod_i = "[" + re.search(r'.*?\[(.*)].*', i).group(1) + "]"
            char_i = i[0] + "[" + re.search(r'.*?\[(.*)].*', i).group(1) + "]"

            # lookup the mz shift for that modification
            mzshift = i[0] + MZSHIFT_DICT[mod_i]

            # replace the aa index in the peptide with the modified aa
            new_peptide = new_peptide[:aa_index] + char_i + new_peptide[(aa_index + 1):]

            pepmz = pepmz[:aa_index] + mzshift + pepmz[(aa_index + 1):]

        # add the results for this mod combination to a list for this peptide
        temp_pepnote.append(pepmod)
        temp_pepmod.append(new_peptide)
        temp_pepmz.append(pepmz)

    prot = mastermod_df[mastermod_df['pep_prot_pair'] == pep_prot_pair]['protein'][0]
    new_df = pd.DataFrame({"Protein Name": ([prot] * len(temp_pepnote)),
                           "PeptideSeq": ([peptide] * len(temp_pepnote)),
                           "Peptide Note": temp_pepnote,
                           "PepModSeq": temp_pepmod,
                           "Peptide Modified Sequence": temp_pepmz})

    out_df = out_df.append(new_df, ignore_index=True)

out_df = out_df.drop_duplicates(keep="first")

sys.stdout.write("Finished finding all possible mod combinations for each peptide.\n")

# parse a Skyline-acceptable "Insert Peptides" csv
skyline_df = out_df[['Peptide Modified Sequence', 'Protein Name', 'Peptide Note']]
skyline_df.to_csv(os.path.join(os.getcwd(), '../../collab_greer/data/skyline_insert_peptides.csv'),
                  index=False, header=False)
out_df.to_csv(os.path.join(os.getcwd(), '../../collab_greer/data/skyline_insert_peptides_annotated.csv'),
              index=False)

out_df.tail(25)
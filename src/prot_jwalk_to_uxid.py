#!/usr/bin/env python

import pandas as pd 
import argparse
import os
from Bio import SeqIO, pairwise2

desc = """Kai Kammer - 2019-09. 
Script to convert jwalk output to uxid. 
Requires a jwalk output file and and file associating uniprot ids with pdb chain ids which
can be acquired through the prot_match_pdb_chains.py script
"""


parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input files separated by spaces. Either .txt (jwalk) or .csv (uniprot to pdb chain id) files (at least one of each required)")
parser.add_argument('-o', '--outname', action="store", dest="outname", default='uxid_distances.csv',
                    help="Name of the output file")
args = parser.parse_args()

def prepare_jwalk_df(df):
    df['Atom1'] = df['Atom1'].str.replace('--', '-')
    df['Atom2'] = df['Atom2'].str.replace('--', '-')
    df['AminoAcid1'], df['Pos1'], df['Chain1'], df['AtomType1'] = df['Atom1'].str.split('-').str
    df['AminoAcid2'], df['Pos2'], df['Chain2'], df['AtomType2'] = df['Atom2'].str.split('-').str
    return df

def merge_uni_with_jwalk(df_uni, df_jwalk):
    df = pd.merge(df_jwalk, df_uni, right_on=['PDBFile', 'PDB_ChainID'], left_on=['Model', 'Chain1'], how='inner')
    df = pd.merge(df, df_uni, right_on=['PDBFile', 'PDB_ChainID'], left_on=['Model', 'Chain2'], how='inner', suffixes=('1', '2'))
    df_no_match = df[df.isna().any(axis=1)] 
    if len(df_no_match) > 0:
        print("No match found for\n", df_no_match)
    return df.dropna()

def get_uxid_df(df):
    df['uxID'] =  df['UniprotID_Full1'] + ':' + df['Pos1'] + ':x:' + df['UniprotID_Full2'] + ':' + df['Pos2']
    df['uxID_inv'] =  df['UniprotID_Full2'] + ':' + df['Pos2'] + ':x:' + df['UniprotID_Full1'] + ':' + df['Pos1']
    return df

def main():
    jwalk_df_list = []
    uni_to_pdb_df_list = []
    for inp in args.input:
        if ".txt" in inp:
            df_tmp = pd.read_csv(inp, delim_whitespace=True)
            jwalk_df_list.append(df_tmp.dropna(axis='columns'))
        elif ".csv" in inp:
            uni_to_pdb_df_list.append(pd.read_csv(inp))
    assert jwalk_df_list and uni_to_pdb_df_list, "Either no txt or csv file found in input; the script needs at least one of each"
    df_jwalk = pd.concat(jwalk_df_list, sort=True)
    #print(df_jwalk)
    df_jwalk_prep = prepare_jwalk_df(df_jwalk.copy())
    #print(df_jwalk)
    df_uni = pd.concat(uni_to_pdb_df_list)
    df_merge = merge_uni_with_jwalk(df_uni, df_jwalk_prep)
    df = get_uxid_df(df_merge)
    df = df.drop(columns=[col for col in df.columns if col not in df_jwalk.columns and col not in ['uxID', 'uxID_inv']])
    df = df.drop(columns=['Index','Atom1','Atom2'])
    df.to_csv(args.outname, index=False)
    print(f"Output written to {args.outname}")

if __name__ == "__main__":
    main()
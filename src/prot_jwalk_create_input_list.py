#!/usr/bin/env python

import pandas as pd 
import argparse
import os
from Bio import SeqIO, pairwise2

desc = """Kai Kammer - 2019-09. 
Script to create input crosslink list for jwalk by combining xquest output and  
combining it with the output of the prot_match_pdb_chains.py script
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input files separated by spaces. \
                        Either from the prot_match_pdb_chains script or xquest output (at least one of each required)")
parser.add_argument('-o', '--outname', action="store", dest="outname", default='jwalk_input',
                    help="Name of the output file")
args = parser.parse_args()

def merge_xq_chains(df_xq, df_chain):
    df = pd.merge(df_xq, df_chain, left_on=['Protein1'], right_on=['UniprotID_Full'], how='inner')
    # the pdb model for the linked proteins has to be the same; let's drop it here to avoid a duplicate column
    # also duplicated UniprotIDs are not necessary
    df = df.drop(columns=['UniprotID','UniprotID_Full','PDBFile'])
    df = pd.merge(df, df_chain, left_on=['Protein2'], right_on=['UniprotID_Full'], how='inner', suffixes=('1', '2'))
    df = df.drop(columns=['UniprotID','UniprotID_Full'])
    return df

def write_out_files(df):
    models = df['PDBFile'].unique()
    df_n = pd.DataFrame()
    df_n['tmp'] = df['AbsPos1'].astype(int).astype(str) + '|' + df['PDB_ChainID1'] + '|' + df['AbsPos2'].astype(int).astype(str) + '|' + df['PDB_ChainID2'] 
    for model in models:
        outname = f"{args.outname}_{model.split('.')[0]}.txt"
        df_n.to_csv(outname, index=False, header=False)
        print(f"Output written to {outname}")

def main():
    chain_df_list = []
    xquest_df_list = []
    for inp in args.input:
        df_tmp = pd.read_csv(inp)
        if 'PDBFile' in df_tmp.columns:
            chain_df_list.append(df_tmp)
        elif 'Type' in df_tmp.columns:
            xquest_df_list.append(df_tmp)
        else:
            print(f"Error: the file type of {inp} was not recognized. Exiting")
            exit(1)
    assert chain_df_list and xquest_df_list, "Either no chain match or xquest file found in input; the script needs at least one of each"
    df_xq = pd.concat(xquest_df_list)
    df_chain = pd.concat(chain_df_list)
    df = merge_xq_chains(df_xq, df_chain)
    write_out_files(df)

if __name__ == "__main__":
    main()
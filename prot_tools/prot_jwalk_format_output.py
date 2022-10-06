#!/usr/bin/env python

import pandas as pd
import argparse
import prot_lib

desc = """Kai Kammer - 2019-09. 
Script to convert jwalk output to uxid. 
Requires a jwalk output file and and file associating uniprot ids with pdb chain ids which
can be acquired through the prot_match_pdb_chains.py script
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input files separated by spaces. Either .txt (jwalk) or .csv (uniprot to pdb chain "
                         "id) files (at least one of each required)")
parser.add_argument('-off', '--offsets_file', action="store", required=True,
                    help="File containing fasta-pdb offsets (required!)")
parser.add_argument('-f', '--format', action="store", dest="format", default="unspecific",
                    help="Output format. xquest creates uxid compatible with xquest output. "
                         "Possible values: unspecific (default), xquest")
parser.add_argument('-o', '--outname', action="store", dest="outname", default='prot_distances.csv',
                    help="Name of the output file")
args = parser.parse_args()


def prepare_jwalk_df(df):
    df['Atom1'] = df['Atom1'].str.replace('--', '-')
    df['Atom2'] = df['Atom2'].str.replace('--', '-')
    df[[prot_lib.COL_AMINO_ACID_1, prot_lib.COL_POS_1, prot_lib.COL_CHAIN_1, 'AtomType1']] = df['Atom1'].str.split('-',
                                                                                                                   expand=True)
    df[[prot_lib.COL_AMINO_ACID_2, prot_lib.COL_POS_2, prot_lib.COL_CHAIN_2, 'AtomType2']] = df['Atom2'].str.split('-',
                                                                                                                   expand=True)
    df[[prot_lib.COL_POS_1, prot_lib.COL_POS_2]] = df[[prot_lib.COL_POS_1, prot_lib.COL_POS_2]].astype(int)
    return df


def merge_uni_with_jwalk(df_uni, df_jwalk):
    df = pd.merge(df_jwalk, df_uni, right_on=[prot_lib.COL_PDB_FILE, prot_lib.COL_PDB_CHAIN_ID],
                  left_on=['Model', prot_lib.COL_CHAIN_1], how='inner')
    df = pd.merge(df, df_uni, right_on=[prot_lib.COL_PDB_FILE, prot_lib.COL_PDB_CHAIN_ID],
                  left_on=['Model', prot_lib.COL_CHAIN_2], how='inner',
                  suffixes=('1', '2'))
    if 'Offset' in df_uni:
        df['Offset1'] = df['Offset1'].astype(int)
        df[prot_lib.COL_POS_1] += df['Offset1']
        df[prot_lib.COL_POS_2] += df['Offset2']
    df_no_match = df[df.isna().any(axis=1)]
    if len(df_no_match) > 0:
        print("No match found for\n", df_no_match)
    return df.dropna()


def get_uxid_df(df):
    df[prot_lib.COL_UXID] = df[prot_lib.COL_UNIPROT_ID_1] + ':' + df[prot_lib.COL_POS_1].astype(str) + ':x:' + df[
        prot_lib.COL_UNIPROT_ID_2] + ':' + df[
                                prot_lib.COL_POS_2].astype(str)
    df[prot_lib.COL_UXID_INV] = df[prot_lib.COL_UNIPROT_ID_2] + ':' + df[prot_lib.COL_POS_2].astype(str) + ':x:' + df[
        prot_lib.COL_UNIPROT_ID_1] + ':' + df[
                                    prot_lib.COL_POS_1].astype(str)
    return df


def main():
    df_offsets = None
    if args.offsets_file:
        df_offsets = pd.read_csv(args.offsets_file)
    jwalk_df_list = []
    for inp in args.input:
        if ".txt" in inp:
            df_tmp = pd.read_csv(inp, delim_whitespace=True)
            jwalk_df_list.append(df_tmp.dropna(axis='columns'))
    assert jwalk_df_list and df_offsets is not None, "Either no jwalk output or offsets file given; the script needs at least one of each"
    df_jwalk = pd.concat(jwalk_df_list, sort=True)
    df_jwalk_prep = prepare_jwalk_df(df_jwalk.copy())
    df_merge = merge_uni_with_jwalk(df_offsets, df_jwalk_prep)
    if args.format == "xquest":
        df = get_uxid_df(df_merge)
    else:
        df = df_merge
    df = df.drop(columns=[col for col in df.columns if
                          col not in df_jwalk.columns and col not in [prot_lib.COL_UXID, prot_lib.COL_UXID_INV,
                                                                      prot_lib.COL_POS_1,
                                                                      prot_lib.COL_POS_2, prot_lib.COL_AMINO_ACID_1,
                                                                      prot_lib.COL_AMINO_ACID_2, prot_lib.COL_CHAIN_1,
                                                                      prot_lib.COL_CHAIN_2, prot_lib.COL_UNIPROT_ID_1,
                                                                      prot_lib.COL_UNIPROT_ID_2]])
    df = df.drop(columns=['Index', 'Atom1', 'Atom2'], errors='ignore')
    df = df.rename(columns={'Model': prot_lib.COL_PDB_FILE})
    df.to_csv(args.outname, index=False)
    print(f"Output written to {args.outname}")


if __name__ == "__main__":
    main()

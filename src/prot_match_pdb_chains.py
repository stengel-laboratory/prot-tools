#!/usr/bin/env python

import pandas as pd
import argparse
import os
from Bio import SeqIO#, pairwise2

desc = """Kai Kammer - 2019-09. 
Script to match chains in pdb files with entries in a fasta file. 
Right now only works if the pdb contains seqres information matching the fasta file. 
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input files separated by spaces. Either .pdb or .fasta files (at least one of each "
                         "required)")
parser.add_argument('-o', '--outname', action="store", dest="outname", default='uni_to_pdb_chain.csv',
                    help="Name of the output file")
args = parser.parse_args()


def get_uni_ids(fasta_file_list):
    uni_dict = {}
    for fasta_file in fasta_file_list:
        uni_short_to_full_dict = {}
        fasta_records = SeqIO.parse(fasta_file, "fasta")
        for seq in fasta_records:
            id_split = seq.id.split('|')
            if len(id_split) > 1:
                uid = id_split[1]
            else:
                uid = seq.id
            uni_short_to_full_dict[uid] = seq.id
        uni_dict[os.path.basename(fasta_file)] = uni_short_to_full_dict
    return uni_dict


def get_pdb_chains(pdb_file_list):
    pdb_dict = {}
    for pdb_file in pdb_file_list:
        chain_id_to_uni_dict = {}
        # converting SeqIO records from generator to list; otherwise would be empty after assert statement
        pdb_records = list(SeqIO.parse(pdb_file, "pdb-seqres"))
        assert list(pdb_records), f"No seqres information found in {pdb_file}. Exiting"
        for seq in pdb_records:
            chain_id_to_uni_dict[seq.id.split(':')[1]] = seq.dbxrefs[0].split(':')[1]
        pdb_dict[os.path.basename(pdb_file)] = chain_id_to_uni_dict
    return pdb_dict


def get_merge_df(uni_dict, pdb_dict):
    df_uni_list = []
    df_pdb_list = []
    for uni_file, uni_short_to_full_dict in uni_dict.items():
        df = pd.DataFrame(
            {'UniprotID': list(uni_short_to_full_dict.keys()), 'UniprotID_Full': list(uni_short_to_full_dict.values())})
        df['FastaFile'] = uni_file
        df_uni_list.append(df)
    for pdb_file, chain_id_to_uni_dict in pdb_dict.items():
        df = pd.DataFrame(
            {'UniprotID': list(chain_id_to_uni_dict.values()), 'PDB_ChainID': list(chain_id_to_uni_dict.keys())})
        df['PDBFile'] = pdb_file
        df_pdb_list.append(df)
    df_uni = pd.concat(df_uni_list)
    df_pdb = pd.concat(df_pdb_list)
    df = pd.merge(df_uni, df_pdb, on='UniprotID', how='outer')
    df_no_match = df[df.isna().any(axis=1)]
    if len(df_no_match) > 0:
        print("No match found for\n", df_no_match)
    return df


def main():
    pdb_file_list = []
    fasta_file_list = []
    for inp in args.input:
        if ".pdb" in inp:
            pdb_file_list.append(inp)
        elif ".fasta" in inp:
            fasta_file_list.append(inp)
    assert pdb_file_list or not fasta_file_list, "Either no fasta or pdb file found in input; the script needs at " \
                                                 "least one of each "
    uni_dict = get_uni_ids(fasta_file_list)
    pdb_dict = get_pdb_chains(pdb_file_list)
    df = get_merge_df(uni_dict, pdb_dict)
    df.to_csv(args.outname, index=False)
    print(f"Output written to {args.outname}")


if __name__ == "__main__":
    main()

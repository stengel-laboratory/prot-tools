#!/usr/bin/env python
# %%
import sys, os
from typing import List

import numpy as np
import prot_lib
import pandas as pd
import argparse


from pKAI.pKAI import pKAI
from pandas import DataFrame

desc = """Kai Kammer - 2022-05. 
Create a csv file containing pKa values of lysines for use in a kinetic model.
Accepts a single or multiple pdb file as input and requires pKAI installed (https://pypi.org/project/pKAI).
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# not sure if just using type=str is better
parser.add_argument('input', action="store", default=None, type=argparse.FileType('r'), nargs='+',
                    help="List of input files separated by spaces. Minimum requirement is one pdb file. ")
parser.add_argument('-off', '--offsets_file', action="store",
                    help="Optional file containing fasta-pdb offsets")
parser.add_argument('-o', '--outname', type=argparse.FileType('w'),
                    default='prot_pka.csv', help="Name of the output file")
args = parser.parse_args()


def get_pka_df(pdb_file, df_offsets=None) -> "pd.DataFrame":
    pdb_file = pdb_file.name
    sys.stdout = open(os.devnull, 'w')  # avoid printing pKAI output
    pks = pKAI(pdb_file)
    sys.stdout = sys.__stdout__
    df = pd.DataFrame(pks, columns=[prot_lib.COL_PDB_CHAIN_ID, prot_lib.COL_POS, prot_lib.COL_RES, prot_lib.COL_PKA])
    df = df[df[prot_lib.COL_RES].str.contains(prot_lib.STR_LYS)]
    df[prot_lib.COL_PDB_FILE] = pdb_file

    if df_offsets is not None:
        pdb_file = df_offsets[prot_lib.COL_PDB_FILE].values[0]
        df[prot_lib.COL_POS] = df.apply(
            lambda x: x[prot_lib.COL_POS] + prot_lib.get_pdb_fasta_offset(pdb_file, df_offsets,
                                                                          x[prot_lib.COL_PDB_CHAIN_ID]),
            axis=1)
    return df


def main():
    df_offsets = None
    if args.offsets_file:
        df_offsets = pd.read_csv(args.offsets_file)
    df_pka_list: List[DataFrame] = []
    for inp in args.input:
        if prot_lib.is_pdb_file(inp):
            df_pka_list.append(get_pka_df(inp, df_offsets))
    assert df_pka_list, "No valid pdb file found in input; stopping"
    df_pka = pd.concat(df_pka_list)
    df_pka.to_csv(args.outname, index=False)


if __name__ == "__main__":
    main()

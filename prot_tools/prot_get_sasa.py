#!/usr/bin/env python
# %%
import sys
from typing import List

import numpy as np
import prot_lib
import pandas as pd
import argparse
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

from pandas import DataFrame

desc = """Kai Kammer - 2022-05. 
Create a csv file containing SASA values of lysines (via ShrakeRupley) for use in a kinetic model.
Accepts a single or multiple pdb file as input and requires biopython installed.
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# not sure if just using type=str is better
parser.add_argument('input', action="store", default=None, type=argparse.FileType('r'), nargs='+',
                    help="List of input files separated by spaces. Minimum requirement is one pdb file")
parser.add_argument('-off', '--offsets_file', action="store",
                    help="Optional file containing fasta-pdb offsets")
parser.add_argument('-o', '--outname', type=argparse.FileType('w'),
                    default='prot_sasa.csv', help="Name of the output file")
args = parser.parse_args()


def get_sasa_df(pdb_file, df_offsets=None) -> pd.DataFrame:
    pdb_file = pdb_file.name
    list_sasa = []
    p = PDBParser(QUIET=1)
    struct = p.get_structure("pdb", pdb_file)
    sr = ShrakeRupley()
    sr.compute(struct, level="R")  # compute SASA on residue level
    for chain in struct.get_chains():
        for res in chain.get_residues():
            if res.get_resname() == 'LYS':  # pre-filter for lys residues
                list_sasa.append((chain.get_id(), res.get_id()[1], res.get_resname(), res.sasa))
    df = pd.DataFrame(list_sasa, columns=[prot_lib.COL_PDB_CHAIN_ID, prot_lib.COL_POS, prot_lib.COL_RES, prot_lib.COL_SASA])
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
    df_sasa_list: List[DataFrame] = []
    for inp in args.input:
        if prot_lib.is_pdb_file(inp):
            df_sasa_list.append(get_sasa_df(inp, df_offsets))
    assert df_sasa_list, "No valid pdb file found in input; stopping"
    df_sasa = pd.concat(df_sasa_list)
    df_sasa.to_csv(args.outname, index=False)


if __name__ == "__main__":
    main()

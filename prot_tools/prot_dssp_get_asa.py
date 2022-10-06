#!/usr/bin/env python
import pandas as pd
import argparse
import prot_lib
import os.path
from collections import defaultdict
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP


desc = """Kai Kammer - 2021-03. 
Run dssp (needs to installed and in path) and extract ASA values
"""
parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", help="PDB file")
parser.add_argument('-off', '--offsets_file', action="store",
                    help="Optional file containing fasta-pdb offsets")
parser.add_argument('-o', '--output', action="store", default='dssp_asa.csv',
                    help="File to save to")
args = parser.parse_args()


def get_asa_df(pdb_file, df_offset=None):
    p = PDBParser()
    structure = p.get_structure("tmp", pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)
    asa_dict = defaultdict(list)
    pdbfile_basename = os.path.basename(pdb_file)
    if df_offset is None:
        print("Assume offset of 0 for asa dict")
    for k in dssp.keys():
        chain_id = k[0]
        if df_offset is not None:
            offset = prot_lib.get_pdb_fasta_offset(pdb_file=pdb_file, pdb_chain_id=chain_id, df_offset=df_offset)
        else:
            offset = 0
        asa_dict[prot_lib.COL_PDB_CHAIN_ID].append(chain_id)
        asa_dict[prot_lib.COL_POS].append(k[1][1] + offset)
        asa_dict[prot_lib.COL_ASA].append(dssp[k][3])
        asa_dict[prot_lib.COL_RES].append(dssp[k][1])
        asa_dict[prot_lib.COL_SEC_STRUC].append(dssp[k][2])
        asa_dict[prot_lib.COL_PDB_FILE].append(pdbfile_basename)
    df = pd.DataFrame(asa_dict)
    return df


def main():
    df_offsets = None
    if args.offsets_file:
        df_offsets = pd.read_csv(args.offsets_file)
    df = get_asa_df(args.input, df_offsets)
    df.to_csv(args.output, index=False)
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()

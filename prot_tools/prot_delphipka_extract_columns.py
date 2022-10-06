#!/usr/bin/env python
import pandas as pd
import argparse
import prot_lib

desc = """Kai Kammer - 2021-03. 
Parse DelphiPKa web server out into multiple columns
"""
parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", help="DelphiPKa output file")
parser.add_argument('-off', '--offsets_file', action="store",
                    help="Optional file containing fasta-pdb offsets. First pdb entry is assumed "
                         "to be the one matching the given dssp file")
parser.add_argument('-o', '--output', action="store", default='delphi_pka_cols.csv',
                    help="File to save to")
args = parser.parse_args()


def main():
    df = pd.read_csv(args.input)
    df[prot_lib.COL_POS] = df['ResName'].str[3:7].astype(int)
    df[prot_lib.COL_PDB_CHAIN_ID] = df['ResName'].str[7:]
    df[prot_lib.COL_RES] = df['ResName'].str[0:3]
    if args.offsets_file:
        df_offsets = pd.read_csv(args.offsets_file)
        pdb_file = df_offsets[prot_lib.COL_PDB_FILE].values[0]
        df[prot_lib.COL_POS] = df.apply(
            lambda x: x[prot_lib.COL_POS] + prot_lib.get_pdb_fasta_offset(pdb_file, df_offsets,
                                                                          x[prot_lib.COL_PDB_CHAIN_ID]),
            axis=1)
    df.to_csv(args.output, index=False)
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()

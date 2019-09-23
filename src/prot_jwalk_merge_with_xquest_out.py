#!/usr/bin/env python

import pandas as pd 
import argparse

desc = """Kai Kammer - 2019-09. 
Script to convert merge jwalk output with xquest output. 
Requires file(s) prepared by the prot_jwalk_to_uxid.py script and xquest output.
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input files separated by spaces. \
                        Either from the prot_jwalk_to_uxid script xquest output (at least one of each required)")
parser.add_argument('-o', '--outname', action="store", dest="outname", default='xquest_jwalk_merge.csv',
                    help="Name of the output file")
parser.add_argument('-l', '--long', action="store_true", dest="long",
                    help="Output complete xquest file merged with distances; otherwise will output compact form\
                         only containing distances")
parser.add_argument('-e', '--exp_name', action="store_true", dest="exp_name",
                    help="Optionally provide experiment names for the corresponding model files")
args = parser.parse_args()


# merges on uxID and inverse uxID since the construction of the uxID can be ambiguous
def merge_xq_jwalk(df_xq, df_jwalk):
    df1 = pd.merge(df_xq, df_jwalk, on='uxID')
    print(f"Matched {len(df1)} uxIDs on first merge")
    df_jwalk = df_jwalk.drop(columns=['uxID'])
    df2 = pd.merge(df_xq, df_jwalk, left_on='uxID', right_on='uxID_inv')
    print(f"Matched {len(df2)} uxIDs on second merge")
    df = pd.concat([df1, df2], sort=False)
    print(f"Combined into {len(df)} entries. xQuest input had {len(df_xq)} and jwalk input had {len(df_jwalk)} entries")
    df = df.drop_duplicates(subset=['uxID', 'Model'])
    df = df.drop(columns='uxID_inv')
    print(f"After dropping duplicates {len(df)} entries are left")
    if not args.long:
        df = df.drop(columns=[col for col in df.columns if col in df_xq and col != 'uxID'])
    df = df.sort_values(['uxID'])
    return df


def assign_exp_names(df):
    exp_dict = {}
    for model in df['Model'].unique():
        sel = input(f"Please type the experiment name for model file {model}:")
        exp_dict[model] = sel
    df['exp_name'] = df['Model'].map(exp_dict)
    df = df.drop(columns='Model')
    return df


def main():
    jwalk_df_list = []
    xquest_df_list = []
    for inp in args.input:
        df_tmp = pd.read_csv(inp)
        if 'SASD' in df_tmp.columns:
            jwalk_df_list.append(df_tmp)
        elif 'Type' in df_tmp.columns:
            xquest_df_list.append(df_tmp)
        else:
            print(f"Error: the file type of {inp} was not recognized. Exiting")
            exit(1)
    assert jwalk_df_list and xquest_df_list, "Either no txt or csv file found in input; the script needs at least one of each"
    df_xq = pd.concat(xquest_df_list)
    df_jwalk = pd.concat(jwalk_df_list)
    df = merge_xq_jwalk(df_xq, df_jwalk)
    if args.exp_name:
        df = assign_exp_names(df)
    df.to_csv(args.outname, index=False)
    print(f"Output written to {args.outname}")

if __name__ == "__main__":
    main()
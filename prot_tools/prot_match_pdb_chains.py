#!/usr/bin/env python

import pandas as pd
import argparse
import os
from Bio import SeqIO, pairwise2
from collections import defaultdict
import prot_lib

desc = """Kai Kammer - 2019-09. 
Script to match chains in pdb files with entries in a fasta file. 
Right now only works if the pdb contains seqres information matching the fasta file. 
Also accounts for offsets between pdb and fasta (can be turned off with -n)
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# not sure if just using type=str is better
parser.add_argument('input', action="store", default=None, type=argparse.FileType('r'), nargs='+',
                    help="List of input files separated by spaces. Either pdb or fasta files (at least one of each "
                         "required). Files are validated with biopython; suffixes do not matter")
parser.add_argument('-n', '--no_compute_offsets', action="store_true", default=False,
                    help="Optionally do not compute offsets between fasta and pdb. Much faster")
parser.add_argument('-o', '--outname', type=argparse.FileType('w'),
                    default='prot_uni_to_pdb_chain.csv', help="Name of the output file")
args = parser.parse_args()


def get_uni_ids(fasta_file_list):
    uni_dict = {}
    file_to_uni_id_dict = defaultdict(list)
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
            file_to_uni_id_dict[fasta_file].append(uid)
        uni_dict[fasta_file] = uni_short_to_full_dict

    return uni_dict, file_to_uni_id_dict


def get_merge_df(uni_dict, pdb_dict, offset_dict=None):
    df_uni_list = []
    df_pdb_list = []
    for uni_file, uni_short_to_full_dict in uni_dict.items():
        df_tmp = pd.DataFrame(
            {prot_lib.COL_UNIPROT_ID: list(uni_short_to_full_dict.keys()), prot_lib.COL_UNIPROT_ID_FULL: list(uni_short_to_full_dict.values())})
        df_tmp[prot_lib.COL_FASTA_FILE] = os.path.basename(uni_file)
        df_uni_list.append(df_tmp)
    for pdb_file, chain_id_to_uni_dict in pdb_dict.items():
        df_tmp = pd.DataFrame(
            {prot_lib.COL_UNIPROT_ID: list(chain_id_to_uni_dict.values()), prot_lib.COL_PDB_CHAIN_ID: list(chain_id_to_uni_dict.keys())})
        df_tmp[prot_lib.COL_PDB_FILE] = os.path.basename(pdb_file)
        df_pdb_list.append(df_tmp)

    df_uni = pd.concat(df_uni_list)
    df_pdb = pd.concat(df_pdb_list)
    df = pd.merge(df_uni, df_pdb, on=prot_lib.COL_UNIPROT_ID, how='outer')
    df_no_match = df[df.isna().any(axis=1)]
    if len(df_no_match) > 0:
        print("No match found for\n", df_no_match)
    if offset_dict:
        df_offset_list = []
        for pdb_file, chain_id_to_offset_dict in offset_dict.items():
            df_tmp = pd.DataFrame(
                {prot_lib.COL_OFFSET: list(chain_id_to_offset_dict.values()), prot_lib.COL_PDB_CHAIN_ID: list(chain_id_to_offset_dict.keys())})
            df_tmp[prot_lib.COL_PDB_FILE] = os.path.basename(pdb_file)
            df_offset_list.append(df_tmp)
        df_offset = pd.concat(df_offset_list)
        df = pd.merge(df, df_offset, on=[prot_lib.COL_PDB_FILE, prot_lib.COL_PDB_CHAIN_ID], how='outer')
        df_no_match = df[df.isna().any(axis=1)]
        if len(df_no_match) > 0:
            print("No offsets found for\n", df_no_match, "\nOffsets will be set to 0")
        df = df.fillna(0)
        df[prot_lib.COL_OFFSET] = df[prot_lib.COL_OFFSET].astype(int)
    return df


def get_offsets(fasta_file_to_uni_id, uni_id_to_file_pdb):
    SCORE = 0  # lowest score to accept a match between UniID and PDB; note that if only a subunit in a big complex is a match, then very low scores (~-1000) are possible
    SEQ_ID = 20  # minmum sequence identity between UniID and PDB sequence matches (0 - 100)
    SEQ_ID_GAPLESS = 90  # minimum gapless sequence id. between UniID and PDB seq. matches (0 - 100)
    offset_dict = {}
    for fasta_file, uni_id_list in fasta_file_to_uni_id.items():
        pdb_chain_dict = {}
        fasta_records = SeqIO.parse(fasta_file, "fasta")
        for uni_id in uni_id_list:
            pdb_file = uni_id_to_file_pdb[uni_id]
            pdb_records = list(SeqIO.parse(pdb_file, "pdb-atom"))
            for fasta_rec in fasta_records:

                for pdb_rec in pdb_records:

                    # print(pdb_rec)
                    aligned = pairwise2.align.localms(fasta_rec.seq, pdb_rec.seq, 5, -4, -10, -0.5,
                                                      penalize_end_gaps=False)
                    aligned_a, aligned_b, score, begin, end = aligned[0]

                    seq_id, g_seq_id = _calculate_identity(aligned_a, aligned_b)
                    if score > SCORE and seq_id > SEQ_ID and g_seq_id > SEQ_ID_GAPLESS:
                        offset = begin - pdb_rec.annotations['start'] + 1
                        chain = pdb_rec.annotations['chain']
                        # print(score, seq_id, g_seq_id)
                        print(pdb_file, fasta_rec.id, pdb_rec.id, chain)
                        print("alignment sequence start-end, chain", (begin, end, chain))
                        print("offset from pdb to fasta index (fasta_pos-offset=pdb_pos)", offset)
                        pdb_chain_dict[chain] = offset
            offset_dict[os.path.basename(pdb_file)] = pdb_chain_dict
    return offset_dict


def _calculate_identity(seqA, seqB):
    """
    Returns the percentage of identical characters between two sequences.
    Assumes the sequences are aligned. Adapted from: https://gist.github.com/JoaoRodrigues/8c2f7d2fc5ae38fc9cb2
    Modified to take the minimum sequence length of A and B instead of just taking A's length
    """

    sa, sb, sl = seqA, seqB, min(len(seqA), len(seqB))
    matches = [sa[i] == sb[i] for i in range(sl)]
    seq_id = (100 * sum(matches)) / sl

    gapless_sl = sum([1 for i in range(sl) if (sa[i] != '-' and sb[i] != '-')])
    gap_id = (100 * sum(matches)) / gapless_sl
    return (seq_id, gap_id)


def main():
    pdb_file_list = []
    fasta_file_list = []
    for inp in args.input:
        print(inp.name)
        if prot_lib.is_pdb_file(inp.name):
            pdb_file_list.append(inp.name)
        elif prot_lib.is_fasta_file(inp.name):
            fasta_file_list.append(inp.name)
    assert pdb_file_list or not fasta_file_list, "Either no fasta or pdb file found in input; the script needs at " \
                                                 "least one of each "
    uni_dict, fasta_file_to_uni_id = get_uni_ids(fasta_file_list)
    pdb_dict, uni_id_to_file_pdb = prot_lib.get_pdb_chains(pdb_file_list)
    offset_dict = None
    if not args.no_compute_offsets:
        offset_dict = get_offsets(fasta_file_to_uni_id, uni_id_to_file_pdb)
    df = get_merge_df(uni_dict, pdb_dict, offset_dict)
    df.to_csv(args.outname.name, index=False)
    print(f"Output written to {args.outname.name}")


if __name__ == "__main__":
    main()

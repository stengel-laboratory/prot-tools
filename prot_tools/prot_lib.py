from Bio.PDB import PDBParser
from Bio import SeqIO, SeqUtils
import pandas as pd

COL_FASTA_FILE = "FastaFile"
COL_PDB_FILE = 'PDBFile'
COL_PDB_CHAIN_ID = 'PDB_ChainID'
COL_OFFSET = 'Offset'
COL_POS = 'Pos'
COL_POS_1 = 'Pos1'
COL_POS_2 = 'Pos2'
COL_AMINO_ACID_1 = 'AminoAcid1'
COL_AMINO_ACID_2 = 'AminoAcid2'
COL_CHAIN_1 = 'Chain1'
COL_CHAIN_2 = 'Chain2'
COL_UNIPROT_ID = "UniprotID"
COL_UNIPROT_ID_1 = COL_UNIPROT_ID + '1'
COL_UNIPROT_ID_2 = COL_UNIPROT_ID + '2'
COL_UNIPROT_ID_FULL = COL_UNIPROT_ID + "_Full"
COL_UNIPROT_ID_FULL_1 = COL_UNIPROT_ID_FULL + "1"
COL_UNIPROT_ID_FULL_2 = COL_UNIPROT_ID_FULL + "2"
COL_RES = 'Res'
COL_ASA = 'ASA'
COL_SASA = 'SASA'
COL_PKA = 'pKa'
COL_SEC_STRUC = 'SecondaryStructure'
COL_DIST_EU = "eucDist"
COL_UXID = "uxID"
COL_UXID_INV = COL_UXID + "_inv"
COL_DIST_SASD = "SASD"
COL_ENERGY = "energy"
COL_PARAM = "Param"
COL_VALUE = "Value"
COL_ORIGIN = 'Origin'
COL_VALUE_NORMALIZED = COL_VALUE + "_Norm"
STR_LYS = "LYS"


def get_pdb_fasta_offset(pdb_file, df_offset, pdb_chain_id=None):
    if not pdb_chain_id:
        print(f"No chain id given for pdb fasta offset; using first chain: {pdb_chain_id}")
    if COL_PDB_FILE in df_offset and COL_PDB_CHAIN_ID in df_offset.columns:
        if pdb_file in df_offset[COL_PDB_FILE].values and pdb_chain_id in df_offset[COL_PDB_CHAIN_ID].values:
            return df_offset[(df_offset[COL_PDB_FILE] == pdb_file) & (df_offset[COL_PDB_CHAIN_ID] == pdb_chain_id)][
                COL_OFFSET].values[0]
        else:
            print(f"Warning: pdb file {pdb_file} or pdb chain {pdb_chain_id} not found in offsets file\n"
                  f"Assuming offset of 0")
            return 0
    else:
        print(f"Warning: column {COL_PDB_FILE} or {COL_PDB_CHAIN_ID} not found in offsets file\n"
              f"Assuming offset of 0")
        return 0


def is_pdb_file(file_path) -> "bool":
    p = PDBParser(QUIET=1)
    num_atoms = len(list(p.get_structure("pdb", file_path)))
    if num_atoms > 0:
        return True
    return False


# from: https://stackoverflow.com/a/44294079
def is_fasta_file(file_path):
    with open(file_path, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)


def is_csv_file(file_path):
    try:
        df = pd.read_csv(file_path)
        return df
    except pd.errors.ParserError:
        return None


def check_cols(file_path, cols):
    df = is_csv_file(file_path)
    if df is not None:
        if all(col in df.columns for col in cols):
            return True
    return False


def is_params_file(file_path):
    return check_cols(file_path=file_path, cols=[COL_PARAM])


def is_fasta_to_pdb_file(file_path):
    return check_cols(file_path=file_path, cols=[COL_UNIPROT_ID, COL_PDB_FILE, COL_OFFSET])


def is_sasa_file(file_path):
    return check_cols(file_path=file_path, cols=[COL_SASA])


def is_pka_file(file_path):
    return check_cols(file_path=file_path, cols=[COL_PKA])


def is_sasd_file(file_path):
    return check_cols(file_path=file_path, cols=[COL_DIST_SASD])


def check_for_cols(file_path, cols):
    with open(file_path, "r") as file:
        first_line = file.readline()
        if all(col in first_line for col in cols):
            return True
    return False


def get_fasta_records(fasta):
    with open(fasta, "r") as handle:
        record = list(SeqIO.parse(handle, "fasta"))
        return record


def get_prot_lys_pos_dict(fasta_records):
    prot_pos_dict = {}
    for record in fasta_records:
        prot_pos_dict[record.id] = [n + 1 for n, c in enumerate(record.seq) if c == 'K']
    return prot_pos_dict

def get_prot_lys_pos_dict_pdb(pdb_file):
    prot_pos_dict = {}
    chain_id_to_uni_dict = {}
    # converting SeqIO records from generator to list; otherwise would be empty after assert statement
    if ".cif" in pdb_file:
        parse_string = "cif-seqres"
    else:
        parse_string = "pdb-seqres"
    pdb_records = list(SeqIO.parse(pdb_file, parse_string))
    assert list(pdb_records), f"No seqres information found in {pdb_file}. Exiting"
    for record in pdb_records:
        uni_id = record.dbxrefs[0].split(':')[1]
        chain_id = record.id.split(':')[1]
        chain_id_to_uni_dict[chain_id] = uni_id
        prot_pos_dict[chain_id] = [n for n, c in enumerate(record.seq) if c == 'K']
    return chain_id_to_uni_dict, prot_pos_dict

def get_pdb_chains(pdb_file_list):
    pdb_dict = {}
    uni_id_to_file_dict = {}
    for pdb_file in pdb_file_list:
        if ".cif" in pdb_file:
            parse_string = "cif-seqres"
        else:
            parse_string = "pdb-seqres"

        chain_id_to_uni_dict = {}
        # converting SeqIO records from generator to list; otherwise would be empty after assert statement
        pdb_records = list(SeqIO.parse(pdb_file, parse_string))
        assert list(pdb_records), f"No seqres information found in {pdb_file}. Exiting"
        for record in pdb_records:
            uni_id = record.dbxrefs[0].split(':')[1]
            chain_id_to_uni_dict[record.id.split(':')[1]] = uni_id
            uni_id_to_file_dict[uni_id] = pdb_file
        pdb_dict[pdb_file] = chain_id_to_uni_dict
    return pdb_dict, uni_id_to_file_dict


def get_pdb_mol_weight(pdb_file):
    mw = 0
    if ".cif" in pdb_file:
        parse_string = "cif-seqres"
    else:
        parse_string = "pdb-seqres"
    pdb_records = list(SeqIO.parse(pdb_file, parse_string))
    for record in pdb_records:
        seq = record.seq
        seq = seq.replace("X", "V") # replace X (any aa) with valine which is close to average weight; otherwise this will not work
        mw += SeqUtils.molecular_weight(seq, seq_type="protein")
    return mw
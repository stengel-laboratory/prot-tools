import pytest
import os.path
import pandas as pd
from pandas.testing import assert_frame_equal
from prot_tools import prot_lib

INPUT = './test/input'
GROUND_TRUTH = os.path.join(INPUT, 'ground_truth')
PDB_FILE = os.path.join(INPUT, 'eg5_3wpn.pdb')
PKA_FILE = os.path.join(GROUND_TRUTH, 'prot_pka.csv')
SASA_FILE = os.path.join(GROUND_TRUTH, 'prot_sasa.csv')


def test_is_pdb():
    assert prot_lib.is_pdb_file(PDB_FILE)


def test_sasa():
    df_sasa = prot_lib.get_sasa_df(PDB_FILE)
    df_sasa_ref = pd.read_csv(SASA_FILE)
    assert_frame_equal(df_sasa, df_sasa_ref, check_like=True)


def test_pka():
    df_pka = prot_lib.get_pka_df(PDB_FILE)
    df_pka_ref = pd.read_csv(PKA_FILE)
    assert_frame_equal(df_pka, df_pka_ref, check_like=True)


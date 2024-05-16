import os
import sys

SCRIPT_DIR = os.path.join(
    os.path.split(os.path.dirname(os.path.abspath(__file__)))[0], "src", "charlie"
)
sys.path.append(SCRIPT_DIR)
from util import smk_base


def test_smk_base():
    assert smk_base("").endswith("charlie/")

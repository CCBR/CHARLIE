import os
import re
import sys

SCRIPT_DIR = os.path.join(
    os.path.split(os.path.dirname(os.path.abspath(__file__)))[0], "src", "charlie"
)
sys.path.append(SCRIPT_DIR)
from util import *


def test_smk_base():
    assert smk_base("").endswith("charlie/")


def test_semantic_version():
    assert bool(
        re.match(
            r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$",
            get_version(),
        )
    )

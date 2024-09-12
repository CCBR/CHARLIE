import os


def smk_base(rel_path):
    basedir = os.path.split(
        os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
    )[0]
    return os.path.join(basedir, rel_path)


def get_version():
    with open(smk_base("VERSION"), "r") as f:
        version = f.readline().strip()
    return version

## functions


def get_file_size(filename):
    filename = filename.strip()
    if check_readaccess(filename):
        return os.stat(filename).st_size


def get_peorse(wildcards):
    return SAMPLESDF["PEorSE"][wildcards.sample]


def check_existence(filename):
    """Checks if file exists on filesystem
    :param filename <str>: Name of file to check
    """
    filename = filename.strip()
    if not os.path.exists(filename):
        sys.exit("File: {} does not exists!".format(filename))
    return True


def check_readaccess(filename):
    """Checks permissions to see if user can read a file
    :param filename <str>: Name of file to check
    """
    filename = filename.strip()
    check_existence(filename)
    if not os.access(filename, os.R_OK):
        sys.exit(
            "File: {} exists, but user cannot read from file due to permissions!".format(
                filename
            )
        )
    return True


def check_writeaccess(filename):
    """Checks permissions to see if user can write to a file
    :param filename <str>: Name of file to check
    """
    filename = filename.strip()
    check_existence(filename)
    if not os.access(filename, os.W_OK):
        sys.exit(
            "File: {} exists, but user cannot write to file due to permissions!".format(
                filename
            )
        )
    return True


def append_files_in_list(flist, ofile):
    if not os.path.exists(ofile):
        print("FILE %s does not exist! Creating it!"%(ofile),flush=True)
        with open(ofile, "w") as outfile:
            for fname in flist:
                with open(fname) as infile:
                    l = infile.read()
                    l = l.strip()
                    outfile.write("%s\n" % (l))
    return True


##### load config and sample sheets #####

# check_readaccess("config/config.yaml")
# configfile: "config/config.yaml"

## set memory limit
## used for sambamba sort, etc
# MEMORYG="100G"


def _is_true(variable):
    if variable == True or variable == "True" or variable == "TRUE":
        return True
    else:
        return False


def _convert_to_int(variable):
    if variable:
        return 1  # True
    if not variable:
        return 0  # False
    return -1  # Unknown


# resouce absolute path
WORKDIR = config["workdir"]
SCRIPTS_DIR = config["scriptsdir"]
RESOURCES_DIR = config["resourcesdir"]
FASTAS_GTFS_DIR = config["fastas_gtfs_dir"]
RUN_CLEAR = _is_true(config["run_clear"])
RUN_DCC = _is_true(config["run_dcc"])
RUN_MAPSPLICE = _is_true(config["run_mapsplice"])
RUN_CIRCRNAFINDER = _is_true(config["run_circRNAFinder"])
RUN_NCLSCAN = _is_true(config["run_nclscan"])
RUN_FINDCIRC = _is_true(config["run_findcirc"])
N_RUN_CLEAR = _convert_to_int(RUN_CLEAR)
N_RUN_DCC = _convert_to_int(RUN_DCC)
N_RUN_MAPSPLICE = _convert_to_int(RUN_MAPSPLICE)
N_RUN_CIRCRNAFINDER = _convert_to_int(RUN_CIRCRNAFINDER)
N_RUN_NCLSCAN = _convert_to_int(RUN_NCLSCAN)
N_RUN_FINDCIRC = _convert_to_int(RUN_FINDCIRC)
MAPSPLICE_MIN_MAP_LEN = config["mapsplice_min_map_len"]
MAPSPLICE_FILTERING = config["mapsplice_filtering"]
FLANKSIZE = config['flanksize']
FIND_CIRC_DIR = config['find_circ_dir']

REF_DIR = join(WORKDIR, "ref")
if not os.path.exists(REF_DIR):
    os.mkdir(REF_DIR)
STAR_INDEX_DIR = join(REF_DIR, "STAR_no_GTF")
if not os.path.exists(STAR_INDEX_DIR):
    os.mkdir(STAR_INDEX_DIR)
# strip trailing slashes if any
for d in [
    WORKDIR,
    SCRIPTS_DIR,
    RESOURCES_DIR,
    FASTAS_GTFS_DIR,
    STAR_INDEX_DIR,
    REF_DIR,
]:
    d = d.strip("r\/")
BWA_INDEX = join(REF_DIR, "ref")

HOST = config["host"]  # hg38 or mm39
ADDITIVES = config["additives"]  # ERCC and/or BAC16Insert
ADDITIVES = ADDITIVES.replace(" ", "")
if ADDITIVES != "":
    HOST_ADDITIVES = HOST + "," + ADDITIVES
else:
    HOST_ADDITIVES = HOST
VIRUSES = config["viruses"]
VIRUSES = VIRUSES.replace(" ", "")
REPEATS_GTF = join(FASTAS_GTFS_DIR, HOST + ".repeats.gtf")

HOST_ADDITIVES_VIRUSES = HOST_ADDITIVES + "," + VIRUSES
HOST_VIRUSES = HOST + "," + VIRUSES
HOST_ADDITIVES_VIRUSES = HOST_ADDITIVES_VIRUSES.split(",")
FASTAS = [join(FASTAS_GTFS_DIR, f + ".fa") for f in HOST_ADDITIVES_VIRUSES]
REGIONS = [join(FASTAS_GTFS_DIR, f + ".fa.regions") for f in HOST_ADDITIVES_VIRUSES]
REGIONS_HOST = [join(FASTAS_GTFS_DIR, f + ".fa.regions") for f in HOST.split(",")]
REGIONS_VIRUSES = [join(FASTAS_GTFS_DIR, f + ".fa.regions") for f in VIRUSES.split(",")]
GTFS = [join(FASTAS_GTFS_DIR, f + ".gtf") for f in HOST_ADDITIVES_VIRUSES]
FASTAS_REGIONS_GTFS = FASTAS.copy()
FASTAS_REGIONS_GTFS.extend(REGIONS)
FASTAS_REGIONS_GTFS.extend(GTFS)

ANNOTATION_LOOKUP = config["annotation_lookups"][HOST]

if not os.path.exists(join(WORKDIR, "fastqs")):
    os.mkdir(join(WORKDIR, "fastqs"))
if not os.path.exists(join(WORKDIR, "results")):
    os.mkdir(join(WORKDIR, "results"))

REQUIRED_FILES = [config[f] for f in ["samples", "tools", "cluster"]]
REQUIRED_FILES.append(ANNOTATION_LOOKUP)
REQUIRED_FILES.extend(FASTAS)
REQUIRED_FILES.extend(REGIONS)
REQUIRED_FILES.extend(GTFS)
for f in REQUIRED_FILES:
    check_readaccess(f)

REF_FA = join(REF_DIR, "ref.fa")
REF_REGIONS = join(REF_DIR, "ref.fa.regions")
REF_REGIONS_HOST = join(REF_DIR, "ref.fa.regions.host")
REF_REGIONS_VIRUSES = join(REF_DIR, "ref.fa.regions.viruses")
REF_GTF = join(REF_DIR, "ref.gtf")
append_files_in_list(FASTAS, REF_FA)
append_files_in_list(REGIONS, REF_REGIONS)
append_files_in_list(REGIONS_HOST, REF_REGIONS_HOST)
append_files_in_list(REGIONS_VIRUSES, REF_REGIONS_VIRUSES)
append_files_in_list(GTFS, REF_GTF)

SAMPLESDF = pd.read_csv(config["samples"], sep="\t", header=0)
# NCLscan does not like hyphens in sample names change them to underscores
# see https://github.com/CCBR/CCBR_circRNA_DAQ/issues/53
SAMPLESDF["sampleName"] = SAMPLESDF["sampleName"].str.replace("-", "_")
SAMPLESDF.set_index(["sampleName"], inplace=True)
SAMPLES = list(SAMPLESDF.index)
SAMPLESDF["R1"] = join(RESOURCES_DIR, "dummy")
SAMPLESDF["R2"] = join(RESOURCES_DIR, "dummy")
SAMPLESDF["PEorSE"] = "PE"

for sample in SAMPLES:
    R1file = SAMPLESDF["path_to_R1_fastq"][sample]
    R2file = SAMPLESDF["path_to_R2_fastq"][sample]
    # print(sample,R1file,R2file)
    check_readaccess(R1file)
    R1filenewname = join(WORKDIR, "fastqs", sample + ".R1.fastq.gz")
    if not os.path.exists(R1filenewname):
        os.symlink(R1file, R1filenewname)
        # os.link(R1file,R1filenewname)
    SAMPLESDF.loc[[sample], "R1"] = R1filenewname
    if str(R2file) != "nan":
        check_readaccess(R2file)
        R2filenewname = join(WORKDIR, "fastqs", sample + ".R2.fastq.gz")
        if not os.path.exists(R2filenewname):
            os.symlink(R2file, R2filenewname)
            # os.link(R2file,R2filenewname)
        SAMPLESDF.loc[[sample], "R2"] = R2filenewname
    else:
        SAMPLESDF.loc[[sample], "PEorSE"] = "SE"
# print(SAMPLESDF)
# sys.exit()

## Load tools from YAML file
with open(config["tools"]) as f:
    TOOLS = yaml.safe_load(f)

## Load cluster.json
with open(config["cluster"]) as json_file:
    CLUSTER = yaml.safe_load(json_file)

## Create lambda functions to allow a way to insert read-in values
## as rule directives
getthreads = (
    lambda rname: int(CLUSTER[rname]["threads"])
    if rname in CLUSTER and "threads" in CLUSTER[rname]
    else int(CLUSTER["__default__"]["threads"])
)
getmemg = (
    lambda rname: CLUSTER[rname]["mem"]
    if rname in CLUSTER and "mem" in CLUSTER[rname]
    else CLUSTER["__default__"]["mem"]
)
getmemG = lambda rname: getmemg(rname).replace("g", "G")

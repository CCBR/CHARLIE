## functions

def get_file_size(filename):
    filename=filename.strip()
    if check_readaccess(filename):
        return os.stat(filename).st_size

def get_peorse(wildcards):
    return SAMPLESDF["PEorSE"][wildcards.sample]

def check_existence(filename):
    """Checks if file exists on filesystem
    :param filename <str>: Name of file to check
    """
    filename=filename.strip()
    if not os.path.exists(filename):
        sys.exit("File: {} does not exists!".format(filename))
    return True


def check_readaccess(filename):
    """Checks permissions to see if user can read a file
    :param filename <str>: Name of file to check
    """
    filename=filename.strip()
    check_existence(filename)
    if not os.access(filename,os.R_OK):
        sys.exit("File: {} exists, but user cannot read from file due to permissions!".format(filename))
    return True


def check_writeaccess(filename):
    """Checks permissions to see if user can write to a file
    :param filename <str>: Name of file to check
    """
    filename=filename.strip()
    check_existence(filename)
    if not os.access(filename,os.W_OK):
        sys.exit("File: {} exists, but user cannot write to file due to permissions!".format(filename))
    return True

def append_files_in_list(flist,ofile):
    if not os.path.exists(ofile):
        with open(ofile, 'w') as outfile:
            for fname in flist:
                with open(fname) as infile:
                    outfile.write(infile.read())
    return True

##### load config and sample sheets #####

# check_readaccess("config/config.yaml")
# configfile: "config/config.yaml"

## set memory limit 
## used for sambamba sort, etc
# MEMORYG="100G"

def _is_true(variable):
    if variable==True or variable=="True" or variable == "TRUE":
        return True
    else:
        return False

#resouce absolute path
WORKDIR=config['workdir']
SCRIPTS_DIR=config['scriptsdir']
RESOURCES_DIR=config['resourcesdir']
FASTAS_GTFS_DIR=config['fastas_gts_dir']
RUN_CLEAR=_is_true(config['run_clear'])
RUN_DCC=_is_true(config['run_dcc'])
RUN_MAPSPLICE=_is_true(config['run_mapsplice'])
RUN_NCLSCAN=_is_true(config['run_nclscan'])


REF_DIR=join(WORKDIR,"ref")
if not os.path.exists(REF_DIR):
    os.mkdir(REF_DIR)
STAR_INDEX_DIR=join(REF_DIR,"STAR_no_GTF")
if not os.path.exists(STAR_INDEX_DIR):
    os.mkdir(STAR_INDEX_DIR)
# strip trailing slashes if any
for d in [WORKDIR,SCRIPTS_DIR,RESOURCES_DIR,FASTAS_GTFS_DIR,STAR_INDEX_DIR,REF_DIR]:
    d=d.strip('r\/')
BWA_INDEX=join(REF_DIR,"ref")

HOST=config['host'] # hg38 or mm39
ADDITIVES=config['additives'] # ERCC and/or BAC16Insert
if ADDITIVES != "":
    HOST_ADDITIVES=HOST+","+ADDITIVES
else:
    HOST_ADDITIVES=HOST
VIRUSES=config['viruses']
REPEATS_GTF=join(FASTAS_GTFS_DIR,HOST+".repeats.gtf")

HOST_ADDITIVES_VIRUSES=HOST_ADDITIVES+","+VIRUSES
HOST_ADDITIVES_VIRUSES=HOST_ADDITIVES_VIRUSES.split(",")
FASTAS=[join(FASTAS_GTFS_DIR,f+".fa") for f in HOST_ADDITIVES_VIRUSES]
REGIONS=[join(FASTAS_GTFS_DIR,f+".fa.regions") for f in HOST_ADDITIVES_VIRUSES]
GTFS=[join(FASTAS_GTFS_DIR,f+".gtf") for f in HOST_ADDITIVES_VIRUSES]
FASTAS_REGIONS_GTFS=FASTAS.copy()
FASTAS_REGIONS_GTFS.extend(REGIONS)
FASTAS_REGIONS_GTFS.extend(GTFS)

ANNOTATION_LOOKUP=config['annotation_lookups'][HOST]

if not os.path.exists(join(WORKDIR,"fastqs")):
    os.mkdir(join(WORKDIR,"fastqs"))
if not os.path.exists(join(WORKDIR,"results")):
    os.mkdir(join(WORKDIR,"results"))

REQUIRED_FILES=[config[f] for f in ["samples", "tools", "cluster"]]
REQUIRED_FILES.append(ANNOTATION_LOOKUP)
REQUIRED_FILES.extend(FASTAS)
REQUIRED_FILES.extend(REGIONS)
REQUIRED_FILES.extend(GTFS)
for f in REQUIRED_FILES:
    check_readaccess(f)

REF_FA=join(REF_DIR,"ref.fa")
REF_REGIONS=join(REF_DIR,"ref.fa.regions")
REF_GTF=join(REF_DIR,"ref.gtf")
append_files_in_list(FASTAS,REF_FA)
append_files_in_list(REGIONS,REF_REGIONS)
append_files_in_list(GTFS,REF_GTF)

SAMPLESDF = pd.read_csv(config["samples"],sep="\t",header=0,index_col="sampleName")
SAMPLES = list(SAMPLESDF.index)
SAMPLESDF["R1"]=join(RESOURCES_DIR,"dummy")
SAMPLESDF["R2"]=join(RESOURCES_DIR,"dummy")
SAMPLESDF["PEorSE"]="PE"

for sample in SAMPLES:
    R1file=SAMPLESDF["path_to_R1_fastq"][sample]
    R2file=SAMPLESDF["path_to_R2_fastq"][sample]
    # print(sample,R1file,R2file)
    check_readaccess(R1file)
    R1filenewname=join(WORKDIR,"fastqs",sample+".R1.fastq.gz")
    if not os.path.exists(R1filenewname):
        os.symlink(R1file,R1filenewname)
        # os.link(R1file,R1filenewname)
    SAMPLESDF.loc[[sample],"R1"]=R1filenewname
    if str(R2file)!='nan':
        check_readaccess(R2file)
        R2filenewname=join(WORKDIR,"fastqs",sample+".R2.fastq.gz")
        if not os.path.exists(R2filenewname):
            os.symlink(R2file,R2filenewname)
            # os.link(R2file,R2filenewname)
        SAMPLESDF.loc[[sample],"R2"]=R2filenewname
    else:
        SAMPLESDF.loc[[sample],"PEorSE"]="SE"
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
getthreads=lambda rname:int(CLUSTER[rname]["threads"]) if rname in CLUSTER and "threads" in CLUSTER[rname] else int(CLUSTER["__default__"]["threads"])
getmemg=lambda rname:CLUSTER[rname]["mem"] if rname in CLUSTER and "mem" in CLUSTER[rname] else CLUSTER["__default__"]["mem"]
getmemG=lambda rname:getmemg(rname).replace("g","G")




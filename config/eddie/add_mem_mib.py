#!/usr/bin/env python

import ruamel.yaml

yaml = ruamel.yaml.YAML()
yaml.preserve_quotes = True
yaml.explicit_start = True

with open("cluster.yaml", "r") as infile:
    data = yaml.load(infile)

for k, v in data.items():
    if "mem" in v:
        data[k]["mem_mib"] = int(v["mem"].rstrip("g")) * 1024

with open("cluster.yaml.2", "w") as outfile:
    yaml.dump(data, outfile)

import json
import pandas as pd
import sys

from pprint import pprint


with open(sys.argv[1], "r") as file:
    data_list = [json.loads(line) for line in file]

for d in data_list:
    if d["clust_name"] == 'YacG':
        pprint(d)

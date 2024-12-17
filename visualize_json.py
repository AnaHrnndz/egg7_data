import json
import pandas as pd

from pprint import pprint


with open("/data/projects/Egg7_data/process_results/output_ogs_with_annotations.jsonl", "r") as file:
    data_list = [json.loads(line) for line in file]
pprint(data_list)

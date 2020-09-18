import re
import pandas as pd

FIELD_DTYPES = {
    "segmentation": int,
    "id": int, "valid":int,
    "pre_pt_position": "object",
    "pre_pt_supervoxel_id": int,
    "pre_pt_root_id": int,
    "ctr_pt_position": "object",
    "post_pt_position": "object",
    "post_pt_supervoxel_id": int,
    "post_pt_root_id": int,
    "size": int, "spine_vol": float,
    "exclude_conn": int
}
COORDS_REGEXP = re.compile(r"\[([0-9]+) +([0-9]+) +([0-9]+)\]")


def read_pyc_graph(filename):
    rawdf = pd.read_csv(filename, dtype=FIELD_DTYPES, index_col=0)
    
    coordinatecols = [
        "pre_pt_position",
        "ctr_pt_position",
        "post_pt_position"
    ]
    
    for colname in coordinatecols:
        rawdf[colname] = [evalcoord(v) for v in rawdf[colname]]

    return rawdf


def evalcoord(valstr):
    return list(COORDS_REGEXP.match(valstr).groups())
    
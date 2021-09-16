"""
Computing compartment path length

Measuring the amount of path length covered by each compartment
"""
import sys
import argparse

import pandas as pd
from tqdm import tqdm

sys.path.append(".")
from lib import skel, u


ENGLISHLABELS = ["soma", "axon", "basal", "apical",
                 "ambiguous dendrite", "ambiguous"]


def main(cellid_filename, output_filename):
    cellids = pd.read_csv(cellid_filename, index_col=0).pt_root_id

    # builds a dataframe by records
    records = list()
    for cellid in tqdm(cellids):
        skl = skel.read_neuron_skel(cellid)
        lbl = skel.read_neuron_complbls(cellid)
        record = dict(cellid=cellid)
        for compartment in range(6):
            englbl = ENGLISHLABELS[compartment]
            length = skel.compartment_pathlength(skl, lbl, compartment)
            record[f"{englbl}_length_um"] = length
        records.append(record)

    pd.DataFrame.from_records(records).to_csv(output_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("cellid_filename")
    parser.add_argument("output_filename")

    main(**vars(parser.parse_args()))

"""
Compute OSi, DSi, and determine significantly tuned cells.
"""
import numpy as np
import pandas as pd


def main():

	id_fname = "../data/pyc_func_id_list.csv"
	id_df = pd.read_csv(id_fname)

	seg_id_list = np.array(id_df["segment_id"])
	scan_list = np.array(id_df["scan"])
	n = seg_id_list.shape[0]


if __name__ == "__main__":

	main()
"""
Compute axon/dendrite length [um] using skeletons.
"""
import numpy as np
import pandas as pd
import argparse

from meshparty import skeleton_io


def main():

	id_fname = "../data/pyc_id_list.csv"
	id_df = pd.read_csv(id_fname)
	seg_id_list = np.array(id_df["segment_id"])
	n = seg_id_list.shape[0]

	dendrite_length = np.zeros(n)
	axon_length = np.zeros(n)
	print("Computing axon/dendrite length...")
	for i in range(n):

		seg_id = seg_id_list[i]

		seg_skel, skel_lab = load_skeleton(seg_id)

		dendrite_length[i] = compute_path_len(seg_skel, skel_lab, "dendrite")
		axon_length[i] = compute_path_len(seg_skel, skel_lab, "axon")

		if (i+1)%100==0 or i==n-1:
			print("{} / {} completed.".format(i+1, n))

	write_data(seg_id_list, dendrite_length, axon_length)
	print("Results saved in ../data/neurite_length.csv")


def load_skeleton(seg_id):
	"""
	Load skeleton.
	"""
	path_skel = "../data/smoothed_skeletons_v185/"

	skel_lab = np.load(path_skel+str(seg_id)+"_skeleton_label.npy")
	seg_skel = skeleton_io.read_skeleton_h5(path_skel+str(seg_id)+"_skeleton.h5")

	return seg_skel, skel_lab


def compute_path_len(skel, skel_label, neurite_type="dendrite"):
	"""
	Compute axon/dendrite length [um].
	"""
	
	g = skel.csgraph

	if neurite_type=="dendrite":
		valid_idx = np.where((skel_label==2)+(skel_label==3)+(skel_label==4))[0]

	elif neurite_type=="axon":
		valid_idx = np.where((skel_label==1))[0]

	g = g[valid_idx][:,valid_idx]

	return np.sum(g)/1000


def write_data(seg_ids, dendrite_length, axon_length):
	"""
	Save axon/dendrite length data.
	"""

	df = pd.DataFrame(data={"seg_id": seg_ids,
													"axon_len": axon_length,
													"dendrite_len": dendrite_length})
	output_fname = "../data/neurite_length.csv"
	df.to_csv(output_fname, index=False)


if __name__ == "__main__":

	main()

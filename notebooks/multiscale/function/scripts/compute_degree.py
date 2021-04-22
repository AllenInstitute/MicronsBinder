"""
Compute in-connection degree, in-synapse degree.
"""
import numpy as np
import pandas as pd


def main():

	id_fname = "../data/pyc_id_list.csv"
	id_df = pd.read_csv(id_fname)
	seg_id_list = np.array(id_df["segment_id"])
	n = seg_id_list.shape[0]

	conn_fname = "../data/pyc_pyc_subgraph.csv"
	conn_df = pd.read_csv(conn_fname)

	syn_arr, conn_arr = get_connectivity_array(conn_df)

	outsyn_deg = np.zeros(n, dtype="int")
	insyn_deg = np.zeros(n, dtype="int")
	outconn_deg = np.zeros(n, dtype="int")
	inconn_deg = np.zeros(n, dtype="int")
	print("Computing degree...")
	for i in range(n):

		seg_id = seg_id_list[i]

		outsyn_deg[i] = compute_deg(syn_arr, seg_id, "out")
		insyn_deg[i] = compute_deg(syn_arr, seg_id, "in")
		outconn_deg[i] = compute_deg(conn_arr, seg_id, "out")
		inconn_deg[i] = compute_deg(conn_arr, seg_id, "in")

		if (i+1)%100==0 or i==n-1:
			print("{} / {} completed.".format(i+1, n))
 	
	write_data(seg_id_list, outsyn_deg, insyn_deg, outconn_deg, inconn_deg)
 	

def get_connectivity_array(conn_df):

	pre_id_list = np.array(conn_df["pre_root_id"])
	post_id_list = np.array(conn_df["post_root_id"])

	n = pre_id_list.shape[0]
	# Connectivity array with all synapses
	synapse_arr = np.zeros((n,2), dtype="uint64")
	synapse_arr[:,0] = pre_id_list
	synapse_arr[:,1] = post_id_list

	# Connectivity array with all cell pairs
	connection_arr = np.unique(synapse_arr, axis=0)

	return synapse_arr, connection_arr


def compute_deg(arr, seg_id, direction="in"):

	if direction=="in":
		deg = np.sum(arr[:,1]==seg_id)

	elif direction=="out":
		deg = np.sum(arr[:,0]==seg_id)

	return deg


def write_data(seg_ids, outsyn_deg, insyn_deg, outconn_deg, inconn_deg):
	"""
	Save out/in-degree data.
	"""

	df = pd.DataFrame(data={"seg_id": seg_ids,
													"outsyn_deg": outsyn_deg,
													"insyn_deg": insyn_deg,
													"outconn_deg": outconn_deg,
													"inconn_deg": inconn_deg})
	output_fname = "../data/degree.csv"
	df.to_csv(output_fname, index=False)
	print("Results saved in {}".format(output_fname))


if __name__ == "__main__":

	main()
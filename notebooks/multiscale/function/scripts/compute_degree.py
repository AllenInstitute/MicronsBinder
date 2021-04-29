"""
Compute in-connection degree, in-synapse degree.
"""
import numpy as np
import pandas as pd
import pickle


def main():

	Neuron = load_dict("../data/Neuron.pkl")
	seg_id_list = Neuron["segment_id"]
	n = seg_id_list.shape[0]

	conn_fname = "../data/pyc_pyc_subgraph.csv"
	conn_all_fname = "../data/pyc_all_synapse.csv"
	
	conn_df = pd.read_csv(conn_fname)
	conn_all_df = pd.read_csv(conn_all_fname)

	conn_post_df = filter_conn_df(conn_all_df, seg_id_list)

	syn_arr, conn_arr = get_connectivity_array(conn_df)
	syn_all_arr, conn_all_arr = get_connectivity_array(conn_post_df)

	outsyn_deg = np.zeros(n, dtype="int")
	insyn_deg = np.zeros(n, dtype="int")
	outconn_deg = np.zeros(n, dtype="int")
	inconn_deg = np.zeros(n, dtype="int")
	total_insyn_deg = np.zeros(n, dtype="int")
	print("Computing degree...")
	for i in range(n):

		seg_id = seg_id_list[i]

		outsyn_deg[i] = compute_deg(syn_arr, seg_id, "out")
		insyn_deg[i] = compute_deg(syn_arr, seg_id, "in")
		outconn_deg[i] = compute_deg(conn_arr, seg_id, "out")
		inconn_deg[i] = compute_deg(conn_arr, seg_id, "in")
		total_insyn_deg[i] = compute_deg(syn_all_arr, seg_id, "in")

		if (i+1)%100==0 or i==n-1:
			print("{} / {} completed.".format(i+1, n))
 	
	write_data(seg_id_list, outsyn_deg, insyn_deg, outconn_deg, inconn_deg, total_insyn_deg)


def get_soma_loc(seg_ids):

	Soma = load_dict("../data/Soma.pkl")

	n = seg_ids.shape[0]
	soma_loc = np.zeros((n,3))
	for i in range(n):

		seg_id = seg_ids[i]
		soma_loc[i,:] = Soma["loc"][Soma["segment_id"]==seg_id][0]

	soma_loc = soma_loc/1000

	return soma_loc


def filter_conn_df(conn_all_df, seg_ids):

	post_id_all = np.array(conn_all_df["post_root_id"], dtype="uint64")
	post_idx = np.where(np.isin(post_id_all, seg_ids))[0]

	conn_post_df = conn_all_df.iloc[post_idx]

	post_centroids = np.zeros((post_idx.shape[0],3))
	for i in range(post_idx.shape[0]):
		post_centroids[i,:] = str2coord(conn_post_df.iloc[i]["ctr_pt_position"])

	soma_loc = get_soma_loc(seg_ids)

	post_id_list = np.array(conn_post_df["post_root_id"])
	post_soma_loc = np.zeros(post_centroids.shape)
	for i in range(post_id_list.shape[0]):
		post_soma_loc[i,:] = soma_loc[seg_ids==post_id_list[i],:]

	res = np.array([0.004, 0.004, 0.04])
	post_centroids = post_centroids*res

	d_post = distance(post_centroids, post_soma_loc)
	valid = (d_post>15)

	return conn_post_df.iloc[valid]


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


def distance(x, y):

	return np.sqrt(np.sum((x-y)**2, axis=1))


def str2coord(str_coord, dtype="int"):
    
  l = str_coord[1:-1].split(" ")
  coord = []
  for j in range(len(l)):
    if l[j] != "":
      if dtype == "int": 
        coord.append(int(l[j]))
      elif dtype == "float":
        coord.append(float(l[j]))
    
  coord = np.array(coord)
    
  return coord


def load_dict(fname):

  with open(fname, 'rb') as handle:
    data = pickle.load(handle)

  return data


def write_data(seg_ids, outsyn_deg, insyn_deg, outconn_deg, inconn_deg, total_insyn_deg):
	"""
	Save out/in-degree data.
	"""

	df = pd.DataFrame(data={"segment_id": seg_ids,
													"outsyn_deg": outsyn_deg,
													"insyn_deg": insyn_deg,
													"outconn_deg": outconn_deg,
													"inconn_deg": inconn_deg,
													"total_insyn_deg": total_insyn_deg})

	output_fname = "../data/SynDegree.csv"
	df.to_csv(output_fname, index=False)
	print("Results saved in {}".format(output_fname))


if __name__ == "__main__":

	main()
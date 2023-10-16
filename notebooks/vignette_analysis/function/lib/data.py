"""
Helper functions for fetching data from data tables.
"""
import numpy as np


# Get soma center coordinates ([um])
def get_soma_loc(data_dict, seg_id):

	seg_ids = data_dict["segment_id"]
	soma_locs = data_dict["loc"]

	return soma_locs[seg_ids==seg_id][0]


# get stimulus label
def get_stim_label(data_dict, scan_id):

	scan_ids = data_dict["scan_id"]
	conditions = data_dict["condition"]

	return np.round(conditions[scan_ids==scan_id][0],1)


# Get calcium trace
def get_trace(data_dict, seg_id, scan_id, trace_type="spike"):

	seg_ids = data_dict["segment_id"]
	scan_ids = data_dict["scan_id"]
	traces = data_dict[trace_type]

	valid = (seg_ids==seg_id)*(scan_ids==scan_id)

	return traces[np.where(valid)[0][0]]


# Get tuning curve
def get_tuning(data_dict, seg_id, scan_id):

	seg_ids = data_dict["segment_id"]
	scan_ids = data_dict["scan_id"]
	tune = data_dict["tune_curve"]

	valid = (seg_ids==seg_id)*(scan_ids==scan_id)

	return tune[np.where(valid)[0][0]]


# Get scan id
def get_scan(data_dict, seg_id):

	seg_ids = data_dict["segment_id"]
	scan_ids = data_dict["scan_id"]

	return scan_ids[seg_ids==seg_id][0]


# Get synapse_conneciton density
def get_density(data_dict, seg_id, dens_type="inconn_dens"):

	seg_ids = data_dict["segment_id"]
	densities = data_dict[dens_type]

	return densities[seg_ids==seg_id][0]
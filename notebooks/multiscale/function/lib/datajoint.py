"""
Helper functions for fetching data from DataJoint database.
"""
import numpy as np


# Get soma center coordinates ([nm])
def get_soma_loc(database, seg_id):
    
  key = {"segmentation": 185, "segment_id": seg_id}
    
  return (database.Soma() & key).fetch1("loc")/1000


# Get stimulus label
def get_stim_label(database, scan_id):
    
  key = {"scan_id": scan_id}
    
  return np.round((database.Stimulus() & key).fetch1("condition"), 1)


# Get calcium trace
def get_trace(database, seg_id, scan_id, trace_type="spike"):
    
  trace = (database.EASETrace() & {"segment_id": seg_id, "scan_id": scan_id}).fetch1(trace_type)
    
  return trace


# Get tuning curve
def get_tuning(database, seg_id, scan_id):

	tune = (database.EASETuning() & {"segment_id": seg_id, "scan_id": scan_id}).fetch1("tune_curve")

	return tune


# Get scan id
def get_scan(database, seg_id):

	scan_id = (database.EASETuning() & {"segment_id": seg_id}).fetch1("scan_id")

	return scan_id


# Get synapse/connection density
def get_density(database, seg_id, dens_type="inconn_dens"):

	density = (database.SynDensity() & {"segment_id": seg_id}).fetch1(dens_type)

	return density
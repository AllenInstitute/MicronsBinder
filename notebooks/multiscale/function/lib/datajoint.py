# DataJoint
# Get soma center coordinates ([nm])
def get_soma_loc(database, seg_id):
    
  key = {"segmentation": 185, "segment_id": seg_id}
    
  return (database.Soma() & key).fetch1("loc")/1000


# Get stimulus label
def get_stim_label(database, scan_id):
    
  key = {"scan_id": scan_id}
    
  return (database.Stimulus() & key).fetch1("condition")


# Get functional responses
def get_trace(database, seg_id, scan_id, trace_type="spike"):
    
  trace = (database.EASETrace() & {"segment_id": seg_id, "scan_id": scan_id}).fetch1(trace_type)
    
  return trace
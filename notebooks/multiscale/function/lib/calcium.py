"""
Visual response related functions.
"""
import numpy as np


def get_section(conditions, angle):
    
  if np.isnan(angle):
    valid = np.isnan(conditions).astype("int")
  else:
    valid = (conditions==angle).astype("int")
    
  valid_diff = np.diff(valid)
  valid_diff = valid_diff.reshape((1,-1))
  
  st_idx = np.where(valid_diff==1)[1] + 1
  end_idx = np.where(valid_diff==-1)[1]
    
  if st_idx.shape[0] > end_idx.shape[0]:
    end_idx = np.concatenate((end_idx,np.array([len(conditions)])))
  elif st_idx.shape[0] < end_idx.shape[0]:
    st_idx = np.concatenate((np.array([0]),st_idx))
  elif st_idx[0] > end_idx[0]:
    st_idx = np.concatenate((np.array([0]),st_idx))
    end_idx = np.concatenate((end_idx,np.array([len(conditions)])))
    
  section_list = []
  for i in range(st_idx.shape[0]):
    section_list.append((st_idx[i], end_idx[i]+1))
    
  return section_list
    

def get_peakamp_tdarray(trace, condition):

  valid = ~np.isnan(condition)
  angle_list = np.unique(condition[valid])

  tdarray = np.zeros((30,16))
  for i in range(angle_list.shape[0]):

    angle = angle_list[i]
    section_list = get_section(condition, angle)

    offset = 0
    for j in range(len(section_list)):

      if len(section_list)!=30:
	      offset = 30-len(section_list)

      s = section_list[j] 
      trace_sect = trace[s[0]:s[1]]

      max_idx = np.argmax(trace_sect)
      max_val = trace_sect[max_idx]
      if (max_idx==0):
        tdarray[j+offset,i] = 0      
      elif (max_idx==trace_sect.shape[0]-1):
        trace_post = trace[s[1]:s[1]+15]
        tdarray[j+offset,i] = np.max(trace_post)
      elif trace_sect[0]>0.5*max_val:
        tdarray[j+offset,i] = 0
      else:
        tdarray[j+offset,i] = max_val

  return tdarray


def get_active_tdarray(spike, condition, thr=3):
    
  spike = spike*(spike>=thr*np.std(spike))
  
  valid = ~np.isnan(condition)
  angle_list = np.unique(condition[valid])
  
  tdarray = np.zeros((30,16))
  for i in range(angle_list.shape[0]):
        
    angle = angle_list[i]
    section_list = get_section(condition, angle)

    for j in range(len(section_list)):

      s = section_list[j]

      spike_section = spike[s[0]:s[1]]
      if np.sum(spike_section) > 0:
        tdarray[j,i] = 1
      else:
        tdarray[j,i] = 0
                
  return tdarray


def compute_intermittency(spike, condition, pref_idx, thr=3):
    
  spike = spike*(spike>=thr*np.std(spike))

  valid = ~np.isnan(condition)
  angle_list = np.unique(condition[valid])

  trial_arr = np.ones((16,30))*np.nan
  for i in range(angle_list.shape[0]):
        
    angle = angle_list[i]
    section_list = get_section(condition, angle)

    for j in range(len(section_list)):

      s = section_list[j]

      spike_section = spike[s[0]:s[1]]
      if np.sum(spike_section) > 0:
        trial_arr[i,j] = 0
      else:
        trial_arr[i,j] = 1
    
  frac_angle = np.nansum(trial_arr, axis=1)/np.sum(~np.isnan(trial_arr), axis=1) 
    
  return np.mean(frac_angle[pref_idx])
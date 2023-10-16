"""
Tuning curve related functions.
"""
import numpy as np
from scipy.stats import percentileofscore


def tuning_curve(tdarray):
    
  u, s, vh = np.linalg.svd(tdarray, full_matrices=False)
  tune = np.abs(vh[0,:])
  
  return tune


def dsi(tune):
    
  tune = tune.reshape((-1,1))
  
  angles = np.linspace(0, 2*np.pi, 17)[:-1]
  vec = np.zeros((1,16), dtype="complex")
  for i in range(16):
    vec[0,i] = np.exp(complex(0,angles[i]))
  K = vec@tune/np.sum(tune)
    
  return np.abs(K)[0]


def osi(tune):
    
  tune = tune.reshape((-1,1))
  angles = np.linspace(0, 4*np.pi, 17)[:-1]
  vec = np.zeros((1,16), dtype="complex")
  for i in range(16):
    vec[0,i] = np.exp(complex(0,angles[i]))
  K = vec@tune/np.sum(tune)
    
  return np.abs(K)[0]


# Permutation test
def shuffle_amparr(amp_arr):
    
  amp_arr_copy = np.copy(amp_arr)
  for i in range(amp_arr.shape[0]):
        
    tmp = amp_arr_copy[i,:]
    np.random.shuffle(tmp)
      
  return amp_arr_copy


def permutation_test(amp_arr, n_shuff, mode="dsi"):
    
  si_shuffled = np.zeros(n_shuff)
    
  if mode=="dsi":
    
    tune_true = tuning_curve(amp_arr)    
    si_true = dsi(tune_true)
    for t in range(n_shuff):

      amp_arr_shuffled = shuffle_amparr(amp_arr)
      tune_shuffled = tuning_curve(amp_arr_shuffled)
      si_shuffled[t] = dsi(tune_shuffled)
        
  elif mode=="osi":
    
    tune_true = tuning_curve(amp_arr)
    si_true = osi(tune_true)
    for t in range(n_shuff):

      amp_arr_shuffled = shuffle_amparr(amp_arr)
      tune_shuffled = tuning_curve(amp_arr_shuffled)
      si_shuffled[t] = osi(tune_shuffled)
            
  p = percentileofscore(si_shuffled, si_true)
  p = 1-p/100
        
  return si_shuffled, p
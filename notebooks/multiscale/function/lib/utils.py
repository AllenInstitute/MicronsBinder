# Utils
import numpy as np


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
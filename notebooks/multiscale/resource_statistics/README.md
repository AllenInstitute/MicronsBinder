# Resource statistics computation

* `neurite_lengths.ipynb` 
  Plotting axon and dendrite branch lengths for each cell with clean compartment labels and a soma in the volume

### Base data

* `data/soma_ids/p100_pyr_soma_IDs_v185.csv`  
  A list of PyC ids with their soma captured by the volume for proofreading version 185
* `data/soma_ids/p100_inh_soma_IDs_v185.csv`  
  A list of interneuron ids with their soma captured by the volume for proofreading version 185
* `data/clean_compartment_ids_v185.csv`  
  A list of ids for which the skeleton compartment labels are sufficiently correct
* `data/smoothed_skels`  
  A set of smoothed skeletons for all cells with their soma captured by the volume. These were computed using MeshParty

### Intermediate data

* `data/pyr_neurite_lengths.csv`
  A dataframe that describes the length of each branch for every PyC
* `data/inh_neurite_lengths.csv`
  A dataframe that describes the length of each branch for every putative interneuron

### Scripts

* `scripts/compute_neurite_lengths.py`
  Computes the length of each branch for each PyC and putative interneuron from the smoothed skeletons

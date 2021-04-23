# Compute neurite lengths for each cell id in data/soma_ids
# -> data/pyr_dist_to_leaves.csv
# -> data/inh_dist_to_leaves.csv
python scripts/compute_dist_to_leaf.py

# Compute total path length covered by each compartment label
# -> data/pyr_compartment_lengths.csv
# -> data/inh_compartment_lengths.csv
python scripts/compute_compartment_lengths.py data/soma_ids/p100_pyr_soma_IDs_v185.csv
python scripts/compute_compartment_lengths.py data/soma_ids/p100_inh_soma_IDs_v185.csv

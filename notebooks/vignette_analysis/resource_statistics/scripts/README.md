# Computing intermediate data

These scripts are designed to be run from the `resource_statistics` directory above this one.

## Computing distance to each leaf of the skeleton
```
$ python scripts/compute_dist_to_leaf.py
```

## Computing the path length covered by each compartment
```
$ python scripts/compute_compartment_lengths.py data/soma_ids/p100_pyr_soma_IDs_v185.csv data/pyr_compartment_lengths.csv
$ python scripts/compute_compartment_lengths.py data/soma_ids/p100_inh_soma_IDs_v185.csv data/inh_compartment_lengths.csv
```

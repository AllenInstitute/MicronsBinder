#!/bin/bash

# Computing basic branch statistics
# NOTE: this requires downloading the meshes
#   (through downloadNeuronMeshes.sh)
python scripts/computebasicbranchstats.py \
    data/clean_and_complete_soma_ids_v185.csv \
    data/basicbranchstats.csv


# Associating mitochondria to skeleton nodes
# NOTE: this requires downloading the mitochondria meshes
#  (through downloadmitomeshes.py)
python scripts/mitotoskel.py \
    data/pni_mito_cellswskel_v185.csv \
    data/temp/assoc_
python scripts/assemblemitotoskel.py \
    data/pni_mito_cellswskel_v185.csv \
    data/pycids_v185.csv data/mitotoskel.h5

# Computing MCI and assigning mitochondria to compartments
# NOTE: this requires downloading the mitochondria meshes
#  (through downloadmitomeshes.py)
python scripts/computeextramitostats.py \
    data/pni_mito_cellswskel_v185.csv \
    data/mitotoskel.h5 \
    data/pni_mito_cellswskel_v185_fullstats.csv


# Computing density statistics
python scripts/computedensitystats.py \
    data/pni_mito_cellswskel_v185_fullstats.csv \
    data/pycids_v185.csv data/neuron_received_synapses_v185.csv \
    data/basicbranchstats.csv \
    data/pni_mito_cell_overlaps_v185.csv \
    data/pni_nucleus_segments.csv \
    data/clean_and_complete_nucleus_lookup.csv \
    data/fullbranchstats.csv \
    data/cellwisestats.csv

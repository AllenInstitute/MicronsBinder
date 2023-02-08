#!/bin/bash

# Computing basic branch statistics
# NOTE: this requires downloading the meshes
#   (through downloadNeuronMeshes.sh)
python scripts/computebasicbranchstats.py \
    data/pyc_analysis_ids.csv \
    data/distancebinstats.csv \
    --branchtype distancebins


# Associating mitochondria to skeleton nodes
# NOTE: this requires downloading the mitochondria meshes
#  (through downloadmitomeshes.py)
python scripts/mitotoskel.py \
    data/pni_mito_cellswskel_v185.csv \
    data/pni_analysis_ids.csv \
    data/temp/assoc_
python scripts/assemblemitotoskel.py \
    data/pni_mito_cellswskel_v185.csv \
    data/pni_analysis_ids.csv \
    data/mitotoskel.h5


# Computing MCI and assigning mitochondria to compartments
# NOTE: this requires downloading the mitochondria meshes
#  (through downloadmitomeshes.py)
python scripts/computeextramitostats.py \
    data/pni_mito_cellswskel_v185.csv \
    data/mitotoskel_v185.h5 \
    data/pni_mito_analysisids_v185_extrastats.csv


# Computing density statistics
python scripts/computedensitystats.py \
    data/pni_mito_analysisids_v185_extrastats.csv \
    data/pni_analysis_ids.csv data/neuron_received_synapses_v185.csv \
    data/distancebinstats.csv \
    data/pni_mito_cell_overlaps_v185.csv \
    data/pni_nucleus_segments.csv \
    data/pni_nucleus_cell_overlaps_v185.csv \
    data/distancebinfullstats.csv \
    data/notused.csv \
    --branchtype distancebins


# Computing branch diameter
python scripts/branchdiameter.py \
    data/pni_analysis_ids.csv \
    data/mitotoskel_v185.h5 \
    data/distancebindiameters.csv \
    --branchtype distancebins


# Computing segment distance to soma
python scripts/branchproximity.py \
    data/pni_analysis_ids.csv \
    data/distancebinproximity.csv \
    --branchtype distancebins


# Computing mitochondrion length
python computemitolengths.py \
    data/pni_mito_analysisids_v185_extrastats.csv \
    data/pni_analysis_ids.csv \
    data/mitotoskel_v185.h5 \
    data/pni_mito_analysisids_v185_fullstats.csv

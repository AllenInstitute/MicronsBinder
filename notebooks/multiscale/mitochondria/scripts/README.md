# Computing intermediate data

These scripts are designed to be run from the `mitochondria` directory above this one.

### Computing basic branch statistics
Here, a branch refers to a connected component of a skeleton that possesses
the same compartment label. This means that branches can also refer to the
soma.

This script computes a data frame with an ID, volume, and surface area
measurement for each branch. The branch ID is the minimum skeleton node index
that is contained within the branch.

NOTE: this requires downloading the meshes (through downloadNeuronMeshes.sh)
```
python scripts/computebasicbranchstats.py \
    data/pycids.csv \
    data/basicbranchstats.csv
```


### Associating mitochondria to skeleton nodes
More precisely, finding the set of skeleton nodes that is the closest to
some mesh vertex of each mitochondrion. Labels each such mitochondrion
with the set of skeleton node indices that is paired with it in this way.

NOTE: this requires downloading the mitochondria meshes (through downloadmitomeshes.py)
```
python scripts/mitotoskel.py \
    data/pni_mito_cellswskel_v185.csv \
    data/temp/assoc_
python scripts/assemblemitotoskel.py \
    data/pni_mito_cellswskel_v185.csv \
    data/pycids.csv data/mitotoskel.h5
```


### Computing MCI and assigning mitochondria to compartments
This adds a surface area measurement, mitochondrial complexity index,
and a compartment assignment to each row of the mitochondria stats
data frame.

NOTE: this requires downloading the mitochondria meshes (through downloadmitomeshes.py)
```
python scripts/computeextramitostats.py \
    data/pni_mito_cellswskel_v185.csv \
    data/mitotoskel.h5 \
    data/fullmitostats.csv
```


### Computing density statistics
For each compartment branch, we count the number of synapses associated
with each branch and compute its total mitochondrial volume. We then
compute synapse and mitochondrial density metrics, and write out
two dataframes as a result - one describing a branch per row, and
another describing a cell per row, where each branch of a compartment
is summed within each cell.

The mitochondrial segmentation was performed in a larger volume than the
cellular segmentation. This means that several mitochondria poke out from
the tips of branches at the ends of the volume. To quantify mitochondrial
density accurately at these points, we only count the volume that overlaps
with a given segmentation ID.
```
python scripts/computedensitystats.py \
    data/fullmitostats.csv \
    data/pycids.csv data/pycinputs.csv \
    data/basicbranchstats.csv \
    data/mitochondriaoverlap.csv \
    data/pni_nucleus_segments.csv \
    data/pni_nucleus_lookup.csv \
    data/fullbranchstats.csv \
    data/cellwisestats.csv
```

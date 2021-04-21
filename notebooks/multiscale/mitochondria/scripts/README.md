# Computing intermediate data

These scripts are written to be run from the `mitochondria` directory above this one.

### Computing basic branch statistics
NOTE: this requires downloading the meshes through downloadNeuronMeshes.sh
```
python scripts/computebasicbranchstats.py \
    data/pycids.csv \
    data/basicbranchstats.csv
```


# Associating mitochondria to skeleton nodes
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
NOTE: this requires downloading the mitochondria meshes (through downloadmitomeshes.py)
```
python scripts/computeextramitostats.py \
    data/pni_mito_cellswskel_v185.csv \
    data/mitotoskel.h5 \
    data/fullmitostats.csv
```


### Computing density statistics
```
python scripts/computedensitystats.py \
    data/fullmitostats.csv \
    data/pycids.csv data/pycinputs.csv \
    data/basicbranchstats.csv \
    data/mitochondriaoverlap.csv \
    data/fullbranchstats.csv \
    data/cellwisestats.csv
```

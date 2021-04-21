# Mitochondria analysis

* `visualize.ipynb`
  A brief notebook that demonstrates one way to view a neuron rendered with its overlapping mitochondria.
* `compartment_comparison.ipynb` 
  An analysis notebook comparing mitochondria across neuronal compartments (Figure 4D)
* `density_analysis.ipynb`
  An analysis notebook relating mitochondrial volume density against synapse surface density


### Base data

* `pni_mito_cellswskel_v185.csv`
  A dataframe that describes each mitochondrion that overlaps with one of the cells with a soma in the volume and a skeleton.
* `pni_nucleus_segments.csv`
  A dataframe that describes each predicted nucleus segment in the volume.


### Intermediate data

* `basicbranchstats.csv`
  Contains branch IDs, volume measurements and surface area measurement for each compartment branch of the analyzed cells. Created by `scripts/computebasicbranchstats.py`
* `mitotoskel.h5`
  Contains associations between mitochondria and the skeleton nodes of each cell. Created by `scripts/mitotoskel.py`
* `fullmitostats.csv`
  A dataframe that adds a surface area measurement, the mitochondrial complexity index, and an inferred compartment label for each mitochondrion in `pni_mito_cellswskel_v185.csv`. Created by `scripts/computeextramitostats.py`
* `fullbranchstats.csv`
  A dataframe that adds synapse count, synapse surface density, mitochondrial volume, and mitochondrial volume density to the branch statistics in `basicbranchstats.csv`. Computed by `scripts/computedensitystats.py`
* `cellwisestats.csv`
  A dataframe that agglomerates the statistics from `fullbranchstats.csv` for each cell's compartments in bulk. Created by `scripts/computedensitystats.py`

### Scripts

See `scripts/README.md` and `scripts/derivedata.sh` for a walkthrough of how to compute the intermediate data.

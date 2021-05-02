[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AllenInstitute/MicronsBinder/master?urlpath=lab) [![](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/download/releases/3.7.0/)

# MicronsBinder
A collection of notebooks to provide examples of using data from [microns-explorer.org](https://microns-explorer.org).  The repository is designed to work with [mybinder.org](https://mybinder.org)


# Contents

### Introductory Notebooks
We've created some [introductory notebooks](notebooks/intro) to demonstrate some potential uses of the data. See:

* [MostSynapsesInAndOut](notebooks/intro/MostSynapsesInAndOut.ipynb)  
This notebook introduces you to reading synapses and the soma valence table.  It creates neuroglancer links that let you explore the inputs and outputs of individual neurons.
* [DashSynapseExplorer](notebooks/intro/DashSynapseExplorer.ipynb)  
This notebook shows you how to create dynamic scatterplots that recreate some of the results about layer 2/3 to layer 2/3 connections that were reported in Dorkenwald et al. 2019.
* [ImageAndSegmentationDownload](notebooks/ImageAndSegmentationDownload.ipynb)  
This notebook shows you how to create figures with overlaid EM and segmentation figures.

The introductory notebooks below are not intended to be run on mybinder for reasons specified below. To run them you should set up a local python environment (see [these instructions](#localenv)).

* [MeshExample](notebooks/intro/MeshExample.ipynb)
This demonstrates basic 3D visualization of meshes and skeletons using vtk, as well as calculating shortest paths along a mesh. This example uses more memory than allocated to most binder instances, and may be killed during processing while using those resources.
* [Render3DScaleBar](notebooks/intro/Render3DScaleBar.ipynb)
This demonstrates two techinques to create 3D scale bars on 3D visualization plots. It requires access to an X windows system to view these plots.

### Multiscale manuscript analyses  
These notebooks walk through some newer analyses studying the Phase I electron microscopy volume.  

These include:
* [Motif analysis](notebooks/multiscale/motifs) of a proofread synaptic connectivity graph between pyramidal cells.
* [Functional analysis](notebooks/multiscale/function) of a subset of the same pyramidal cells based on two-photon calcium imaging. The analysis relates local connectivity to function.
* [Mitochondria analysis](notebooks/multiscale/mitochondria) comparing mitochondria across neuronal compartments (axon, dendrite, soma), as well as a relation between mitochondrial density and synapse density.
* [Resource statistics](notebooks/multiscale/resource_statistics) that summarize neurite branch length and the expected completeness of the dendritic arbors in the volume.

See each [directory](notebooks/multiscale) and our [biorxiv manuscript](https://www.biorxiv.org/content/10.1101/2020.10.14.338681v3) for more details.


# <a name="localenv"></a> Local environment
A local environment for running the intermediate code generation scripts can be installed using the Anaconda environment installed within the binder and the `postBuild` script
```
conda env create -f environment.yml
bash postBuild
```
This installs the required python packages for running the basic code and the jupyter extensions for any plots and visualizations.

If you'd like to use these notebooks as part of your general jupyter environment. You'll likely need to install this environment into your ipython kernels.
```
conda activate micronsbinder
python -m ipykernel install --user --name=micronsbinder
```
You can then select this python environment when opening the relevant notebooks.


# Related projects
The notebooks contained here make heavy use of standard python tools, but also tools built as part of the collaboration between the Connectomics group at Allen Institute, the Seung Lab at Princeton, and the Tolias lab at Baylor, along with neuroglancer (developed
by Jeremy Maitin-Shepard from the Connectomics group at Google).

* [neuroglancer](https://www.github.com/google/neuroglancer)  
This is the main neuroglancer repository developed by Jeremy Maitin-Shepard. 
* [neuroglancer Seung-lab](https://www.github.com/seung-lab/neuroglancer)  
This is the Seung lab's fork of neuroglancer that has some alternative features added by many different Seung lab members. 
* [NeuroglancerAnnotationUI (nglui)](https://www.github.com/seung-lab/NeuroglancerAnnotationUI)  
This is a package principally developed by Casey Schneider-Mizell from the Allen Institute.  The package helps to create a pipeline that connects [pandas](https://pandas.pydata.org/) dataframes to neuroglancer links that visualize the contained data.
* [CloudVolume](https://www.github.com/seung-lab/cloud-volume)  
This is a python library developed principally by Will Silversmith from the Seung Lab for reading and writing volumetric data (e.g. EM images, segmentation), meshes, and skeletons to a variety of storage locations (e.g. cloud buckets, chunked files).
* [MeshParty](https://www.github.com/sdorkenw/MeshParty)  
This is a package developed by Sven Dorkenwald (Princeton), Forrest Collman (Allen), and Casey Schneider-Mizell (Allen) to make downloading meshes (via [CloudVolume](https://www.github.com/seung-lab/cloud-volume)), performing analysis (with tools like [trimesh](https://github.com/mikedh/trimesh), and [scipy](https://www.scipy.org/)) and visualization (via [vtk](https://pypi.org/project/vtk/)) of neuronal meshes easier.  There are also some tools for helping make dynamic movies of these data.
* [DashDataFrame](https://www.github.com/AllenInstitute/DashDataFrame)  
This is a package developed by Leila Elabaddy, Melissa Hendershott, and Forrest Collman at the Allen Institute.  It simplifies constructing dynamic visualization from pandas dataframes using [Dash](https://www.github.com/plotly/dash), including making dynamic links out to external services.  In this case, we use this to make dynamic scatterplots that allow you to select variables to plot, select and filter data points, and construct neuroglancer views of the specific locations in the dataset of those data points.


# Level of Support
We are releasing this repository as-is, and plan to update it without a fixed schedule.
It is intended to be a teaching tool to start interacting with the data. Community involvement is encouraged through both issues and pull requests.


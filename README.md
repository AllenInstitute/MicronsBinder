[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AllenInstitute/MicronsBinder/master?urlpath=lab)

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

These notebooks are not intended to be run on mybinder, because they require use of VTK, which interfaces with your graphics card to make 3d graphics.  To run them you should setup a local python environment with the requirements outlined in [environment.yml](environment.yml)

* [MeshExample](notebooks/intro/MeshExample.ipynb)
This demonstrates some basic 3d visualization of a meshes and skeletons using vtk, as well as calculating shortest paths along a mesh.
* [Render3dScaleBar](notebooks/intro/Render3dScaleBar.ipynb)
This demonstrates two techinques to create 3d scale bars on 3d visualization plots.

### Multiscale manuscript analyses  
These notebooks walk through some newer analyses. See each [directory](notebooks/multiscale) to see the types of analysis.


# Local environment
A local environment for running the intermediate code generation scripts can be installed using the Anaconda environment installed within the binder and the `postBuild` script
```
conda env create -f environment.yml
bash postBuild
```
This installs the required python packages for running the basic code and the jupyter extensions for any plots and visualizations.


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


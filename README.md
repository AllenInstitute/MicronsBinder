[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/AllenInstitute/MicronsBinder/master?urlpath=lab)

# MicronsBinder
A collection of notebooks to provide examples of using data from [microns-explorer.org](https://microns-explorer.org).  The repository is designed to work with [mybinder.org](https://mybinder.org)

# Reccomended Notebook Order

### [MostSynapsesInAndOut](notebooks/MostSynapsesInAndOut,ipynb)
This notebook introduces you to reading synapses and the soma valence table.  It creates neuroglancer links that let you explore the inputs and outputs of individual neurons.

###  [DashSynapseExplorer](notebooks/DashSynapseExplorer.ipynb)
This notebook shows you how to create dynamic scatterplots that recreate some of the results about layer 2/3 to layer 2/3 connections that were reported in Dorkenwald et al. 2019.

### [ImageAndSegmentationDownload](notebooks/ImageAndSegmentationDownload.ipynb)
This notebook shows you how to create figures with overlaid EM and segmentation figures.

## Local Notebooks
These notebooks are not intended to be run on mybinder, because they require use of VTK, which interfaces with your graphics card to make 3d graphics.  To run them you should setup a local python environment with the requirements outlined in [environment.yml](environment.yml)

### [MeshExample](notebooks/MeshExample.ipynb)
This demonstrates some basic 3d visualization of a meshes and skeletons using vtk, as well as calculating shortest paths along a mesh.

### [Render3dScaleBar](notebooks/Render3dScaleBar.ipynb)
This demonstrates two techinques to create 3d scale bars on 3d visualization plots.

# Related projects
The notebooks contained here make heavy use of standard python tools, but also tools built as part of the collaboration between the Connectomics group at Allen Institute, the Seung Lab at Princeton, and the Tolias lab at Baylor, along with neuroglancer whose development
by Jeremy Maiten-Shepard from the Connectomics group at Google.

Here is a list of some of the relevant repositories to know about

## [neuroglancer](https://www.github.com/google/neuroglancer)
This is the main neuroglancer repository developed by Jeremy Maiten-Shepard. 

## [neuroglancer Seung-lab](https://www.github.com/seung-lab/neuroglancer)
This is the Seung lab's fork of neuroglancer that has some alternative features
added by many different Seung lab members. 

## [NeuroglancerAnnotationUI (nglui)](https://www.github.com/seung-lab/NeuroglancerAnnotationUI)
This is a package principled developed by Casey Schneider-Mizell from the Allen Institute,
which helps you create a pipeline of pandas dataframes to neuroglancer links.

## [CloudVolume](https://www.github.com/seung-lab/cloud-volume)
This is a python library developed principly by Will Silversmith from the Seung Lab for reading and writing image data, segmentation data, meshes and skeletons to a variety of storage locations (cloud buckets, files, etc).

## [MeshParty](https://www.github.com/sdorkenw/MeshParty)
This is a package developed by Sven Dorkenwald (Princeton), Forrest Collman (Allen),
and Casey Schneider-Mizell (Allen) to make downloading meshes (via cloud-volume),
performing analysis (with tools like trimesh, and scipy) and visualization (via vtk) of neuronal meshes easier.  There are also some tools for helping make dynamic movies of these data.

## [DashDataFrame](https://www.github.com/AllenInstitute/DashDataFrame)
This is a package developed by Leila Elabaddy, Melissa Hendershott, and Forrest Collman at the Allen Institute.  It simplifies constructing dynamic visualization from pandas dataframes using [Dash](https://www.github.com/plotly/dash), including making dynamic links out to external services.  In this case, we use this to make dynamic scatterplots that allow you to select variables to plot, select and filter data points, and construct neuroglancer views of the specific locations in the dataset of those data points.

# Level of Support
We are releasing this repository as is, and plan on update it without a fixed schedule.
It is intended to be a teaching tool to start interacting with the data. Community involvement is encouraged through both issues and pull requests.


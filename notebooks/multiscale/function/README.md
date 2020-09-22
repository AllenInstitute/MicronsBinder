# Functional data
2-photon calcium imaging data acquired by Baylor College of Medicine.
![](figures/function_info_fig.png)
### Raw data
#### Calcium video
- Calcium videos can be downloaded from [here](https://drive.google.com/drive/folders/1nL0_asZkqiWrgkE-tpXIswf84tEdBwq_?usp=sharing)
- Each video is tiff file with size 256 (x) x 256 (y) x 27300 (time).
#### Visual stimulus
- Stimulus videos and the labels can be downloaded from [here](https://drive.google.com/drive/folders/1-hLrXYclGwQmCX0VhjyrqJ8rpLsDSLgK?usp=sharing)
##### Stimulus video
- Each video is tiff file with size 90 (x) x 160 (y) x 27300 (time).
##### Stimulus labels
- Each file has length 27300 (time).
- The value indicates the angle of the directional stimulus at that time frame.
- If it's empty (NaN), it means that the noise stimulus is shown at that time frame.

## Skeletons for pyramidal cells
Skeletons for pyramidal cells can be downloaded from [here](https://drive.google.com/drive/folders/1_6jVwOx0YQE9all7cnf75xYmyWCYsiUk?usp=sharing).  
`skeletons`: Skeleton files for pyramidal cells named by its cell id.  
`labels`: Skeleton compartment labels (n_vertices x 1). Each value indicates which compartment each skeleton node belongs to.  
> 0: Soma  
> 1: Axon  
> 2: Basal dendrite  
> 3: Apical dendrite  
> 4: Ambiguous dendrite  
> 5: Ambiguous  

#### How to load skeletons
1. Install [MeshParty](https://meshparty.readthedocs.io/en/latest/includeme.html#meshparty)  
2. Load skeletons using [skeleton_io](https://meshparty.readthedocs.io/en/latest/guide/skeletons.html)  
```
from meshparty import skeleton_io

skeleton = skeleton_io.read_skeleton_h5("648518346349537297_skeleton.h5")
```
- `skeleton.vertices` : Coordinates of vertices [nm].  
- `skeleton.csgraph` : Skeleton graph where each element refers to distance from node i to j [nm].
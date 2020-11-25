# Structure-function analysis (Figure 6)

## DataJoint database
This public database contains extracted structural and functional data such as visual responses and connectivity data.

### Registration
You need to be registered to access the database.  
Please fill out the [registration form](https://forms.gle/6SeDGRT8zoLqpbfU9) to receive user id and password.

### Installation
DataJoint for Python requires Python 3.4 or above to function properly.
```
pip3 install datajoint
```
For more information, please checkout the [DataJoint Tutorials](https://tutorials.datajoint.io/setting-up/datajoint-python.html).  

### Database configuration
- HOST: datajoint.ninai.org
- USER: Given after registration
- PASSWORD: Given after registration

### Accessing the database
```python3
import datajoint as dj

# Datajoint credentials
dj.config["database.host"] = "datajoint.ninai.org"

dj.conn() # Then it will ask for your net id and password

pinky = dj.create_virtual_module("seung_pinky", "seung_pinky")
```
*Pinky is the nickname for this dataset named after the American animated television series, *Pinky and the Brain*.

### Fetching data from the database
Example extracting visual response trace from the database.
```python3
key = {"segment_id": 648518346349539895, "scan_id": 2}
trace = (pinky.EASETrace() & key).fetch1("trace")
```
Please refer to the [DataJoint Tutorials](https://tutorials.datajoint.io/setting-up/datajoint-python.html) for additional help.  

## Raw data

2-photon calcium imaging data acquired by Baylor College of Medicine.
![](figures/function_info_fig.png)

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

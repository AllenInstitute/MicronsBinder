# Structure-function analysis

## Outline
### scripts
- `compute_neurite_length.py`: Compute axon/dendrite length using skeletons.
- `compute_degree.py`: Compute out/in-degree from the connectivity graph.
```python3
python [script_name.py]
```
The results will be saved in the `data` folder.
### notebooks
- `structure_function_analysis.ipynb` (Figure 7): Analysis testing whether visual response (response strength, intermittency, OSi) is correlated with local connectivity. 
- `orientation_direction_tuning.ipynb`: Compute tuning curves and determine significantly tuned cells with permutation test.
- `spatial_organization.ipynb` (Figure S7): Spatial organization of in-connection density in the volume.

## DataJoint database
This public database which contains extracted structural and functional data such as visual responses and connectivity data.  
  
![](figures/pinky_schema.png)  
*Refer to [SeungLab DataJoint repository](https://github.com/seung-lab/datajoint_seung) for detailed information of the tables in the schema.

### Registration
You need to be registered to access the database.  
Please fill out the [registration form](https://forms.gle/6SeDGRT8zoLqpbfU9) to receive user id and password.
*Currently, it is temporarily available without registration.

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

dj.conn() # Then it will ask for your id and password

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

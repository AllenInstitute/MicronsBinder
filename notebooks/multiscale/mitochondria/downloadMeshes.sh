mkdir -p data/neuronmeshes
https://drive.google.com/file/d//view?usp=sharing
gdown --id 1lPWicWMBsIiXur68Sqf_FMS5ZzHoSqLj -O data/neuronmeshes/neuron_meshes_v185.tgz \
&& pushd data/neuronmeshes \
&& tar -xzf neuron_meshes_v185.tgz \
&& popd

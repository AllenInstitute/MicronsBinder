mkdir -p data/neuronmeshes
gdown --id 1lPWicWMBsIiXur68Sqf_FMS5ZzHoSqLj -O data/neuronmeshes/neuron_meshes_v185.tgz \
&& pushd data/neuronmeshes \
&& tar -xzf neuron_meshes_v185.tgz \
&& popd

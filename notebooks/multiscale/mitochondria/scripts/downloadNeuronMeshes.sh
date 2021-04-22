mkdir -p ../data/neuron_meshes_v185
pushd ../data/neuron_meshes_v185 \
&& gdown --id 1lPWicWMBsIiXur68Sqf_FMS5ZzHoSqLj -O neuron_meshes_v185.tgz \
&& tar -xzf neuron_meshes_v185.tgz \
&& popd

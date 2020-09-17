mkdir -p data

gdown --id 1wLwIUQtIxO2HvfbPm9YG0DC0gRrtEvMb -O data/multiscale_resource_intermediate_data.tgz \
&& pushd data \
&& tar -xzf multiscale_resource_intermediate_data.tgz \
&& popd

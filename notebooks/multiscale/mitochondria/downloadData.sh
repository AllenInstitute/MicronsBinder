mkdir -p data
gdown --id 1vg8eRuLm9OJsfGt5CotJ2ML6k5zgmz-f -O data/multiscale_mito_intermediate_data.tgz \
&& pushd data/ \
&& tar -xzf multiscale_mito_intermediate_data.tgz \
&& popd

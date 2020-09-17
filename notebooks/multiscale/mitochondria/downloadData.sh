mkdir -p data
gdown --id 1OIUpHr9gVYKVagDqpom9Sx9qb5iAX5Km -O data/multiscale_mito_intermediate_data.tgz \
&& pushd data/ \
&& tar -xzf multiscale_mito_intermediate_data.tgz \
&& popd

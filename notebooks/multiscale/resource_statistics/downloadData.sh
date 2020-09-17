mkdir -p data

gdown --id 1Roj4WFBofLsRQVY8FHTE_0fq9EJPKGRt -O data/multiscale_resource_intermediate_data.tgz \
&& pushd data \
&& tar -xzf multiscale_resource_intermediate_data.tgz \
&& popd

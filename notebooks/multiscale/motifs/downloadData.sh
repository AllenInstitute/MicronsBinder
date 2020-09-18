mkdir -p data
gdown --id 1SIRrG8KuobgjuByBHxRt11qkM_Fj7JSN -O data/multiscale_motif_intermediate_data.tgz \
&& pushd data/ \
&& tar -xzf multiscale_motif_intermediate_data.tgz \
&& popd

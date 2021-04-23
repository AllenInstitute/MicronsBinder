mkdir -p data/mitomeshes
gdown --id XXX -O data/mitomeshes/mito_meshes.tgz \
&& pushd data/mitomeshes \
&& tar -xzf mito_meshes.tgz \
&& popd

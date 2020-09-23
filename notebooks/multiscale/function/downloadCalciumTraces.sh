mkdir -p data
gdown --id 1S2qE6-oNP_KMvkobRxswkH0DqCJDbPnu -O data/calcium_trace.tgz \
&& pushd data/ \
&& tar -xzf calcium_trace.tgz \
&& popd

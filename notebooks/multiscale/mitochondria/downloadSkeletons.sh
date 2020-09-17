mkdir -p data
gdown --id 1Tpa-p7rwOZ65rIwDPbm4famUFtqhbqOK -O data/smoothed_skeletons_v185.tgz \
&& pushd data/ \
&& tar -xzf smoothed_skeletons_v185.tgz \
&& popd

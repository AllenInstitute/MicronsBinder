mkdir -p data/smoothed_skels

gdown --id 1Tpa-p7rwOZ65rIwDPbm4famUFtqhbqOK -O data/smoothed_skels/smoothed_skeletons_v185.tgz \
&& pushd data/smoothed_skels \
&& tar -xzf smoothed_skeletons_v185.tgz \
&& popd

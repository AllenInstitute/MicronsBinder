mkdir -p data
mkdir -p data/soma_ids
mkdir -p data/smoothed_skels
mkdir -p data/resource_vis_info

gdown --id 1_Xciy8oPcpegRMuM3AUNU7tpHgHtrLlv -O data/soma_ids/p100_pyr_soma_IDs_v185.csv
gdown --id 1k3XtA3LPcK3MN3BwUIRKcqJEU9PTi3Gz -O data/soma_ids/p100_inh_soma_IDs_v185.csv

gdown --id 1Tpa-p7rwOZ65rIwDPbm4famUFtqhbqOK -O data/smoothed_skels/smoothed_skeletons_v185.tgz \
&& pushd data/smoothed_skels \
&& tar -xzf smoothed_skeletons_v185.tgz \
&& popd

gdown --id 1thwNHGh8Boa6iLMK7o_aOnA59QeCz2j2 -O data/resource_vis_info/resource_vis_info.tgz \
&& pushd data/resource_vis_info \
&& tar -xzf resource_vis_info.tgz \
&& popd



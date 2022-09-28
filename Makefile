pwd= $( pwd )

all: jupyter_extensions intro \
	vignette_analysis/data \
	vignette_analysis/function vignette_analysis/mitochondria \
	vignette_analysis/motifs vignette_analysis/resource_statistics

clean: clean-intro_datadir \
	clean-vignette_analysis/data clean-vignette_analysis/mitochondria \
	clean-vignette_analysis/motifs clean-vignette_analysis/resource_statistics \
	clean-vignette_analysis/function


jupyter_extensions:
	printf "\033[1;33m%s\033[0;39m\n" \
		"Installing jupyter extensions"
	mkdir -p ${HOME}/.jupyter
	echo "c.NotebookApp.iopub_data_rate_limit=1e22" >> \
		${HOME}/.jupyter/jupyter_notebook_config.py
	jupyter nbextension install --py --sys-prefix itkwidgets
	jupyter nbextension enable --py --sys-prefix itkwidgets
	cd /tmp/
	jupyter labextension install @jupyter-widgets/jupyterlab-manager --no-build
	jupyter labextension install jupyter-matplotlib jupyterlab-datawidgets \
		jupyter-webrtc itkwidgets --no-build
	jupyter labextension install jupyterlab-dash --no-build
	jupyter lab build --minimize=False  # minimize uses too much memory for binder
	cd ${pwd}


intro_datadir= notebooks/intro/data
intro: jupyter_extensions
	printf "\033[1;33m%s\033[0;39m\n" \
		"Downloading intro notebook data"
	mkdir -p ${intro_datadir}
	wget -q -O ${intro_datadir}/pni_synapses_v185.csv \
		"https://zenodo.org/record/3710459/files/pni_synapses_v185.csv?download=1"
	wget -q -O ${intro_datadir}/soma_valence_v185.csv \
		"https://zenodo.org/record/3710459/files/soma_valence_v185.csv?download=1"
	wget -q -O ${intro_datadir}/soma_subgraph_synapses_spines_v185.csv \
		"https://zenodo.org/record/3710459/files/soma_subgraph_synapses_spines_v185.csv?download=1"

clean-intro_datadir: 
	if [ -d ${intro_datadir} ]; then rm -r ${intro_datadir}; fi


mm3_intro: jupyter_extensions


vignette_datadir= notebooks/vignette_analysis/data
vignette_tarfile= ${vignette_datadir}/211019_smoothed_skeletons.tgz 
vignette_url= https://zenodo.org/record/6363348/files/211019_smoothed_skeletons.tgz
vignette_analysis/data:
	mkdir -p ${vignette_datadir}
	if [ ! -f ${vignette_tarfile} ]; \
		then printf "\033[1;33m%s\033[0;39m\n" \
			"Downloading base vignette analysis data"; \
		wget -q -O ${vignette_tarfile} ${vignette_url} \
		&& tar -C ${vignette_datadir} -xzf ${vignette_tarfile}; \
	fi

clean-vignette_analysis/data:
	if [ -d ${vignette_datadir} ]; then rm -r ${vignette_datadir}; fi


func_datadir= notebooks/vignette_analysis/function/data
func_tarfile= ${func_datadir}/211019_vignette_functional_analysis_data.tgz
func_datatablefile= ${func_datadir}/function_data_tables.tgz
func_url= https://zenodo.org/record/6363348/files/211019_vignette_functional_analysis_data.tgz
vignette_analysis/function: vignette_analysis/data
	mkdir -p ${func_datadir}
	if [ ! -f ${func_tarfile} ]; \
		then printf "\033[1;33m%s\033[0;39m\n" \
			"Downloading functional vignette analysis data"; \
		wget -q -O ${func_tarfile} ${func_url} \
		&& tar -C ${func_datadir} -xzf ${func_tarfile} \
		&& tar -C ${func_datadir} -xzf ${func_datatablefile}; \
	fi

clean-vignette_analysis/function:
	if [ -d ${func_datadir} ]; then rm -r ${func_datadir}; fi


mito_datadir= notebooks/vignette_analysis/mitochondria/data
mito_tarfile= ${mito_datadir}/vignette_mito_analysis_data.tgz
vignette_analysis/mitochondria: vignette_analysis/data
	mkdir -p ${mito_datadir}
	if [ ! -f ${mito_tarfile} ]; \
		then printf "\033[1;33m%s\033[0;39m\n" \
			"Downloading mitochondria vignette analysis data"; \
		gdown --id 1MruDsi7q77G5Mq0b1t95TDFwVoHDK3WB -O ${mito_tarfile} \
		&& tar -C ${mito_datadir} -xzf ${mito_tarfile}; \
	fi

clean-vignette_analysis/mitochondria:
	if [ -d ${mito_datadir} ]; then rm -r ${mito_datadir}; fi


motif_dir= notebooks/vignette_analysis/motifs
motif_tarfile= ${motif_dir}/211019_vignette_motif_analysis_data.tgz
motif_url= https://zenodo.org/record/6363348/files/211019_vignette_motif_analysis_data.tgz
vignette_analysis/motifs: vignette_analysis/data
	if [ ! -f ${motif_tarfile} ]; \
		then printf "\033[1;33m%s\033[0;39m\n" \
			"Downloading motif vignette analysis data"; \
		wget -O ${motif_tarfile} -q ${motif_url} \
		&& tar -C ${motif_dir} -xzf ${motif_tarfile} \
		&& tar -C ${motif_dir} -xzf ${motif_dir}/vignette_motif_base_data.tgz \
		&& tar -C ${motif_dir} -xzf ${motif_dir}/vignette_motif_intermediate_data.tgz; \
	fi

clean-vignette_analysis/motifs:
	if [ -d ${motif_dir}/data ]; then rm -r ${motif_dir}/data; fi
	if [ -d ${motif_dir}/saved ]; then rm -r ${motif_dir}/saved; fi
	if [ -f ${motif_tarfile} ]; then rm ${motif_tarfile}; fi
	if [ -f ${motif_dir}/vignette_motif_base_data.tgz ]; \
		then rm ${motif_dir}/vignette_motif_base_data.tgz; fi
	if [ -f ${motif_dir}/vignette_motif_intermediate_data.tgz ]; \
		then rm ${motif_dir}/vignette_motif_intermediate_data.tgz; fi


resource_datadir= notebooks/vignette_analysis/resource_statistics/data
resource_tarfile= ${resource_datadir}/vignette_resourcestats_analysis_data.tgz
vignette_analysis/resource_statistics: vignette_analysis/data
	mkdir -p ${resource_datadir}
	if [ ! -f ${resource_tarfile} ]; \
		then printf "\033[1;33m%s\033[0;39m\n" \
			"Downloading vignette resource statistics data"; \
		gdown --id 1JzpHNK0wzj7JfmdM4fECR4yPZT6BE9ZB \
			-O ${resource_tarfile} \
		&& tar -C ${resource_datadir} -xzf ${resource_tarfile}; \
	fi

clean-vignette_analysis/resource_statistics:
	if [ -d ${resource_datadir} ]; then rm -r ${resource_datadir}; fi

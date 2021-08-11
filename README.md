# Installation

MLPack dependencies are listed here: https://github.com/mlpack/mlpack#3-dependencies

To install all dependencies (hopefully):

	sudo apt -y install gcc g++ perl pkgconf make cmake liblapack-dev libblas-dev libboost-all-dev libarmadillo-dev libensmallen-dev

To build everything, do the following:

	git clone https://github.com/msanderford/mlpack_sg_lasso/
	cd mlpack_sg_lasso
	bash build_script.sh

# Components

## preprocess:

	parameter 1: response matrix file
	parameter 2: file containing list of alignment file paths
	parameter 3: basename of output files
	optional parameters (must be specified after parameters 1-3):
		n: "normalize" feature weights by column
		ct {N}: ignore mutations observed fewer than {N} times (must be an integer)
sample usage:

	cd sample_files
	../mlpack-3.2.2/build/bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns.txt angiosperm_input
	mv angiosperm_input ..
	cd ..

Notes and Caveats

The paths for the alignment files must be relative to the directory that preprocess is run from.

Any species with a response value of 0 in the predictions file (first argument to preprocess executable) will be left out of the features file entirely.

The preprocessing script is not very well tested. The easiest way to check if it worked correctly is to inspect the feature mapping file - if something went wrong, it will usually have nonsense characters as the alleles for some feature columns.

The list of alignments given to the preprocess program can also generate overlapping groups of input for the overlapping SGLasso algorithm, by specifying multiple comma-separated files/genes on a single line.

For example, for a non-overlapping set of input groups, the file might look like:

	aln_dir/gene1.fas
	aln_dir/gene2.fas
	aln_dir/gene3.fas
	aln_dir/gene4.fas

While input where gene2 shares properties with both gene1 and gene3 might look like:

	aln_dir/gene1.fas,aln_dir/gene2.fas
	aln_dir/gene3.fas,aln_dir/gene2.fas
	aln_dir/gene4.fas


## mlpack_sg_lasso_leastr

	required inputs:
	  --features_file (-f)          Matrix containing feature set A.
	  --opts_ind_file (-n)          Matrix of indices defining non-overlapping group information.
	  --responses_file (-r)         Vector containing responses y.
	optional inputs:
	  --feature_weights_file (-w)   Output file to write learned feature weights to.
	  --lambda1 (-z)                Feature regularization parameter (z1 >=0). Default value 0.
	  --lambda2 (-y)                Group regularization parameter (z2 >=0). Default value 0.

sample usage:

	mlpack-3.2.2/build/bin/mlpack_sg_lasso_leastr -v -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights.xml


## mlpack_overlapping_sg_lasso_leastr

	required inputs:
	  --features_file (-f)          Matrix containing feature set A.
	  --opts_ind_file (-n)          Matrix of indices defining overlapping group information.
	  --field_file (-g)             Vector of feature indices for overlapping groups.
	  --responses_file (-r)         Vector containing responses y.
	optional inputs:
	  --feature_weights_file (-w)   Output file to write learned feature weights to.
	  --lambda1 (-z)                Feature regularization parameter (z1 >=0). Default value 0.
	  --lambda2 (-y)                Group regularization parameter (z2 >=0). Default value 0.

sample usage:

	mlpack-3.2.2/build/bin/mlpack_overlapping_sg_lasso_leastr -v -f angiosperm_input/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n angiosperm_input/group_indices_angiosperm_input.txt -g angiosperm_input/field_angiosperm_input.txt -r angiosperm_input/response_angiosperm_input.txt -w angiosperm_out_feature_weights.xml


# Parsing outputs

The weights in the output XML file are in the same order as the lines in the feature mapping file.
Simple commands to remove the XML formatting, merge the two, and remove all features with a weight of zero:

	grep -P "<item>.*</item>" angiosperm_out_feature_weights.xml | sed -re "s/.*<item>(.*)<\/item>.*/\1/" > temp_angiosperm_out_feature_weights.txt
	paste <(sed -e "1d" angiosperm_input/feature_mapping_angiosperm_input.txt) temp_angiosperm_out_feature_weights.txt | grep -v "0.00000000000000000e+00" > angiosperm_out_feature_weights.txt

# Phylogeny Testing Pipeline

The pipeline is implemented as a python script which takes a set of gene alignments and a newick tree with at least one named internal node then, for each named internal node, creates a model predicting the chance that a given sequence descends from that node, and finally applies each predictive model to each input sequence and generates a table of predictive values grouped by gene for each sequence.

In order to make a local installation of python 3.9.5 on the cluster and install the required python packages, use the following commands:

	wget https://www.python.org/ftp/python/3.9.5/Python-3.9.5.tgz
	tar xvf Python-3.9.5.tgz
	cd Python-3.9.5
	./configure --enable-optimizations --with-ensurepip=install --prefix=$HOME
	make -j 8
	make altinstall
	python3.9 -m pip install biopython matplotlib

sample usage:

	python3.9 mmlp_pipeline.py sample_files/mmlp_test.nwk sample_files/angiosperm_100_sample_alns.txt -o sample_output

Additionally, the phylogeny testing pipeline can be run in ensemble mode, which will randomly partition the alignment files into N partitions (--ensemble_parts) and build a model with each partition separately, then repeat this process M times (--ensemble_coverage), aggregating the results to analyze the importance of each gene.

sample usage:

	python3.9 mmlp_pipeline.py sample_files/mmlp_test.nwk sample_files/angiosperm_100_sample_alns.txt -o sample_output --ensemble_parts 5 --ensemble_coverage 3

Finally, when running in ensemble mode, the --sparsify argument may be added, which will cause the pipeline to run the analysis with an increasing group sparsity parameter until the number of genes assigned non-zero weights will fit into a single partition (i.e. if --ensemble_parts is set to 5, the pipeline will remove 80% of the genes). Once the set of selected genes reaches the threshold, it will perform ridge regression to assign a final set of weights to features in those genes.

	python3.9 mmlp_pipeline.py sample_files/mmlp_test.nwk sample_files/angiosperm_100_sample_alns.txt -o sample_output --ensemble_parts 5 --ensemble_coverage 3 --sparsify


#SGLASSO Settings

Additional parameters for the SGLASSO program can be set by creating an options file like the included slep_opts_example.txt file.

To provide your options file to the pipeline, use the argument "--slep_opts slep_opts_example.txt".

If using the mlpack_sg_lasso executables directly, your options file can be specified by adding the argument "-s slep_opts_example.txt".

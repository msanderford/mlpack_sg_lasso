Dependencies are listed here: https://github.com/mlpack/mlpack#3-dependencies

Command to install all dependencies (hopefully):

	sudo apt -y install gcc g++ perl pkgconf make cmake liblapack-dev libblas-dev libboost-all-dev libarmadillo-dev libensmallen-dev

To build everything, just do:

	bash build_script.sh


# preprocess program:

	parameter 1: response matrix file
	parameter 2: file containing list of alignment file paths
	parameter 3: basename of output files
	optional parameters (must be specified after parameters 1-3):
		n: "normalize" feature weights by column
		is: "ignore singletons" gives a weight of zero to any feature that occurs in only one species

sample usage:

	cd sample_files
	../mlpack-3.2.2/build/bin/preprocess angiosperm_20spec_pred.txt angiosperm_100_sample_alns.txt angiosperm_input [n] [is]

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


# mlpack_sg_lasso_leastr

	required inputs:
	  --features_file (-f)          Matrix containing feature set A.
	  --opts_ind_file (-n)          Matrix of indices defining non-overlapping group information.
	  --responses_file (-r)         Vector containing responses y.
	optional inputs:
	  --feature_weights_file (-w)   Output file to write learned feature weights to.
	  --lambda1 (-z)                Feature regularization parameter (z1 >=0). Default value 0.
	  --lambda2 (-y)                Group regularization parameter (z2 >=0). Default value 0.

sample usage:

	mlpack-3.2.2/build/bin/mlpack_sg_lasso_leastr -v -f sample_files/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n sample_files/group_indices_angiosperm_input.txt -r sample_files/response_angiosperm_input.txt -w angiosperm_out_feature_weights.xml


# mlpack_overlapping_sg_lasso_leastr

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

	mlpack-3.2.2/build/bin/mlpack_overlapping_sg_lasso_leastr -v -f sample_files/feature_angiosperm_input.txt -z 0.1 -y 0.5 -n sample_files/group_indices_angiosperm_input.txt -g sample_files/field_angiosperm_input.txt -r sample_files/response_angiosperm_input.txt -w angiosperm_out_feature_weights.xml


# Parsing output

The weights in the output XML file are in the same order as the lines in the feature mapping file.
Simple commands to remove the XML formatting, merge the two, and remove all features with a weight of zero:

	grep -P "<item>.*</item>" angiosperm_out_feature_weights.xml | sed -re "s/.*<item>(.*)<\/item>.*/\1/" > temp_angiosperm_out_feature_weights.txt
	paste <(sed -e "1d" sample_files/feature_mapping_angiosperm_input.txt) temp_angiosperm_out_feature_weights.txt | grep -v "0.00000000000000000e+00" > angiosperm_out_feature_weights.txt


Dependencies are listed here: https://github.com/mlpack/mlpack#3-dependencies

To build everything, just do:
	bash build_script.sh



Notes and Caveats

Any species with a response value of 0 in the predictions file (first argument to preprocess executable) will be left out of the features file entirely.

Every alignment file must have a sequence for every non-zero species in the predictions file. Be especially careful with this one, because the preprocess executable won't give an error or a warning, it will just silently give you messed up output. The easiest way to check this is to inspect the feature mapping file - if a species was missing from one of the alignments, it will have nonsense characters as the alleles for some feature columns.

The weights in the output XML file are in the same order as the lines in the feature mapping file.
Simple commands to remove the XML formatting, merge the two, and remove all features with a weight of zero:
	grep -P "<item>.*</item>" test_file_out_feature_weights.xml | sed -re "s/.*<item>(.*)<\/item>.*/\1/" > temp_test_file_out_feature_weights.txt
	paste feature_mapping_test_file.txt temp_test_file_out_feature_weights.txt | grep -v "0.00000000000000000e+00" > test_file_out_feature_weights.txt

Sample usage commands are in sample_run_commands.txt

-z parameter is feature-level sparsity coefficient (larger = more sparse)
-y parameter is group-level sparsity coefficient (larger = more sparse)


preprocess program:
	parameter 1: response matrix file
	parameter 2: file containing list of alignment file paths

mlpack_sg_lasso program:
	parameter 1:

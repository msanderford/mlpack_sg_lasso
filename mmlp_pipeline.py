import argparse
import shutil
import pipeline_funcs as pf


def main(args):
	hypothesis_file_list = pf.generate_hypothesis_set(args.tree, args.nodelist)
	features_filename, groups_filename, response_filename_list = pf.generate_input_matrices(args.aln_list, hypothesis_file_list, args.output)
	weights_file_list = pf.run_mlp(features_filename, groups_filename, response_filename_list)
	pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename, features_filename)
	for hypothesis_filename in hypothesis_file_list:
		shutil.move(hypothesis_filename, args.output)
		shutil.move(hypothesis_filename.replace(".txt","_out_feature_weights.xml"), args.output)
		shutil.move(hypothesis_filename.replace("hypothesis.txt","gene_predictions.txt"), args.output)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Phylogenetic hypothesis tester.")
	parser.add_argument("tree", help="Input phylogeny to perform testing on.", type=str)
	parser.add_argument("aln_list", help="List of alignment files to extract features from.", type=str)
	parser.add_argument("--nodelist", help="File containing list of named internal nodes to test. If no file is specified, each named internal node of input phylogeny will be tested.", type=str, default=None)
	parser.add_argument("-z", "--lambda1", help="Feature sparsity parameter.", type=float, default=0.1)
	parser.add_argument("-y", "--lambda2", help="Group sparsity parameter.", type=float, default=0.1)
	parser.add_argument("-o", "--output", help="Output directory.", type=str, default="output")
	args = parser.parse_args()
	main(args)
	# generate_gene_prediction_table(weights_filename, responses_filename, groups_filename, features_filename, output_filename)
	# if False:
	# 	pf.generate_gene_prediction_table("angiosperm_out_feature_weights.xml", "angiosperm_20spec_pred.txt", "angiosperm_input/group_indices_angiosperm_input.txt", "angiosperm_input/feature_angiosperm_input.txt", "test_output.txt")
	# if False:
	# 	hypothesis_file_list = pf.generate_hypothesis_set("testtree2.nwk", "nodelist.txt")
	# 	features_filename, groups_filename, response_filename_list = pf.generate_input_matrices("sample_files/angiosperm_100_sample_alns.txt", hypothesis_file_list)
	# 	weights_file_list = pf.run_mlp(features_filename, groups_filename, response_filename_list)
	# 	pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename, features_filename)




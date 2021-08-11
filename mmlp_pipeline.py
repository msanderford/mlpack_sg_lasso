import os
import argparse
import shutil
import math
import copy
import pipeline_funcs as pf
import gene_contribution_visualizer as gcv


def main(args):
	hypothesis_file_list = pf.generate_hypothesis_set(args.tree, args.nodelist, args.response)
	HSS = {}
	missing_seqs = set()
	merged_rep_predictions_files = {}
	if args.ensemble_parts is not None and args.ensemble_parts >= 1:
		tempdir_list = []
		merged_parts_prediction_files = {}
		for i in range(0, args.ensemble_coverage):
			partitioned_aln_lists = pf.split_gene_list(args.aln_list, args.ensemble_parts)
			j = 0
			gene_prediction_files = {}
			for part_aln_list in partitioned_aln_lists:
				tempdir = "{}_rep{}_part{}".format(args.output, i+1, j+1)
				tempdir_list.append(tempdir)
				features_filename_list, groups_filename_list, response_filename_list, gene_list = pf.generate_input_matrices(part_aln_list, hypothesis_file_list, args)
				with open(os.path.join(args.output, "missing_seqs_" + args.output + ".txt"), "r") as file:
					for line in file:
						data = line.strip().split("\t")
						missing_seqs.add((data[1], os.path.splitext(os.path.basename(data[0]))[0]))
				weights_file_list = pf.run_mlp(features_filename_list, groups_filename_list, response_filename_list, args.lambda1, args.lambda2, args.method, args.slep_opts)
				pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, HSS, missing_seqs)
				for hypothesis_filename in hypothesis_file_list:
					if i+1 == args.ensemble_coverage and j+1 == args.ensemble_parts and not args.sparsify:
						shutil.move(hypothesis_filename, args.output)
					else:
						shutil.copy(hypothesis_filename, args.output)
					shutil.move(hypothesis_filename.replace(".txt","_out_feature_weights.xml"), args.output)
					shutil.move(hypothesis_filename.replace("hypothesis.txt","gene_predictions.txt"), args.output)
					gene_prediction_files[hypothesis_filename].append(os.path.join(tempdir, hypothesis_filename.replace("hypothesis.txt","gene_predictions.txt")))
					shutil.move(hypothesis_filename.replace("hypothesis.txt", "mapped_feature_weights.txt"), args.output)
					shutil.move(hypothesis_filename.replace("hypothesis.txt", "PSS.txt"), args.output)
					shutil.move(hypothesis_filename.replace("hypothesis.txt", "GSS.txt"), args.output)
				shutil.move(args.output, tempdir)
				j += 1
			for hypothesis_filename in hypothesis_file_list:
				merged_parts_prediction_files[hypothesis_filename].append(gcv.merge_predictions(gene_prediction_files[hypothesis_filename],hypothesis_filename.replace("hypothesis.txt","merged_gene_predictions_rep{}.txt".format(i))))
		for hypothesis_filename in hypothesis_file_list:
			merged_rep_predictions_files[hypothesis_filename] = gcv.merge_predictions(merged_parts_prediction_files,hypothesis_filename.replace("hypothesis.txt","merged_gene_predictions_final.txt"))
			gcv.main(merged_rep_predictions_files[hypothesis_filename])
		os.mkdir(args.output)
		for tempdir in tempdir_list:
			shutil.move(tempdir, args.output)
		with open(os.path.join(args.output, "HSS.txt"), 'w') as file:
			file.write("{}\t{}\n".format("Hypothesis", "HSS"))
			for hypothesis_filename in hypothesis_file_list:
				file.write("{}\t{}\n".format(hypothesis_filename.replace("_hypothesis.txt", ""), HSS[hypothesis_filename]))
		result_files_list = pf.find_result_files(args, hypothesis_file_list)
		weights = pf.parse_result_files(args, result_files_list)
		return pf.analyze_ensemble_weights(args, weights)
	else:
		features_filename_list, groups_filename_list, response_filename_list, gene_list = pf.generate_input_matrices(args.aln_list, hypothesis_file_list, args)
		with open(os.path.join(args.output, "missing_seqs_" + args.output + ".txt"), "r") as file:
			for line in file:
				data = line.strip().split("\t")
				missing_seqs.add((data[1], os.path.splitext(os.path.basename(data[0]))[0]))
		weights_file_list = pf.run_mlp(features_filename_list, groups_filename_list, response_filename_list, args.lambda1, args.lambda2, args.method, args.slep_opts)
		pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, HSS, missing_seqs)
		for hypothesis_filename in hypothesis_file_list:
			shutil.move(hypothesis_filename, args.output)
			shutil.move(hypothesis_filename.replace(".txt","_out_feature_weights.xml"), args.output)
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "gene_predictions.txt"), args.output)
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "mapped_feature_weights.txt"), args.output)
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "PSS.txt"), args.output)
			shutil.move(hypothesis_filename.replace("hypothesis.txt", "GSS.txt"), args.output)
		with open(os.path.join(args.output, "HSS.txt"), 'w') as file:
			file.write("{}\t{}\n".format("Hypothesis", "HSS"))
			for hypothesis_filename in hypothesis_file_list:
				file.write("{}\t{}\n".format(hypothesis_filename.replace("_hypothesis.txt", ""), HSS[hypothesis_filename]))
		return None


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Phylogenetic hypothesis tester.")
	parser.add_argument("tree", help="Input phylogeny to perform testing on.", type=str)
	parser.add_argument("aln_list", help="List of alignment files to extract features from.", type=str)
	parser.add_argument("--nodelist", help="File containing list of named internal nodes to test. If no file is specified, each named internal node of input phylogeny will be tested.", type=str, default=None)
	parser.add_argument("--response", help="File containing list of named node/response value pairs.", type=str, default=None)
	parser.add_argument("-z", "--lambda1", help="Feature sparsity parameter.", type=float, default=0.1)
	parser.add_argument("-y", "--lambda2", help="Group sparsity parameter.", type=float, default=0.1)
	parser.add_argument("-o", "--output", help="Output directory.", type=str, default="output")
	parser.add_argument("--upsample_balance", help="Balance positive and negative response sets by upsampling the underpopulated set.", action='store_true', default=False)
	parser.add_argument("--downsample_balance", help="Balance positive and negative response sets by downsampling the overpopulated set.", action='store_true', default=False)
	parser.add_argument("--ensemble_parts", help="Build gene-wise ensemble models, splitting the set of genes into N partitions for each run.", type=int, default=None)
	parser.add_argument("--ensemble_coverage", help="Number of ensemble models to build. Each gene will be included in this many individual models.", type=int, default=5)
	parser.add_argument("--sparsify", help="Iteratively increase sparsity until selected set of genes fits in one partition.", action='store_true', default=False)
	parser.add_argument("--method", help="SGLasso type to use. Options are \"leastr\" or \"logistic\". Defaults to \"leastr\".", type=str, default="leastr")
	parser.add_argument("--slep_opts", help="File of tab-separated name-value pairs (one per line) to specify SLEP options.", type=str, default=None)
	args = parser.parse_args()
	score_tables = main(args)
	gene_target = None
	if args.sparsify and score_tables is not None:
		args_original = copy.deepcopy(args)
		output_folder = args.output
		aln_list_dir = os.path.dirname(args.aln_list)
		aln_list_basename = os.path.splitext(os.path.basename(args.aln_list))[0]
		aln_file_list = {}
		with open(args.aln_list, 'r') as file:
			for line in file:
				basename = os.path.splitext(os.path.basename(line.strip()))[0]
				if basename in aln_file_list:
					raise Exception("Found multiple alignment files with identical basename {}.".format(basename))
				aln_file_list[basename] = line.strip()
		for hypothesis in score_tables.keys():
			args = copy.deepcopy(args_original)
			if gene_target is None:
				gene_target = math.ceil(len(score_tables[hypothesis])/args.ensemble_parts)
			elif gene_target != math.ceil(len(score_tables[hypothesis])/args.ensemble_parts):
				raise Exception("Mismatched gene counts between hypotheses.")
			selected_count = sum([1 for val in score_tables[hypothesis].values() if val[0] > 0])
			counter = 1
			final_round = False
			shutil.move("{}_hypothesis.txt".format(hypothesis), "{}.txt".format(hypothesis))
			while selected_count > gene_target or final_round:
				# Generate new aln_list file and point args.aln_list at it
				aln_list_filename = os.path.join(aln_list_dir, "{}_{}_{}.txt".format(aln_list_basename, hypothesis, counter))
				with open(aln_list_filename, 'w') as file:
					for gene in score_tables[hypothesis].keys():
						if score_tables[hypothesis][gene][0] > 0:
							file.write("{}\n".format(aln_file_list[gene]))
				args.aln_list = aln_list_filename
				# Also point --response at this hypothesis' response file
				args.response = "{}.txt".format(hypothesis)
				# args.output = os.path.join(output_folder, "{}_sparsify_round{}".format(hypothesis, counter))
				args.output = "{}_sparsify_round{}".format(hypothesis, counter)
				# temp_score_tables = main(args)
				if final_round:
					args.output = "{}_sparsify_final".format(hypothesis)
					args.lambda1 = 10**-6
					args.lambda2 = 10**-6
					args.ensemble_parts = 1
					args.ensemble_coverage = 1
					args.sparsify = False
				score_tables.update(main(args))
				args.sparsify = True
				shutil.move(args.output, output_folder)
				new_selected_count = sum([1 for val in score_tables[hypothesis].values() if val[0] > 0])
				if not final_round and new_selected_count == selected_count:
					# This round didn't increase sparsity, so increase group sparsity parameter for the next round
					print("The number of selected genes was unchanged, increasing group sparsity parameter from {} to {}.".format(args.lambda2, args.lambda2 * 1.2))
					args.lambda2 = args.lambda2 * 1.2
				else:
					selected_count = new_selected_count
				if selected_count <= gene_target and not final_round:
					final_round = True
				else:
					final_round = False
				counter += 1
			if os.path.exists("{}.txt".format(hypothesis)):
				os.remove("{}.txt".format(hypothesis))

	# generate_gene_prediction_table(weights_filename, responses_filename, groups_filename, features_filename, output_filename)
	# if False:
	# 	pf.generate_gene_prediction_table("angiosperm_out_feature_weights.xml", "angiosperm_20spec_pred.txt", "angiosperm_input/group_indices_angiosperm_input.txt", "angiosperm_input/feature_angiosperm_input.txt", "test_output.txt")
	# if False:
	# 	hypothesis_file_list = pf.generate_hypothesis_set("testtree2.nwk", "nodelist.txt")
	# 	features_filename, groups_filename, response_filename_list = pf.generate_input_matrices("sample_files/angiosperm_100_sample_alns.txt", hypothesis_file_list)
	# 	weights_file_list = pf.run_mlp(features_filename, groups_filename, response_filename_list)
	# 	pf.process_weights(weights_file_list, hypothesis_file_list, groups_filename, features_filename)




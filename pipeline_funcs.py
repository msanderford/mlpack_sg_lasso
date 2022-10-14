import sys
import os
import subprocess
import copy
import time
import random
import shutil
import numpy
import xml.etree.ElementTree as ET
from Bio import Phylo


def find_result_files(args, hypothesis_file_list):
	if args.ensemble_parts is not None and args.ensemble_parts >= 1:
		# Do analysis for ensemble model directory structure
		pass
		result_files = {}
		for hypothesis_file in hypothesis_file_list:
			hypothesis = os.path.splitext(os.path.basename(hypothesis_file))[0].replace("_hypothesis", "")
			result_files[hypothesis] = {}
			for i in range(1, args.ensemble_coverage+1):
				result_files[hypothesis][i] = {}
				for j in range(1, args.ensemble_parts+1):
					result_files[hypothesis][i][j] = {}
					result_files[hypothesis][i][j]["weights"] = os.path.join(args.output,
										     "{}_rep{}_part{}".format(args.output, i, j),
										     "{}_mapped_feature_weights.txt".format(hypothesis))
					result_files[hypothesis][i][j]["predictions"] = os.path.join(args.output,
										     "{}_rep{}_part{}".format(args.output, i, j),
										     "{}_gene_predictions.txt".format(hypothesis))
									     
	else:
		# Do analysis for singular model directory structure
		pass
		result_files = {}
		for hypothesis_file in hypothesis_file_list:
			hypothesis = os.path.splitext(os.path.basename(hypothesis_file))[0].replace("_hypothesis", "")
			result_files[hypothesis] = {}
			result_files[hypothesis]["weights"] = os.path.join(args.output,
									       "{}_mapped_feature_weights.txt".format(hypothesis))
			result_files[hypothesis]["predictions"] = os.path.join(args.output,
									       "{}_gene_predictions.txt".format(hypothesis))
	return result_files


def parse_result_files(args, file_dict):
	weights = {}
	for hypothesis in file_dict.keys():
		temp_weights = {}
		last_gene = ""
		for i in range(1, args.ensemble_coverage+1):
			for j in range(1, args.ensemble_parts+1):
				with open(file_dict[hypothesis][i][j]["weights"], 'r') as file:
					for line in file:
						data = line.strip().split("\t")
						rowname = data[0].split("_")
						feature = "_".join(rowname[-2:])
						gene = data[0].replace("_{}".format(feature), "")
						if i == 1:
							if gene != last_gene:
								last_gene = gene
								temp_weights[gene] = {}
							temp_weights[gene][feature] = []
						temp_weights[gene][feature].append(float(data[1]))
		weights[hypothesis] = temp_weights
	return weights


def analyze_ensemble_weights(args, weights):
	score_tables = {}
	for hypothesis in weights.keys():
		counts = {}
		totals = {}
		scores = {}
		outfile = os.path.join(args.output, "{}_weight_analysis.txt".format(hypothesis))
		for gene in weights[hypothesis].keys():
			counts[gene] = {}
			for feature in weights[hypothesis][gene].keys():
				counts[gene][feature] = sum([1 for x in weights[hypothesis][gene][feature] if abs(x)>0])
			totals[gene] = sum(counts[gene].values())
			scores[gene] = totals[gene]/len(counts[gene])
		with open(outfile, 'w') as file:
			file.write("{}\t{}\t{}\n".format("Gene","Total Non-zero","Score"))
			for gene in weights[hypothesis].keys():
				file.write("{}\t{}\t{}\n".format(gene, totals[gene], scores[gene]))
		score_tables[hypothesis] = {gene:(totals[gene], scores[gene]) for gene in weights[hypothesis].keys()}
	return score_tables


def generate_gene_prediction_table(weights_filename, responses_filename, groups_filename, features_filename, output_filename, gene_list, missing_seqs, group_list, field_filename=None):
	# Read weights, responses, and group indices files
	model = xml_model_to_dict(weights_filename)
	seqlist = []
	responses = {}
	with open(responses_filename, 'r') as file:
		for line in file:
			data = line.strip().split("\t")
			seqlist.append(data[0])
			responses[data[0]] = data[1]
	with open(groups_filename, 'r') as file:
		line1_data = [int(x) for x in file.readline().strip().split(",")]
		line2_data = [int(x) for x in file.readline().strip().split(",")]
		line3_data = [float(x) for x in file.readline().strip().split(",")]
	field = list(range(0, len(model["weight_list"])))
	if field_filename is not None:
		with open(field_filename, 'r') as file:
			field = [int(x)-1 for x in file.readline().strip().split(",") if len(x) > 0]
	group_indices = list(zip(line1_data, line2_data, line3_data))
	group_weights = []
	#print(len(model["weight_list"]))
	#print(len(field))
	for group in group_indices:
		#print(group)
		group_weights.append(numpy.asarray([model["weight_list"][field[x]] for x in range(group[0]-1, group[1])]))
	# Open features file and process 1 row at a time,
	#  can be modified to process in chunks that fit in memory later.
	group_sums = []
	predictions = []
	with open(features_filename, 'r') as file:
		for line in file:
			data = numpy.asarray([float(x) for x in line.strip().split(",")])
			sums = []
			for (index, weight) in zip(group_indices, group_weights):
				#print(index)
				#print(len(data[index[0]-1:index[1]]))
				sums.append(sum(numpy.asarray([data[x] for x in field[index[0]-1:index[1]]]) * weight))
			group_sums.append(sums)
			predictions.append(sum((data * numpy.asarray(model["weight_list"]))) + model["intercept"])
	# Write gene predictions table
	with open(output_filename, 'w') as file:
		file.write("SeqID\tResponse\tPrediction\tIntercept\t{}\n".format("\t".join([",".join(x) for x in group_list])))
		for (seqid, gene_sums, prediction) in zip(seqlist, group_sums, predictions):
			for i in range(0, len(gene_sums)):
				if sum([1 for x in group_list[i] if (seqid, x) in missing_seqs]) == len(group_list[i]):
					gene_sums[i] = "N/A"
			file.write("{}\t{}\t{}\t{}\t{}\n".format(seqid, responses[seqid], prediction, model["intercept"], "\t".join([str(x) for x in gene_sums])))
	with open(str(output_filename).replace("_gene_predictions.txt", "_GSS.txt"), 'w') as file:
		file.write("{}\t{}\n".format("Gene","GSS"))
		for (gene, weights) in zip(gene_list, group_weights):
			file.write("{}\t{}\n".format(gene, str(sum(numpy.abs(weights)))))


# Takes a list of alignment files and splits it into randomly selected subsets of equal size and returns the filenames
def split_gene_list(aln_list_filename, partitions):
	partition_files = []
	dirname = os.path.dirname(aln_list_filename)
	basename = "{}_{}".format(os.path.splitext(os.path.basename(aln_list_filename))[0], random.randrange(100000,999999))	
	with open(aln_list_filename, 'r') as file:
		aln_file_list = [] #List
		partitioned_file_list = [] #List of lists
		for line in file:
			aln_file_list.append(line.strip())
	random.shuffle(aln_file_list)
	partition_sizes = [0 for i in range(0, partitions)]
	for i in range(0, len(aln_file_list)):
		partition_sizes[i % partitions] += 1;
	for i in range(0, partitions):
		filename = "{}_partition_{}.txt".format(basename, i+1)
		filename = os.path.join(dirname, filename)
		with open(filename, 'w') as file:
			for j in range(0, partition_sizes[i]):
				file.write("{}\n".format(aln_file_list.pop()))
		partition_files.append(filename)
	return partition_files


def xml_model_to_dict(model_filename):
	params = {}
	params["intercept"] = 0
	params["lambda1"] = 0
	params["weight_list"] = []
	# Read weights and responses file
	xml_weights = ET.parse(model_filename)
	xml_tree = xml_weights.getroot()
	if xml_tree.tag != "model":
		raise Exception("Unexpected model XML format")
	for child1 in xml_tree:
		if child1.tag == "intercept_value":
			params["intercept"] = float(child1.text)
		elif child1.tag == "lambda1":
			params["lambda1"] = float(child1.text)
		elif child1.tag == "parameters":
			for child2 in child1:
				if child2.tag == "item":
					params["weight_list"].append(float(child2.text))
				else:
					params[child2.tag] = child2.text
	return params


def lookup_by_names(tree):
	names = {}
	for clade in tree.find_clades():
		if clade.name:
			if clade.name in names:
				raise ValueError("Duplicate key: %s" % clade.name)
			names[clade.name] = clade
	return names


def generate_hypothesis_set(newick_filename, nodelist_filename=None, response_filename=None, auto_name_nodes=False, cladesize_cutoff=0, auto_name_length=5):
	tree = Phylo.parse(newick_filename, 'newick').__next__()
	taxa_list = [x.name for x in tree.get_terminals()]
	taxa_list.reverse()
	auto_names = []
	if auto_name_nodes:
		i = 0
		for clade in tree.find_clades():
			if not clade.name:
				new_name = "{}_{}".format(clade[0].get_terminals()[0].name[0:auto_name_length], clade[1].get_terminals()[0].name[0:auto_name_length])
				if new_name in auto_names:
					raise ValueError("Duplicate auto generated name: {}\nIncrease size of auto_name_length parameter and try again.".format(new_name))
				else:
					clade.name = "{}_{}".format(i, new_name)
					auto_names += [new_name]
					i += 1
	nodes = lookup_by_names(tree)
	# print(tree)
	if nodelist_filename is None:
		nodelist = [key for key in nodes if key not in taxa_list]
	else:
		with open(nodelist_filename, 'r') as file:
			nodelist = [line.strip() for line in file]
	nodelist = [x for x in nodelist if len(nodes[x].get_terminals()) >= cladesize_cutoff and len(nodes[x].get_terminals()) < len(taxa_list) ]
	# print(nodelist)
	responses = {}
	if response_filename is None:
		for nodename in nodelist:
			responses[nodename] = {x: -1 for x in taxa_list}
			for terminal in nodes[nodename].get_terminals():
				responses[nodename][terminal.name] = 1
	else:
		with open(response_filename, 'r') as file:
			basename = os.path.splitext(os.path.basename(response_filename))[0]
			responses[basename] = {x: None for x in taxa_list}
			custom_responses = [tuple(line.strip().split("\t")) for line in file]
			for response in custom_responses:
				for terminal in nodes[response[0]].get_terminals():
					if responses[basename][terminal.name] is None:
						responses[basename][terminal.name] = response[1]
					else:
						raise Exception("Response value of sequence {} specified more than once".format(terminal.name))
			for key in responses[basename].keys():
				if responses[basename][key] is None:
					responses[basename][key] = "0"
	for nodename in responses.keys():
		with open("{}_hypothesis.txt".format(nodename), 'w') as file:
			for taxa in taxa_list:
				if responses[nodename][taxa] != "0":
					file.write("{}\t{}\n".format(taxa, responses[nodename][taxa]))
	return ["{}_hypothesis.txt".format(nodename) for nodename in responses.keys()]


def split_path(path):
	path_list = []
	while 1:
		path, dir = os.path.split(path)
		path_list.append(dir)
		if path == "":
			break
	path_list.reverse()
	return path_list


def generate_input_matrices(alnlist_filename, hypothesis_filename_list, args):
	output_basename = args.output
	options = ""
	modified_response = False
	if args.upsample_balance:
		options = "{} {}".format(options.strip(),"ub")
		modified_response = True
	elif args.downsample_balance:
		options = "{} {}".format(options.strip(),"db")
		modified_response = True
	response_file_list = []
	group_indices_file_list = []
	features_file_list = []
	aln_file_list = {}
	gene_list = []
	group_list = []
	field_file_list = []
	# Generate gene list from alignment list file
	with open(alnlist_filename) as file:
		for line in file:
			group = []
			for aln_filename in line.strip().split(","):
				basename = os.path.splitext(os.path.basename(aln_filename.strip()))[0]
				group.append(basename)
				if basename in aln_file_list.keys():
					if aln_file_list[basename] != aln_filename.strip():
						raise Exception("Found multiple alignment files with identical basename {}.".format(basename))
				else:
					aln_file_list[basename] = aln_filename.strip()
			gene_list.extend(group)
			group_list.append(group)
	#preprocess_exe = "/home/tuf79348/git/pipeline/mlpack-3.2.2/build/bin/preprocess"
	#preprocess_exe = os.path.join(os.getcwd(), "mlpack-3.2.2", "build", "bin", "preprocess")
	preprocess_exe = os.path.join(os.getcwd(), "bin", "preprocess")
	preprocess_cwd, alnlist_filename = os.path.split(alnlist_filename)
	if preprocess_cwd == "":
			preprocess_cwd = "."
	if not modified_response:
		# Construct preprocessing command for first hypothesis file
		preprocess_cmd = "{} {} {} {} {}".format(preprocess_exe, os.path.join(os.getcwd(), hypothesis_filename_list[0]), alnlist_filename, output_basename, options)
		print(preprocess_cmd)
		subprocess.call(preprocess_cmd.split(" "), stderr=subprocess.STDOUT, cwd=preprocess_cwd)
		# Move generated inputs to top level directory
		if preprocess_cwd != ".":
			shutil.move(os.path.join(preprocess_cwd, output_basename), ".")
		# Substitute gene/group penalties
		if args.gene_penalties is not None:
			gene_penalties = {}
			penalty_list = []
			with open(os.path.join(output_basename, "group_indices_" + output_basename + ".txt"), 'r') as groups_file:
				lines = [val.strip() for val in groups_file.readlines()]
			with open(args.gene_penalties, 'r') as penalties_file:
				for line in penalties_file.readlines():
					data = line.strip().split('\t')
					gene_penalties[data[0]] = float(data[1])
			with open(os.path.join(preprocess_cwd, alnlist_filename), 'r') as file:
				for line in file.readlines():
					try:
						penalty_list.append(sum([gene_penalties[os.path.splitext(os.path.basename(filename))[0]] for filename in line.strip().split(',')]))
					except:
						raise Exception("Could not find gene penalty for one of the following genes: {}".format(','.join([os.path.splitext(os.path.basename(filename))[0] for filename in line.strip().split(',')])))
			with open(os.path.join(output_basename, "group_indices_" + output_basename + ".txt"), 'w') as groups_file:
				groups_file.write("{}\n".format(lines[0]))
				groups_file.write("{}\n".format(lines[1]))
				groups_file.write("{}\n".format('\t'.join([str(x) for x in penalty_list])))
		# Construct response input file for each additional hypothesis file
		for filename in hypothesis_filename_list:
			with open(filename, 'r') as infile:
				temp_fname = os.path.join(output_basename, "response_" + os.path.splitext(os.path.basename(filename))[0] + ".txt")
				with open(temp_fname, 'w') as outfile:
					for line in infile:
						outfile.write("{}\n".format(line.strip().split("\t")[1]))
					response_file_list.append(temp_fname)
					group_indices_file_list.append(os.path.join(output_basename, "group_indices_" + output_basename + ".txt"))
					features_file_list.append(os.path.join(output_basename, "feature_" + output_basename + ".txt"))
					field_file_list.append(os.path.join(output_basename, "field_" + output_basename + ".txt"))
		#return [os.path.join(output_basename, "feature_" + output_basename + ".txt"), os.path.join(output_basename, "group_indices_" + output_basename + ".txt"), response_file_list, gene_list]
		return [features_file_list, group_indices_file_list, response_file_list, gene_list, field_file_list, group_list]
	else:
		for filename in hypothesis_filename_list:
			# Construct preprocessing command
			preprocess_cmd = "{} {} {} {} {}".format(preprocess_exe, os.path.join(os.getcwd(), filename), alnlist_filename, output_basename, options)
			# print(preprocess_cmd)
			subprocess.call(preprocess_cmd.split(" "), stderr=subprocess.STDOUT, cwd=preprocess_cwd)
			if preprocess_cwd != ".":
				if os.path.exists(output_basename):
					shutil.copytree(os.path.join(preprocess_cwd, output_basename), output_basename, dirs_exist_ok=True)
					shutil.rmtree(os.path.join(preprocess_cwd, output_basename))
				else:
					shutil.move(os.path.join(preprocess_cwd, output_basename), output_basename)
			hypothesis_basename = os.path.splitext(os.path.basename(filename))[0]
			shutil.move(os.path.join(output_basename, "feature_" + output_basename + ".txt"), os.path.join(output_basename, "feature_" + hypothesis_basename + ".txt"))
			shutil.move(os.path.join(output_basename, "group_indices_" + output_basename + ".txt"), os.path.join(output_basename, "group_indices_" + hypothesis_basename + ".txt"))
			shutil.move(os.path.join(output_basename, "response_" + output_basename + ".txt"), os.path.join(output_basename, "response_" + hypothesis_basename + ".txt"))
			shutil.move(os.path.join(output_basename, "field_" + output_basename + ".txt"), os.path.join(output_basename, "field_" + hypothesis_basename + ".txt"))
			shutil.move(os.path.join(output_basename, "feature_mapping_" + output_basename + ".txt"), os.path.join(output_basename, "feature_mapping_" + hypothesis_basename + ".txt"))
			shutil.move(os.path.join(output_basename, "resampled_" + output_basename + ".txt"), os.path.join(output_basename, "resampled_" + hypothesis_basename + ".txt"))
			response_file_list.append(os.path.join(output_basename, "response_" + hypothesis_basename + ".txt"))
			group_indices_file_list.append(os.path.join(output_basename, "group_indices_" + hypothesis_basename + ".txt"))
			field_file_list.append(os.path.join(output_basename, "field_" + hypothesis_basename + ".txt"))
			features_file_list.append(os.path.join(output_basename, "feature_" + hypothesis_basename + ".txt"))
		return [features_file_list, group_indices_file_list, response_file_list, gene_list, field_file_list, group_list]


def run_mlp(features_filename_list, groups_filename_list, response_filename_list, field_filename_list, sparsity, group_sparsity, method, slep_opts):
	if method == "leastr":
		method = "sg_lasso_leastr"
	elif method == "logistic":
		method = "sg_lasso"
	elif method == "ol_leastr":
		method = "overlapping_sg_lasso_leastr"
	else:
		raise Exception("Provided method name not recognized, please provide a valid method name.")
	weights_file_list = []
	mlp_exe = os.path.join(os.getcwd(), "bin", method)
	# Run sg_lasso for each response file in response_filename_list
	for response_filename, features_filename, groups_filename, field_filename in zip(response_filename_list, features_filename_list, groups_filename_list, field_filename_list):
		basename = str(os.path.splitext(os.path.basename(response_filename))[0]).replace("response_","")
		if slep_opts is None:
			mlp_cmd = "{} -f {} -z {} -y {} -n {} -r {} -w {}".format(mlp_exe, features_filename, sparsity, group_sparsity, groups_filename, response_filename, basename + "_out_feature_weights")
		else:
			mlp_cmd = "{} -f {} -z {} -y {} -n {} -r {} -s {} -w {}".format(mlp_exe, features_filename, sparsity, group_sparsity, groups_filename, response_filename, slep_opts, basename + "_out_feature_weights")
		if method == "overlapping_sg_lasso_leastr":
			mlp_cmd = mlp_cmd + " -g {}".format(field_filename)
		print(mlp_cmd)
		subprocess.call(mlp_cmd.split(" "), stderr=subprocess.STDOUT)
		weights_file_list.append(basename + "_out_feature_weights.xml")
	return weights_file_list


def process_weights(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list, gene_list, HSS, missing_seqs, group_list):
	for (weights_filename, hypothesis_filename, groups_filename, features_filename) in zip(weights_file_list, hypothesis_file_list, groups_filename_list, features_filename_list):
		generate_gene_prediction_table(weights_filename, hypothesis_filename, groups_filename, features_filename, hypothesis_filename.replace("_hypothesis.txt", "_gene_predictions.txt"), gene_list, missing_seqs, group_list, groups_filename.replace("group_indices_", "field_"))
		total_significance = generate_mapped_weights_file(weights_filename, groups_filename.replace("group_indices_", "feature_mapping_"))
		HSS[hypothesis_filename] = HSS.get(hypothesis_filename, 0) + total_significance

def generate_mapped_weights_file(weights_filename, feature_map_filename):
	# Read weights and feature mapping files
	PSS = {}
	posname_list = []
	last_posname = ""
	model = xml_model_to_dict(weights_filename)
	feature_map = {}
	output_filename = str(weights_filename).replace("_hypothesis_out_feature_weights.xml", "_mapped_feature_weights.txt")
	with open(feature_map_filename, 'r') as file:
		for line in file:
			data = line.strip().split("\t")
			if len(data) == 2:
				feature_map[int(data[0])] = data[1]
	with open(output_filename, 'w') as file:
		for i in range(0, len(model["weight_list"])):
			file.write("{}\t{}\n".format(feature_map[i+1], model["weight_list"][i]))
			posname = feature_map[i + 1][0:-2]
			if posname != last_posname:
				posname_list.append(posname)
				last_posname = posname
			PSS[posname] = PSS.get(posname, 0.0) + abs(model["weight_list"][i])
	with open(str(output_filename).replace("_mapped_feature_weights.txt", "_PSS.txt"), 'w') as file:
		file.write("{}\t{}\n".format("Position Name", "PSS"))
		for posname in posname_list:
			file.write("{}\t{}\n".format(posname, PSS[posname]))
	# Return sum of all position significance scores
	return sum(PSS.values())


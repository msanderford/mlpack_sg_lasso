import sys
import os
import subprocess
import copy
import time
import shutil
import numpy
import xml.etree.ElementTree as ET
from Bio import Phylo


def generate_gene_prediction_table(weights_filename, responses_filename, groups_filename, features_filename, output_filename, gene_list):
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
		line1_data = [int(x) for x in file.readline().strip().split("\t")]
		line2_data = [int(x) for x in file.readline().strip().split("\t")]
		line3_data = [float(x) for x in file.readline().strip().split("\t")]
	group_indices = list(zip(line1_data, line2_data, line3_data))
	group_weights = []
	# print(len(model["weight_list"]))
	for group in group_indices:
		# print(group)
		group_weights.append(numpy.asarray([model["weight_list"][x] for x in range(group[0]-1, group[1])]))
	# Open features file and process 1 row at a time,
	#  can be modified to process in chunks that fit in memory later.
	group_sums = []
	with open(features_filename, 'r') as file:
		for line in file:
			data = numpy.asarray([float(x) for x in line.strip().split("\t")])
			sums = []
			for (index, weight) in zip(group_indices, group_weights):
				sums.append(sum(data[index[0]-1:index[1]] * weight))
			group_sums.append(sums)
	# Write gene predictions table
	with open(output_filename, 'w') as file:
		file.write("SeqID\tResponse\tPrediction\t{}\n".format("\t".join(gene_list)))
		for (seqid, gene_sums) in zip(seqlist, group_sums):
			file.write("{}\t{}\t{}\t{}\n".format(seqid, responses[seqid], sum(gene_sums) + model["intercept"], "\t".join([str(x) for x in gene_sums])))


def xml_model_to_dict(model_filename):
	params = {}
	params["intercept"] = 0
	params["lambda1"] = 0
	params["weight_list"] = []
	# Read weights and responses file
	xml_weights = ET.parse(model_filename)
	xml_tree = xml_weights.getroot()
	if xml_tree[0].tag != "model":
		raise Exception("Unexpected model XML format")
	for child1 in xml_tree[0]:
		if child1.tag == "intercept":
			params["intercept"] = float(child1.text)
		elif child1.tag == "lambda1":
			params["lambda1"] = float(child1.text)
		elif child1.tag == "parameters":
			for child2 in child1:
				if child2.tag == "item":
					params["weight_list"].append(float(child2.text))
				else:
					params[child2.tag] = child2.text
	params["intercept"] = 0
	return params


def lookup_by_names(tree):
	names = {}
	for clade in tree.find_clades():
		if clade.name:
			if clade.name in names:
				raise ValueError("Duplicate key: %s" % clade.name)
			names[clade.name] = clade
	return names


def generate_hypothesis_set(newick_filename, nodelist_filename=None, response_filename=None):
	tree = Phylo.parse(newick_filename, 'newick').__next__()
	nodes = lookup_by_names(tree)
	taxa_list = [x.name for x in tree.get_terminals()]
	taxa_list.reverse()
	# print(tree)
	if nodelist_filename is None:
		nodelist = [key for key in nodes if key not in taxa_list]
	else:
		with open(nodelist_filename, 'r') as file:
			nodelist = [line.strip() for line in file]
	# print(nodelist)
	responses = {}
	if response_filename is None:
		for nodename in nodelist:
			responses[nodename] = {x: -1 for x in taxa_list}
			for terminal in nodes[nodename].get_terminals():
				responses[nodename][terminal.name] = 1
	else:
		with open(response_filename, 'r') as file:
			responses["custom"] = {x: None for x in taxa_list}
			custom_responses = [tuple(line.strip().split("\t")) for line in file]
			for response in custom_responses:
				for terminal in nodes[response[0]].get_terminals():
					if responses["custom"][terminal.name] is None:
						responses["custom"][terminal.name] = response[1]
					else:
						raise Exception("Response value of sequence {} specified more than once".format(terminal.name))
			for key in responses["custom"].keys():
				if responses["custom"][key] is None:
					responses["custom"][key] = 0
	for nodename in responses.keys():
		with open("{}_hypothesis.txt".format(nodename), 'w') as file:
			for taxa in taxa_list:
				if responses[nodename][taxa] != 0:
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


def generate_input_matrices(alnlist_filename, hypothesis_filename_list, output_basename):
	response_file_list = []
	gene_list = []
	# Generate gene list from alignment list file
	with open(alnlist_filename) as file:
		for line in file:
			gene_list.append(os.path.splitext(os.path.basename(line.strip()))[0])
	# Construct preprocessing command for first hypothesis file
	preprocess_exe = "/home/tuf79348/git/pipeline/mlpack-3.2.2/build/bin/preprocess"
	preprocess_exe = os.path.join(os.getcwd(), "mlpack-3.2.2", "build", "bin", "preprocess")
	preprocess_cwd, alnlist_filename = os.path.split(alnlist_filename)
	options = ""
	preprocess_cmd = "{} {} {} {} {}".format(preprocess_exe, os.path.join(os.getcwd(), hypothesis_filename_list[0]), alnlist_filename, output_basename, options)
	# print(preprocess_cmd)
	if preprocess_cwd == "":
		preprocess_cwd = "."
	subprocess.call(preprocess_cmd.split(" "), stderr=subprocess.STDOUT, cwd=preprocess_cwd)
	# Move generated inputs to top level directory
	if preprocess_cwd != ".":
		shutil.move(os.path.join(preprocess_cwd, output_basename), ".")
	# Construct response input file for each additional hypothesis file
	for filename in hypothesis_filename_list:
		with open(filename, 'r') as infile:
			temp_fname = os.path.join(output_basename, "response_" + os.path.splitext(os.path.basename(filename))[0] + ".txt")
			with open(temp_fname, 'w') as outfile:
				for line in infile:
					outfile.write("{}\n".format(line.strip().split("\t")[1]))
				response_file_list.append(temp_fname)
	return [os.path.join(output_basename, "feature_" + output_basename + ".txt"), os.path.join(output_basename, "group_indices_" + output_basename + ".txt"), response_file_list, gene_list]


def run_mlp(features_filename, groups_filename, response_filename_list):
	weights_file_list = []
	mlp_exe = os.path.join(os.getcwd(), "mlpack-3.2.2", "build", "bin", "mlpack_sg_lasso_leastr")
	# Run sg_lasso for each response file in response_filename_list
	for response_filename in response_filename_list:
		basename = str(os.path.splitext(os.path.basename(response_filename))[0]).replace("response_","")
		mlp_cmd = "{} -v -f {} -z 0.1 -y 0.5 -n {} -r {} -w {}".format(mlp_exe, features_filename, groups_filename, response_filename, basename + "_out_feature_weights.xml")
		print(mlp_cmd)
		subprocess.call(mlp_cmd.split(" "), stderr=subprocess.STDOUT)
		weights_file_list.append(basename + "_out_feature_weights.xml")
	return weights_file_list


def process_weights(weights_file_list, hypothesis_file_list, groups_filename, features_filename, gene_list):
	for (weights_filename, hypothesis_filename) in zip(weights_file_list, hypothesis_file_list):
		generate_gene_prediction_table(weights_filename, hypothesis_filename, groups_filename, features_filename, str(hypothesis_filename).replace("_hypothesis.txt", "_gene_predictions.txt"), gene_list)
		generate_mapped_weights_file(weights_filename, str(groups_filename).replace("group_indices_", "feature_mapping_"))

def generate_mapped_weights_file(weights_filename, feature_map_filename):
	# Read weights and feature mapping files
	model = xml_model_to_dict(weights_filename)
	feature_map = {}
	output_filename = str(weights_filename).replace("_hypothesis_out_feature_weights.xml", "_mapped_feature_weights.txt")
	with open(feature_map_filename, 'r') as file:
		for line in file:
			data = line.strip().split("\t")
			if len(data) == 2:
				feature_map[int(data[0])] = data[1]
	with open(output_filename, 'w') as file:
		for i in range(1, len(model["weight_list"])):
			file.write("{}\t{}\n".format(feature_map[i], model["weight_list"][i]))


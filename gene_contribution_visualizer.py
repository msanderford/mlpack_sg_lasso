import os
import argparse
import shutil
import math
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import TwoSlopeNorm
import numpy as np


def main(predictions_table, lead_cols=4, response_idx=2, prediction_idx=3, output=None):
	sample_size = 24
	sample_size = None
	with open(predictions_table, 'r') as file:
		header = file.readline().split("\t")
		num_cols = len(header)
		file.seek(0)
		if sample_size is None:
			data = np.genfromtxt(file, dtype=None, skip_header=1, missing_values="N/A", encoding=None)
		else:
			data = np.genfromtxt(file, dtype=None, skip_header=1, missing_values="N/A", encoding=None, usecols=range(0, lead_cols + sample_size))
			header = header[0: lead_cols + sample_size]
		species = [(val)[0] for val in data]
		data = np.asarray([list(val)[1:] for val in data], dtype="float")
	num_rows, num_cols = data.shape
	print("\nRows: {}\nCols: {}\n".format(num_rows, num_cols))
	#data = np.random.rand(10, 10) * 2 - 1


	# Sort columns by sum of contributions
	ssq_scores = list(zip(range(3, num_cols), [np.sum(gene_col**2) for gene_col in data.transpose()[3:]]))
	ssq_scores.sort(key=lambda tup: tup[1], reverse=True)
	sum_scores = list(zip(range(3, num_cols), [np.sum(gene_col) for gene_col in data.transpose()[3:]]))
	sum_scores.sort(key=lambda tup: tup[1], reverse=True)
	data = data[:, list(range(0, 3)) + [val[0] for val in ssq_scores]]


	# Sort rows by ground truth and predicted value
	sorted_rows = list(zip(range(0, num_rows), data[:, 0], data[:, 1]))
	sorted_rows.sort(key=lambda tup: (tup[1], tup[2]), reverse=True)
	data = data[[val[0] for val in sorted_rows]]


	# draw gridlines
	xtick_width = 20
	ytick_width = 20
	#fig, ax = plt.subplots()
	fig, ax = plt.subplots(figsize=(num_cols * 0.1, num_rows * 0.1))
	#fig.set_size_inches()
	xlabel_size = 3
	ylabel_size = 3
	cell_label_size = 2.0
	DPI = 600
	ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=0.2)
	# Make red-yellow-green color map with NaN=blue
	cmap = copy.copy(mpl.cm.get_cmap("RdYlGn"))
	cmap.set_bad(color='lightskyblue')
	# Apply color map to first three columns, and then separately to the rest of the columns
	ax.imshow(data[:, 0:3], cmap=cmap, norm=TwoSlopeNorm(0), extent=(0, (lead_cols - 1) * xtick_width, num_rows * ytick_width, 0))
	ax.imshow(data[:, 3:], cmap=cmap, norm=TwoSlopeNorm(0), extent=((lead_cols - 1)*xtick_width, num_cols*xtick_width, num_rows*ytick_width, 0))
	#ax.set_xticks(np.arange(-.5, 100, 10));
	ax.set_xticks(np.arange(0, num_cols * xtick_width, xtick_width));
	ax.set_xticklabels(header[1:], rotation=90, ha='left', size=xlabel_size)
	#ax.set_yticks(np.arange(-.5, 10, 1));
	ax.set_yticks(np.arange(0, num_rows * ytick_width, ytick_width));
	ax.set_yticklabels(species, va="top", size=ylabel_size)
	for (i, j), z in np.ndenumerate(data):
		ax.text((j+0.5) * xtick_width, (i+0.5)*ytick_width, '{:0.2f}'.format(z), ha='center', va='center', size=cell_label_size)

	if output is None:
		output = "{}.png".format(os.path.splitext(os.path.basename(predictions_table))[0])
	plt.savefig(output, dpi=DPI, bbox_inches='tight')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Gene contribution visualizer for ESL pipeline.")
	parser.add_argument("predictions_table", help="Table of gene prediction contribution values to sort and color.", type=str)
	parser.add_argument("--lead_cols", help="Number of columns ahead of gene contribution columns.", type=int, default=4)
	parser.add_argument("--response_idx", help="1-based index of response column.", type=int, default=2)
	parser.add_argument("--prediction_idx", help="1-based index of prediction column.", type=int, default=3)
	parser.add_argument("--output", help="Output image file.", type=str, default=None)
	args = parser.parse_args()
	main(args.predictions_table, args.lead_cols, args.response_idx, args.prediction_idx, args.output)
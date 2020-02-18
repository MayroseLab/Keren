from ete3 import Tree
import re, os

#### auxiliary functions for converting tree string to INDELible's readable format ####

def reformat_num(number):
	num = float(number)
	reformatted_num = str("{0:.10f}".format(num))
	return reformatted_num

def reformat_branch(branch):
	fixed_branch = branch[0] + ":" + branch[1:]
	return fixed_branch


#### utils in alterations and learining of trees ####

# measures the total branch length of the tree
def get_tree_size(tree):
	dem = Tree()
	if type(tree) != type(dem):
		tree = Tree(tree)
	size = 0
	for node in tree.traverse():
		size+= node.dist
	return size

# shrinks or blow the tree by multiplying its branches by a constant
def rescale_tree(orig_tree, size, new_tree_path):
	tree = Tree(orig_tree) # orig_tree could be a path to a tree file of a string of a tree
	orig_size = get_tree_size(tree)
	factor = size / orig_size
	for node in tree.traverse():
		node.dist = node.dist * factor
	tree.write(outfile=new_tree_path)
	return 0

# converts a tree string to an INDElible readable format
def adjust_tree(tree):

	tree_file = open(tree, 'r')
	tree_str = tree_file.read()
	tree_file.close()

	# if there is ROOT string in the end of the tree -> remove it
	tree_str = tree_str.replace('ROOT', '')

	# remove internal node names
	internal_nodes_pattern = 'N[0-9]+'
	pattern = re.compile(internal_nodes_pattern)
	tree_str = pattern.sub('', tree_str)

	# reformat small numbers
	reformatted_number = re.compile(r'\d+[\.\d+]*e\-\d+')

	for number in reformatted_number.findall(tree_str):
		reformatted_num = reformat_num(number)
		tree_str = tree_str.replace(number, reformatted_num, 1)

	comb_branch_pattern = re.compile('\)[0-9]+\.[0-9]+')
	for branch in re.findall(comb_branch_pattern, tree_str):
		reformatted_branch = reformat_branch(branch)
		tree_str = re.sub(re.escape(branch), reformatted_branch, tree_str)

	#remove bootsrap numbers
	bootsrap_num_expr = re.compile(r"-{0,1}[0-9]+\.[0-9]+:")
	tree_str = bootsrap_num_expr.sub("", tree_str)

	return tree_str

# auxiliary function to reduce a tree with a given subset of taxa
def reduce_tree(input_tree, output_tree, subset):
	tree = Tree(input_tree, format=1)
	if "NA" in subset:
		subset.remove("NA")
	tree.prune(subset, preserve_branch_length=True)
	tree.write(outfile=output_tree, format=5)  # do not write to the file the internal nodes names
	return 0

# auxiliary function to remove branches length from a tree
def remove_branches_length(in_tree):
	out_tree = in_tree + ".no_branches"
	tree = Tree(in_tree)
	tree.write(format=9, outfile=out_tree)
	res = os.system("rm -r " + in_tree)
	res = os.rename(out_tree, in_tree)
	return 0

def rename_nodes(input_tree_path, output_tree_path, prefix="S"):
	tree = Tree(input_tree_path)
	for node in tree.traverse():
		node.name = prefix + node.name
	tree.write(outfile=output_tree_path, format=5)
	return 0

import os, socket
from ete3 import Tree



#### nj auxiliary functions ####

# extent sequences names in a phylip file to 11 chars exactly, and writes the updated data to output
def add_sapces(tree_name, msa_path, output_dir):
    output = output_dir + tree_name + "_protdist_compatible_msa.phy"
    name_expr = re.compile(r"N\d+[ ]+")
    with open(msa_path, "r") as in_file:
        content = in_file.read()
    for name in name_expr.findall(content):
        curr_len = len(name)
        new_name = name
        while curr_len < 10:
            new_name = new_name + " "
            curr_len = curr_len + 1
        content = content.replace(name, new_name)
    with open(output, "w") as out_file:
        out_file.write(content)
    return output

# create distance matrix
def generate_dist_mat(tree_name, msa_file, output_dir):
    protdist_file = output_dir + tree_name + "_distance_matrix"
    res = os.system('echo "' + msa_file + '\nF\n'+ protdist_file +'\nY\n" | protdist')  #for my file
    return protdist_file


#### generating trees functions ####

# create nj tree
def create_nj_tree(tree_name, msa_path, nj_garbage_dir, nj_tree_path):
	# set garbage directories
	msa_files_dir = nj_garbage_dir + "protdist_compatible_msa/"
	if not os.path.exists(msa_files_dir):
		res = os.system("mkdir " + msa_files_dir)
	dist_mat_dir = nj_garbage_dir + "protdist_distance_matrices/"
	if not os.path.exists(dist_mat_dir):
		res = os.system("mkdir " + dist_mat_dir)

	# generate helper files
	msa_file = add_sapces(tree_name, msa_path, msa_files_dir)
	distance_matrix_file = generate_dist_mat(tree_name, msa_file, dist_mat_dir)

	# create a folder for PHYLIP to generate the tree in
	tree_folder = nj_garbage_dir + tree_name + "/"
	if not os.path.isdir(tree_folder):
		res = os.system("mkdir " + tree_folder)

	# call PHYLIP to generate the tree
	output1 = tree_folder + "nj_data_" + tree_name
	output2 = tree_folder + "nj_tree_" + tree_name
	res = os.system("cd " + tree_folder)
	res = os.system('echo "' + distance_matrix_file + '\nF\n' + output1 + '\nY\nF\n' + output2 + '\n" | neighbor')

	# get the output to the nj_tree_path
	res = os.rename(output2, nj_tree_path)

	# remove the tree folder
	res = os.system("rm -r " + tree_folder)

	return 0

# create FastTree tree
def create_fasttree_tree(msa_path, tree_path):
    cmd = 'FastTree < ' + msa_path + ' > ' + tree_path
    res = os.system(cmd)
    if res != 0:
        return -1
    # convert to format r4s likes
    tree = Tree(tree_path, format=1)
    tree.write(outfile=tree_path, format=5)
    return 0

# create RAxML tree
def create_raxml_tree(msa_path, output_dir, output_name):
    tree_path = output_dir+ "/RAxML_bestTree." + output_name
    if not(os.path.exists(tree_path)):
        res = os.system("raxmlHPC -m PROTCAT -s " + msa_path + " -w " + output_dir + " -n " + output_name + " -p 1")
        if res != 0:
            return -1
    return tree_path

# create EAxML tree
def create_eaxml_tree(msa_path, base_tree, output_dir, output_name):
    tree_path = output_dir+ "/ExaML_result." + output_name
    if not (os.path.exists(tree)):
        bin_name = "bin_" + output_name
        res = os.system("parse-examl -s " + msa_path + " -m PROT -n " + bin_name)
        if res != 0:
            return -1
        res = os.system("examl -s " + output_dir + bin_name + ".binary -t " + base_tree + " -m " + model.rstrip() + " -w " + output_dir + " -n " + output_name)
        if res != 0:
            return -1
    return tree_path

# create UPGMA tree
def create_upgma_tree(fasta_path, output_path):
    msa_path = output_path + ".aln"
    res = os.system("mafft --treeout --auto " + fasta_path + " > " + msa_path)
    if mafft_result != 0:
        return -1
    tree_path = fasta_path + ".tree"
    with open(tree_path, "r") as tree_file:
        tree_str = tree_file.read()
    redundant_patterns = ['\([0-9]+_+', ',[0-9]+_+', '_']
    for pattern in redundant_patterns:
        tree_str = re.sub(pattern, '(', tree_str)
    with open(output_path, "w") as output_file:
        output_file.write(tree_str)
    return 0

# create BIONJ tree
def create_bionj_tree(msa_path, output_path):
    hostname = socket.gethostname()
    print("hostname="+hostname)
    if hostname == 'jekyl.tau.ac.il':
        print("error, phyml is absent from jekyl")
        return -1
    else:
        res = os.system("phyml -i " + msa_path + " -d aa -n 1 -m JTT -o n")
        if res != 0:
            return -1
        tree_path = msa_path + "_phyml_tree.txt"
        res = os.system("cp " + tree_path + " " + output_path)
        return 0
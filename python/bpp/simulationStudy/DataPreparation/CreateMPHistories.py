import re, os, argparse
from time import sleep
from ete3 import Tree
import numpy as np

def size(tree):
    s = 0
    for node in tree.traverse():
        s += node.dist
    return s

def reroot(history_path, tree_path):
    history = Tree(history_path, format=1)
    tree = Tree(tree_path)
    if len(tree.get_children()) > 2:
        print("original tree is not rooted, no need to fix history")
        return
    child1 = history.get_children()[0].detach()
    child2 = history.get_children()[1].detach()
    missing_length = tree.get_children()[0].dist
    history.get_tree_root().add_child(name="missing_node{0}", dist=missing_length)
    history.get_children()[-1].add_child(child1)
    history.get_children()[-1].add_child(child2)
    # make sure that now the tree and the history are of the same length
    if not size(history) == size(tree):
        print("Error! failed to fix history tree")
        print("size(history) = ", size(history) , "\n size(tree) = ", size(tree))
        exit(1)
    # make sure that all the lengths were written to tree string
    history_str = history.write(outfile=None, format=1)
    bl_re = re.compile("\:(\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    bls = [float(match.group(1)) for match in bl_re.finditer(history_str)]
    if abs(np.sum(bls)-size(tree)) > 0.00001:
        print("Error! failed to fix history tree newick format")
        print("size(history) = ", size(history) , "\n size(tree) = ", np.sum(bls))
        exit(1)
    history.write(outfile=history_path, format=1)

if __name__ == '__main__':
    # process input from command line
    parser = argparse.ArgumentParser(
        description='Creates maximum parsimony histories and parameters files')
    parser.add_argument('--trees_dir', '-t', help='path to the input trees for the mp solution program', required=True)
    parser.add_argument('--character_data_dir', '-c',
                        help='directory that holds trees and character data', required=True)
    parser.add_argument('--parameter_files_dir', '-o',
                        help='directory to write the parameter files to',
                        required=True)
    parser.add_argument('--relax_parameter_files_dir', '-p', help='directory of relax parameter files', required=True)
    parser.add_argument('--use_mp', '-u', help='1 for mp histories, 0 for for ml histories', required=False, default=1)

    args = parser.parse_args()
    trees_dir = args.trees_dir
    character_data_dir = args.character_data_dir
    relax_parameter_files_dir = args.relax_parameter_files_dir
    parameter_files_dir = args.parameter_files_dir
    if not os.path.exists(parameter_files_dir):
        res = os.system("mkdir -p " + parameter_files_dir)
    use_mp = int(args.use_mp)

    # # get the tree height
    # tbl_regex = re.compile("tbl_(\d*\.?\d*)")
    # tbl = tbl_regex.search(character_data_dir).group(1)

    replicate_regex = re.compile("replicate_(\d*)")
    for path in os.listdir(character_data_dir):
        if "replicate" in path and os.path.isdir(character_data_dir+path):
            replicate = replicate_regex.search(path).group(1)
            
            hist_prefix = "mp"
            if use_mp == 0:
                hist_prefix = "ml"

            # create input parameters file for the program
            character_data_path = character_data_dir + path + "/character_data/character_data.fas"
            if not os.path.exists(character_data_path):
                continue
            tree_path = trees_dir + replicate + ".nwk"
            hist_dir = character_data_dir + path + "/" + hist_prefix + "_data/"
            if not os.path.exists(hist_dir):
                res = os.system("mkdir -p " + hist_dir)
            history_path = hist_dir + hist_prefix + "_history.nwk"
            lables_path = hist_dir + "mp_partition.txt"
            with open(hist_dir + hist_prefix + "_solution_parameters.bpp", "w") as outfile:
                outfile.write("input.character.file = " + character_data_path + "\n")
                outfile.write("init.tree = user\ninput.tree.file = " + tree_path + "\n")
                outfile.write("input.tree.format = Newick\ninit.brlen.method = Input\n")
                outfile.write("output.tree.file = " + history_path + "\n")
                outfile.write("output.mp.partition = " + lables_path + "\n")

            # execute the program on the parameters file
            script_path = "/groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite/writemphistory"
            if use_mp == 0:
                script_path = "/groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite/writeasrsolution"
            res = os.system(script_path + " param=" + hist_dir + hist_prefix + "_solution_parameters.bpp")
            while not os.path.exists(lables_path):
               sleep(2)

            # re-rot the tree according to the input tree
            reroot(history_path, tree_path)

            # open mp history and write history tree (should be exactly like the original tree) to a file
            history = Tree(history_path, format=1)
            history_str = history.write(outfile=None, format=5)
            history_str = history_str.replace("{0}", "")
            history_str = history_str.replace("{1}", "")
            history_tree_path = history_path.replace("mp_history.nwk", "mp_tree.nwk")
            with open(history_tree_path, "w") as history_tree_file:
               history_tree_file.write(history_str)
            tree_path_info = "input.tree.file = " + tree_path + "\n"

            # parse the labels of the history
            with open(lables_path, "r") as infile:
                labels_data = infile.readlines()
                BG_labels = "model1.nodes_id = " + labels_data[0]
                FG_labels = "model2.nodes_id = " + labels_data[1]

            # now create the parameters file
            BG_labels_regex = re.compile("model1.nodes_id = .*\n")
            FG_labels_regex = re.compile("model2.nodes_id = .*\n")
            tree_path_regex = re.compile("input.tree.file = .*?\n")
            relax_parameter_files_path = relax_parameter_files_dir + replicate + ".bpp"
            with open(relax_parameter_files_path, "r") as infile:
                content = infile.read()
            # replace the labels of the true history with labels of the mp history
            # content = content.replace("/mu_", "/tbl_" + tbl + "_mu_")
            content = re.sub("model1.nodes_id = .*\n", BG_labels, content)
            content = re.sub("model2.nodes_id = .*\n*", FG_labels, content)
            content = re.sub("input.tree.file = .*?\n", tree_path_info, content)
            # replace the directory of the history tree with the one of the tree
            with open(parameter_files_dir + replicate + ".bpp", "w") as outfile:
                outfile.write(content)





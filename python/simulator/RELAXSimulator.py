import argparse, os, re
from ete3 import Tree

''' how to use rhe simulator:
# cd to the directory to create the sequence data in (because indelible will write to the directory from which the script is executed
cd "/groups/itay_mayrose/halabikeren/TraitRELAX/bpp/simulations/indeliable_simulations/TraitRELAX/mu_10_pi0_0.5_kappa_1_omega0_0.1_omega1_1_omega2_2_theta1_0.5_theta2_0.8_k_1/seq_data/"

# execute the script
python "/groups/itay_mayrose/halabikeren/myScripts/python/simulator/TraitRELAXSimulator.py" -o /groups/itay_mayrose/halabikeren/TraitRELAX/bpp/simulations/indeliable_simulations/TraitRELAX/mu_10_pi0_0.5_kappa_1_omega0_0.1_omega1_1_omega2_2_theta1_0.5_theta2_0.8_k_1/ -t /groups/itay_mayrose/halabikeren/TraitRELAX/bpp/simulations/indeliable_simulations/TraitRELAX/base_tree.nwk -k 1

# remove spaces from the names of the taxa in the simulated MSAs because Bio++ doesn't like spaces
python "/groups/itay_mayrose/halabikeren/myScripts/python/simulator/removeSpacesFromFiles.py" -i "/groups/itay_mayrose/halabikeren/TraitRELAX/bpp/simulations/indeliable_simulations/TraitRELAX/mu_10_pi0_0.5_kappa_1_omega0_0.1_omega1_1_omega2_2_theta1_0.5_theta2_0.8_k_1/seq_data/"

'''


def fix_tree_format(tree_path):
    tree = Tree(tree_path, format=1)
    tree_str = tree.write(outfile=None, format=1)
    bad_format_numbers_regex = re.compile(":(\d*\.?\d*e-\d*)", re.MULTILINE | re.DOTALL)
    for match in bad_format_numbers_regex.finditer(tree_str):
        bad_format = match.group(1)
        good_format = "{:.10f}".format(float(bad_format), 10)
        tree_str = tree_str.replace(bad_format, good_format)
    with open(tree_path, "w") as tree_file:
        tree_file.write(tree_str)

def fix_tree_str_format(tree_str):
    bad_format_numbers_regex = re.compile(":(\d*\.?\d*e-\d*)", re.MULTILINE | re.DOTALL)
    for match in bad_format_numbers_regex.finditer(tree_str):
        bad_format = match.group(1)
        good_format = "{:.10f}".format(float(bad_format), 10)
        tree_str = tree_str.replace(bad_format, good_format)
    return tree_str

def simulate_sequence_data(kappa, omega0, omega1, omega2, omega0_weight, omega1_weight, selection_intensity_parameter,
                           true_history_path, output_dir, num_of_replicates, aln_len):
    # prepare directory for simulation output
    sequence_output_dir = output_dir  # _kappa_" + str(kappa) + "_omega0_" + str(omega0) + "_omega1_" + str(omega1) + "_omega2_" + str(omega2) + "_theta1_" + str(omega0_weight) + "_theta2_" + str(omega1_weight/(1-omega0_weight)) + "/"
    if not os.path.exists(sequence_output_dir):
        res = os.system("mkdir -p " + sequence_output_dir)
    control_file_path = output_dir + "control.txt"
    bpp_labels_path = output_dir + "labels.bpp"

    # read the character history and derive from it a tree and a labeling in INDELible control file compatible format
    true_history = Tree(true_history_path, format=1)
    label_regex = re.compile("\{(.*?)\}", re.DOTALL)
    node_to_branch_length_expression = dict()  # will contain the node name is well if the node is internal, as internal node names should be omitted from the tree string
    node_to_label = dict()
    node_to_branch_index = dict()  # will help in creation of branch labeling in relax parameters file
    i = 0
    for node in true_history.traverse("postorder"):
        if not node.is_root():
            node_label = label_regex.search(node.name).group(1)
            node_name = node.name.replace("{" + node_label + "}", "")
            node.name = node_name  # update the name of the node to exclude the label
            if node_label == "0":
                node_to_label[node_name] = "BG"
            else:
                node_to_label[node_name] = "FG"
            node_to_branch_index[node_name] = i
            i += 1
            node_dist = node.dist
            if int(node.dist) == node.dist:
                node_dist = int(node.dist)
            node_to_branch_length_expression[node_name] = node_name + ":" + str(node_dist)
    tree_str = true_history.write(outfile=None,
                                  format=5)  # get a newick representation of the tree without internal nodes names
    # fix tree str number formatting
    tree_str = fix_tree_str_format(tree_str)
    tree_labels_str = true_history.write(outfile=None,
                                         format=1)  # get a newick representation of the tree with internal nodes names
    bpp_bg_labels_str = "model1.nodes_id = "
    bpp_fg_labels_str = "model2.nodes_id = "
    for node_name in node_to_label.keys():
        before = tree_labels_str
        node = true_history.search_nodes(name=node_name)[0]
        if node.is_leaf():
            tree_labels_str = tree_labels_str.replace(node_to_branch_length_expression[node_name],
                                                      node_name + " #" + node_to_label[node_name])
        else:
            tree_labels_str = tree_labels_str.replace(node_to_branch_length_expression[node_name],
                                                      " #" + node_to_label[node_name])
        after = tree_labels_str
        if (before == after):
            print("\nfailed to replace expression")
            print("node name: ", node_name)
            print("branch length expression: ", node_to_branch_length_expression[node_name])
            print("node label expression: ", " #" + node_to_label[node_name])
        if node_to_label[node_name] == "BG":
            bpp_bg_labels_str = bpp_bg_labels_str + str(node_to_branch_index[node_name]) + ","
        else:
            bpp_fg_labels_str = bpp_fg_labels_str + str(node_to_branch_index[node_name]) + ","
    # set the label of the root as the label of one of its immediate sons
    root = true_history.get_tree_root()
    root_son = root.get_children()[0]
    son_label = node_to_label[root_son.name]
    tree_labels_str = tree_labels_str.replace(";", "#" + son_label + ";")
    # write the labels into the labels path
    with open(bpp_labels_path, "w") as bpp_labels_file:
        bpp_labels_file.write(bpp_bg_labels_str[:-1] + "\n")
        bpp_labels_file.write(bpp_fg_labels_str[:-1])

    # create control file
    control_file_template = '''[TYPE] CODON 1                    

[SETTINGS]
    [ancestralprint]    FALSE
	[fastaextension]    fas
	[output]          	FASTA
	[fileperrep]      	TRUE
    [printrates]        TRUE

  /* Notice that the number of classes and proportions do not change below. */

                   // Kap p0  p1  w0  w1  w2     (p2=1-p1-p0=0.5) 
[MODEL] BG [submodel] <kappa> <omega0_weight> <omega1_weight> <bg_omega0> <bg_omega1> <bg_omega2>
[MODEL] FG [submodel] <kappa> <omega0_weight> <omega1_weight> <fg_omega0> <fg_omega1> <fg_omega2>

  /* 
     Like before, to get a correctly formatted [BRANCHES] block from a [TREE] block
     simply cut and paste the tree and change the branch lengths to model names.
     The stationary frequencies of the model at the root are used to generate the
     root sequence and this model defines the number of site categories used.
  */


[TREE]     t1  <tree_str>

[BRANCHES] b1  <tree_labels_str>


[PARTITIONS] Pname  [t1 b1 <aln_len>]    // tree t1, branchclass b1, root length 1000

[EVOLVE]     Pname  <num_of_replicates>  sequence_data  // 10 replicates generated from partition Pname'''
    control_file_content = control_file_template.replace("<kappa>", str(kappa))
    control_file_content = control_file_content.replace("<bg_omega0>", str(omega0))
    control_file_content = control_file_content.replace("<bg_omega1>", str(omega1))
    control_file_content = control_file_content.replace("<bg_omega2>", str(omega2))
    control_file_content = control_file_content.replace("<omega0_weight>", str(omega0_weight))
    control_file_content = control_file_content.replace("<omega1_weight>", str(omega1_weight))
    control_file_content = control_file_content.replace("<fg_omega0>", str(omega0 ** selection_intensity_parameter))
    control_file_content = control_file_content.replace("<fg_omega1>", str(omega1 ** selection_intensity_parameter))
    control_file_content = control_file_content.replace("<fg_omega2>", str(omega2 ** selection_intensity_parameter))
    control_file_content = control_file_content.replace("<tree_str>", tree_str)
    control_file_content = control_file_content.replace("<tree_labels_str>", tree_labels_str)
    control_file_content = control_file_content.replace("<num_of_replicates>", str(num_of_replicates))
    control_file_content = control_file_content.replace("<aln_len>", str(aln_len))
    with open(control_file_path, "w") as control_file:
        control_file.write(control_file_content)

    # execute INDELible
    res = os.system(
        "(echo " + control_file_path + " && cat) | /groups/itay_mayrose/halabikeren/indelible/INDELibleV1.03/src/indelible")

    return 0


if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='simulates alignments and character history under the TraitRELAX null and alternative models using Bio++ and INDELible')
    parser.add_argument('--output_dir', '-o', help='directory that holds the simulation output',
                        required=True)
    parser.add_argument('--true_history_path', '-t', help='path to the input history to simulate data over',
                        required=True)
    parser.add_argument('--kappa', '-kappa', help='nucleotide substitution rate parameter kappa', required=False,
                        default=2)
    parser.add_argument('--omega0', '-omega0', help='purifying selection omega to simulate with', required=False,
                        default=0.1)
    parser.add_argument('--omega1', '-omega1', help='neutral selection omega to simulate with', required=False,
                        default=1)
    parser.add_argument('--omega2', '-omega2', help='positive selection omega to simulate with', required=False,
                        default=2)
    parser.add_argument('--omega0_weight', '-p0', help='probability of purifying selection omega to simulate with',
                        required=False, default=0.5)
    parser.add_argument('--omega1_weight', '-p1', help='probability of neutral selection omega to simulate with',
                        required=False, default=0.4)
    parser.add_argument('--selection_intensity_parameter', '-k',
                        help='selection intensity parameter value under the alternative model', required=False,
                        default=0.5)
    parser.add_argument('--num_of_replicates', '-rep', help="number of codon MSAs to simulate", required=False,
                        default=100)
    parser.add_argument('--aln_len', '-al', help='number of codon positions to simulated', required=False, default=1000)

    args = parser.parse_args()
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        res = os.system(
            "mkdir -p " + output_dir)  # -p allows recusive mkdir in case one of the upper directories doesn't exist
    true_history_path = args.true_history_path
    if not os.path.exists(true_history_path):
        print("tree path " + true_history_path + " does not exist")
        exit(1)
    fix_tree_format(true_history_path)
    kappa = float(args.kappa)
    omega0 = float(args.omega0)
    omega1 = float(args.omega1)
    omega2 = float(args.omega2)
    omega0_weight = float(args.omega0_weight)
    omega1_weight = float(args.omega1_weight)
    selection_intensity_parameter = float(args.selection_intensity_parameter)
    num_of_replicates = int(args.num_of_replicates)
    aln_len = int(args.aln_len)


    # simulate sequence data
    simulate_sequence_data(kappa, omega0, omega1, omega2, omega0_weight, omega1_weight, selection_intensity_parameter,
                           true_history_path, output_dir, num_of_replicates, aln_len)
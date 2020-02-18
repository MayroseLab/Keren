import argparse, os, re
from ete3 import Tree
from time import sleep

def convert_history_to_simmap(tree_path, output_path):
    label_regex = re.compile("{(.*?)}")
    history = Tree(tree_path, format=1)
    #tree_str = history.write(format=8, outfile=None) # get a string of the tree newick only with names
    visited_nodes = [history.get_tree_root()]
    for node in history.traverse("postorder"):
        if not node in visited_nodes:
            # print("node name: ", node.name)
            node_label = label_regex.search(node.name).group(1)
            node_name = node.name.replace(label_regex.search(node.name).group(0),"")
            node_expression = ""
            if node.is_leaf():
                node_expression = node_expression + node_name
            node_expression = node_expression + ":{" + node_label + "," + str(node.dist)
            visited_nodes.append(node)
            if not node == history.get_tree_root():
                try:
                    curNode = node.up
                    while "mapping" in curNode.name:
                        node_label = label_regex.search(curNode.name).group(1)
                        node_expression = node_expression + ":" + node_label + "," + str(curNode.dist)
                        visited_nodes.append(curNode)
                        curNode = curNode.up
                except:
                    pass
            node_expression = node_expression + "}"
            node.name = node_expression
            #tree_str = tree_str.replace(node.name, node_expression)
    # remove mapping nodes
    for node in history.traverse("postorder"):
        if "mapping" in node.name:
            node.delete()
    # write tree on your own because ete3 writer has a bug
    tree_str = getTreeStr(history.get_tree_root())
    # print("tree_str: ", tree_str)
    with open(output_path, "w") as output_file:
        output_file.write(tree_str + ";")

# writes tree with nodes names only in newick format
def getTreeStr(root):
    if root.is_leaf():
        return root.name

    tree_str = "("
    for child in root.get_children():
        tree_str = tree_str + getTreeStr(child) + ","
    tree_str = tree_str[:-1] + ")" + root.name
    return tree_str

def fix_tree_format(tree_path):
    tree = Tree(tree_path,format=1)
    tree_str = tree.write(outfile=None,format=1)
    bad_format_numbers_regex = re.compile(":(\d*\.?\d*e-\d*)", re.MULTILINE | re.DOTALL)
    for match in bad_format_numbers_regex.finditer(tree_str):
        bad_format = match.group(1)
        good_format = "{:.10f}".format(float(bad_format),10)
        tree_str = tree_str.replace(bad_format, good_format)
    with open(tree_path, "w") as tree_file:
        tree_file.write(tree_str)

def remove_spaces(filepath):
    with open(filepath, "r") as input:
        file_content = input.read()
        file_content = file_content.replace(" ", "")
    with open(filepath, "w") as output:
        output.write(file_content)

def fix_tree_str_format(tree_str):
    bad_format_numbers_regex = re.compile(":(\d*\.?\d*e-\d*)", re.MULTILINE | re.DOTALL)
    for match in bad_format_numbers_regex.finditer(tree_str):
        bad_format = match.group(1)
        good_format = "{:.10f}".format(float(bad_format),10)
        tree_str = tree_str.replace(bad_format, good_format)
    return tree_str

# don't stop until you get a simulation with both types of states
def simulate_character_data(character_model_mu, character_model_pi0, tree_path, output_dir):

    return_res = False
    while not return_res:
        character_output_dir = output_dir + "character_data/"
        if not os.path.exists(character_output_dir):
            res = os.system("mkdir -p " + character_output_dir)
        true_history_path = character_output_dir + "true_history.nwk"
        character_data_path = character_output_dir + "character_data.fas"
        history_tree_path = character_output_dir + "history_tree.nwk"

        # create parameters file
        character_simulation_parameters_path = character_output_dir + "characterParams.bpp"
        with open(character_simulation_parameters_path, "w") as char_param_file:
            char_param_file.write("init.tree = user\ninput.tree.format = Newick\ninput.tree.file = " + tree_path + "\n")
            char_param_file.write("character_model.mu = " + str(character_model_mu) + "\n")
            char_param_file.write("character_model.pi0 = " + str(character_model_pi0) + "\n")
            char_param_file.write("character_history.seq_path = " + str(character_data_path) + "\n")
            char_param_file.write("character_history.tree_path = " + str(true_history_path) + "\n")

        # call to simulator with parameters file
        res = os.system("/groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite/traitsimulator param=" + character_simulation_parameters_path)
        fix_tree_format(true_history_path)

        # write the tree without labels into a file
        true_history = Tree(true_history_path,format=1)
        for node in true_history.traverse():
            node_name = re.sub("\{.*?\}","",node.name)
            node.name = node_name
        true_history.write(outfile=history_tree_path, format=5)
        true_history_str = true_history.write(outfile=None, format=5)
        # #insufficient condition - need the character data to consist of a mix of 0s and 1s
        # with open(true_history_path, "r") as true_history_file:
        #     true_history_str = true_history_file.read()
        #     if "{1}" in true_history_str and "{0}" in true_history_str:
        #         return_res = True
        true_history = Tree(true_history_path, format=1)
        availableStates = []
        for leaf in true_history:
            if "{0}" in leaf.name:
                availableStates.append(0)
            elif "{1}" in leaf.name:
                availableStates.append(1)
        # print("availableStates: ", availableStates) # debug
        if 0 in availableStates and 1 in availableStates:
            return_res = True

        # write true history in simmap format
        true_history_simmap_path = true_history_path.replace(".nwk", "_simmap.nwk")
        convert_history_to_simmap(true_history_path, true_history_simmap_path)

        # # plot mapping in R
        # true_history_visual_path = true_history_path.replace(".nwk", ".pdf")
        # res=os.system("Rscript --vanilla /groups/itay_mayrose/halabikeren/myScripts/R/plot_history.R " + true_history_simmap_path + " " + true_history_visual_path)

    return true_history_path, character_data_path, history_tree_path

def set_relax_param_file(output_path, sequence_data_path, tree_path, kappa, omega0, omega1, omega2, omega0_weight, omega1_weight, selection_intensity_parameter, nuc1_theta, nuc1_theta1, nuc1_theta2, nuc2_theta, nuc2_theta1, nuc2_theta2, nuc3_theta, nuc3_theta1, nuc3_theta2, labels):

    param_template = '''# Global variables:
verbose = 1

# ----------------------------------------------------------------------------------------
#                                     Input alignment file
# ----------------------------------------------------------------------------------------

alphabet=Codon(letter=DNA)
genetic_code=Standard
input.sequence.file = <sequence_data_path>
input.sequence.format = Fasta
input.sequence.sites_to_use = all
input.sequence.max_gap_allowed = 100%
input.sequence.remove_stop_codons = yes

# ----------------------------------------------------------------------------------------
#                                     Input tree file
# ----------------------------------------------------------------------------------------

init.tree = user
input.tree.file = <tree_path>
input.tree.format = Newick
init.brlen.method = Input

# ----------------------------------------------------------------------------------------
#                                    Model specification
# ----------------------------------------------------------------------------------------

model1 = RELAX(kappa=<kappa>,p=<p>,omega1=<omega1>,omega2=<omega2>,k=1,theta1=<theta1>,theta2=<theta2>,frequencies=F3X4, 1_Full.theta=<nuc1_theta>,1_Full.theta1=<nuc1_theta1>, 1_Full.theta2=<nuc1_theta2>, 2_Full.theta=<nuc2_theta>,2_Full.theta1=<nuc2_theta1>, 2_Full.theta2=<nuc2_theta2>,3_Full.theta=<nuc3_theta>,3_Full.theta1=<nuc3_theta1>, 3_Full.theta2=<nuc3_theta2>)
model2 = RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,frequencies=F3X4,1_Full.theta=RELAX.1_Full.theta_1,1_Full.theta1=RELAX.1_Full.theta1_1,1_Full.theta2=RELAX.1_Full.theta2_1,2_Full.theta=RELAX.2_Full.theta_1,2_Full.theta1=RELAX.2_Full.theta1_1,2_Full.theta2=RELAX.2_Full.theta2_1,3_Full.theta=RELAX.3_Full.theta_1,3_Full.theta1=RELAX.3_Full.theta1_1,3_Full.theta2=RELAX.3_Full.theta2_1,k=<selection_intensity_parameter>)
nonhomogeneous = general
nonhomogeneous.number_of_models = 2
nonhomogeneous.stationarity = yes
site.number_of_paths = 2
site.path1 = model1[YN98.omega_1]&model2[YN98.omega_1]
site.path2 = model1[YN98.omega_2]&model2[YN98.omega_2]
rate_distribution = Constant() //Gamma(n=4, alpha=0.358)
likelihood.recursion = simple
likelihood.recursion_simple.compression = recursive

# ----------------------------------------------------------------------------------------
#                                    optimization parameters
# ----------------------------------------------------------------------------------------

optimization.tolerance = 0.000001
optimization.max_number_f_eval = 10000
optimization = FullD(derivatives=Newton,nstep=10)
optimization.final = powell

# ----------------------------------------------------------------------------------------
#                                    branches partition
# ----------------------------------------------------------------------------------------

<labels>
'''

    # translate simulator input to bio++ parameters
    p = omega0 / omega1
    theta1 = omega0_weight
    theta2 = omega1_weight / (1-omega0_weight)

    param_content = param_template.replace("<sequence_data_path>", sequence_data_path)
    param_content = param_content.replace("<tree_path>", tree_path)
    param_content = param_content.replace("<kappa>", str(kappa))
    param_content = param_content.replace("<p>", str(p))
    param_content = param_content.replace("<omega1>", str(omega1))
    param_content = param_content.replace("<omega2>", str(omega2))
    param_content = param_content.replace("<theta1>", str(theta1))
    param_content = param_content.replace("<theta2>", str(theta2))
    param_content = param_content.replace("<nuc1_theta>", str(nuc1_theta))
    param_content = param_content.replace("<nuc1_theta1>", str(nuc1_theta1))
    param_content = param_content.replace("<nuc1_theta2>", str(nuc1_theta2))
    param_content = param_content.replace("<nuc2_theta>", str(nuc2_theta))
    param_content = param_content.replace("<nuc2_theta1>", str(nuc2_theta1))
    param_content = param_content.replace("<nuc2_theta2>", str(nuc2_theta2))
    param_content = param_content.replace("<nuc3_theta>", str(nuc3_theta))
    param_content = param_content.replace("<nuc3_theta1>", str(nuc3_theta1))
    param_content = param_content.replace("<nuc3_theta2>", str(nuc3_theta2))
    param_content = param_content.replace("<selection_intensity_parameter>", str(selection_intensity_parameter))
    param_content = param_content.replace("<labels>", labels)

    with open(output_path, "w") as output_file:
        output_file.write(param_content)

def set_traitrelax_param_file(output_dir, output_path, sequence_data_path, tree_path, character_data_path, character_model_mu, character_model_pi0, kappa, omega0, omega1, omega2, omega0_weight, omega1_weight, selection_intensity_parameter, history_tree_path, nuc1_theta, nuc1_theta1, nuc1_theta2, nuc2_theta, nuc2_theta1, nuc2_theta2, nuc3_theta, nuc3_theta1, nuc3_theta2, labels_str):
    param_template = '''# Global variables:
verbose = 1

# ----------------------------------------------------------------------------------------
#                                     Input character file
# ----------------------------------------------------------------------------------------

input.character.file = <character_data_path>

# ----------------------------------------------------------------------------------------
#                                     Input alignment file
# ----------------------------------------------------------------------------------------

alphabet=Codon(letter=DNA)
genetic_code=Standard
input.sequence.file = <sequence_data_path>
input.sequence.format = Fasta
input.sequence.sites_to_use = all
input.sequence.max_gap_allowed = 100%
input.sequence.remove_stop_codons = yes

# ----------------------------------------------------------------------------------------
#                                     Input tree file
# ----------------------------------------------------------------------------------------

init.tree = user
input.tree.file = <tree_path>
input.tree.format = Newick
init.brlen.method = Input

# ----------------------------------------------------------------------------------------
#                                     Character Model specification
# ----------------------------------------------------------------------------------------

character_model.set_initial_parameters = true
character_model.mu = <mu>
character_model.pi0 = <pi0>

# ----------------------------------------------------------------------------------------
#                                     Sequence Model specification
# ----------------------------------------------------------------------------------------

sequence_model.set_initial_parameters = true
model1 = RELAX(kappa=<kappa>,p=<p>,omega1=<omega1>,omega2=<omega2>,k=1,theta1=<theta1>,theta2=<theta2>,frequencies=F3X4, 1_Full.theta=<nuc1_theta>,1_Full.theta1=<nuc1_theta1>, 1_Full.theta2=<nuc1_theta2>, 2_Full.theta=<nuc2_theta>,2_Full.theta1=<nuc2_theta1>, 2_Full.theta2=<nuc2_theta2>,3_Full.theta=<nuc3_theta>,3_Full.theta1=<nuc3_theta1>, 3_Full.theta2=<nuc3_theta2>)
model2 = RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,frequencies=F3X4,1_Full.theta=RELAX.1_Full.theta_1,1_Full.theta1=RELAX.1_Full.theta1_1,1_Full.theta2=RELAX.1_Full.theta2_1,2_Full.theta=RELAX.2_Full.theta_1,2_Full.theta1=RELAX.2_Full.theta1_1,2_Full.theta2=RELAX.2_Full.theta2_1,3_Full.theta=RELAX.3_Full.theta_1,3_Full.theta1=RELAX.3_Full.theta1_1,3_Full.theta2=RELAX.3_Full.theta2_1,k=<selection_intensity_parameter>)

# ----------------------------------------------------------------------------------------
#                                    optimization parameters
# ----------------------------------------------------------------------------------------

optimization.tolerance = 0.000001
optimization.max_number_f_eval = 10000
optimization = FullD(derivatives=Newton,nstep=10)
optimization.final = powell
#optimization.message_handler = std
optimization.ignore_parameters = RELAX.k_1,BrLen

# ----------------------------------------------------------------------------------------
#                                    output files data                                    
# ----------------------------------------------------------------------------------------

optimization.profiler = <log_path>
output.tree.file = <expected_history_path>
output.debug.dir = <debug_dir>
output.tree.format = Newick
true_history.tree.file = <history_tree_path>
'''

    # translate simulator input to bio++ parameters
    if not os.path.exists(output_dir):
        res = os.system("mkdir -p " + output_dir)
    expected_history_path = output_dir + "expected_history.nwk"
    debug_dir = output_dir + "histories_evaluation/"
    if not os.path.exists(debug_dir):
        res = os.system("mkdir -p " + debug_dir)
    log_path = output_dir + "traitrelax_optimization.log"
    p = omega0 / omega1
    theta1 = omega0_weight
    theta2 = omega1_weight / (1-omega0_weight)

    param_content = param_template.replace("<character_data_path>", character_data_path)
    param_content = param_content.replace("<sequence_data_path>", sequence_data_path)
    param_content = param_content.replace("<tree_path>", tree_path)
    param_content = param_content.replace("<mu>", str(character_model_mu))
    param_content = param_content.replace("<pi0>", str(character_model_pi0))
    param_content = param_content.replace("<kappa>", str(kappa))
    param_content = param_content.replace("<p>", str(p))
    param_content = param_content.replace("<omega1>", str(omega1))
    param_content = param_content.replace("<omega2>", str(omega2))
    param_content = param_content.replace("<theta1>", str(theta1))
    param_content = param_content.replace("<theta2>", str(theta2))
    param_content = param_content.replace("<nuc1_theta>", str(nuc1_theta))
    param_content = param_content.replace("<nuc1_theta1>", str(nuc1_theta1))
    param_content = param_content.replace("<nuc1_theta2>", str(nuc1_theta2))
    param_content = param_content.replace("<nuc2_theta>", str(nuc2_theta))
    param_content = param_content.replace("<nuc2_theta1>", str(nuc2_theta1))
    param_content = param_content.replace("<nuc2_theta2>", str(nuc2_theta2))
    param_content = param_content.replace("<nuc3_theta>", str(nuc3_theta))
    param_content = param_content.replace("<nuc3_theta1>", str(nuc3_theta1))
    param_content = param_content.replace("<nuc3_theta2>", str(nuc3_theta2))
    param_content = param_content.replace("<selection_intensity_parameter>", str(selection_intensity_parameter))
    param_content = param_content.replace("<log_path>", log_path)
    param_content = param_content.replace("<expected_history_path>", expected_history_path)
    param_content = param_content.replace("<debug_dir>", debug_dir)
    param_content = param_content.replace("<history_tree_path>", history_tree_path)
    adjusted_labels_str = labels_str.replace("model1.nodes_id", "true_history.model1.nodes_id")
    adjusted_labels_str = adjusted_labels_str.replace("model2.nodes_id", "true_history.model2.nodes_id")
    param_content = param_content + adjusted_labels_str

    with open(output_path, "w") as output_file:
        output_file.write(param_content)

def get_labels(true_history_path):

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
    labels_str = bpp_bg_labels_str[:-1] + "\n" + bpp_fg_labels_str[:-1]

    return labels_str


if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(description='simulates character histories Bio++')
    parser.add_argument('--output_dir', '-o', help='directory that holds the simulation output',required=True)
    parser.add_argument('--tree_path', '-t', help='path to the hold the simulated input trees',required=True)
    parser.add_argument('--character_model_mu', '-mu', help='character substitution rate to set as initial parameter', required=False, default=10.)
    parser.add_argument('--character_model_pi0', '-pi0', help='character state 0 frequency ot simulate with', required=False, default=0.5)
    parser.add_argument('--initial_character_model_mu', '-imu', help='character substitution rate to set as initial parameter', required=False, default=10.)
    parser.add_argument('--initial_character_model_pi0', '-ipi0', help='character state 0 frequency ot simulate with', required=False, default=0.5)
    parser.add_argument('--initial_kappa', '-ikappa', help='nucleotide substitution rate parameter kappa', required=False,
                        default=1)
    parser.add_argument('--initial_omega0', '-iomega0', help='purifying selection omega to set as initial parameter', required=False,
                        default=0.1)
    parser.add_argument('--initial_omega1', '-iomega1', help='neutral selection omega to set as initial parameter', required=False,
                        default=1)
    parser.add_argument('--initial_omega2', '-iomega2', help='positive selection omega to set as initial parameter', required=False,
                        default=2)
    parser.add_argument('--initial_omega0_weight', '-ip0', help='probability of purifying selection omega to set as initial parameter',
                        required=False, default=0.5)
    parser.add_argument('--initial_omega1_weight', '-ip1', help='probability of neutral selection omega to set as initial parameter',
                        required=False, default=0.4)
    parser.add_argument('--initial_selection_intensity_parameter', '-ik',
                        help='selection intensity parameter value under the alternative model', required=False,
                        default=0.5)
    parser.add_argument('--num_of_replicates', '-rep', help="number of codon MSAs to simulate", required=False, default=50)
    parser.add_argument('--sequence_data_path', '-s', help="path ot the sequence data file", required=True)

    parser.add_argument('--initial_nuc1_theta', '-in1t', help='Bio++ parameter 1_Full.theta from which simulation values of codon frequencies will be derived', required=False, default=0.5)
    parser.add_argument('--initial_nuc1_theta1', '-in1t1', help='Bio++ parameter 1_Full.theta1 from which simulation values of codon frequencies will be derived',required=False, default=0.5)
    parser.add_argument('--initial_nuc1_theta2', '-in1t2', help='Bio++ parameter 1_Full.theta from which simulation values of codon frequencies will be derived', required=False, default=0.5)
    parser.add_argument('--initial_nuc2_theta', '-in2t', help='Bio++ parameter 2_Full.theta from which simulation values of codon frequencies will be derived', required=False, default=0.5)
    parser.add_argument('--initial_nuc2_theta1', '-in2t1', help='Bio++ parameter 2_Full.theta1 from which simulation values of codon frequencies will be derived',required=False, default=0.5)
    parser.add_argument('--initial_nuc2_theta2', '-in2t2', help='Bio++ parameter 2_Full.theta from which simulation values of codon frequencies will be derived', required=False, default=0.5)
    parser.add_argument('--initial_nuc3_theta', '-in3t', help='Bio++ parameter 3_Full.theta from which simulation values of codon frequencies will be derived', required=False, default=0.5)
    parser.add_argument('--initial_nuc3_theta1', '-in3t1', help='Bio++ parameter 3_Full.theta1 from which simulation values of codon frequencies will be derived',required=False, default=0.5)
    parser.add_argument('--initial_nuc3_theta2', '-in3t2', help='Bio++ parameter 3_Full.theta from which simulation values of codon frequencies will be derived', required=False, default=0.5)

    args = parser.parse_args()
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        res = os.system("mkdir -p " + output_dir) # -p allows recusive mkdir in case one of the upper directories doesn't exist
    relax_param_dir = output_dir + "relax_param/"
    if not os.path.exists(relax_param_dir):
        res = os.system("mkdir -p " + relax_param_dir)
    traitrelax_param_dir = output_dir + "traitrelax_param/"
    if not os.path.exists(traitrelax_param_dir):
        res = os.system("mkdir -p " + traitrelax_param_dir)
    num_of_replicates = int(args.num_of_replicates)
    tree_path = args.tree_path
    character_model_mu = float(args.character_model_mu)
    character_model_pi0 = float(args.character_model_pi0)
    initial_kappa = float(args.initial_kappa)
    initial_omega0 = float(args.initial_omega0)
    initial_omega1 = float(args.initial_omega1)
    initial_omega2 = float(args.initial_omega2)
    initial_omega0_weight = float(args.initial_omega0_weight)
    initial_omega1_weight = float(args.initial_omega1_weight)
    initial_selection_intensity_parameter = float(args.initial_selection_intensity_parameter)
    sequence_data_path = args.sequence_data_path

    initial_nuc1_theta = float(args.initial_nuc1_theta)
    initial_nuc1_theta1 = float(args.initial_nuc1_theta1)
    initial_nuc1_theta2 = float(args.initial_nuc1_theta2)
    initial_nuc2_theta = float(args.initial_nuc2_theta)
    initial_nuc2_theta1 = float(args.initial_nuc2_theta1)
    initial_nuc2_theta2 = float(args.initial_nuc2_theta2)
    initial_nuc3_theta = float(args.initial_nuc3_theta)
    initial_nuc3_theta1 = float(args.initial_nuc3_theta1)
    initial_nuc3_theta2 = float(args.initial_nuc3_theta2)

    for rep in range(num_of_replicates):
        fix_tree_format(tree_path)
        print("**** simulating replicate " + str(rep) + " ****")
        # set simulation output directory
        simulation_output_dir = output_dir + "replicate_" + str(rep) + "/"
        if not os.path.exists(simulation_output_dir):
            res = os.system("mkdir -p " + simulation_output_dir)
        # simulate character data
        true_history_path, character_data_path, history_tree_path = simulate_character_data(character_model_mu, character_model_pi0, tree_path, simulation_output_dir)
        # simulate sequence data
        labels_str = get_labels(true_history_path)
        # set the parameters file for RELAX
        set_relax_param_file(relax_param_dir + str(rep) + ".bpp", sequence_data_path, history_tree_path, initial_kappa, initial_omega0, initial_omega1, initial_omega2, initial_omega0_weight, initial_omega1_weight, initial_selection_intensity_parameter, initial_nuc1_theta, initial_nuc1_theta1, initial_nuc1_theta2, initial_nuc2_theta, initial_nuc2_theta1, initial_nuc2_theta2, initial_nuc3_theta, initial_nuc3_theta1, initial_nuc3_theta2, labels_str)
        # set parameters file for TraitRELAX
        set_traitrelax_param_file(simulation_output_dir+"traitrelax_result/", traitrelax_param_dir + str(rep) + ".bpp", sequence_data_path, tree_path, character_data_path, character_model_mu, character_model_pi0, initial_kappa, initial_omega0, initial_omega1, initial_omega2, initial_omega0_weight, initial_omega1_weight, initial_selection_intensity_parameter, history_tree_path, initial_nuc1_theta, initial_nuc1_theta1, initial_nuc1_theta2, initial_nuc2_theta, initial_nuc2_theta1, initial_nuc2_theta2, initial_nuc3_theta, initial_nuc3_theta1, initial_nuc3_theta2, labels_str)

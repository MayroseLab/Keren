import argparse, os, re
from ete3 import Tree
from time import sleep


def convert_history_to_simmap(tree_path, output_path):
    label_regex = re.compile("{(.*?)}")
    history = Tree(tree_path, format=1)
    # tree_str = history.write(format=8, outfile=None) # get a string of the tree newick only with names
    visited_nodes = [history.get_tree_root()]
    for node in history.traverse("postorder"):
        if not node in visited_nodes:
            node_label = label_regex.search(node.name).group(1)
            node_name = node.name.replace(label_regex.search(node.name).group(0), "")
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
            # tree_str = tree_str.replace(node.name, node_expression)
    # remove mapping nodes
    for node in history.traverse("postorder"):
        if "mapping" in node.name:
            node.delete()
    # write tree on your own because ete3 writer has a bug
    tree_str = getTreeStr(history.get_tree_root())
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


def fix_tree_format(tree_path, base_tree=False):
    format=1
    if base_tree:
        format=5
    tree = Tree(tree_path, format=1)
    tree_str = tree.write(outfile=None, format=format)
    bad_format_numbers_regex = re.compile(":(\d*\.?\d*e-\d*)", re.MULTILINE | re.DOTALL)
    for match in bad_format_numbers_regex.finditer(tree_str):
        bad_format = match.group(1)
        good_format = "{:.10f}".format(float(bad_format), 10)
        tree_str = tree_str.replace(bad_format, good_format)
    with open(tree_path, "w") as tree_file:
        tree_file.write(tree_str)


def scale_tree(input_tree_path, output_tree_path, scaling_factor=1.0):
    tree = Tree(input_tree_path, format=1)
    for node in tree.traverse():
        node.dist = node.dist * scaling_factor
    tree.write(outfile=output_tree_path, format=1)


def remove_spaces(filepath):
    with open(filepath, "r") as input:
        file_content = input.read()
        file_content = file_content.replace(" ", "")
    with open(filepath, "w") as output:
        output.write(file_content)

# writes parameters file to be given to the simulator
def write_simulation_parameters(tree_path, model_parameters, sequence_data_path, character_data_path,
                                          history_tree_path, true_history_path, aln_len, simulation_parameters_path):
    parameters_str_template = '''# ----------------------------------------------------------------------------------------
#                                     Input tree file
# ----------------------------------------------------------------------------------------

init.tree = user
input.tree.file = <tree_path>
input.tree.format = Newick
init.brlen.method = Input

# ----------------------------------------------------------------------------------------
#                             Character simulation parameters 
# ----------------------------------------------------------------------------------------

character_model.mu = <mu>
character_model.pi0 = <pi0>

# ----------------------------------------------------------------------------------------
#                             Sequence simulation parameters 
# ----------------------------------------------------------------------------------------

sequence.num_of_sites = <aln_len>
sequence_model.kappa = <kappa>
sequence_model.omega0 = <omega0>
sequence_model.p0 = <omega0_weight>
sequence_model.omega1 = <omega1>
sequence_model.p1 = <omega1_weight>
sequence_model.omega2 = <omega2>
sequence_model.k = <selection_intensity_parameter>
sequence_model.frequencies = F3X4
sequence_model.frequencies.codon_pos1 = <nuc1_theta>, <nuc1_theta1>, <nuc1_theta2>
sequence_model.frequencies.codon_pos2 = <nuc2_theta>, <nuc2_theta1>, <nuc2_theta2>
sequence_model.frequencies.codon_pos3 = <nuc3_theta>, <nuc3_theta1>, <nuc3_theta2>

# ----------------------------------------------------------------------------------------
#                                Output settings ans specs
# ----------------------------------------------------------------------------------------

character.data_path = <character_data_path>
character.unlabeled_history_path = <history_tree_path>
character.labeled_history_path = <true_history_path>
sequence.data_path = <sequence_data_path>
output.sequence.format = Fasta
output.tree.format = Newick
'''
    parameters_str = parameters_str_template.replace("<aln_len>", str(aln_len))
    parameters_str = parameters_str.replace("<character_data_path>", character_data_path)
    parameters_str = parameters_str.replace("<sequence_data_path>", sequence_data_path)
    parameters_str = parameters_str.replace("<history_tree_path>", history_tree_path)
    parameters_str = parameters_str.replace("<true_history_path>", true_history_path)
    parameters_str = parameters_str.replace("<tree_path>", tree_path)

    for parameter_name in model_parameters:
        parameters_str = parameters_str.replace("<"+parameter_name+">", str(model_parameters[parameter_name]))
    with open(simulation_parameters_path, "w") as outfile:
        outfile.write(parameters_str)

def fix_tree_str_format(tree_str):
    bad_format_numbers_regex = re.compile(":(\d*\.?\d*e-\d*)", re.MULTILINE | re.DOTALL)
    for match in bad_format_numbers_regex.finditer(tree_str):
        bad_format = match.group(1)
        good_format = "{:.10f}".format(float(bad_format), 10)
        tree_str = tree_str.replace(bad_format, good_format)
    return tree_str


def set_relax_param_file(output_path, sequence_data_path, tree_path, kappa, omega0, omega1, omega2, omega0_weight,
                         omega1_weight, selection_intensity_parameter, nuc1_theta, nuc1_theta1, nuc1_theta2, nuc2_theta,
                         nuc2_theta1, nuc2_theta2, nuc3_theta, nuc3_theta1, nuc3_theta2, labels):
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

model1 = RELAX(kappa=<kappa>,p=<p>,omega1=<omega1>,omega2=<omega2>,k=1,theta1=<theta1>,theta2=<theta2>,frequencies=F3X4,1_Full.theta=<nuc1_theta>,1_Full.theta1=<nuc1_theta1>, 1_Full.theta2=<nuc1_theta2>, 2_Full.theta=<nuc2_theta>,2_Full.theta1=<nuc2_theta1>, 2_Full.theta2=<nuc2_theta2>,3_Full.theta=<nuc3_theta>,3_Full.theta1=<nuc3_theta1>, 3_Full.theta2=<nuc3_theta2>)
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
    theta2 = omega1_weight / (1 - omega0_weight)

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

def set_traitrelax_param_file(output_dir, output_path, sequence_data_path, tree_path, character_data_path, kappa,
                              omega0, omega1, omega2, omega0_weight,
                              omega1_weight, selection_intensity_parameter, history_tree_path, nuc1_theta, nuc1_theta1,
                              nuc1_theta2, nuc2_theta, nuc2_theta1, nuc2_theta2, nuc3_theta, nuc3_theta1, nuc3_theta2,
                              labels_str, character_model_mu=None, character_model_pi0=None):
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
    theta2 = omega1_weight / (1 - omega0_weight)

    param_content = param_template.replace("<character_data_path>", character_data_path)
    param_content = param_content.replace("<sequence_data_path>", sequence_data_path)
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
    param_content = param_content.replace("<log_path>", log_path)
    param_content = param_content.replace("<expected_history_path>", expected_history_path)
    param_content = param_content.replace("<debug_dir>", debug_dir)
    param_content = param_content.replace("<history_tree_path>", history_tree_path)
    adjusted_labels_str = labels_str.replace("model1.nodes_id", "true_history.model1.nodes_id")
    adjusted_labels_str = adjusted_labels_str.replace("model2.nodes_id", "true_history.model2.nodes_id")
    param_content = param_content + adjusted_labels_str

    if character_model_mu != None and character_model_pi0 != None:
        param_content = param_content + '''# ----------------------------------------------------------------------------------------
#                                     Character Model specification
# ----------------------------------------------------------------------------------------

character_model.set_initial_parameters = true
character_model.mu = <mu>
character_model.pi0 = <pi0>'''

        param_content = param_content.replace("<mu>", str(character_model_mu))
        param_content = param_content.replace("<pi0>", str(character_model_pi0))

    with open(output_path, "w") as output_file:
        output_file.write(param_content)

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='simulates alignments and character history under the TraitRELAX null and alternative models using Bio++ and INDELible')
    parser.add_argument('--simulator_path', '-sp', help='full path to the simulator program', required=False, default="/groups/itay_mayrose/halabikeren/programs/Bio++/TraitRELAX/TraitRELAX/simulator")
    parser.add_argument('--output_dir', '-o', help='directory that holds the simulation output',
                        required=True)
    # parser.add_argument('--tree_path', '-t', help='path to the input tree to simulate data over',
    #                     required=True)
    parser.add_argument('--tree_dir', '-t', help='path to the hold the simulated input trees',
                        required=True)
    parser.add_argument('--character_model_mu', '-mu', help='character substitution rate to simulate with',
                        required=False, default=10.)
    parser.add_argument('--character_model_pi0', '-pi0', help='character state 0 frequency ot simulate with',
                        required=False, default=0.5)
    parser.add_argument('--kappa', '-kappa', help='nucleotide substitution rate parameter kappa', required=False,
                        default=2)
    parser.add_argument('--omega0', '-omega0', help='purifying selection omega to simulate with', required=False,
                        default=0.1)
    parser.add_argument('--omega1', '-omega1', help='neutral selection omega to simulate with', required=False,
                        default=0.8)
    parser.add_argument('--omega2', '-omega2', help='positive selection omega to simulate with', required=False,
                        default=2)
    parser.add_argument('--omega0_weight', '-p0', help='probability of purifying selection omega to simulate with',
                        required=False, default=0.5)
    parser.add_argument('--omega1_weight', '-p1', help='probability of neutral selection omega to simulate with',
                        required=False, default=0.4)
    parser.add_argument('--selection_intensity_parameter', '-k',
                        help='selection intensity parameter value under the alternative model', required=False,
                        default=1)
    parser.add_argument('--num_of_replicates', '-rep', help="number of codon MSAs to simulate", required=False,
                        default=30)
    parser.add_argument('--aln_len', '-al', help='number of codon positions to simulated', required=False, default=1000)

    parser.add_argument('--scaling_factor', '-sc',
                        help='scaling factor of the tree with which simulations are conducted and inference is executed',
                        required=False, default=1)

    parser.add_argument('--initial_mu', '-imu',
                        help='character substitution rate to initialize traitrelax parameters file with',
                        required=False,
                        default=None)
    parser.add_argument('--initial_pi0', '-ipi0',
                        help='character state 0 frequency to initialize traitrelax parameters file with',
                        required=False,
                        default=None)
    parser.add_argument('--initial_kappa', '-ikappa',
                        help='nucleotide substitution rate parameter kappa to initialize parameters file with',
                        required=False,
                        default=2)
    parser.add_argument('--initial_omega0', '-iomega0',
                        help='purifying selection omega to initialize parameters file with', required=False,
                        default=0.1)
    parser.add_argument('--initial_omega1', '-iomega1',
                        help='neutral selection omega to initialize parameters file with', required=False,
                        default=0.8)
    parser.add_argument('--initial_omega2', '-iomega2',
                        help='positive selection omega to initialize parameters file with', required=False,
                        default=2)
    parser.add_argument('--initial_omega0_weight', '-ip0',
                        help='probability of purifying selection omega to initialize parameters file with',
                        required=False, default=0.5)
    parser.add_argument('--initial_omega1_weight', '-ip1',
                        help='probability of neutral selection omega to initialize parameters file with',
                        required=False, default=0.4)
    parser.add_argument('--initial_selection_intensity_parameter', '-ik',
                        help='selection intensity parameter value under the alternative model to initialize parameters file with',
                        required=False,
                        default=1)

    parser.add_argument('--character_data_path', '-cd',
                        help='path to the character data file to be given as input to traitrelax', required=False,
                        default="")
    parser.add_argument('--nuc1_theta', '-n1t',
                        help='Bio++ parameter 1_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--nuc1_theta1', '-n1t1',
                        help='Bio++ parameter 1_Full.theta1 from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--nuc1_theta2', '-n1t2',
                        help='Bio++ parameter 1_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--nuc2_theta', '-n2t',
                        help='Bio++ parameter 2_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--nuc2_theta1', '-n2t1',
                        help='Bio++ parameter 2_Full.theta1 from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--nuc2_theta2', '-n2t2',
                        help='Bio++ parameter 2_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--nuc3_theta', '-n3t',
                        help='Bio++ parameter 3_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--nuc3_theta1', '-n3t1',
                        help='Bio++ parameter 3_Full.theta1 from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--nuc3_theta2', '-n3t2',
                        help='Bio++ parameter 3_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)

    parser.add_argument('--initial_nuc1_theta', '-in1t',
                        help='Bio++ parameter 1_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--initial_nuc1_theta1', '-in1t1',
                        help='Bio++ parameter 1_Full.theta1 from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--initial_nuc1_theta2', '-in1t2',
                        help='Bio++ parameter 1_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--initial_nuc2_theta', '-in2t',
                        help='Bio++ parameter 2_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--initial_nuc2_theta1', '-in2t1',
                        help='Bio++ parameter 2_Full.theta1 from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--initial_nuc2_theta2', '-in2t2',
                        help='Bio++ parameter 2_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--initial_nuc3_theta', '-in3t',
                        help='Bio++ parameter 3_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--initial_nuc3_theta1', '-in3t1',
                        help='Bio++ parameter 3_Full.theta1 from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)
    parser.add_argument('--initial_nuc3_theta2', '-in3t2',
                        help='Bio++ parameter 3_Full.theta from which simulation values of codon frequencies will be derived',
                        required=False, default=0.5)

    args = parser.parse_args()
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        res = os.system(
            "mkdir -p " + output_dir)  # -p allows recusive mkdir in case one of the upper directories doesn't exist
    relax_param_dir = output_dir + "relax_param/"
    if not os.path.exists(relax_param_dir):
        res = os.system("mkdir -p " + relax_param_dir)
    traitrelax_param_dir = output_dir + "traitrelax_param/"
    if not os.path.exists(traitrelax_param_dir):
        res = os.system("mkdir -p " + traitrelax_param_dir)

    simulator_path = args.simulator_path
    trees_dir = args.tree_dir
    trees_paths = [path for path in os.listdir(trees_dir) if not ".nfs" in path]

    # character model parameters
    character_model_mu = float(args.character_model_mu)
    character_model_pi0 = float(args.character_model_pi0)

    # sequence model parameters
    kappa = float(args.kappa)
    omega0 = float(args.omega0)
    omega1 = float(args.omega1)
    omega2 = float(args.omega2)
    omega0_weight = float(args.omega0_weight)
    omega1_weight = float(args.omega1_weight)
    nuc1_theta = float(args.nuc1_theta)
    nuc1_theta1 = float(args.nuc1_theta1)
    nuc1_theta2 = float(args.nuc1_theta2)
    nuc2_theta = float(args.nuc2_theta)
    nuc2_theta1 = float(args.nuc2_theta1)
    nuc2_theta2 = float(args.nuc2_theta2)
    nuc3_theta = float(args.nuc3_theta)
    nuc3_theta1 = float(args.nuc3_theta1)
    nuc3_theta2 = float(args.nuc3_theta2)
    selection_intensity_parameter = float(args.selection_intensity_parameter)

    # general simulation parameters
    scaling_factor = float(args.scaling_factor)
    aln_len = int(args.aln_len)
    num_of_replicates = int(args.num_of_replicates)

    # parmeters to initalize model with via parameters file
    initial_kappa = float(args.initial_kappa)
    initial_omega0 = float(args.initial_omega0)
    initial_omega1 = float(args.initial_omega1)
    initial_omega2 = float(args.initial_omega2)
    initial_omega0_weight = float(args.initial_omega0_weight)
    initial_omega1_weight = float(args.initial_omega1_weight)
    initial_selection_intensity_parameter = float(args.initial_selection_intensity_parameter)
    initial_character_model_mu = args.initial_mu
    if initial_character_model_mu != None:
        initial_character_model_mu = float(initial_character_model_mu)
    initial_character_model_pi0 = args.initial_pi0
    if initial_character_model_pi0 != None:
        initial_character_model_pi0 = float(initial_character_model_pi0)
    orig_character_data_path = args.character_data_path
    use_scaled_as_hist = False
    if orig_character_data_path != "":
        use_scaled_as_hist = True
    initial_nuc1_theta = float(args.initial_nuc1_theta)
    initial_nuc1_theta1 = float(args.initial_nuc1_theta1)
    initial_nuc1_theta2 = float(args.initial_nuc1_theta2)
    initial_nuc2_theta = float(args.initial_nuc2_theta)
    initial_nuc2_theta1 = float(args.initial_nuc2_theta1)
    initial_nuc2_theta2 = float(args.initial_nuc2_theta2)
    initial_nuc3_theta = float(args.initial_nuc3_theta)
    initial_nuc3_theta1 = float(args.initial_nuc3_theta1)
    initial_nuc3_theta2 = float(args.initial_nuc3_theta2)

    rep_regex = re.compile("(\d*).nwk", re.MULTILINE | re.DOTALL)
    for tree_filepath in os.listdir(trees_dir):
        tree_path = trees_dir + tree_filepath
        if not ".nwk" in tree_path:
            continue
        rep = int(rep_regex.search(tree_filepath).group(1))
        if int(rep) >= num_of_replicates:
            continue
        fix_tree_format(tree_path, base_tree=True)
        if scaling_factor == 1:
            scaled_tree_path = tree_path
        else:
            scaled_tree_path = tree_path.replace(".nwk", "_scaled_by_" + str(scaling_factor) + ".nwk")
            scale_tree(tree_path, scaled_tree_path, scaling_factor=scaling_factor)
        input_tree_path = ""
        if use_scaled_as_hist:  # if character data was given
            input_tree_path = scaled_tree_path
        print("**** simulating replicate " + str(rep) + " ****")

        # simulate data with Bio++
        replicate_output_dir = output_dir + "replicate_" + str(rep) + "/"
        if not os.path.exists(replicate_output_dir):
            res = os.system("mkdir -p " + replicate_output_dir)
        simulation_parameters_path = replicate_output_dir + "simulation_params.bpp"
        model_parameters = {"mu": character_model_mu,
                           "pi0": character_model_pi0,
                            "kappa": kappa,
                            "omega0": omega0,
                            "omega1": omega1,
                            "omega2": omega2,
                            "omega0_weight": omega0_weight,
                            "omega1_weight": omega1_weight,
                            "nuc1_theta": nuc1_theta,
                            "nuc1_theta1": nuc1_theta1,
                            "nuc1_theta2": nuc1_theta2,
                            "nuc2_theta": nuc2_theta,
                            "nuc2_theta1": nuc2_theta1,
                            "nuc2_theta2": nuc2_theta2,
                            "nuc3_theta": nuc3_theta,
                            "nuc3_theta1": nuc3_theta1,
                            "nuc3_theta2": nuc3_theta2,
                            "selection_intensity_parameter": selection_intensity_parameter}
        sequence_data_path = replicate_output_dir + "sequence_data.fas"
        character_data_path = replicate_output_dir + "character_data.fas"
        history_tree_path = replicate_output_dir + "unlabeled_trait_history.nwk"
        true_history_path = replicate_output_dir + "labeled_trait_history.nwk"
        write_simulation_parameters(tree_path, model_parameters, sequence_data_path, character_data_path, history_tree_path, true_history_path, aln_len, simulation_parameters_path)
        simulation_output_log = replicate_output_dir + "simulator_log.txt"
        res = os.system(simulator_path + " param=" + simulation_parameters_path + " > " + simulation_output_log)

        # extract the labeling of nodes in the trait history for relax parameters and traitrelax debugging
        labels_str_regex = re.compile("\#* assignment of nodes to model in the simulated trait history is: \#*\n\n(.*?)\n\n\#*", re.MULTILINE | re.DOTALL)
        with open(simulation_output_log, "r") as infile:
            content = infile.read()
        labels_str = labels_str_regex.search(content).group(1)

        # set the parameters file for RELAX
        set_relax_param_file(relax_param_dir + str(rep) + ".bpp", sequence_data_path, scaled_tree_path,
                             initial_kappa, initial_omega0, initial_omega1, initial_omega2, initial_omega0_weight,
                             initial_omega1_weight, initial_selection_intensity_parameter, initial_nuc1_theta,
                             initial_nuc1_theta1, initial_nuc1_theta2, initial_nuc2_theta, initial_nuc2_theta1,
                             initial_nuc2_theta2, initial_nuc3_theta, initial_nuc3_theta1, initial_nuc3_theta2,
                             labels_str)

        # set parameters file for TraitRELAX
        set_traitrelax_param_file(replicate_output_dir + "traitrelax_result/",
                                  traitrelax_param_dir + str(rep) + ".bpp", sequence_data_path, scaled_tree_path,
                                  character_data_path, initial_kappa, initial_omega0, initial_omega1,
                                  initial_omega2, initial_omega0_weight, initial_omega1_weight,
                                  initial_selection_intensity_parameter, history_tree_path,
                                  initial_nuc1_theta, initial_nuc1_theta1, initial_nuc1_theta2, initial_nuc2_theta,
                                  initial_nuc2_theta1, initial_nuc2_theta2, initial_nuc3_theta, initial_nuc3_theta1,
                                  initial_nuc3_theta2, labels_str, initial_character_model_mu,
                                  initial_character_model_pi0)

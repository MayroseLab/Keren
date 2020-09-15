import os, re, json, argparse
from ete3 import Tree

def write_traitrelax_parameters(json_data, char_data_path, msa_path, tree_path, debug_dir, output_path):
    template = '''# Global variables:
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
model1 = RELAX(kappa=<kappa>,p=<p>,omega1=<omega1>,omega2=<omega2>,k=1,theta1=<theta1>,theta2=<theta2>,frequencies=F3X4)
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

output.debug.dir = <debug_dir>
output.tree.format = Newick
    '''

    omega0 = json_data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["0"]["omega"]
    omega1 = json_data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["1"]["omega"]
    omega2 = json_data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["2"]["omega"]
    p0 = json_data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["0"]["proportion"]
    p1 = json_data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["1"]["proportion"]
    k = json_data["test results"]["relaxation or intensification parameter"]

    param = template.replace("<character_data_path>", char_data_path)
    param = param.replace("<sequence_data_path>", msa_path)
    param = param.replace("<tree_path>", tree_path)
    param = param.replace("<kappa>", "2")
    param = param.replace("<p>", str(omega0/ omega1))
    param = param.replace("<omega1>", str(omega1))
    param = param.replace("<omega2>", str(omega2))
    param = param.replace("<theta1>", str(p0))
    param = param.replace("<theta2>", str((1-p1)/p0))
    param = param.replace("<selection_intensity_parameter>", str(k))
    param = param.replace("<debug_dir>", debug_dir)

    with open(output_path, "w") as outfile:
        outfile.write(param)


def write_relax_parameters(json_data, msa_path, tree_path,  node_id_to_label, output_path):
    template = '''# Global variables:
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

model1 = RELAX(kappa=<kappa>,p=<p>,omega1=<omega1>,omega2=<omega2>,k=1,theta1=<theta1>,theta2=<theta2>,frequencies=F3X4)
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

    omega0 = json_data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["0"]["omega"]
    omega1 = json_data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["1"]["omega"]
    omega2 = json_data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["2"]["omega"]
    p0 = json_data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["0"]["proportion"]
    p1 = json_data["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]["1"]["proportion"]
    k = json_data["test results"]["relaxation or intensification parameter"]

    param = template.replace("<sequence_data_path>", msa_path)
    param = param.replace("<tree_path>", tree_path)
    param = param.replace("<kappa>", "2")
    param = param.replace("<p>", str(omega0/ omega1))
    param = param.replace("<omega1>", str(omega1))
    param = param.replace("<omega2>", str(omega2))
    param = param.replace("<theta1>", str(p0))
    param = param.replace("<theta2>", str((1-p1)/p0))
    param = param.replace("<selection_intensity_parameter>", str(k))

    state_0_nodes = [str(id) for id in node_id_to_label.keys() if node_id_to_label[id] == 1]
    state_1_nodes = [str(id) for id in node_id_to_label.keys() if node_id_to_label[id] == 2]
    labels = "model1.nodes_id = " + ','.join(state_0_nodes) + "\nmodel2.nodes_id = " + ','.join(state_1_nodes)
    param = param.replace("<labels>", labels)

    with open(output_path, "w") as outfile:
        outfile.write(param)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='generates based on relax json output file input files for Bio++ execution')
    parser.add_argument('--json_dir', '-i', help='directory holding the json output files', required=True)
    parser.add_argument('--MSAs_dir', '-m', help='directory of the input MSAs given to hyphy', required=True)
    parser.add_argument('--output_dir', '-o', help='directory to hold the input for Bio++', required=True)
    args = parser.parse_args()
    json_dir = args.json_dir
    MSAs_dir = args.MSAs_dir
    output_dir = args.output_dir
    trees_dir = output_dir + "/trees/" # directory to hold the trees with the branch lengths inferred by hyphy
    os.system("mkdir -p " + trees_dir)
    output_MSAs_dir = output_dir + "/MSAs/"
    os.system("mkdir -p " + output_MSAs_dir) # directory to hold the input MSAs for Bio++
    character_data_dir = output_dir + "/charData/"
    os.system("mkdir -p " + character_data_dir)  # directory to hold the input trait data for Bio++
    traitrelax_params_dir = output_dir + "/traitrelax_parameters/"
    os.system("mkdir -p " + traitrelax_params_dir)
    relax_params_dir = output_dir + "/relax_parameters/"
    os.system("mkdir -p " + relax_params_dir)

    for path in os.listdir(json_dir):
        if "json" in path:
            with open(json_dir + "/" + path, 'r') as j:
                json_data = json.load(j)

        # write the MSA
        dataset_name_regex = re.compile(".*/([^/]*)")
        dataset_name = dataset_name_regex.search(json_data["input"]["file name"]).group(1)
        input_MSA_path = MSAs_dir + "/" + dataset_name
        output_MSA_path = output_MSAs_dir + dataset_name + ".fasta"
        with open(input_MSA_path, "r") as infile:
            content = infile.read()
        with open(input_MSA_path, "r") as infile:
            tree = infile.readlines()[-1]
        msa = content.replace(tree, "")
        with open(output_MSA_path, "w") as outfile:
            outfile.write(msa)

        # write the tree
        tree_str_with_internals = json_data["input"]["trees"]["0"]
        tree_with_internals = Tree(tree_str_with_internals + ";", format=8)

        nodes = [node for node in tree_with_internals.traverse("postorder")]

        label_regex = re.compile("{(.*)}")
        node_id_to_label = dict()
        for i in range(len(nodes)):
            node = nodes[i]
            if node.name != "":
                label = json_data["tested"]["0"][node.name]
                if label == "Unclassified" or label == "Reference":
                    node_id_to_label[i] = 1
                else:
                    node_id_to_label[i] = 2
                bl = json_data["branch attributes"]["0"][node.name]["RELAX alternative"]
                node.dist = float(bl)
        tree_path = trees_dir + dataset_name + ".nwk"
        tree_with_internals.write(outfile=tree_path, format=5)

        # write the character data
        character_data_path = character_data_dir + dataset_name + ".fas"
        with open(character_data_path, "w") as outfile:
            for leaf in tree_with_internals.get_leaves():
                label = json_data["tested"]["0"][leaf.name]
                if label == "Unclassified" or label == "Reference":
                    label = "0"
                else:
                    label = "1"
                outfile.write(">" + leaf.name + "\n" + label + "\n")


        # write the parameters file
        debug_dir = output_dir + "/traitrelax_result/" + dataset_name
        res = os.system("mkdir -p " + debug_dir)
        write_traitrelax_parameters(json_data, character_data_path, output_MSA_path, tree_path, debug_dir, traitrelax_params_dir + dataset_name + ".bpp")
        write_relax_parameters(json_data, output_MSA_path, tree_path, node_id_to_label, relax_params_dir + dataset_name + ".bpp")

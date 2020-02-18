import re, os
from ete3 import Tree
from utils.alterAlignments import convert_formats


def extract_tree(input_path, biopp_ibput_aln_path, biopp_base_tree_output_path, biopp_relax_param_output_path, biopp_char_data_output_path, biopp_traitrelax_param_output_path,  hyphy_base_tree_output_path, hyphy_labeled_tree_output_path):
    tree_regex = re.compile("begin trees;\n\s*tree .*?\s*=\s*[^\(]*(.*?;)\nend;")
    label_regex = re.compile("\{(.*?)\}")

    # process tree from input file
    with open(input_path, "r") as input_file:
        content = input_file.read()

    # write the labeled tree as hyphy input
    labeled_tree_str = tree_regex.search(content).group(1)
    labeled_tree_str = re.sub("\[&label=\{\"\wG\"\},!color=\#.*?\]", "", labeled_tree_str) # clear color labels
    with open(hyphy_labeled_tree_output_path, "w") as hyphy_labeled_tree_output_file:
        hyphy_labeled_tree_output_file.write(labeled_tree_str)

    # save the labels and while at it, clean the tree from labeles and write the base tree as both hyphy an bio++ input
    labeled_tree = Tree(hyphy_labeled_tree_output_path, format=1)
    biopp_lables_dict = {0: [], 1: []}
    node_name_to_node_label = dict()
    biopp_leaf_to_label = {}
    node_id = 0
    # give names to internal nodes in post-order
    ni = 1

    for node in labeled_tree.traverse("postorder"):
        if node != labeled_tree.get_tree_root():
            if node.name == "{R}" or node.name == "{T}" or node.name == "":
                node_name = "Node" + str(ni) + node.name
                node.name = node_name
            ni += 1

    for node in labeled_tree.traverse("postorder"):
        if node != labeled_tree.get_tree_root():
            try:
                text_label = label_regex.search(node.name).group(1)
            except:
                print("node: ", node.name, " is missing a label in the supplementary data")
                # add label based on the majority rule on the labels of the children
                support_values = {0: 0, 1: 0}
                for child in node.get_children():
                    child_label = node_name_to_node_label[child.name]
                    support_values[child_label] += 1
                if support_values[0] >= support_values[1]:
                    text_label = "R"
                else:
                    text_label = "T"
                node.name = node.name + "{" + text_label + "}"

            numeric_label = 0
            if text_label == "T":
                numeric_label = 1

            biopp_lables_dict[numeric_label].append(node_id)

            if node.is_leaf():
                biopp_leaf_to_label[node.name] = numeric_label

            node_name_to_node_label[node.name] = numeric_label
            node_id += 1

    # remove internal names
    for node in labeled_tree.traverse():
        if not node.is_leaf() and node != labeled_tree.get_tree_root():
            text_label = label_regex.search(node.name).group(1)
            node.name = "{" + text_label + "}"
    # write the labeled tree to hyphy labeled tree output path
    labeled_tree.write(outfile=hyphy_labeled_tree_output_path, format=1)

    # now remove the lables from the nodes names and then write the base tree
    for node in labeled_tree.traverse("postorder"):
        node_name = re.sub("\{.*?\}", "", node.name)
        node.name = node_name

    labeled_tree.write(outfile=biopp_base_tree_output_path, format=5)
    labeled_tree.write(outfile=hyphy_base_tree_output_path, format=5)

    # write traitRELAX character data
    char_str = ""
    for leaf_name in biopp_leaf_to_label.keys():
        char_str = char_str + ">" + leaf_name + "\n" + str(biopp_leaf_to_label[leaf_name]) + "\n"
    with open(biopp_char_data_output_path, "w") as char_data_file:
        char_data_file.write(char_str)

    # write bio++ parameters file for RELAX and TraitRELAX
    relax_param_template = '''# Global variables:
verbose = 1

# ----------------------------------------------------------------------------------------
#                                     Input alignment file
# ----------------------------------------------------------------------------------------

alphabet=Codon(letter=DNA)
genetic_code=Standard
input.sequence.file = <SEQ_DATA_PATH>
input.sequence.format = Fasta
input.sequence.sites_to_use = all
input.sequence.max_gap_allowed = 100%
input.sequence.remove_stop_codons = yes

# ----------------------------------------------------------------------------------------
#                                     Input tree file
# ----------------------------------------------------------------------------------------

init.tree = user
input.tree.file = <TREE_PATH>
input.tree.format = Newick
init.brlen.method = Input

# ----------------------------------------------------------------------------------------
#                                    Model specification
# ----------------------------------------------------------------------------------------

model1 = RELAX(kappa=2.0,p=0.5,omega1=1.0,omega2=1.0,k=1,theta1=0.5,theta2=0.8,frequencies=F0)
model2 = RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,frequencies=F0,k=1.0)
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

<LABELS_STR>
'''

    # generate labels string
    labels_str = "model1.nodes_id = "
    for node_id in biopp_lables_dict[0]:
        labels_str = labels_str + str(node_id) + ","
    labels_str = labels_str[:len(labels_str)-1] + "\nmodel2.nodes_id = "
    for node_id in biopp_lables_dict[1]:
        labels_str = labels_str + str(node_id) + ","
    labels_str = labels_str[:len(labels_str) - 1]

    server_biopp_ibput_aln_path = biopp_ibput_aln_path.replace("C:/Users/ItayMNB7/Google Drive/biopp_data/", "/groups/itay_mayrose/halabikeren/TraitRELAX/")
    relax_param_content = relax_param_template.replace("<SEQ_DATA_PATH>", server_biopp_ibput_aln_path)
    server_biopp_base_tree_output_path = biopp_base_tree_output_path.replace("C:/Users/ItayMNB7/Google Drive/biopp_data/", "/groups/itay_mayrose/halabikeren/TraitRELAX/")
    relax_param_content = relax_param_content.replace("<TREE_PATH>", server_biopp_base_tree_output_path)
    relax_param_content = relax_param_content.replace("<LABELS_STR>", labels_str)

    with open(biopp_relax_param_output_path, "w") as biopp_relax_param_input_file:
        biopp_relax_param_input_file.write(relax_param_content)

    # .replace("<CHAR_DATA_PATH>", biopp_char_data_output_path)

    traitrelax_param_template = '''# Global variables:
verbose = 1

# ----------------------------------------------------------------------------------------
#                                     Input character file
# ----------------------------------------------------------------------------------------

input.character.file = <CHAR_DATA_PATH>

# ----------------------------------------------------------------------------------------
#                                     Input alignment file
# ----------------------------------------------------------------------------------------

alphabet=Codon(letter=DNA)
genetic_code=Standard
input.sequence.file = <SEQ_DATA_PATH>
input.sequence.format = Fasta
input.sequence.sites_to_use = all
input.sequence.max_gap_allowed = 100%
input.sequence.remove_stop_codons = yes

# ----------------------------------------------------------------------------------------
#                                     Input tree file
# ----------------------------------------------------------------------------------------

init.tree = user
input.tree.file = <TREE_PATH>
input.tree.format = Newick
init.brlen.method = Input

# ----------------------------------------------------------------------------------------
#                                     Character Model specification
# ----------------------------------------------------------------------------------------

character_model.set_initial_parameters = true
character_model.mu = 2.0
character_model.pi0 = 0.5

# ----------------------------------------------------------------------------------------
#                                     Sequence Model specification
# ----------------------------------------------------------------------------------------

sequence_model.set_initial_parameters = true
model1 = RELAX(kappa=2.0,p=0.5,omega1=1.0,omega2=1.0,k=1,theta1=0.5,theta2=0.8,frequencies=F0)
model2 = RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,frequencies=F0,k=1.0)

# ----------------------------------------------------------------------------------------
#                                    optimization parameters
# ----------------------------------------------------------------------------------------

optimization.tolerance = 0.000001
optimization.max_number_f_eval = 10000
optimization = FullD(derivatives=Newton,nstep=10)
optimization.final = powell
#optimization.message_handler = std
optimization.ignore_parameters = RELAX.k_1,BrLen
'''

    traitrelax_param_content = traitrelax_param_template.replace("<SEQ_DATA_PATH>", server_biopp_ibput_aln_path)
    traitrelax_param_content = traitrelax_param_content.replace("<TREE_PATH>", server_biopp_base_tree_output_path)
    server_biopp_char_data_output_path = biopp_char_data_output_path.replace("C:/Users/ItayMNB7/Google Drive/biopp_data/", "/groups/itay_mayrose/halabikeren/TraitRELAX/")
    traitrelax_param_content = traitrelax_param_content.replace("<CHAR_DATA_PATH>", server_biopp_char_data_output_path)

    with open(biopp_traitrelax_param_output_path, "w") as biopp_traitrelax_param_output_file:
        biopp_traitrelax_param_output_file.write(traitrelax_param_content)

    return 0


if __name__ == '__main__':

    # process input
    input_dir = "C:/Users/ItayMNB7/Downloads/SupportingInformation/"
    output_dir = "C:/Users/ItayMNB7/Google Drive/real_data/RELAX_real_data/"

    # for each file in the input directory, set a folder with the relevant name
    input_paths = [input_dir + input_path for input_path in os.listdir(input_dir)]
    input_name_regex = re.compile("Fig(.*?)\.nex")
    for input_path in input_paths:
        if ".nex" in input_path:
            input_name = input_name_regex.search(input_path).group(1)
            print("input name: ", input_name)
            spec_output_dir = output_dir + input_name + "/"
            if not os.path.exists(spec_output_dir):
                os.mkdir(spec_output_dir)
            spec_hyphy_output_dir = spec_output_dir + "hyphy/"
            if not os.path.exists(spec_hyphy_output_dir):
                os.mkdir(spec_hyphy_output_dir)
            spec_biopp_output_dir = spec_output_dir + "biopp/"
            if not os.path.exists(spec_biopp_output_dir):
                os.mkdir(spec_biopp_output_dir)

            # process alignment
            convert_formats(input_path, spec_hyphy_output_dir+"aln.fas", "nexus", "fasta")
            convert_formats(input_path, spec_biopp_output_dir + "aln.fas", "nexus", "fasta")
            with open(spec_biopp_output_dir + "aln.fas", "r") as outfile:
                content = outfile.read()
            content = content.replace("\r\n", "\n")
            with open(spec_biopp_output_dir + "aln.fas", "w") as outfile:
                outfile.write(content)

            # process tree
            extract_tree(input_path, spec_biopp_output_dir + "aln.fas", spec_biopp_output_dir+"tree.nwk", spec_biopp_output_dir+"relax_param.bpp", spec_biopp_output_dir+"char_data.fas", spec_biopp_output_dir+"traitrelax_param.bpp", spec_hyphy_output_dir+"tree.nwk", spec_hyphy_output_dir+"labeled_tree.nwk")

    print("extraction of real data is complete")
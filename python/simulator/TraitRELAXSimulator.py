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
    print("tree_str: ", tree_str)
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
    tree = Tree(tree_path, format=1)
    tree_str = tree.write(outfile=None, format=1)
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


def compute_codon_frequencies(nuc1_theta, nuc1_theta1, nuc1_theta2, nuc2_theta, nuc2_theta1, nuc2_theta2, nuc3_theta,
                              nuc3_theta1, nuc3_theta2):
    codon_to_frequency = dict()
    nucleotides = ["A", "C", "G", "T"]

    # / *codon frequencies parameterization using F3X4:
    # for each _Full.theta, corresponding to a a codon position over {0, 1, 2}:
    #     getFreq_(0) = theta1 * (1. - theta);
    # getFreq_(1) = (1 - theta2) * theta;
    # getFreq_(2) = theta2 * theta;
    # getFreq_(3) = (1 - theta1) * (1. - theta); * /

    pos_to_frequencies = []
    pos_to_frequencies.append(dict())
    pos_to_frequencies[0]["A"] = nuc1_theta1 * (1 - nuc1_theta)
    pos_to_frequencies[0]["C"] = nuc1_theta * (1 - nuc1_theta2)
    pos_to_frequencies[0]["G"] = nuc1_theta2 * nuc1_theta
    pos_to_frequencies[0]["T"] = (1 - nuc1_theta1) * (1 - nuc1_theta)

    pos_to_frequencies.append(dict())
    pos_to_frequencies[1]["A"] = nuc2_theta1 * (1 - nuc2_theta)
    pos_to_frequencies[1]["C"] = nuc2_theta * (1 - nuc2_theta2)
    pos_to_frequencies[1]["G"] = nuc2_theta2 * nuc2_theta
    pos_to_frequencies[1]["T"] = (1 - nuc2_theta1) * (1 - nuc2_theta)

    pos_to_frequencies.append(dict())
    pos_to_frequencies[2]["A"] = nuc3_theta * (1 - nuc3_theta)
    pos_to_frequencies[2]["C"] = nuc3_theta * (1 - nuc3_theta2)
    pos_to_frequencies[2]["G"] = nuc3_theta2 * nuc3_theta
    pos_to_frequencies[2]["T"] = (1 - nuc3_theta1) * (1 - nuc3_theta)

    # compute frequencies of stop codons and normalize the other frequencies by the left sum
    stop_codons = ["TAG", "TAA", "TGA"]
    stop_codons_frequencies = 0
    for codon in stop_codons:
        stop_codons_frequencies += pos_to_frequencies[0][codon[0]] * pos_to_frequencies[1][codon[1]] * \
                                   pos_to_frequencies[2][codon[2]]

    for i in nucleotides:
        freq1 = pos_to_frequencies[0][i]
        for j in nucleotides:
            freq2 = pos_to_frequencies[0][j]
            for k in nucleotides:
                freq3 = pos_to_frequencies[0][k]
                if i + j + k in stop_codons:
                    codon_to_frequency[i + j + k] = 0
                else:
                    codon_to_frequency[i + j + k] = freq1 * freq2 * freq3 / (1 - stop_codons_frequencies)

    return codon_to_frequency

def fix_tree_str_format(tree_str):
    bad_format_numbers_regex = re.compile(":(\d*\.?\d*e-\d*)", re.MULTILINE | re.DOTALL)
    for match in bad_format_numbers_regex.finditer(tree_str):
        bad_format = match.group(1)
        good_format = "{:.10f}".format(float(bad_format), 10)
        tree_str = tree_str.replace(bad_format, good_format)
    return tree_str

# don't stop until you get a simulation with both types of dats
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
        res = os.system(
            "/groups/itay_mayrose/halabikeren/biopp_test/bppsuite/build/bppSuite/traitsimulator param=" + character_simulation_parameters_path)
        fix_tree_format(true_history_path)

        # write the tree without labels into a file
        true_history = Tree(true_history_path, format=1)
        for node in true_history.traverse():
            node_name = re.sub("\{.*?\}", "", node.name)
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
        print("availableStates: ", availableStates)  # debug
        if 0 in availableStates and 1 in availableStates:
            return_res = True

        # write true history in simmap format
        true_history_simmap_path = true_history_path.replace(".nwk", "_simmap.nwk")
        convert_history_to_simmap(true_history_path, true_history_simmap_path)

        # # plot mapping in R
        # true_history_visual_path = true_history_path.replace(".nwk", ".pdf")
        # res=os.system("Rscript --vanilla /groups/itay_mayrose/halabikeren/myScripts/R/plot_history.R " + true_history_simmap_path + " " + true_history_visual_path)

    return true_history_path, character_data_path, history_tree_path

# indelible's branch-site simulation feature does not allow selective regime to swtich upon transitions between bramch classes, just like in our implementation
# one can force selective regimes transitions along the tree using a pre-given labeling, as can be seem at the bottom of the documentation in:
# http://abacus.gene.ucl.ac.uk/software/indelible/tutorial/branchsite.shtml
def simulate_sequence_data(kappa, omega0, omega1, omega2, omega0_weight, omega1_weight, selection_intensity_parameter,
                           true_history_path, output_dir, num_of_replicates, aln_len, nuc1_theta, nuc1_theta1,
                           nuc1_theta2, nuc2_theta, nuc2_theta1, nuc2_theta2, nuc3_theta, nuc3_theta1, nuc3_theta2):
    # prepare directory for simulation output
    sequence_output_dir = output_dir + "sequence_data/"  # _kappa_" + str(kappa) + "_omega0_" + str(omega0) + "_omega1_" + str(omega1) + "_omega2_" + str(omega2) + "_theta1_" + str(omega0_weight) + "_theta2_" + str(omega1_weight/(1-omega0_weight)) + "/"
    if not os.path.exists(sequence_output_dir):
        res = os.system("mkdir -p " + sequence_output_dir)
    control_file_path = sequence_output_dir + "control.txt"

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
    labels_str = bpp_bg_labels_str[:-1] + "\n" + bpp_fg_labels_str[:-1]

    # compute codon frequencies
    codon_to_frequency = compute_codon_frequencies(nuc1_theta, nuc1_theta1, nuc1_theta2, nuc2_theta, nuc2_theta1,
                                                   nuc2_theta2, nuc3_theta, nuc3_theta1, nuc3_theta2)

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
    [statefreq]  
          <TTT_freq> <TTC_freq> <TTA_freq> <TTG_freq>    //  TTT  TTC  TTA  TTG
          <TCT_freq> <TCC_freq> <TCA_freq> <TCG_freq>    //  TCT  TCC  TCA  TCG 
          <TAT_freq> <TAC_freq> <TAA_freq> <TAG_freq>    //  TAT  TAC  TAA  TAG
          <TGT_freq> <TGC_freq> <TGA_freq> <TGG_freq>    //  TGT  TGC  TGA  TGG

          <CTT_freq> <CTC_freq> <CTA_freq> <CTG_freq>    //  CTT  CTC  CTA  CTG 
          <CCT_freq> <CCC_freq> <CCA_freq> <CCG_freq>    //  CCT  CCC  CCA  CCG 
          <CAT_freq> <CAC_freq> <CAA_freq> <CAG_freq>    //  CAT  CAC  CAA  CAG 
          <CGT_freq> <CGC_freq> <CGA_freq> <CGG_freq>    //  CGT  CGC  CGA  CGG 

          <ATT_freq> <ATC_freq> <ATA_freq> <ATG_freq>    //  ATT  ATC  ATA  ATG  
          <ACT_freq> <ACC_freq> <ACA_freq> <ACG_freq>    //  ACT  ACC  ACA  ACG  
          <AAT_freq> <AAC_freq> <AAA_freq> <AAG_freq>    //  AAT  AAC  AAA  AAG  
          <AGT_freq> <AGC_freq> <AGA_freq> <AGG_freq>    //  AGT  AGC  AGA  AGG 

          <GTT_freq> <GTC_freq> <GTA_freq> <GTG_freq>    //  GTT  GTC  GTA  GTG 
          <GCT_freq> <GCC_freq> <GCA_freq> <GCG_freq>    //  GCT  GCC  GCA  GCG  
          <GAT_freq> <GAC_freq> <GAA_freq> <GAG_freq>    //  GAT  GAC  GAA  GAG  
          <GGT_freq> <GGC_freq> <GGA_freq> <GGG_freq>    //  GGT  GGC  GGA  GGG
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
    control_file_content = control_file_template.replace("<kappa>", str("%.15f" % kappa))
    control_file_content = control_file_content.replace("<bg_omega0>", str("%.15f" % omega0))
    control_file_content = control_file_content.replace("<bg_omega1>", str("%.15f" % omega1))
    control_file_content = control_file_content.replace("<bg_omega2>", str("%.15f" % omega2))
    control_file_content = control_file_content.replace("<omega0_weight>", str("%.15f" % omega0_weight))
    control_file_content = control_file_content.replace("<omega1_weight>", str("%.15f" % omega1_weight))
    control_file_content = control_file_content.replace("<fg_omega0>",
                                                        str("%.15f" % omega0 ** selection_intensity_parameter))
    control_file_content = control_file_content.replace("<fg_omega1>",
                                                        str("%.15f" % omega1 ** selection_intensity_parameter))
    control_file_content = control_file_content.replace("<fg_omega2>",
                                                        str("%.15f" % omega2 ** selection_intensity_parameter))
    control_file_content = control_file_content.replace("<tree_str>", tree_str)
    control_file_content = control_file_content.replace("<tree_labels_str>", tree_labels_str)
    control_file_content = control_file_content.replace("<num_of_replicates>", str(num_of_replicates))
    control_file_content = control_file_content.replace("<aln_len>", str(aln_len))
    control_file_content = control_file_content.replace("<TTT_freq>", str(codon_to_frequency["TTT"]))
    control_file_content = control_file_content.replace("<TTC_freq>", str(codon_to_frequency["TTC"]))
    control_file_content = control_file_content.replace("<TTA_freq>", str(codon_to_frequency["TTA"]))
    control_file_content = control_file_content.replace("<TTG_freq>", str(codon_to_frequency["TTG"]))
    control_file_content = control_file_content.replace("<TCT_freq>", str(codon_to_frequency["TCT"]))
    control_file_content = control_file_content.replace("<TCC_freq>", str(codon_to_frequency["TCC"]))
    control_file_content = control_file_content.replace("<TCA_freq>", str(codon_to_frequency["TCA"]))
    control_file_content = control_file_content.replace("<TCG_freq>", str(codon_to_frequency["TCG"]))
    control_file_content = control_file_content.replace("<TAT_freq>", str(codon_to_frequency["TAT"]))
    control_file_content = control_file_content.replace("<TAC_freq>", str(codon_to_frequency["TAC"]))
    control_file_content = control_file_content.replace("<TAA_freq>", str(codon_to_frequency["TAA"]))
    control_file_content = control_file_content.replace("<TAG_freq>", str(codon_to_frequency["TAG"]))
    control_file_content = control_file_content.replace("<TGT_freq>", str(codon_to_frequency["TGT"]))
    control_file_content = control_file_content.replace("<TGC_freq>", str(codon_to_frequency["TGC"]))
    control_file_content = control_file_content.replace("<TGA_freq>", str(codon_to_frequency["TGA"]))
    control_file_content = control_file_content.replace("<TGG_freq>", str(codon_to_frequency["TGG"]))
    control_file_content = control_file_content.replace("<CTT_freq>", str(codon_to_frequency["CTT"]))
    control_file_content = control_file_content.replace("<CTC_freq>", str(codon_to_frequency["CTC"]))
    control_file_content = control_file_content.replace("<CTA_freq>", str(codon_to_frequency["CTA"]))
    control_file_content = control_file_content.replace("<CTG_freq>", str(codon_to_frequency["CTG"]))
    control_file_content = control_file_content.replace("<CCT_freq>", str(codon_to_frequency["CCT"]))
    control_file_content = control_file_content.replace("<CCC_freq>", str(codon_to_frequency["CCC"]))
    control_file_content = control_file_content.replace("<CCA_freq>", str(codon_to_frequency["CCA"]))
    control_file_content = control_file_content.replace("<CCG_freq>", str(codon_to_frequency["CCG"]))
    control_file_content = control_file_content.replace("<CAT_freq>", str(codon_to_frequency["CAT"]))
    control_file_content = control_file_content.replace("<CAC_freq>", str(codon_to_frequency["CAC"]))
    control_file_content = control_file_content.replace("<CAA_freq>", str(codon_to_frequency["CAA"]))
    control_file_content = control_file_content.replace("<CAG_freq>", str(codon_to_frequency["CAG"]))
    control_file_content = control_file_content.replace("<CGT_freq>", str(codon_to_frequency["CGT"]))
    control_file_content = control_file_content.replace("<CGC_freq>", str(codon_to_frequency["CGC"]))
    control_file_content = control_file_content.replace("<CGA_freq>", str(codon_to_frequency["CGA"]))
    control_file_content = control_file_content.replace("<CGG_freq>", str(codon_to_frequency["CGG"]))
    control_file_content = control_file_content.replace("<ATT_freq>", str(codon_to_frequency["ATT"]))
    control_file_content = control_file_content.replace("<ATC_freq>", str(codon_to_frequency["ATC"]))
    control_file_content = control_file_content.replace("<ATA_freq>", str(codon_to_frequency["ATA"]))
    control_file_content = control_file_content.replace("<ATG_freq>", str(codon_to_frequency["ATG"]))
    control_file_content = control_file_content.replace("<ACT_freq>", str(codon_to_frequency["ACT"]))
    control_file_content = control_file_content.replace("<ACC_freq>", str(codon_to_frequency["ACC"]))
    control_file_content = control_file_content.replace("<ACA_freq>", str(codon_to_frequency["ACA"]))
    control_file_content = control_file_content.replace("<ACG_freq>", str(codon_to_frequency["ACG"]))
    control_file_content = control_file_content.replace("<AAT_freq>", str(codon_to_frequency["AAT"]))
    control_file_content = control_file_content.replace("<AAC_freq>", str(codon_to_frequency["AAC"]))
    control_file_content = control_file_content.replace("<AAA_freq>", str(codon_to_frequency["AAA"]))
    control_file_content = control_file_content.replace("<AAG_freq>", str(codon_to_frequency["AAG"]))
    control_file_content = control_file_content.replace("<AGT_freq>", str(codon_to_frequency["AGT"]))
    control_file_content = control_file_content.replace("<AGC_freq>", str(codon_to_frequency["AGC"]))
    control_file_content = control_file_content.replace("<AGA_freq>", str(codon_to_frequency["AGA"]))
    control_file_content = control_file_content.replace("<AGG_freq>", str(codon_to_frequency["AGG"]))
    control_file_content = control_file_content.replace("<GTT_freq>", str(codon_to_frequency["GTT"]))
    control_file_content = control_file_content.replace("<GTC_freq>", str(codon_to_frequency["GTC"]))
    control_file_content = control_file_content.replace("<GTA_freq>", str(codon_to_frequency["GTA"]))
    control_file_content = control_file_content.replace("<GTG_freq>", str(codon_to_frequency["GTG"]))
    control_file_content = control_file_content.replace("<GCT_freq>", str(codon_to_frequency["GCT"]))
    control_file_content = control_file_content.replace("<GCC_freq>", str(codon_to_frequency["GCC"]))
    control_file_content = control_file_content.replace("<GCA_freq>", str(codon_to_frequency["GCA"]))
    control_file_content = control_file_content.replace("<GCG_freq>", str(codon_to_frequency["GCG"]))
    control_file_content = control_file_content.replace("<GAT_freq>", str(codon_to_frequency["GAT"]))
    control_file_content = control_file_content.replace("<GAC_freq>", str(codon_to_frequency["GAC"]))
    control_file_content = control_file_content.replace("<GAA_freq>", str(codon_to_frequency["GAA"]))
    control_file_content = control_file_content.replace("<GAG_freq>", str(codon_to_frequency["GAG"]))
    control_file_content = control_file_content.replace("<GGT_freq>", str(codon_to_frequency["GGT"]))
    control_file_content = control_file_content.replace("<GGC_freq>", str(codon_to_frequency["GGC"]))
    control_file_content = control_file_content.replace("<GGA_freq>", str(codon_to_frequency["GGA"]))
    control_file_content = control_file_content.replace("<GGG_freq>", str(codon_to_frequency["GGG"]))
    with open(control_file_path, "w") as control_file:
        control_file.write(control_file_content)

    # execute INDELible
    res = os.chdir(sequence_output_dir)
    res = os.chdir(sequence_output_dir)
    res = os.system("/groups/itay_mayrose/halabikeren/programs/indelible/INDELibleV1.03/src/indelible")

    # check if the simulation is done, and sleep until done
    while not os.path.exists(sequence_output_dir + "sequence_data_1.fas"):
        sleep(3)

    # remove spaces from file names and rename files
    res = os.system("rm -r " + sequence_output_dir + "LOG.txt")
    os.chdir(sequence_output_dir)
    res = os.system("rm -r " + sequence_output_dir + "sequence_data_TRUE_1.fas")
    remove_spaces(sequence_output_dir + "sequence_data_1.fas")

    return sequence_output_dir + "sequence_data_1.fas", labels_str

def set_relax_param_file(output_path, sequence_data_path, tree_path, kappa, omega0, omega1, omega2, omega0_weight,
                         omega1_weight, selection_intensity_parameter, labels):
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

model1 = RELAX(kappa=<kappa>,p=<p>,omega1=<omega1>,omega2=<omega1>,k=1,theta1=<theta1>,theta2=<theta2>,frequencies=F0)
model2 = RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,frequencies=F0,k=<selection_intensity_parameter>)
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
    param_content = param_content.replace("<selection_intensity_parameter>", str(selection_intensity_parameter))
    param_content = param_content.replace("<labels>", labels)

    with open(output_path, "w") as output_file:
        output_file.write(param_content)

def set_traitrelax_param_file(output_dir, output_path, sequence_data_path, tree_path, character_data_path,
                              kappa, omega0, omega1, omega2, omega0_weight,
                              omega1_weight, selection_intensity_parameter, history_tree_path, labels_str, character_model_mu = None, character_model_pi0 = None):
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
model1 = RELAX(kappa=<kappa>,p=<p>,omega1=<omega1>,omega2=<omega2>,k=1,theta1=<theta1>,theta2=<theta2>,frequencies=F0)
model2 = RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,frequencies=F0,k=<selection_intensity_parameter>)

# ----------------------------------------------------------------------------------------
#                                    optimization parameters
# ----------------------------------------------------------------------------------------

optimization.tolerance = 0.000001
optimization.max_number_f_eval = 10000
optimization = FullD(derivatives=Newton,nstep=10)
optimization.final = powell
#optimization.message_handler = std
optimization.ignore_parameters = RELAX.k_1,BrLen
optimization.scale.tree = 0
optimization.advanced = 1


# ----------------------------------------------------------------------------------------
#                                    output files data                                    
# ----------------------------------------------------------------------------------------

joint.likelihood.debug = 0
optimization.profiler = <log_path>
output.tree.file = <expected_history_path>
output.debug.dir = <debug_dir>
optimization.backup.file = <backup_dir>
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
    backup_dir = output_dir + "inferred_values"
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
    param_content = param_content.replace("<selection_intensity_parameter>", str(selection_intensity_parameter))
    param_content = param_content.replace("<log_path>", log_path)
    param_content = param_content.replace("<expected_history_path>", expected_history_path)
    param_content = param_content.replace("<debug_dir>", debug_dir)
    param_content = param_content.replace("<backup_dir>", backup_dir)
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
character_model.pi0 = <pi0>
'''
        param_content = param_content.replace("<mu>", str(character_model_mu))
        param_content = param_content.replace("<pi0>", str(character_model_pi0))

    with open(output_path, "w") as output_file:
        output_file.write(param_content)

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='simulates alignments and character history under the TraitRELAX null and alternative models using Bio++ and INDELible')
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
                        default=1)
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
                        default=1)
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
                        default=0.5)

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
    num_of_replicates = int(args.num_of_replicates)
    # tree_path = args.tree_path
    # fix_tree_format(tree_path)
    # if not os.path.exists(tree_path):
    #     print("tree path " + tree_path + " does not exist")
    #     exit(1)
    trees_dir = args.tree_dir
    trees_paths = [path for path in os.listdir(trees_dir) if not ".nfs" in path]
    character_model_mu = float(args.character_model_mu)
    character_model_pi0 = float(args.character_model_pi0)
    kappa = float(args.kappa)
    omega0 = float(args.omega0)
    omega1 = float(args.omega1)
    omega2 = float(args.omega2)
    omega0_weight = float(args.omega0_weight)
    omega1_weight = float(args.omega1_weight)
    selection_intensity_parameter = float(args.selection_intensity_parameter)
    aln_len = int(args.aln_len)

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

    nuc1_theta = float(args.nuc1_theta)
    nuc1_theta1 = float(args.nuc1_theta1)
    nuc1_theta2 = float(args.nuc1_theta2)
    nuc2_theta = float(args.nuc2_theta)
    nuc2_theta1 = float(args.nuc2_theta1)
    nuc2_theta2 = float(args.nuc2_theta2)
    nuc3_theta = float(args.nuc3_theta)
    nuc3_theta1 = float(args.nuc3_theta1)
    nuc3_theta2 = float(args.nuc3_theta2)

    initial_nuc1_theta = float(args.initial_nuc1_theta)
    initial_nuc1_theta1 = float(args.initial_nuc1_theta1)
    initial_nuc1_theta2 = float(args.initial_nuc1_theta2)
    initial_nuc2_theta = float(args.initial_nuc2_theta)
    initial_nuc2_theta1 = float(args.initial_nuc2_theta1)
    initial_nuc2_theta2 = float(args.initial_nuc2_theta2)
    initial_nuc3_theta = float(args.initial_nuc3_theta)
    initial_nuc3_theta1 = float(args.initial_nuc3_theta1)
    initial_nuc3_theta2 = float(args.initial_nuc3_theta2)

    scaling_factor = float(args.scaling_factor)

    rep_regex = re.compile("(\d*).nwk", re.MULTILINE | re.DOTALL)
    for tree_filepath in os.listdir(trees_dir):
        tree_path = trees_dir + tree_filepath
        if not ".nwk" in tree_path:
            continue
        rep = int(rep_regex.search(tree_filepath).group(1))
        if int(rep) >= num_of_replicates:
            continue
        fix_tree_format(tree_path)
        scaled_tree_path = tree_path.replace(".nwk", "_scaled_by_" + str(scaling_factor) + ".nwk")
        scale_tree(tree_path, scaled_tree_path, scaling_factor=scaling_factor)
        print("**** simulating replicate " + str(rep) + " ****")
        # set simulation output directory
        simulation_output_dir = output_dir + "replicate_" + str(rep) + "/"
        if not os.path.exists(simulation_output_dir):
            res = os.system("mkdir -p " + simulation_output_dir)
        # simulate character data
        true_history_path, character_data_path, history_tree_path = simulate_character_data(character_model_mu,
                                                                                            character_model_pi0,
                                                                                            tree_path,
                                                                                            simulation_output_dir)
        # simulate sequence data
        sequence_data_path, labels_str = simulate_sequence_data(kappa, omega0, omega1, omega2, omega0_weight,
                                                                omega1_weight, selection_intensity_parameter,
                                                                scaled_tree_path, simulation_output_dir, 1, aln_len,
                                                                nuc1_theta, nuc1_theta1, nuc1_theta2, nuc2_theta,
                                                                nuc2_theta1, nuc2_theta2, nuc3_theta, nuc3_theta1,
                                                                nuc3_theta2)
        # set the parameters file for RELAX
        set_relax_param_file(relax_param_dir + str(rep) + ".bpp", sequence_data_path, scaled_tree_path, initial_kappa,
                             initial_omega0, initial_omega1, initial_omega2, initial_omega0_weight,
                             initial_omega1_weight, initial_selection_intensity_parameter, labels_str)
        # set parameters file for TraitRELAX
        set_traitrelax_param_file(simulation_output_dir + "traitrelax_result/",
                                  traitrelax_param_dir + str(rep) + ".bpp", sequence_data_path, scaled_tree_path,
                                  character_data_path,
                                  initial_kappa, initial_omega0, initial_omega1, initial_omega2, initial_omega0_weight,
                                  initial_omega1_weight, initial_selection_intensity_parameter, history_tree_path,
                                  labels_str, initial_character_model_mu, initial_character_model_pi0,)

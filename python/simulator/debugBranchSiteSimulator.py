import os, re, argparse
from Bio import AlignIO
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ete3 import Tree
plt.switch_backend('agg')
sns.set_style('whitegrid')
import numpy as np



def extract_omegas_assignment(omegas_assignments_filepath, simulation_source):
    omega_category_to_positions = dict()
    already_assigned = []
    if simulation_source == "INDELible":
        assignment_regex = re.compile("(\d+)\t(\d)\t\d\n")
        with open(omegas_assignments_filepath, "r") as omegas_assignments_file:
            omegas_assignments = omegas_assignments_file.read()
        for match in assignment_regex.finditer(omegas_assignments):
            position = int(match.group(1))
            omega_category = int(match.group(2))
            if not omega_category in omega_category_to_positions and not position in already_assigned:
                omega_category_to_positions[omega_category] = [position]
                already_assigned.append(position)
            elif not position in already_assigned:
                omega_category_to_positions[omega_category].append(position)
                already_assigned.append(position)
    else:
        with open(omegas_assignments_filepath, "r") as omegas_assignments_file:
            omegas_assignments_content = omegas_assignments_file.read()
        num_of_positions_regex = re.compile("num_of_pos=(\d*)")
        omega0_weight_regex = re.compile("p0=(\d*\.?\d*)")
        omega1_weight_regex = re.compile("p1=(\d*\.?\d*)")
        # process omegas weights and number of positions
        num_of_positions = int(num_of_positions_regex.search(omegas_assignments_content).group(1))
        omega0_weight = float(omega0_weight_regex.search(omegas_assignments_content).group(1))
        omega1_weight = float(omega1_weight_regex.search(omegas_assignments_content).group(1))
        omega2_weight = 1-omega0_weight-omega1_weight
        # divide the number of positions to partitions with increasing omega assignment per bin
        omega_category_to_positions[0] = []
        omega_category_to_positions[1] = []
        omega_category_to_positions[2] = []
        for position in range(1,num_of_positions+1):
            if position < int(num_of_positions*omega0_weight):
                omega_category_to_positions[0].append(position)
            elif position < int(num_of_positions*omega0_weight+num_of_positions*omega1_weight):
                omega_category_to_positions[1].append(position)
            else:
                omega_category_to_positions[2].append(position)
    return omega_category_to_positions


def extract_node_name(index, tree):
    i = 0
    for node in tree.traverse():
        if i == int(index):
            return node.name
        i += 1
    return 0


def extract_nodes_assignments(params_file_path, tree_path):
    nodes_assignment = dict()
    tree = Tree(tree_path,format=1)
    # since the tree in this case is always a star tree, each node is a leaf, and has a well defines classification into either BG or FG branch category
    with open(params_file_path, "r") as params_file:
        content = params_file.read()
    BG_nodes_regex = re.compile("model1.nodes_id\s*=\s*([\d*\,]*)", re.MULTILINE | re.DOTALL)
    FG_nodes_regex = re.compile("model2.nodes_id\s*=\s*([\d*\,]*)", re.MULTILINE | re.DOTALL)
    BG_evolution_size = 0
    FG_evolution_size = 0
    BG_nodes_indices = (BG_nodes_regex.search(content).group(1)).split(",")
    FG_nodes_indices = (FG_nodes_regex.search(content).group(1)).split(",")
    for index in BG_nodes_indices:
        node_name = extract_node_name(int(index)+1, tree) # not the correct index - 0 is root and is the last to be encountered - go by the params file
        nodes_assignment[node_name] = "BG"
        try:
            BG_evolution_size += (tree.search_nodes(name=node_name)[0]).dist
        except:
            if not node_name == tree.get_tree_root().name:
                print("node_name: ", node_name)
                exit(1)
    for index in FG_nodes_indices:
        node_name = extract_node_name(int(index)+1, tree) # not the correct index - 0 is root and is the last to be encountered - go by the params file
        nodes_assignment[node_name] = "FG"
        try:
            FG_evolution_size += (tree.search_nodes(name=node_name)[0]).dist
        except:
            if not node_name == tree.get_tree_root().name:
                print("node_name: ", node_name)
                exit(1)
    return nodes_assignment, BG_evolution_size, FG_evolution_size


if __name__ == '__main__':

    # test values
    omega0 = 0.1
    omega1 = 1
    omega2 = 2
    global_prefix = "C:/Users/ItayMNB7/Desktop/debugSimulator/star_tree_1000_positions_bl_0.5/"
    output_path = global_prefix + "dNdS_fig.jpg"
    tree_path = global_prefix + "treeWithInternals.nwk"
    params_filepath = global_prefix + "RELAX_params_template.bpp"
    simulation_source = "INDELible"
    Ks = [0.2, 0.5, 1, 2]
    subplots_order = [(0, 0), (0, 1), (1, 0), (1, 1)]
    fig, axis = plt.subplots(2, 2, sharex="row", sharey="col")
    boxprops = dict(linestyle='-', linewidth=1, color='k')
    medianprops = dict(linestyle='-', linewidth=0, color='k')
    meanprops = dict(linestyle='-', linewidth=1, color='b')
    color = dict(boxes='k', whiskers='k', medians='w', caps='k', means='b')
    fig.add_subplot(111, frameon=False)
    plt.grid(False)

    for k in Ks:
        true_values = [omega0, omega0 ** k, omega1, omega1 ** k, omega2, omega2 ** k]
        prefix = global_prefix + "k_" + str(k) + "/"
        k_regex = re.compile("k_(\d*\.?\d*)")
        k = float(k_regex.search(prefix).group(1))
        msa_path = prefix+"sequence_data_1.fas"
        omegas_assignments_filepath = prefix+"sequence_data_RATES.txt"
        ancestral_sequences_path = prefix+"sequence_data_ANCESTRAL_1.fas"

        omega_category_to_positions = extract_omegas_assignment(omegas_assignments_filepath, simulation_source)
        branch_to_class, BG_evolution_size, FG_evolution_size = extract_nodes_assignments(params_filepath, tree_path)

        # process sequences
        branch_to_sequence = dict()
        # internal sequences
        with open(ancestral_sequences_path, "r") as root_sequence_file:
            root_sequence_content = AlignIO.read(root_sequence_file, "fasta")
            for record in root_sequence_content:
                branch_to_sequence[record.name] = record.seq.translate()
        # external sequences
        with open(msa_path, "r") as msa_file:
            msa = AlignIO.read(msa_file, "fasta")
            num_of_positions = int(msa.get_alignment_length() / 3)
            for record in msa:
                branch_to_sequence[record.name] = record.seq.translate()

        # reset dictionaries that map to position number of non-synonymous and synonymous transitions
        BG_position_to_dN_count = dict()
        BG_position_to_dS_count = dict()
        FG_position_to_dN_count = dict()
        FG_position_to_dS_count = dict()
        for position in range(1, num_of_positions+1):
            BG_position_to_dN_count[position] = 0
            BG_position_to_dS_count[position] = 0
            FG_position_to_dN_count[position] = 0
            FG_position_to_dS_count[position] = 0

            # for each branch and position, count synonymous and non-synonymous transitions
            tree = Tree(tree_path, format=1)
            for node in tree.traverse():
                if not node is tree.get_tree_root():
                    child_sequence = branch_to_sequence[node.name]
                    parent_sequence = branch_to_sequence[node.up.name]
                    category = branch_to_class[node.name]
                    if child_sequence[position-1] != parent_sequence[position-1]:
                        if category == "BG":
                            BG_position_to_dN_count[position] = BG_position_to_dN_count[position] + 1
                        else:
                            FG_position_to_dN_count[position] = FG_position_to_dN_count[position] + 1
                    else:
                        if category == "BG":
                            BG_position_to_dS_count[position] = BG_position_to_dS_count[position] + 1
                        else:
                            FG_position_to_dS_count[position] = FG_position_to_dS_count[position] + 1

        # for each position, convert the dN and dS counts to a dN/dS ratio per branch category
        BG_position_to_dNdS = dict()
        FG_position_to_dNdS = dict()
        for position in range(1, num_of_positions+1):
            if BG_position_to_dS_count[position] > 0:
                BG_position_to_dNdS[position] = BG_position_to_dN_count[position] / BG_position_to_dS_count[position]
            if FG_position_to_dS_count[position] > 0:
                FG_position_to_dNdS[position] = FG_position_to_dN_count[position] / FG_position_to_dS_count[position]

        # compute for each selective regime the dN/dS values under each branch category (FG or BG) ofr the different relevant positions
        # print(omega_category_to_positions)
        dNdS_df = pd.DataFrame(columns=["branch_category", "omega_category", "dN/dS", "position"])
        for omega_category in [0, 1, 2]:
            for branch_category in ["BG", "FG"]:
                for position in omega_category_to_positions[omega_category]:
                    if position in BG_position_to_dNdS.keys() and position in FG_position_to_dNdS.keys():
                        df_record = dict()
                        df_record["branch_category"] = branch_category
                        df_record["omega_category"] = omega_category
                        df_record["position"] = position
                        if branch_category == "BG":
                            df_record["dN/dS"] = BG_position_to_dNdS[position]
                        else:
                            df_record["dN/dS"] = FG_position_to_dNdS[position]
                        dNdS_df = dNdS_df.append(df_record, ignore_index=True)

        # plot a boxplot of the values, groups by (level 1: omega category, level 2: branch category)
        index = Ks.index(k)
        bp = dNdS_df.boxplot(column="dN/dS", by=["omega_category", "branch_category"], grid=False, showmeans=True, meanline=True, boxprops=boxprops,
                               medianprops=medianprops,
                               meanprops=meanprops, ax=axis[subplots_order[index][0],subplots_order[index][1]], fontsize=8)
        sfig = axis[subplots_order[index][0], subplots_order[index][1]].get_figure()
        sfig.suptitle('')
        axis[subplots_order[index][0], subplots_order[index][1]].set_xlabel('')
        axis[subplots_order[index][0], subplots_order[index][1]].scatter(
            x=axis[subplots_order[index][0], subplots_order[index][1]].get_xticks(), y=np.array(true_values),
            label="simulated value", edgecolors=None, color="g", marker="s")
        axis[subplots_order[index][0], subplots_order[index][1]].set_ylim([0, max(dNdS_df["dN/dS"]) + 0.5])
        axis[subplots_order[index][0], subplots_order[index][1]].title.set_fontsize(8)
        axis[subplots_order[index][0], subplots_order[index][1]].title.set_text("simulated k = " + str(k))


    handles, labels = axis[subplots_order[0][0], subplots_order[0][1]].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center')
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel('simulated values')
    plt.ylabel('estimated values')
    plt.subplots_adjust()
    plt.tight_layout()
    plt.savefig(output_path)

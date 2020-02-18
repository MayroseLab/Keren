import re, os, argparse
from ete3 import Tree
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd


#### functions for histories in R format (simmap) ####

# convert history to tree in simmap format
def convert_history_to_simmap(history, base_tree, output_path):
    history_str = base_tree.write(format=1)  # initiate string to be the basic newick format (with no branch lengths) of the base tree
    for original_node in base_tree.traverse():
        if not original_node.is_root():
            original_parent = original_node.up
            history_node = history.search_nodes(name=original_node.name)[0]
            history_parent = history.search_nodes(name=original_parent.name)[0]
            branch_history = []  # list of history steps in reverse order. each history step is a tuple of (state, time_spent)
            current_node = get_next_history_step(history_node, history_parent)
            while current_node != None and current_node.name != history_node.name:
                label = current_node.label
                if label == "FG":
                    state = 1
                else:
                    state = 0
                time_spent = "{0:.6f}".format(current_node.dist)
                history_step = (state, time_spent)
                branch_history.append(history_step)
                history_parent = current_node
                current_node = get_next_history_step(history_node, history_parent)
            # document the first step in the history
            label = history_node.label
            if label == "FG":
                state = 1
            else:
                state = 0
            time_spent = "{0:.6f}".format(history_node.dist)
            branch_history.append((state, time_spent))
            # convert branch_history to simmap format and concat it to the original node name
            '''((A:{0,0.3:1,0.4:0,0.3},B:{0,1.0}):{0,1.0},(C:{0,0.25:1,0.75},D:{0,1.0}):{0,1.0});
                        branch lengths have been replaced by { curly brackets } according to the following convention: { state i, time spent in state i : state j, time spent in state j : etc. }.'''
            branch_history_str = ""
            for i in range(len(branch_history) - 1, 0, -1):
                history_step = branch_history[i]
                branch_history_str = branch_history_str + str(history_step[0]) + "," + str(history_step[1]) + ":"
            last_history_step = branch_history[0]
            branch_history_str = "{" + branch_history_str + str(last_history_step[0]) + "," + str(
                last_history_step[1]) + "}"
            bl = str(original_node.dist) if not str(original_node.dist).endswith(".0") else str(int(original_node.dist))
            original_node_name = re.compile(original_node.name + ":" + bl)
            if original_node.is_leaf():
                new_node_name = original_node.name + ":" + branch_history_str
            else:  # if the node is internal, exclude its name
                new_node_name = ":" + branch_history_str
            history_str = original_node_name.sub(new_node_name, history_str) # bug here - replacement isn't working

    with open(output_path, "w") as output:
        output.write(history_str)
    return 0


# parse branch history in simmap format
def parse_branch_history(tree_segment, node_history, labels_translator, internals_counter):
    first = 1
    for step in node_history:
        last = 1 if node_history.index(step) == len(node_history) - 1 else 0
        step_components = step.split(",")
        label = labels_translator[step_components[0]]
        length = step_components[1]
        if first:
            if last:
                tree_segment = tree_segment + "{" + label + "}:" + length
            else:
                tree_segment = "(" + tree_segment + "{" + label + "}:" + length + ")"
            first = 0
        else:
            if last:
                tree_segment = tree_segment + "internal" + str(internals_counter) + "{" + label + "}:" + length
            else:
                tree_segment = "(" + tree_segment + "internal" + str(
                    internals_counter) + "{" + label + "}:" + length + ")"
            internals_counter += 1
    return tree_segment, internals_counter


# parse simmap tree into a ete3 tree instance with internal branches
def parse_simmap_history(input_path, base_tree_path):

    labels_translator = {"a": "BG", "b": "FG"}
    with open(input_path, "r") as simmap_file:
        simmap_tree = simmap_file.read()

    # process the base tree to get the internal nodes names
    with open(base_tree_path, "r") as base_tree_file:
        base_tree = base_tree_file.read().rstrip()

    # first, catch internal nodes with no name, and give them names
    node_with_no_name_regex = re.compile("\):{")
    node_name_regex = re.compile("\)(.*?):")
    for simmap_match, base_tree_match in zip(node_with_no_name_regex.finditer(simmap_tree), node_name_regex.finditer(base_tree)):
        simmap_match_str = simmap_match.group(0)
        node_name = base_tree_match.group(1)
        replacement = ")" + node_name + ":{"
        simmap_tree = simmap_tree.replace(simmap_match_str, replacement, 1)

    # start by incorporating the history of internal branches into the tree
    # solve iteratively by repeatedly handing the innermost "()" expressions and exchanging them by some non-"()" expression mapped to them
    # internal_env_regex = re.compile("\([^\(]*?\)")
    internal_branch_history_regex = re.compile("\(([^\)|\(]*?)\)(Node\d+):{(.*?)}") # problem is caused by extra () in internals processing
    match = internal_branch_history_regex.search(simmap_tree)
    replacements_map = dict()
    internals_counter = 0
    replacements_counter = 0
    while match != None:
        tree_segment = match.group(0)
        node_name = match.group(2)
        tree_segment = tree_segment.replace(node_name+":{"+match.group(3)+"}", node_name)
        node_history = match.group(3).split(":")
        tree_segment, internals_counter = parse_branch_history(tree_segment, node_history, labels_translator, internals_counter)
        replacements_map["replacement_"+str(replacements_counter)+"_"] = tree_segment
        simmap_tree = simmap_tree.replace(match.group(0), "replacement_"+str(replacements_counter)+"_")
        replacements_counter += 1
        match = internal_branch_history_regex.search(simmap_tree)

    # replace replacements with their corresponding tree segments
    while "replacement" in simmap_tree:
        for replacement in replacements_map.keys():
            simmap_tree = simmap_tree.replace(replacement, replacements_map[replacement])

    # now, incorporate history of external branches of the tree
    branch_history_regex = re.compile("([^\(|\)|,]*?):({([a|b],.*?[:|}])*)")
    for match in branch_history_regex.finditer(simmap_tree):
        tree_segment = match.group(0)
        tree_segment = tree_segment.replace(":" + match.group(2), "")
        history = match.group(2).replace("}","").replace("{","")
        node_history = history.split(":")
        tree_segment, internals_counter = parse_branch_history(tree_segment, node_history, labels_translator,
                                                               internals_counter)
        simmap_tree = simmap_tree.replace(match.group(0), tree_segment)

    # now, parse the history from the processed string
    history = Tree(simmap_tree, format=1)
    label_regex = re.compile("{(.*?)}")
    for node in history.traverse():
        if node.is_root():
            # determine the node label by the label of its single child
            node_child = node.get_children()[0]
            node_label = label_regex.search(node_child.name).group(1)
        else:
            node_label = label_regex.search(node.name).group(1)
        node.add_feature("label", node_label)
        node.name = node.name.replace("{" + node_label + "}", "")
    return history # missing label for the root


#### functions for histories in Bio++ format ####

# set the name of the root node
def set_root_name(tree):
    nodes = []
    for node in tree.traverse("postorder"):
        nodes.append(node)
    root = tree.get_tree_root()
    root.name = "_baseInternal_" + str(nodes.index(root))


# parse history in biopp format
def parse_biopp_history(history_path):
    node_data_regex = re.compile("([^(|)]*?)\{(\d)\}")
    # read the tree from the file
    history = Tree(history_path, format=1)
    for node in history.traverse():
        if node != history.get_tree_root():
            node_name = (node_data_regex.search(node.name)).group(1)
            node_state = (node_data_regex.search(node.name)).group(2)
            node.name = node_name
            if node_state == "0":
                node.add_feature("label", "BG")
            else:
                node.add_feature("label", "FG")
        else:
            node.add_feature("label", "BG") # root is always BG
        history.get_tree_root().name = "_baseInternal_30"
    return history


#### general functions for histories evaluation ####


# returns tree with union of nodes from history_1 and history_2, where each node has two features: hist_1_label, hist_2_label
def parse_union_tree(history_1, history_2, base_tree_path, debug=False):
    base_tree = Tree(base_tree_path, format=1)
    # add for debugging
    base_tree.get_tree_root().name = "_baseInternal_30"
    united_tree = Tree()
    united_tree.dist = 0  # initialize distance to 0
    united_tree.get_tree_root().name = history_1.get_tree_root().name # set the name of the root
    united_tree.add_feature("history_1_label", history_1.get_tree_root().label)
    united_tree.add_feature("history_2_label", history_2.get_tree_root().label)
    union_nodes_number = 0
    for original_node in base_tree.traverse("preorder"): # traverse the tree in pre-order to assure that for any visited node, its parent from the base branch is already in the united tree
        original_parent = original_node.up
        if original_parent != None:  # will be none only in the case the original node is the root
            if debug:
                print("handled branch: (", original_node.name, ",", original_parent.name, ")")
            curr_union_parent = united_tree.search_nodes(name=original_parent.name)[0]
            hist_1_done = True
            hist_1_curr_child = None
            hist_1_parent = history_1.search_nodes(name=original_parent.name)[0]  # need to check names consistency across the 3 trees
            for child in hist_1_parent.children:
                if len(base_tree.search_nodes(name=child.name)) == 0 and len(child.search_nodes(
                        name=original_node.name)) > 0:  # if the child is a root in a tree that holds the original child node, then this child must be on the branch of interest
                    hist_1_curr_child = child
                    hist_1_done = False
                    break
            if hist_1_done:
                hist_1_curr_child = history_1.search_nodes(name=original_node.name)[0]
            hist_1_current_label = hist_1_curr_child.label

            hist_2_done = True
            hist_2_curr_child = None
            hist_2_parent = history_2.search_nodes(name=original_parent.name)[0]  # need to check names consistency across the 3 trees
            for child in hist_2_parent.children:
                if len(base_tree.search_nodes(name=child.name)) == 0 and len(child.search_nodes(
                        name=original_node.name)) > 0:  # if the child is a root in a tree that holds the original child node, then this child must be on the branch of interest
                    hist_2_curr_child = child
                    hist_2_done = False
                    break
            if hist_2_done:
                hist_2_curr_child = history_2.search_nodes(name=original_node.name)[0]
            hist_2_current_label = hist_2_curr_child.label

            while not hist_1_done or not hist_2_done:

                hist_1_dist = float("inf")
                hist_2_dist = float("inf")
                if not hist_1_done:  # if there is a node closer to the original node in history 1 -> add it to the united tree first
                    hist_1_dist = hist_1_curr_child.get_distance(original_parent.name) - curr_union_parent.get_distance(
                        original_parent.name)
                if not hist_2_done:
                    hist_2_dist = hist_2_curr_child.get_distance(original_parent.name) - curr_union_parent.get_distance(
                        original_parent.name)

                if debug:
                    if not hist_1_done:
                        print("history 1 has current child of ", original_parent.name, ": ", hist_1_curr_child.name,
                              " with label: ", hist_1_current_label, " and distance from parent is: ", hist_1_dist)
                    if not hist_2_done:
                        print("history 2 has current child of ", original_parent.name, ": ", hist_2_curr_child.name,
                              " with label: ", hist_2_current_label, " and distance from parent is: ", hist_2_dist)


                # first, check if now the two current children have the same name, and if this name is in the base tree - exit
                if hist_1_curr_child.name == hist_2_curr_child.name and len(base_tree.search_nodes(name=hist_1_curr_child.name)) > 0:
                    break

                # else, at least one of the histories has more than one step to go before reaching the bottom of the branch
                if hist_1_dist < hist_2_dist:  # add the node from history 1 and travel down to the next node in history 1
                    if debug:
                        print("adding child from history 1 which precedes to the one from history 2")
                        print("the label of the added node in history 1 is: ", hist_1_curr_child.label)
                        print("the label of the added node in histroy 2 remains like papa: ", hist_2_current_label)
                    curr_union_parent = curr_union_parent.add_child(child=None,
                                                                    name="internal_" + str(union_nodes_number),
                                                                    dist=hist_1_dist,
                                                                    support=None)
                    curr_union_parent.add_feature("history_1_label", hist_1_curr_child.label)
                    curr_union_parent.add_feature("history_2_label", hist_2_current_label)
                    hist_1_parent = hist_1_curr_child
                    if len(hist_1_parent.children) == 1:
                        hist_1_curr_child = hist_1_parent.children[0]
                    else:
                        hist_1_done = True
                    if debug:
                        print("united tree is now: \n", united_tree)
                        if hist_1_done:
                            print("history 1 on the handled branch is complete")
                        else:
                            print("history 1 on the handled branch isn't complete yet")

                else:  # add the node from history 2 and travel down to the next node in history 2
                    if debug:
                        print("adding child from history 2 which precedes to the one from history 1")
                        print("the label of the added node in history 2 is: ", hist_2_curr_child.label)
                        print("the label of the added node in history 1 remains like papa: ", hist_1_current_label)
                    curr_union_parent = curr_union_parent.add_child(child=None,
                                                                    name="internal_" + str(union_nodes_number),
                                                                    dist=hist_2_dist)  # added as a new branch
                    curr_union_parent.add_feature("history_1_label", hist_1_current_label)
                    curr_union_parent.add_feature("history_2_label", hist_2_curr_child.label)
                    hist_2_parent = hist_2_curr_child
                    if len(hist_2_parent.children) == 1:
                        hist_2_curr_child = hist_2_parent.children[0]
                    else:
                        hist_2_done = True
                    if debug:
                        print("united tree is now: \n", united_tree)
                        if hist_2_done:
                            print("history 2 on the handled branch is complete")
                        else:
                            print("history 2 on the handled branch isn't complete yet")
                union_nodes_number += 1

            # now add the original node as the child of the current parent
            original_dist = original_node.dist
            residual = original_dist - curr_union_parent.get_distance(
                united_tree.search_nodes(name=original_parent.name)[0])
            curr_union_parent = curr_union_parent.add_child(child=None, name=original_node.name, dist=residual)
            curr_union_parent.add_feature("history_1_label", history_1.search_nodes(name=original_node.name)[0].label)
            curr_union_parent.add_feature("history_2_label", history_2.search_nodes(name=original_node.name)[0].label)

    return united_tree


# compute the level of disagreement between the histories
def compute_distance(history_1, history_2, base_tree_path, debug=False):
    # 4 categories: 1) A = in FG in history_1 and in history_2
    #               2) B = in FG in history_1 and in BG in history_2
    #               3) C = in BG in history_1 and in FG in history_2
    #               4) D = in BG in history_1 and in history_2
    #               we compute the distance as 1-((A+D)/(A+B+C+D))=1-((fraction of tree in disagreement)/(total branch length))
    united_tree = parse_union_tree(history_1, history_2, base_tree_path, debug=debug)
    # resent departments sizes
    A = B = C = D = 0
    for node in united_tree.traverse():
        if debug:
            print("handled node: ", node.name)
            print("label in history 1: ", node.history_1_label)
            print("label in history 2: ", node.history_2_label)
            print("branch length: ", node.dist)
        if node.history_1_label == "FG" and node.history_2_label == "FG":
            if debug:
                print("added to category A: FG in history 1, FG in history 2")
            A += node.dist
        elif node.history_1_label == "FG" and node.history_2_label == "BG":
            if debug:
                print("added to category B: FG in history 1, BG in history 2")
            B += node.dist
        elif node.history_1_label == "BG" and node.history_2_label == "FG":
            if debug:
                print("added to category C: BG in history 1, FG in history 2")
            C += node.dist
        else:
            if debug:
                print("added to category D: BG in history 1, BG in history 2")
            D += node.dist
    if debug:
        print("absolute value of disagreement: ", (B + C))
        print("relative value of disagreement: ", (B + C) / (A + B + C + D))
    return (B + C) / (A + B + C + D)


# compare distance as the level of disagreement on each branch in the base tree between the two provided histories
def compute_distance_by_branch(history_1, history_2, base_tree_path):
    united_tree = parse_union_tree(history_1, history_2, base_tree_path)
    base_tree = Tree(base_tree_path, format=1)
    distance_per_branch = dict()
    for node in base_tree.traverse():
        if not node.is_root():
            original_node_name = node.name
            original_parent_name = node.up.name
            node_in_united_tree = united_tree.search_nodes(name=original_node_name)[0]
            parent_in_united_tree = united_tree.search_nodes(name=original_parent_name)[0]
            # compute the disagreement between the two histories in the range [original_node_name,original_parent_name]
            A = B = C = D = 0
            while node_in_united_tree.name != parent_in_united_tree.name: # while you haven't reached the end of the branch yet
                if node_in_united_tree.history_1_label == "FG" and node_in_united_tree.history_2_label == "FG":
                    A += node.dist
                elif node_in_united_tree.history_1_label == "FG" and node_in_united_tree.history_2_label == "BG":
                    B += node.dist
                elif node_in_united_tree.history_1_label == "BG" and node_in_united_tree.history_2_label == "FG":
                    C += node.dist
                else:
                    D += node.dist
                node_in_united_tree = node_in_united_tree.up
            distance_per_branch[original_node_name] = (B + C) / (A + B + C + D)

    return distance_per_branch


# plot a density plot of the histories distances from the true history, and place the distance of the expected history from the true history on it
def plot_histories_evaluation(distances_vector, output_path):
    # sns.distplot(distances_vector, rug=False, kde=True, hist=False) the kde is weird
    weights = np.ones_like(distances_vector)/float(len(distances_vector))
    plt.hist(distances_vector, weights=weights)
    plt.savefig(output_path)


# compute total branch length
def compute_tree_size(tree):
    size = 0
    for node in tree.traverse():
        size += node.dist
    return size


# get the next step in a history of a branch
def get_next_history_step(node, parent):
    children = parent.children
    for child in children:
        if len(child.search_nodes(name=node.name)) > 0:
            return child
    return None


# debug union tree properties
def debug_union_tree(union_tree):
    for node in union_tree.traverse():
        print("node name: ", node.name)
        if not node.is_root():
            print("parent name: ", node.up.name)
        print("branch length: ", node.dist)
        print("label in true history: ", node.history_1_label)
        print("label in expected history: ", node.history_2_label)
        print("")


# count the number of transitions in a history
def count_transitions(history):
    transitions_num = 0
    for node in history.traverse():
        node_children = node.get_children()
        for child in node_children:
            if node.label != child.label:
                transitions_num += 1
    return transitions_num


def analyse_histories(input_dir, true_history_path, base_tree_path):
    histories_dir = input_dir
    # process all the histories
    base_tree = Tree(base_tree_path, format=1)
    # print("base tree size: ", compute_tree_size(base_tree))
    true_history = parse_biopp_history(true_history_path)
    # print("true history size: ", compute_tree_size(true_history))
    if os.path.exists(input_dir+"distances.csv"):
        distances = pd.read_csv(input_dir+"distances.csv")
    else:
        distances = pd.DataFrame(columns=["distance"])
        for history_filename in os.listdir(histories_dir):
            history_path = histories_dir + history_filename
            history = parse_biopp_history(history_path)
            # print("history size: ", compute_tree_size(history))
            values = {"distance": compute_distance(history, true_history, base_tree_path)}
            distances = distances.append(values, ignore_index=True)
        distances.to_csv(input_dir +"distances.csv")
    print("mean distance: ", np.mean(distances))
    print("std of distances: ", np.std(distances))
    return distances

if __name__ == '__main__':

    print("**********************************************************")
    print("* analysis of the expected histories given by TraitRELAX *")
    print("**********************************************************")

    colors = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]

    # for debugging
    global_input_dir = "C:/Users/ItayMNB7/Desktop/TraitRELAXEvaluation/mu_10_pi0_0.5_kappa_2_p_0.1_omega1_1_omega2_2_theta1_0.5_theta2_0.8/"
    output_path = global_input_dir + "_expected_histories_analysis.jpg"
    taxa_num_options = [16, 32, 64]
    positions_num_options = [300, 600, 1000]
    k_values_options = [0.2, 0.5, 1, 1.6, 2]
    combo_to_df = dict()
    for taxa_num in taxa_num_options:
        input_dir = global_input_dir + str(taxa_num) + "_taxa/"
        true_history_path = input_dir + "true_history.nwk"
        true_history_simmap_path = input_dir + "true_history_simmap.nwk"
        base_tree_path = input_dir + "baseTree.nwk"
        true_history = parse_biopp_history(true_history_path)
        base_tree = Tree(base_tree_path, format=1)
        for positions_num in positions_num_options:
            # now join the distances achieved by all the k values executions
            dataframes = []
            for k in k_values_options:
                expected_histories_dir = input_dir + str(positions_num) + "_codons/k_" + str(k) + "/expected_histories/"
                # convert_history_to_simmap(true_history, base_tree, true_history_simmap_path)
                dataframes.append(analyse_histories(expected_histories_dir, true_history_path, base_tree_path))
            combo_to_df[(taxa_num, positions_num)] = pd.concat(dataframes)

    # plot a histogram where: taxa num is the x axis, there is a bar for each positions num option, and the bar shows the mean distance
    labels = [str(positions_num) + " positions" for positions_num in positions_num_options] # label per taxa num
    hist_colors = [colors[i] for i in range(len(positions_num_options))]
    data = [] # divided by the value of k
    for taxa_num in taxa_num_options:
        for positions_num in positions_num_options:
            mean_distances = []
            df = combo_to_df[(taxa_num, positions_num)]
            mean_distances.append(df["distances"].mean())
            data.append(np.array(mean_distances))

    # plot the histogram
    fig, axis = plt.subplots(nrows=1, ncols=1)
    axis.hist(np.array(data), density=True, histtype='bar', color=hist_colors, label=labels)
    axis.legend(prop={'size': 10})
    axis.set_xlabel("Number of species")
    axis.set_xticks([16, 32, 64])
    axis.set_ylabel("Mean distance of expected history from the true history")
    plt.savefig(output_path)
    plt.clf()
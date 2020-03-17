import re, argparse, os
from ete3 import Tree
import pandas as pd

# parse history in biopp format
# history path can also be string of the tree
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
                node.add_feature("label", "0")
            else:
                node.add_feature("label", "1")
        else:
            node.add_feature("label", "0") # root is always BG
        # history.get_tree_root().name = "_baseInternal_30"
    return history

# returns tree with union of nodes from history_1 and history_2, where each node has two features: hist_1_label, hist_2_label
def parse_union_tree(history_1, history_2, base_tree_path, debug=False):
    base_tree = Tree(base_tree_path, format=1)
    # add for debugging
    # base_tree.get_tree_root().name = "_baseInternal_30"
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

# debug mapping function
def debug_mapping(mapping, base_tree_path):
    base_tree = Tree(base_tree_path, format=1)
    for node in base_tree.traverse("postorder"):
        parent = node.up
        node_in_mapping = mapping.search_nodes(name=node.name)[0]
        parent_in_mapping = mapping.search_nodes(name=parent.name)[0]
        dist = node.dist
        dist_in_mapping = node_in_mapping.get_distance(parent_in_mapping.name)
        if dist != dist_in_mapping:
            print("bug in stochastic mapping product at branch (", parent.name, ", ", node.name, ")")
            print("distance in base tree: ", dist)
            print("distance in mapping: ", dist_in_mapping)

# compute the level of disagreement between the histories
def compute_distance(history_1, history_2, base_tree_path, debug=False):
    # 4 categories: 1) A = in FG in history_1 and in history_2
    #               2) B = in FG in history_1 and in BG in history_2
    #               3) C = in BG in history_1 and in FG in history_2
    #               4) D = in BG in history_1 and in history_2
    #               we compute the distance as 1-((A+D)/(A+B+C+D))=1-((fraction of tree in disagreement)/(total branch length))
    try:
        united_tree = parse_union_tree(history_1, history_2, base_tree_path, debug=debug)
    except:
        print("failed to parse united tree")
        exit(1)
    # resent departments sizes
    A = B = C = D = 0
    for node in united_tree.traverse():
        if debug:
            print("handled node: ", node.name)
            print("label in history 1: ", node.history_1_label)
            print("label in history 2: ", node.history_2_label)
            print("branch length: ", node.dist)
        if node.history_1_label == "1" and node.history_2_label == "1":
            if debug:
                print("added to category A: FG in history 1, FG in history 2")
            A += node.dist
        elif node.history_1_label == "1" and node.history_2_label == "0":
            if debug:
                print("added to category B: FG in history 1, BG in history 2")
            B += node.dist
        elif node.history_1_label == "0" and node.history_2_label == "1":
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
                if node_in_united_tree.history_1_label == "1" and node_in_united_tree.history_2_label == "1":
                    A += node.dist
                elif node_in_united_tree.history_1_label == "1" and node_in_united_tree.history_2_label == "0":
                    B += node.dist
                elif node_in_united_tree.history_1_label == "0" and node_in_united_tree.history_2_label == "1":
                    C += node.dist
                else:
                    D += node.dist
                node_in_united_tree = node_in_united_tree.up
            distance_per_branch[original_node_name] = (B + C) / (A + B + C + D)

    return distance_per_branch

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

# compute the number of ancestral asssignments wihch are in disagreement among all the nodes of the base tree
def compute_asr_disagreement(history_1, history_2, base_tree_path, debug=False):
    base_tree = Tree(base_tree_path, format=1)
    asr_disagreement_count = 0
    nodes_count = 0
    for base_node in base_tree.traverse("postorder"):
        try:
            nodes_count += 1
            history_1_node = history_1.search_nodes(name=base_node.name)[0]
            history_2_node = history_2.search_nodes(name=base_node.name)[0]
            if history_1_node.label != history_2_node.label:
                if debug:
                    print("ASR disagreement in node: ", base_node.name)
                asr_disagreement_count += 1
        except Exception as e:
            print("failed to find asr for node: ", base_node.name)
            print("error: ", e)
            exit(1)
    if debug:
        print("asr_disagreement_count: ", asr_disagreement_count)
        print("nodes_count: ", nodes_count)
    return asr_disagreement_count / nodes_count



if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
    description='computes distances of mappings from true history')
    parser.add_argument('--mu', '-mu', help='mu', required=True)
    parser.add_argument('--taxa_num', '-tn', help='taxa_num', required=False, default=32)
    parser.add_argument('--positions_num', '-pn', help='positions_num', required=False, default=150)
    parser.add_argument('--k', '-k', help='k', required=False, default=2)
    parser.add_argument('--replicate', '-replicate', help='replicate', required=True)
    parser.add_argument('--mappings_num', '-mappings_num', help='mappings_num', required=True)
    parser.add_argument('--mappings_file_path', '-i', help='path of the mappings whose distances should be computed',required=True)
    parser.add_argument('--sampling_based_expected_mapping_path', '-se', help='path that holds the sampling based approximation expected mapping', required=Tree)
    parser.add_argument('--analytic_expected_mapping_path', '-ae', help='path that holds the analytic approximation expected mapping', required=Tree)
    parser.add_argument('--base_tree_path', '-b', help='path to base tree file', required=True)
    parser.add_argument('--true_history_path', '-t', help='path to true history file', required=True)
    parser.add_argument('--output_path', '-o', help='path of file to write the df in csv format to', required=True)


    args = parser.parse_args()
    mu = float(args.mu)
    taxa_num = int(args.taxa_num)
    positions_num = int(args.positions_num)
    k = float(args.k)
    replicate = int(args.replicate)
    mappings_num = int(args.mappings_num)
    mappings_file_path = args.mappings_file_path
    sampling_based_expected_mapping_path = args.sampling_based_expected_mapping_path
    analytic_expected_mapping_path = args.analytic_expected_mapping_path
    base_tree_path = args.base_tree_path
    true_history_path = args.true_history_path
    output_path = args.output_path

    # if not os.path.exists(output_path):
    # reset dataframe
    df = pd.DataFrame(
        columns=["mu", "taxa_num", "positions_num", "k", "replicate", "mappings_num", "sampling_based_expected_history_distance",
                  "analytic_expected_history_distance", "mapping_order", "mapping_distance", "analytic_vs_sampling_based_expected_mapping_distance", "analytic_vs_sampling_based_expected_mapping_asr_disagreement"])
    values = {"mu": mu, "taxa_num": taxa_num, "positions_num": positions_num, "k": k, "replicate": replicate, "mappings_num": mappings_num}
    true_history = parse_biopp_history(true_history_path)
    with open(mappings_file_path, "r") as mappings_input_file:
        mappings_reads = mappings_input_file.readlines()
        mappings = [mapping for mapping in mappings_reads if mapping != "\n"]
    if len(mappings) < mappings_num:
        print("not enough mappings in output at: ", mappings_file_path)
    values["mappings_num"] = mappings_num
    sampling_based_expected_history = parse_biopp_history(sampling_based_expected_mapping_path)
    analytic_expected_history = parse_biopp_history(analytic_expected_mapping_path)
    values["sampling_based_expected_history_distance"] = compute_distance(sampling_based_expected_history, true_history, base_tree_path)
    values["analytic_expected_history_distance"] = compute_distance(analytic_expected_history, true_history, base_tree_path)
    for mapping_order in range(mappings_num):
        values["mapping_order"] = mapping_order
        mapping_str = mappings[mapping_order]
        mapping = parse_biopp_history(mapping_str)
        values["mapping_distance"] = compute_distance(mapping, true_history, base_tree_path)
        # compute the distance between the sampling based and analyic expected mappings
        values["analytic_vs_sampling_based_expected_mapping_distance"] = compute_distance(sampling_based_expected_history, analytic_expected_history, base_tree_path)
        values["analytic_vs_sampling_based_expected_mapping_asr_disagreement"] = compute_asr_disagreement(sampling_based_expected_history,analytic_expected_history,base_tree_path)
        df = df.append(values, ignore_index=True)

    print("distance between the analytic and sampling based expected mappings: ", values["analytic_vs_sampling_based_expected_mapping_distance"])
    print("fraction of disagreement between the analytic and sampling based expected mappings with respect to ASR: ", values["analytic_vs_sampling_based_expected_mapping_asr_disagreement"])
    df.to_csv(output_path)




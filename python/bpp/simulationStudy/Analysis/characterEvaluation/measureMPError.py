import re, argparse, os
from ete3 import Tree
import pandas as pd
import numpy as np

# add internal node names to base tree
def create_base_tree(input_tree_path, output_path):
    tree = Tree(input_tree_path, format=1)
    for node in tree.traverse():
        if "mapping" in node.name:
            node.get_children()[0].dist += node.dist # added the length of the removed branch to its single child
            node.delete()
        else:
            node_name = node.name
            node_name = node_name.replace("{0}", "")
            node_name = node_name.replace("{1}", "")
            node.name = node_name
    # give names to internal nodes
    tree.get_tree_root().name = "root"
    tree.write(outfile=output_path, format=1)
    return output_path

# parse history in biopp format
# history path can also be string of the tree
def parse_biopp_history(history_path, base_tree_path):
    base_tree = Tree(base_tree_path, format=1)
    node_data_regex = re.compile("([^(|)]*?)\{(\d)\}")
    # read the tree from the file
    history = Tree(history_path, format=1)
    for node in history.traverse("postorder"):
        if node != history.get_tree_root():
            node_name = (node_data_regex.search(node.name)).group(1)
            node_state = (node_data_regex.search(node.name)).group(2)
            if "missing_node" in node_name:
                parent = base_tree.search_nodes(name=node.get_children()[0].name)[0].up
                node_name = parent.name
            node.name = node_name
            if node_state == "0":
                node.add_feature("label", "0")
            else:
                node.add_feature("label", "1")
        else:
            # check if the root has more than 2 children, and if yes, reroot the tree in accordance with the base tree
            if len(node.get_children()) > 2:
                original_children = base_tree.get_tree_root().get_children()
                current_children = history.get_children()
                missing_node = original_children[0]
                missing_children = []
                in_curr = False
                for orig_node in original_children:
                    for curr_node in current_children:
                        if orig_node.name in curr_node.name:
                            in_curr = True
                    if not in_curr:
                        missing_node = orig_node
                    in_curr = False
                apparent_node = [node for node in original_children if not node == missing_node][0]
                missing_children = [node for node in current_children if not apparent_node.name in node.name]
                new = history.get_tree_root().add_child(child=None, name=missing_node.name, dist=missing_node.dist, support=None)
                # now make the two extra children of the current root the children of the new node
                for child in missing_children:
                    child.detach()
                    new.add_child(child=child)
                # set the label as pwe the mp solution
                if new.get_children()[0].label == new.get_children()[1].label:
                    new.add_feature("label", new.get_children()[0].label)
                else:
                    apparent_node_in_history = [node for node in history.get_tree_root().get_children() if apparent_node.name in node.name][0]
                    new.add_feature("label", apparent_node_in_history.label)
                node.add_feature("label", new.label)  # root is always BG
            else:
                node.add_feature("label", node.get_children()[0].label)
    history.get_tree_root().name = "root"
    return history

# returns tree with union of nodes from history_1 and history_2, where each node has two features: hist_1_label, hist_2_label
def parse_union_tree(history_1, history_2, base_tree_path, debug=False):
    base_tree = Tree(base_tree_path, format=1)
    base_tree.get_tree_root().name = "root"
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
            curr_union_parent = united_tree.search_nodes(name=original_parent.name.rstrip())[0]
            hist_1_done = True
            hist_1_curr_child = None
            hist_1_parent = history_1.search_nodes(name=original_parent.name.rstrip())[0]  # need to check names consistency across the 3 trees
            for child in hist_1_parent.children:
                if len(base_tree.search_nodes(name=child.name)) == 0 and len(child.get_children()) == 1 and child.get_children()[0].name == original_node.name:  # if the child does not exist in the base tree, it represents a mapping node that was created out of breaking a branch in the original tree
                    hist_1_curr_child = child
                    hist_1_done = False
                    break
            if hist_1_done:
                hist_1_curr_child = history_1.search_nodes(name=original_node.name.rstrip())[0]
            hist_1_current_label = hist_1_curr_child.label

            hist_2_done = True
            hist_2_curr_child = None
            hist_2_parent = history_2.search_nodes(name=original_parent.name.rstrip())[0]  # need to check names consistency across the 3 trees
            for child in hist_2_parent.children:
                if len(base_tree.search_nodes(name=child.name)) == 0 and len(child.get_children()) == 1 and child.get_children()[0].name == original_node.name:  #:  # if the child is a root in a tree that holds the original child node, then this child must be on the branch of interest
                    hist_2_curr_child = child
                    hist_2_done = False
                    break
            if hist_2_done:
                hist_2_curr_child = history_2.search_nodes(name=original_node.name.rstrip())[0]
            hist_2_current_label = hist_2_curr_child.label

            while not hist_1_done or not hist_2_done:

                if hist_1_curr_child.name == hist_2_curr_child.name: # both have reached the original child
                    print("error! original child wasn't recognized in the end of the loop")
                    exit(1)

                hist_1_dist = history_1.search_nodes(name=original_node.name.rstrip())[0].dist
                hist_2_dist = history_2.search_nodes(name=original_node.name.rstrip())[0].dist
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
                        print("the label of the added node in history 2 remains like papa: ", hist_2_current_label)
                    curr_union_parent = curr_union_parent.add_child(child=None,
                                                                    name="internal_" + str(union_nodes_number),
                                                                    dist=hist_1_dist,
                                                                    support=None)
                    curr_union_parent.add_feature("history_1_label", hist_1_curr_child.label)
                    curr_union_parent.add_feature("history_2_label", hist_2_current_label)
                    hist_1_parent = hist_1_curr_child
                    if len(hist_1_parent.children) == 1:
                        hist_1_curr_child = hist_1_parent.children[0]
                        if hist_1_curr_child.name == original_node.name:
                            hist_1_done = True
                    else:
                        hist_1_done = True
                    if debug:
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
                        if hist_2_curr_child.name == original_node.name:
                            hist_2_done = True
                    else:
                        hist_2_done = True
                    if debug:
                        if hist_2_done:
                            print("history 2 on the handled branch is complete")
                        else:
                            print("history 2 on the handled branch isn't complete yet")
                union_nodes_number += 1

            # now add the original node as the child of the current parent
            original_dist = original_node.dist
            residual = original_dist - curr_union_parent.get_distance(
                united_tree.search_nodes(name=original_parent.name.rstrip())[0])
            if residual < 0:
                print("error on residual computation for branch leading to ", original_node.name)
                print("original_dist: ", original_dist)
                print("curr_union_parent.get_distance(united_tree.search_nodes(name=original_parent.name.rstrip())[0]): ", curr_union_parent.get_distance(
                united_tree.search_nodes(name=original_parent.name.rstrip())[0]))
                exit(1)
            curr_union_parent = curr_union_parent.add_child(child=None, name=original_node.name, dist=residual)
            curr_union_parent.add_feature("history_1_label", history_1.search_nodes(name=original_parent.name.rstrip())[0].label)
            curr_union_parent.add_feature("history_2_label", history_2.search_nodes(name=original_parent.name.rstrip())[0].label)

    if debug:
        for node in united_tree.traverse("postorder"):
            print("node=", node.name)
            print("label in hist1=", node.history_1_label)
            print("label in hist2=", node.history_2_label)
            print("branch length=", node.dist)
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
    united_tree = parse_union_tree(history_1, history_2, base_tree_path, debug=debug)
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

# count the number of transitions in a history
def count_transitions(history):
    transitions_num = 0
    for node in history.traverse():
        node_children = node.get_children()
        for child in node_children:
            if node.label != child.label:
                transitions_num += 1
    return transitions_num

# main program
if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
    description='Analyses the estimated expected histories based on the multiple histories approximation and the true history based computation')
    parser.add_argument('--input_dir', '-i', help='directory that holds the the maximum parsimony solution based histories and the troe histories', required=True)
    parser.add_argument('--base_trees_dir', '-t', help='directory to the tree on which character evolution was simulated', required=True)
    parser.add_argument('--output_path', '-o', help='path that will hold a csv file with the distances between the true and mp history per replicate', required=True)

    args = parser.parse_args()
    input_dir = args.input_dir
    base_trees_dir = args.base_trees_dir
    output_path = args.output_path

    # initialize an output dataframe
    df = pd.DataFrame(columns=["tree_length", "mu", "#taxa", "k", "replicate", "expected(#transitions)", "simulated(#transitions)",  "mp(#transitions)", "distance(true_history,mp_history)"])
    tree_length_regex = re.compile("tbl_(.*?)_")
    mu_regex = re.compile("mu_(.*?)_")
    taxanum_regex = re.compile("([^/]*?)_taxa")
    k_regex = re.compile("k_(.*?)\/")
    replicate_regex = re.compile("replicate_(.*?)\/")

    tree_length = int(tree_length_regex.search(input_dir).group(1))
    mu = int(mu_regex.search(input_dir).group(1))
    expected_transitions_num = tree_length * mu
    taxa_num = int(taxanum_regex.search(input_dir).group(1))
    k = float(k_regex.search(input_dir).group(1))

    for path in os.listdir(input_dir):
        if "replicate" in path:
            record = {"tree_length": tree_length, "mu": mu, "#taxa": taxa_num, "k": k, "expected(#transitions)": expected_transitions_num}
            full_path = input_dir + path + "/"
            record["replicate"] = int(replicate_regex.search(full_path).group(1))
            if not os.path.exists(base_trees_dir):
                res = os.system("mkdir -p " + base_trees_dir)
            base_tree_path = create_base_tree(full_path + "character_data/true_history.nwk", base_trees_dir + str(record["replicate"]) + ".nwk")
            true_history = parse_biopp_history(full_path + "character_data/true_history.nwk", base_tree_path)
            mp_history = parse_biopp_history(full_path + "mp_data/mp_history.nwk", base_tree_path)
            record["simulated(#transitions)"] = count_transitions(true_history)
            record["mp(#transitions)"] = count_transitions(mp_history)
            record["distance(true_history,mp_history)"] = compute_distance(true_history, mp_history, base_tree_path)
            df = df.append(record, ignore_index=True)

    df.to_csv(output_path)

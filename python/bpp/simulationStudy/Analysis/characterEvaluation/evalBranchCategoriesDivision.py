import argparse, re, os
from ete3 import Tree
import numpy as np

#############################################################

mu = 8
tbl = 4
taxa_num = 32
pos_num = 300
pi0_options = [0.1, 0.3, 0.5, 0.9]
k_options = [0.5]

#############################################################

# count the number of transitions in a history
def count_transitions(history):
    transitions_num = 0
    for node in history.traverse():
        node_children = node.get_children()
        for child in node_children:
            if node.label != child.label:
                transitions_num += 1
    return transitions_num

# parse history in biopp format
def parse_biopp_history(history_path):
    node_data_regex = re.compile("([^(|)]*?)\{(\d)\}")
    # read the tree from the file
    print("history_path: ", history_path)
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
        # history.get_tree_root().name = "_baseInternal_30"
    return history

# computes relative frequencies of character states in the leaves
def eval_fraction_of_tip_states(history):
    tip_number = 0
    BG_tips = 0
    FG_tips = 0

    for leaf in history:
        tip_number += 1
        try:
            if leaf.label == "BG":
                BG_tips += 1
            elif leaf.label == "FG":
                FG_tips += 1
            else:
                print("Error: no category for leaf ", leaf.name)
                exit(1)
        except:
            print("leaf ", leaf.name, " has no label attribute")
            exit(1)

    BG_relative_tips = BG_tips / tip_number
    FG_relative_tips = FG_tips / tip_number
    if not BG_relative_tips + FG_relative_tips == 1:
        print("Error in evaluation: BG and FG fractions add up to ", str(BG_relative_tips + FG_relative_tips), " instead of 1")
        exit(1)
    return BG_relative_tips

# evaluate fractions of evolutions under BG and FG categories
def eval_division_ratio(history):

    overall_evolution = 0
    BG_evolution = 0
    FG_evolution = 0

    for node in history.traverse():
        evolution = node.dist
        try:
            category = node.label
            overall_evolution += evolution
            if category == "BG":
                BG_evolution += evolution
            elif category == "FG":
                FG_evolution += evolution
            else:
                print("Error: no category for branch ", node.name)
                exit(1)
        except:
            print("node ", node.name, " has no label attribute")
            exit(1)

    BG_relative_evolution = BG_evolution / overall_evolution
    FG_relative_evolution = FG_evolution / overall_evolution
    if not BG_relative_evolution + FG_relative_evolution > 0.98:
        print("Error in evaluation: BG and FG fractions add up to ", str(BG_relative_evolution + FG_relative_evolution)," instead of 1")
        exit(1)

    return BG_relative_evolution

#############################################################

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
    description='Analyses the estimated expected histories based on the multiple histories approximation and the true history based computation')
    parser.add_argument('--input_dir', '-i', help='directory that holds the true histories', required=False, default="/groups/itay_mayrose/halabikeren/TraitRELAX/simulations/newSimulationStudy/")


    args = parser.parse_args()
    input_dir = args.input_dir

    # traverse all pi0 combos and for each one - gather the measurements
    pi0_to_measurements = dict()
    for pi0 in pi0_options:
        pi0_to_measurements[pi0] = dict()
        (pi0_to_measurements[pi0])["tips_division"] = []
        (pi0_to_measurements[pi0])["evolution_division"] = []
        (pi0_to_measurements[pi0])["number_of_transitions"] = []
        for k in k_options:
            histories_input_dir = input_dir + "tbl_" + str(tbl) + "_mu_" + str(mu) + "_pi0_" + str(pi0) + "_kappa_2_p_0.1_omega1_1_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa_num) + "_taxa/" + str(pos_num) + "_codons/k_" + str(k) + "/"
            for i in range(50):
                history_path = histories_input_dir + "replicate_" + str(i) + "/character_data/true_history.nwk"
                history = parse_biopp_history(history_path)
                (pi0_to_measurements[pi0]["tips_division"]).append(eval_fraction_of_tip_states(history))
                (pi0_to_measurements[pi0]["evolution_division"]).append(eval_division_ratio(history))
                (pi0_to_measurements[pi0]["number_of_transitions"]).append(count_transitions(history))

    # now, report the mean and std of measurements for each pi0
    for pi0 in pi0_options:
        print("********** pi0 = " + str(pi0) + " **********")
        tips_division_data = pi0_to_measurements[pi0]["tips_division"]
        mean = np.mean(tips_division_data)
        std = np.std(tips_division_data)
        print("mean fraction of tip taxa in BG state: ", mean)
        print("std of fraction of tip taxa in BG state:", std)
        print("\n")
        evolution_division_data = pi0_to_measurements[pi0]["evolution_division"]
        mean = np.mean(evolution_division_data)
        std = np.std(evolution_division_data)
        print("mean fraction of evolution under BG state: ", mean)
        print("std of fraction of evolution under BG state:", std)
        print("\n")
        transitons_number_data = pi0_to_measurements[pi0]["number_of_transitions"]
        mean = np.mean(transitons_number_data)
        std = np.std(transitons_number_data)
        print("mean number of transitions: ", mean)
        print("std of number of transitions:", std)
        print("\n ************************* \n")
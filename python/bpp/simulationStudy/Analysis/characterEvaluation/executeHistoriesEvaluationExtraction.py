# likelihood based plot - read txt. ith logl values - the last corresponds to the expected history + fix the one of the true history based on the single output
import re, os, argparse, sys, math
from ete3 import Tree
import numpy as np
import seaborn as sns
sns.set_style('whitegrid')
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd
sys.path.append("/groups/itay_mayrose/halabikeren/myScripts/python/")
from utils.createJobFile import set_job_env, create_job_file

# add internal node names to base tree
def create_base_tree(input_tree_path, output_path):
    print('input_path: ', input_tree_path)
    tree = Tree(input_tree_path, format=1)
    for node in tree.traverse():
        if "mapping" in node.name:
            node.delete()
        else:
            node_name = node.name
            node_name = node_name.replace("{0}", "")
            node_name = node_name.replace("{1}", "")
            node.name = node_name
    tree.write(outfile=output_path, format=1)
    return output_path

# extract from execution output files the information
def extract_distance_values(mu, taxa_num, positions_num, k, input_dir, num_of_replicates):

    script_path = "/groups/itay_mayrose/halabikeren/myScripts/python/bpp/SimulationStudy/Analysis/characterEvaluation/compute_distances.py"
    for rep in range(num_of_replicates):
        distances_input_dir = input_dir + "replicate_" + str(rep) + "/"
        true_history_path = distances_input_dir + "/character_data/true_history.nwk"
        base_trees_dir = "/groups/itay_mayrose/halabikeren/TraitRELAX/simulations/historiesEvaluation/tbl_4_mu_" + str(mu) + "_pi0_0.5_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa_num) + "_taxa/trees_with_internal_names/"
        if not os.path.exists(base_trees_dir):
            res=os.system("mkdir -p " + base_trees_dir)
        base_tree_path = create_base_tree(true_history_path, base_trees_dir+str(rep) + ".nwk")
        error_files_dir = job_files_dir = distances_input_dir + "traitrelax_result/histories_evaluation/"
        for num_of_mappings in num_of_mappings_options:
            mappings_input_path = distances_input_dir + "traitrelax_result/histories_evaluation/" + str(num_of_mappings) + "_mappings_in_nwk.txt"
            sampling_based_expected_mapping_path = distances_input_dir + "traitrelax_result/histories_evaluation/" + str(num_of_mappings) + "_sampling_based_expected_mapping.nwk"
            analytic_expected_mapping_path = distances_input_dir + "traitrelax_result/histories_evaluation/" + str(
                num_of_mappings) + "_analytic_expected_mapping.nwk"
            output_path = distances_input_dir + "traitrelax_result/histories_evaluation/" + str(num_of_mappings) + "_mappings_distances.csv"
            # if not os.path.exists(output_path):
            cmd = "python " + script_path + " -mu " + str(mu) + " -tn " + str(taxa_num) + " -pn " + str(positions_num) + " -k " + str(k) + " -replicate " + str(rep) + " -mappings_num " + str(num_of_mappings) + " -i " + mappings_input_path + " -se " + sampling_based_expected_mapping_path + " -ae "  + analytic_expected_mapping_path + " -b " + base_tree_path + " -t " + true_history_path + " -o " + output_path
            # create job file with cmd
            job_name = "distances_computation_rep" + str(rep) + "_" + str(num_of_mappings) + "_mappings"
            file_name = job_name + ".sh"
            commands = ["module load python/python-anaconda3.6.5-michaldrori", cmd]  # ["export OMP_NUM_THREADS=4", cmd]
            touch_file_path = job_name + "_flag_done"
            full_job = create_job_file(job_name, commands, file_name, error_files_dir, job_files_dir,
                                           0, 1,
                                           touch_file_path, limit_nodes=False, python=False,
                                           openmpi=False, language="bash", queue="itaymr")
            res = os.system("qsub " + full_job)

# extract from job output files the log likelihood computations info
def extract_logl_values(mu, taxa_num, positions_num, k, data_dir, output_path):

    df = pd.DataFrame(
        columns=["mu", "taxa_num", "positions_num", "k", "replicate", "true_history_logl", "mappings_num", "exhaustive_computation_logl",
                 "sampling_based_expected_history_logl", "analytic_expected_history_logl", "mapping_order", "mapping_logl"])

    values = {"mu": mu, "taxa_num": taxa_num, "positions_num": positions_num, "k": k}
    replicate_regex = re.compile("replicate_(\d*)", re.MULTILINE | re.DOTALL)
    mappings_data_regex_str = "Analysis based on MAPPINGS_NUM mappings.*?\*{39}"
    character_logl_regex = re.compile("Character model log likelihood\:\s*\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    mapping_sequence_logl_regex = re.compile("Initializing data structure.*?\n(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    sampling_based_expected_sequence_logl_regex = re.compile("Sampling based analysis.*?Computing log likelihood based on the expected history approximation.*?Sequence Log likelihood\.*\:\s(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    analytic_expected_sequence_logl_regex = re.compile("Analytic rewards based analysis.*?Sequence Log likelihood\.*\:\s(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    true_history_joint_logl_regex = re.compile("Joint model log likelihood\:\s*\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)

    for path in os.listdir(data_dir):
        if ".OU" in path:
            fullpath = data_dir+path
            with open(fullpath, "r") as infile:
                content = infile.read()
                rep = int(replicate_regex.search(content).group(1))
                values["replicate"] = rep
                try:
                    # get the true history logl
                    values["true_history_logl"] = float(true_history_joint_logl_regex.search(content).group(1))
                    character_logl = float(character_logl_regex.search(content).group(1))
                    for num_of_mappings in num_of_mappings_options:
                        regex_str = mappings_data_regex_str.replace("MAPPINGS_NUM", str(num_of_mappings))
                        regex = re.compile(regex_str, re.MULTILINE | re.DOTALL)
                        mappings_data = regex.findall(content)[0]
                        # print("sp2")
                        values["mappings_num"] = num_of_mappings
                        logls = []
                        for match in mapping_sequence_logl_regex.finditer(mappings_data):
                            logls.append(float(match.group(1)))
                        order = 1
                        if len(logls) == 0:
                            print("no logl of mapping caught")
                            exit(1)
                        values["sampling_based_expected_history_logl"] = character_logl + float(sampling_based_expected_sequence_logl_regex.search(mappings_data).group(1))
                        values["analytic_expected_history_logl"] = character_logl + float(analytic_expected_sequence_logl_regex.search(mappings_data).group(1))
                        # compute the exhaustive based logl
                        max_sequence_logl = np.max(logls)
                        exhaustive_sequence_computation_logl = max_sequence_logl + math.log(np.sum([math.exp(mapping_logl - max_sequence_logl) for mapping_logl in logls])) - math.log(len(logls))
                        values["exhaustive_computation_logl"] = character_logl + exhaustive_sequence_computation_logl
                        if values["exhaustive_computation_logl"] == float('nan'):
                            continue
                        for logl in logls:
                            values["mapping_order"] = order
                            order += 1
                            values["mapping_logl"] = logl
                            df = df.append(values, ignore_index=True)

                except Exception as e:
                        print("failed to extract logl data data for mu=", mu, ", replicate=", rep)
                        print("input path: ", fullpath)
                        print("error: ", e)
                        exit(1)
                        continue
    # print("output_path: ", output_path)
    df.to_csv(output_path)
    return df


if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
    description='Analyses the estimated expected histories based on the multiple histories approximation and the true history based computation')
    parser.add_argument('--input_dir', '-i', help='directory that holds the stdout of the histories evaluation jobs', required=False, default="/groups/itay_mayrose/halabikeren/TraitRELAX/simulations/newSimulationStudy/historiesEvaluation/")
    parser.add_argument('--logl_data_dir', '-l', help='directory that holds the outputfiles of bio++ program debugexpectedhistory', required=True)
    parser.add_argument('--mu_options', '-mu', help='list of values of mu to include in the analysis', required=False, default=[1, 4, 8])
    parser.add_argument('--taxa_num_options', '-tn', help='list of taxa number values to include in the analysis', required=False, default=[32])
    parser.add_argument('--positions_num_options', '-pn', help='list of positions number values to include in the analysis', required=False, default=[300])
    parser.add_argument('--k_options', '-ko', help='list of k values to include in the analysis', required=False, default=[0.5])
    parser.add_argument('--num_of_replicates', '-rn', help='number of replicates to execute jobs on per combo', required=False, default=50)
    parser.add_argument('--num_of_mappings_options', '-nm', help='number of mappings combos used in the evaluation', required=False, default=[1000])
    parser.add_argument('--distance_analysis', '-da', help='1 is distance analysis s required, 0 else (will perform logl analysis', required=False, default=0)

    args = parser.parse_args()
    input_dir = args.input_dir
    logl_data_dir = args.logl_data_dir
    mu_options = args.mu_options
    taxa_num_options = args.taxa_num_options
    positions_num_options = args.positions_num_options
    k_options = args.k_options
    num_of_replicates = int(args.num_of_replicates)
    num_of_mappings_options = args.num_of_mappings_options
    distance_analysis = int(args.distance_analysis)

    #### distances based analysis ####

    # gather the data #
    if distance_analysis == 1:
        for mu in mu_options:
            for taxa_num in taxa_num_options:
                for positions_num in positions_num_options:
                    for k in k_options:
                        data_input_dir = input_dir + "tbl_4_mu_" + str(mu) + "_pi0_0.5_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa_num) + "_taxa/" + str(positions_num) + "_codons/k_" + str(k) + "/"
                        print("data_input_dir: ", data_input_dir)
                        extract_distance_values(mu, taxa_num, positions_num, k, data_input_dir, num_of_replicates)


    ### log likelihood based analysis ####

    # gather the data #
    else:
        print("analysis of logl data")
        for mu in mu_options:
            for taxa_num in taxa_num_options:
                for positions_num in positions_num_options:
                    for k in k_options:
                        data_dir = logl_data_dir + "tbl_4_mu_" + str(mu) + "_pi0_0.5_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa_num) + "_taxa/" + str(positions_num) + "_codons/k_" + str(k) + "/"
                        data_input_dir = input_dir + "tbl_4_mu_" + str(
                            mu) + "_pi0_0.5_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa_num) + "_taxa/" + str(positions_num) + "_codons/k_" + str(k) + "/"
                        print("data_input_dir: ", data_input_dir)
                        extract_logl_values(mu, taxa_num, positions_num, k, data_dir, data_input_dir+"logl_analysis.csv")


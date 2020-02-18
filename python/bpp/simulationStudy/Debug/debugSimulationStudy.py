import os, re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
plt.switch_backend('agg')
from scipy.stats.distributions import chi2
from mpl_toolkits.mplot3d import axes3d, Axes3D #<-- Note the capitalization!
import matplotlib.gridspec as gridspec
sns.set_style('whitegrid')
import math
from matplotlib.ticker import FormatStrFormatter, FuncFormatter

# hardcoded data
colors = ["steelblue", "darkorange", "g", "slateblue", "lightcoral"]
MAX_RECORDS_TO_PROCESS = 50
tbl_options = [1, 4, 8, 16, 32]
mu_options = [1, 2, 4, 8, 16, 32] # there is also 2 and 2 but not for FPR values
pi0_options = [0.1, 0.3, 0.5, 0.7, 0.9]
taxa_num_options = [16, 32, 64]
codon_positions_num_options = [150, 300, 600]
k_values_options = [0.2, 0.5, 0.8, 1, 1.2, 1.6, 2]
simulation_study_output_dir = "C:/Users/ItayMNB7/Google Drive/PhD/TraitRELAX/Results/Simulation"
TR_simulation_study_output_dir = simulation_study_output_dir + "/TraitRELAX"
R_simulation_study_output_dir = simulation_study_output_dir + "/RELAXWithTrueHistory"
R_MP_simulation_study_output_dir = simulation_study_output_dir + "/RELAXWithMPHistory"
output_dir = simulation_study_output_dir
grid_data_path = output_dir + "/grid_data.csv"

#####################################################################################

def convertToPercent(x, pos=0):
    return '%1.0f%%'%(100*x)

def doLRT(null_logl, alternative_logl, df=1):
    LR = 2 * (alternative_logl - null_logl)
    pvalue = chi2.sf(LR, df)
    return pvalue


def report_theoretical_test_result(combo_to_df, tbl, mu, pi0, taxa_num, positions_num):
    print("**** Theoretical analysis ****")
    df = combo_to_df[(tbl, mu, pi0, taxa_num, positions_num, 1)]
    print("number of replicates: ", df.shape[0])
    if df.shape[0] > 0:
        print("False positive rate: ", df[df.significant == True].shape[0] / df.shape[0])
    for k in k_values_options:
        if k != 1:
            try:
                df = combo_to_df[(tbl, mu, pi0, taxa_num, positions_num, k)]
                if df.shape[0] > 0:
                    print("For k = ", str(k), ": power = ", df[df.significant == True].shape[0] / df.shape[0])
            except:
                continue
    print("\n")


def find_empirical_LR_threshold(combo_to_df, tbl, mu, pi0, taxa_num, positions_num):

    print("**** Empirical analysis ****")

    # extract the data for which k=1
    null_df = combo_to_df[(tbl, mu, pi0, taxa_num, positions_num, 1)]
    if null_df.shape[0] == 0:
        return 0
    null_LRs = list(null_df["LR"])

    # check the value of the LR from the top in the location under which you have 5% of the samples - it will define the new LR threshold - reprot FDR as 5%
    null_LRs.sort()
    LR_cutoff = int(5/100*len(null_LRs))
    if LR_cutoff == 0:
        LR_cutoff += 1
    LR_threshold = null_LRs[len(null_LRs)-LR_cutoff] # any LR >= to LR_threshold will be considered significant
    FPR = null_df[null_df.LR >= LR_threshold].shape[0] / null_df.shape[0]
    print("LR threshold: ", LR_threshold, "\n")
    print("For k = 1: FPR = " + str(FPR))

    for k in k_values_options:
        if k != 1:
            try:
                alternative_df = combo_to_df[(tbl, mu, pi0, taxa_num, positions_num, k)]
                if alternative_df.shape[0] > 0:
                    power = alternative_df[alternative_df.LR >= LR_threshold].shape[0] / alternative_df.shape[0]
                    print("For k = " + str(k) + ": power = " + str(power))
            except:
                continue
    print("\n")
    return LR_threshold

#####################################################################################

# extract data from TraitRELAX execution output files in input_dir return a dataframe
def extract_TraitRELAX_data(input_dir):

    parameters = ["mu", "pi0", "kappa", "p", "omega1", "omega2", "theta1", "theta2", "k"]
    # set regex expressions
    days_regex = re.compile("Total execution time\: (\d*\.?\d*)d", re.MULTILINE | re.DOTALL)
    hours_regex = re.compile("Total execution time\:.*?(\d*\.?\d*)h", re.MULTILINE | re.DOTALL)
    minutes_regex = re.compile("Total execution time\:.*?(\d*\.?\d*)m", re.MULTILINE | re.DOTALL)
    seconds_regex = re.compile("Total execution time\:.*?(\d*\.?\d*)s", re.MULTILINE | re.DOTALL)
    replicate_regex = re.compile("Parsing.*?replicate_(\d*).*?fas", re.DOTALL | re.MULTILINE)
    initial_character_logl_regex = re.compile("Character Log likelihood\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    initial_sequence_logl_regex = re.compile("Sequence Log likelihood\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    initial_joint_logl_regex = re.compile("Overall Log likelihood\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    null_character_logl_regex = re.compile("Null model fitting.*?Character Log likelihood\.*\:\s*(-\d*\.?\d*)",
                                           re.MULTILINE | re.DOTALL)
    null_sequence_logl_regex = re.compile("Null model fitting.*?Sequence Log likelihood\.*\:\s*(-\d*\.?\d*)",
                                          re.MULTILINE | re.DOTALL)
    null_joint_logl_regex = re.compile("Null model fitting.*?Overall Log likelihood\.*\:\s*(-\d*\.?\d*)",
                                       re.MULTILINE | re.DOTALL)
    alternative_character_logl_regex = re.compile("Alternative model likelihood after optimization.*?Character Log likelihood\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    alternative_sequence_logl_regex = re.compile("Alternative model likelihood after optimization.*?Sequence Log likelihood\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL) # bug in some output files - covert by actually computing the overall-char
    alternative_joint_logl_regex = re.compile("Alternative model likelihood after optimization.*?Overall Log likelihood\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    num_of_opt_cycles_regex = re.compile("Number of optimization cycles\.*\:\s*(\d)")
    non_parametric_regex_t_colname = {"replicate": replicate_regex,
                                      "initial_character_logl": initial_character_logl_regex,
                                      "initial_sequence_logl": initial_sequence_logl_regex,
                                      "initial_joint_logl": initial_joint_logl_regex,
                                      "null_character_logl": null_character_logl_regex,
                                      "null_sequence_logl": null_sequence_logl_regex,
                                      "null_joint_logl": null_joint_logl_regex,
                                      "alternative_character_logl": alternative_character_logl_regex,
                                      "alternative_sequence_logl": alternative_sequence_logl_regex,
                                      "alternative_joint_logl": alternative_joint_logl_regex,
                                      "num_of_cycles": num_of_opt_cycles_regex}

    simulated_parameter_regex_str = "Parsing file\s*.*PARAMETER_(\d*\.*\d*).*?for options"
    null_parameter_regex_str = "Null model fitting.*?\.PARAMETER\.*\:\s*(\d*\.?\d*)"
    alternative_parameter_regex_str = "Alternative model likelihood after optimization.*?\w*\.PARAMETER(_\d)?\.*\:\s*(\d*\.?\d*)"

    # initialize dataframe
    colnames = list(non_parametric_regex_t_colname.keys())

    for parameter in parameters:
        colnames.append("simulated_"+parameter)
        colnames.append("null_"+parameter)
        colnames.append("alternative_" + parameter)
    # manually computad parameters
    manual_parameters = ["omega0", "p0", "p1"]
    for parameter in manual_parameters:
        colnames.append("simulated_"+parameter)
        colnames.append("null_"+parameter)
        colnames.append("alternative_" + parameter)
    colnames = colnames + ["LR", "pvalue", "significant", "correctDirection", "duration(hours)"]
    df = pd.DataFrame(columns=colnames)
    clean_df = pd.DataFrame(columns=colnames)
    records_num = 0

    # fill in the dataframe
    for file in os.listdir(input_dir):
        if ".OU" in file and ((not "k_1/" in input_dir and records_num < MAX_RECORDS_TO_PROCESS) or "k_1/" in input_dir):
            # print("analysing ", input_dir, file)
            with open(input_dir+file, "r") as input_file:
                content=input_file.read()
                values = dict()

                # extract the values for the non-parametric columns
                for colname in non_parametric_regex_t_colname.keys():
                    try:
                        if "logl" in colname:
                            values["alternative_character_logl"] = float(alternative_character_logl_regex.search(content).group(1))
                            values["alternative_joint_logl"] = float(alternative_joint_logl_regex.search(content).group(1))
                            values["alternative_sequence_logl"] = values["alternative_joint_logl"] -  values["alternative_character_logl"]
                    except:
                        print("unable to get log likelihood of alternative model in file ", input_dir+file)
                        exit(1)

                for colname in non_parametric_regex_t_colname.keys():
                    try:
                        regex = non_parametric_regex_t_colname[colname]
                        if not colname in values:
                            values[colname] = float(regex.search(content).group(1))
                    except:
                        print("failed to find " + colname + " match in " + input_dir+file + "-> setting as Nan")
                        values[colname] = float('nan')
                        exit(1)
                        continue

                # extract the parameters
                for parameter in parameters:

                    # try to extract the simulated
                    try:
                        regex_str = simulated_parameter_regex_str.replace("PARAMETER", parameter)
                        regex = re.compile(regex_str, re.MULTILINE | re.DOTALL)

                        values["simulated_" + parameter] = float(regex.search(content).group(1))
                    except:
                        print("failed to find simulated " + parameter + " match in " + input_dir + file + "-> setting as Nan")
                        values["simulated_" + parameter] = float('nan')
                        exit(1)

                    # try to extract the null
                    try:
                        regex_str = null_parameter_regex_str.replace("PARAMETER", parameter)
                        regex = re.compile(regex_str, re.MULTILINE | re.DOTALL)
                        values["null_" + parameter] = float(regex.search(content).group(1))
                    except:
                        print("failed to find null " + parameter + " match in " + input_dir + file + "-> setting as Nan")
                        values["null_" + parameter] = float('nan')
                        exit(1)

                    # try to extract the alternative
                    try:
                        regex_str = alternative_parameter_regex_str.replace("PARAMETER", parameter)
                        regex = re.compile(regex_str, re.MULTILINE | re.DOTALL)
                        values["alternative_" + parameter] = float(regex.search(content).group(2))
                    except:
                        print("failed to find " + parameter + " match in " + input_dir + file + "-> setting as Nan")
                        print("regex:\n", regex.pattern)
                        values["alternative_" + parameter] = float('nan')
                        exit(1)

                # add manual translation of parameters: omega0, p0, p1

                # omega0
                values["simulated_omega0"] = values["simulated_p"] * values["simulated_omega1"]
                values["null_omega0"] = values["null_p"] * values["null_omega1"]
                values["alternative_omega0"] = values["alternative_p"] * values["alternative_omega1"]

                # p0
                values["simulated_p0"] = values["simulated_theta1"]
                values["null_p0"] = values["null_theta1"]
                values["alternative_p0"] = values["alternative_theta1"]

                # p1
                values["simulated_p1"] = values["simulated_theta2"] * (1 - values["simulated_theta1"])
                values["null_p1"] = values["null_theta2"] * (1 - values["null_theta1"])
                values["alternative_p1"] = values["alternative_theta2"] * (1 - values["alternative_theta1"])

                # LR and pvalue
                values["LR"] = 2 * (values["alternative_joint_logl"] - values["null_joint_logl"])
                values["pvalue"] = chi2.sf(values["LR"], 1) # 1 degree of freedom is the diff in num of parameters between alternative and null models

                # significant: True of False
                if values["pvalue"] < 0.05:
                    values["significant"] = 1
                else:
                    values["significant"] = 0

                # correctDirection: compare to simulated k -> if both < 1 or > 1 True, else False
                if (values["simulated_k"] < 1 and values["alternative_k"] < 1) or (values["simulated_k"] > 1 and values["alternative_k"] >1):
                    values["correctDirection"] = 1
                else:
                    values["correctDirection"] = 0

                # add the overall duration
                duration = 0
                for match in days_regex.finditer(content):
                    duration += float(match.group(1))*24
                for match in hours_regex.finditer(content):
                    duration += float(match.group(1))
                for match in minutes_regex.finditer(content):
                    duration += float(match.group(1))*(1/60)
                for match in seconds_regex.finditer(content):
                    duration += float(match.group(1))*(1/60)*(1/60)
                values["duration(hours)"] = duration

                # add the record to the dataframe
                add = True
                for key in list(values.keys()):
                    if math.isnan(values[key]):
                        add = False
                # add the record to the dataframe
                if add:
                    if (values['null_joint_logl'] - values['alternative_joint_logl']) > 0.5:
                        print("optimization bug expressed in " + input_dir+file)
                        print("exhibit logl values:\n")
                        print("values['null_joint_logl']: ", values['null_joint_logl'])
                        print("values['alternative_joint_logl']: ", values['alternative_joint_logl'])
                        print("\n")
                        continue
                    else:
                        clean_df = clean_df.append(values, ignore_index=True)
                        records_num += 1
                    df = df.append(values, ignore_index=True)

    df.drop_duplicates(subset='replicate', keep="last")
    clean_df.drop_duplicates(subset='replicate', keep="last")
    if df.shape[0] * 0.5 > clean_df.shape[1]:
        print("PROBLEM! More than half the data lost in cleaning process")

    return clean_df, df

# extract data from RELAX execution output files in input_dir return a dataframe
def extract_RELAX_data(input_dir, tbl, mu, pi0, taxa_num, positions_num, k):

    parameters = ["kappa", "p", "omega1", "omega2", "theta1", "theta2", "k"]
    # set regex expressions
    days_regex = re.compile("Total execution time\: (\d*\.?\d*)d", re.MULTILINE | re.DOTALL)
    hours_regex = re.compile("Total execution time\:.*?(\d*\.?\d*)h", re.MULTILINE | re.DOTALL)
    minutes_regex = re.compile("Total execution time\:.*?(\d*\.?\d*)m", re.MULTILINE | re.DOTALL)
    seconds_regex = re.compile("Total execution time\:.*?(\d*\.?\d*)s", re.MULTILINE | re.DOTALL)
    replicate_regex = re.compile("Parsing.*?replicate_(\d*).*?fas", re.DOTALL | re.MULTILINE)
    initial_logl_regex = re.compile("Computing intial log likelihood.*?Log likelihood\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    null_logl_regex = re.compile("Fitting the null model.*?Log likelihood\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    alternative_logl_regex_str = "Fitting the alternative model.*?Log likelihood\.*\:\s*(-\d*\.?\d*)"
    non_parametric_regex_t_colname = {"replicate": replicate_regex,
                                      "initial_logl": initial_logl_regex,
                                      "null_logl": null_logl_regex,
                                      "alternative_logl": alternative_logl_regex_str}
    simulated_parameter_regex_str = "Parsing file\s*.*[_|\/]PARAMETER_(\d*\.*\d*).*?for options"
    null_parameter_regex_str = "Fitting the null model.*?model 2.*?.RELAX.PARAMETER\.*\:\s*(\d*\.?\d*)"
    alternative_parameter_regex_str = "Fitting the alternative model.*?model 2.*?.RELAX.PARAMETER\.*\:\s*(\d*\.?\d*)"
    # intialize dataframe
    colnames = ["tbl", "mu", "taxa_num", "pi0", "positions_num", "k"]  # keren - 4.7.19
    colnames = colnames + list(non_parametric_regex_t_colname.keys())

    for parameter in parameters:
        colnames.append("simulated_"+parameter)
        colnames.append("null_"+parameter)
        colnames.append("alternative_" + parameter)
    # manually computad parameters
    manual_parameters = ["omega0", "p0", "p1"]
    for parameter in manual_parameters:
        colnames.append("simulated_"+parameter)
        colnames.append("null_"+parameter)
        colnames.append("alternative_" + parameter)
    colnames = colnames + ["LR", "pvalue", "significant", "correctDirection", "duration(hours)"]
    # colnames.append("winning_starting_point") # keren - 3.7.19 - two sp execution
    df = pd.DataFrame(columns=colnames)
    clean_df = pd.DataFrame(columns=colnames)
    records_num = 0

    # fill in the dataframe
    for file in os.listdir(input_dir):
        if ".OU" in file and ((not "k_1/" in input_dir and records_num < MAX_RECORDS_TO_PROCESS) or "k_1/" in input_dir):

            with open(input_dir+file, "r") as input_file:
                content=input_file.read()
                if "Line minimization failed!" in content:
                    print("Bio++ optimization bug expressed in " + input_dir + file)
                    continue
                if "no model associated" in content:
                    continue
                values = dict()

                # keren - 4.7.19
                values["tbl"] = tbl
                values["mu"] = mu
                values["pi0"] = pi0
                values["taxa_num"] = taxa_num
                values["positions_num"] = positions_num
                values["k"] = k

                # extract the winning sp in the alternative model

                non_parametric_regex_t_colname["alternative_logl"] = re.compile(alternative_logl_regex_str, re.MULTILINE | re.DOTALL)

                # extract the values for the non-parametric columns
                for colname in non_parametric_regex_t_colname.keys():
                    try:
                        regex = non_parametric_regex_t_colname[colname]
                        values[colname] = float(regex.search(content).group(1))
                    except:
                        # print("failed to find " + colname + " match in " + input_dir+file + "-> setting as Nan")
                        values[colname] = float('nan')

                # extract the parameters
                for parameter in parameters:

                    # try to extract the simulated
                    try:
                        regex_str = simulated_parameter_regex_str.replace("PARAMETER", parameter)
                        regex = re.compile(regex_str, re.MULTILINE | re.DOTALL)

                        values["simulated_" + parameter] = float(regex.search(content).group(1))
                    except:
                        # print("failed to find simulated " + parameter + " match in " + input_dir + file + "-> setting as Nan")
                        values["simulated_" + parameter] = float('nan')

                    # try to extract the null
                    try:
                        regex_str = null_parameter_regex_str.replace("PARAMETER", parameter)
                        regex = re.compile(regex_str, re.MULTILINE | re.DOTALL)
                        values["null_" + parameter] = float(regex.search(content).group(1))
                    except:
                        # print("failed to find null " + parameter + " match in " + input_dir + file + "-> setting as Nan")
                        values["null_" + parameter] = float('nan')

                    # try to extract the alternative
                    try:
                        regex_str = alternative_parameter_regex_str.replace("PARAMETER", parameter)
                        regex = re.compile(regex_str, re.MULTILINE | re.DOTALL)
                        values["alternative_" + parameter] = float(regex.search(content).group(1))
                    except:
                        # print("failed to find alternative " + parameter + " match in " + input_dir + file + "-> setting as Nan")
                        values["alternative_" + parameter] = float('nan')

                # add manual translation of parameters: omega0, p0, p1

                # omega0
                values["simulated_omega0"] = values["simulated_p"] * values["simulated_omega1"]
                values["null_omega0"] = values["null_p"] * values["null_omega1"]
                values["alternative_omega0"] = values["alternative_p"] * values["alternative_omega1"]

                # p0
                values["simulated_p0"] = values["simulated_theta1"]
                values["null_p0"] = values["null_theta1"]
                values["alternative_p0"] = values["alternative_theta1"]

                # p1
                values["simulated_p1"] = values["simulated_theta2"] * (1 - values["simulated_theta1"])
                values["null_p1"] = values["null_theta2"] * (1 - values["null_theta1"])
                values["alternative_p1"] = values["alternative_theta2"] * (1 - values["alternative_theta1"])

                # LR and pvalue
                values["LR"] = 2 * (values["alternative_logl"] - values["null_logl"])
                values["pvalue"] = chi2.sf(values["LR"], 1) # 1 degree of freedom is the diff in num of parameters between alternative and null models

                # significant: True of False
                if values["pvalue"] < 0.05:
                    values["significant"] = 1
                else:
                    values["significant"] = 0

                # correctDirection: compare to simulated k -> if both < 1 or > 1 True, else False
                if (values["simulated_k"] < 1 and values["alternative_k"] < 1) or (values["simulated_k"] > 1 and values["alternative_k"] >1):
                    values["correctDirection"] = 1
                else:
                    values["correctDirection"] = 0

                # add the overall duration
                duration = 0
                for match in days_regex.finditer(content):
                    duration += float(match.group(1))*24
                for match in hours_regex.finditer(content):
                    duration += float(match.group(1))
                for match in minutes_regex.finditer(content):
                    duration += float(match.group(1))*(1/60)
                for match in seconds_regex.finditer(content):
                    duration += float(match.group(1))*(1/60)*(1/60)
                values["duration(hours)"] = duration



                add = True
                # add the record to the dataframe
                if add:
                    if (values['null_logl'] - values['alternative_logl']) > 0.5:
                        print("optimization bug expressed in " + input_dir + file)
                        print("values['null_logl']: ", values['null_logl'])
                        print("values['alternative_logl']: ", values['alternative_logl'])
                        # exit(1)
                    else:
                        clean_df = clean_df.append(values, ignore_index=True)
                    df = df.append(values, ignore_index=True)
                    records_num += 1

    df.drop_duplicates(subset='replicate', keep="last")
    clean_df.drop_duplicates(subset='replicate', keep="last")
    if df.shape[0] * 0.5 > clean_df.shape[0]:
        print("PROBLEM! More than half the data lost in cleaning process")

    return clean_df, df

#####################################################################################

def compare_null_likelihoods(TraitRELAX_DF, RELAX_DF):

        nullCompToCombo = dict()
        print("combo\tmean_diff\t#TraitRELAX_wins\t#RELAX_wins")
        tbl_to_mu = {1:32, 4:1, 4:4, 4:8, 8:4, 16:2, 32:1}
        for tbl in tbl_to_mu:
            mu = tbl_to_mu[tbl]
            for taxa_num in taxa_num_options:
                for codon_positions_num in codon_positions_num_options:
                    if (tbl, mu, 0.5, taxa_num, codon_positions_num, 1) in TraitRELAX_DF and (tbl, mu, 0.5, taxa_num, codon_positions_num, 1) in RELAX_DF:
                        TraitRELAX_data = TraitRELAX_DF[(tbl, mu, 0.5, taxa_num, codon_positions_num, 1)]
                        RELAX_data = RELAX_DF[(tbl, mu, 0.5, taxa_num, codon_positions_num, 1)]
                        TraitRELAX_dataset_to_null_likelihood = TraitRELAX_data[["replicate", "null_sequence_logl"]]
                        RELAX_dataset_to_null_likelihood = RELAX_data[["replicate", "null_logl"]]
                        merged_data = pd.merge(TraitRELAX_dataset_to_null_likelihood, RELAX_dataset_to_null_likelihood, on="replicate")
                        diffs = merged_data["null_sequence_logl"] - merged_data["null_logl"]
                        mean_diff = np.mean(diffs)
                        num_of_traitrelax_wins = len([diff for diff in diffs if diff > 0])
                        num_of_relax_wins = len([diff for diff in diffs if diff < 0])
                        nullCompToCombo[(tbl, mu, pi0, taxa_num, codon_positions_num, k_value)] = (mean_diff, num_of_traitrelax_wins, num_of_relax_wins)
                        print(",".join([str(tbl), str(mu), str(0.5), str(taxa_num), str(codon_positions_num), str(1)]), "\t", mean_diff, "\t", num_of_traitrelax_wins, "\t", num_of_relax_wins)
        return nullCompToCombo


def find_path(method, replicate, tbl, mu, pi0, taxa_num, codon_positions_num, k_value):
    replicate_regex = re.compile("Parsing.*?replicate_(\d*).*?fas", re.DOTALL | re.MULTILINE)
    data_dir = TR_simulation_study_output_dir + "/tbl_" + str(tbl) + "_mu_" + str(mu) + "_pi0_" + str(pi0) + "_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa_num) + "_taxa/" + str(codon_positions_num) + "_codons/k_" + str(k_value) + "/"
    if method == "RELAX":
        data_dir = R_simulation_study_output_dir + "/tbl_" + str(tbl) + "_mu_" + str(mu) + "_pi0_" + str(pi0) + "_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa_num) + "_taxa/" + str(codon_positions_num) + "_codons/k_" + str(k_value) + "/"
    for path in os.listdir(data_dir):
        if ".OU" in path:
            with open(data_dir + path, "r") as file:
                content = file.read()
            file_replicate = replicate_regex.search(content).group(1)
            if float(file_replicate) == replicate:
                return data_dir + path
    return "None"


def report_convergence_issue_examples(DF_dict, method="RELAX"):

    print("diff(alternative,null)\talternative_k\treplicate\tpath")
    tbl_to_mu = {1: 32, 4: 1, 4: 4, 4: 8, 8: 4, 16: 2, 32: 1}
    for tbl in tbl_to_mu:
        mu = tbl_to_mu[tbl]
        for taxa_num in taxa_num_options:
            for codon_positions_num in codon_positions_num_options:
                if (tbl, mu, 0.5, taxa_num, codon_positions_num, 1) in DF_dict:
                    DF = DF_dict[(tbl, mu, 0.5, taxa_num, codon_positions_num, 1)]
                    if method == "RELAX":
                        data = DF[["replicate", "null_logl", "alternative_logl", "alternative_k"]]
                        for index, row in data.iterrows():
                            diff = row["alternative_logl"] - row["null_logl"]
                            if diff > 1 and abs(row["alternative_k"]-1) < 0.01:
                                print(str(diff), "\t", str(row["alternative_k"]), "\t", str(row["replicate"]), "\t", find_path(method, row["replicate"], tbl, mu, 0.5, taxa_num, codon_positions_num, 1))
                    elif method == "TraitRELAX":
                        data = DF[["replicate", "null_sequence_logl", "alternative_sequence_logl", "alternative_k"]]
                        for index, row in data.iterrows():
                            diff = row["alternative_sequence_logl"] - row["null_sequence_logl"]
                            if diff > 1 and abs(row["alternative_k"] - 1) < 0.01:
                                print(str(diff), "\t", str(row["alternative_k"]), "\t", str(row["replicate"]), "\t",
                                      find_path(method, row["replicate"], tbl, mu, 0.5, taxa_num, codon_positions_num, 1))


if __name__ == '__main__':

    # extract the data from TraitRELAX executions
    TraitRELAXComboToDf = dict()
    EmpiricalTraitRELAXLRThresholds = dict()
    for tbl in tbl_options:
        for mu in mu_options:
            for pi0 in pi0_options:
                for taxa_num in taxa_num_options:
                    for codon_positions_num in codon_positions_num_options:
                        skip_combo = False
                        for k_value in k_values_options:

                            # declare the visited combo
                            input_dir = TR_simulation_study_output_dir + "/tbl_" + str(tbl) + "_mu_" + str(mu) + "_pi0_" + str(pi0) + "_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa_num) + "_taxa/" + str(codon_positions_num) + "_codons/k_" + str(k_value) + "/"
                            if not os.path.exists(input_dir):
                                if k_value == 1:
                                    skip_combo = True
                                continue
                            # extract the results
                            if not os.path.exists(input_dir + "full_analysis.csv") or not os.path.exists(input_dir + "clean_analysis.csv"):
                                df, full_df = extract_TraitRELAX_data(input_dir)
                                full_df.to_csv(input_dir + "full_analysis.csv")
                                df.to_csv(input_dir + "clean_analysis.csv")
                            else:
                                full_df = pd.read_csv(input_dir + "full_analysis.csv")
                                df = pd.read_csv(input_dir + "clean_analysis.csv")

                            # insert the dataframe into a dictionary
                            TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, codon_positions_num, k_value)] = df

    # extract the data fron RELAX + simulated character histories executions
    RELAXComboToDf = dict()
    EmpiricalRELAXLRThresholds = dict()
    for tbl in tbl_options:
        for mu in mu_options:
            for pi0 in pi0_options:
                for taxa_num in taxa_num_options:
                    for codon_positions_num in codon_positions_num_options:
                        skip_combo = False
                        for k_value in k_values_options:
                            # declare the visited combo
                            input_dir = R_simulation_study_output_dir + "/tbl_" + str(tbl) + "_mu_" + str(mu) + "_pi0_" + str(pi0) + "_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa_num) + "_taxa/" + str(codon_positions_num) + "_codons/k_" + str(k_value) + "/"
                            if not os.path.exists(input_dir):
                                if k_value == 1:
                                    skip_combo = True
                                continue

                            # extract the results
                            if not os.path.exists(input_dir + "full_analysis.csv") or not os.path.exists(input_dir + "clean_analysis.csv"):
                                df, full_df = extract_RELAX_data(input_dir, tbl, mu, pi0, taxa_num, codon_positions_num, k_value)
                                full_df.to_csv(input_dir + "full_analysis.csv")
                                df.to_csv(input_dir + "clean_analysis.csv")
                            else:
                                full_df = pd.read_csv(input_dir + "full_analysis.csv")
                                df = pd.read_csv(input_dir + "clean_analysis.csv")

                            # insert the dataframe into a dictionary
                            # print("processed combo: (tbl=", tbl, ", mu=", mu, ", pi0=", pi0, ", #taxa=", taxa_num, ", #pos=", codon_positions_num, ", k=", k_value, ")")
                            RELAXComboToDf[(tbl, mu, pi0, taxa_num, codon_positions_num, k_value)] = df

    print("*** Comparison of RELAX and TraitRELAX null likelihoods ****")
    nullCompToCombo = compare_null_likelihoods(TraitRELAXComboToDf, RELAXComboToDf)
    print("\n\n")
    print("*** Looking for RELAX examples with convergence issue ****")
    report_convergence_issue_examples(RELAXComboToDf, method="RELAX")
    print("\n\n")
    print("*** Looking for TraitRELAX examples with convergence issue ****")
    report_convergence_issue_examples(TraitRELAXComboToDf, method="TraitRELAX")
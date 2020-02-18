import os, re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

plt.switch_backend('agg')
from scipy.stats.distributions import chi2
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import axes3d, Axes3D  # <-- Note the capitalization!
import matplotlib.gridspec as gridspec

sns.set_style('whitegrid')
import math
from matplotlib.ticker import FormatStrFormatter, FuncFormatter, MaxNLocator

# hardcoded data
colors = ["lightgrey", "grey", "k"]
alternative_colors = ["mediumblue", "darkslateblue", "blueviolet"]
markers = ["o", "^", "s"]
MAX_RECORDS_TO_PROCESS = 50
tbl_options = [1, 4, 8, 16, 32]
mu_options = [1, 2, 4, 8, 16, 32]  # there is also 2 and 2 but not for FPR values
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
    return '%1.0f%%' % (100 * x)


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
    LR_cutoff = int(5 / 100 * len(null_LRs))
    if LR_cutoff == 0:
        LR_cutoff += 1
    LR_threshold = null_LRs[len(null_LRs) - LR_cutoff]  # any LR >= to LR_threshold will be considered significant
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
    alternative_character_logl_regex = re.compile(
        "Alternative model likelihood after optimization.*?Character Log likelihood\.*\:\s*(-\d*\.?\d*)",
        re.MULTILINE | re.DOTALL)
    alternative_sequence_logl_regex = re.compile(
        "Alternative model likelihood after optimization.*?Sequence Log likelihood\.*\:\s*(-\d*\.?\d*)",
        re.MULTILINE | re.DOTALL)  # bug in some output files - covert by actually computing the overall-char
    alternative_joint_logl_regex = re.compile(
        "Alternative model likelihood after optimization.*?Overall Log likelihood\.*\:\s*(-\d*\.?\d*)",
        re.MULTILINE | re.DOTALL)
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
        colnames.append("simulated_" + parameter)
        colnames.append("null_" + parameter)
        colnames.append("alternative_" + parameter)
    # manually computad parameters
    manual_parameters = ["omega0", "p0", "p1"]
    for parameter in manual_parameters:
        colnames.append("simulated_" + parameter)
        colnames.append("null_" + parameter)
        colnames.append("alternative_" + parameter)
    colnames = colnames + ["LR", "pvalue", "significant", "correctDirection", "duration(hours)"]
    df = pd.DataFrame(columns=colnames)
    clean_df = pd.DataFrame(columns=colnames)
    records_num = 0

    # fill in the dataframe
    for file in os.listdir(input_dir):
        if ".OU" in file and (
            (not "k_1/" in input_dir and records_num < MAX_RECORDS_TO_PROCESS) or "k_1/" in input_dir):
            # print("analysing ", input_dir, file)
            with open(input_dir + file, "r") as input_file:
                content = input_file.read()
                values = dict()

                # extract the values for the non-parametric columns
                for colname in non_parametric_regex_t_colname.keys():
                    try:
                        if "logl" in colname:
                            values["alternative_character_logl"] = float(
                                alternative_character_logl_regex.search(content).group(1))
                            values["alternative_joint_logl"] = float(
                                alternative_joint_logl_regex.search(content).group(1))
                            values["alternative_sequence_logl"] = values["alternative_joint_logl"] - values[
                                "alternative_character_logl"]
                    except:
                        print("unable to get log likelihood of alternative model in file ", input_dir + file)
                        exit(1)

                for colname in non_parametric_regex_t_colname.keys():
                    try:
                        regex = non_parametric_regex_t_colname[colname]
                        if not colname in values:
                            values[colname] = float(regex.search(content).group(1))
                    except:
                        print("failed to find " + colname + " match in " + input_dir + file + "-> setting as Nan")
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
                        print(
                            "failed to find simulated " + parameter + " match in " + input_dir + file + "-> setting as Nan")
                        values["simulated_" + parameter] = float('nan')
                        exit(1)

                    # try to extract the null
                    try:
                        regex_str = null_parameter_regex_str.replace("PARAMETER", parameter)
                        regex = re.compile(regex_str, re.MULTILINE | re.DOTALL)
                        values["null_" + parameter] = float(regex.search(content).group(1))
                    except:
                        print(
                            "failed to find null " + parameter + " match in " + input_dir + file + "-> setting as Nan")
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
                values["pvalue"] = chi2.sf(values["LR"], 1)  # 1 degree of freedom is the diff in num of parameters between alternative and null models

                # significant: True of False
                if values["pvalue"] < 0.05:
                    values["significant"] = 1
                else:
                    values["significant"] = 0

                # correctDirection: compare to simulated k -> if both < 1 or > 1 True, else False
                if (values["simulated_k"] < 1 and values["alternative_k"] < 1) or (
                        values["simulated_k"] > 1 and values["alternative_k"] > 1):
                    values["correctDirection"] = 1
                else:
                    values["correctDirection"] = 0

                # add the overall duration
                duration = 0
                for match in days_regex.finditer(content):
                    duration += float(match.group(1)) * 24
                for match in hours_regex.finditer(content):
                    duration += float(match.group(1))
                for match in minutes_regex.finditer(content):
                    duration += float(match.group(1)) * (1 / 60)
                for match in seconds_regex.finditer(content):
                    duration += float(match.group(1)) * (1 / 60) * (1 / 60)
                values["duration(hours)"] = duration

                # add the record to the dataframe
                add = True
                for key in list(values.keys()):
                    if math.isnan(values[key]):
                        add = False
                # add the record to the dataframe
                if add:
                    if (values['null_joint_logl'] - values['alternative_joint_logl']) > 0.5:
                        print("optimization bug expressed in " + input_dir + file)
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
    initial_logl_regex = re.compile("Computing intial log likelihood.*?Log likelihood\.*\:\s*(-\d*\.?\d*)",
                                    re.MULTILINE | re.DOTALL)
    null_logl_regex = re.compile("Fitting the null model.*?Log likelihood\.*\:\s*(-\d*\.?\d*)",
                                 re.MULTILINE | re.DOTALL)
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
        colnames.append("simulated_" + parameter)
        colnames.append("null_" + parameter)
        colnames.append("alternative_" + parameter)
    # manually computad parameters
    manual_parameters = ["omega0", "p0", "p1"]
    for parameter in manual_parameters:
        colnames.append("simulated_" + parameter)
        colnames.append("null_" + parameter)
        colnames.append("alternative_" + parameter)
    colnames = colnames + ["LR", "pvalue", "significant", "correctDirection", "duration(hours)"]
    # colnames.append("winning_starting_point") # keren - 3.7.19 - two sp execution
    df = pd.DataFrame(columns=colnames)
    clean_df = pd.DataFrame(columns=colnames)
    records_num = 0

    # fill in the dataframe
    for file in os.listdir(input_dir):
        if ".OU" in file and (
            (not "k_1/" in input_dir and records_num < MAX_RECORDS_TO_PROCESS) or "k_1/" in input_dir):

            with open(input_dir + file, "r") as input_file:
                content = input_file.read()
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

                non_parametric_regex_t_colname["alternative_logl"] = re.compile(alternative_logl_regex_str,
                                                                                re.MULTILINE | re.DOTALL)

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
                values["pvalue"] = chi2.sf(values["LR"], 1)  # 1 degree of freedom is the diff in num of parameters between alternative and null models

                # significant: True of False
                if values["pvalue"] < 0.05:
                    values["significant"] = 1
                else:
                    values["significant"] = 0

                # correctDirection: compare to simulated k -> if both < 1 or > 1 True, else False
                if (values["simulated_k"] < 1 and values["alternative_k"] < 1) or (
                        values["simulated_k"] > 1 and values["alternative_k"] > 1):
                    values["correctDirection"] = 1
                else:
                    values["correctDirection"] = 0

                # add the overall duration
                duration = 0
                for match in days_regex.finditer(content):
                    duration += float(match.group(1)) * 24
                for match in hours_regex.finditer(content):
                    duration += float(match.group(1))
                for match in minutes_regex.finditer(content):
                    duration += float(match.group(1)) * (1 / 60)
                for match in seconds_regex.finditer(content):
                    duration += float(match.group(1)) * (1 / 60) * (1 / 60)
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

# plot info: fixed to tbl=4, pi0=0.5, mu=8, taxa_num=32, codons=300
# x axis: simulated value of k
# y axis: %significant (over all alternative datasets)
def plot_power_vs_k_across_posnum(taxa_num, ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, title,
                                  empirical=True, add_legend=False, plot_labels=True, add_fpr_line=False):
    # specify the data to use
    tbl = 4
    mu = 8
    pi0 = 0.5

    ax.grid(False)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    # for each positions num, plot two curves: full for TraitRELAX, dashed for TraitRELAX with the true history
    for positions_num in codon_positions_num_options:
        # gather the data
        y_full_values = []
        for k in k_values_options:
            full_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
            if empirical:
                full_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
                full_power = full_df[full_df.LR >= full_threshold].shape[0] / full_df.shape[0]
                y_full_values.append(full_power)
            else:
                full_power = full_df[full_df.pvalue <= 0.05].shape[0] / full_df.shape[0]
                y_full_values.append(full_power)
        # plot the data in ax
        ax.plot(k_values_options, y_full_values, color=colors[codon_positions_num_options.index(positions_num)],
                marker=markers[taxa_num_options.index(taxa_num)], lw=2.5, label=str(positions_num) + "P",
                linestyle="dashed", markersize=16)
        if plot_labels:
            ax.set_xlabel(r"$k$", fontdict={'size': 30})
            ax.set_ylabel("Null rejected", fontdict={'size': 30})

    vals = [0, 0.2, 0.4, 0.6, 0.8, 1]
    ax.set_yticks(vals, ['{:,.0%}'.format(x) for x in vals])
    ax.yaxis.set_major_formatter(FuncFormatter(convertToPercent))
    ax.set_ylim(0, 1)
    ax.set_xticks(k_values_options)

    if add_fpr_line:
        ax.axhline(y=0.05, color='red', linestyle='dashed', label="5%")

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='upper center', prop={'size': 30}, frameon=False)


# plot info: fixed to tbl=4, pi0=0.5, mu=8, taxa_num=32, codons=300
# x axis: simulated value of k
# y axis: %significant (over all alternative datasets)
def plot_power_vs_k_across_taxanum(positions_num, ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, title,
                                   empirical=True, add_legend=False, plot_labels=True):
    # specify the data to use
    tbl = 4
    mu = 8
    pi0 = 0.5

    ax.grid(False)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    # for each positions num, plot two curves: full for TraitRELAX, dashed for TraitRELAX with the true history
    for taxa_num in taxa_num_options:
        # gather the data
        y_full_values = []
        for k in k_values_options:
            full_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
            if empirical:
                full_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
                full_power = full_df[full_df.LR >= full_threshold].shape[0] / full_df.shape[0]
                y_full_values.append(full_power)
            else:
                full_power = full_df[full_df.pvalue <= 0.05].shape[0] / full_df.shape[0]
                y_full_values.append(full_power)
        # plot the data in ax
        ax.plot(k_values_options, y_full_values, color=colors[codon_positions_num_options.index(positions_num)],
                marker=markers[taxa_num_options.index(taxa_num)], lw=2.5, label=str(taxa_num) + "S", linestyle="dashed",
                markersize=16)
        if plot_labels:
            ax.set_xlabel(r"$k$", fontdict={'size': 30})
            # ax.set_ylabel("Null rejected", fontdict={'size': 30})

    vals = [0, 0.2, 0.4, 0.6, 0.8, 1]
    ax.set_yticks(vals, ['{:,.0%}'.format(x) for x in vals])
    ax.yaxis.set_major_formatter(FuncFormatter(convertToPercent))
    ax.set_ylim(0, 1)
    ax.set_xticks(k_values_options)

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='upper center', prop={'size': 30}, frameon=False)


# plot info: fixed to tbl=4, pi0=0.5, taxa_num=32, codons=300, k=0.5
# x axis: mu
# y axis: power
# 2 curves: one for TraitRELAX standard execution (full), another for TraitRELAX given the true history (dashed)
def plot_power_vs_mu(ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, RELAXComboToDf,
                     EmpiricalRELAXLRThresholds, MPComboToDf, EmpiricalRELAXMPLRTThesholds, title):
    # gather the data
    tbl = 4
    pi0 = 0.5
    taxa_num = 32
    positions_num = 300
    k = 0.5
    sig_fraction_standard_execution = []
    sig_fraction_given_true_history = []
    sig_fraction_given_mp_history = []
    mu_local_options = [1, 2, 4, 8, 16]

    ax.set_xticks(mu_local_options)
    ax.set_ylim(0, 1)
    vals = [0, 0.2, 0.4, 0.6, 0.8, 1]
    ax.set_yticks(vals, ['{:,.0%}'.format(x) for x in vals])

    for mu in mu_local_options:
        standard_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        given_true_history_df = RELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        given_mp_history_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        LR_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
        sig_fraction_standard_execution.append(
            standard_df[standard_df.LR >= LR_threshold].shape[0] / standard_df.shape[0])
        LR_threshold = EmpiricalRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
        sig_fraction_given_true_history.append(
            given_true_history_df[given_true_history_df.LR >= LR_threshold].shape[0] / given_true_history_df.shape[0])
        LR_threshold = EmpiricalRELAXMPLRTThesholds[(tbl, mu, pi0, taxa_num, positions_num)]
        sig_fraction_given_mp_history.append(
            given_mp_history_df[given_mp_history_df.LR >= LR_threshold].shape[0] / given_mp_history_df.shape[0])

    ax.grid(False)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    # ax.plot(mu_local_options, sig_fraction_given_mp_history, linestyle=":", color='black', label="RELAX with maximum parsimony partition", lw=2.5)
    ax.plot(mu_local_options, sig_fraction_given_true_history, color=colors[codon_positions_num_options.index(300)],
            label="RELAX with true history", linestyle='dotted', marker=markers[taxa_num_options.index(32)], lw=2.5,
            markersize=16)
    ax.plot(mu_local_options, sig_fraction_standard_execution, label="TraitRELAX",
            color=colors[codon_positions_num_options.index(300)], linestyle='dashed',
            marker=markers[taxa_num_options.index(32)], lw=2.5, markersize=16)

    ax.set_xlabel(r"$\mu$", fontdict={'size': 30})
    # ax.set_ylabel("Null rejected", fontdict={'size': 30})

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, prop={'size': 30}, loc='best')  # , bbox_to_anchor=(0.5, -0.02))


# plot info: fixed to mu=0.5, pi0=0.5, taxa_num=32, codons=300, k=0.5
# x axis: tree length
# y axis: Power
# 2 curves: one for TraitRELAX standard execution (full), another for TraitRELAX given the true history (dashed)
def plot_power_vs_tbl(ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, RELAXComboToDf,
                      EmpiricalRELAXLRThresholds, MPComboToDf, EmpiricalRELAXMPLRTThesholds, title):
    # gather the data
    pi0 = 0.5
    taxa_num = 32
    positions_num = 300
    k = 0.5

    tbl_to_mu = {1: 32, 4: 8, 8: 4, 16: 2, 32: 1}
    ax.grid(False)
    sig_fraction_standard_execution = []
    sig_fraction_given_true_history = []
    sig_fraction_given_mp_history = []
    for tbl in tbl_options:
        mu = tbl_to_mu[tbl]
        standard_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        given_true_history_df = RELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        given_mp_history_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        try:
            LR_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            sig_fraction_standard_execution.append(
                standard_df[standard_df.LR >= LR_threshold].shape[0] / standard_df.shape[0])
        except:
            print("No standard TraitRELAX data is available for total branches lengths: ", tbl)
            sig_fraction_standard_execution.append(float('nan'))
        try:
            LR_threshold = EmpiricalRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            sig_fraction_given_true_history.append(
                given_true_history_df[given_true_history_df.LR >= LR_threshold].shape[0] / given_true_history_df.shape[
                    0])
        except:
            print("No RELAX data given true history is available for total branches lengths: ", tbl)
            sig_fraction_given_true_history.append(float('nan'))
        try:
            LR_threshold = EmpiricalRELAXMPLRTThesholds[(tbl, mu, pi0, taxa_num, positions_num)]
            sig_fraction_given_mp_history.append(
                given_mp_history_df[given_mp_history_df.LR >= LR_threshold].shape[0] / given_mp_history_df.shape[0])
        except Exception as e:
            print("No RELAX data given MP history is available for total branches lengths: ", tbl)
            sig_fraction_given_mp_history.append(float('nan'))

    # ax.plot(tbl_options, sig_fraction_given_mp_history, linestyle=":", color='black', label="RELAX with maximum parsimony partition", lw=2.5)
    ax.plot(tbl_options, sig_fraction_given_true_history, linestyle="--", color='black',
            label="RELAX with true history", lw=2.5)
    ax.plot(tbl_options, sig_fraction_standard_execution, label="TraitRELAX", linestyle='dashed', color='black', lw=2.5)

    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.set_xticks(tbl_options)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])

    ax.set_xlabel("Tree length", fontdict={'size': 30})
    ax.set_ylabel("Power", fontdict={'size': 30})

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, prop={'size': 30}, loc='best')


# plot info: fixed to tbl=4, pi0=0.5, taxa_num=32, codons=300, k=0.5
# x axis: mu
# y axis: empirical LR
# 2 curves: one for TraitRELAX standard execution (full), another for TraitRELAX given the true history (dashed)
def plot_empirical_LRs_vs_mu(ax, EmpiricalTraitRELAXLRThresholds, EmpiricalRELAXLRThresholds,
                             EmpiricalRELAXMPLRTThesholds, title):
    # gather the data
    pi0 = 0.5
    tbl = 4
    taxa_num = 32
    positions_num = 300

    ax.grid(False)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    LR_standard_execution = []
    LR_given_true_history = []
    LR_given_mp_history = []
    mu_local_options = [1, 2, 4, 8, 16]
    for mu in mu_local_options:
        # divide by 2 to achieve LR threshold instead of 2*LR threshold
        LR_standard_execution.append(EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)] / 2)
        try:
            # divide by 2 to achieve LR threshold instead of 2*LR threshold
            LR_given_true_history.append(EmpiricalRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)] / 2)
        except:
            print("no real data for combo: ", (tbl, mu, pi0, taxa_num, positions_num), " yet")
            LR_given_true_history.append(float("nan"))
            continue
        try:
            # divide by 2 to achieve LR threshold instead of 2*LR threshold
            LR_given_mp_history.append(EmpiricalRELAXMPLRTThesholds[(tbl, mu, pi0, taxa_num, positions_num)] / 2)
        except:
            print("no real data for combo: ", (tbl, mu, pi0, taxa_num, positions_num), " yet")
            LR_given_mp_history.append(float("nan"))
            continue
    ax.plot(mu_local_options, LR_standard_execution, label="TraitRELAX", linestyle='dashed', color='black', lw=2.5)
    # ax.plot(mu_local_options, LR_given_mp_history, linestyle=":", color='black', label="RELAX with maximum parsimony partition", lw=2.5)
    ax.plot(mu_local_options, LR_given_true_history, linestyle="--", color='black', label="RELAX with true history",
            lw=2.5)

    ax.set_yticks([1, 2, 3, 4, 5])
    ax.set_xticks(mu_local_options)

    ax.set_xlabel(r"$\mu$", fontdict={'size': 30})
    ax.set_ylabel("Computed empirical LR threshold", fontdict={'size': 30})
    ax.axhline(y=3.841, color='red', linestyle='dashed')


# plot info: fixed to: pi0=0.5, taxa_num=32, codons=300, k=0.5
# x axis: tbl
# y axis: empirical LR
# 2 curves: one for TraitRELAX standard execution (full), another for TraitRELAX given the true history (dashed)
def plot_empirical_LRs_vs_tbl(ax, EmpiricalTraitRELAXLRThresholds, EmpiricalRELAXLRThresholds,
                              EmpiricalRELAXMPLRTThesholds, title):
    # gather the data
    pi0 = 0.5
    taxa_num = 32
    positions_num = 300

    # plot b: vs. tbl
    tbl_to_mu = {1: 32, 4: 8, 8: 4, 16: 2, 32: 1}
    ax.grid(False)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    LR_standard_execution = []
    LR_given_true_history = []
    LR_given_mp_history = []
    for tbl in tbl_to_mu.keys():
        # divide by 2 to achieve LR threshold instead of 2*LR threshold
        LR_standard_execution.append(
            EmpiricalTraitRELAXLRThresholds[(tbl, tbl_to_mu[tbl], pi0, taxa_num, positions_num)] / 2)
        try:
            # divide by 2 to achieve LR threshold instead of 2*LR threshold
            LR_given_true_history.append(
                EmpiricalRELAXLRThresholds[(tbl, tbl_to_mu[tbl], pi0, taxa_num, positions_num)] / 2)
        except:
            print("no real data for combo: ", (tbl, tbl_to_mu[tbl], pi0, taxa_num, positions_num), " yet")
            LR_given_true_history.append(float("nan"))
            continue
        try:
            # divide by 2 to achieve LR threshold instead of 2*LR threshold
            LR_given_mp_history.append(
                EmpiricalRELAXMPLRTThesholds[(tbl, tbl_to_mu[tbl], pi0, taxa_num, positions_num)] / 2)
        except:
            print("no real data for combo: ", (tbl, tbl_to_mu[tbl], pi0, taxa_num, positions_num), " yet")
            LR_given_mp_history.append(float("nan"))
            continue

    ax.plot(list(tbl_to_mu.keys()), LR_standard_execution, label="TraitRELAX", linestyle='dashed', color='black',
            lw=2.5)
    # ax.plot(list(tbl_to_mu.keys()), LR_given_mp_history, linestyle=":", color='black', label="RELAX with maximum parsimony partition", lw=2.5)
    ax.plot(list(tbl_to_mu.keys()), LR_given_true_history, linestyle="--", color='black',
            label="RELAX with true history", lw=2.5)

    ax.set_xticks(tbl_options)
    ax.set_yticks([1, 2, 3, 4, 5])

    ax.set_xlabel("Tree length", fontdict={'size': 30})
    ax.set_ylabel("Computed empirical LR threshold", fontdict={'size': 30})
    ax.axhline(y=3.841, color='red', linestyle='dashed')


# power and FPR assessment
# figure 0 for results:
## plot 1: power based on empirical threshold: simulated k in x axis, significant(5) in y axis - for 32 taxa
## plot 2: power based on theoretical threshold: simulated k in x axis, significant(5) in y axis
## plot 3: power vs. mu
# figures for supp materials
# figure 1: same as plots 1 + 2 above for 16 and 64 taxa + power vs. tbl and empirical LR vs. tbl
def plot_power_and_FPR_assessment(TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, RELAXComboToDf,
                                  EmpiricalRELAXLRThresholds, MPComboToDf, EmpiricalMPRELAXLRThresholds,
                                  output_path_res, output_path_supp_1, output_path_supp_2):
    # figure 0 (for results)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=1, ncols=3, sharex="none", sharey="none", figsize=[3 * 8.2 + 2, 7.58])
    plot_power_vs_k_across_posnum(32, axis[0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "A\n",
                                  empirical=True, add_legend=True)
    # plot_power_vs_k_across_posnum(32, axis[1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "B\n", empirical=False, add_legend=True)
    plot_power_vs_k_across_taxanum(300, axis[1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "B\n",
                                   empirical=True, add_legend=True, plot_labels=True)
    plot_power_vs_mu(axis[2], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, RELAXComboToDf,
                     EmpiricalRELAXLRThresholds, MPComboToDf, EmpiricalMPRELAXLRThresholds, "C\n")
    # fig.text(0.95, 0.17, r"$k=0.5$", ha='center', fontdict={'size': 30})
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(output_path_res, bbox_inches='tight', transparent=True)
    plt.clf()

    # figure 1 (for supp materials)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=2, ncols=3, sharey='all', sharex='all', figsize=[3 * 8.2 + 2, 2 * 7.58])
    plot_power_vs_k_across_posnum(16, axis[0][0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "A\n",
                                  empirical=True, add_legend=True, plot_labels=False, add_fpr_line=True)
    plot_power_vs_k_across_posnum(32, axis[0][1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "B\n",
                                  empirical=True, add_legend=False, plot_labels=False, add_fpr_line=True)
    plot_power_vs_k_across_posnum(64, axis[0][2], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "C\n",
                                  empirical=True, add_legend=False, plot_labels=False, add_fpr_line=True)
    plot_power_vs_k_across_posnum(16, axis[1][0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "D\n",
                                  empirical=False, add_legend=False, plot_labels=False)
    plot_power_vs_k_across_posnum(32, axis[1][1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "E\n",
                                  empirical=False, add_legend=False, plot_labels=False)
    plot_power_vs_k_across_posnum(64, axis[1][2], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "F\n",
                                  empirical=False, add_legend=False, plot_labels=False)
    fig.text(0.53, -0.02, r"$k$", ha='center', fontdict={'size': 30})
    fig.text(-0.015, 0.75, 'Null rejected (5% FPR)', va='center', rotation='vertical', fontdict={'size': 30})
    fig.text(-0.015, 0.25, 'Null rejected (> 1.92 LR)', va='center', rotation='vertical', fontdict={'size': 30})
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(output_path_supp_1, bbox_inches='tight', transparent=True)
    plt.clf()

    # figure 2 (for supp materials)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=1, ncols=3, sharex="none", sharey="none", figsize=[3 * 8.2 + 2, 7.58])
    plot_power_vs_tbl(axis[0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, RELAXComboToDf,
                      EmpiricalRELAXLRThresholds, MPComboToDf, EmpiricalMPRELAXLRThresholds, "A\n")
    plot_empirical_LRs_vs_tbl(axis[1], EmpiricalTraitRELAXLRThresholds, EmpiricalRELAXLRThresholds,
                              EmpiricalMPRELAXLRThresholds, "B\n")
    plot_empirical_LRs_vs_mu(axis[2], EmpiricalTraitRELAXLRThresholds, EmpiricalRELAXLRThresholds,
                             EmpiricalMPRELAXLRThresholds, "C\n")
    fig.subplots_adjust()
    fig.tight_layout(pad=0.5)
    plt.savefig(output_path_supp_2, bbox_inches='tight', transparent=True)
    plt.clf()


#####################################################################################

# fixed parameters: tbl=4, mu=8, pi0=0.5
# plot assessment of inference with respect to k
# x axis: simulated value of k
# y axis: mean value inferred value of k
# add errorbars if not too noisy
# a panel with plot for each taxa num
# In each plot, 6 curves in 3 colors (color per positions number) - full for TraitRELAX, dashed for TraitRELAX with true history
# in green - scatter the simulated value of k (with legend)
def plot_accuracy_vs_k(ax, TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, title, add_legend=False,
                       include_reference=False):
    # declare fixed parameters
    tbl = 4
    mu = 8
    pi0 = 0.5
    taxa_num = 32
    positions_num = 300

    ax.grid(False)
    rel_error_sim_vs_standard = []
    rel_error_given_true_history_vs_standard = []
    rel_error_given_mp_history_vs_standard = []
    mu_local_options = [1, 2, 4, 8, 16]
    for k in k_values_options:
        standard_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        standard_df = standard_df.loc[(standard_df["alternative_k"] > 0)]
        # given_true_history_df = RELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        # given_true_history_df = given_true_history_df.loc[(given_true_history_df["alternative_k"] > 0)]
        given_mp_history_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        given_mp_history_df = given_mp_history_df.loc[(given_mp_history_df["alternative_k"] > 0)]
        rel_error_sim_vs_standard.append((abs(
            np.log(standard_df["alternative_k"] + 0.000001) - np.log(standard_df["simulated_k"] + 0.000001))).mean())
        # rel_error_given_true_history_vs_standard.append(abs(np.log(given_true_history_df["alternative_k"]+0.000001)-np.log(given_true_history_df["simulated_k"]+0.000001)).mean())
        rel_error_given_mp_history_vs_standard.append(abs(
            np.log(given_mp_history_df["alternative_k"] + 0.000001) - np.log(
                given_mp_history_df["simulated_k"] + 0.000001)).mean())

    ax.plot(k_values_options, rel_error_sim_vs_standard, color='black', linestyle="dashed", label="TraitRELAX", lw=2.5)
    ax.plot(k_values_options, rel_error_given_mp_history_vs_standard, color='black', linestyle=":",
            label="RELAX with maximum\nparsimony partition", lw=2.5)
    # if include_reference:
    #     ax.plot(k_values_options, rel_error_given_true_history_vs_standard, color='black', linestyle="--" , label="RELAX with true history")

    ax.set_xticks(k_values_options)
    ax.set_xlabel(r"$k$", fontdict={'size': 30})
    ax.set_ylabel("Mean error (" + r"$\^k$" + ")", fontdict={'size': 30})
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1])

    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])  # , 1.2, 1.4])

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best', prop={'size': 30}, frameon=False)


# plot relative error: dashed - compared to best possible (based on true history), full - based on simulated
# plot info: fixed to tbl=1, pi0=0.5, taxa_num=32, codons=600
# x axis: mu
# y axis: mean relative error of k = abs(simulated_k-inferred_k)/simulated_k (over all alternative datasets)
# 2 curves: one for TraitRELAX standard execution (full),
#           another for TraitRELAX given the true history (dashed),
def plot_accuracy_vs_mu(ax, TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, title, add_legend=False, add_ylabel=False,
                        include_reference=False):
    # gather the data
    tbl = 4
    pi0 = 0.5
    taxa_num = 32
    positions_num = 300

    ax.grid(False)
    rel_error_sim_vs_standard = []
    rel_error_given_true_history_vs_standard = []
    rel_error_given_mp_history_vs_standard = []
    mu_local_options = [1, 2, 4, 8, 16]
    for mu in mu_local_options:
        standard_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        standard_df = standard_df.loc[(standard_df["alternative_k"] > 0)]
        given_true_history_df = RELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        given_true_history_df = given_true_history_df.loc[(given_true_history_df["alternative_k"] > 0)]
        given_mp_history_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        given_mp_history_df = given_mp_history_df.loc[(given_mp_history_df["alternative_k"] > 0)]
        rel_error_sim_vs_standard.append((abs(
            np.log(standard_df["alternative_k"] + 0.000001) - np.log(standard_df["simulated_k"] + 0.000001))).mean())
        rel_error_given_true_history_vs_standard.append(abs(
            np.log(given_true_history_df["alternative_k"] + 0.000001) - np.log(
                given_true_history_df["simulated_k"] + 0.000001)).mean())
        rel_error_given_mp_history_vs_standard.append(abs(
            np.log(given_mp_history_df["alternative_k"] + 0.000001) - np.log(
                given_mp_history_df["simulated_k"] + 0.000001)).mean())

    ax.plot(mu_local_options, rel_error_sim_vs_standard, color='black', label="TraitRELAX", linestyle='dashed', lw=2.5)
    ax.plot(mu_local_options, rel_error_given_mp_history_vs_standard, color='black', linestyle=":",
            label="RELAX with maximum\nparsimony partition", lw=2.5)
    if include_reference:
        ax.plot(mu_local_options, rel_error_given_true_history_vs_standard, color='black', linestyle="--",
                label="RELAX with true history")

    ax.set_xticks(mu_local_options)
    ax.set_xlabel(r"$\mu$", fontdict={'size': 30})
    if add_ylabel:
        ax.set_ylabel("Mean error(" + r"$\^k$" + ")", fontdict={'size': 30})

    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best', prop={'size': 30}, frameon=False)


# plot relative error: dashed - compared to best possible (based on true history), full - based on simulated
# plot info: fixed to tbl=1, mu=0.5, taxa_num=32, codons=600
# x axis: mu
# y axis:relative error of k = abs(simulated_k-inferred_k)/simulated_k (over all alternative datasets)
# 2 curves: one for TraitRELAX standard execution (full),
#           another for TraitRELAX given the true history (dashed)
def plot_accuracy_vs_pi0(ax, TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, title, add_legend=False,
                         add_ylabel=False, include_reference=False):
    # gather the data
    tbl = 4
    mu = 8
    taxa_num = 32
    positions_num = 300

    ax.grid(False)
    rel_error_sim_vs_standard = []
    rel_error_given_true_history_vs_standard = []
    rel_error_given_mp_history_vs_standard = []
    for pi0 in pi0_options:
        standard_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        rel_error_sim_vs_standard.append((abs(
            np.log(standard_df["alternative_k"] + 0.000001) - np.log(standard_df["simulated_k"] + 0.000001))).mean())
        given_true_history_df = RELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        rel_error_given_true_history_vs_standard.append((abs(
            np.log(given_true_history_df["alternative_k"] + 0.000001) - np.log(
                given_true_history_df["simulated_k"] + 0.000001))).mean())
        given_mp_history_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        rel_error_given_mp_history_vs_standard.append((abs(
            np.log(given_mp_history_df["alternative_k"] + 0.000001) - np.log(
                given_mp_history_df["simulated_k"] + 0.000001))).mean())

    ax.plot(pi0_options, rel_error_sim_vs_standard, color='black', linestyle='dashed', label="TraitRELAX", lw=2.5)
    ax.plot(pi0_options, rel_error_given_mp_history_vs_standard, color='black', linestyle=":",
            label="RELAX with maximum\nparsimony partition", lw=2.5)
    if include_reference:
        ax.plot(pi0_options, rel_error_given_true_history_vs_standard, color='black', linestyle="--",
                label="RELAX with true history")

    ax.set_xticks(pi0_options)
    ax.set_xlabel(r"$\pi_0$", fontdict={'size': 30})
    if add_ylabel:
        ax.set_ylabel("Mean error (" + r"$\^k$" + ")", fontdict={'size': 30})

    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_xticks(pi0_options)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_ylabel("Mean error (" + r"$\^k$" + ")", fontdict={'size': 30, 'color': 'white'})

    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best', prop={'size': 30}, frameon=False)


# fixed parameters: tbl=4, mu=8, pi0=0.5
# plot assessment of inference with respect to k
# x axis: simulated value of k
# y axis: mean value inferred value of k
# add errorbars if not too noisy
# a panel with plot for each taxa num
# In each plot, 6 curves in 3 colors (color per positions number) - full for TraitRELAX, dashed for TraitRELAX with true history
# in green - scatter the simulated value of k (with legend)
def plot_accuracy_vs_k_across_positions_num(taxa_num, ax, TraitRELAXComboToDf, title, add_legend=False,
                                            plot_labels=True):
    # declare fixed parameters
    tbl = 4
    mu = 8
    pi0 = 0.5

    ax.grid(False)
    # for each positions num, plot two curves: full for TraitRELAX, dashed for TraitRELAX with the true history
    for positions_num in codon_positions_num_options:
        # gather the data
        y_full_values = []
        for k in k_values_options:
            full_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
            alternative_k_values = list(full_df["alternative_k"])
            for i in range(len(alternative_k_values)):
                if alternative_k_values[i] == 0:
                    alternative_k_values[i] = 0.000001
            k_error_values = (abs(np.log(np.array(alternative_k_values)) - np.log(full_df["simulated_k"])))
            y_full_values.append(np.mean(k_error_values))
        # plot the data in ax
        ax.plot(k_values_options, y_full_values, color=colors[codon_positions_num_options.index(positions_num)],
                label=str(positions_num) + "P", lw=2.5)

    if plot_labels:
        ax.set_xlabel(r"$k$", fontdict={'size': 30})
        ax.set_ylabel("Mean error (" + r"$\^k$" + ")", fontdict={'size': 30})

    ax.set_xticks(k_values_options)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.set_yticks([0, 0.4, 0.8, 1.2, 1.6, 2])

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best', prop={'size': 30}, frameon=False)


# fixed parameters: tbl=4, mu=8, pi0=0.5
# plot assessment of inference with respect to k
# x axis: simulated value of k
# y axis: mean value inferred value of k
# add errorbars if not too noisy
# a panel with plot for each taxa num
# In each plot, 6 curves in 3 colors (color per positions number) - full for TraitRELAX, dashed for TraitRELAX with true history
# in green - scatter the simulated value of k (with legend)
def plot_accuracy_vs_k_across_taxa_num(positions_num, ax, TraitRELAXComboToDf, title, add_legend=False,
                                       plot_labels=True):
    # declare fixed parameters
    tbl = 4
    mu = 8
    pi0 = 0.5

    ax.grid(False)
    # for each positions num, plot two curves: full for TraitRELAX, dashed for TraitRELAX with the true history
    for taxa_num in taxa_num_options:
        # gather the data
        y_full_values = []
        for k in k_values_options:
            full_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
            alternative_k_values = list(full_df["alternative_k"])
            for i in range(len(alternative_k_values)):
                if alternative_k_values[i] == 0:
                    alternative_k_values[i] = 0.000001
            k_error_values = (abs(np.log(np.array(alternative_k_values)) - np.log(full_df["simulated_k"])))
            y_full_values.append(np.mean(k_error_values))
        # plot the data in ax
        ax.plot(k_values_options, y_full_values, color=colors[taxa_num_options.index(taxa_num)],
                label=str(taxa_num) + "S", lw=2.5)

    if plot_labels:
        ax.set_xlabel(r"$k$", fontdict={'size': 30})
        # ax.set_ylabel("Mean error (" + r"$\^k$" + ")", fontdict={'size': 30})

    ax.set_xticks(k_values_options)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])  # , 1.2, 1.4])

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best', prop={'size': 30}, frameon=False)


def plot_2d_grid(ax, grid_data_path, title, dist):
    # collect data
    grid_data = pd.read_csv(grid_data_path, header=1)
    param1_values = grid_data.iloc[:, 0]
    param2_values = grid_data.iloc[:, 1]
    logl_values = grid_data.iloc[:, 2]

    # plot grid
    ax.grid(False)
    ax.plot_trisurf(param1_values, param2_values, logl_values, linewidth=1, antialiased=True,
                    cmap='viridis')  # rstride=1, cstride=1, cmap='viridis', edgecolor='none', alpha=.8)

    # scatter the ML point
    mlX, mlY = np.meshgrid([0.890321], [2.12883])
    mlZ = np.array([-15170.1873250066])
    ax.plot(mlX, mlY, mlZ, color='red', marker='o', label="MLE")

    # scatter the simulated point
    mlX, mlY = np.meshgrid([0.8], [2])
    mlZ = np.array([-15190.836444])
    ax.plot(mlX, mlY, mlZ, color='green', marker='o', label="True")

    ax.set_xticks([0, 0.5, 1, 1.5, 2])
    ax.set_xlabel("\n\n" + r"$k$", fontdict={'size': 30}, labelpad=8)
    ax.set_yticks([1, 2, 3, 4])
    ax.set_ylabel("\n\n" + r"$\omega_2$", fontdict={'size': 30}, labelpad=8)
    ax.set_zlabel("\n\nlog likelihood", fontdict={'size': 30}, labelpad=30)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, prop={'size': 30}, frameon=False, bbox_to_anchor=(0.9, 1))
    ax.dist = dist
    ax.tick_params(axis='z', which='major', pad=24)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')


def format_fn(tick_val, tick_pos):
    return "  " + str(tick_val) + "    "


def plot_inferred_vs_simulated_k_across_positions_num(taxa_num, ax, combo_to_df, title, plot_labels=False):
    # declare fixed parameters
    tbl = 4
    mu = 8
    pi0 = 0.5

    ax.grid(False)
    positions_to_combo = {(0.2, 150): 0, (0.2, 300): 0.3, (0.2, 600): 0.6, (0.5, 150): 1.2, (0.5, 300): 1.5,
                          (0.5, 600): 1.8, (0.8, 150): 2.4, (0.8, 300): 2.7, (0.8, 600): 3, (1, 150): 3.6,
                          (1, 300): 3.9, (1, 600): 4.2, (1.2, 150): 4.8, (1.2, 300): 5.1, (1.2, 600): 5.4,
                          (1.6, 150): 6, (1.6, 300): 6.3, (1.6, 600): 6.6, (2, 150): 7.2, (2, 300): 7.5, (2, 600): 7.8}

    # plot the data
    zorder_index = 1
    for k in k_values_options:
        boxplots = []
        # inferred_k_values = []
        # inferred_k_errors = []
        for positions_num in codon_positions_num_options:
            df = combo_to_df[(tbl, mu, pi0, taxa_num, positions_num, k)]
            # inferred_k_values.append(df["alternative_k"].mean())
            # inferred_k_errors.append(df["alternative_k"].std())
            filtered_df = df.loc[(df["alternative_k"] < 3)]
            boxplots.append(
                ax.boxplot(filtered_df["alternative_k"], widths=0.2, whis=[5, 95], showfliers=True, patch_artist=True,
                           positions=[positions_to_combo[(k, positions_num)]],
                           zorder=zorder_index))  # , color=colors[codon_positions_num_options.index[positions_num]]
            zorder_index += 1
        for boxplot in boxplots:
            color = colors[boxplots.index(boxplot)]
            patch = boxplot['boxes'][0]
            patch.set_facecolor(color)

        # ax.scatter(k_values_options, inferred_k_values, color=colors[codon_positions_num_options.index(positions_num)]) #, label=str(positions_num)+" positions")
        # ax.errorbar(k_values_options, inferred_k_values, yerr=inferred_k_errors, color=colors[codon_positions_num_options.index(positions_num)])
        x_for_scatter = [0, 0.3, 0.6, 1.2, 1.5, 1.8, 2.4, 2.7, 3, 3.6, 3.9, 4.2, 4.8, 5.1, 5.4, 6, 6.3, 6.6, 7.2, 7.5,
                         7.8]
        y_for_scatter = [0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.8, 0.8, 0.8, 1, 1, 1, 1.2, 1.2, 1.2, 1.6, 1.6, 1.6, 2, 2, 2]
        ax.scatter(x_for_scatter, y_for_scatter, label=" simulated values", color='green', zorder=zorder_index)

    if plot_labels:
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$^k$")

    ax.set_xticklabels(
        ["", "0.2", "", "", "0.5", "", "", "0.8", "", "", "1", "", "", "1.2", "", "", "1.6", "", "", "2", ""])
    ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')


# plot relative error: dashed - compared to best possible (based on true history), full - based on simulated
# plot info: fixed to tbl=1, pi0=0.5, taxa_num=32, codons=600
# x axis: mu
# y axis: mean relative error of k = abs(simulated_k-inferred_k)/simulated_k (over all alternative datasets)
# 2 curves: one for TraitRELAX standard execution (full),
#           another for TraitRELAX given the true history (dashed),
def plot_accuracy_vs_tbl(ax, TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, title, add_legend=False,
                         add_ylabel=False, include_reference=False):
    # gather the data
    tbl_to_mu = {1: 32, 4: 8, 8: 4, 16: 2, 32: 1}
    pi0 = 0.5
    taxa_num = 32
    positions_num = 300

    ax.grid(False)
    ax.grid(False)
    rel_error_sim_vs_standard = []
    rel_error_given_true_history_vs_standard = []
    rel_error_given_mp_history_vs_standard = []
    tbls = list(tbl_to_mu.keys())
    tbls.sort()
    for tbl in tbls:
        mu = tbl_to_mu[tbl]
        standard_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        standard_df = standard_df.loc[(standard_df["alternative_k"] > 0)]
        given_true_history_df = RELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        given_true_history_df = given_true_history_df.loc[(given_true_history_df["alternative_k"] > 0)]
        given_mp_history_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        given_mp_history_df = given_mp_history_df.loc[(given_mp_history_df["alternative_k"] > 0)]
        rel_error_sim_vs_standard.append((abs(
            np.log(standard_df["alternative_k"] + 0.000001) - np.log(standard_df["simulated_k"] + 0.000001))).mean())
        rel_error_given_true_history_vs_standard.append(abs(
            np.log(given_true_history_df["alternative_k"] + 0.000001) - np.log(
                given_true_history_df["simulated_k"] + 0.000001)).mean())
        rel_error_given_mp_history_vs_standard.append(abs(
            np.log(given_mp_history_df["alternative_k"] + 0.000001) - np.log(
                given_mp_history_df["simulated_k"] + 0.000001)).mean())

    ax.plot(tbls, rel_error_sim_vs_standard, color='black', linestyle='dashed', label="TraitRELAX", lw=2.5)
    ax.plot(tbls, rel_error_given_mp_history_vs_standard, color='black', linestyle=":",
            label="RELAX with maximum\nparsimony partition", lw=2.5)
    if include_reference:
        ax.plot(tbls, rel_error_given_true_history_vs_standard, color='black', linestyle="--",
                label="RELAX with true history", lw=2.5)

    ax.set_xticks(tbls)
    ax.set_xlabel("Tree length", fontdict={'size': 30})
    if add_ylabel:
        ax.set_ylabel("Mean error (" + r"$\^k$" + ")", fontdict={'size': 30})

    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best', prop={'size': 30}, frameon=False)


# fix: tbl=1, pi0=0.5, mu=0.5, taxa_num=32, codons=600
# a panel of 5 plots: one per k.
# each plot has 4 boxplots: k, omega0, omega1, omega2
# on each boxplot, scatter simulated value in green
def plot_k_to_omegas_accuracy(ax, comboToDf, title):
    # declare fixed parameters
    tbl = 4
    mu = 8
    pi0 = 0.5
    taxa_num = 32
    positions_num = 300
    fig_parameter_names = {"omega0": r"$\omega_0$",
                           "omega1": r"$\omega_1$",
                           "omega2": r"$\omega_2$",
                           "k": r"$k$"}

    # collect the  data
    colnames = ["simulated value of k"]
    for parameter in fig_parameter_names.keys():
        colnames.append(fig_parameter_names[parameter])
    df = pd.DataFrame(columns=colnames)
    true_values = []
    k_values_options = [0.2, 0.5, 1, 1.6, 2]
    for k in k_values_options:
        values = {"simulated value of k": k}
        source_df = comboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        # filter out rows with value of omega2 which is higher than 10
        source_df = source_df[(source_df.alternative_omega2 < 4) & (source_df.alternative_k < 4)]
        for index, row in source_df.iterrows():
            for parameter in fig_parameter_names.keys():
                parameter_value = row["alternative_" + parameter]
                values[fig_parameter_names[parameter]] = parameter_value
            df = df.append(values, ignore_index=True)
        for parameter in fig_parameter_names.keys():
            true_values.append(source_df["simulated_" + parameter].unique()[0])
    df.set_index("simulated value of k")
    # melt the dataframe by simulated_k
    melted_df = pd.melt(df, id_vars="simulated value of k")

    # plot the boxplots
    ax.grid(False)
    sns.boxplot(x="simulated value of k", y='value', hue='variable', data=melted_df, ax=ax, palette="Set3")
    ax.legend().set_title('')

    # scatter the simulated values
    x_for_scatter = [-0.3, -0.1, 0.1, 0.3, 0.7, 0.9, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.7, 2.9, 3.1, 3.3, 3.7, 3.9, 4.1,
                     4.3]
    ax.scatter(x=np.array(x_for_scatter), y=np.array(true_values), label="Simulated value", edgecolors=None,
               color="green", marker="s")

    ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4])
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.set_ylabel(r"$MLEs$", fontdict={'size': 30})
    ax.set_xlabel(r"$k$", fontdict={'size': 30})


# fix: tbl=1, pi0=0.5, mu=0.5, taxa_num=32, codons=600
# a panel of 5 plots: one per k.
# each plot has 4 boxplots: k, omega0, omega1, omega2
# on each boxplot, scatter simulated value in green
def plot_k_distribution_across_k_and_posnum(ax, comboToDf, title):
    # declare fixed parameters
    tbl = 4
    mu = 8
    pi0 = 0.5
    taxa_num = 32
    zorder_index = 1

    # collect the  data
    k_values_to_plot = [0.2, 1, 2]
    positions_to_combo = {(0.2, 150): 0, (0.2, 300): 0.3, (0.2, 600): 0.6, (1, 150): 1.2, (1, 300): 1.5, (1, 600): 1.8,
                          (2, 150): 2.4, (2, 300): 2.7, (2, 600): 3}
    for k in k_values_to_plot:
        boxplots = []
        for positions_num in codon_positions_num_options:
            df = comboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
            filetred_df = df.loc[(df["alternative_k"] < 3)]
            k_values = filetred_df["alternative_k"]
            boxplots.append(ax.boxplot(k_values, widths=0.2, whis=[5, 95], showfliers=True, patch_artist=True,
                                       positions=[positions_to_combo[(k, positions_num)]],
                                       zorder=zorder_index))  # , color=colors[codon_positions_num_options.index[positions_num]]
            zorder_index += 1
        for boxplot in boxplots:
            color = colors[boxplots.index(boxplot)]
            patch = boxplot['boxes'][0]
            patch.set_facecolor(color)

    x_for_scatter = [0, 0.3, 0.6, 1.2, 1.5, 1.8, 2.4, 2.7, 3]
    ax.scatter(x=x_for_scatter, y=[0.2, 0.2, 0.2, 1, 1, 1, 2, 2, 2], label="Simulated value", edgecolors=None,
               color="green", marker="s", zorder=zorder_index)

    # add legend
    # create custom legend
    custom_lines = []
    custom_names = []
    for positions_num in codon_positions_num_options:
        custom_lines.append(Line2D([0], [0], color=colors[codon_positions_num_options.index(positions_num)], lw=6))
        custom_names.append(str(positions_num) + "P")
    ax.legend(custom_lines, custom_names, loc='best', prop={'size': 30}, frameon=False)

    ax.grid(False)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.set_xticklabels(["", "0.2", "", "", "1", "", "", "2", ""])
    ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5, 3])
    ax.set_xlabel(r"$k$", fontdict={'size': 30})
    # ax.set_ylabel(r"$^k$", fontdict={'size': 30})


def plot_accuracy_analysis(TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, grid_data_path, res_output_path,
                           supp_output_path_1, supp_output_path_2, supp_output_path_3):
    # figure 0
    plt.grid(False)
    fig = plt.figure(figsize=[3 * 8.2 + 4, 7.58 + 2])
    gs = fig.add_gridspec(ncols=3, nrows=1)  # , width_ratios=[1, 1, 1.5])
    ax1 = fig.add_subplot(gs[0, 0])
    plot_k_distribution_across_k_and_posnum(ax1, TraitRELAXComboToDf, "A\n")
    ax2 = fig.add_subplot(gs[0, 1])
    plot_accuracy_vs_k(ax2, TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, "B\n", add_legend=True)
    ax3 = fig.add_subplot(gs[0, 2], projection='3d')
    plot_2d_grid(ax3, grid_data_path, "C\n", 14)
    fig.text(-0.01, 0.5, "" + r"$\^k$", va='center', rotation='vertical', fontdict={'size': 30})
    fig.tight_layout(pad=1)
    fig.subplots_adjust()
    plt.savefig(res_output_path, bbox_inches='tight', transparent=True)
    plt.clf()

    # # figure 0 (for results)
    # plt.grid(False)
    # fig = plt.figure(figsize=[3 * 8.2 + 2, 2 * 7.58]) #, constrained_layout=True)
    # gs = fig.add_gridspec(ncols=3, nrows=2)
    # ax1 = fig.add_subplot(gs[0,0])
    # plot_accuracy_vs_k(ax1, TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, "A\n")
    # ax2 = fig.add_subplot(gs[0,1])
    # plot_accuracy_vs_mu(ax2, TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, "B\n")
    # ax3 = fig.add_subplot(gs[0,2])
    # plot_accuracy_vs_pi0(ax3, TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, "C\n", add_legend=True)
    # ax4 = fig.add_subplot(gs[1,0])
    # plot_accuracy_vs_k_across_positions_num(32, ax4, TraitRELAXComboToDf, "D\n", add_legend=True, plot_labels=True)
    # ax5 = fig.add_subplot(gs[1,1])
    # plot_accuracy_vs_k_across_taxa_num(300, ax5, TraitRELAXComboToDf, "E\n", add_legend=True, plot_labels=True)
    # ax6 = fig.add_subplot(gs[1,2], projection='3d')
    # plot_2d_grid(ax6, grid_data_path, "F\n", 14)
    # fig.text(0.5, 0.95, r"$k=0.5$", ha='center', fontdict={'size': 30})
    # fig.text(0.95, 0.95, r"$k=0.5$", ha='center', fontdict={'size': 30})
    # fig.subplots_adjust()
    # fig.tight_layout(pad=2)
    # plt.savefig(res_output_path, bbox_inches='tight', transparent=True)
    # plt.clf()

    # figure 1 (for supp materials)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=1, ncols=3, sharey='none', sharex='none', figsize=[3 * 8.2 + 2, 1 * 7.58 + 2])
    plot_accuracy_vs_mu(axis[0], TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, "A\n", add_ylabel=True)
    plot_accuracy_vs_pi0(axis[1], TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, "B\n")
    plot_accuracy_vs_tbl(axis[2], TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, "C\n", add_legend=True,
                         add_ylabel=False)
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(supp_output_path_1, bbox_inches='tight', transparent=True)
    plt.clf()

    # figure 2 (for supp materials)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=2, ncols=3, sharey='none', sharex='none', figsize=[3 * 8.2 + 2, 2 * 7.58 + 2])
    plot_accuracy_vs_k_across_positions_num(16, axis[0][0], TraitRELAXComboToDf, "A\n", add_legend=True,
                                            plot_labels=False)
    plot_accuracy_vs_k_across_positions_num(32, axis[0][1], TraitRELAXComboToDf, "B\n", add_legend=False,
                                            plot_labels=False)
    plot_accuracy_vs_k_across_positions_num(64, axis[0][2], TraitRELAXComboToDf, "C\n", add_legend=False,
                                            plot_labels=False)
    plot_inferred_vs_simulated_k_across_positions_num(16, axis[1][0], TraitRELAXComboToDf, "D\n", plot_labels=False)
    plot_inferred_vs_simulated_k_across_positions_num(32, axis[1][1], TraitRELAXComboToDf, "E\n", plot_labels=False)
    plot_inferred_vs_simulated_k_across_positions_num(64, axis[1][2], TraitRELAXComboToDf, "F\n", plot_labels=False)
    fig.text(0.53, -0.02, r"$k$", ha='center', fontdict={'size': 30})
    fig.text(-0.02, 0.75, "Mean error (" + r"$\^k$" + ")", va='center', rotation='vertical', fontdict={'size': 30})
    fig.text(-0.02, 0.25, "Mean(" + r"$\^k$" + ")", va='center', rotation='vertical', fontdict={'size': 30})
    fig.subplots_adjust()
    fig.tight_layout(pad=1)
    plt.savefig(supp_output_path_2, bbox_inches='tight', transparent=True)
    plt.clf()

    # figure 3 (for supp materials)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=1, ncols=1, sharex="none", sharey="none", figsize=[1 * 8.2 + 2, 7.58])
    plot_k_to_omegas_accuracy(axis, TraitRELAXComboToDf, "")
    fig.subplots_adjust()
    fig.tight_layout(pad=1)
    plt.savefig(supp_output_path_3, bbox_inches='tight', transparent=True)
    plt.clf()


#####################################################################################

# fixed parameters: tbl=1, mu=0.5, pi0=0.5
# plot duration:
# 3 colors, 1 per taxa num
# x axis: positions number
# y axis: mean duration
def plot_duration(TraitRELAXComboToDf, output_path):
    # declare fixed parameters
    tbl = 4
    mu = 8
    pi0 = 0.5

    # gather the data
    duration_data_by_taxa = dict()
    for taxa_num in taxa_num_options:
        duration_data_by_taxa[taxa_num] = []
        for positions_num in codon_positions_num_options:
            dataframes = []
            for k in k_values_options:
                dataframes.append(TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)])
            joint_df = pd.concat(dataframes)
            duration_data_by_taxa[taxa_num].append(joint_df["duration(hours)"].mean())

    print("duration_data_by_taxa: ", duration_data_by_taxa)

    # plot a line for each taxa_num
    fig, axis = plt.subplots(nrows=1, ncols=1)
    plt.grid(False)
    for taxa_num in taxa_num_options:
        axis.plot(codon_positions_num_options, duration_data_by_taxa[taxa_num],
                  color=colors[taxa_num_options.index(taxa_num)], label=str(taxa_num) + " species")
    axis.set_xticks(codon_positions_num_options)
    axis.set_xlabel("Number of codon positions", fontdict={'size': 30})
    axis.set_ylabel("Mean duration (hours)", fontdict={'size': 30})
    handles, labels = axis.get_legend_handles_labels()
    lgd = fig.legend(handles, labels, loc='right', bbox_to_anchor=(1.2, 0.5), prop={'size': 30})
    # axis.set_yticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    plt.tight_layout(w_pad=0.5)
    fig.subplots_adjust(right=0.8)
    plt.savefig(output_path, bbox_extra_artists=(lgd,), bbox_inches='tight', figsize=[11, 7.58], transparent=True)
    plt.clf()


#####################################################################################


if __name__ == '__main__':

    matplotlib.rc('xtick', labelsize=22)
    matplotlib.rc('ytick', labelsize=22)
    matplotlib.rc('legend', fontsize=16)

    # create a directory of figures
    figures_dir = output_dir + "/figures/"
    if not os.path.exists(figures_dir):
        res = os.makedirs(figures_dir)

    # extract the data from TraitRELAX executions
    print("**** Processing TraitRELAX standard executions ****")
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
                            input_dir = TR_simulation_study_output_dir + "/tbl_" + str(tbl) + "_mu_" + str(
                                mu) + "_pi0_" + str(
                                pi0) + "_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(
                                taxa_num) + "_taxa/" + str(codon_positions_num) + "_codons/k_" + str(k_value) + "/"
                            if not os.path.exists(input_dir):
                                if k_value == 1:
                                    skip_combo = True
                                continue
                            # extract the results
                            if not os.path.exists(input_dir + "full_analysis.csv") or not os.path.exists(
                                            input_dir + "clean_analysis.csv"):
                                df, full_df = extract_TraitRELAX_data(input_dir)
                                full_df.to_csv(input_dir + "full_analysis.csv")
                                df.to_csv(input_dir + "clean_analysis.csv")
                            else:
                                full_df = pd.read_csv(input_dir + "full_analysis.csv")
                                df = pd.read_csv(input_dir + "clean_analysis.csv")

                            # insert the dataframe into a dictionary
                            # print("processed combo: (tbl=", tbl, ", mu=", mu, ", pi0=", pi0,
                            #       ", #taxa=", taxa_num, ", #pos=", codon_positions_num, ", k=", k_value, ")")
                            TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, codon_positions_num, k_value)] = df

                        # report theoretical and empirical results
                        if not skip_combo:
                            print("combo: (tbl=", tbl, ", mu=", mu, ", pi0=", pi0, ", taxa_num=",
                                  taxa_num, ", positions_num=", codon_positions_num, ")")
                            report_theoretical_test_result(TraitRELAXComboToDf, tbl, mu, pi0, taxa_num,
                                                           codon_positions_num)
                            EmpiricalTraitRELAXLRThresholds[
                                (tbl, mu, pi0, taxa_num, codon_positions_num)] = find_empirical_LR_threshold(
                                TraitRELAXComboToDf, tbl, mu, pi0, taxa_num, codon_positions_num)

    # perform individual analysis of traitrelax
    figures_subdir = output_dir + "/figures/TraitRELAX/"
    if not os.path.exists(figures_subdir):
        res = os.makedirs(figures_subdir)

    print("***********************************************\n\n")

    # extract the data fron RELAX + simulated character histories executions
    print("**** Processing TraitRELAX executions given the true history ****")
    # input_path = R_simulation_study_output_dir + "/RELAX_united.csv"
    # united_df = pd.read_csv(input_path)
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
                            input_dir = R_simulation_study_output_dir + "/tbl_" + str(tbl) + "_mu_" + str(
                                mu) + "_pi0_" + str(
                                pi0) + "_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(
                                taxa_num) + "_taxa/" + str(codon_positions_num) + "_codons/k_" + str(k_value) + "/"
                            # input_dir = R_simulation_study_output_dir + "/tbl_" + str(tbl) + "_mu_" + str(mu) + "_pi0_" + str(pi0) + "_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.6/" + str(taxa_num) + "_taxa/" + str(codon_positions_num) + "_codons/k_" + str(k_value) + "/"

                            if not os.path.exists(input_dir):
                                if k_value == 1:
                                    skip_combo = True
                                continue

                            # extract the results
                            if not os.path.exists(input_dir + "full_analysis.csv") or not os.path.exists(
                                            input_dir + "clean_analysis.csv"):
                                df, full_df = extract_RELAX_data(input_dir, tbl, mu, pi0, taxa_num, codon_positions_num,
                                                                 k_value)
                                full_df.to_csv(input_dir + "full_analysis.csv")
                                df.to_csv(input_dir + "clean_analysis.csv")
                            else:
                                full_df = pd.read_csv(input_dir + "full_analysis.csv")
                                df = pd.read_csv(input_dir + "clean_analysis.csv")

                            # df = united_df.loc[(united_df["tbl"] == tbl) & (united_df["mu"] == mu) & (united_df["pi0"] == pi0) & (united_df["taxa_num"] == taxa_num) & (united_df["positions_num"] == codon_positions_num) & (united_df["k"] == k_value)]

                            # # report the replicates for which the alternative_k is 10
                            # print("#taxa=", taxa_num, ", #positions=", codon_positions_num, ", k=", k_value)
                            # bads = df.loc[(df["alternative_k"] == 10)]
                            # print(list(bads["replicate"]))
                            # print("")
                            # print("df shape: ", df.shape[0])

                            # insert the dataframe into a dictionary
                            print("processed combo: (tbl=", tbl, ", mu=", mu, ", pi0=", pi0, ", #taxa=", taxa_num,
                                  ", #pos=", codon_positions_num, ", k=", k_value, ")")
                            RELAXComboToDf[(tbl, mu, pi0, taxa_num, codon_positions_num, k_value)] = df
                            kZeroInf = df.loc[(df["alternative_k"] == 0)]
                            print("replicates where estimated k=0: ", list(kZeroInf["replicate"]))

                        # report theoretical and empirical results
                        if not skip_combo:
                            print("combo: (tbl=", tbl, ", mu=", mu, ", pi0=", pi0, ", taxa_num=", taxa_num,
                                  ", positions_num=", codon_positions_num, ")")
                            report_theoretical_test_result(RELAXComboToDf, tbl, mu, pi0, taxa_num, codon_positions_num)
                            EmpiricalRELAXLRThresholds[
                                (tbl, mu, pi0, taxa_num, codon_positions_num)] = find_empirical_LR_threshold(
                                RELAXComboToDf, tbl, mu, pi0, taxa_num, codon_positions_num)

    # perform individual analysis of traitrelax
    figures_subdir = output_dir + "/figures/RELAXWithTrueHistory/"
    if not os.path.exists(figures_subdir):
        res = os.makedirs(figures_subdir)

    print("***********************************************\n\n")

    # extract the data from RELAX + simulated character histories executions
    print("**** Processing TraitRELAX executions given the maximum parsimony history ****")
    MPComboToDf = dict()
    EmpiricalMPRELAXLRThresholds = dict()
    for tbl in tbl_options:
        for mu in mu_options:
            for pi0 in pi0_options:
                for taxa_num in taxa_num_options:
                    for codon_positions_num in codon_positions_num_options:
                        skip_combo = False
                        for k_value in k_values_options:
                            # declare the visited combo
                            input_dir = R_MP_simulation_study_output_dir + "/tbl_" + str(tbl) + "_mu_" + str(
                                mu) + "_pi0_" + str(
                                pi0) + "_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(
                                taxa_num) + "_taxa/" + str(codon_positions_num) + "_codons/k_" + str(k_value) + "/"
                            if not os.path.exists(input_dir):
                                if k_value == 1:
                                    skip_combo = True
                                continue

                            # extract the results
                            if not os.path.exists(input_dir + "full_analysis.csv") or not os.path.exists(
                                            input_dir + "clean_analysis.csv"):
                                # df, full_df = extract_RELAX_data(input_dir, tbl, mu, pi0, taxa_num, codon_positions_num, k_value)
                                df, full_df = extract_RELAX_data(input_dir, tbl, mu, pi0, taxa_num, codon_positions_num,
                                                                 k_value)
                                full_df.to_csv(input_dir + "full_analysis.csv")
                                df.to_csv(input_dir + "clean_analysis.csv")
                            else:
                                full_df = pd.read_csv(input_dir + "full_analysis.csv")
                                df = pd.read_csv(input_dir + "clean_analysis.csv")

                            # insert the dataframe into a dictionary
                            print("processed combo: (tbl=", tbl, ", mu=", mu, ", pi0=", pi0,
                                  ", #taxa=", taxa_num, ", #pos=", codon_positions_num, ", k=", k_value, ")")
                            MPComboToDf[(tbl, mu, pi0, taxa_num, codon_positions_num, k_value)] = df

                        # report theoretical and empirical results
                        if not skip_combo:
                            print("combo: (tbl=", tbl, ", mu=", mu, ", pi0=", pi0, ", taxa_num=",
                                  taxa_num, ", positions_num=", codon_positions_num, ")")
                            report_theoretical_test_result(MPComboToDf, tbl, mu, pi0, taxa_num, codon_positions_num)
                            EmpiricalMPRELAXLRThresholds[
                                (tbl, mu, pi0, taxa_num, codon_positions_num)] = find_empirical_LR_threshold(
                                MPComboToDf, tbl, mu, pi0, taxa_num, codon_positions_num)

    # perform individual analysis of traitrelax
    figures_subdir = output_dir + "/figures/RELAXWithMPHistory/"
    if not os.path.exists(figures_subdir):
        res = os.makedirs(figures_subdir)

    print("***********************************************\n\n")

    plot_power_and_FPR_assessment(TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, RELAXComboToDf,
                                  EmpiricalRELAXLRThresholds, MPComboToDf, EmpiricalMPRELAXLRThresholds,
                                  figures_dir + "power_and_fpr_assessment_for_results.svg",
                                  figures_dir + "power_and_fpr_assessment_for_supp_1.svg",
                                  figures_dir + "power_and_fpr_assessment_for_supp_2.svg")

    # plot_accuracy_analysis(TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, grid_data_path,
    #                        figures_dir + "accuracy_assessment_for_results.svg",
    #                        figures_dir + "accuracy_assessment_for_supp_1.svg",
    #                        figures_dir + "accuracy_assessment_for_supp_2.svg",
    #                        figures_dir + "accuracy_assessment_for_supp_3.svg")
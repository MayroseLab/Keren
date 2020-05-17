from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import os, re, matplotlib, math
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats.distributions import chi2
from matplotlib.lines import Line2D
sns.set_style('whitegrid') #make the plots appear on white grid
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
import matplotlib.cm as cm #make the plots look pretty
sns.set()

# hardcoded data
colors = ["lightgrey", "grey", "k"]
alternative_colors = ["mediumblue", "darkslateblue", "blueviolet"]
pretty_colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']
markers = ["o", "^", "s"]
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
R_HyPhy_MP_simulation_study_output_dir = simulation_study_output_dir + "/RELAX_HyPhy"
output_dir = simulation_study_output_dir
grid_data_path = output_dir + "/grid_data.csv"


#####################################################################################

def format_fn(tick_val, tick_pos):
    return "  " + str(tick_val) + "    "


def rel_error(simulated_values, inferred_values):
    return abs((simulated_values-inferred_values)/simulated_values)


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
                    power = alternative_df.loc[(alternative_df["LR"] >= LR_threshold) & (alternative_df["correctDirection"] == 1) ].shape[0] / alternative_df.shape[0]
                    print("For k = " + str(k) + ": power = " + str(power))
            except:
                continue
    print("\n")
    return LR_threshold


#####################################################################################

def extract_hyphy_data(input_dir):

    colname_to_regex = {"replicate": re.compile("replicate_(\d*)", re.MULTILINE | re.DOTALL),
                                      "null_logl": re.compile("Fitting the null \(K \:= 1\) model.*?Log\(L\)\s*=\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "alternative_logl": re.compile("Fitting the alternative model.*?Log\(L\)\s*=\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "pvalue" : re.compile("Likelihood ratio test\s*\*\*p\s*=\s*(\d*\.?\d*)\*\*", re.MULTILINE | re.DOTALL),
                                      "simulated_omega0" : re.compile("_p_(\d*\.?\d*)_omega1_(\d*\.?\d*)_", re.MULTILINE | re.DOTALL),
                                      "simulated_omega1" : re.compile("_omega1_(\d*\.?\d*)_", re.MULTILINE | re.DOTALL),
                                      "simulated_omega2": re.compile("_omega2_(\d*\.?\d*)_", re.MULTILINE | re.DOTALL),
                                      "simulated_p0" : re.compile("_theta1_(\d*\.?\d*)_", re.MULTILINE | re.DOTALL),
                                      "simulated_p1" : re.compile("_theta1_(\d*\.?\d*)_theta2_(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "simulated_k" : re.compile("\/k_(\d*\.?\d*)\/", re.MULTILINE | re.DOTALL),
                                      "null_omega0": re.compile("Fitting the null(.*?\|){12}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "null_omega1": re.compile("Fitting the null(.*?\|){17}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "null_omega2": re.compile("Fitting the null(.*?\|){22}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "null_p0": re.compile("Fitting the null(.*?\|){13}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "null_p1": re.compile("Fitting the null(.*?\|){18}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "alternative_omega0": re.compile("Fitting the alternative model.*?The following rate distribution was inferred for \*\*reference\*\* branches(.*?\|){12}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "alternative_omega1": re.compile("Fitting the alternative model.*?The following rate distribution was inferred for \*\*reference\*\* branches(.*?\|){17}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "alternative_omega2": re.compile("Fitting the alternative model.*?The following rate distribution was inferred for \*\*reference\*\* branches(.*?\|){22}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "alternative_p0": re.compile("Fitting the alternative model.*?The following rate distribution was inferred for \*\*reference\*\* branches(.*?\|){13}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "alternative_p1": re.compile("Fitting the alternative model.*?The following rate distribution was inferred for \*\*reference\*\* branches(.*?\|){18}\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "alternative_k": re.compile("Relaxation\/intensification parameter\s*\(K\)\s*=\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                                      "job_id" : re.compile("(\d*)\.power8", re.MULTILINE | re.DOTALL)}
    colnames = list(colname_to_regex.keys()) + ["LR", "significant", "correctDirection", "job_id"]
    df = pd.DataFrame(columns=colnames)

    for path in os.listdir(input_dir):
        if not "power8" in path:
            continue
        record = dict()
        record["job_id"] = colname_to_regex["job_id"].search(path).group(1)
        with open(input_dir+path, "r") as infile:
            content = infile.read()
        for colname in list(colname_to_regex.keys()):
            if colname == "simulated_omega0":
                match = colname_to_regex[colname].search(content)
                p = match.group(1)
                omega1 = match.group(2)
                record[colname] = float(p) * float(omega1)
            elif colname == "simulated_p1":
                match = colname_to_regex[colname].search(content)
                theta1 = match.group(1)
                theta2 = match.group(2)
                record[colname] = (1-float(theta1))*float(theta2)
            elif not "simulated" in colname and "omega" in colname:
                record[colname] = float(colname_to_regex[colname].search(content).group(2))
            elif not "simulated" in colname and ("p0" in colname or "p1" in colname):
                record[colname] = float(colname_to_regex[colname].search(content).group(2)) / 100
            elif colname != "job_id":
                # print("colname: ", colname, "\nregex: ", colname_to_regex[colname].pattern, "\npath: ", input_dir+path)
                record[colname] = float(colname_to_regex[colname].search(content).group(1))

        record["LR"] = 2 * (record["alternative_logl"] - record["null_logl"])
        pvalue = doLRT(record["null_logl"], record["alternative_logl"])
        if abs(pvalue - record["pvalue"]) > 0.05:
            # record["pvalue"] = pvalue
            # print("error! computed pvalue is different than the one reported by hyhphy!")
            # print("input path: ", input_dir + path)
            print("replicate: ", record["replicate"])
            # print("null logl: ", record["null_logl"])
            # print("alternative logl: ", record["alternative_logl"])
            print("computed pvalue: ", pvalue)
            print("reported pvalue: ", record["pvalue"], "\n")
            # exit(1)
        if pvalue < 0.05:
            record["significant"] = 1
        else:
            record["significant"] = 0
        if (record["simulated_k"] < 1 and record["alternative_k"] < 1) or (record["simulated_k"] > 1 and record["alternative_k"] > 1) or (record["simulated_k"] == 1 and record["alternative_k"] == 1):
            record["correctDirection"] = 1
        else:
            record["correctDirection"] = 0

        df = df.append(record, ignore_index=True)
    return df


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
            not "k_1/" in input_dir or "k_1/" in input_dir):
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
                        values["simulated_k"] > 1 and values["alternative_k"] > 1) or (values["simulated_k"] == 1 and values["alternative_k"] == 1):
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
    null_parameter_regex_str = "Fitting the null model.*?iterative.*?model 2.*?.RELAX.PARAMETER\.*\:\s*(\d*\.?\d*)"
    alternative_parameter_regex_str = "Fitting the alternative model.*?iterative.*?model 2.*?.RELAX.PARAMETER\.*\:\s*(\d*\.?\d*)"
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
            not "k_1/" in input_dir or "k_1/" in input_dir):

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
                        values["simulated_k"] > 1 and values["alternative_k"] > 1) or (values["simulated_k"] == 1 and values["alternative_k"] == 1):
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

# get the percentage of null reject in the write and wrong direction
def get_div_for_null_rejected(df, empricial_threshold, theoretical=0, out_of_all=False):


    if theoretical == 1:
        significants_df = df.loc[(df["significant"] == 1)]
    else:
        significants_df = df.loc[(df["LR"] >= empricial_threshold)]
    num_of_significants = significants_df.shape[0]
    if num_of_significants == 0:
        return 0, 0

    num_of_replicates = df.shape[0]
    correct_direction = significants_df.loc[(significants_df["correctDirection"] == 1)]
    num_of_correct_direction = correct_direction.shape[0]
    fraction_correct_direction = num_of_correct_direction / num_of_significants
    fraction_wrong_direction = 1 - fraction_correct_direction
    if out_of_all:
        fraction_correct_direction = num_of_correct_direction / num_of_replicates
        fraction_wrong_direction = (num_of_significants-num_of_correct_direction) / num_of_replicates

    return fraction_correct_direction, fraction_wrong_direction


# plot info: fixed to tbl=4, pi0=0.5, mu=8, taxa_num=32, codons=300
# x axis: simulated value of k
# y axis: %significant (over all alternative datasets)
def plot_power_vs_k_across_posnum(taxa_num, ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, title,
                                  empirical=True, add_legend=False, plot_labels=True, add_fpr_line=False, linestyle='solid'):
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
                full_power = full_df.loc[(full_df["LR"] >= full_threshold)].shape[0] / full_df.shape[0]
                y_full_values.append(full_power)
            else:
                full_power = full_df.loc[(full_df["pvalue"] <= 0.05)].shape[0] / full_df.shape[0]
                y_full_values.append(full_power)
        # plot the data in ax
        ax.plot(k_values_options, y_full_values, color=colors[codon_positions_num_options.index(positions_num)],
                marker=markers[taxa_num_options.index(taxa_num)], lw=2.5, label=str(positions_num) + " sites",
                linestyle=linestyle, markersize=16)
        if plot_labels:
            ax.set_xlabel(r"$k$", fontdict={'size': 30})
            ax.set_ylabel("Null rejected", fontdict={'size': 30})

    if add_fpr_line:
        ax.axhline(y=0.05, color='red', linestyle='solid', label="Required FPR")

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='upper center', prop={'size': 30}, frameon=False)

    vals = [0, 0.05, 0.2, 0.4, 0.6, 0.8, 1]
    ax.set_yticks(vals, ['{:,.0%}'.format(x) for x in vals])
    ax.yaxis.set_major_formatter(FuncFormatter(convertToPercent))
    ax.set_ylim(0, 1)
    ax.set_xticks(k_values_options)

# plot info: fixed to tbl=4, pi0=0.5, mu=8, taxa_num=32, codons=300
# x axis: simulated value of k
# y axis: %significant (over all alternative datasets)
def plot_power_vs_k_across_taxanum(positions_num, ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, title,
                                   empirical=True, add_legend=False, plot_labels=True, linestyle='solid'):
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
                full_power = full_df.loc[(full_df["LR"] >= full_threshold)].shape[0] / full_df.shape[0]
                y_full_values.append(full_power)
            else:
                full_power = full_df.loc[(full_df["pvalue"] <= 0.05)].shape[0] / full_df.shape[0]
                y_full_values.append(full_power)
        # plot the data in ax
        ax.plot(k_values_options, y_full_values, color=colors[codon_positions_num_options.index(positions_num)],
                marker=markers[taxa_num_options.index(taxa_num)], lw=2.5, label=str(taxa_num) + " taxa", linestyle=linestyle,
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
                     EmpiricalRELAXLRThresholds, MPComboToDf, EmpiricalMPThresholds, title):
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
    ax.yaxis.set_major_formatter(FuncFormatter(convertToPercent))

    for mu in mu_local_options:
        standard_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        given_true_history_df = RELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        LR_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
        sig_fraction_standard_execution.append(
            standard_df.loc[(standard_df["LR"] >= LR_threshold)].shape[0] / standard_df.shape[0])
        LR_threshold = EmpiricalRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
        sig_fraction_given_true_history.append(
            given_true_history_df.loc[(given_true_history_df["LR"] >= LR_threshold)].shape[0] / given_true_history_df.shape[0])
        given_mp_history_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        #LR_threshold = EmpiricalMPThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
        sig_fraction_given_mp_history.append(given_mp_history_df.loc[(given_mp_history_df["pvalue"] <= 0.05)].shape[0] / given_mp_history_df.shape[0])


    ax.grid(False)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.plot(mu_local_options, sig_fraction_standard_execution, label="TraitRELAX",
            color=colors[codon_positions_num_options.index(300)], linestyle='solid',
            marker=markers[taxa_num_options.index(32)], lw=2.5, markersize=16)
    ax.plot(mu_local_options, sig_fraction_given_true_history, color=colors[codon_positions_num_options.index(300)],
            label="TraitRELAX+TH", linestyle='dotted', marker=markers[taxa_num_options.index(32)], lw=2.5,
            markersize=16)
    ax.plot(mu_local_options, sig_fraction_given_mp_history, label="RELAX",
            color=colors[codon_positions_num_options.index(300)], linestyle='dashed',
            marker=markers[taxa_num_options.index(32)], lw=2.5, markersize=16)


    ax.set_xlabel(r"$\mu$", fontdict={'size': 30})
    # ax.set_ylabel("Null rejected", fontdict={'size': 30})

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, prop={'size': 30}, loc='upper right', frameon=False, bbox_to_anchor=(1.1, 0.88))


# plot info: fixed to mu=0.5, pi0=0.5, taxa_num=32, codons=300, k=0.5
# x axis: tree length
# y axis: Power
# 2 curves: one for TraitRELAX standard execution (full), another for TraitRELAX given the true history (dashed)
def plot_power_vs_tbl(ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, RELAXComboToDf,
                      EmpiricalRELAXLRThresholds, title):
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
        try:
            LR_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            sig_fraction_standard_execution.append(
                standard_df.loc[(standard_df["LR"] >= LR_threshold) & (standard_df["correctDirection"] == 1)].shape[0] / standard_df.shape[0])
        except:
            print("No standard TraitRELAX data is available for total branches lengths: ", tbl)
            sig_fraction_standard_execution.append(float('nan'))
        try:
            LR_threshold = EmpiricalRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            sig_fraction_given_true_history.append(
                given_true_history_df.loc[(given_true_history_df["LR"] >= LR_threshold) & (given_true_history_df["correctDirection"] == 1)].shape[0] / given_true_history_df.shape[
                    0])
        except:
            print("No RELAX data given true history is available for total branches lengths: ", tbl)
            sig_fraction_given_true_history.append(float('nan'))

    ax.plot(tbl_options, sig_fraction_given_true_history, linestyle='dotted', color='black',
            label="TraitRELAX with true history", lw=2.5)
    ax.plot(tbl_options, sig_fraction_standard_execution, label="TraitRELAX", linestyle='solid', color='black', lw=2.5)

    ax.plot([],[], linestyle='dashed', color='red', label="Theoretical threshold")

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
def plot_empirical_LRs_vs_mu(ax, EmpiricalTraitRELAXLRThresholds, EmpiricalRELAXLRThresholds, title):
    # gather the data
    pi0 = 0.5
    tbl = 4
    taxa_num = 32
    positions_num = 300

    ax.grid(False)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    LR_standard_execution = []
    LR_given_true_history = []
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
    ax.plot(mu_local_options, LR_standard_execution, label="TraitRELAX", linestyle='solid', color='black', lw=2.5)
    #ax.plot(mu_local_options, LR_given_true_history, linestyle='dotted', color='black', label="TraitRELAX with true history",
            # lw=2.5)

    ax.set_yticks([1, 2, 3, 4, 5])
    ax.set_xticks(mu_local_options)

    ax.set_xlabel(r"$\mu$", fontdict={'size': 30})
    ax.set_ylabel("Computed empirical LR threshold", fontdict={'size': 30})
    ax.axhline(y=3.841/2, color='red', linestyle='solid', label="Theoretical threshold")


# plot info: fixed to: pi0=0.5, taxa_num=32, codons=300, k=0.5
# x axis: tbl
# y axis: empirical LR
# 2 curves: one for TraitRELAX standard execution (full), another for TraitRELAX given the true history (dashed)
def plot_empirical_LRs_vs_tbl(ax, EmpiricalTraitRELAXLRThresholds, EmpiricalRELAXLRThresholds, title):
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

    ax.plot(list(tbl_to_mu.keys()), LR_standard_execution, label="TraitRELAX", linestyle='solid', color='black',
            lw=2.5)
    # ax.plot(list(tbl_to_mu.keys()), LR_given_mp_history, linestyle=":", color='black', label="RELAX with maximum parsimony partition", lw=2.5)
    ax.plot(list(tbl_to_mu.keys()), LR_given_true_history, linestyle='dotted', color='black',
            label="TraitRELAX with true history", lw=2.5)

    ax.set_xticks(tbl_options)
    ax.set_yticks([1, 2, 3, 4, 5])

    ax.set_xlabel("Tree length", fontdict={'size': 30})
    ax.set_ylabel("Computed empirical LR threshold", fontdict={'size': 30})
    ax.axhline(y=3.841, color='red', linestyle='solid')


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
    fig, axis = plt.subplots(nrows=1, ncols=3, sharex="none", sharey="none", figsize=[3 * 8.2 + 2, 7.58], frameon=True)
    plot_power_vs_k_across_posnum(32, axis[0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "a\n",
                                  empirical=True, add_legend=True)
    plot_power_vs_k_across_taxanum(300, axis[1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "b\n",
                                   empirical=True, add_legend=True, plot_labels=True)
    plot_power_vs_mu(axis[2], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, RELAXComboToDf, EmpiricalRELAXLRThresholds, MPComboToDf, EmpiricalMPRELAXLRThresholds, "c\n")
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(output_path_res, bbox_inches='tight') #, transparent=True)
    plt.clf()

    # figure 1 (for supp materials)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=2, ncols=3, sharey='all', sharex='all', figsize=[3 * 8.2 + 2, 2 * 7.58], frameon=True)
    plot_power_vs_k_across_posnum(16, axis[0][0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "a\n",
                                  empirical=True, add_legend=True, plot_labels=False, add_fpr_line=True)
    plot_power_vs_k_across_posnum(32, axis[0][1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "b\n",
                                  empirical=True, add_legend=False, plot_labels=False, add_fpr_line=True)
    plot_power_vs_k_across_posnum(64, axis[0][2], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "c\n",
                                  empirical=True, add_legend=False, plot_labels=False, add_fpr_line=True)
    plot_power_vs_k_across_posnum(16, axis[1][0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "d\n",
                                  empirical=False, add_legend=False, plot_labels=False)
    plot_power_vs_k_across_posnum(32, axis[1][1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "e\n",
                                  empirical=False, add_legend=False, plot_labels=False)
    plot_power_vs_k_across_posnum(64, axis[1][2], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "f\n",
                                  empirical=False, add_legend=False, plot_labels=False)
    fig.text(0.53, -0.02, r"$k$", ha='center', fontdict={'size': 30})
    fig.text(-0.015, 0.75, 'Null rejected (5% FPR)', va='center', rotation='vertical', fontdict={'size': 30})
    fig.text(-0.015, 0.25, 'Null rejected (> 1.92 LR)', va='center', rotation='vertical', fontdict={'size': 30})
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(output_path_supp_1, bbox_inches='tight') #, transparent=True)
    plt.clf()

    # figure 2 (for supp materials)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=1, ncols=3, sharex="none", sharey="none", figsize=[3 * 8.2 + 2, 7.58], frameon=True)
    plot_power_vs_tbl(axis[0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, RELAXComboToDf,
                      EmpiricalRELAXLRThresholds, "a\n")
    plot_empirical_LRs_vs_tbl(axis[1], EmpiricalTraitRELAXLRThresholds, EmpiricalRELAXLRThresholds, "b\n")
    plot_empirical_LRs_vs_mu(axis[2], EmpiricalTraitRELAXLRThresholds, EmpiricalRELAXLRThresholds, "c\n")
    fig.subplots_adjust()
    fig.tight_layout(pad=0.5)
    plt.savefig(output_path_supp_2, bbox_inches='tight') #, transparent=True)
    plt.clf()


#####################################################################################

# fixed parameters: tbl=4, mu=8, pi0=0.5
# plot assessment of inference with respect to k
# x axis: simulated value of k
# y axis: mean value inferred value of k
# add errorbars if not too noisy
# a panel with plot for each taxa num
# In each plot, 6 curves in 3 colors (color per positions number) - full for TraitRELAX, dashed for TraitTraitRELAX with true history
# in green - scatter the simulated value of k (with legend)

def global_accuracy(df, debug=False):

    p0 = df["simulated_p0"][0]
    p1 = df["simulated_p1"][0]
    p2 = 1-p1-p0

    bg_omega0_values = df["alternative_omega0"]
    simulated_bg_omega0_values = df["simulated_omega0"]
    bg_omega1_values = df["alternative_omega1"]
    simulated_bg_omega1_values = df["simulated_omega1"]
    bg_omega2_values = df["alternative_omega2"]
    simulated_bg_omega2_values = df["simulated_omega2"]

    k_values = df["alternative_k"]
    simulated_k_values = df["simulated_k"]

    fg_omega0_values = bg_omega0_values ** k_values
    simulated_fg_omega0_values = simulated_bg_omega0_values ** simulated_k_values
    fg_omega1_values = bg_omega1_values ** k_values
    simulated_fg_omega1_values = simulated_bg_omega1_values ** simulated_k_values
    fg_omega2_values = bg_omega2_values ** k_values
    simulated_fg_omega2_values = simulated_bg_omega2_values ** simulated_k_values

    error = p0*((rel_error(simulated_bg_omega0_values, bg_omega0_values)+rel_error(simulated_fg_omega0_values, fg_omega0_values))/2) + p1*((rel_error(simulated_bg_omega1_values, bg_omega1_values)+rel_error(simulated_fg_omega1_values, fg_omega1_values))/2) + p2*((rel_error(simulated_bg_omega2_values, bg_omega2_values)+rel_error(simulated_fg_omega2_values, fg_omega2_values))/2)

    if debug:
        print("errors: ", error)
        print("error turned as nan")
        print("p0=", p0, ",p1=", p1)
        print("bg_omega0_values: ", bg_omega0_values)
        print("bg_omega1_values: ", bg_omega1_values)
        print("bg_omega2_values: ", bg_omega2_values)
        print("k_values: ", k_values)
        print("simulated_k_values: ", simulated_k_values)
        print("fg_omega0_values: ", fg_omega0_values)
        print("fg_omega1_values: ", fg_omega1_values)
        print("fg_omega2_values: ", fg_omega2_values)
        exit()

    return error

def plot_accuracy_vs_k(ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, MPComboToDf, EmpiricalMPLRThresholds, title, add_legend=False, add_ylabel=True, use_global_accuracy=False, only_of_significants=False, use_empirical=False):
    # declare fixed parameters
    tbl = 4
    mu = 8 # 4 - just a test
    pi0 = 0.5
    taxa_num = 32
    positions_num = 300

    ax.grid(False)
    rel_error_sim_vs_standard = []
    rel_error_given_mp_history_vs_standard = []
    for k in k_values_options:
        standard_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        if only_of_significants and use_empirical:
            empirical_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]  
            standard_df = standard_df.loc[(standard_df["LR"] >= empirical_threshold)]
        elif only_of_significants:
            standard_df = standard_df.loc[(standard_df["significant"] == 1)]
        given_mp_history_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
        if only_of_significants and use_empirical:
            empirical_threshold = EmpiricalMPLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            given_mp_history_df = given_mp_history_df.loc[(given_mp_history_df["LR"] >= empirical_threshold)]
        elif only_of_significants:
            given_mp_history_df = given_mp_history_df.loc[(given_mp_history_df["significant"] == 1)]
        if use_global_accuracy:
            rel_error_sim_vs_standard.append(np.mean(global_accuracy(standard_df)))
            rel_error_given_mp_history_vs_standard.append(np.mean(global_accuracy(given_mp_history_df))) #, debug=True)))
        else:
            rel_error_sim_vs_standard.append((abs(
                np.log(standard_df["alternative_k"] + 0.000001) - np.log(standard_df["simulated_k"] + 0.000001))).mean())
            rel_error_given_mp_history_vs_standard.append(abs(
                np.log(given_mp_history_df["alternative_k"] + 0.000001) - np.log(
                    given_mp_history_df["simulated_k"] + 0.000001)).mean())

    ax.plot(k_values_options, rel_error_sim_vs_standard, color='black', linestyle='solid', label="TraitRELAX", lw=2.5)
    ax.plot(k_values_options, rel_error_given_mp_history_vs_standard, color='black', linestyle="dashed",
            label="RELAX", lw=2.5)

    ax.set_xticks(k_values_options)
    ax.set_xlabel(r"$k$", fontdict={'size': 30})
    if add_ylabel:
        if not use_global_accuracy:
            ax.set_ylabel("Mean error(" + r"$\^k$" + ")", fontdict={'size': 30})
        else:
            ax.set_ylabel("Mean error(global)", fontdict={'size': 30})
    if use_global_accuracy:
        ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
    else:
        ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5])

    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best', prop={'size': 30}, frameon=False)


# plot relative error: dashed - compared to best possible (based on true history), full - based on simulated
# plot info: fixed to tbl=1, pi0=0.5, taxa_num=32, codons=600
# x axis: mu
# y axis: mean relative error of k = abs(simulated_k-inferred_k)/simulated_k (over all alternative datasets)
# 2 curves: one for TraitRELAX standard execution (full),
#           another for TraitRELAX given the true history (dashed),
def plot_accuracy_vs_mu(ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, MPComboToDf, EmpiricalMPLRThresholds, title, add_legend=False, add_ylabel=False,
                        use_global_accuracy=False, only_of_significants=False, use_empirical=False, aggregate_k=False, add_mp=False, use_boxplot=True):
    # gather the data
    tbl = 4
    pi0 = 0.5
    taxa_num = 32
    positions_num = 300
    mu_local_options = [1, 2, 4, 8, 16]

    ax.grid(False)
    if not use_boxplot:
        rel_error_sim_vs_standard = []
        rel_error_given_mp_history_vs_standard = []
        for mu in mu_local_options:
            standard_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
            if aggregate_k:
                standard_df = pd.concat([TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)] for k in k_values_options if (tbl, mu, pi0, taxa_num, positions_num, k) in TraitRELAXComboToDf])
            if only_of_significants and use_empirical:
                empirical_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
                standard_df = standard_df.loc[(standard_df["LR"] >= empirical_threshold)]
            elif only_of_significants:
                standard_df = standard_df.loc[(standard_df["significant"] == 1)]
            given_mp_history_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
            if aggregate_k:
                given_mp_history_df = pd.concat([MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)] for k in k_values_options  if (tbl, mu, pi0, taxa_num, positions_num, k) in MPComboToDf])
            if only_of_significants and use_empirical:
                empirical_threshold = EmpiricalMPLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
                given_mp_history_df = given_mp_history_df.loc[(given_mp_history_df["LR"] >= empirical_threshold)]
            elif only_of_significants:
                given_mp_history_df = given_mp_history_df.loc[(given_mp_history_df["significant"] == 1)]

            if use_global_accuracy:
                rel_error_sim_vs_standard.append(np.mean(global_accuracy(standard_df)))
                rel_error_given_mp_history_vs_standard.append(np.mean(global_accuracy(given_mp_history_df)))
            else:
                rel_error_sim_vs_standard.append((abs(
                    np.log(standard_df["alternative_k"] + 0.000001) - np.log(standard_df["simulated_k"] + 0.000001))).mean())
                rel_error_given_mp_history_vs_standard.append(abs(
                    np.log(given_mp_history_df["alternative_k"] + 0.000001) - np.log(
                        given_mp_history_df["simulated_k"] + 0.000001)).mean())

        ax.plot(mu_local_options, rel_error_sim_vs_standard, color='black', label="TraitRELAX", linestyle='solid', lw=2.5)
        if add_mp:
            ax.plot(mu_local_options, rel_error_given_mp_history_vs_standard, color='black', linestyle="dashed",
                label="RELAX", lw=2.5)

    else:
        mu_and_method_to_combo = {(1, 'TraitRELAX'): 0, (1, 'RELAX'): 1.6, (2, 'TraitRELAX'): 4.8, (2, 'RELAX'): 6.4, (4, 'TraitRELAX'): 8.0,
                                  (4, 'RELAX'): 9.6, (8, 'TraitRELAX'): 11.2, (8, 'RELAX'): 12.8, (16, 'TraitRELAX'): 14.4, (16, 'RELAX'): 16.0}
        # plot the data
        boxplots = []
        zorder_index = 1
        for mu in mu_local_options:
            traitrelax_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
            if only_of_significants and use_empirical:
                empirical_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
                traitrelax_df = traitrelax_df.loc[(traitrelax_df["LR"] >= empirical_threshold)]
            elif only_of_significants:
                traitrelax_df = traitrelax_df.loc[(traitrelax_df["significant"] == 1)]
            boxplots.append(
                ax.boxplot(traitrelax_df["alternative_k"], widths=0.5, whis=[5, 95], showfliers=False, patch_artist=True,
                           positions=[mu_and_method_to_combo[(mu, 'TraitRELAX')]],
                           zorder=zorder_index, boxprops = dict(facecolor=pretty_colors[5]), medianprops=dict(linewidth=6, color='k'))) # labels='TraitRELAX'
            zorder_index += 1
            relax_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
            if only_of_significants and use_empirical:
                empirical_threshold = EmpiricalMPLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
                relax_df = relax_df.loc[(relax_df["LR"] >= empirical_threshold)]
            elif only_of_significants:
                relax_df = relax_df.loc[(relax_df["significant"] == 1)]
            boxplots.append(
                ax.boxplot(relax_df["alternative_k"], widths=0.5, whis=[5, 95], showfliers=False, patch_artist=True,
                           positions=[mu_and_method_to_combo[(mu, 'RELAX')]],
                           zorder=zorder_index, boxprops = dict(facecolor=pretty_colors[7]), medianprops=dict(linewidth=6, color='k'))) # labels='RELAX'
            zorder_index += 1
        if add_ylabel:
            ax.set_ylabel(r"$^k$")

    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    if not use_boxplot:
        ax.set_xticks(mu_local_options)
    else:
        ax.set_xticks([0.8, 5.6, 8.8, 12, 15.2])
        ax.set_xticklabels([str(mu) for mu in mu_local_options])
    ax.set_xlabel(r"$\mu$", fontdict={'size': 30})
    if use_global_accuracy:
        ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
    else:
        if use_boxplot:
            ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        else:
            ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')

    if add_ylabel:
        if not use_boxplot:
            if not use_global_accuracy:
                ax.set_ylabel("Mean error(" + r"$\^k$" + ")", fontdict={'size': 30})
            else:
                ax.set_ylabel("Mean error(global)", fontdict={'size': 30})
        else:
            ax.set_ylabel(r"$\^k$", fontdict={'size': 30})

    if add_legend:
        if use_boxplot:
            circ1 = mpatches.Patch(facecolor=pretty_colors[5], label='TraitRELAX')
            circ2 = mpatches.Patch(facecolor=pretty_colors[7], label='RELAX')
            ax.legend(handles=[circ1, circ2], loc='best', prop={'size': 30}, frameon=False)
        else:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc='best', prop={'size': 30}, frameon=False)


# plot relative error: dashed - compared to best possible (based on true history), full - based on simulated
# plot info: fixed to tbl=1, mu=0.5, taxa_num=32, codons=600
# x axis: mu
# y axis:relative error of k = abs(simulated_k-inferred_k)/simulated_k (over all alternative datasets)
# 2 curves: one for TraitRELAX standard execution (full),
#           another for TraitRELAX given the true history (dashed)
def plot_accuracy_vs_pi0(ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, MPComboToDf, EmpiricalMPLRThresholds, title, add_legend=False,
                         add_ylabel=False, include_mp=False, use_global_accuracy=False, only_of_significants=False, use_empirical=False, aggregate_k=False):
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
        if aggregate_k:
            standard_df = pd.concat(
                [TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)] for k in k_values_options if (tbl, mu, pi0, taxa_num, positions_num, k) in TraitRELAXComboToDf])

        if only_of_significants and use_empirical:
            empirical_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            standard_df = standard_df.loc[(standard_df["LR"] >= empirical_threshold)]
        elif only_of_significants:
            standard_df = standard_df.loc[(standard_df["significant"] == 1)]
        given_mp_history_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        if aggregate_k:
            given_mp_history_df = pd.concat(
                [MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)] for k in k_values_options if (tbl, mu, pi0, taxa_num, positions_num, k) in MPComboToDf])
        if only_of_significants and use_empirical:
            empirical_threshold = EmpiricalMPLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            given_mp_history_df = given_mp_history_df.loc[(given_mp_history_df["LR"] >= empirical_threshold)]
        elif only_of_significants:
            given_mp_history_df = given_mp_history_df.loc[(given_mp_history_df["significant"] == 1)]
        if use_global_accuracy:
            rel_error_sim_vs_standard.append(np.mean(global_accuracy(standard_df)))
            rel_error_given_mp_history_vs_standard.append(np.mean(global_accuracy(given_mp_history_df)))
        else:
            rel_error_sim_vs_standard.append((abs(
                np.log(standard_df["alternative_k"] + 0.000001) - np.log(standard_df["simulated_k"] + 0.000001))).mean())
            rel_error_given_mp_history_vs_standard.append((abs(
                np.log(given_mp_history_df["alternative_k"] + 0.000001) - np.log(
                    given_mp_history_df["simulated_k"] + 0.000001))).mean())

    ax.plot(pi0_options, rel_error_sim_vs_standard, color='black', linestyle='solid', label="TraitRELAX", lw=2.5)
    ax.plot(pi0_options, rel_error_given_mp_history_vs_standard, color='black', linestyle="dashed",
            label="RELAX", lw=2.5)
    if include_mp:
        ax.plot(pi0_options, rel_error_given_true_history_vs_standard, color='black', linestyle="--",
                label="TraitRELAX with true history")

    ax.set_xticks(pi0_options)
    ax.set_xlabel(r"$\pi_0$", fontdict={'size': 30})
    if add_ylabel:
        if not use_global_accuracy:
            ax.set_ylabel("Mean error(" + r"$\^k$" + ")", fontdict={'size': 30})
        else:
            ax.set_ylabel("Mean error(global)", fontdict={'size': 30})

    if use_global_accuracy:
        ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
    else:
        ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_xticks(pi0_options)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_ylabel("Mean error (" + r"$\^k$" + ")", fontdict={'size': 30, 'color': 'white'})

    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best', prop={'size': 30}, frameon=False)

# plot relative error: dashed - compared to best possible (based on true history), full - based on simulated
# plot info: fixed to tbl=1, pi0=0.5, taxa_num=32, codons=600
# x axis: mu
# y axis: mean relative error of k = abs(simulated_k-inferred_k)/simulated_k (over all alternative datasets)
# 2 curves: one for TraitRELAX standard execution (full),
#           another for TraitRELAX given the true history (dashed),
def plot_accuracy_vs_tbl(ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, MPComboToDf, EmpiricalMPLRThresholds, title, add_legend=False,
                         add_ylabel=False, include_mp=False, use_global_accuracy=False, only_of_significants=False, use_empirical=False, aggregate_k=False):
    # gather the data
    tbl_to_mu = {1: 32, 4: 8, 8: 4, 16: 2, 32: 1}
    pi0 = 0.5
    taxa_num = 32
    positions_num = 300

    ax.grid(False)
    ax.grid(False)
    rel_error_sim_vs_standard = []
    rel_error_given_mp_history_vs_standard = []
    tbls = list(tbl_to_mu.keys())
    tbls.sort()
    for tbl in tbls:
        mu = tbl_to_mu[tbl]
        standard_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        if aggregate_k:
            standard_df = pd.concat([TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)] for k in k_values_options if (tbl, mu, pi0, taxa_num, positions_num, k) in TraitRELAXComboToDf])
        if only_of_significants and use_empirical:
            empirical_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            standard_df = standard_df.loc[(standard_df["LR"] >= empirical_threshold)]
        elif only_of_significants:
            standard_df = standard_df.loc[(standard_df["significant"] == 1)]
        given_mp_history_df = MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 0.5)]
        if aggregate_k:
            given_mp_history_df = pd.concat(
                [MPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)] for k in k_values_options if (tbl, mu, pi0, taxa_num, positions_num, k) in MPComboToDf])
        if only_of_significants and use_empirical:
            empirical_threshold = EmpiricalMPLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            given_mp_history_df = given_mp_history_df.loc[(given_mp_history_df["LR"] >= empirical_threshold)]
        elif only_of_significants:
            given_mp_history_df = given_mp_history_df.loc[(given_mp_history_df["significant"] == 1)]

        if use_global_accuracy:
            rel_error_sim_vs_standard.append(np.mean(global_accuracy(standard_df)))
            rel_error_given_mp_history_vs_standard.append(np.mean(global_accuracy(given_mp_history_df)))
        else:
            rel_error_sim_vs_standard.append((abs(
                np.log(standard_df["alternative_k"] + 0.000001) - np.log(
                    standard_df["simulated_k"] + 0.000001))).mean())
            rel_error_given_mp_history_vs_standard.append(abs(
                np.log(given_mp_history_df["alternative_k"] + 0.000001) - np.log(
                    given_mp_history_df["simulated_k"] + 0.000001)).mean())

    ax.plot(tbls, rel_error_sim_vs_standard, color='black', linestyle='solid', label="TraitRELAX", lw=2.5)
    if include_mp:
        ax.plot(tbls, rel_error_given_mp_history_vs_standard, color='black', linestyle="dashed",
                label="RELAX", lw=2.5)

    ax.set_xticks(tbls)
    ax.set_xlabel("Tree length", fontdict={'size': 30})
    if add_ylabel:
        if not use_global_accuracy:
            ax.set_ylabel("Mean error(" + r"$\^k$" + ")", fontdict={'size': 30})
        else:
            ax.set_ylabel("Mean error(global)", fontdict={'size': 30})

    if use_global_accuracy:
        ax.set_yticks([0, 1, 2, 3, 4, 5, 6])
    else:
        ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
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
# In each plot, 6 curves in 3 colors (color per positions number) - full for TraitRELAX, dashed for TraitTraitRELAX with true history
# in green - scatter the simulated value of k (with legend)
def plot_accuracy_vs_k_across_positions_num(taxa_num, ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, title, add_legend=False,
                                            plot_labels=True, only_of_significants=False, use_empirical=False, linestyle='solid'):
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
            if only_of_significants and use_empirical:
                empirical_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
                full_df = full_df.loc[(full_df["LR"] >= empirical_threshold)]
            elif only_of_significants:
                full_df = full_df.loc[(full_df["significant"] == 1)]
            alternative_k_values = list(full_df["alternative_k"])
            for i in range(len(alternative_k_values)):
                if alternative_k_values[i] == 0:
                    alternative_k_values[i] = 0.000001
            k_error_values = (abs(np.log(np.array(alternative_k_values)) - np.log(full_df["simulated_k"])))
            y_full_values.append(np.mean(k_error_values))
        # plot the data in ax
        ax.plot(k_values_options, y_full_values, color=colors[codon_positions_num_options.index(positions_num)],
                label=str(positions_num) + " sites", lw=2.5, linestyle=linestyle)

    if plot_labels:
        ax.set_xlabel(r"$k$", fontdict={'size': 30})
        ax.set_ylabel("Mean error (" + r"$\^k$" + ")", fontdict={'size': 30})

    ax.set_xticks(k_values_options)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5])

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best', prop={'size': 30}, frameon=False)


# fixed parameters: tbl=4, mu=8, pi0=0.5
# plot assessment of inference with respect to k
# x axis: simulated value of k
# y axis: mean value inferred value of k
# add errorbars if not too noisy
# a panel with plot for each taxa num
# In each plot, 6 curves in 3 colors (color per positions number) - full for TraitRELAX, dashed for TraitTraitRELAX with true history
# in green - scatter the simulated value of k (with legend)
def plot_accuracy_vs_k_across_taxa_num(positions_num, ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, title, add_legend=False,
                                       plot_labels=True, only_of_significants=False, use_empirical=False, linestyle='solid'):
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
            if only_of_significants and use_empirical:
                empirical_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
                full_df = full_df.loc[(full_df["LR"] >= empirical_threshold)]
            elif only_of_significants:
                full_df = full_df.loc[(full_df["significant"] == 1)]
            alternative_k_values = list(full_df["alternative_k"])
            for i in range(len(alternative_k_values)):
                if alternative_k_values[i] == 0:
                    alternative_k_values[i] = 0.000001
            k_error_values = (abs(np.log(np.array(alternative_k_values)) - np.log(full_df["simulated_k"])))
            y_full_values.append(np.mean(k_error_values))
        # plot the data in ax
        ax.plot(k_values_options, y_full_values, color=colors[taxa_num_options.index(taxa_num)],
                label=str(taxa_num) + " taxa", lw=2.5, linestyle=linestyle)

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
    # ax.grid(False)
    ax.grid(color='white')
    ax.plot_trisurf(param1_values, param2_values, logl_values, linewidth=1, antialiased=True,
                    cmap='viridis')  # rstride=1, cstride=1, cmap='viridis', edgecolor='none', alpha=.8)

    # scatter the ML point
    mlX, mlY = np.meshgrid([0.890321], [2.12883])
    mlZ = np.array([-15170.1873250066])
    ax.plot(mlX, mlY, mlZ, color='red', marker='o', label="MLE", markersize=14)

    # scatter the simulated point
    mlX, mlY = np.meshgrid([0.8], [2])
    mlZ = np.array([-15190.836444])
    ax.plot(mlX, mlY, mlZ, color='k', marker='o', label="True", markersize=14)

    ax.set_xticks([0, 0.5, 1, 1.5, 2])
    ax.set_xlabel("\n\n" + r"$k$", fontdict={'size': 30}, labelpad=8)
    ax.set_yticks([1, 2, 3, 4])
    ax.set_ylabel("\n\n" + r"$\omega_2$", fontdict={'size': 30}, labelpad=8)
    # ax.set_zlabel("\n\nlog likelihood", fontdict={'size': 30}, labelpad=30)

    legend_elements = [Line2D([0], [0], marker='o', color='w', label='MLE',
                          markerfacecolor='red', markersize=14), Line2D([0], [0], marker='o', color='w', label='True',
                          markerfacecolor='k', markersize=14)]
    ax.legend(handles=legend_elements, prop={'size': 30}, frameon=False, bbox_to_anchor=(0.95, 1))
    ax.dist = dist
    ax.tick_params(axis='z', which='major', pad=24)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')


def plot_inferred_vs_simulated_k_across_positions_num(taxa_num, ax, combo_to_df, empirical_thresholds, title, plot_labels=False, only_of_significants=False, use_empirical=False):
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
            if only_of_significants and use_empirical:
                empirical_threshold = empirical_thresholds[(tbl, mu, pi0, taxa_num, positions_num)]
                df = df.loc[(df["LR"] >= empirical_threshold)]
            elif only_of_significants:
                df = df.loc[(df["significant"] == 1)]
            boxplots.append(
                ax.boxplot(df["alternative_k"], widths=0.2, whis=[5, 95], showfliers=False, patch_artist=True,
                           positions=[positions_to_combo[(k, positions_num)]],
                           zorder=zorder_index))  # , color=colors[codon_positions_num_options.index[positions_num]]
            zorder_index += 1
        for boxplot in boxplots:
            color = colors[boxplots.index(boxplot)]
            patch = boxplot['boxes'][0]
            patch.set_facecolor(color)
            # patch.set_linewidth(2)

        x_for_scatter = [0, 0.3, 0.6, 1.2, 1.5, 1.8, 2.4, 2.7, 3, 3.6, 3.9, 4.2, 4.8, 5.1, 5.4, 6, 6.3, 6.6, 7.2, 7.5,
                         7.8]
        y_for_scatter = [0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.8, 0.8, 0.8, 1, 1, 1, 1.2, 1.2, 1.2, 1.6, 1.6, 1.6, 2, 2, 2]
        ax.scatter(x_for_scatter, y_for_scatter, label=" simulated values", color='green', zorder=zorder_index)

    if plot_labels:
        ax.set_xlabel(r"$k$")
        ax.set_ylabel(r"$^k$")

    ax.set_xticklabels(
        ["", "0.2", "", "", "0.5", "", "", "0.8", "", "", "1", "", "", "1.2", "", "", "1.6", "", "", "2", ""])
    ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5])
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')


# fix: tbl=1, pi0=0.5, mu=0.5, taxa_num=32, codons=600
# a panel of 5 plots: one per k.
# each plot has 4 boxplots: k, omega0, omega1, omega2
# on each boxplot, scatter simulated value in green
def plot_k_to_omegas_accuracy(ax, comboToDf, empricial_thresholds, title, only_of_significants=False, use_empirical=False):
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
        if only_of_significants and use_empirical:
            empirical_threshold = empricial_thresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            source_df = source_df.loc[(source_df["LR"] >= empirical_threshold)]
        elif only_of_significants:
            source_df = source_df.loc[(source_df["significant"] == 1)]

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
    matplotlib.rc('legend', fontsize=16)
    matplotlib.rc('legend', frameon=False)
    ax.grid(False)
    sns.boxplot(x="simulated value of k", y='value', hue='variable', data=melted_df, ax=ax, palette="Set3", showfliers=False)
    ax.legend().set_title('')
    # ax.legend().set_loc('upper right')
    # ax.legend().set_fontsize()

    # scatter the simulated values
    x_for_scatter = [-0.3, -0.1, 0.1, 0.3, 0.7, 0.9, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.7, 2.9, 3.1, 3.3, 3.7, 3.9, 4.1,
                     4.3]
    ax.scatter(x=np.array(x_for_scatter), y=np.array(true_values), label="Simulated value", edgecolors=None,
               color="k", marker="s")

    ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4])
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.set_ylabel(r"Parameter value", fontdict={'size': 30})
    ax.set_xlabel(r"$k$", fontdict={'size': 30})

    matplotlib.rc('legend', frameon=False)
    matplotlib.rc('legend', fontsize=30)

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
            boxplots.append(ax.boxplot(k_values, widths=0.2, whis=[5, 95], showfliers=False, patch_artist=True,
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

    ax.set_ylabel("", fontdict={'size': 30})

    # add legend
    # create custom legend
    custom_lines = []
    custom_names = []
    for positions_num in codon_positions_num_options:
        custom_lines.append(Line2D([0], [0], color=colors[codon_positions_num_options.index(positions_num)], lw=6))
        custom_names.append(str(positions_num) + " sites")
    ax.legend(custom_lines, custom_names, loc='best', prop={'size': 30}, frameon=False)

    ax.grid(False)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.set_xticklabels(["", "0.2", "", "", "1", "", "", "2", ""])
    ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5])


def plot_k_distribution_across_k_and_taxanum(positions_num, ax, comboToDf, title):
    # declare fixed parameters
    tbl = 4
    mu = 8
    pi0 = 0.5
    zorder_index = 1

    # collect the  data
    k_values_to_plot = [0.2, 1, 2]
    taxanum_to_combo = {(0.2, 16): 0, (0.2, 32): 0.3, (0.2, 64): 0.6, (1, 16): 1.2, (1, 32): 1.5, (1, 64): 1.8,
                          (2, 16): 2.4, (2, 32): 2.7, (2, 64): 3}
    for k in k_values_to_plot:
        boxplots = []
        for taxa_num in taxa_num_options:
            df = comboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
            filetred_df = df.loc[(df["alternative_k"] < 3)]
            k_values = filetred_df["alternative_k"]
            boxplots.append(ax.boxplot(k_values, widths=0.2, whis=[5, 95], showfliers=False, patch_artist=True,
                                       positions=[taxanum_to_combo[(k, taxa_num)]],
                                       zorder=zorder_index, medianprops=dict(linewidth=6, color='k')))  # , color=colors[codon_positions_num_options.index[positions_num]]
            zorder_index += 1
        for boxplot in boxplots:
            box = boxplot['boxes'][0]
            box.set(edgecolor = 'white')
            box.set(linewidth = 0.75)
            box.set(facecolor  = pretty_colors[boxplots.index(boxplot)])
            median = boxplot['medians'][0]
            median.set(color='k', linewidth=2,)

    x_for_scatter = [0, 0.3, 0.6, 1.2, 1.5, 1.8, 2.4, 2.7, 3]
    ax.scatter(x=x_for_scatter, y=[0.2, 0.2, 0.2, 1, 1, 1, 2, 2, 2], label="Simulated value", edgecolors=None,
               color="k", marker="s", zorder=zorder_index, linewidths=6)

    # add legend
    # create custom legend
    # custom_lines = []
    # custom_names = []
    # for taxa_num in taxa_num_options:
    #     custom_lines.append(Line2D([0], [0], color=colors[taxa_num_options.index(taxa_num)], lw=6))
    #     custom_names.append(str(taxa_num) + " taxa")
    # ax.legend(custom_lines, custom_names, loc='best', prop={'size': 30}, frameon=False)
    a_val = 0.6
    circ1 = mpatches.Patch(facecolor=pretty_colors[0], label='16 taxa') # ,alpha=a_val)
    circ2 = mpatches.Patch(facecolor=pretty_colors[1], label='32 taxa') # ,alpha=a_val)
    circ3 = mpatches.Patch(facecolor=pretty_colors[2], label='64 taxa') # ,alpha=a_val)
    ax.legend(handles=[circ1, circ2, circ3], loc='best', prop={'size': 30}, frameon=False)

    ax.grid(False)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.set_xticklabels(["", "0.2", "", "", "1", "", "", "2", ""])
    ax.set_yticks([0, 0.5, 1, 1.5, 2, 2.5])
    ax.set_xlabel(r"$k$", fontdict={'size': 30})
    ax.set_ylabel(r"$\^k$", fontdict={'size': 30})


def plot_accuracy_analysis(TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, grid_data_path, res_output_path,
                           supp_output_path_1, supp_output_path_2, supp_output_path_3):


    # alternative figure 0 (for results)
    plt.grid(False)
    fig = plt.figure(figsize=[3 * 8.5 + 2, 1 * 7.58 + 2], frameon=True)
    ax1 = fig.add_subplot(1, 3, 1)
    plot_k_distribution_across_k_and_taxanum(600, ax1, TraitRELAXComboToDf, "a\n")
    ax2 = fig.add_subplot(1, 3, 2)
    plot_accuracy_vs_mu(ax2, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "b\n", add_ylabel=True, add_legend=True, only_of_significants=True, use_empirical=True, use_boxplot=True)
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    plot_2d_grid(ax3, grid_data_path, "c\n", 9)
    # Rotate angle
    l2 = np.array((1, 1))
    trans_angle = plt.gca().transData.transform_angles(np.array((86,)),
                                                       l2.reshape((1, 2)))[0]
    fig.text(1.04, 0.4, 'Log likelihood', va='center', rotation_mode='anchor', fontdict={'size': 30}, rotation=86) #=trans_angle)
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(res_output_path, bbox_inches='tight') #, transparent=True)
    plt.clf()

    # figure 1 (for supp materials)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=2, ncols=2, sharey='none', sharex='none', figsize=[2 * 8.2 + 2, 2 * 7.58 + 2], frameon=True)
    plot_accuracy_vs_mu(axis[0][0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, MPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "a\n", add_ylabel=False, use_global_accuracy=False) #, only_of_significants=True, use_empirical=True)
    plot_accuracy_vs_pi0(axis[0][1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, MPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "b\n", add_ylabel=False, use_global_accuracy=False) #, only_of_significants=True, use_empirical=True)
    plot_accuracy_vs_tbl(axis[1][0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, MPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "c\n", add_legend=False,
                         add_ylabel=False, use_global_accuracy=False) #, only_of_significants=True, use_empirical=True)
    plot_accuracy_vs_k(axis[1][1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, MPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "d\n", add_legend=True, add_ylabel=False, use_global_accuracy=False) #, only_of_significants=True, use_empirical=True)
    fig.text(-0.02, 0.5, "Mean error(k)", va='center', rotation='vertical', fontdict={'size': 30})
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(supp_output_path_1, bbox_inches='tight') #, transparent=True)
    plt.clf()

    # figure 2 (for supp materials)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=3, ncols=2, sharey='none', sharex='none', figsize=[2 * 12 + 2, 3 * 7.58 + 2], frameon=True)
    plot_inferred_vs_simulated_k_across_positions_num(16, axis[0][0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "a\n", plot_labels=False) #, only_of_significants=True, use_empirical=True)
    plot_accuracy_vs_k_across_positions_num(16, axis[0][1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "b\n", add_legend=False,
                                            plot_labels=False) #, only_of_significants=True, use_empirical=True)
    plot_inferred_vs_simulated_k_across_positions_num(32, axis[1][0], TraitRELAXComboToDf,
                                                      EmpiricalTraitRELAXLRThresholds, "c\n",
                                                      plot_labels=True)  # , only_of_significants=True, use_empirical=True)
    plot_accuracy_vs_k_across_positions_num(32, axis[1][1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "d\n", add_legend=False,
                                            plot_labels=True) #, only_of_significants=True, use_empirical=True)
    plot_inferred_vs_simulated_k_across_positions_num(64, axis[2][0], TraitRELAXComboToDf,
                                                      EmpiricalTraitRELAXLRThresholds, "e\n",
                                                      plot_labels=False)  # , only_of_significants=True, use_empirical=True)
    plot_accuracy_vs_k_across_positions_num(64, axis[2][1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "f\n", add_legend=True,
                                            plot_labels=False) #, only_of_significants=True, use_empirical=True)

    # fig.text(0.53, -0.02, r"$k$", ha='center', fontdict={'size': 30})
    fig.text(-0.02, 0.5, "Mean(" + r"$\^k$" + ")", va='center', rotation='vertical', fontdict={'size': 30})
    fig.text(0.5, -0.02, r"$\^k$" , va='center', rotation='horizontal', fontdict={'size': 30})
    fig.subplots_adjust()
    fig.tight_layout(pad=1)
    plt.savefig(supp_output_path_2, bbox_inches='tight') #, transparent=True)
    plt.clf()

    # figure 3 (for supp materials)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=1, ncols=1, sharex="none", sharey="none", figsize=[1 * 8.2 + 2, 7.58], frameon=True)
    plot_k_to_omegas_accuracy(axis, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "", only_of_significants=True, use_empirical=True)
    fig.subplots_adjust()
    fig.tight_layout(pad=1)
    plt.savefig(supp_output_path_3, bbox_inches='tight') #, transparent=True)
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

    # plot a line for each taxa_num
    fig, axis = plt.subplots(nrows=1, ncols=1, frameon=True)
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
    plt.savefig(output_path, bbox_extra_artists=(lgd,), bbox_inches='tight', figsize=[11, 7.58]) #, transparent=True)
    plt.clf()


#####################################################################################

# plot info: fixed to tbl=4, pi0=0.5, mu=8, taxa_num=32, codons=300
# x axis: simulated value of k
# y axis: %significant (over all alternative datasets)
def plot_direction_vs_k_across_posnum(taxa_num, ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, title,
                                  empirical=True, add_legend=False, plot_labels=True, direction='wrong', linestyle='solid'):
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
            if k == 1:
                continue
            full_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
            if empirical:
                full_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
                full_direction = 0
                if direction == 'wrong':
                    full_direction = full_df.loc[(full_df["LR"] >= full_threshold) & (full_df["correctDirection"] == 0)].shape[0] / full_df.shape[0]
                else:
                    full_direction = full_df.loc[(full_df["LR"] >= full_threshold) & (full_df["correctDirection"] == 1)].shape[0] / full_df.shape[0]
                y_full_values.append(full_direction)
            else:
                if direction == 'wrong':
                    full_direction = full_df.loc[(full_df["pvalue"] <= 0.05) & (full_df["correctDirection"] == 0)].shape[0] / full_df.shape[0]
                else:
                    full_direction = full_df.loc[(full_df["pvalue"] <= 0.05) & (full_df["correctDirection"] == 1)].shape[0] / full_df.shape[0]
                y_full_values.append(full_direction)
        # plot the data in ax
        ax.plot([k for k in k_values_options if k!=1], y_full_values, color=colors[codon_positions_num_options.index(positions_num)],
                marker=markers[taxa_num_options.index(taxa_num)], lw=2.5, label=str(positions_num) + " sites",
                linestyle=linestyle, markersize=16)
        if plot_labels:
            ax.set_xlabel(r"$k$", fontdict={'size': 30})
            ax.set_ylabel("Null rejected in " + direction + "\nselection intensity direction", fontdict={'size': 30})

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='upper center', prop={'size': 30}, frameon=False)

    # vals = [0, 0.05, 0.2, 0.4, 0.6, 0.8, 1]
    # ax.set_ylim(0,1)
    # ax.set_yticks(vals, ['{:,.0%}'.format(x) for x in vals])
    ax.yaxis.set_major_formatter(FuncFormatter(convertToPercent))
    ax.set_xticks([k for k in k_values_options if k!=1])

# plot info: fixed to tbl=4, pi0=0.5, mu=8, taxa_num=32, codons=300
# x axis: simulated value of k
# y axis: %significant (over all alternative datasets)
def plot_direction_vs_k_across_taxanum(positions_num, ax, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, title,
                                   empirical=True, add_legend=False, plot_labels=True, direction='wrong', linestyle='solid'):
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
            if k == 1:
                continue
            full_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
            if empirical:
                full_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
                if direction == 'wrong':
                    full_direction = full_df.loc[(full_df["LR"] >= full_threshold) & (full_df["correctDirection"] == 0)].shape[0] / full_df.shape[0]
                else:
                    full_direction = full_df.loc[(full_df["LR"] >= full_threshold) & (full_df["correctDirection"] == 1)].shape[0] / full_df.shape[0]
                y_full_values.append(full_direction)
            else:
                if direction == 'wrong':
                    full_direction = full_df.loc[(full_df["pvalue"] <= 0.05) & (full_df["correctDirection"] == 0)].shape[0] / full_df.shape[0]
                else:
                    full_direction = full_df.loc[(full_df["pvalue"] <= 0.05) & (full_df["correctDirection"] == 1)].shape[0] / full_df.shape[0]
                y_full_values.append(full_direction)
        # plot the data in ax
        ax.plot([k for k in k_values_options if k!=1], y_full_values, color=colors[codon_positions_num_options.index(positions_num)],
                marker=markers[taxa_num_options.index(taxa_num)], lw=2.5, label=str(taxa_num) + " taxa", linestyle=linestyle,
                markersize=16)
        if plot_labels:
            ax.set_xlabel(r"$k$", fontdict={'size': 30})

    # vals = [0, 0.2, 0.4, 0.6, 0.8, 1]
    # ax.set_yticks(vals, ['{:,.0%}'.format(x) for x in vals])
    # ax.yaxis.set_major_formatter(FuncFormatter(convertToPercent))
    ax.set_xticks([k for k in k_values_options if k!=1])

    if add_legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='upper center', prop={'size': 30}, frameon=False)

def report_power_fpr_comparison(TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds):

    tbl = 4
    mu = 8
    pi0 = 0.5

    # first, TraitRELAX
    print("******************************** TraitRELAX ********************************")
    print("taxa_num\tpositions_num\tk\ttheoretical_power\ttheoretical_FPR\ttheoretical_correct_direction\ttheoretical_wrong_direction\tempirical_power\tempirical_FPR\tempirical_correct_direction\tempirical_wrong_direction")
    for taxa_num in taxa_num_options:
        for positions_num in codon_positions_num_options:
            fpr_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 1)]
            theoretical_fpr = (fpr_df.loc[(fpr_df["significant"] == 1)].shape[0]) / fpr_df.shape[0]
            empirical_threshold = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            empirical_fpr = fpr_df[fpr_df.LR >= empirical_threshold].shape[0] / fpr_df.shape[0]
            for k in k_values_options:
                if k != 1:
                    df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
                    theoretical_power = (df.loc[(df["significant"] == 1)].shape[0]) / df.shape[0]
                    theoretical_correct_dir, theoretical_wrong_dir = get_div_for_null_rejected(df, empirical_threshold, theoretical=1)
                    empirical_power = (df.loc[(df["LR"] >= empirical_threshold)].shape[0]) / df.shape[0]
                    empirical_correct_dir, empirical_wrong_dir = get_div_for_null_rejected(df, empirical_threshold, theoretical=0)
                    print(taxa_num, "\t", positions_num, "\t", k, "\t", theoretical_power, "\t", theoretical_fpr, "\t", theoretical_correct_dir, "\t", theoretical_wrong_dir, "\t", empirical_power, "\t", empirical_fpr, "\t", empirical_correct_dir, "\t", empirical_wrong_dir)


    # Then, HyPhy
    print("******************************** RELAX_HyPhy ********************************")
    print("taxa_num\tpositions_num\tk\ttheoretical_power\ttheoretical_FPR\ttheoretical_correct_direction\ttheoretical_wrong_direction\tempirical_power\tempirical_FPR\tempirical_correct_direction\tempirical_wrong_direction")
    for taxa_num in taxa_num_options:
        for positions_num in codon_positions_num_options:
            fpr_df = HyPhyMPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, 1)]
            theoretical_fpr = (fpr_df.loc[(fpr_df["significant"] == 1)].shape[0]) / fpr_df.shape[0]
            empirical_threshold = EmpiricalHyPhyMPRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            empirical_fpr = fpr_df[fpr_df.LR >= empirical_threshold].shape[0] / fpr_df.shape[0]
            for k in k_values_options:
                if k != 1:
                    df = HyPhyMPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k)]
                    theoretical_power = (df.loc[(df["significant"] == 1)].shape[0]) / df.shape[0]
                    theoretical_correct_dir, theoretical_wrong_dir = get_div_for_null_rejected(df, empirical_threshold, theoretical=1)
                    empirical_power = (df.loc[(df["LR"] >= empirical_threshold)].shape[0]) / df.shape[0]
                    empirical_correct_dir, empirical_wrong_dir = get_div_for_null_rejected(df, empirical_threshold, theoretical=0)
                    print(taxa_num, "\t", positions_num, "\t", k, "\t", theoretical_power, "\t", theoretical_fpr, "\t", theoretical_correct_dir, "\t", theoretical_wrong_dir, "\t", empirical_power, "\t", empirical_fpr, "\t", empirical_correct_dir, "\t", empirical_wrong_dir)

def plot_correct_direction_distribution(ax, title, TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, by_taxa=1, k_options=(0.2, 0.5, 0.8, 1.2, 1.6, 2), add_legend=False, out_of_all=0):

    tbl = 4
    mu = 8
    pi0 = 0.5

    traitrelax_data = []
    hyphy_data = []

    # create a data structure of your result
    if by_taxa:
        positions_num = 300
        for taxa_num in taxa_num_options:
            # collect TraitRELAX data
            TraitRELAX_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k_options[0])]
            for i in range(1,len(k_options)):
                TraitRELAX_df.append(TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k_options[i])])
            traitrelax_empirical_cutoff = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            traitrelax_correct_dir, traitrelax_wrong_dir = get_div_for_null_rejected(TraitRELAX_df, traitrelax_empirical_cutoff, theoretical=0, out_of_all=out_of_all)
            traitrelax_data.append(traitrelax_correct_dir)

            # collect hyphy data
            hyphy_df = HyPhyMPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k_options[0])]
            for i in range(1, len(k_options)):
                hyphy_df.append(HyPhyMPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k_options[i])])
            hyphy_empirical_cutoff = EmpiricalHyPhyMPRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            hyphy_correct_dir, hyphy_wrong_dir = get_div_for_null_rejected(hyphy_df, hyphy_empirical_cutoff,
                                                                                     theoretical=1, out_of_all=out_of_all)
            hyphy_data.append(hyphy_correct_dir)
    else:
        taxa_num = 32
        for positions_num in codon_positions_num_options:
            # collect TraitRELAX data
            TraitRELAX_df = TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k_options[0])]
            for i in range(1,len(k_options)):
                TraitRELAX_df.append(TraitRELAXComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k_options[i])])
            traitrelax_empirical_cutoff = EmpiricalTraitRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            traitrelax_correct_dir, traitrelax_wrong_dir = get_div_for_null_rejected(TraitRELAX_df, traitrelax_empirical_cutoff, theoretical=0)
            traitrelax_data.append(traitrelax_correct_dir)

            # collect hyphy data
            hyphy_df = HyPhyMPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k_options[0])]
            for i in range(1, len(k_options)):
                hyphy_df.append(HyPhyMPComboToDf[(tbl, mu, pi0, taxa_num, positions_num, k_options[i])])
            hyphy_empirical_cutoff = EmpiricalHyPhyMPRELAXLRThresholds[(tbl, mu, pi0, taxa_num, positions_num)]
            hyphy_correct_dir, hyphy_wrong_dir = get_div_for_null_rejected(hyphy_df, hyphy_empirical_cutoff,
                                                                                     theoretical=1)
            hyphy_data.append(hyphy_correct_dir)

    ax.grid(False)
    if by_taxa == 1:
        index = ['16 taxa', '32 taxa', '64 taxa']
    else:
        index = ['150 sites', '300 sites', '600 sites']

    df = pd.DataFrame({'TraitRELAX': traitrelax_data,
                       'RELAX': hyphy_data}, index=index)

    # set a greyscale colormap
    greyColorMap = ListedColormap(colors)
    df.plot.bar(ax=ax, rot=0, legend=False, colormap=greyColorMap, figsize=(8.2, 6)) #5.58))
    vals = [0, 0.5, 1]
    ax.set_yticks(vals, ['{:,.0%}'.format(x) for x in vals])
    ax.yaxis.set_major_formatter(FuncFormatter(convertToPercent))
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')

def plot_comparison_to_hyphy(TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, output_path_1, output_path_2, output_path_3, output_path_4, output_path_5):

    # figure 0 (for results)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=2, ncols=2, sharex="all", sharey="all", figsize=[2 * 8.2 + 2, 2* 7.58], frameon=True)
    plot_direction_vs_k_across_posnum(32, axis[0][0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "a\n",
                                  empirical=True, add_legend=True, plot_labels=False, direction='right')
    plot_direction_vs_k_across_taxanum(300, axis[0][1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "b\n",
                                   empirical=True, add_legend=True, plot_labels=False, direction='right')
    plot_direction_vs_k_across_posnum(32, axis[1][0], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "c\n",
                                  empirical=False, add_legend=False, plot_labels=False, direction='right', linestyle='dashed')
    plot_direction_vs_k_across_taxanum(300, axis[1][1], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "d\n",
                                   empirical=False, add_legend=False, plot_labels=False, direction='right', linestyle='dashed')
    fig.text(-0.015, 0.5, 'Null rejected in correct direction', va='center', rotation='vertical', fontdict={'size': 30})
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(output_path_1, bbox_inches='tight') #, transparent=True)
    plt.clf()

    # figure 1: power plots for TraitRELAX by empirical test (A,B,C) and for HyPhy by theoretical test (D,E,F) with addition of required FPR line
    plt.grid(False)
    fig, axis = plt.subplots(nrows=2, ncols=2, sharey='all', sharex='all', figsize=[2 * 8.2 + 2, 2 * 7.58], frameon=True)
    plot_direction_vs_k_across_posnum(32, axis[0][0], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "a\n",
                                  empirical=True, add_legend=True, plot_labels=False, direction='wrong')
    plot_direction_vs_k_across_taxanum(300, axis[0][1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, "b\n",
                                  empirical=True, add_legend=True, plot_labels=False, direction='wrong')
    plot_direction_vs_k_across_posnum(32, axis[1][0], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "c\n",
                                  empirical=False, add_legend=False, plot_labels=False, direction='wrong', linestyle='dashed')
    plot_direction_vs_k_across_taxanum(300, axis[1][1], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "d\n",
                                  empirical=False, add_legend=False, plot_labels=False, direction='wrong', linestyle='dashed')
    fig.text(-0.015, 0.5, 'Null rejected in wrong direction', va='center', rotation='vertical', fontdict={'size': 30})
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(output_path_2, bbox_inches='tight') #, transparent=True)
    plt.clf()

    # # figure 2: 2 histograms: in one the x-axis is the #taxa for 300 positions and in another it is #positions for 32 taxa
    # #           in each histogram, we plot the fraction of correct directions of TraitREAX and HyPhy
    # #           two options are enabled: aggregation across all the k values (except 1) and focus on a single k. This is done by passing a list argument of "k options"
    # plt.grid(False)
    # fig, axis = plt.subplots(nrows=2, ncols=1, sharey='none', sharex='none', figsize=[2 * 8.2 + 4, 1 * 7.58], frameon=True)
    # plot_correct_direction_distribution(axis[0], "a\n", TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, by_taxa=1, add_legend=True, out_of_all=1)
    # plot_correct_direction_distribution(axis[1], "b\n", TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, by_taxa=0, add_legend=False, out_of_all=1)
    # # plot custom legend
    # circ1 = mpatches.Patch(facecolor='lightgrey', label='TraitRELAX')
    # circ2 = mpatches.Patch(facecolor='k', label='RELAX')
    # axis[0].legend(handles=[circ1, circ2], loc='right', prop={'size': 26}, frameon=False, bbox_to_anchor=(1, 1.35))
    #
    # fig.text(-0.015, 0.5, 'Null rejected in correct direction', va='center', rotation='vertical', fontdict={'size': 30})
    # fig.subplots_adjust()
    # fig.tight_layout()
    # plt.savefig(output_path_2, bbox_inches='tight', transparent=True)
    # plt.clf()

    # figure 3: accuracy comparison between TraitRELAX and HyPhy
    plt.grid(False)
    fig, axis = plt.subplots(nrows=1, ncols=3, sharey='none', sharex='none', figsize=[3 * 8.2 + 4, 1 * 7.58], frameon=True)
    plot_k_distribution_across_k_and_taxanum(600, axis[0], HyPhyMPComboToDf, "a\n")
    plot_accuracy_vs_k(axis[1], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "b\n", add_legend=False, add_ylabel=True, use_global_accuracy=False, only_of_significants=True, use_empirical=True)
    plot_accuracy_vs_mu(axis[2], TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "c\n", add_legend=True, add_ylabel=True, use_global_accuracy=False, only_of_significants=False, use_empirical=False, aggregate_k=False, add_mp=True, use_boxplot=False)
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(output_path_3, bbox_inches='tight') #, transparent=True)
    plt.clf()


    # figure 2 (for supp materials)
    plt.grid(False)
    fig, axis = plt.subplots(nrows=3, ncols=2, sharey='none', sharex='none', figsize=[2 * 12 + 2, 3 * 7.58 + 2], frameon=True)
    plot_inferred_vs_simulated_k_across_positions_num(16, axis[0][0], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "a\n", plot_labels=False) #, only_of_significants=True, use_empirical=True)
    plot_accuracy_vs_k_across_positions_num(16, axis[0][1], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "b\n",
                                            add_legend=False, plot_labels=False)  # , only_of_significants=True, use_empirical=True)
    plot_inferred_vs_simulated_k_across_positions_num(32, axis[1][0], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "c\n", plot_labels=True) #, only_of_significants=True, use_empirical=True)
    plot_accuracy_vs_k_across_positions_num(32, axis[1][1], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "d\n",
                                            add_legend=False, plot_labels=True)  # , only_of_significants=True, use_empirical=True)
    plot_inferred_vs_simulated_k_across_positions_num(64, axis[2][0], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "e\n", plot_labels=False) #, only_of_significants=True, use_empirical=True)
    plot_accuracy_vs_k_across_positions_num(64, axis[2][1], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "f\n", add_legend=True,
                                            plot_labels=False) #, only_of_significants=True, use_empirical=True)

    # fig.text(0.53, -0.02, r"$k$", ha='center', fontdict={'size': 30})
    # fig.text(-0.02, 0.5, "Mean(" + r"$\^k$" + ")", va='center', rotation='vertical', fontdict={'size': 30})
    # fig.text(, 0.25, "Mean error (" + r"$\^k$" + ")", va='center', rotation='vertical', fontdict={'size': 30})
    fig.subplots_adjust()
    fig.tight_layout(pad=1)
    plt.savefig(output_path_4, bbox_inches='tight') #, transparent=True)
    plt.clf()

    # figure s4 - power of hyphy
    plt.grid(False)
    fig, axis = plt.subplots(nrows=1, ncols=2, sharex="all", sharey="all", figsize=[2 * 8.2 + 2, 7.58], frameon=True)
    plot_power_vs_k_across_posnum(32, axis[0], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "a\n",
                                  empirical=False, add_legend=True, linestyle='dashed')
    plot_power_vs_k_across_taxanum(300, axis[1], HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, "b\n",
                                   empirical=False, add_legend=True, plot_labels=True, linestyle='dashed')
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(output_path_5, bbox_inches='tight') #, transparent=True)
    plt.clf()

#####################################################################################


if __name__ == '__main__':

    matplotlib.rc('xtick', labelsize=22)
    matplotlib.rc('ytick', labelsize=22)
    matplotlib.rc('legend', frameon=False)
    matplotlib.rc('legend', fontsize=30)
    matplotlib.rc('figure', facecolor='white')
    matplotlib.rc('axes', facecolor='white')
    matplotlib.rc('axes', grid=False)
    matplotlib.rc('figure', edgecolor='k')
    matplotlib.rc('figure', frameon=True)
    plt.rcParams["axes.edgecolor"] = "lightgrey"
    plt.rcParams["axes.linewidth"] = 1.25

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
                            RELAXComboToDf[(tbl, mu, pi0, taxa_num, codon_positions_num, k_value)] = df
                        # report theoretical and empirical results
                        if not skip_combo:
                            print("combo: (tbl=", tbl, ", mu=", mu, ", pi0=", pi0, ", taxa_num=", taxa_num,
                                  ", positions_num=", codon_positions_num, ")")
                            report_theoretical_test_result(RELAXComboToDf, tbl, mu, pi0, taxa_num, codon_positions_num)
                            EmpiricalRELAXLRThresholds[
                                (tbl, mu, pi0, taxa_num, codon_positions_num)] = find_empirical_LR_threshold(
                                RELAXComboToDf, tbl, mu, pi0, taxa_num, codon_positions_num)

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
                                df, full_df = extract_RELAX_data(input_dir, tbl, mu, pi0, taxa_num, codon_positions_num,
                                                                 k_value)
                                full_df.to_csv(input_dir + "full_analysis.csv")
                                df.to_csv(input_dir + "clean_analysis.csv")
                            else:
                                full_df = pd.read_csv(input_dir + "full_analysis.csv")
                                df = pd.read_csv(input_dir + "clean_analysis.csv")

                            # insert the dataframe into a dictionary
                            MPComboToDf[(tbl, mu, pi0, taxa_num, codon_positions_num, k_value)] = df

                        # report theoretical and empirical results
                        if not skip_combo:
                            print("combo: (tbl=", tbl, ", mu=", mu, ", pi0=", pi0, ", taxa_num=",
                                  taxa_num, ", positions_num=", codon_positions_num, ")")
                            report_theoretical_test_result(MPComboToDf, tbl, mu, pi0, taxa_num, codon_positions_num)
                            EmpiricalMPRELAXLRThresholds[
                                (tbl, mu, pi0, taxa_num, codon_positions_num)] = find_empirical_LR_threshold(
                                MPComboToDf, tbl, mu, pi0, taxa_num, codon_positions_num)

    print("***********************************************\n\n")

    # extract the data from RELAX + simulated character histories executions
    print("**** Processing RELAX (HyPhy) executions given the maximum parsimony history ****")
    HyPhyMPComboToDf = dict()
    EmpiricalHyPhyMPRELAXLRThresholds = dict()
    for tbl in tbl_options:
        for mu in mu_options:
            for pi0 in pi0_options:
                for taxa_num in taxa_num_options:
                    for codon_positions_num in codon_positions_num_options:
                        skip_combo = False
                        for k_value in k_values_options:
                            # declare the visited combo
                            input_dir = R_HyPhy_MP_simulation_study_output_dir + "/tbl_" + str(tbl) + "_mu_" + str(
                                mu) + "_pi0_" + str(
                                pi0) + "_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(
                                taxa_num) + "_taxa/" + str(codon_positions_num) + "_codons/k_" + str(k_value) + "/"
                            if not os.path.exists(input_dir):
                                if k_value == 1:
                                    skip_combo = True
                                continue

                            # extract the results
                            if not os.path.exists(input_dir + "clean_analysis.csv"):
                                df = extract_hyphy_data(input_dir)
                                df.to_csv(input_dir + "clean_analysis.csv")
                            else:
                                df = pd.read_csv(input_dir + "clean_analysis.csv")
                            # df = extract_hyphy_data(input_dir)
                            # df.to_csv(input_dir + "clean_analysis.csv")

                            # insert the dataframe into a dictionary
                            HyPhyMPComboToDf[(tbl, mu, pi0, taxa_num, codon_positions_num, k_value)] = df

                        # report theoretical and empirical results
                        if not skip_combo:
                            print("combo: (tbl=", tbl, ", mu=", mu, ", pi0=", pi0, ", taxa_num=",
                                  taxa_num, ", positions_num=", codon_positions_num, ")")
                            report_theoretical_test_result(HyPhyMPComboToDf, tbl, mu, pi0, taxa_num, codon_positions_num)
                            EmpiricalHyPhyMPRELAXLRThresholds[ (tbl, mu, pi0, taxa_num, codon_positions_num)] = find_empirical_LR_threshold(
                                HyPhyMPComboToDf, tbl, mu, pi0, taxa_num, codon_positions_num)

    print("***********************************************\n\n")

    plot_power_and_FPR_assessment(TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, RELAXComboToDf,
                                  EmpiricalRELAXLRThresholds, HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds,
                                  figures_dir + "power_and_fpr_assessment_for_results.svg",
                                  figures_dir + "power_and_fpr_assessment_for_supp_1.svg",
                                  figures_dir + "power_and_fpr_assessment_for_supp_2.svg")

    plot_accuracy_analysis(TraitRELAXComboToDf, RELAXComboToDf, MPComboToDf, grid_data_path,
                           figures_dir + "accuracy_assessment_for_results.svg",
                           figures_dir + "accuracy_assessment_for_supp_1.svg",
                           figures_dir + "accuracy_assessment_for_supp_2.svg",
                           figures_dir + "accuracy_assessment_for_supp_3.svg")

    # plot_duration(TraitRELAXComboToDf, figures_dir + "duration_analysis.svg")

    # report_power_fpr_comparison(TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds)

    plot_comparison_to_hyphy(TraitRELAXComboToDf, EmpiricalTraitRELAXLRThresholds, HyPhyMPComboToDf, EmpiricalHyPhyMPRELAXLRThresholds, figures_dir + "hyphy_comparison_1.svg", figures_dir + "hyphy_comparison_2.svg", figures_dir + "hyphy_comparison_3.svg", figures_dir + "hyphy_comparison_4.svg", figures_dir + "hyphy_comparison_5.svg")
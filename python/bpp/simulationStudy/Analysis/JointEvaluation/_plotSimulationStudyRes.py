import os, re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats.distributions import chi2
sns.set_style('whitegrid')
from sklearn import metrics

# hardcoded data
colors = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]

def extract_data(input_dir):

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
        "Alternative model fitting.*?Character Log likelihood\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    alternative_sequence_logl_regex = re.compile(
        "Alternative model fitting.*?Sequence Log likelihood\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
    alternative_joint_logl_regex = re.compile("Alternative model fitting.*?Overall Log likelihood\.*\:\s*(-\d*\.?\d*)",
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
    null_k_regex = re.compile("Null model fitting.*?model 2.*?RELAX.k\.*\:\s*(\d*\.?\d*)",
                                     re.MULTILINE | re.DOTALL)
    alternative_parameter_regex_str = "Alternative model fitting.*?\.PARAMETER\.*\:\s*(\d*\.?\d*)"
    alternative_k_regex = re.compile("Alternative model fitting.*model 2.*?RELAX.k\.*\:\s*(\d*\.?\d*)",
                                     re.MULTILINE | re.DOTALL)

    # intialize dataframe
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

    # fill in the dataframe
    for file in os.listdir(input_dir):
        if ".OU" in file:
            with open(input_dir+file, "r") as input_file:
                content=input_file.read()
                values = dict()

                # extract the values for the non-parametric columns
                for colname in non_parametric_regex_t_colname.keys():
                    try:
                        regex = non_parametric_regex_t_colname[colname]
                        values[colname] = float(regex.search(content).group(1))
                    except:
                        print("failed to find " + colname + " match in " + input_dir+file + "-> setting as Nan")
                        values[colname] = float('nan')

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

                    if parameter != "k":

                        # try to extract the null
                        try:
                            regex_str = null_parameter_regex_str.replace("PARAMETER", parameter)
                            regex = re.compile(regex_str, re.MULTILINE | re.DOTALL)
                            values["null_" + parameter] = float(regex.search(content).group(1))
                        except:
                            print("failed to find null " + parameter + " match in " + input_dir + file + "-> setting as Nan")
                            values["null_" + parameter] = float('nan')

                        # try to extract the alternative
                        try:
                            regex_str = alternative_parameter_regex_str.replace("PARAMETER", parameter)
                            regex = re.compile(regex_str, re.MULTILINE | re.DOTALL)
                            for match in regex.finditer(content):
                                pass
                            values["alternative_" + parameter] = float(match.group(1))
                        except:
                            print("failed to find alternative " + parameter + " match in " + input_dir + file + "-> setting as Nan")
                            values["alternative_" + parameter] = float('nan')

                    else: # different handling for k

                        # try to extract null
                        try:
                            values["null_" + parameter] = float(null_k_regex.search(content).group(1))
                        except:
                            print("failed to find null " + parameter + " match in " + input_dir + file + "-> setting as Nan")
                            values["null_" + parameter] = float('nan')


                        # try to extract alternative
                        try:
                            values["alternative_" + parameter] = float(alternative_k_regex.search(content).group(1))
                        except:
                            print(
                                "failed to find alternative " + parameter + " match in " + input_dir + file + "-> setting as Nan")
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
                if values['initial_sequence_logl'] > values['alternative_sequence_logl']:
                    print("optimization bug expressed in " + input_dir + "for replicate " + values[
                        "replicate"] + ", therefore, dataset will not be included in the analysis")
                else:
                    clean_df = clean_df.append(values, ignore_index=True)

                df = df.append(values, ignore_index=True)

    if df.shape[0] * 0.5 > clean_df.shape[1]:
        print("PROBLEM! More than half the data lost in cleaning process")

    return clean_df, df


def plot_inference_quality(combo_to_df, taxa_num, codon_positions_num, output_path):
    # parameters to evaluate the inference quality of
    fig_parameter_names = {"mu": r"$\mu$",
                           "pi0": r"$\pi_0$",
                           "kappa": r"$\kappa$",
                           "omega0": r"$\omega_0$",
                           "omega1": r"$\omega_1$",
                           "omega2": r"$\omega_2$",
                           "p0": r"$p_0$",
                           "p1": r"$p_1$",
                           "k": "k"}
    boxprops = dict(linestyle='-', linewidth=1, color='k')
    medianprops = dict(linestyle='-', linewidth=0, color='r')
    meanprops = dict(linestyle='-', linewidth=1, color='b')
    fig, axis = plt.subplots(1, 5, sharex="row", sharey="col", figsize=[32,4.8])
    fig.add_subplot(111, frameon=False)
    plt.grid(False)

    # plot a figure with all the remaining parameters, for each one a boxplot of the inferred values in opt1 cycle
    for k in k_values_options:
        df = combo_to_df[(taxa_num, codon_positions_num, k)]
        subdf = df[["alternative_" + parameter for parameter in fig_parameter_names.keys()]]
        index = k_values_options.index(k)
        bp = subdf.boxplot(grid=False, showmeans=True, meanline=True, boxprops=boxprops,
                           medianprops=medianprops,
                           meanprops=meanprops, ax=axis[index], fontsize=8)
        true_values = []
        inferred_values = []
        for parameter in list(fig_parameter_names.keys()):
            true_values.append(df["simulated_" + parameter].unique()[0])
            inferred_values.append(df["alternative_" + parameter].unique()[0])

        axis[index].scatter(x=axis[index].get_xticks(), y=np.array(true_values), label="true value", edgecolors=None,
                            color="g", marker="s")
        axis[index].set_xticklabels([fig_parameter_names[param] for param in list(fig_parameter_names.keys())])
        axis[index].set_ylim([0, 2.5])
        axis[index].title.set_fontsize(8)
        axis[index].title.set_text("simulated k = " + str(k))


    handles, labels = axis[index].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.04, 1))
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel('Simulated values')
    plt.ylabel('Inferred values')
    plt.subplots_adjust()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.clf()


def plot_power_vs_k(combo_to_df, taxa_num, output_path):

    fig, axis = plt.subplots(nrows=1, ncols=1)
    plt.grid(False)

    # collect data
    positions_num_to_power = dict()
    k_values = [k for k in k_values_options if k != 1]
    for positions_num in codon_positions_num_options:
        positions_num_to_power[positions_num] = []
        for k in k_values:
            if k != 1:
                df = combo_to_df[(taxa_num, positions_num, k)]
                positions_num_to_power[positions_num].append(df[df.significant == True].shape[0] / df.shape[0])

    # plot a line for each positions_num option
    for positions_num in codon_positions_num_options:
        axis.plot(k_values, positions_num_to_power[positions_num], color=colors[codon_positions_num_options.index(positions_num)], label=str(positions_num) + " positions")
    axis.set_xlabel("Simulated value of k")
    axis.set_ylabel("Power")
    handles, labels = axis.get_legend_handles_labels()
    fig.legend(handles, labels, loc='best')
    plt.savefig(output_path)
    plt.clf()


def plot_inferred_vs_simulated_k(combo_to_df, taxa_num, output_path):
    fig, axis = plt.subplots(nrows=1, ncols=1)
    plt.grid(False)

    # plot the data
    for positions_num in codon_positions_num_options:
        inferred_k_values = []
        inferred_k_errors = []
        for k in k_values_options:
            df = combo_to_df[(taxa_num, positions_num, k)]
            inferred_k_values.append(df["alternative_k"].mean())
            inferred_k_errors.append(df["alternative_k"].std())
        axis.scatter(k_values_options, inferred_k_values, label=str(positions_num)+" positions", color=colors[codon_positions_num_options.index(positions_num)])
        axis.errorbar(k_values_options, inferred_k_values, yerr=inferred_k_errors, color=colors[codon_positions_num_options.index(positions_num)])

    axis.set_xlabel("Simulated value of k")
    axis.set_ylabel("Inferred value of k")
    handles, labels = axis.get_legend_handles_labels()
    fig.legend(handles, labels, loc='best')
    plt.savefig(output_path)
    plt.clf()


def plot_inferred_k_error(combo_to_df, taxa_num, output_path):
    fig, axis = plt.subplots(nrows=1, ncols=1)
    plt.grid(False)

    # plot the data
    for positions_num in codon_positions_num_options:
        inferred_k_errors = []
        for k in k_values_options:
            df = combo_to_df[(taxa_num, positions_num, k)]
            inferred_k_errors.append((abs(df["alternative_k"]-k)/k).mean()) # compute the mean of relative errors (regardless of direction of error
        axis.plot(k_values_options, inferred_k_errors, label=str(positions_num)+" positions", color=colors[codon_positions_num_options.index(positions_num)])

    axis.set_xlabel("Simulated value of k")
    axis.set_ylabel("Error")
    handles, labels = axis.get_legend_handles_labels()
    fig.legend(handles, labels, loc='best')
    plt.savefig(output_path)
    plt.clf()


def plot_power_vs_k_over_taxa_num(combo_to_df, positions_num, output_path):

    # get the power per alternative k for the set of positions_num
    labels = [str(taxa_num) + " species" for taxa_num in taxa_num_options] # label per taxa num
    hist_colors = [colors[i] for i in range(len(taxa_num_options))]
    data = [] # divided by the value of k
    for taxa_num in taxa_num_options:
        power_values = []
        for k in alternative_k_values:
            df = combo_to_df[(taxa_num, positions_num, k)]
            power_values.append(df[df.significant == True].shape[0] / df.shape[0])
        data.append(np.array(power_values))

    # plot the histogram
    fig, axis = plt.subplots(nrows=1, ncols=1)
    axis.hist(np.array(data), density=True, histtype='bar', color=hist_colors, label=labels)
    axis.legend(prop={'size': 10})
    axis.set_xlabel("Simulated k")
    axis.set_ylabel("Power")
    plt.savefig(output_path)
    plt.clf()


def get_ROC_data(combo_to_df, taxa_num, positions_num):

    # create a joint dataframe for all the k values with concatanation
    dataframes = [combo_to_df[(taxa_num, positions_num, k)] for k in k_values_options]
    joint_df = pd.concat(dataframes)

    # get a binary vector of the true labels according to the simulated values of k (1 for k != 1, else 0)
    true_labels = joint_df["simulated_k"] != 1

    # get a binary vector of the actual labels according to the significant records regardless of direction (1 if siginificant, else 0)
    inferred_labels = joint_df["significant"] == True

    # plot the ROC curve and return the AUC for the combination pf taxa num and positions num
    fpr, tpr, thresholds = metrics.roc_curve(true_labels, inferred_labels)
    roc_auc = metrics.roc_auc_score(true_labels, inferred_labels)

    return fpr, tpr, thresholds, roc_auc


def plot_ROC_curve(combo_to_df, taxa_num, positions_num, output_path):

    fpr, tpr, thresholds, roc_auc = get_ROC_data(combo_to_df, taxa_num, positions_num)
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.savefig(output_path)
    plt.clf()


def plot_AUC_by_taxa_num(combo_to_df, output_path):

    # collect the auc data
    labels = [str(positions_num) + " positions" for positions_num in codon_positions_num_options] # label per positions number
    hist_colors = [colors[i] for i in range(len(codon_positions_num_options))]
    data = []  # divided by taxa num over x-axis
    for taxa_num in taxa_num_options:
        auc_values = []
        for positions_num in codon_positions_num_options:
            fpr, tpr, thresholds, roc_auc = get_ROC_data(combo_to_df, taxa_num, positions_num)
            auc_values.append(roc_auc)
        data.append(np.array(auc_values))

    # plot the histogram
    fig, axis = plt.subplots(nrows=1, ncols=1)
    axis.hist(np.array(data), density=True, histtype='bar', color=hist_colors, label=labels)
    axis.legend(prop={'size': 10})
    axis.set_xlabel("Number of species")
    axis.set_ylabel("AUC")
    plt.savefig(output_path)
    plt.clf()



if __name__ == '__main__':

    # simulation_study_output_dir = ""
    # taxa_num_options = [16, 32, 64]
    # codon_positions_num_options = [300, 600, 1000]
    # k_values_options = [0.2, 0.5, 1, 1.6, 2]

    # process input from cmd
    parser = argparse.ArgumentParser(description='analyses output of RELAX simulation study')
    parser.add_argument('--simulation_study_output_dir', '-i', help='directiory that holds the simulation study output',
                        required=True)
    parser.add_argument('--taxa_num_options', '-tn', help='list of taxa numbers combinations used in the simulation study',
                        required=False, default=[16, 32, 64])
    parser.add_argument('--codon_positions_num_options', '-pn', help='list of the codon positions number combinations used in the simulaton study',
                        required=False, default=[300, 600, 1000])
    parser.add_argument('--k_values_options', '-k', help='list of k values combinations used in the simulation study', required=False, default=[0.2, 0.5, 1, 1.6, 2])

    args = parser.parse_args()
    simulation_study_output_dir = args.simulation_study_output_dir
    taxa_num_options = args.taxa_num_options
    if not type(taxa_num_options) == list:
        taxa_num_options = list(args.taxa_num_options.split(","))
    codon_positions_num_options = args.codon_positions_num_options
    if not type(codon_positions_num_options) == list:
        codon_positions_num_options = list(args.codon_positions_num_options.split(","))
    k_values_options = args.k_values_options
    if not type(k_values_options) == list:
        k_values_options = list(args.k_values_options.split(","))
    alternative_k_values = [k for k in k_values_options if k != 1]

    # process the data for each info combination
    combo_to_df = dict()
    for taxa_num in taxa_num_options:
        for codon_positions_num in codon_positions_num_options:
            for k_value in k_values_options:

                # declare the visited combo
                input_dir = simulation_study_output_dir + "/" + str(taxa_num) + "_taxa/" + str(codon_positions_num) + "_codons/k_" + str(k_value) + "/"
                if not os.path.exists(input_dir):
                    continue
                print("\nanalysing TraitRELAX output on the following simulations combo:\ntaxa_num=" + str(
                    taxa_num) + "\ncodon_positions_num=" + str(codon_positions_num) + "\nsimulated_k=" + str(
                    k_value) + "\n")

                # extract the results
                if not os.path.exists(input_dir + "full_analysis.csv") or not os.path.exists(input_dir + "clean_analysis.csv"):
                    df, full_df = extract_data(input_dir)
                    full_df.to_csv(input_dir + "full_analysis.csv")
                    df.to_csv(input_dir + "clean_analysis.csv")
                else:
                    full_df = pd.read_csv(input_dir + "full_analysis.csv")
                    df = pd.read_csv(input_dir + "clean_analysis.csv")

                # insert the dataframe into a dictionary
                combo_to_df[(taxa_num, codon_positions_num, k_value)] = df

                # declare various statistics of the analysed output
                print("mean number of optimization cycles: ", df["num_of_cycles"].mean())
                print("mean duration(hours): ", df["duration(hours)"].mean())

                # in case of null simulation, report the FDR
                if k_value == 1:
                    print("False positive rate: ", df[df.significant == True].shape[0] / df.shape[0])
                # in case of alternative simulation, report the overall power and the % of significant results in the correct direction
                else:
                    print("power: ", df[df.significant == True].shape[0] / df.shape[0])
                    print("power in the correct direction of selection intensity: ", df[df.correctDirection == 0][df.significant == True].shape[0]) / (df[df.significant == True].shape[0])

    # plot the figures for the simulation study
    figures_dir = simulation_study_output_dir + "figures/"
    if not os.path.exists(figures_dir):
        res = os.system("mkdir -p " + figures_dir)

    # plot a boxplot of inference quality per combo
    inference_quality_figures_dir = figures_dir + "inference_quality/"
    if not os.path.exists(inference_quality_figures_dir):
        res = os.system("mkdir -p " + inference_quality_figures_dir)
    for taxa_num in taxa_num_options:
        for codon_positions_num in codon_positions_num_options:
            output_path = inference_quality_figures_dir + str(taxa_num) + "_taxa_" + str(codon_positions_num) + "_codons.jpg"
            plot_inference_quality(combo_to_df, taxa_num, codon_positions_num, output_path)

    # plot the power vs. the simulated k, over different codon positions num, per taxa num
    power_vs_k_figures_dir = figures_dir + "power_vs_k/"
    if not os.path.exists(power_vs_k_figures_dir):
        res = os.system("mkdir -p " + power_vs_k_figures_dir)
    for taxa_num in taxa_num_options:
        output_path = power_vs_k_figures_dir + str(taxa_num) + "taxa.jpg"
        plot_power_vs_k(combo_to_df, taxa_num, output_path)

    # plot the mean inferred k (with std bar) vs. the simulated k over different codon positions num, per taxa num
    inferred_k_vs_simulated_k_figures_dir = figures_dir + "inferred_vs_simulated_k/"
    if not os.path.exists(inferred_k_vs_simulated_k_figures_dir):
        res = os.system("mkdir -p " + inferred_k_vs_simulated_k_figures_dir)
    for taxa_num in taxa_num_options:
        output_path = inferred_k_vs_simulated_k_figures_dir + str(taxa_num) + "_taxa.jpg"
        plot_inferred_vs_simulated_k(combo_to_df, taxa_num, output_path)

    # plot the mean relative error (with std bar) of the estimated k vs. the simulated k over different codon positions num, per taxa num
    inferred_k_error_figures_dir = figures_dir + "inferred_k_error/"
    if not os.path.exists(inferred_k_error_figures_dir):
        res = os.system("mkdir -p " + inferred_k_error_figures_dir)
    for taxa_num in taxa_num_options:
        output_path = inferred_k_error_figures_dir + str(taxa_num) + "_taxa.jpg"
        plot_inferred_k_error(combo_to_df, taxa_num, output_path)

    # plot a histogram of the power vs. the simulated k over the different taxa num, per codon_posotions_num option
    power_vs_k_over_taxa_num_figures_dir = figures_dir + "power_vs_k_over_taxa_num/"
    if not os.path.exists(power_vs_k_over_taxa_num_figures_dir):
        res = os.system("mkdir -p " + power_vs_k_over_taxa_num_figures_dir)
    for positions_num in codon_positions_num_options:
        output_path = power_vs_k_over_taxa_num_figures_dir + str(positions_num) + "_codons.jpg"
        plot_power_vs_k_over_taxa_num(combo_to_df, positions_num, output_path)

    # plot a ROC curve for each combo of taxa num and positions num over all the k values at once, add the AUC to each df while at it
    roc_figures_dir = figures_dir + "ROC/"
    if not os.path.exists(roc_figures_dir):
        res = os.system("mkdir -p " + roc_figures_dir)
    for taxa_num in taxa_num_options:
        for codon_positions_num in codon_positions_num_options:
            output_path = roc_figures_dir + str(taxa_num) + "_taxa_" + str(codon_positions_num) + "_codons.jpg"
            plot_ROC_curve(combo_to_df, taxa_num, codon_positions_num, output_path)

    # plot a histogram of the AUC vs. the taxa num over the different codon positions num options
    output_path = figures_dir + "AUC_vs_taxa_num_over_taxa_num.jpg"
    plot_AUC_by_taxa_num(combo_to_df, output_path)



# likelihood based plot - read txt. ith logl values - the last corresponds to the expected history + fix the one of the true history based on the single output
import re, os, argparse
import numpy as np
import seaborn as sns
sns.set_style('whitegrid')
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd
from matplotlib.lines import Line2D


###################################################################

# plot the mean overlap vs. mu
def plot_mappings_distance_vs_mu(ax, df, title):

    ax.set_frame_on(True)
    mappings_num = 1000
    ax.grid(False)

    # gather the data
    mappings_y_values = []
    mappings_yerr_values = []
    sampling_based_expected_y_values = []
    sampling_based_expected_yerr_values = []
    analytic_expected_y_values = []
    analytic_expected_yerr_values = []
    for mu in mu_options:
        data = df.loc[(df["mu"] == mu) & (df["mappings_num"] == mappings_num)]
        # mean
        mappings_y_values.append(np.mean(np.array(data["mapping_distance"])))
        sampling_based_expected_y_values.append(np.mean(np.array(data["sampling_based_expected_history_distance"])))
        analytic_expected_y_values.append(np.mean(np.array(data["analytic_expected_history_distance"])))

        # std
        mappings_yerr_values.append(np.std(np.array(data["mapping_distance"])))
        sampling_based_expected_yerr_values.append(np.std(np.array(data["sampling_based_expected_history_distance"])))
        analytic_expected_yerr_values.append(np.std(np.array(data["analytic_expected_history_distance"])))

    # plot the data
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.errorbar(x=mu_options, y=mappings_y_values, yerr=mappings_yerr_values, label='Exhaustive', lw=1, marker="s", color="black", markersize=16, linestyle='dashed')
    ax.errorbar(x=mu_options, y=sampling_based_expected_y_values, yerr=sampling_based_expected_yerr_values, label='Expected history (sampled)', lw=1, marker="^", color="black", markersize=16, linestyle='dashed')
    ax.errorbar(x=mu_options, y=analytic_expected_y_values, yerr=analytic_expected_yerr_values, label='Expected history (analytic)', lw=1, marker="o", color="black", markersize=16, linestyle='dashed')

    # add plot settings
    ax.set_xticks(mu_options)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_ylabel("Fraction of histories disagreement", fontdict={'family': 'sans-serif', 'size': 20})
    # ax.update_datalim()

# plot the percentile of expected overlap vs. mu for different number of mappings
def plot_expected_history_rank_among_mappings(ax, distances_df, title):

    ax.set_frame_on(True)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    ax.grid(False)

    x_values = mu_options
    sampling_based_y_values = dict()
    sampling_based_y_errors = dict()
    analytic_y_values = dict()
    analytic_y_errors = dict()

    # gather the data
    for mappings_num in num_of_mappings_options:
        sampling_based_y_values[mappings_num] = []
        sampling_based_y_errors[mappings_num] = []
        analytic_y_values[mappings_num] = []
        analytic_y_errors[mappings_num] = []
        for mu in mu_options:
            data = distances_df.loc[(distances_df["mu"] == mu) & (distances_df["mappings_num"] == mappings_num) & (distances_df["taxa_num"] == 32) & (distances_df["positions_num"] == 300)]
            replicates = data.replicate.unique()
            sampling_based_precentiles_data = []
            analytic_precentiles_data = []
            for replicate in replicates:
                replicate_data = data.loc[(data["replicate"] == replicate)]
                mappings_distances = list(replicate_data["mapping_distance"])
                mappings_overlaps = [1-distance for distance in mappings_distances]
                mappings_overlaps.sort()
                sampling_based_expected_mapping_distance = list(replicate_data["sampling_based_expected_history_distance"])[0]
                sampling_based_expected_mapping_overlap = 1 - sampling_based_expected_mapping_distance
                analytic_expected_mapping_distance = list(replicate_data["analytic_expected_history_distance"])[0]
                analytic_expected_mapping_overlap = 1 - analytic_expected_mapping_distance
                # find the percentile in which the expected_mapping_distance lies within the mappings distances
                sampling_based_expected_percentile = (len([overlap for overlap in mappings_overlaps if overlap < sampling_based_expected_mapping_overlap]))
                sampling_based_precentiles_data.append(sampling_based_expected_percentile)
                analytic_expected_percentile = (len([overlap for overlap in mappings_overlaps if overlap < analytic_expected_mapping_overlap]))
                analytic_precentiles_data.append(analytic_expected_percentile)
            sampling_based_y_values[mappings_num].append(np.mean(sampling_based_precentiles_data))
            sampling_based_y_errors[mappings_num].append(np.std(sampling_based_precentiles_data))
            analytic_y_values[mappings_num].append(np.mean(analytic_precentiles_data))
            analytic_y_errors[mappings_num].append(np.std(analytic_precentiles_data))

    # plot the data only for 1000 mappings
    ax.errorbar(x_values, sampling_based_y_values[1000], yerr=sampling_based_y_errors[1000], marker="^", lw=1, label="Expected history (sampled)", color="black", markersize=16, linestyle='dashed')
    ax.errorbar(x_values, analytic_y_values[1000], yerr=analytic_y_errors[1000], marker="o", lw=1, label="Expected history (analytic)", color="black", markersize=16, linestyle='dashed')

    ax.set_xticks(mu_options)
    ax.set_yticks([0, 200, 400, 600, 800, 1001])
    ax.set_yticklabels(["0", "200", "400", "600", "800", "1000"])
    ax.set_ylabel("rank (fraction of histories agreement) ", fontdict={'family': 'sans-serif', 'size': 20})
    # ax.update_datalim()

###################################################################

# extract from execution output files
def extract_logl_values(mu, taxa_num, positions_num, k, input_dir, jobs_dir):
    if os.path.exists(input_dir + "logls.csv"):
        df = pd.read_csv(input_dir + "logls.csv")
        if df.shape[0] < 2:
            print("no data in csv at: ", input_dir)
            exit(1)
    else:
        df = pd.DataFrame(
            columns=["mu", "taxa_num", "positions_num", "k", "replicate", "true_history_logl", "mappings_num",
                     "exhaustive_computation_logl"
                     "sampling_based_expected_history_logl", "sampling_based_expected_history_max_logl_diff",
                     "analytic_expected_history_logl", "mapping_order", "mapping_logl"])

        values = {"mu": mu, "taxa_num": taxa_num, "positions_num": positions_num, "k": k}
        replicate_regex = re.compile("replicate_(\d*)", re.MULTILINE | re.DOTALL)
        mappings_data_regex_str = "Analysis based on MAPPINGS_NUM mappings.*?Sequence log likelihood\.*\:\s*-\d*\.?\d*"
        character_logl_regex = re.compile("Character model log likelihood\:\s*\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
        mapping_sequence_logl_regex = re.compile("Initializing data structure.*?\n(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
        expected_sequence_logl_regex = re.compile("Computing log likelihood based on the expected history approximation.*?Sequence log likelihood\.*\:\s*(-\d*\.?\d*)", re.MULTILINE | re.DOTALL)
        for path in os.listdir(jobs_dir):
            if ".OU" in path:
                fullpath = jobs_dir+path
                with open(fullpath, "r") as infile:
                    content = infile.read()
                    rep = int(replicate_regex.search(content).group(1))
                    values["replicate"] = rep
                    try:
                        character_logl = float(character_logl_regex.search(content).group(1))
                        # get the true history logl
                        true_history_logl_path = input_dir + "replicate_" + str(rep) + "/traitrelax_result/histories_evaluation/true_history_likelihood.txt"
                        with open(true_history_logl_path, "r") as true_history_logl_file:
                            logl_str = true_history_logl_file.read()
                            logl_str = logl_str.replace("\n", "")
                            values["true_history_logl"] = -1 * float(logl_str) - character_logl
                        for num_of_mappings in num_of_mappings_options:
                            regex_str = mappings_data_regex_str.replace("MAPPINGS_NUM", str(num_of_mappings))
                            regex = re.compile(regex_str, re.MULTILINE | re.DOTALL)
                            mappings_data = regex.findall(content)[0]
                            values["mappings_num"] = num_of_mappings
                            logls = []
                            for match in mapping_sequence_logl_regex.finditer(mappings_data):
                                logls.append(match.group(1))
                            order = 1
                            values["sampling_based_expected_history_logl"] = float(expected_sequence_logl_regex.search(mappings_data).group(1))
                            for logl in logls:
                                values["mapping_order"] = order
                                order += 1
                                values["mapping_logl"] = logl
                                df = df.append(values, ignore_index=True)
                    except:
                        print("failed to extract logl data data for mu=", mu, ", replicate=", rep)
                        print("input path: ", fullpath)
                        exit(0)
                        # continue
                print("processed: ", fullpath)
                df.to_csv(input_dir + "logls.csv")
    return df

def plot_logl_distance_vs_mu(ax, df, title):

    ax.set_frame_on(True)

    # gather data: this plot is based on a single replicate
    mu_to_approximation_diffs = dict()
    for mu in mu_options:
        all_mappings_data = df.loc[(df["mu"] == mu)]

        # compute the diff from true history logl for each replicate
        replicates = all_mappings_data.replicate.unique()
        sampling_based_expected_histories_diffs = []
        exhaustive_diffs = []
        analytic_expected_histories_diffs = []
        for replicate in replicates:
            relevant_df = all_mappings_data.loc[(all_mappings_data["replicate"] == replicate)]
            relevant_record = relevant_df["true_history_logl"]
            true_history_logl = float(list(relevant_record)[0])
            sampling_based_expected_history_logl = list(relevant_df["sampling_based_expected_history_logl"])[0]
            sampling_based_expected_history_logl_diff = true_history_logl - sampling_based_expected_history_logl
            sampling_based_expected_histories_diffs.append(sampling_based_expected_history_logl_diff)
            analytic_expected_history_logl = list(relevant_df["analytic_expected_history_logl"])[0]
            analytic_expected_history_logl_diff = true_history_logl - analytic_expected_history_logl
            analytic_expected_histories_diffs.append(analytic_expected_history_logl_diff)
            exhaustive_computation_logl = list(relevant_df["exhaustive_computation_logl"])[0]
            exhaustive_computation_logl_diff = true_history_logl - exhaustive_computation_logl
            exhaustive_diffs.append(exhaustive_computation_logl_diff)
        mu_to_approximation_diffs[mu] = {"exhaustive": exhaustive_diffs, "sampling_based_expected": sampling_based_expected_histories_diffs, "analytic_expected": analytic_expected_histories_diffs}

    # normalize the distances
    max_distances_by_mu = []
    for mu in mu_options:
        max_dist_1 = np.max(mu_to_approximation_diffs[mu]["exhaustive"])
        max_dist_2 = np.max(mu_to_approximation_diffs[mu]["sampling_based_expected"])
        max_dist_3 = np.max(mu_to_approximation_diffs[mu]["analytic_expected"])
        max_distances_by_mu.append(np.max([max_dist_1, max_dist_2, max_dist_3]))
    max_dist = np.max(max_distances_by_mu)
    for mu in mu_options:
        for i in range(len(mu_to_approximation_diffs[mu]["exhaustive"])):
            mu_to_approximation_diffs[mu]["exhaustive"][i] = mu_to_approximation_diffs[mu]["exhaustive"][i] / max_dist
        for i in range(len(mu_to_approximation_diffs[mu]["sampling_based_expected"])):
            mu_to_approximation_diffs[mu]["sampling_based_expected"][i] = mu_to_approximation_diffs[mu]["sampling_based_expected"][i] / max_dist
        for i in range(len(mu_to_approximation_diffs[mu]["analytic_expected"])):
            mu_to_approximation_diffs[mu]["analytic_expected"][i] = mu_to_approximation_diffs[mu]["analytic_expected"][i] / max_dist

    # now, plot the results in error bars
    ax.grid(False)
    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 30}, loc='left')
    # plot results for exhaustive approximation
    y = []
    yerr = []
    for mu in mu_options:
        y.append(np.mean(mu_to_approximation_diffs[mu]["exhaustive"]))
        yerr.append(np.std(mu_to_approximation_diffs[mu]["exhaustive"]))
    ax.errorbar(mu_options, y, yerr=yerr, label='Exhaustive', lw=1, marker="s", color="black", markersize=16, linestyle='dashed')

    # plot results for sampling_based_expected approximation
    y = []
    yerr = []
    for mu in mu_options:
        y.append(np.mean(mu_to_approximation_diffs[mu]["sampling_based_expected"]))
        yerr.append(np.std(mu_to_approximation_diffs[mu]["sampling_based_expected"]))
    ax.errorbar(mu_options, y, yerr=yerr, label='Stochastic single history', lw=1, marker="^", color="black", markersize=16)

    # plot results for analytic_expected approximation
    y = []
    yerr = []
    for mu in mu_options:
        y.append(np.mean(mu_to_approximation_diffs[mu]["analytic_expected"]))
        yerr.append(np.std(mu_to_approximation_diffs[mu]["analytic_expected"]))
    ax.errorbar(mu_options, y, yerr=yerr, label='Analytic single history', lw=1, marker="o", color="black", markersize=16, linestyle='dashed')

    ax.set_xticks(mu_options)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_ylabel("Relative difference in log likelihood", fontdict={'family': 'sans-serif', 'size': 20})
    # ax.update_datalim()

def plot_mappings_logl_diff_from_true_distribution(ax, df, title, mu, mu_options, use_analytic=True):

    ax.set_title(title + r"                                    $\mu$=" + str(mu), fontdict={'family': 'sans-serif', 'size': 20}, loc='left')
    ax.grid(False)

    # gather data: this plot is based on a single replicate
    all_mappings_data = df.loc[(df["mu"] == mu)]

    # compute the diff from true history logl for each replicate
    replicates = all_mappings_data.replicate.unique()
    num_of_mappings_to_sampling_based_expected_histories_diffs = []
    num_of_mappings_to_histories_diffs = []
    num_of_mappings_to_analytic_expected_histories_diffs = []
    for replicate in replicates:
        relevant_df = all_mappings_data.loc[(all_mappings_data["replicate"] == replicate)]
        relevant_record = relevant_df["true_history_logl"]
        true_history_logl = float(list(relevant_record)[0])
        sampling_based_expected_history_logl = list(relevant_df["sampling_based_expected_history_logl"])[0]
        sampling_based_expected_history_logl_diff = true_history_logl - sampling_based_expected_history_logl
        num_of_mappings_to_sampling_based_expected_histories_diffs.append(sampling_based_expected_history_logl_diff)
        analytic_expected_history_logl = list(relevant_df["analytic_expected_history_logl"])[0]
        analytic_expected_history_logl_diff = true_history_logl - analytic_expected_history_logl
        num_of_mappings_to_analytic_expected_histories_diffs.append(analytic_expected_history_logl_diff)
        exhaustive_computation_logl = list(relevant_df["exhaustive_computation_logl"])[0]
        exhaustive_computation_logl_diff = true_history_logl - exhaustive_computation_logl
        num_of_mappings_to_histories_diffs.append(exhaustive_computation_logl_diff)

    # plot distribution of real model likelihood computation diffs
    sns.kdeplot(num_of_mappings_to_histories_diffs, ax=ax, color="steelblue", lw=3, shade=True)

    # plot distribution of expected mappings diffs
    if not use_analytic:
        sns.kdeplot(num_of_mappings_to_sampling_based_expected_histories_diffs, ax=ax, color="darkorange", lw=3, shade=True)

    # plot distribution of expected mappings diffs
    if use_analytic:
        sns.kdeplot(num_of_mappings_to_analytic_expected_histories_diffs, ax=ax, color="g", lw=3, shade=True)

    ax.set_xticks([-10, 0, 10, 20, 30])
    ax.set_yticks([0, 0.2, 0.4, 0.6])

    if mu == mu_options[0]:
        ax.set_ylabel('Frequency', fontdict={'family': 'sans-serif', 'size': 20})
    if mu == mu_options[1]:
        ax.set_xlabel('Difference in log likelihood\n(true-approximated)', fontdict={'family': 'sans-serif', 'size': 20})


def plot_mappings_logl_approx_diff_distribution(ax, df, title, mu, mu_options, use_analytic=True):

    ax.set_title(title, fontdict={'family': 'sans-serif', 'size': 20}, loc='left')
    ax.grid(False)

    # gather data: this plot is based on a single replicate
    all_mappings_data = df.loc[(df["mu"] == mu)]

    # compute the diff from true history logl for each replicate
    replicates = all_mappings_data.replicate.unique()
    approximations_logl_diff = []
    for replicate in replicates:
        relevant_df = all_mappings_data.loc[(all_mappings_data["replicate"] == replicate)]
        exhaustive_computation_logl = list(relevant_df["exhaustive_computation_logl"])[0]
        if not use_analytic:
            sampling_based_expected_history_logl = list(relevant_df["sampling_based_expected_history_logl"])[0]
            approximations_logl_diff.append(exhaustive_computation_logl - sampling_based_expected_history_logl)
        else:
            analytic_expected_history_logl = list(relevant_df["analytic_expected_history_logl"])[0]
            approximations_logl_diff.append(exhaustive_computation_logl - analytic_expected_history_logl)

    # plot distribution of approximations' likelihood computation diffs
    sns.kdeplot(approximations_logl_diff, ax=ax, color="grey", lw=3, shade=True)

    ax.set_xticks([-30, -20, -10, 0])
    ax.set_yticks([0, 0.1, 0.2])

    if mu == mu_options[1]:
        ax.set_xlabel('Difference in log likelihood\n(exhaustive-analytic)', fontdict={'family': 'sans-serif', 'size': 20})
    if mu == mu_options[0]:
        ax.set_ylabel('Frequency', fontdict={'family': 'sans-serif', 'size': 20})

###################################################################

def plot_full_figure(distances_df, logl_df, mu_options, output_path):

    # initialize a figure with 3 subplots
    plt.grid(False)
    fig, axis = plt.subplots(nrows=1, ncols=len(mu_options), sharex='all', frameon=True, figsize=[len(mu_options)*8+2, 7.58])
    plt.xticks(mu_options)

    # (a) - mean distance (logl) vs. mu
    plot_mappings_distance_vs_mu(axis[0], distances_df, "A\n")

    # (b) - mean distance (partition) vs. mu
    plot_logl_distance_vs_mu(axis[1], logl_df, "B\n")

    # (c) - rank (expected history) vs. mu
    plot_expected_history_rank_among_mappings(axis[2], distances_df, "C\n")

    # plot properties: shared xlabel (mu) and xticks (mu_options), shared legend: title: approximations, categories: exhaustive, expected history (sampled), expected history (analytic)
    axis[1].set_xlabel(r"$\mu$", fontdict={'family': 'sans-serif', 'size': 20})
    handles, labels = axis[0].get_legend_handles_labels()
    lgd = axis[0].legend(handles, labels, loc='upper left', prop={'size': 20}, frameon=False)
    fig.subplots_adjust()
    fig.tight_layout()
    plt.draw()  # necessary to render figure before saving
    plt.savefig(output_path, bbox_extra_artists=(lgd,), bbox_inches='tight', transparent=True)
    plt.clf()


def plot_alternative_figure(distances_df, logl_df, mu_options, output_path):

    # initialize a figure with 3 subplots
    plt.grid(False)
    fig, axis = plt.subplots(nrows=1, ncols=len(mu_options), figsize=[len(mu_options)*8+2, 1*7.58], frameon=True)
    plt.xticks(mu_options)

    # a
    plot_mappings_logl_diff_from_true_distribution(axis[0], logl_df, "A\n", mu_options[0], mu_options)

    # b
    plot_mappings_logl_diff_from_true_distribution(axis[1], logl_df, "B\n", mu_options[1], mu_options)

    # c
    plot_mappings_logl_diff_from_true_distribution(axis[2], logl_df, "C\n", mu_options[2], mu_options)

    custom_lines = []
    custom_names = []
    custom_lines.append(Line2D([0], [0], color="steelblue", lw=3))
    custom_names.append("Exhaustive approximation")
    custom_lines.append(Line2D([0], [0], color='g', lw=3))
    custom_names.append("Analytic approximation")
    lgd = axis[0][2].legend(custom_lines, custom_names, loc='right', bbox_to_anchor=(1, 0.9), prop={'size': 20}, frameon=False)
    fig.subplots_adjust()
    fig.tight_layout()
    plt.draw()  # necessary to render figure before saving
    plt.savefig(output_path, bbox_inches='tight', transparent=True, bbox_extra_artists=(lgd,))
    plt.clf()


###################################################################

def get_duration(str):
    days_regex = re.compile("Total execution time\: (\d*\.?\d*)d", re.MULTILINE | re.DOTALL)
    hours_regex = re.compile("Total execution time\:.*?(\d*\.?\d*)h", re.MULTILINE | re.DOTALL)
    minutes_regex = re.compile("Total execution time\:.*?(\d*\.?\d*)m", re.MULTILINE | re.DOTALL)
    seconds_regex = re.compile("Total execution time\:.*?(\d*\.?\d*)s", re.MULTILINE | re.DOTALL)
    duration = 0
    for match in days_regex.finditer(str):
        duration += float(match.group(1)) * 24
    for match in hours_regex.finditer(str):
        duration += float(match.group(1))
    for match in minutes_regex.finditer(str):
        duration += float(match.group(1)) * (1 / 60)
    for match in seconds_regex.finditer(str):
        duration += float(match.group(1)) * (1 / 60) * (1 / 60)
    return duration

def report_durations(dir):
    exhaustive_approximation_data_regex = re.compile(
        "Computing sequence log likelihoods given the different mappings.*?Total execution time.*?\n",
        re.MULTILINE | re.DOTALL)
    expected_history_approximation_regex = re.compile(
        "Computing log likelihood based on the expected history approximation.*?Total execution time.*?\n",
        re.MULTILINE | re.DOTALL)
    paths = []
    for path, subdirs, files in os.walk(dir):
        for name in files:
            paths.append(os.path.join(path, name))
    exhaustive_durations = []
    sampled_expected_durations = []
    analytic_expected_durations = []
    for path in paths:
        if ".OU" in path:
            with open(path, "r") as infile:
                content = infile.read()
            exhaustive_durations.append(get_duration(exhaustive_approximation_data_regex.search(content).group(0)))
            matches = expected_history_approximation_regex.findall(content)
            sampled_expected_durations.append(get_duration(matches[0]))
            analytic_expected_durations.append(get_duration(matches[1]))

    # for each approximation approach, report its mean, median, std, min and max
    print("approximation\tmean\tmedian\tstd\tmin\tmax")

    # exhaustive
    mean = np.mean(exhaustive_durations)
    median = np.median(exhaustive_durations)
    std = np.std(exhaustive_durations)
    min = np.min(exhaustive_durations)
    max = np.max(exhaustive_durations)
    print("exhaustive\t", mean, "\t", median, "\t", std, "\t", min, "\t", max)

    # sampled expected
    mean = np.mean(sampled_expected_durations)
    median = np.median(sampled_expected_durations)
    std = np.std(sampled_expected_durations)
    min = np.min(sampled_expected_durations)
    max = np.max(sampled_expected_durations)
    print("sampled_expected\t", mean, "\t", median, "\t", std, "\t", min, "\t", max)

    # analytic expected
    mean = np.mean(analytic_expected_durations)
    median = np.median(analytic_expected_durations)
    std = np.std(analytic_expected_durations)
    min = np.min(analytic_expected_durations)
    max = np.max(analytic_expected_durations)
    print("sampled_expected\t", mean, "\t", median, "\t", std, "\t", min, "\t", max)

###################################################################

if __name__ == '__main__':

    matplotlib.rc('xtick', labelsize=22)
    matplotlib.rc('ytick', labelsize=22)

    # process input from command line
    parser = argparse.ArgumentParser(
    description='Analyses the estimated expected histories based on the multiple histories approximation and the true history based computation')
    parser.add_argument('--input_dir', '-i', help='directory that holds the stdout of the histories evaluation jobs', required=True)
    parser.add_argument('--mu_options', '-mu', help='list of values of mu to include in the analysis', required=False, default=[1, 2, 4]) #, 8])
    parser.add_argument('--taxa_num_options', '-tn', help='list of taxa number values to include in the analysis', required=False, default=[32])
    parser.add_argument('--positions_num_options', '-pn', help='list of positions number values to include in the analysis', required=False, default=[300])
    parser.add_argument('--k_options', '-ko', help='list of k values to include in the analysis', required=False, default=[0.5])
    parser.add_argument('--num_of_replicates', '-rn', help='number of replicates to execute jobs on per combo', required=False, default=50)
    parser.add_argument('--num_of_mappings_options', '-nm', help='number of mappings combos used in the evaluation', required=False, default=[1000])
    parser.add_argument('--jobs_output_dir', '-err', help='directory of computations of likelihoods given histories ', required=True)
    parser.add_argument('--output_path', '-o', help='full pah to png output file with the approximation analysis figure', required=True)

    args = parser.parse_args()
    input_dir = args.input_dir

    mu_options = args.mu_options
    if not type(mu_options) == list:
        mu_options = mu_options.split(",")
        for i in range(len(mu_options)):
            try:
                mu_options[i] = int(mu_options[i])
            except Exception as e:
                mu_options[i] = float(mu_options[i])

    taxa_num_options = args.taxa_num_options
    if not type(taxa_num_options) == list:
        taxa_num_options = [int(tn) for tn in taxa_num_options.split(",")]

    positions_num_options = args.positions_num_options
    if not type(positions_num_options) == list:
        positions_num_options = [int(pn) for pn in positions_num_options.split(",")]

    k_options = args.k_options
    if not type(k_options) == list:
        k_options = k_options.split(",")
        for i in range(len(k_options)):
            try:
                k_options[i] = int(k_options[i])
            except Exception as e:
                k_options[i] = float(k_options[i])

    num_of_replicates = int(args.num_of_replicates)
    num_of_mappings_options = args.num_of_mappings_options
    jobs_output_dir = args.jobs_output_dir
    output_path = args.output_path

    # gather the histories data #
    if os.path.exists(input_dir + "distance_data.csv"):
        distances_df = pd.read_csv(input_dir + "distance_data.csv")
    else:
        dfs = []
        for mu in mu_options:
            for taxa_num in taxa_num_options:
                for positions_num in positions_num_options:
                    for k in k_options:
                        for rep in range(num_of_replicates):
                            for mappings_num in num_of_mappings_options:
                                data_input_dir = input_dir + "tbl_4_mu_" + str(mu) + "_pi0_0.5_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa_num) + "_taxa/" + str(positions_num) + "_codons/k_" + str(k) + "/replicate_" + str(rep) + "/traitrelax_result/histories_evaluation/"
                                csv_path = input_dir + "tbl_4_mu_" + str(mu) + "_pi0_0.5_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa_num) + "_taxa/" + str(positions_num) + "_codons/k_" + str(k) + "/replicate_" + str(rep) + "/traitrelax_result/histories_evaluation/" + str(mappings_num) + "_mappings_distances.csv"
                                try:
                                    dfs.append(pd.read_csv(csv_path))
                                except:
                                    print("failed to extract distance data for mu=", mu, ", replicate=", rep, "mappings num=", mappings_num)
                                    print("csv_path: ", csv_path)
                                    print("\n")
                                    continue
        distances_df = pd.concat(dfs)
        distances_df.to_csv(input_dir + "distance_data.csv")

    print("finished distance data processing")

    # gather the logl data #

    dfs = []
    for mu in mu_options:
        for taxa_num in taxa_num_options:
            for positions_num in positions_num_options:
                for k in k_options:
                    data_input_path = input_dir + "tbl_4_mu_" + str(
                        mu) + "_pi0_0.5_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(
                        taxa_num) + "_taxa/" + str(positions_num) + "_codons/k_" + str(k) + "/logl_analysis.csv"
                    df = pd.read_csv(data_input_path)
                    df.dropna()
                    dfs.append(df)
    logls_df = pd.concat(dfs)
    logls_df.drop_duplicates(keep="last")
    logls_df.to_csv(input_dir + "logl_data.csv")

    print("finished logl data processing")

    # plot the results
    #plot_full_figure(distances_df, logls_df, mu_options, output_path)
    plot_alternative_figure(distances_df, logls_df, mu_options, output_path)
    print("figure plotted to: ", output_path)

    # get_durations_statistics
    # report_durations(input_dir)





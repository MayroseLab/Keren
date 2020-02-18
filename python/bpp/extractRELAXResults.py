import re, os, argparse, math
import pandas as pd
from scipy.stats.distributions import chi2
import matplotlib

matplotlib.use('agg')
import seaborn as sns

sns.set_style('whitegrid')


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

def doLRT(null_logl, alternative_logl, df=1):
    LR = 2 * (alternative_logl - null_logl)
    pvalue = chi2.sf(LR, df)
    return LR, pvalue

def extract_data_by_hypothesis(str, hypothesis, dictionary):

    regex_strings = {hypothesis + "_p": "Fitting the " + hypothesis + " model.*?model 2.*?RELAX\.p\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete",
                     hypothesis + "_omega1": "Fitting the " + hypothesis + " model.*?model 2.*?RELAX\.omega1\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete",
                     hypothesis + "_omega2": "Fitting the " + hypothesis + " model.*?model 2.*?RELAX\.omega2\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete",
                     hypothesis + "_theta1": "Fitting the " + hypothesis + " model.*?model 2.*?RELAX\.theta1\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete",
                     hypothesis + "_theta2": "Fitting the " + hypothesis + " model.*?model 2.*?RELAX\.theta2\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete",
                     hypothesis + "_k": "Fitting the " + hypothesis + " model.*?model 2.*?RELAX\.k\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete"}

    if hypothesis == "null":
        regex_strings["null_logl"] = "Fitting the null model.*Current log likelihood\.*\:\s*(-\d*\.?\d*).*?Fitting the alternative model"
    else:
        regex_strings["alternative_logl"] = "Fitting the alternative model.*Current log likelihood\.*\:\s*(-\d*\.?\d*).*?iteraive optimzation complete"

    # extract the basic field
    for field in regex_strings.keys():
        if not (hypothesis == "null" and "_k" in field):
            regex = re.compile(regex_strings[field].replace("HYPOTHESIS", hypothesis), re.MULTILINE | re.DOTALL)
            try:
                dictionary[field] = float(regex.search(str).group(1))
            except:
                print("failed to extract ", field, " for dataset ", dictionary["dataset_id"])
                print("regex: ", regex.pattern)
                return 1
    # compute the induced parameters
    dictionary[hypothesis + "_omega0"] = dictionary[hypothesis + "_p"] * dictionary[hypothesis + "_omega1"]
    dictionary[hypothesis + "_p0"] = dictionary[hypothesis + "_theta1"]
    dictionary[hypothesis + "_p1"] = dictionary[hypothesis + "_theta2"] * (1 - dictionary[hypothesis + "_theta1"])
    return 0

def extract_real_data(input_dir):
    colnames = ["dataset_id", "null_p", "null_omega1", "null_omega2", "null_omega0", "null_theta1", "null_theta2",
                "null_p0", "null_p1", "null_logl", "alternative_p", "alternative_omega1", "alternative_omega2",
                "alternative_omega0", "alternative_theta1", "alternative_theta2", "alternative_p0",
                "alternative_p1", "alternative_logl", "LRT_statistic", "pvalue", "alternative_k"]

    df = pd.DataFrame(columns=colnames)
    for path in os.listdir(input_dir):
        input_path = input_dir + path
        if ".OU" in input_path:
            dictionary = dict()
            with open(input_path, "r") as input_file:
                content = input_file.read()

            # catch the dataset id
            dataset_id_regex = re.compile("Parsing file .*?([^\/]*?).bpp for options", re.MULTILINE | re.DOTALL)
            try:
                dictionary["dataset_id"] = dataset_id_regex.search(content).group(1)
                print("input_path: ", input_path)
                print("dataset id: ", dictionary["dataset_id"], "\n")
            except:
                print("failed to extract dataset_id for path ", input_path)
                return None

            # catch the MLEs of the null optimization
            extract_data_by_hypothesis(content, "null", dictionary)

            # catch the MLEs of the alternative optimization
            extract_data_by_hypothesis(content, "alternative", dictionary)

            # perform LRT and add to dictionary the reslut
            LRT_statistic, pvalue = doLRT(dictionary["null_logl"], dictionary["alternative_logl"])
            dictionary["LRT_statistic"] = LRT_statistic
            dictionary["pvalue"] = pvalue

            df = df.append(dictionary, ignore_index=True)

    # write data to output file
    return df


if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='Extract the results of RELAX analysis in Bio++ on real datasets')
    parser.add_argument('--input_dir', '-i', help='directory that holds the output files by relax to be analyzed',
                        required=True)
    parser.add_argument('--output_dir', '-o', help='directory that will hold the csv output of the analysis',
                        required=True)

    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir

    # process data
    df = extract_real_data(input_dir)

    df.to_csv(output_dir + "res.csv")


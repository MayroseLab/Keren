import re, os, argparse
import pandas as pd
from scipy.stats.distributions import chi2
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')
# HARDCODED - the paraemters to extract
colnames = ['dataset_id','null_logl','alternative_logl','null_mu','alternative_mu','null_pi0','alternative_pi0','null_kappa','alternative_kappa','null_p','alternative_p','null_omega1','alternative_omega1','null_omega2','alternative_omega2','null_theta1','alternative_theta1','null_theta2','alternative_theta2','null_k','alternative_k','null_omega0','alternative_omega0','null_p0','alternative_p0','null_p1','alternative_p1',"LRT_statistic",'pvalue']

def doLRT(null_logl, alternative_logl, df=1):
    LR = 2 * (alternative_logl - null_logl)
    pvalue = chi2.sf(LR, df)
    return LR, pvalue

def extract_traitrelax_parameters(input_path):

    with open(input_path, "r") as infile:
        content = infile.read()

    dictionary = dict()
    dataset_id_regex = re.compile("Parsing file .*?([^\/]*?).bpp for options", re.MULTILINE | re.DOTALL)
    dictionary["dataset_id"] = dataset_id_regex.search(content).group(1)

    print("input_path: ", input_path)
    print("dataset: ", dictionary["dataset_id"])

    regex_strings = {"null_logl": "Null model fitting.*?Overall Log likelihood\.*?\s*\:\s*(-\d*\.?\d*)",
                     "null_kappa": "Null model fitting.*?RELAX.kappa_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "null_p": "Null model fitting.*?RELAX.p_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "null_omega1": "Null model fitting.*?RELAX.omega1_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "null_omega2": "Null model fitting.*?RELAX.omega2_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "null_theta1": "Null model fitting.*?RELAX.theta1_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "null_theta2": "Null model fitting.*?RELAX.theta2_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "null_k": "Null model fitting.*?RELAX.k_2\.*?\s*\:\s*(\d*\.?\d*)",
                     "null_mu": "Null model fitting.*?TwoParameterBinary\.mu\.*?\s*\:\s*(\d*\.?\d*)",
                     "null_pi0": "Null model fitting.*?TwoParameterBinary\.pi0\.*?\s*\:\s*(\d*\.?\d*)",
                     "alternative_logl": "Alternative model fitting.*Overall Log likelihood\.*?\s*\:\s*(-\d*\.?\d*)",
                     "alternative_kappa": "Alternative model fitting.*RELAX.kappa_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "alternative_p": "Null model fitting.*?RELAX.p_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "alternative_omega1": "Alternative model fitting.*RELAX.omega1_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "alternative_omega2": "Alternative model fitting.*RELAX.omega2_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "alternative_theta1": "Alternative model fitting.*RELAX.theta1_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "alternative_theta2": "Alternative model fitting.*RELAX.theta2_1\.*?\s*\:\s*(\d*\.?\d*)",
                     "alternative_k": "Alternative model fitting.*RELAX.k_2\.*?\s*\:\s*(\d*\.?\d*)",
                     "alternative_mu": "Alternative model fitting.*TwoParameterBinary\.mu\.*?\s*\:\s*(\d*\.?\d*)",
                     "alternative_pi0": "Alternative model fitting.*TwoParameterBinary\.pi0\.*?\s*\:\s*(\d*\.?\d*)"}
    # extract the basic field
    for field in regex_strings.keys():
        try:
            regex = re.compile(regex_strings[field], re.MULTILINE | re.DOTALL)
            dictionary[field] = regex.search(content).group(1)
        except:
            print("failed to extract ", field, " for dataset ", dictionary["dataset_id"])
            print("regex: ", regex_strings[field])
            print("function: extract_traitrelax_parameters")
            exit(1)
    # compute the induced parameters
    dictionary["null_omega0"] = str(float(dictionary["null_p"]) * float(dictionary["null_omega1"]))
    dictionary["null_p0"] = dictionary["null_theta1"]
    dictionary["null_p1"] = str(float(dictionary["null_theta2"]) * (1 - float(dictionary["null_theta1"])))
    dictionary["alternative_omega0"] = str(float(dictionary["alternative_p"]) * float(dictionary["alternative_omega1"]))
    dictionary["alternative_p0"] = dictionary["alternative_theta1"]
    dictionary["alternative_p1"] = str(float(dictionary["alternative_theta2"]) * (1 - float(dictionary["alternative_theta1"])))

    # LR and pvalue
    dictionary["LRT_statistic"] = 2 * (float(dictionary["alternative_logl"]) - float(dictionary["null_logl"]))
    dictionary["pvalue"] = chi2.sf(dictionary["LRT_statistic"], 1)  # 1 degree of freedom is the diff in num of parameters between alternative and null models

    return dictionary


if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='Extract the results of RELAX analysis in Bio++ on real datasets')
    parser.add_argument('--input_dir', '-i', help='directory that holds the output files by relax to be analyzed',
                        required=True)
    parser.add_argument('--output_dir', '-o', help='directory that will hold the csv output of the analysis', required=True)

    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir

    # process data
    df = pd.DataFrame(columns=colnames)
    for input_path in os.listdir(input_dir):
        if ".OU" in input_path:
            data = extract_traitrelax_parameters(input_dir + input_path)
            df = df.append(data, ignore_index=True)
    # write data to output file
    df.to_csv(output_dir+"res.csv")


import argparse, os, re
import pandas as pd
import numpy as np


def extract_data(input_path):

    with open(input_path, "r") as infile:
        content = infile.read()
    values = dict()

    dataset_id_regex = re.compile("Parsing file.*?([^/]*?)\.bpp", re.M | re.DOTALL)
    values["dataset_id"] = dataset_id_regex.search(content).group(1)

    parameters = ["1_Full.theta", "1_Full.theta1", "1_Full.theta2", "2_Full.theta", "2_Full.theta1", "2_Full.theta2",
               "3_Full.theta", "3_Full.theta1", "3_Full.theta2", "kappa", "p", "omega1", "omega2", "theta1", "theta2", "k"]

    regex_str = "Fitting the alternative model.*model 2.*?RELAX\.<parameter>\.*\:\s*(\d*\.?\d*)"

    for parameter in parameters:
        model_num = 1
        if "_2" in parameter:
            model_num = 2
        parameter_name = re.sub("_\d", "", parameter)
        parameter_regex_str = regex_str.replace("<model_num>", str(model_num))
        parameter_regex_str = parameter_regex_str.replace("<parameter>", parameter_name)
        parameter_regex = re.compile(parameter_regex_str, re.MULTILINE | re.DOTALL)

        try:
            parameter_value = float(parameter_regex.search(content).group(1))
            values[parameter] = parameter_value
        except:
            print("input path: ", input_path)
            print("failed to extract parameter ", parameter_name, " in model ", model_num)
            print("with regex: ", parameter_regex.pattern)
            return 1

    return values


def estimate_relax_parameters(record):

    template = '''model1 = RELAX(kappa=<kappa>,p=<p>,omega1=<omega1>,omega2=<omega2>,k=1,theta1=<theta1>,theta2=<theta2>,frequencies=F3X4,1_Full.theta=<1_Full.theta>,1_Full.theta1=<1_Full.theta1>,1_Full.theta2=<1_Full.theta2>,2_Full.theta=<2_Full.theta>,2_Full.theta1=<2_Full.theta1>,2_Full.theta2=<2_Full.theta2>,3_Full.theta=<3_Full.theta>,3_Full.theta1=<3_Full.theta1>,3_Full.theta2=<3_Full.theta2>)
model2 = RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,frequencies=F3X4,1_Full.theta=RELAX.1_Full.theta_1,1_Full.theta1=RELAX.1_Full.theta1_1,1_Full.theta2=RELAX.1_Full.theta2_1,2_Full.theta=RELAX.2_Full.theta_1,2_Full.theta1=RELAX.2_Full.theta1_1,2_Full.theta2=RELAX.2_Full.theta2_1,3_Full.theta=RELAX.3_Full.theta_1,3_Full.theta1=RELAX.3_Full.theta1_1,3_Full.theta2=RELAX.3_Full.theta2_1,k=<k>)
'''
    # first handle parameters which translate directroy to RELAX parameters
    relax_info = template.replace("<1_Full.theta>", str(float(record["1_Full.theta"])))
    relax_info = relax_info.replace("<1_Full.theta1>", str(float(record["1_Full.theta1"])))
    relax_info = relax_info.replace("<1_Full.theta2>", str(float(record["1_Full.theta2"])))
    relax_info= relax_info.replace("<2_Full.theta>", str(float(record["2_Full.theta"])))
    relax_info= relax_info.replace("<2_Full.theta1>", str(float(record["2_Full.theta1"])))
    relax_info= relax_info.replace("<2_Full.theta2>", str(float(record["2_Full.theta2"])))
    relax_info= relax_info.replace("<3_Full.theta>", str(float(record["3_Full.theta"])))
    relax_info= relax_info.replace("<3_Full.theta1>", str(float(record["3_Full.theta1"])))
    relax_info= relax_info.replace("<3_Full.theta2>", str(float(record["3_Full.theta2"])))
    relax_info = relax_info.replace("<kappa>", str(float(record["kappa"])))
    relax_info = relax_info.replace("<p>", str(float(record["p"])))
    relax_info = relax_info.replace("<omega1>", str(float(record["omega1"])))
    relax_info = relax_info.replace("<omega2>", str(float(record["omega2"])))
    relax_info = relax_info.replace("<theta1>", str(float(record["theta1"])))
    relax_info = relax_info.replace("<theta2>", str(float(record["theta2"])))
    relax_info = relax_info.replace("<k>", str(float(record["k"])))
    return relax_info


if __name__ == '__main__':

    # parse input from cmd
    parser = argparse.ArgumentParser(
        description='Extract the results of RELAX analysis in Bio++ on real datasets')
    parser.add_argument('--input_dir', '-i', help='directory that holds the output files of the RELAX fitting',
                        required=True)
    parser.add_argument('--traitrelax_param_dir', '-r', help='directory to the traitrelax parameter files that will be used as templates for the new parameter files', required=True)
    parser.add_argument('--output_dir', '-o', help='directory that will hold the TraitRELAX parameter files that were inspired by the RELAX fitting', required=True)

    args = parser.parse_args()
    input_dir = args.input_dir
    traitrelax_param_dir = args.traitrelax_param_dir
    output_dir = args.output_dir

    # create a dataframe with the RELAX fitting result
    df = pd.DataFrame(columns=["dataset_id", "1_Full.theta", "1_Full.theta1", "1_Full.theta2", "2_Full.theta", "2_Full.theta1", "2_Full.theta2", "3_Full.theta", "3_Full.theta1", "3_Full.theta2", "kappa", "p", "omega1", "omega2", "theta1", "theta2", "k"])
    for path in os.listdir(input_dir):
        if ".OU" in path:
            values = extract_data(input_dir+path)
            if values != 1:
                df = df.append(values, ignore_index=True)

    dataset_id_regex = re.compile("(.*?).bpp")
    for path in os.listdir(traitrelax_param_dir):
        dataset_id = dataset_id_regex.search(path).group(1)
        row = df.loc[df["dataset_id"] == dataset_id]
        relax_parameters_str = ""
        if row.shape[0] > 0:
            relax_parameters_str = estimate_relax_parameters(row)
        orig_relax_param_path = traitrelax_param_dir + dataset_id + ".bpp"
        orig_init_template = re.compile("model1 = RELAX.*model2 = RELAX\(.*?\)", re.MULTILINE | re.DOTALL)
        new_relax_param_path = output_dir + dataset_id + ".bpp"
        with open(orig_relax_param_path, "r") as infile:
            content = infile.read()
        if relax_parameters_str != "":
            content = re.sub(orig_init_template, relax_parameters_str, content)
        with open(new_relax_param_path, "w") as outfile:
            outfile.write(content)



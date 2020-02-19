import argparse, os, re, sys
from Bio import AlignIO

sys.path.append("/groups/itay_mayrose/halabikeren/myScripts/python/")
from utils.createJobFile import set_job_env, create_job_file


def extract_traitrelax_parameters(content, hypothesis, dictionary):
    regex_strings = dict()
    if hypothesis == "alternative":
        regex_strings[
            "logl"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?Overall Log likelihood\.*\:\s*(-\d*\.?\d*)"
        regex_strings[
            "1_Full.theta"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.1_Full.theta_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "1_Full.theta1"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.1_Full.theta1_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "1_Full.theta2"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.1_Full.theta2_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "2_Full.theta"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.2_Full.theta_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "2_Full.theta1"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.2_Full.theta1_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "2_Full.theta2"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.2_Full.theta2_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "3_Full.theta"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.3_Full.theta_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "3_Full.theta1"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.3_Full.theta1_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "3_Full.theta2"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.3_Full.theta2_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "kappa"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.kappa_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "p"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.p_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "omega1"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.omega1_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "omega2"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.omega2_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "theta1"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.theta1_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "theta2"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.theta2_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "k"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.k_2\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "mu"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.mu\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "pi0"] = "\*\*\*\*\s*Alternative model likelihood after optimization\s*\*\*\*\*.*?\.pi0\.*\:\s*(\d*\.?\d*)"
    else:
        regex_strings[
            "logl"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\* .*?Overall Log likelihood\.*\:\s*(-\d*\.?\d*)"
        regex_strings[
            "1_Full.theta"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?\.1_Full.theta_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "1_Full.theta1"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?\.1_Full.theta1_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "1_Full.theta2"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?\.1_Full.theta2_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "2_Full.theta"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?\.2_Full.theta_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "2_Full.theta1"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?\.2_Full.theta1_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "2_Full.theta2"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?\.2_Full.theta2_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "3_Full.theta"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?\.3_Full.theta_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "3_Full.theta1"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?\.3_Full.theta1_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "3_Full.theta2"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?\.3_Full.theta2_1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "kappa"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?model 2.*?\.kappa\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "p"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?model 2.*?\.p\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "omega1"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?model 2.*?\.omega1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "omega2"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?model 2.*?\.omega2\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "theta1"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?model 2.*?\.theta1\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "theta2"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?model 2.*?\.theta2\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "k"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\* .*?model 2.*?\.k\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "mu"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?\.mu\.*\:\s*(\d*\.?\d*)"
        regex_strings[
            "pi0"] = "\*\*\*\*\s*Null model fitting\s*\*\*\*\*.*?\*\*\s*Model parameters after sequence model optimizaiton\s*\*\*.*?\.pi0\.*\:\s*(\d*\.?\d*)"
    regex_strings["scaling_factor"] = "Tree scaled by\.*\: \s*(\d*\.?\d*)"

    null_section_regex = re.compile("Null model fitting(.*?)Alternative model fitting", re.MULTILINE | re.DOTALL)
    alternative_section_regex = re.compile("Alternative model fitting(.*?)", re.MULTILINE | re.DOTALL)
    # extract the basic field
    for field in regex_strings.keys():
        regex = re.compile(regex_strings[field], re.MULTILINE | re.DOTALL)
        if field != "scaling_factor":
            try:
                dictionary[field] = regex.search(content).group(1)
            except:
                print("failed to extract ", field, " for dataset ", dictionary["dataset_id"])
                print("regex: ", regex_strings[field])
                print("function: extract_traitrelax_parameters")
                return 1
        else:
            section = null_section_regex.search(content).group(1)
            if hypothesis == "alternative":
                section = alternative_section_regex.search(content).group(1)
            scaling_factor = 1
            for match in regex.finditer(section):
                scaling_factor = scaling_factor * float(match.group(1))
            dictionary[field] = str(scaling_factor)
    # compute the induced parameters
    dictionary["omega0"] = str(float(dictionary["p"]) * float(dictionary["omega1"]))
    dictionary["p0"] = dictionary["theta1"]
    dictionary["p1"] = str(float(dictionary["theta2"]) * (1 - float(dictionary["theta1"])))
    return 0


def extract_relax_parameters(content, dictionary):
    regex_strings = {"logl": "Fitting the null model.*?Log likelihood\.*?\:\s*(-\d*\.?\d*)",
                     "1_Full.theta": "Fitting the null model.*model 2.*?\.1_Full\.theta\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "1_Full.theta1": "Fitting the null model.*model 2.*?\.1_Full\.theta1\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "1_Full.theta2": "Fitting the null model.*model 2.*?\.1_Full\.theta2\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "2_Full.theta": "Fitting the null model.*model 2.*?\.2_Full\.theta\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "2_Full.theta1": "Fitting the null model.*model 2.*?\.2_Full\.theta1\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "2_Full.theta2": "Fitting the null model.*model 2.*?\.2_Full\.theta2\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "3_Full.theta": "Fitting the null model.*model 2.*?\.3_Full\.theta\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "3_Full.theta1": "Fitting the null model.*model 2.*?\.3_Full\.theta1\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "3_Full.theta2": "Fitting the null model.*model 2.*?\.3_Full\.theta2\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "kappa": "Fitting the null model.*model 2.*?\.kappa\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "p": "Fitting the null model.*model 2.*?\.p\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "omega1": "Fitting the null model.*model 2.*?\.omega1\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "omega2": "Fitting the null model.*model 2.*?\.omega2\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "theta1": "Fitting the null model.*model 2.*?\.theta1\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "theta2": "Fitting the null model.*model 2.*?\.theta2\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "k": "Fitting the null model.*model 2.*?\.k\.*\:\s*(\d*\.?\d*).*?iteraive optimzation complete.*?Fitting the alternative model",
                     "scaling_factor": "Tree scaled by\.*\: \s*(\d*\.?\d*)"}

    null_section_regex = re.compile("Fitting the null model(.*?)Fitting the alternative model",
                                    re.MULTILINE | re.DOTALL)
    # extract the basic field
    for field in regex_strings.keys():
        if field != "scaling_factor":
            regex = re.compile(regex_strings[field], re.MULTILINE | re.DOTALL)
            try:
                dictionary[field] = regex.search(content).group(1)
            except:
                print("failed to extract ", field, " for dataset ", dictionary["dataset_id"])
                print("regex: ", regex.pattern)
                print("function: extract_relax_parameters")
                return 1
        else:
            section = null_section_regex.search(content).group(1)
            scaling_factor = 1
            scaling_regex = re.compile(regex_strings[field], re.MULTILINE | re.DOTALL)
            for match in scaling_regex.finditer(section):
                scaling_factor = scaling_factor * float(match.group(1))
            dictionary[field] = str(scaling_factor)
    # compute the induced parameters
    dictionary["omega0"] = str(float(dictionary["p"]) * float(dictionary["omega1"]))
    dictionary["p0"] = dictionary["theta1"]
    dictionary["p1"] = str(float(dictionary["theta2"]) * (1 - float(dictionary["theta1"])))
    return 0


def parse_results(input_path, method="traitrelax", hypothesis_to_use="null"):
    dictionary = dict()
    with open(input_path, "r") as input_file:
        content = input_file.read()

    # catch the dataset id
    dataset_id_regex = re.compile("Parsing file .*?([^\/]*?).bpp for options", re.MULTILINE | re.DOTALL)
    try:
        dictionary["dataset_id"] = dataset_id_regex.search(content).group(1)
    except:
        print("failed to extract dataset_id for path ", input_path)
        return None

    # catch the MLEs of the null optimization
    if method == "relax":
        extract_relax_parameters(content, dictionary)
    else:
        extract_traitrelax_parameters(content, hypothesis_to_use, dictionary)

    # write data to output file
    return dictionary


def get_initial_parameters(input_path, method="traitrelax"):
    regex_data = {"1_Full.theta": re.compile("1_Full.theta=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "1_Full.theta1": re.compile("1_Full.theta1=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "1_Full.theta2": re.compile("1_Full.theta2=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "2_Full.theta": re.compile("2_Full.theta=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "2_Full.theta1": re.compile("2_Full.theta1=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "2_Full.theta2": re.compile("2_Full.theta2=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "3_Full.theta": re.compile("3_Full.theta=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "3_Full.theta1": re.compile("3_Full.theta1=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "3_Full.theta2": re.compile("3_Full.theta2=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "kappa": re.compile("kappa=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "p": re.compile("p=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "omega1": re.compile("omega1=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "omega2": re.compile("omega2=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "theta1": re.compile("theta1=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "theta2": re.compile("theta2=(\d*\.?\d*)", re.MULTILINE | re.DOTALL),
                  "k": re.compile("model2.*?k=(\d*\.?\d*)", re.MULTILINE | re.DOTALL)}

    if method == "traitrelax":
        regex_data["mu"] = re.compile("character_model.mu\s*=\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL)
        regex_data["pi0"] = re.compile("character_model.pi0\s*=\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL)

    with open(input_path, "r") as infile:
        content = infile.read()

    parameters = dict()
    for field in regex_data:
        try:
            parameters[field] = regex_data[field].search(content).group(1)
        except:
            print("failed to extract ", field)
            print("regex: ", regex_data[field].pattern)
            print("input_path: ", input_path)
            print("function: get_initial_parameters")
            return 1
    parameters["omega0"] = str(float(parameters["p"]) * float(parameters["omega1"]))
    parameters["p0"] = parameters["theta1"]
    parameters["p1"] = str(float(parameters["theta2"]) * (1 - float(parameters["theta1"])))

    return parameters


if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='Creates simulations based on inferred parameters and executes inference program on them. This pipeline is comparible both with RELAX and with TraitRELAX')
    parser.add_argument('--program_path', '-p', help='path to the program file', required=False,
                        default="/groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite/traitrelax")
    parser.add_argument('--input_path', '-i',
                        help='path to the parameter file given as input to the program to create the inference result. Form here, we will extract the initial parameters for the parameter files of the simulations',
                        required=True)
    parser.add_argument('--inference_path', '-r', help='path to the output file of the program on the data',
                        required=True)
    parser.add_argument('--simulations_dir', '-s', help='directory to create simulations environment in', required=True)
    parser.add_argument('--data_path', '-d', help='path to the sequence alignment of the data', required=True)
    parser.add_argument('--tree_path', '-t', help='path to the tree path of the data in data_path', required=True)
    parser.add_argument('--history_path', '-hist',
                        help='path to the labeled history. use in case of parametric bootstrapping on RELAX',
                        required=False, default="")
    parser.add_argument('--bootstrap_on_char_only', '-c',
                        help='boolean indicating weather only character data should be simulated or also sequence data',
                        required=False, default=False)
    parser.add_argument('--set_advanced_parameters', '-sc',
                        help='set argument in the parameter files that allows scaling of the tree by sequence model',
                        required=False, default=True)
    parser.add_argument('--num_of_replicates', '-rep', help='number of replicates to simulate and infer results on',
                        required=False, default=200)
    parser.add_argument('--queue', '-q', help='name of the queue to send inference jobs to', required=False,
                        default="itaym1")
    parser.add_argument('--priority', '-pr', help='priority of inference jobs', required=False, default="0")
    parser.add_argument('--without_trait_simulations', '-nt',
                        help='boolean indicating weather trait data should be simulated or not', required=False,
                        default=True)
    parser.add_argument('--character_data_path', '-cd',
                        help='path to the character data file to be given as input to traitrelax parameter file in case trait data is not simulated',
                        required=False, default="")
    parser.add_argument('--run_inference', '-ri', help='indicator of weather inference should be executed now or later',
                        required=False, default=False)

    args = parser.parse_args()
    program_path = args.program_path
    input_path = args.input_path
    inference_path = args.inference_path
    simulations_dir = args.simulations_dir
    data_path = args.data_path
    bootstrap_on_char_only = bool(int(args.bootstrap_on_char_only))
    set_advanced_parameters = bool(int(args.set_advanced_parameters))
    tree_path = args.tree_path
    history_path = args.history_path
    num_of_replicates = int(args.num_of_replicates)
    queue = args.queue
    priority = args.priority
    without_trait_simulations = bool(int(args.without_trait_simulations))
    character_data_path = args.character_data_path
    run_inference = bool(int(args.run_inference))

    # extract the program name from the program path
    method = "relax"
    if "traitrelax" in program_path:
        method = "traitrelax"

    # parse the inference results from the inference_path
    inferred_parameters = parse_results(inference_path, method=method, hypothesis_to_use="null")

    # parse data size from data_path (assumes the input file is in fasta format)
    with open(data_path, "r") as datafile:
        alignment = AlignIO.read(datafile, "fasta")
    alignment_length = alignment.get_alignment_length()
    sequences_names = []
    for record in alignment:
        sequences_names.append(record.id)
    sequences_num = len(sequences_names)

    # extract the initial parameters for the simulation parameter files - to ensure similar initial conditions for the original dataset and the simulated ones
    initial_parameters = get_initial_parameters(input_path, method=method)

    # create environment for simulations generation and creation
    if not os.path.exists(simulations_dir + "data/"):
        res = os.system("mkdir -p " + simulations_dir + "data/")
    if not os.path.exists(simulations_dir + "jobs/"):
        res = os.system("mkdir -p " + simulations_dir + "jobs/")
    if not os.path.exists(simulations_dir + "output/"):
        res = os.system("mkdir -p " + simulations_dir + "output/")

    # create simulations
    if method == "traitrelax" and bootstrap_on_char_only:  # simulate only character data, and create parameter files that use the data_path
        cmd = "python /groups/itay_mayrose/halabikeren/myScripts/python/simulator/TraitSimulator.py -o " + simulations_dir + "data/ -t " + tree_path + " -mu " + \
              inferred_parameters["mu"] + " -pi0 " + inferred_parameters["pi0"] + " -imu " + initial_parameters[
                  "mu"] + " -ipi0 " + initial_parameters["pi0"] + " -ikappa " + initial_parameters[
                  "kappa"] + " -iomega0 " + initial_parameters["omega0"] + " -iomega1 " + initial_parameters[
                  "omega1"] + " -iomega2 " + initial_parameters["omega2"] + " -ip0 " + initial_parameters[
                  "p0"] + " -ip1 " + initial_parameters["p1"] + " -ik " + initial_parameters["k"] + " -rep " + str(
            num_of_replicates) + " -s " + data_path + " -rep " + str(num_of_replicates) + " -sc " + inferred_parameters[
                  "sc"] + " -in1t " + initial_parameters["1_Full.theta"] + " -n1t1 " + inferred_parameters[
                  "1_Full.theta1"] + " -in1t1 " + initial_parameters["1_Full.theta1"] + " -n1t2 " + inferred_parameters[
                  "1_Full.theta2"] + " -in1t2 " + initial_parameters["1_Full.theta2"] + " -n2t " + inferred_parameters[
                  "2_Full.theta"] + " -in2t " + initial_parameters["2_Full.theta"] + " -n2t1 " + inferred_parameters[
                  "2_Full.theta1"] + " -in2t1 " + initial_parameters["2_Full.theta1"] + " -n2t2 " + inferred_parameters[
                  "2_Full.theta2"] + " -in2t2 " + initial_parameters["2_Full.theta2"] + " -n3t " + inferred_parameters[
                  "3_Full.theta"] + " -in3t " + initial_parameters["3_Full.theta"] + " -n3t1 " + inferred_parameters[
                  "3_Full.theta1"] + " -in3t1 " + initial_parameters["3_Full.theta1"] + " -n3t2 " + inferred_parameters[
                  "3_Full.theta2"] + " -in3t2 " + initial_parameters["3_Full.theta2"] + " -imu " + initial_parameters[
                  "mu"] + " -ipi0 " + initial_parameters["pi0"]
    if method == "traitrelax" and without_trait_simulations:
        cmd = "python /groups/itay_mayrose/halabikeren/myScripts/python/simulator/RELAXSimulator.py -o " + simulations_dir + "data/ -t " + history_path + " -kappa " + \
              inferred_parameters["kappa"] + " -omega0 " + inferred_parameters["omega0"] + " -omega1 " + \
              inferred_parameters["omega1"] + " -omega2 " + inferred_parameters["omega2"] + " -p0 " + \
              inferred_parameters["p0"] + " -p1 " + inferred_parameters["p1"] + " -k 1 -ikappa " + initial_parameters[
                  "kappa"] + " -iomega0 " + initial_parameters["omega0"] + " -iomega1 " + initial_parameters[
                  "omega1"] + " -iomega2 " + initial_parameters["omega2"] + " -ip0 " + initial_parameters[
                  "p0"] + " -ip1 " + initial_parameters["p1"] + " -ik " + initial_parameters["k"] + " -al " + str(
            alignment_length) + " -rep " + str(num_of_replicates) + " -n1t " + inferred_parameters[
                  "1_Full.theta"] + " -in1t " + initial_parameters["1_Full.theta"] + " -n1t1 " + inferred_parameters[
                  "1_Full.theta1"] + " -in1t1 " + initial_parameters["1_Full.theta1"] + " -n1t2 " + inferred_parameters[
                  "1_Full.theta2"] + " -in1t2 " + initial_parameters["1_Full.theta2"] + " -n2t " + inferred_parameters[
                  "2_Full.theta"] + " -in2t " + initial_parameters["2_Full.theta"] + " -n2t1 " + inferred_parameters[
                  "2_Full.theta1"] + " -in2t1 " + initial_parameters["2_Full.theta1"] + " -n2t2 " + inferred_parameters[
                  "2_Full.theta2"] + " -in2t2 " + initial_parameters["2_Full.theta2"] + " -n3t " + inferred_parameters[
                  "3_Full.theta"] + " -in3t " + initial_parameters["3_Full.theta"] + " -n3t1 " + inferred_parameters[
                  "3_Full.theta1"] + " -in3t1 " + initial_parameters["3_Full.theta1"] + " -n3t2 " + inferred_parameters[
                  "3_Full.theta2"] + " -in3t2 " + initial_parameters["3_Full.theta2"] + " -imu " + initial_parameters[
                  "mu"] + " -ipi0 " + initial_parameters["pi0"] + " -cd " + character_data_path + " -sc " + \
              inferred_parameters["sc"]
    elif method == "traitrelax" and not without_trait_simulations:
        cmd = "python /groups/itay_mayrose/halabikeren/myScripts/python/simulator/TraitRELAXSimulator.py -o " + simulations_dir + "data/ -t " + tree_path + " -mu " + \
              inferred_parameters["mu"] + " -pi0 " + inferred_parameters["pi0"] + " -kappa " + inferred_parameters[
                  "kappa"] + " -omega0 " + inferred_parameters["omega0"] + " -omega1 " + inferred_parameters[
                  "omega1"] + " -omega2 " + inferred_parameters["omega2"] + " -p0 " + inferred_parameters[
                  "p0"] + " -p1 " + inferred_parameters["p1"] + " -k 1 -imu " + initial_parameters["mu"] + " -ipi0 " + \
              initial_parameters["pi0"] + " -ikappa " + initial_parameters["kappa"] + " -iomega0 " + initial_parameters[
                  "omega0"] + " -iomega1 " + initial_parameters["omega1"] + " -iomega2 " + initial_parameters[
                  "omega2"] + " -ip0 " + initial_parameters["p0"] + " -ip1 " + initial_parameters["p1"] + " -ik " + \
              initial_parameters["k"] + " -al " + str(alignment_length) + " -rep " + str(num_of_replicates) + " -n1t " + \
              inferred_parameters["1_Full.theta"] + " -in1t " + initial_parameters["1_Full.theta"] + " -n1t1 " + \
              inferred_parameters["1_Full.theta1"] + " -in1t1 " + initial_parameters["1_Full.theta1"] + " -n1t2 " + \
              inferred_parameters["1_Full.theta2"] + " -in1t2 " + initial_parameters["1_Full.theta2"] + " -n2t " + \
              inferred_parameters["2_Full.theta"] + " -in2t " + initial_parameters["2_Full.theta"] + " -n2t1 " + \
              inferred_parameters["2_Full.theta1"] + " -in2t1 " + initial_parameters["2_Full.theta1"] + " -n2t2 " + \
              inferred_parameters["2_Full.theta2"] + " -in2t2 " + initial_parameters["2_Full.theta2"] + " -n3t " + \
              inferred_parameters["3_Full.theta"] + " -in3t " + initial_parameters["3_Full.theta"] + " -n3t1 " + \
              inferred_parameters["3_Full.theta1"] + " -in3t1 " + initial_parameters["3_Full.theta1"] + " -n3t2 " + \
              inferred_parameters["3_Full.theta2"] + " -in3t2 " + initial_parameters["3_Full.theta2"] + " -sc " + \
              inferred_parameters["sc"] + " -in1t " + initial_parameters["1_Full.theta"] + " -n1t1 " + \
              inferred_parameters[
                  "1_Full.theta1"] + " -in1t1 " + initial_parameters["1_Full.theta1"] + " -n1t2 " + inferred_parameters[
                  "1_Full.theta2"] + " -in1t2 " + initial_parameters["1_Full.theta2"] + " -n2t " + inferred_parameters[
                  "2_Full.theta"] + " -in2t " + initial_parameters["2_Full.theta"] + " -n2t1 " + inferred_parameters[
                  "2_Full.theta1"] + " -in2t1 " + initial_parameters["2_Full.theta1"] + " -n2t2 " + inferred_parameters[
                  "2_Full.theta2"] + " -in2t2 " + initial_parameters["2_Full.theta2"] + " -n3t " + inferred_parameters[
                  "3_Full.theta"] + " -in3t " + initial_parameters["3_Full.theta"] + " -n3t1 " + inferred_parameters[
                  "3_Full.theta1"] + " -in3t1 " + initial_parameters["3_Full.theta1"] + " -n3t2 " + inferred_parameters[
                  "3_Full.theta2"] + " -in3t2 " + initial_parameters["3_Full.theta2"] + " -imu " + initial_parameters[
                  "mu"] + " -ipi0 " + initial_parameters["pi0"]
    else:  # use relax simulator
        cmd = "python /groups/itay_mayrose/halabikeren/myScripts/python/simulator/RELAXSimulator.py -o " + simulations_dir + "data/ -t " + history_path + " -kappa " + \
              inferred_parameters["kappa"] + " -omega0 " + inferred_parameters["omega0"] + " -omega1 " + \
              inferred_parameters["omega1"] + " -omega2 " + inferred_parameters["omega2"] + " -p0 " + \
              inferred_parameters["p0"] + " -p1 " + inferred_parameters["p1"] + " -k 1 -ikappa " + initial_parameters[
                  "kappa"] + " -iomega0 " + initial_parameters["omega0"] + " -iomega1 " + initial_parameters[
                  "omega1"] + " -iomega2 " + initial_parameters["omega2"] + " -ip0 " + initial_parameters[
                  "p0"] + " -ip1 " + initial_parameters["p1"] + " -ik " + initial_parameters["k"] + " -al " + str(
            alignment_length) + " -rep " + str(num_of_replicates) + " -n1t " + inferred_parameters[
                  "1_Full.theta"] + " -in1t " + initial_parameters["1_Full.theta"] + " -n1t1 " + inferred_parameters[
                  "1_Full.theta1"] + " -in1t1 " + initial_parameters["1_Full.theta1"] + " -n1t2 " + inferred_parameters[
                  "1_Full.theta2"] + " -in1t2 " + initial_parameters["1_Full.theta2"] + " -n2t " + inferred_parameters[
                  "2_Full.theta"] + " -in2t " + initial_parameters["2_Full.theta"] + " -n2t1 " + inferred_parameters[
                  "2_Full.theta1"] + " -in2t1 " + initial_parameters["2_Full.theta1"] + " -n2t2 " + inferred_parameters[
                  "2_Full.theta2"] + " -in2t2 " + initial_parameters["2_Full.theta2"] + " -n3t " + inferred_parameters[
                  "3_Full.theta"] + " -in3t " + initial_parameters["3_Full.theta"] + " -n3t1 " + inferred_parameters[
                  "3_Full.theta1"] + " -in3t1 " + initial_parameters["3_Full.theta1"] + " -n3t2 " + inferred_parameters[
                  "3_Full.theta2"] + " -in3t2 " + initial_parameters["3_Full.theta2"] + " -sc " + inferred_parameters[
                  "sc"]
    res = os.system(cmd)

    # add scaling to the parameter files, if needed - as of 19.2.20, scaling optimization in parametric bootstrapping is disabled until the convergence issue related to scalingo ptimization is resolved
    parameters_dir = simulations_dir + "data/"
    traitrelax_param_path = parameters_dir + "traitrelax_param/"
    relax_param_path = parameters_dir + "relax_param/"
    if set_advanced_parameters:
        for path in os.listdir(traitrelax_param_path):
            with open(traitrelax_param_path + path, "a") as infile:
                infile.write("\ncharacter.use_analytic_mapping=1")
                # infile.write("\noptimization.scale.tree=1")
                infile.write("\noptimization.advanced=1")
                # for path in os.listdir(relax_param_path):
                #     with open(relax_param_path + path, "a") as infile:
                #         infile.write("\noptimization.scale.tree=1")

    # submit the jobs to execute the program on the simulations
    parameters_dir = simulations_dir + "data/" + method + "_param/"
    cmd = "python /groups/itay_mayrose/halabikeren/myScripts/python/bpp/runBppProgram.py -pp " + program_path + " -pd " + parameters_dir + " -jd " + simulations_dir + "jobs/" + " -err " + simulations_dir + "output/" + " -q " + queue + " -rn " + str(
        num_of_replicates)
    if run_inference:
        res = os.system(cmd)
    else:
        with open(simulations_dir + "inference_ exec_cmd.txt", "w") as outfile:
            outfile.write(cmd)

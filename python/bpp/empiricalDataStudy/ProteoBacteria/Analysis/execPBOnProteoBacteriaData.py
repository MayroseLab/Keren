import argparse, os, re

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='Creates simulations based on inferred parameters and executes inference program on them. This pipeline is compatible both with RELAX and with TraitRELAX')
    parser.add_argument('--program_path', '-p', help='path to the program file', required=False,
                        default="/groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite/traitrelax")
    parser.add_argument('--input_dir', '-i',
                        help='path to the directory with the parameter files given as input to the program to create the inference result. Form here, we will extract the initial parameters for the parameter files of the simulations',
                        required=True)
    parser.add_argument('--inference_dir', '-r', help='directory of the output files of the program on the data',
                        required=True)
    parser.add_argument('--simulations_dir', '-s', help='directory to create simulations environment in per dataset',
                        required=True)
    parser.add_argument('--data_dir', '-d', help='path to the sequence alignment of the data', required=True)
    parser.add_argument('--tree_dir', '-t', help='directory to the tree paths of the data in data_path', required=True)
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
                        required=False, default=100)
    parser.add_argument('--queue', '-q', help='name of the queue to send inference jobs to', required=False,
                        default="itaym1")
    parser.add_argument('--priority', '-pr', help='priority of inference jobs', required=False, default="0")
    parser.add_argument('--without_trait_simulations', '-nt',
                        help='boolean indicating weather trait data should be simulated or not', required=False,
                        default=True)
    parser.add_argument('--character_data_path', '-cd',
                        help='path to the character data file to be given as input to traitrelax parameter file in case trait data is not simulated',
                        required=False, default="")

    args = parser.parse_args()
    program_path = args.program_path
    input_dir = args.input_dir
    inference_dir = args.inference_dir
    simulations_dir = args.simulations_dir
    data_dir = args.data_dir
    bootstrap_on_char_only = bool(int(args.bootstrap_on_char_only))
    set_advanced_parameters = bool(int(args.set_advanced_parameters))
    tree_dir = args.tree_dir
    history_path = args.history_path
    num_of_replicates = int(args.num_of_replicates)
    queue = args.queue
    priority = args.priority
    without_trait_simulations = bool(int(args.without_trait_simulations))
    character_data_path = args.character_data_path

    dataset_id_regex = re.compile("Parsing file .*?([^\/]*?).bpp for options", re.MULTILINE | re.DOTALL)
    for path in os.listdir(inference_dir):
        try:
            # get the dataset id from the file
            with open(inference_dir + path, "r") as infile:
                content = infile.read()
            dataset_id = dataset_id_regex.search(content).group(1)

            # get tree path
            tree_path = tree_dir + dataset_id + "/optimized_tree.nwk"

            # get input_path
            input_path = input_dir + dataset_id + ".bpp"

            # get data_path
            data_path = data_dir + dataset_id + ".fa"

            # set simulation dir
            dataset_simulation_dir = simulations_dir + dataset_id + "/"
            if not os.path.exists(dataset_simulation_dir):
                res = os.system("mkdir -p " + dataset_simulation_dir)
            else:
                continue

            # set the command
            cmd = 'python /groups/itay_mayrose/halabikeren/myScripts/python/bpp/empiricalDataStudy/parametricBootstrap.py -s ' + dataset_simulation_dir + ' -r ' + inference_dir + path + ' -d ' + data_path + ' -t ' + tree_path + ' -hist ' + history_path + ' -i ' + input_path + ' -cd ' + character_data_path + ' -sc ' + str(
                int(set_advanced_parameters)) + ' -rep ' + str(
                num_of_replicates) + ' -q ' + queue + ' -p ' + program_path + ' -nt ' + str(
                int(without_trait_simulations)) + ' -c ' + str(int(bootstrap_on_char_only)) + ' -ri 0'
            res = os.system(cmd)
        except Exception as e:
            print("FAILED TO EXECUTE PARAMETERIC BOOSTRAPPING ON PATH: ", path)

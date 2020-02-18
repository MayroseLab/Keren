import argparse, os

if __name__ == '__main__':

    # process input from cmd
    parser = argparse.ArgumentParser(description='creates parameter files based on a parameter template files')
    parser.add_argument('--param_files_dir', '-pd', help='directory to create the parameters files in',
                        required=True)
    parser.add_argument('--param_template_path', '-pt', help='path to the parameter template file',
                        required=False, default=os.getcwd() + "/jobs/")
    parser.add_argument('--sequence_files_dir', '-sd', help='path of the input sequence files',
                        required=False, default=os.getcwd() + "/jobs_output/")
    parser.add_argument('--expected_histories_dir', '-hd', help='directory to hold the final evaluated expected history per execution. Relevant only in case of TraitRELAX', required=False, default="")
    parser.add_argument('--optimization_log_dir', '-ld', help='directory to hold the log file of each execution', required=True, default=os.getcwd() + "/log/")

    args = parser.parse_args()
    param_output_dir = args.param_files_dir
    if not os.path.exists(param_output_dir):
        res = os.system(param_output_dir)
    sequence_data_dir = args.sequence_files_dir
    param_template_path = args.param_template_path
    expected_histories_dir = args.expected_histories_dir
    if len(expected_histories_dir) > 0 and not os.path.exists(expected_histories_dir):
        res = os.system("mkdir -p " + expected_histories_dir)
    optimization_log_dir = args.optimization_log_dir
    if not os.path.exists(optimization_log_dir):
        res = os.system("mkdir -p " + optimization_log_dir)

    with open(param_template_path, "r") as param_template_file:
            param_template_content = param_template_file.read()

    for sequence_data_path in os.listdir(sequence_data_dir):
        full_seq_data_path = sequence_data_dir + sequence_data_path
        optimization_log_path = optimization_log_dir + sequence_data_path.replace("sequence_data", "optimzation_output")
        optimization_log_path = optimization_log_path.replace(".fas", ".log")
        expected_hitory_path = expected_histories_dir + sequence_data_path.replace("sequence_data", "expected_history")
        expected_hitory_path = expected_hitory_path.replace(".fas", ".nwk")
        param_content = param_template_content.replace("SEQ_DATA_PATH", full_seq_data_path)
        param_content = param_content.replace("LOG_PATH", optimization_log_path)
        param_content = param_content.replace("HISTORY_PATH", expected_hitory_path)
        param_name = sequence_data_path + ".bpp"
        param_output_path = param_output_dir + "param_" + param_name
        with open(param_output_path, "w") as param_output_file:
            param_output_file.write(param_content)

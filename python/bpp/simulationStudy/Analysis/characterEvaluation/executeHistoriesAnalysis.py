import os, argparse, re, sys
sys.path.append("/groups/itay_mayrose/halabikeren/myScripts/python/")
from utils.createJobFile import set_job_env, create_job_file


if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
    description='simulates alignments and character history under the TraitRELAX null and alternative models using Bio++ and INDELible')
    parser.add_argument('--input_dir', '-i', help='directory that holds the input parameters files', required=True)
    parser.add_argument('--jobs_dir', '-jd', help='directory that will hold the .sh files of the executed jobs', required=True)
    parser.add_argument('--error_dir', '-err', help='directory that will hold the stdout and stderr files of the jobs', required=True)
    parser.add_argument('--mu_options', '-mu', help='list of values of mu to include in the analysis', required=False, default=[1, 4, 8])
    parser.add_argument('--taxa_num_options', '-tn', help='list of taxa number values to include in the analysis', required=False, default=[32])
    parser.add_argument('--positions_num_options', '-pn', help='list of positions number values to include in the analysis', required=False, default=[300])
    parser.add_argument('--k_options', '-ko', help='list of k values to include in the analysis', required=False, default=[0.5])
    parser.add_argument('--replicates_number', '-rn', help='number of replicates to execute jobs on per combo', required=False, default=50)
    parser.add_argument('--queue', '-q', help='name of queue to send jobs to', required=False, default="itaym")

    args = parser.parse_args()
    input_dir = args.input_dir
    jobs_dir = args.jobs_dir
    if not os.path.exists(jobs_dir):
        res = os.system("mkdir -p " + jobs_dir)
    error_dir = args.error_dir
    if not os.path.exists(error_dir):
        res = os.system("mkdir -p " + error_dir)

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

    replicates_number = int(args.replicates_number)
    queue = args.queue
    evaluation_program_path = "/groups/itay_mayrose/halabikeren/biopp/bppsuite/build/bppSuite/debugexpectedhistory"
    id_regex = re.compile("(\d*)\.bpp")
    # for each combo, execute <replicates_number> calls to the evaluation program, for the first replicates
    for mu in mu_options:
        for taxa in taxa_options:
            for pos in positions_options:
                for k in k_options:
                    print("executing jobs for combo: (mu=", mu, ", #taxa=", taxa, ", #pos=", pos, ", k=", k, ")")
                    parameter_files_dir = input_dir + "tbl_4_mu_" + str(int(mu)) + "_pi0_0.5_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa) + "_taxa/" + str(pos) + "_codons/k_" + str(k) + "/traitrelax_param/"
                    print("corresponding directory: ", parameter_files_dir)
                    job_files_dir = jobs_dir + "tbl_4_mu_" + str(mu) + "_pi0_0.5_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa) + "_taxa/" + str(pos) + "_codons/k_" + str(k)
                    if not os.path.exists(job_files_dir):
                        res = os.system("mkdir -p " + job_files_dir)
                    error_files_dir = error_dir + "tbl_4_mu_" + str(mu) + "_pi0_0.5_kappa_2_p_0.125_omega1_0.8_omega2_2_theta1_0.5_theta2_0.8/" + str(taxa) + "_taxa/" + str(pos) + "_codons/k_" + str(k)
                    if not os.path.exists(error_files_dir):
                        res = os.system("mkdir -p " + error_files_dir)
                    if os.path.exists(parameter_files_dir):
                        parameter_files = os.listdir(parameter_files_dir)
                        for i in range(min(len(parameter_files), replicates_number)):
                            parameters_filepath = parameter_files[i]

                            id = int(id_regex.search(parameters_filepath).group(1))
                            job_name = "ExpectedHistoryEval_" + parameters_filepath.replace(".bpp", "").replace("param_simulated_data.", "")
                            file_name = job_name + ".sh"
                            cmd = evaluation_program_path + " param=" + parameter_files_dir + parameters_filepath
                            commands = [cmd]  # ["export OMP_NUM_THREADS=4", cmd]
                            touch_file_path = job_name + "_flag_done"
                            full_job = create_job_file(job_name, commands, file_name, error_files_dir, job_files_dir,
                                                       0, 1,
                                                       touch_file_path, limit_nodes=False, python=False, language="bash", queue=queue)
                            res = os.system("qsub " + full_job)




import sys, re, argparse, os
sys.path.append("/groups/itay_mayrose/halabikeren/myScripts/python/")
from utils.createJobFile import set_job_env, create_job_file

if __name__ == '__main__':

    # process input from cmd
    parser = argparse.ArgumentParser(description='executes RELAX in Bio++ (random effects of restrictive approach) on parameter files in a designated directory')
    parser.add_argument('--param_files_dir', '-pd', help='directory that holds the parameters files',
                        required=True)
    parser.add_argument('--job_files_path', '-jd', help='path to hold the RELAX inference jobs and their output',
                        required=False, default=os.getcwd() + "/jobs/")
    parser.add_argument('--error_files_path', '-err', help='path to hold the job error and output files',
                        required=False, default=os.getcwd() + "/jobs_output/")
    parser.add_argument('--priority', '-pr', help='priority for the submitted jobs', required=False, default=0)
    parser.add_argument('--program_path', '-pp', help='path to the program executable', required=False, default=0)
    parser.add_argument('--queue', '-q', help='the name of the queue to submit jobs to', required=False, default="itaymr")
    parser.add_argument('--replicates_num', '-rn', help='number of replicates to run jobs on', required=False, default=50)
    parser.add_argument('--replicates_range', '-rr', help='min id of replicate and max id of replicates. If provided, the script will submit jobs for the ids in between only', required=False, default=None)
    parser.add_argument('--memory_alloc', '-mem', help='memory allocation per execution', required=False, default="4gb")

    args = parser.parse_args()
    param_files_dir = args.param_files_dir
    job_files_path = args.job_files_path
    if not os.path.exists(job_files_path):
        res = os.system("mkdir -p " + job_files_path)
    error_files_path = args.error_files_path
    if not os.path.exists(error_files_path):
        res = os.system("mkdir -p " + error_files_path)
    priority = int(args.priority)
    queue = args.queue
    program_path = args.program_path
    replicates_num = int(args.replicates_num)
    replicates_range = args.replicates_range
    if not replicates_range == None and not type(replicates_range) == list:
        str_lst = replicates_range.split(",")
        replicates_range = [int(i) for i in str_lst]
    memory_alloc = args.memory_alloc

    print("\n\nreplicates num: ", replicates_num, "\n\n")
    set_job_env(job_files_path, error_files_path)
    jobs_counter = 0
    id_regex = re.compile("(.*?)\.bpp", re.DOTALL)

    if not replicates_range == None:
        for i in range(replicates_range[0], replicates_range[1]):
            param_file = str(i) + ".bpp"
            job_name = "Bpp_" + param_file.replace(".bpp", "").replace("param_simulated_data.", "")
            file_name = job_name + ".sh"
            cmd = program_path + " param=" + param_files_dir + param_file
            commands = [cmd]  # ["export OMP_NUM_THREADS=4", cmd]
            touch_file_path = job_name + "_flag_done"
            full_job = create_job_file(job_name, commands, file_name, error_files_path, job_files_path, priority, 1,
                                       touch_file_path,
                                       allowed_nodes=["compute-0-246", "compute-0-247", "compute-0-248",
                                                      "compute-0-249"], limit_nodes=False, python=False,
                                       openmpi=False, language="bash", queue=queue)
            res = os.system("qsub -p " + str(priority) + " " + full_job)

    else:
        for param_file in os.listdir(param_files_dir):
            if ".bpp" in param_file and jobs_counter < replicates_num:
                id = id_regex.search(param_file).group(1)
                job_name = "Bpp_" + param_file.replace(".bpp", "").replace("param_simulated_data.", "")
                file_name = job_name + ".sh"
                cmd = program_path + " param=" + param_files_dir+param_file
                commands = [cmd] #["export OMP_NUM_THREADS=4", cmd]
                touch_file_path =  job_name + "_flag_done"
                full_job = create_job_file(job_name, commands, file_name, error_files_path, job_files_path, priority, 1,
                                                   touch_file_path, allowed_nodes=["compute-0-246", "compute-0-247", "compute-0-248", "compute-0-249"], limit_nodes=False, python=False,
                                                   openmpi=False, language="bash", queue=queue, mem_alloc=memory_alloc)
                res=os.system("qsub -p " + str(priority) + " " + full_job)
                jobs_counter += 1


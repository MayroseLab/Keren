import sys, re, argparse, os
from ete3 import Tree
sys.path.append("/groups/itay_mayrose/halabikeren/myScripts/python/")
from utils.createJobFile import set_job_env, create_job_file

def get_first_state(tree_path):
    with open (tree_path, "r") as infile:
        tree_str = infile.read()
    state_0_appearance = tree_str.find("{0}")
    state_1_appearance = tree_str.find("{1}")
    if state_0_appearance < state_1_appearance:
        return 0
    return 1

if __name__ == '__main__':

    # process input from cmd
    parser = argparse.ArgumentParser(description='executes RELAX in Bio++ (random effects of restrictive approach) on parameter files in a designated directory using either a maximum parsimony based solution or a simulated history')
    parser.add_argument('--input_data_path', '-i', help='directory that holds the input files', required=True)
    parser.add_argument('--trees_dir', '-t', help='directory holding the trees of the data in case that MP histories need to be constructed', required=True)
    parser.add_argument('--job_files_path', '-jd', help='path to hold the RELAX inference jobs and their output',
                        required=False, default=os.getcwd() + "/jobs/")
    parser.add_argument('--error_files_path', '-err', help='path to hold the job error and output files',
                        required=False, default=os.getcwd() + "/jobs_output/")
    parser.add_argument('--priority', '-pr', help='priority for the submitted jobs', required=False, default=0)
    parser.add_argument('--queue', '-q', help='the name of the queue to submit jobs to', required=False, default="itaymr")
    parser.add_argument('--replicates_num', '-rn', help='number of replicates to run jobs on', required=False, default=50)
    parser.add_argument('--use_true_histories', '-ut', help='boolean indicsting weather true (i.e., simulated) trait histories should be used, or maximum parsimony based solutions. Defalut in maximum parsimony,', required=False, default=0)

    args = parser.parse_args()
    input_data_path = args.input_data_path
    trees_dir = args.trees_dir
    job_files_path = args.job_files_path
    error_files_path = args.error_files_path
    priority = int(args.priority)
    queue = args.queue
    replicates_num = int(args.replicates_num)
    use_true_histories = bool(int(args.use_true_histories))


    set_job_env(job_files_path, error_files_path)

    jobs_counter = 0
    id_regex = re.compile("(.*?)\.bpp", re.DOTALL)

    if not use_true_histories:
        # verify that MP histories exist and if not - create them
        res=os.system(" python /groups/itay_mayrose/halabikeren/myScripts/python/bpp/simulationStudy/DataPreparation/CreateMPHistories.py -t " + trees_dir + ' -c ' + input_data_path + ' -o ' + input_data_path + "mp_param/ -p " + input_data_path + "relax_param/ -u 1")

    for replicate in range(replicates_num):
        replicate_data_path = input_data_path + "replicate_" + str(replicate) + "/"
        sequence_data_path = replicate_data_path + "sequence_data/sequence_data_1.fas"
        labeled_tree_path = replicate_data_path + "mp_data/mp_history.nwk"
        if use_true_histories:
            labeled_tree_path = replicate_data_path + "character_data/true_history.nwk"
        # test_category_index = 1 # compatible with hyphy 2.3.720171228beta(MP) but in 2.5.8 the first state is always 0
        # first_category = get_first_state(labeled_tree_path)
        # if first_category == 0:
        #     test_category_index = 2
        test_category_index = 2
        job_name = "HyPhy_" + str(replicate)
        file_name = job_name + ".sh"
        cmds = ['source ~/.bashrc' , 'conda activate hyphy' , '(printf "1\\\\n7\\\\n1\\\\n' + sequence_data_path + '\\\\n' + labeled_tree_path + '\\\\n' + str(test_category_index) + '\\\\n2\\\\n" && cat) | hyphy'] # switched with HYPHYMP
        touch_file_path = job_name + "_flag_done"
        full_job = create_job_file(job_name, cmds, file_name, error_files_path, job_files_path, priority, 1,
                                                       touch_file_path, limit_nodes=False, python=False,
                                                       openmpi=True, language="bash", queue=queue)
        res=os.system("qsub -q " + queue + " -p " + str(priority) + " " + full_job)




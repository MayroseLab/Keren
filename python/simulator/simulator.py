import os, sys, argparse
from time import sleep
from create_job_file import create_job_file

output_names = {"leafs_msa": "output_TRUE_1.phy", "all_nodes_msa": "output_ANCESTRAL_1.phy", "fasta": "output_1.fas", "rates_and_omegas": "output_RATES.txt"
				,"randomized_tree": "trees.txt", "control": "control"}

orig_to_new_names = {"output_TRUE_1.phy": "leafs_msa.phy", "output_ANCESTRAL_1.phy": "all_nodes_msa.phy", "output_1.fas": "fasta.fas", "output_RATES.txt": "rates_and_omegas.txt"}

job_files_path = "/groups/itay_mayrose/halabikeren/jobs/simulations/"
error_files_path = "/groups/itay_mayrose/halabikeren/error_files/"


def set_sim_env(dir):
	folder_names = output_names.keys()
	for folder_name in folder_names:
		folder_dir = dir + folder_name
		if not os.path.exists(folder_dir):
			res = os.system("mkdir " + folder_dir)


def const_job_per_tree(output_dir, trees_dir, priority):

    # submit a job foreach required dataset
    print("**** submitting jobs for datasets simulation ****")
    jobs_counter = 0
    for tree in os.listdir(trees_dir):
        tree_path = trees_dir + tree
        if jobs_counter == 42:
            sleep(42)
            jobs_counter = 0
        dataset_id = tree.replace(".tree", "")
        job_name = "generate_alignments_" + dataset_id
        file_name = job_name + ".sh"
        command = "python /groups/itay_mayrose/halabikeren/python_scripts/simulator/simulate_constant.py " + dataset_id + " " + output_dir + " " + tree_path + " " + priority
        create_job_file(job_name, command, file_name, error_files_path, job_files_path)
        cmd = "qsub " + job_files_path + file_name
        os.system(cmd)
        jobs_counter += 1


if __name__ == "__main__":

    # process input from cmd
    parser = argparse.ArgumentParser(description='simulate alignment prone to positive selection using INDELible')
    parser.add_argument("--output_dir", "-o", help="direcotry which will hold the simulation's output", required=True)
    parser.add_argument("--simulation_type", "-s", help="the type of the simulation could either be random - in which random trees are generated, or constant - in which the trees for the simulation are provided in advance", required=True)
    parser.add_argument("--trees_dir", "-t", help="directory to the trees based on which the simulations will generate alignments - in case of random trees, the trees will be generated in it, and in case of constant - the trees will be provided in it", required=True)
    parser.add_argument("--error_files_path", "-e", help="directory in which the output logs of each job will be generated. default value: /groups/itay_mayrose/halabikeren/error_files/", required=False, default="/groups/itay_mayrose/halabikeren/error_files/")
    parser.add_argument("--job_files_path", "-j", help="directory in which each job file will be generated. default value: /groups/itay_mayrose/halabikeren/jobs/simulations/", required=False, default="/groups/itay_mayrose/halabikeren/jobs/simulations/")
    parser.add_argument("--priority", "-p", help="priority for the submitted jobs", required=True)
    args = parser.parse_args()
    output_dir = args.output_dir
    simulation_type = args.simulation_type
    trees_dir = args.trees_dir
    if not os.path.exists(trees_dir):
        res = os.system("mkdir " + trees_dir)
    error_files_path = args.error_files_path
    if not os.path.exists(error_files_path):
        res = os.system("mkdir " + error_files_path)
    job_files_path = args.job_files_path
    if not os.path.exists(job_files_path):
        res = os.system("mkdir " + job_files_path)
    priority = args.priority

    # set the simulation environment
    set_sim_env(output_dir)

    if simulation_type == "constant":
        const_job_per_tree(output_dir, trees_dir, priority)
    # elif simulation_type == "random":
    #     continue


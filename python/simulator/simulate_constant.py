import sys, os
from create_control_file import create_control_file
from create_job_file import create_indelible_job_file
from time import sleep

output_names = {"leafs_msa": "output_TRUE_1.phy", "all_nodes_msa": "output_ANCESTRAL_1.phy", "fasta": "output_1.fas", "rates_and_omegas": "output_RATES.txt"
    ,"randomized_tree": "trees.txt", "control": "control"}

orig_to_new_names = {"output_TRUE_1.phy": "leafs_msa.phy", "output_ANCESTRAL_1.phy": "all_nodes_msa.phy", "output_1.fas": "fasta.fas", "output_RATES.txt": "rates_and_omegas.txt", "trees.txt": "randomized_tree_info.txt"}

def set_sim_env(dir):
    folder_names = output_names.keys()
    for folder_name in folder_names:
        folder_dir = dir + folder_name
        if not os.path.exists(folder_dir):
            res = os.system("mkdir " + folder_dir)



def generate_alignments(dataset_id, output_dir, tree_path, priority):

    # set the simulation environment
    set_sim_env(output_dir)

    control_dir = output_dir + "control/" + dataset_id + "/"
    res = os.system("mkdir " + control_dir)
    dataset_output_dir = output_dir + dataset_id + "/"
    res = os.system("mkdir " + dataset_output_dir)
    control_txt = control_dir + "control.txt"
    create_control_file(control_dir, numberOfTaxa, 1, submodel, model_file, dest, tree, "const", indels_rate)  # create a control file that generates an alignment according to a randomized tree
    job_files_path = "/groups/itay_mayrose/halabikeren/jobs/simulations/indelible/"
    error_files_path = "/groups/itay_mayrose/halabikeren/error_files/"
    job_name = "run_indelible_" + dataset_id
    file_name = job_name + ".sh"
    # command style in power8:
    # cd OUTPUT_PATH
    # (echo CONTROL_FILE_PATH & & cat) | /groups/itay_mayrose/halabikeren/indelible/INDELibleV1.03/src/indelible
    create_indelible_job_file(job_name, dataset_output_dir, control_txt, file_name, error_files_path, job_files_path)
    full_job = job_files_path + file_name
    res = os.system("qsub -p " + priority + " " + full_job)
    while not os.path.exists(dataset_output_dir + "flag_indelible_done"):
        sleep(60)

    # while job not done - use dot file now!
    for folder_name in output_names.keys():
        if folder_name != "control":
            orig_file_name = dataset_output_dir + output_names[folder_name]
            new_file_name = output_dir + folder_name + "/" + dataset_id + "_" + orig_to_new_names[output_names[folder_name]]
            os.rename(orig_file_name, new_file_name)

    res = os.system("rm -rf " + dataset_output_dir)
    print('**** alignment generation is complete ****')



if __name__ == "__main__":
    dataset_id = sys.argv[1].rstrip()
    output_dir = sys.argv[2].rstrip()
    tree_path = sys.argv[3].rstrip()
    priority = sys.argv[4].rstrip()

    generate_alignments(dataset_id, output_dir, tree_path, priority)




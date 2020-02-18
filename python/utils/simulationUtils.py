import os, sys
from time import sleep
from ete3 import Tree
sys.path.append("/groups/itay_mayrose/halabikeren/myScripts/python/utils/")
from createJobFile import create_job_file, set_job_env

output_names = {"output_TRUE_1.fas": "leafs_msa", "output_ANCESTRAL_1.fas": "all_nodes_msa", "output_1.fas": "fasta", "output_RATES.txt": "rates_and_omegas"
    ,"trees.txt": "random_tree", "control": "control", "jobs": "jobs", "errors": "errors"}
orig_to_new_names = {"output_TRUE_1.fas": "leafs_msa.fas", "output_ANCESTRAL_1.fas": "all_nodes_msa.fas", "output_1.fas": "fasta.fas", "output_RATES.txt": "rates_and_omegas.txt", "trees.txt": "randomized_tree_info.txt"}
simulation_type_translator = {"AA": "AMINOACID 2", "CODON": "CODON 1", "NUC": "NUCLEOTIDE 1", "sparta": "AMINOACID 2"}
models_translator = {"AA": "JTT", "CODON": "2", "NUC": "HKY", "sparta": "WAG"}


# auxiliary function to set simulation environment
def set_sim_env(dir):
    if not os.path.exists(dir):
        res = os.system("mkdir " + dir)
    folder_names = output_names.values()
    for folder_name in folder_names:
        folder_dir = dir + folder_name
        if not os.path.exists(folder_dir):
            res = os.system("mkdir " + folder_dir)
    return 0

#### create INDELible control file ####

def create_control_file(control_dir, source, tree=None, ultrametric=False, simulation_type="AA", indels_rate=0.01, a_shape=1.7, max_indel_length=50, root_seq_len=100, pinv=0, alpha=0.5, ngamcat=16, numberOfTaxa=None, numberOfTrees=1, random_tree_size=10, positive_selection_flag=False, ps_model="alternative"):
    '''generates control.txt file in indelDir according to the given parameters'''
    control = control_dir + "control.txt"
    with open(control, 'w+') as handle:
        handle.write('[TYPE] ' + simulation_type_translator[simulation_type] + '\n')
        handle.write('\n')
        handle.write('[SETTINGS]\n')
        handle.write('\t[ancestralprint]	NEW\n')
        handle.write('\t[phylipextension] phy\n')
        handle.write('\t[nexusextension] nex\n')
        handle.write('\t[fastaextension] fas\n')
        handle.write('\t[output]          FASTA\n')
        handle.write('\t[fileperrep]      TRUE\n')
        handle.write('\t[printrates] TRUE     // FALSE or TRUE\n')
        handle.write('\n')
        handle.write('[MODEL] tree_model\n')
        if positive_selection_flag:
            handle.write('\t[submodel]  2\n') # set kappa to be 2
            if ps_model == "alternative":
                last_omegas = [1.117919318, 1.484136503, 2.534366734]
            else:
                last_omegas = [1,1,1]
            handle.write('\t         0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02\n')  # set 50 bins of equal area (=probability) of the omegas distribution
            handle.write('\t         0.009485947	0.016198557	0.021505013	0.026346363	0.030985163	0.035541286	0.040083290	0.044656451	0.049294022	0.054022565	0.058864804	0.063841293	0.068971480	0.074274440	0.079769421	0.085476289	0.091415913	0.097610532	0.104084119	0.110862769	0.117975121	0.125452840	0.133331166	0.141649568	0.150452515	0.159790413	0.169720740	0.180309453	0.191632725	0.203779147	0.216852509	0.230975385	0.246293807	0.262983421	0.281257761	0.301379518	0.323676215	0.348562470	0.376572407	0.408408213	0.445015329	0.487703553	0.538351547	0.599772865	0.676420924	0.775883263	0.912492899\t' + str(last_omegas[0]) + '\t' + str(last_omegas[1]) + '\t' + str(last_omegas[2]) + '\n') # fraction of positive selection classes: 0.06
        else:
            handle.write('\t[submodel] ' + models_translator[simulation_type] + '\n')
            handle.write('\t[rates] ' + str(pinv) + ' ' + str(alpha) + ' ' + str(ngamcat) + '      //   16 category discrete gamma with alpha=2\n')
        handle.write('\t[indelmodel] POW ' + str(a_shape) + ' ' + str(max_indel_length) + '\n')  #need to ask about the a=1.5 parameter
        handle.write('\t[indelrate] ' + str(indels_rate) + ' \n')
        handle.write('\n')
        if (source == "random"):
            handle.write('[TREE] tree\n')
            handle.write('\t[rooted] ' + str(numberOfTaxa) + ' 6.7 2.5 0.234 0.31\n')
            handle.write('\t[seed] 1 \n')
            if ultrametric:
                handle.write('\t[branchlengths] ULTRAMETRIC\n')
            handle.write('\t[treelength] ' + str(random_tree_size) + '\n')
        elif (source == "constant"):
            tree_str = adjust_tree(tree)
            if "e-" in tree_str:
                print("error in adjusting tree format!")
                print("tree: ", tree)
                print("source: ", source)
                exit(1)
            handle.write('[TREE] tree ' + tree_str + "\n")
            # if positive_selection_flag:
            #     handle.write('\t[branchlengths] NON-ULTRAMETRIC\n')
        handle.write('\n')
        handle.write('\n')
        handle.write('[PARTITIONS] job\n')
        handle.write('\t[tree tree_model ' + str(root_seq_len) + ']\n')
        handle.write('\n')
        handle.write('[EVOLVE] job ' + str(numberOfTrees) + ' output\n')
    return control


def reformat_num(number):
    # print("number's type: ", type(number))
    num = float(number)
    reformatted_num = str("{0:.10f}".format(num))
    return reformatted_num


def reformat_branch(branch):
    fixed_branch = branch[0] + ":" + branch[1:]
    return fixed_branch


def adjust_tree(tree):
    # Haim's version - used due to a bug in my script
    orig_tree_path = tree
    reformatted_tree_path = orig_tree_path + ".reformatted"
    try:
        res = os.system("perl /groups/pupko/haim/Scripts/validate_tree.pl " + orig_tree_path + " " + reformatted_tree_path)
        if os.path.exists(reformatted_tree_path):
            print("tree format was faulty")
            with open(reformatted_tree_path, "r") as tree_file:
                tree_str = tree_file.read()
            res = os.system("rm -r " + orig_tree_path)
            res = os.rename(reformatted_tree_path, orig_tree_path)
        else:
            with open(orig_tree_path, "r") as tree_file:
                tree_str = tree_file.read()
    except:
        with open(orig_tree_path, "r") as tree_file:
            tree_str = tree_file.read()
    return tree_str
    # # my script - missing some unknown fixing
    # with open(tree, "r") as tree_file:
    #     tree_str = tree_file.read()
    #     # if there is ROOT string in the end of the tree -> remove it
    #     tree_str = tree_str.replace('ROOT', '')
    #     # remove internal node names
    #     internal_nodes_pattern = 'N[0-9]+'
    #     pattern = re.compile(internal_nodes_pattern)
    #     tree_str = pattern.sub('', tree_str)
    #     # reformat small numbers
    #     reformatted_number = re.compile(r'\d*[\.\d+]*e\-\d+', re.MULTILINE|re.DOTALL)
    #     numbers = reformatted_number.findall(tree_str)
    #     for number in numbers:
    #         reformatted_num = reformat_num(number)
    #         tree_str = tree_str.replace(number, reformatted_num, 1)
    #     # reformat the branches
    #     comb_branch_pattern = re.compile(r'\)[0-9]+\.[0-9]+', re.MULTILINE|re.DOTALL)
    #     for branch in re.findall(comb_branch_pattern, tree_str):
    #         reformatted_branch = reformat_branch(branch)
    #         tree_str = re.sub(re.escape(branch), reformatted_branch, tree_str)
    #     #remove bootsrap numbers
    #     bootsrap_num_expr = re.compile(r'\)-{0,1}[0-9]+\.*[0-9]*:', re.MULTILINE|re.DOTALL)
    #     tree_str = bootsrap_num_expr.sub(')', tree_str)
    # return tree_str

# extract tree from simulation output
# auxiliary function to extract the rooted tree from the debug data simulation file
def extract_rooted_tree(data_path, output_path, with_internals=False):
    with open(data_path, "r") as data_file:
        for i in range(8):
            data_file.readline() # get rid of the first 8 lines
        line = data_file.readline()
        line_content = line.split("\t")
        tree_str = line_content[8]
    tree = Tree(tree_str, format=1)
    if with_internals: # if the internal nodes are required, simply return the tree
        tree.write(outfile=output_path, format=1)
    else:
        tree.write(outfile=output_path, format=5) # write the rooted tree with no internal nodes names
    return 0


# generate alignments with INDELible based on given trees
def generate_alignment_and_tree(dataset_id, output_dir, priority, simulation_type, IR, A, RL, numberOfTaxa=100, tree_size=10, ultrametric=False, job_files_path="/groups/itay_mayrose/halabikeren/jobs/simulations/random_simulations/", error_files_path="/groups/itay_mayrose/halabikeren/error_files/random_simulations/", positive_selection_flag=False, positive_selection_model="alternative"):
    # set the simulation environment
    set_sim_env(output_dir)
    # generate a control file
    control_dir = output_dir + "control/" + dataset_id + "/"
    res = os.system("mkdir " + control_dir)
    dataset_output_dir = output_dir + dataset_id + "/"
    if not os.path.exists(dataset_output_dir):
        res = os.system("mkdir " + dataset_output_dir)
    control_txt = control_dir + "control.txt"
    create_control_file(control_dir, "random", ultrametric=ultrametric, simulation_type=simulation_type, indels_rate=IR, a_shape=A, root_seq_len=RL, numberOfTaxa=numberOfTaxa, random_tree_size=tree_size, positive_selection_flag=positive_selection_flag, ps_model=positive_selection_model)  # create a control file that generates an alignment according to a randomized tree
    # set the .sh file
    set_job_env(job_files_path, error_files_path)
    job_name = "run_indelible_" + dataset_id
    file_name = job_name + ".sh"
    commands = ["cd " + dataset_output_dir, "echo " + control_txt + " | /groups/pupko/haim/Programs/indelible/INDELibleV1.03/src/indelible"]
    touch_file_path = dataset_output_dir + "flag_indelible_done"
    full_job = create_job_file(job_name, commands, file_name, error_files_path, job_files_path, priority, 1,touch_file_path)
    res = os.system("qsub -p " + str(priority) + " " + full_job)
    return 0


# generate alignment alone with random trees with INDELible
def generate_alignment_from_tree(dataset_id, output_dir, tree_path, priority, simulation_type, IR, A, max_indel_length, RL, pinv, alpha, ngamcat, job_files_path="/groups/itay_mayrose/halabikeren/jobs/simulations/const_simulations/", error_files_path="/groups/itay_mayrose/halabikeren/error_files/const_simulations/", positive_selection_flag=False, positive_selection_model="alternative"):
    # set the simulation environment
    set_sim_env(output_dir)
    # generate a control file
    control_dir = output_dir + "control/" + dataset_id + "/"
    res = os.system("mkdir " + control_dir)
    dataset_output_dir = output_dir + dataset_id + "/"
    if not os.path.exists(dataset_output_dir):
        res = os.system("mkdir " + dataset_output_dir)
    control_txt = control_dir + "control.txt"
    create_control_file(control_dir, "constant", tree_path, simulation_type=simulation_type, indels_rate=IR, a_shape=A, max_indel_length=max_indel_length, root_seq_len=RL, pinv=pinv, alpha=alpha, ngamcat=ngamcat, positive_selection_flag=positive_selection_flag, ps_model=positive_selection_model)  # create a control file that generates an alignment according to a randomized tree
    # run INDELible via .sh file
    set_job_env(job_files_path, error_files_path)
    job_name = "run_indelible_" + dataset_id
    file_name = job_name + ".sh"
    commands = ["cd " + dataset_output_dir, "echo " + control_txt + " | /groups/pupko/haim/Programs/indelible/INDELibleV1.03/src/indelible"]
    touch_file_path = dataset_output_dir + "flag_indelible_done"
    full_job = create_job_file(job_name, commands, file_name, error_files_path, job_files_path, priority, 1, touch_file_path)
    res = os.system("qsub -p " + priority + " " + full_job)
    return 0


# generate alignment based on a given control file
def generate_sparta_based_alignment(dataset_id, dataset_output_dir, control_path, priority, job_files_path="/groups/itay_mayrose/halabikeren/jobs/simulations/const_simulations/", error_files_path="/groups/itay_mayrose/halabikeren/error_files/const_simulations/"):
    # run INDELible via .sh file
    set_job_env(job_files_path, error_files_path)
    job_name = "run_indelible_sparta_based" + dataset_id
    file_name = job_name + ".sh"
    commands = ["cd " + dataset_output_dir + "\n", "echo " + control_path + " | /groups/pupko/haim/Programs/indelible/INDELibleV1.03/src/indelible\n"]
    for commad in commands:
        res = os.system(commad)
    return 0


def rearrange_simulation_data(dir):
    for folder in os.listdir(dir):
        if folder not in output_names.values() and folder != "trees":
            dataset_output_dir = dir + folder + "/"
            # once the job is complete - change the output file names accordingly
            for file in os.listdir(dataset_output_dir):
                if file != "LOG.txt":
                    orig_file_name = dataset_output_dir + file
                    while not os.path.exists(orig_file_name):
                        print("orig_file_name: ", orig_file_name)
                        sleep(60)
                    new_file_name = dir + output_names[file] + "/" + folder + "_" + orig_to_new_names[file]
                    res = os.rename(orig_file_name, new_file_name)
            res = os.system("rm -rf " + dataset_output_dir)
    return 0




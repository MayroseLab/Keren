import os, sys
sys.path.append('/groups/itay_mayrose/halabikeren/myScripts/python/')
from utils.createJobFile import create_job_file

unaligned_genes_dir = "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/unaligned_genes/"
tree_path = "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/tree.nwk"
char_data_path = "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/char_data.fas"
alignment_jobs_dir = "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/aln_jobs/"
aligned_genes_dir = "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/aligned_genes/"
relax_param_dir = "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/relax_param/"
traitrelax_param_dir = "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/traitrelax_param/"

# generate a prank alignment (best for codons)
def run_prank(input_path, output_path, jobs_dir, job_name, guide_tree_path=None, codon=True, show_tree=True):
    cmd = "prank -d=" + input_path + " -o=" + output_path + " -F"
    if codon:
        cmd = cmd + " -codon"
    if guide_tree_path != None:
        cmd = cmd + " -t=" + guide_tree_path
    if show_tree:
        cmd = cmd + " -showtree"
    if not os.path.exists(output_path):
        output_dir = os.path.dirname(output_path)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

    file_name = job_name + ".sh"
    touch_file_path = job_name + "_flag_done"
    full_job = create_job_file(job_name, [cmd], file_name, jobs_dir, jobs_dir, 0, 1,
                               touch_file_path, queue="itaym1", limit_nodes=False, python=False, openmpi=False, language="bash")
    res = os.system("qsub " + full_job)

if __name__ == '__main__':

    if not os.path.exists(aligned_genes_dir):
        os.makedirs(aligned_genes_dir, exist_ok=True)
    if not os.path.exists(relax_param_dir):
        os.makedirs(relax_param_dir, exist_ok=True)
    if not os.path.exists(traitrelax_param_dir):
        os.makedirs(traitrelax_param_dir, exist_ok=True)
    if not os.path.exists(alignment_jobs_dir):
        os.makedirs(alignment_jobs_dir, exist_ok=True)

    # for each gene
    for filepath in os.listdir(unaligned_genes_dir):

        # replace "._" with "_" in the unaligned gene file to achieve match with names in the hyhyp tree
        with open(unaligned_genes_dir + filepath, "r") as infile:
            content = infile.read()
        content = content.replace("._", "_")
        with open(unaligned_genes_dir + filepath, "w") as outfile:
            outfile.write(content)

        gene_name = filepath.replace(".fa", "")

        # align the gene with prank via a job file, while using the species tree as a guide tree
        input_path = unaligned_genes_dir + filepath
        output_path = aligned_genes_dir + filepath
        # run_prank(input_path, output_path, alignment_jobs_dir, gene_name+"_prank")

        # create a parameters file for relax execution
        relax_template = '''# Global variables:
verbose = 1

# ----------------------------------------------------------------------------------------
#                                     Input alignment file
# ----------------------------------------------------------------------------------------

alphabet=Codon(letter=DNA)
genetic_code=Standard
input.sequence.file = <SEQ_DATA_PATH>
input.sequence.format = Fasta
input.sequence.sites_to_use = all
input.sequence.max_gap_allowed = 100%
input.sequence.remove_stop_codons = yes

# ----------------------------------------------------------------------------------------
#                                     Input tree file
# ----------------------------------------------------------------------------------------

init.tree = user
input.tree.file = <TREE_PATH>
input.tree.format = Newick
init.brlen.method = Input

# ----------------------------------------------------------------------------------------
#                                    Model specification
# ----------------------------------------------------------------------------------------

model1 = RELAX(kappa=2.0,p=0.5,omega1=1.0,omega2=2.0,k=1,theta1=0.5,theta2=0.8,frequencies=F3X4)
model2 = RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,frequencies=F3X4,1_Full.theta=RELAX.1_Full.theta_1,1_Full.theta1=RELAX.1_Full.theta1_1,1_Full.theta2=RELAX.1_Full.theta2_1,2_Full.theta=RELAX.2_Full.theta_1,2_Full.theta1=RELAX.2_Full.theta1_1,2_Full.theta2=RELAX.2_Full.theta2_1,3_Full.theta=RELAX.3_Full.theta_1,3_Full.theta1=RELAX.3_Full.theta1_1,3_Full.theta2=RELAX.3_Full.theta2_1,k=1.0)
nonhomogeneous = general
nonhomogeneous.number_of_models = 2
nonhomogeneous.stationarity = yes
site.number_of_paths = 2
site.path1 = model1[YN98.omega_1]&model2[YN98.omega_1]
site.path2 = model1[YN98.omega_2]&model2[YN98.omega_2]
rate_distribution = Constant() //Gamma(n=4, alpha=0.358)
likelihood.recursion = simple
likelihood.recursion_simple.compression = recursive

# ----------------------------------------------------------------------------------------
#                                    optimization parameters
# ----------------------------------------------------------------------------------------

optimization.tolerance = 0.000001
optimization.max_number_f_eval = 10000
optimization = FullD(derivatives=Newton,nstep=10)
optimization.final = powell

# ----------------------------------------------------------------------------------------
#                                    branches partition
# ----------------------------------------------------------------------------------------

model1.nodes_id = 0,1,2,3,4,5,6,9,10,11,12,13,14,15,16,17,18,19,20,23,24,25,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,88,89,90,91,92,93,94
model2.nodes_id = 7,8,21,22,26,27,28,29,30,31,32,33,78,79,80,81,82,83,84,85,86,87
'''
        relax_content = relax_template.replace("<SEQ_DATA_PATH>", output_path)
        relax_content = relax_content.replace("<TREE_PATH>", tree_path)
        with open(relax_param_dir + gene_name + ".bpp", "w") as relax_param:
            relax_param.write(relax_content)

        # create a parameters file for traitrelax execution
        traitrelax_template = '''# Global variables:
verbose = 1

# ----------------------------------------------------------------------------------------
#                                     Input character file
# ----------------------------------------------------------------------------------------

input.character.file = <CHAR_DATA_PATH>

# ----------------------------------------------------------------------------------------
#                                     Input alignment file
# ----------------------------------------------------------------------------------------

alphabet=Codon(letter=DNA)
genetic_code=Standard
input.sequence.file = <SEQ_DATA_PATH>
input.sequence.format = Fasta
input.sequence.sites_to_use = all
input.sequence.max_gap_allowed = 100%
input.sequence.remove_stop_codons = yes

# ----------------------------------------------------------------------------------------
#                                     Input tree file
# ----------------------------------------------------------------------------------------

init.tree = user
input.tree.file = <TREE_PATH>
input.tree.format = Newick
init.brlen.method = Input

# ----------------------------------------------------------------------------------------
#                                     Character Model specification
# ----------------------------------------------------------------------------------------

character_model.set_initial_parameters = true
character_model.mu = 2.0
character_model.pi0 = 0.5

# ----------------------------------------------------------------------------------------
#                                     Sequence Model specification
# ----------------------------------------------------------------------------------------

sequence_model.set_initial_parameters = true
model1 = RELAX(kappa=2.0,p=0.5,omega1=1.0,omega2=2.0,k=1,theta1=0.5,theta2=0.8,frequencies=F3X4)
model2 = RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,frequencies=F3X4,1_Full.theta=RELAX.1_Full.theta_1,1_Full.theta1=RELAX.1_Full.theta1_1,1_Full.theta2=RELAX.1_Full.theta2_1,2_Full.theta=RELAX.2_Full.theta_1,2_Full.theta1=RELAX.2_Full.theta1_1,2_Full.theta2=RELAX.2_Full.theta2_1,3_Full.theta=RELAX.3_Full.theta_1,3_Full.theta1=RELAX.3_Full.theta1_1,3_Full.theta2=RELAX.3_Full.theta2_1,k=1.0)

# ----------------------------------------------------------------------------------------
#                                    optimization parameters
# ----------------------------------------------------------------------------------------

optimization.tolerance = 0.000001
optimization.max_number_f_eval = 10000
optimization = FullD(derivatives=Newton,nstep=10)
optimization.final = powell
#optimization.message_handler = std
optimization.ignore_parameters = RELAX.k_1,BrLen
'''
        traitrelax_content = traitrelax_template.replace("<CHAR_DATA_PATH>", char_data_path)
        traitrelax_content = traitrelax_content.replace("<SEQ_DATA_PATH>", output_path)
        traitrelax_content = traitrelax_content.replace("<TREE_PATH>", tree_path)
        with open(traitrelax_param_dir + gene_name + ".bpp", "w") as traitrelax_param:
            traitrelax_param.write(traitrelax_content)


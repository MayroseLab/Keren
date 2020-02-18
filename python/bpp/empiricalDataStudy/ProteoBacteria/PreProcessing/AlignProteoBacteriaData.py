import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def translate(input_path, output_path):
    translated_records = []
    with open (input_path, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            translated_records.append(SeqRecord(record.translate(to_stop=True).seq, id=record.id, name=record.name, description=record.description))
    SeqIO.write(translated_records, output_path, "fasta")

def run_mafft(input_path, output_path):
    if not os.path.exists(output_path):
        res = os.system("mafft --localpair --maxiterate 1000 " + input_path + " > " + output_path)
        if res != 0:
            return -1
    return 0

def reverse_translate(codon_data_path, protein_data_path, output_path):
    cmd = "python /groups/itay_mayrose/halabikeren/programs/RevTrans-1.4/revtrans.py " + codon_data_path + " " + protein_data_path + " -Idna fasta -Ipep fasta > " + output_path
    res = os.system(cmd)
    with open(output_path, "r") as outfile:
        content = outfile.read()
    content = content.replace("DNA File type: fasta\n", "")
    content = content.replace("Pep File type: fasta\n", "")
    with open(output_path, "w") as outfile:
        outfile.write(content)

def filter_unrequired_pos(unfiltered_aligned_gene_path, unfiltered_aligned_protein_path, aligned_gene_path, aligned_protein_path, filtering_maps_path, cutoff=0.1):
    cmd = "python /groups/itay_mayrose/halabikeren/myScripts/python/utils/trimAlignment.py -ic " + unfiltered_aligned_gene_path + " -ia " + unfiltered_aligned_protein_path + " -oc " + aligned_gene_path + " -oa " + aligned_protein_path + " -om " + filtering_maps_path + " -c " + str(cutoff)
    res = os.system(cmd)

def set_relax_parameter_file(aligned_gene_path, tree_path, output_path):
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
optimization.scale.tree = 1

# ----------------------------------------------------------------------------------------
#                                    branches partition
# ----------------------------------------------------------------------------------------

model1.nodes_id = 1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,22,23,24,25,26,27,28,29,30,31,32,33,34,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,91,92,93,94,95,96,0
model2.nodes_id = 7,19,20,21,35,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90
'''
    relax_content = relax_template.replace("<SEQ_DATA_PATH>", aligned_gene_path)
    relax_content = relax_content.replace("<TREE_PATH>", tree_path)
    with open(output_path, "w") as relax_param:
        relax_param.write(relax_content)

def set_traitrelax_parameter_file(char_data_path, aligned_gene_path, tree_path, output_path):
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
optimization.scale.tree = 1
character.use_analytic_mapping = 1
optimization.advanced = 1
    '''
    traitrelax_content = traitrelax_template.replace("<CHAR_DATA_PATH>", char_data_path)
    traitrelax_content = traitrelax_content.replace("<SEQ_DATA_PATH>", aligned_gene_path)
    traitrelax_content = traitrelax_content.replace("<TREE_PATH>", tree_path)
    with open(output_path, "w") as traitrelax_param:
        traitrelax_param.write(traitrelax_content)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='aligns proteobacteria empirical data and creates the parameter files for relax and traitrelax')
    parser.add_argument('--unaligned_genes_dir', '-i', help='directory that holds the simulation output', required=True)
    parser.add_argument('--aligned_genes_dir', '-o', help='directory to hold the aligned codon sequences of the genes', required=True)
    parser.add_argument('--help_dir', '-a', help='directory that holds the by products of the script', required=True)
    parser.add_argument('--char_data_path', '-c', help='path to the character trait data', required=True)
    parser.add_argument('--tree_path', '-t', help='path to the input tree', required=True)
    parser.add_argument('--relax_parameters_dir', '-r', help='directory to hold the parameter files for relax', required=True)
    parser.add_argument('--traitrelax_parameters_dir', '-tr', help='directory to hold the parameter files for traitrelax', required=True)

    args = parser.parse_args()
    unaligned_genes_dir = args.unaligned_genes_dir # "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/unaligned_genes/"
    aligned_genes_dir = args.aligned_genes_dir
    if not os.path.exists(aligned_genes_dir):
        os.makedirs(aligned_genes_dir, exist_ok=True)
    help_dir = args.help_dir
    if not os.path.exists(help_dir):
        os.makedirs(help_dir, exist_ok=True)
    unaligned_protein_dir = help_dir + "unaligned_proteins/"
    if not os.path.exists(unaligned_protein_dir):
        os.makedirs(unaligned_protein_dir, exist_ok=True)
    unfiltered_aligned_protein_dir = help_dir + "unfiltered_aligned_proteins/"
    if not os.path.exists(unfiltered_aligned_protein_dir):
        os.makedirs(unfiltered_aligned_protein_dir, exist_ok=True)
    unfiltered_aligned_genes_dir = help_dir + "unfiltered_aligned_genes/"
    if not os.path.exists(unfiltered_aligned_genes_dir):
        os.makedirs(unfiltered_aligned_genes_dir, exist_ok=True)
    filtered_aligned_protein_dir = help_dir + "filtered_aligned_proteins/"
    if not os.path.exists(filtered_aligned_protein_dir):
        os.makedirs(filtered_aligned_protein_dir, exist_ok=True)
    filtering_maps_dir = help_dir + "filtering_maps/"
    if not os.path.exists(filtering_maps_dir):
        os.makedirs(filtering_maps_dir, exist_ok=True)
    tree_path = args.tree_path # "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/tree.nwk"
    char_data_path = args.char_data_path # "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/char_data.fas"

    relax_param_dir = args.relax_parameters_dir # "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/relax_param/"
    if not os.path.exists(relax_param_dir):
        os.makedirs(relax_param_dir, exist_ok=True)
    traitrelax_param_dir = args.traitrelax_parameters_dir # "/groups/itay_mayrose/halabikeren/TraitRELAX/real_data/ProteoBacteria/traitrelax_param/"
    if not os.path.exists(traitrelax_param_dir):
        os.makedirs(traitrelax_param_dir, exist_ok=True)

    # for each gene
    for filepath in os.listdir(unaligned_genes_dir):

        # replace "._" with "_" in the unaligned gene file to achieve match with names in the hyphy tree
        with open(unaligned_genes_dir + filepath, "r") as infile:
            content = infile.read()
        content = content.replace("._", "_")
        with open(unaligned_genes_dir + filepath, "w") as outfile:
            outfile.write(content)

        gene_name = filepath.replace(".fa", "")
        print("#################################################")
        print("gene: ", gene_name)

        # translate the coding sequences to protein sequences
        print("\n\ntranslating gene to protein....\n")
        unaligned_gene_path = unaligned_genes_dir + filepath
        unaligned_protein_path = unaligned_protein_dir + filepath
        translate(unaligned_gene_path, unaligned_protein_path)

        # align the protein sequences
        print("\n\naligning protein sequences....\n")
        unfiltered_aligned_protein_path = unfiltered_aligned_protein_dir + filepath
        run_mafft(unaligned_protein_path, unfiltered_aligned_protein_path)

        # reverse translate back to alignment of coding sequences
        print("\n\nreverse translating alignments....\n")
        unfiltered_aligned_gene_path = unfiltered_aligned_genes_dir + filepath
        reverse_translate(unaligned_gene_path, unfiltered_aligned_protein_path, unfiltered_aligned_gene_path)

        # filter out positions with
        print("\n\nfiltering out non informative positions....\n")
        aligned_protein_path = filtered_aligned_protein_dir + filepath
        aligned_gene_path = aligned_genes_dir + filepath
        filtering_maps_path = filtering_maps_dir + gene_name
        filter_unrequired_pos(unfiltered_aligned_gene_path, unfiltered_aligned_protein_path, aligned_gene_path, aligned_protein_path, filtering_maps_path)

        # create a parameters file for relax execution
        print("\n\ncreating relax parameter file....\n")
        output_path = relax_param_dir + gene_name + ".bpp"
        set_relax_parameter_file(aligned_gene_path, tree_path, output_path)

        # create a parameters file for traitrelax execution
        print("\n\ncreating traitrelax parameter file....\n")
        output_path = traitrelax_param_dir + gene_name + ".bpp"
        set_traitrelax_parameter_file(char_data_path, aligned_gene_path, tree_path, output_path)

        print("#################################################")



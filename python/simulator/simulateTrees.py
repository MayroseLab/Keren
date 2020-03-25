import argparse, os, re
from ete3 import Tree
from time import sleep

def create_control_file(taxa_num, birth_rate, death_rate, sampling_fraction, mutation_rate, rooted, replicates_num, outpath, ultrametric=True, tree_depth=1, tree_template=""):
    ultrametric_template = '''[TYPE] NUCLEOTIDE 1	

[MODEL]    modelname
  [submodel]     JC
  [indelmodel]   NB  0.4 1
  [insertrate]   0.08
  [deleterate]   0.12
  
[TREE] tree
  [<type>] <taxa_num> <birth_rate> <death_rate> <sampling_fraction> <mutation_rate>     // ntaxa birth death sample mut
  [treedepth] <tree_depth>
  
[PARTITIONS] partitionname
  [tree modelname 1]

[EVOLVE] partitionname <replicates_num> out
'''

    non_ultrametric_template = '''[TYPE] NUCLEOTIDE 1	

[MODEL]    modelname
  [submodel]     JC
  [indelmodel]   NB  0.4 1
  [insertrate]   0.08
  [deleterate]   0.12
  
[TREE] tree
  // No branch lengths need to be provided  
  <tree_template>

  [branchlengths] NON-ULTRAMETRIC // All branch lengths will be equal
  
[PARTITIONS] partitionname
  [tree modelname 1]

[EVOLVE] partitionname <replicates_num> out
'''

    type = "unrooted"
    if rooted:
        type = "rooted"

    # set values
    template = ultrametric_template
    if not ultrametric:
        template = non_ultrametric_template

    control_content = template.replace("<type>", type)
    control_content = control_content.replace("<taxa_num>", str(taxa_num))
    control_content = control_content.replace("<birth_rate>", str(birth_rate))
    control_content = control_content.replace("<death_rate>", str(death_rate))
    control_content = control_content.replace("<sampling_fraction>", str(sampling_fraction))
    control_content = control_content.replace("<mutation_rate>", str(mutation_rate))
    control_content = control_content.replace("<tree_depth>", str(tree_depth))
    control_content = control_content.replace("<replicates_num>", str(replicates_num))
    control_content = control_content.replace("<tree_template>", tree_template)
    with open(outpath, "w") as out:
        out.write(control_content)

def extract_trees(output_dir, input_path, tree_size, custom_name=""):
    tree_regex = re.compile("out\s*tree\s*[^\(]*([^\n]*?)\n", re.MULTILINE | re.DOTALL)

    with open(input_path, "r") as input_file:
        input = input_file.read()

        trees = tree_regex.finditer(input)
        trees_count = 0
        for tree in trees:
            output_path = output_dir + str(trees_count) + ".nwk"
            if custom_name != "":
                output_path = output_dir + custom_name
            with open(output_path, "w") as output:
                output.write(tree.group(1))
            trees_count += 1

    res = os.system("rm -r " + input_path)

    # set node names with "S"
    for tree_path in os.listdir(output_dir):
        if (custom_name == "" and ".nwk" in tree_path) or (custom_name != "" and tree_path == custom_name):
            tree = Tree(output_dir + tree_path)
            orig_tree_size = 0
            for node in tree.traverse():
                orig_tree_size += node.dist
                if node.is_leaf():
                    node_name = "S" + node.name
                    node.name = node_name

            print("rescaling factor: ", tree_size / orig_tree_size)
            for node in tree.traverse():
                node.dist = node.dist * tree_size / orig_tree_size

            tree.write(outfile=output_dir + tree_path, format=5)


def simulate_non_ultrametric_trees(taxa_num, birth_rate, death_rate, sampling_fraction, mutation_rate, rooted,
                        replicates_num, output_dir, template_trees_path, tree_size):

    print("simulate_non_ultrametric_trees called ")
    # collect tree templates (without branch lengths)
    tree_templates = []
    tree_regex = re.compile("out\s*tree\s*[^\(]*([^\n]*?)\n", re.MULTILINE | re.DOTALL)

    with open(template_trees_path, "r") as input_file:
        input = input_file.read()

        trees = tree_regex.finditer(input)
        for tree in trees:
            full_tree_str = tree.group(1)
            tree = Tree(full_tree_str)
            topology_str = tree.write(outfile=None, format=9)
            tree_templates.append(topology_str)

    # for each tree template, simulate a single tree and extract it. then, rename tree extracted file to maintain numbering
    trees_counter = 0
    for tree_template in tree_templates:

        if trees_counter > replicates_num:
            return

        # create control file
        create_control_file(taxa_num, birth_rate, death_rate, sampling_fraction, mutation_rate, rooted, 1, output_dir+"control.txt", ultrametric=False, tree_depth=1, tree_template=tree_template)

        # execute INDELible
        res = os.chdir(output_dir)
        res = os.system("/groups/itay_mayrose/halabikeren/indelible/INDELibleV1.03/src/indelible")

        # # delete redundant files
        # res = os.system("rm -r " + output_dir + "LOG.txt")
        # res = os.system("rm -r " + output_dir + "out_TRUE.phy")
        # res = os.system("rm -r " + output_dir + "out.fas")
        # res = os.system("rm -r " + output_dir + "control.txt")

        # extract trees from the tree file
        extract_trees(output_dir, output_dir+"trees.txt", tree_size, custom_name="temp_tree.nwk")

        # rename single tree file
        res = os.rename(output_dir+"temp_tree.nwk", output_dir+str(trees_counter)+".nwk")
        trees_counter += 1

    return

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(description='simulates trees with INDELible')
    parser.add_argument('--output_dir', '-o', help='directory that holds the simulation output',
                        required=True)
    parser.add_argument('--taxa_num', '-tn', help='number of taxa to simulate in tree', required=True)
    parser.add_argument('--birth_rate', '-br', help='birth rate to use', required=False, default=0.3)
    parser.add_argument('--death_rate', '-dr', help='death rate to use', required=False, default=0.1)
    parser.add_argument('--mutation_rate', '-mr', help='mutation rate to use', required=True)
    parser.add_argument('--sampling_fraction', '-sf', help='death rate to use', required=False, default=0.25)
    parser.add_argument('--rooted', '-r', help='true is trees need to be rooted, else false', required=False, default=True)
    parser.add_argument('--replicates_num', '-n', help='number of trees to generate', required=False, default=100)
    parser.add_argument('--ultrametric', '-u', help='boolean: if True simulates ultrametric trees, else, non-utlrametric', required=False, default=True)
    parser.add_argument('--tree_size', '-ts', help='used in case of non-ultrametric trees as the final tree size', required=False, default=0)
    args = parser.parse_args()
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        res = os.system("mkdir -p " + output_dir)
    taxa_num = int(args.taxa_num)
    birth_rate = float(args.birth_rate)
    death_rate = float(args.death_rate)
    mutation_rate = float(args.mutation_rate)
    sampling_fraction = float(args.sampling_fraction)
    rooted = bool(args.rooted)
    replicates_num = int(args.replicates_num)
    ultrametric = bool(int(args.ultrametric))
    print("ultrametric: ", ultrametric)
    tree_size = float(args.tree_size)
    print("tree_size: ", tree_size)
    # create the control file in the output dir
    create_control_file(taxa_num, birth_rate, death_rate, sampling_fraction, mutation_rate, rooted,
                        replicates_num, output_dir+"control.txt")
    # execute INDELible
    res = os.chdir(output_dir)
    res = os.system("/groups/itay_mayrose/halabikeren/programs/indelible/INDELibleV1.03/src/indelible")

    # delete redundant files
    res = os.system("rm -r " + output_dir + "LOG.txt")
    res = os.system("rm -r " + output_dir + "out_TRUE.phy")
    res = os.system("rm -r " + output_dir + "out.fas")
    res = os.system("rm -r " + output_dir + "control.txt")

    # if not ultrametric, read each file without branch length and simulate not ultrametric tree based on it
    if not ultrametric:
        simulate_non_ultrametric_trees(taxa_num, birth_rate, death_rate, sampling_fraction, mutation_rate, rooted, replicates_num, output_dir, output_dir+"trees.txt", tree_size)
    else:
        # extract trees from the tree file
        extract_trees(output_dir, output_dir+"trees.txt", tree_size)



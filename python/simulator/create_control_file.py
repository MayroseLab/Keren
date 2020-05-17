import re

def create_control_file(control_dir, numberOfTaxa, numberOfTrees, submodel, dest, tree, source, indels_rate):
	'''generates conrol.txt file in indelDir according to the given parameters'''
	control = control_dir + "control.txt"
	with open(control, 'w+') as handle:
		handle.write('[TYPE] AMINOACID 2\n')
		handle.write('\n')
		handle.write('[SETTINGS]\n')
		handle.write('\t[ancestralprint]	NEW\n')
		handle.write('\t[phylipextension] phy\n')
		handle.write('\t[output]          PHYLIP\n')
		handle.write('\t[fileperrep]      TRUE\n')
		handle.write('\t[printrates] TRUE     // FALSE or TRUE\n')
		handle.write('\n')
		handle.write('[MODEL] tree_model\n')
		handle.write('\t[submodel] ' + submodel + '\n')
		handle.write('\t[rates] 0 0.5 16      //   16 category discrete gamma with alpha=2\n')
		handle.write('\t[indelmodel] POW 1.7 30\n')  #need to ask about the a=1.5 parameter
		handle.write('\t[indelrate] ' + indels_rate + ' \n')
		handle.write('\n')
		if (dest == 'to_trees'):
			handle.write('[TREE] tree\n')
			handle.write('\t[unrooted] ' + numberOfTaxa + ' 6.7 2.5 0.234 0.31\n')
			handle.write('\t[seed] 1000 \n')
			handle.write('\t[treelength] 100\n')
		elif (dest == 'to_alignments'):
			tree_str = adjust_tree(tree, source)
			handle.write('[TREE] tree ' + tree_str)
			handle.write('\n')
		handle.write('\n')
		handle.write('[PARTITIONS] job\n')
		handle.write('\t[tree tree_model 350]\n')
		handle.write('\n')
		handle.write('[EVOLVE] job ' + str(numberOfTrees) + ' output\n')
	return



def reformat_num(number):
	num = float(number)
	reformatted_num = str("{0:.10f}".format(num))
	return reformatted_num

def reformat_branch(branch):
	fixed_branch = branch[0] + ":" + branch[1:]
	return fixed_branch

def adjust_tree(tree, source):

    tree_file = open(tree, 'r')
    tree_str = tree_file.read()
    tree_file.close()

    if source == "random": # if the source of the tree is random (created by INDELiable)
        # if there is ROOT string in the end of the tree -> remove it
        tree_str = tree_str.replace('ROOT', '')

        # remove internal node names
        internal_nodes_pattern = 'N[0-9]+'
        pattern = re.compile(internal_nodes_pattern)
        tree_str = pattern.sub('', tree_str)

        # reformat small numbers
        reformatted_number = re.compile(r'\d+[\.\d+]*e\-\d+')
        for number in reformatted_number.findall(tree_str):
            reformatted_num = reformat_num(number)
            tree_str = tree_str.replace(number, reformatted_num, 1)

    elif source == "const": # if the source of the tree is FastTree Tests

        # reformat the branches
        comb_branch = '\)[0-9]+\.[0-9]+'
        comb_branch_pattern = re.compile(comb_branch)
        for branch in re.findall(comb_branch_pattern, tree_str):
            reformatted_branch = reformat_branch(branch)
            tree_str = re.sub(re.escape(branch), reformatted_branch, tree_str)

        # reformat small numbers
        reformatted_number = re.compile(r'\d+[\.\d+]*e\-\d+')
        for number in reformatted_number.findall(tree_str):
            print("number: " + number) #debug!!
            reformatted_num = reformat_num(number)
        tree_str = tree_str.replace(number, reformatted_num, 1)


        #remove bootsrap numbers
        bootsrap_num_expr = re.compile(r"-{0,1}[0-9]+\.[0-9]+:")
        tree_str = bootsrap_num_expr.sub("", tree_str)

    return tree_str



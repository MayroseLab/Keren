/* ____________ imports ____________ */

LoadFunctionLibrary("libv3/all-terms.bf"); 				  // must be loaded before CF3x4
LoadFunctionLibrary("libv3/UtilityFunctions.bf"); 		  // namespace 'utility' for convenience functions
LoadFunctionLibrary("libv3/IOFunctions.bf");			  // namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary ("libv3/models/codon/MG_REV.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/tasks/ancestral.bf");
LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");
LoadFunctionLibrary("libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");
LoadFunctionLibrary("libv3/models/binary/charbinary.bf"); 			  			// character model settings 
LoadFunctionLibrary("libv3/models/binary/empirical.bf"); 			  			// character model settings - used to allow setting frequencies as MLEs 
LoadFunctionLibrary("libv3/models/binary.bf");
LoadFunctionLibrary("/scratch300/halabikeren/myScripts/HBL/modules/custom_functions.bf");
LoadFunctionLibrary("libv3/models/model_functions.bf");
LoadFunctionLibrary("libv3/convenience/regexp.bf");
LoadFunctionLibrary("libv3/all-terms.bf");
utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
ExecuteCommands ('LoadFunctionLibrary("chooseGeneticCode");',
					 {"0" : "Universal"});	
utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);	
LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");

/* ____________ hardcoded parameters ____________ */

histories_num = 100;
max_simulations_num = 10000;

/* ____________ main ____________ */

/* process input */
trait_relax.json = { terms.json.input: {},
					  terms.json.fits : {},
					  terms.json.timers : {},
					  terms.json.test_results : {}
					  };

trait_relax.process_input();

/* compute the tree's size */
tree_size = compute_tree_size(trait_relax.trait_trees[0]);

// declare character model parameters
terms.trait_relax.mu = "character transition rate parameter";
mp_data = trait_relax.set_mp_partition(trait_relax.trait_trees["0"], &trait_relax.trait_dataset, trait_relax.char_0.model_branches_name, trait_relax.char_1.model_branches_name);
mp_score = mp_data["0"];
fprintf(stdout, "parsimoniuos number of transitions: ", mp_score, "\n"); // debug
mp_model_map = {"trait_relax.char_1.model" : utility.Filter (mp_data["1"], '_value_', '_value_ == trait_relax.char_1.model_branches_name'),
				"trait_relax.char_0.model" : utility.Filter (mp_data["1"], '_value_', '_value_ == trait_relax.char_0.model_branches_name')};
trait_relax.initial_mu = mp_score/tree_size;
terms.trait_relax.mu_range = {
		terms.lower_bound: (mp_score/tree_size), // shoudn't be the score - need to divide it by the size of the tree - no?
		terms.upper_bound: (2*mp_score/tree_size)
	};
	
terms.trait_relax.pi_0 = "frequency of charcter state 0";
HarvestFrequencies(trait_relax.char_empirical_frequencies, ^trait_relax.trait_filter_names["0"], 1, 1, 1);
trait_relax.initial_pi_0 = trait_relax.char_empirical_frequencies[0][0];
terms.trait_relax.pi_0_range = {
		terms.lower_bound: 0.01, // do not allow value 0 in order to avoid observing state in the character transition matrix
		terms.upper_bound: 0.99  // do not allow value 1 in order to avoid observing state in the character transition matrix 
	};
	
trait_relax.mu = trait_relax.initial_mu;
fprintf(stdout, "trait_relax.mu: ", trait_relax.mu, "\n");
trait_relax.pi_0 = trait_relax.initial_pi_0;
fprintf(stdout, "trait_relax.pi_0: ", trait_relax.pi_0, "\n\n");


/* debug - print the expected number of transitions */
expected_transitions_num = tree_size * trait_relax.mu;
expected_01_transitions_num = trait_relax.pi_0 * expected_transitions_num;
expected_10_transitions_num = expected_transitions_num - expected_01_transitions_num;
fprintf(stdout, "expected_transitions_num: ", expected_transitions_num, "\n");
fprintf(stdout, "expected_01_transitions_num: ", expected_01_transitions_num, "\n");
fprintf(stdout, "expected_10_transitions_num: ", expected_10_transitions_num, "\n\n");

trait_relax_histories_paths = {};
histories_durations_dict = {};

transitions_counter = {"mean_transitions_num" : 0,
					   "mean_01_transitions_num" : 0,
					   "mean_10_transitions_num" : 0};

for (hist = 0; hist < histories_num; hist = hist+1) {  
	 
	// generate a stochastic mapping and put it in the path labelledtree_path	
	history_data = custom_functions.simulate_history(hist, &histories_durations_dict, trait_relax.trait_trees[0], trait_relax.pi_0, trait_relax.mu, max_simulations_num, {"write_histories": TRUE, "histories_dir": histories_dir, "debug_mode": TRUE});
	history = history_data["history"];
	root_label = history_data["root_label"];
	count_transitions(history, root_label, transitions_counter);
}

transitions_counter["mean_transitions_num"] = transitions_counter["mean_transitions_num"] / histories_num;
transitions_counter["mean_01_transitions_num"] = transitions_counter["mean_01_transitions_num"] / histories_num;
transitions_counter["mean_10_transitions_num"] = transitions_counter["mean_10_transitions_num"] / histories_num;
fprintf(stdout, "mean_transitions_num: ", transitions_counter["mean_transitions_num"], "\n");
fprintf(stdout, "mean_01_transitions_num: ", transitions_counter["mean_01_transitions_num"], "\n");
fprintf(stdout, "mean_10_transitions_num: ", transitions_counter["mean_10_transitions_num"], "\n");

/* ____________ auxiliary functions ____________ */

lfunction count_transitions(history, root_label, transitions_counter)
{
	labels_dict = history[utility.getGlobalValue("terms.trees.model_map")];
	node_to_parent = custom_functions.mapNodeToParent(history);
	node_to_children = custom_functions.mapNodeToChildren(history, node_to_parent);
	nodes = trees.BranchNames(history);
	for (i=0; i < Columns(nodes); i += 1) {
		node_label = labels_dict[nodes[i]];
		if (nodes[i] % "Node0") {
			node_label = root_label;
		}
		is_leaf = (history[utility.getGlobalValue("terms.trees.partitioned")])[nodes[i]] % "leaf";
		if (!is_leaf) {
			children_num = (node_to_children[nodes[i]])["childrenNum"];
			children_names = (node_to_children[nodes[i]])["childrenNames"];
			for (j=0; j < children_num; j += 1) {
				child_name = children_names[j];
				child_label = labels_dict[child_name];
				if (node_label % child_label == 0) {
					transitions_counter["mean_transitions_num"] += 1;
					if (node_label % "BG" && child_label % "FG") {
						transitions_counter["mean_01_transitions_num"] += 1;
					} else {
						transitions_counter["mean_10_transitions_num"] += 1;
					}
				}
			}
		}
	}
}

lfunction compute_tree_size(tree) 
{
	tree_size = 0;
	branch_lengths_dict = tree[utility.getGlobalValue("terms.branch_length")];
	branches_names = utility.Keys(branch_lengths_dict);
	for (i=0; i < Abs(branch_lengths_dict); i += 1) {
		branch_name = branches_names[i];
		tree_size += branch_lengths_dict[branch_name];
	}
	return tree_size;
}

lfunction trait_relax.set_base_partition(partition_info) 
{
	tree = ((partition_info["0"])[utility.getGlobalValue("terms.data.tree")]); // [utility.getGlobalValue("terms.trees.newick_with_lengths")];
	branch_names = trees.BranchNames(tree);
	model_map = {};
	char_0_dict = {};
	for (i=0; i<Columns(branch_names)-1; i=i+1) {
		branch = branch_names[i];
		char_0_dict[branch] = "char_0";
	}
	return {"0": char_0_dict};
}	

lfunction trait_relax.set_mp_partition(tree, trait_data, char_0_branches_names, char_1_branches_names) 
{	
	state_to_label = {"0": char_0_branches_names , "1": char_1_branches_names};
	mp_data = custom_functions.getMPTreePartition(tree, state_to_label, trait_data); // move to to he null section and change the branches selector accordingly
	mp_score = mp_data[0];
	mp_selected_branches = mp_data[1];
	return {"0": mp_score, "1": mp_selected_branches};
}

lfunction trait_relax.generate_character_likelihood_function(mu, pi_0, lf_prefix, data_filter, tree) 
{
	character_initial_parameters = {
	 "global":{
	   "Substitution rate from character 0 to character 1":
	   {utility.getGlobalValue("terms.fit.MLE"): mu}
	  },
	 utility.getGlobalValue("terms.efv_estimate"):{
	   utility.getGlobalValue("terms.model"):{
		{pi_0}
		{1-pi_0}
		}
	  }
	};
	
	/* for improvement: since the empirical frequencies are estimates by default, the set_freq option runs over the empirical frequencies 
	with the values provided in the initial values */
	trait_likelihood_function = custom_functions.CreateLFObject_FixedBLs (lf_prefix, data_filter, tree, "models.binaryML.ModelDescription", character_initial_parameters, {utility.getGlobalValue("terms.run_options.retain_lf_object"): TRUE, "set_freq": TRUE, "return_lf": TRUE});
	
	return trait_likelihood_function
}

function trait_relax.process_input() 
{
	
	// accept input paths from the user
	trait_relax.trait_data_path = io.PromptUserForString("\n>Provide full path to the trait character alignment");
	trait_relax.tree_path =  io.PromptUserForString("\n>Provide full path to phylogenetic tree in newick format");
	trait_relax.seq_data_path = io.PromptUserForString("\n>Provide full path to codon sequence alignment");
	trait_relax.json_path = outDir + "simulation_output.json";
	
	trait_relax.branch_selector_function = "trait_relax.set_base_partition";
	namespace trait_relax { 
		// process the tree input files (sequence alignment, tree and trait data)
		LoadFunctionLibrary ("SelectionAnalyses/modules/TraitRELAX_aux.bf"); // used for loading function the below function inside the scope of namespace "trait_relax" (-> all the instances generated within will be inside the namespace)
		load_input ({utility.getGlobalValue("terms.prefix"): "trait_relax", utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : branch_selector_function}}, json_path, seq_data_path, tree_path, trait_data_path); // this function also validates that the species in the tree are exactly the same as the ones in the alignment
	}
	
	// validate name compatibility in character data and codon data. The validation with the names in the tree is performed earlier (function: load_input)
	orig_seq_names = Rows (trait_relax.name_mapping);
	orig_char_names = Rows (trait_relax.trait_name_mapping);
	len_res = Columns(orig_seq_names) == Columns(orig_char_names);
	io.CheckAssertion("len_res==1","number of species in the trait alignment doesn't match the number of species in the codon alignment");
	for (i = 0; i < Columns (orig_seq_names); i += 1) {
		converted_seq_name = trait_relax.name_mapping[orig_seq_names[i]];
		converted_char_name = trait_relax.trait_name_mapping[orig_seq_names[i]];
		comp_res = converted_seq_name%converted_char_name;
		io.CheckAssertion("comp_res==1","inconsistent species names in the trait alignment and codon alignment");
	}

}

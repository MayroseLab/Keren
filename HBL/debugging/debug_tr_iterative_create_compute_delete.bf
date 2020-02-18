RequireVersion("2.3.3");

/* _______________________________ IMPORTS _______________________________ */

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
LoadFunctionLibrary("modules/custom_functions.bf");


/* ________________________ ENVIRONMENT VARIABLES ________________________ */

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);	// defined based on RELAX.bf
utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);	// defined based on RELAX.bf
utility.SetEnvVariable ("REPLACE_TREE_STRUCTURE", TRUE); 	// assures initialization of the tree every time it is re-loaded - required for loading histories repeatedly, with the same name
utility.SetEnvVariable("ACCEPT_ROOTED_TREES", TRUE);		// allows to preserve rooted trees in the generation of a HyPhy tree instance

/* ________________________ HARDCODED PARAMETERS ________________________ */

/* ___________________ HARDCODED PARAMETERS FOR DEBUGGING _______________ */
mu_interval = 0.1;
pie_0_interval = 0.1;
histories_dir = "/groups/itay_mayrose/halabikeren/myScripts/HBL/debug_mem_output/";
// character histories sampling properties
maxSimulationsNum = 100000;

// numeric seetings
almostZero = 1.0e-10;

/* ________________________________ MAIN ________________________________ */

/* initiate script parameters */
trait_relax.init();

/* process user input */
trait_relax.process_input();

/* create the sequence model components (test and reference) and perform preliminary fitting for initialization of the model parameters */
trait_relax.set_model_components();

// generate a basic model map that maps all branches into the refrence category, for the independent case tree
base_model_map = trait_relax.get_base_model_map(trait_relax.trait_trees["0"]); 

trait_relax.create_compute_delete();

/* ___________________ MAIN AUXILIARY FUNCTIONS ___________________ */

		
/**
 * @name trait_relax.init
 * @usage sets the initial parameters for the pipeline
 */
function trait_relax.init() 
{
		
	/* initialization of documentation json file */
	trait_relax.json = { terms.json.input: {},
						  terms.json.fits : {},
						  terms.json.timers : {},
						  terms.json.test_results : {}
						  };
	
	// declare model parameters for json documentation
	trait_relaxation_parameter = "trait_relax.K";
	trait_dependent_site_proportion = "trait_relax.p";
	trait_relax.rate_classes = 3;
	trait_relax.MG94_name = terms.json.mg94xrev_sep_rates;
	trait_relax.alternative.name = "TraitRELAX alternative";
	trait_relax.null.name = "TraitRELAX null";
	trait_relax.trait_charachter_substitution_rate = "trait_relax.trait_substitution_rate"; // == mu
	trait_relax.null.codon_substitution_rate = "trait_relax.null.substitution_rate";		// == trait_relax.scaler before alternative optimization
	trait_relax.alternative.codon_substitution_rate = "trait_relax.null.substitution_rate"; // == trait_relax.scaler after alternative optimization


	// set the information display order in the json file
	trait_relax.display_orders = {terms.original_name: -1,
							terms.json.nucleotide_gtr: 0,
							trait_relax.MG94_name: 1,
							trait_relax.alternative.name: 2,
							trait_relax.null.name: 3,
							trait_relax.trait_charachter_substitution_rate: 4,
							trait_relax.null.codon_substitution_rate: 5, 
							trait_relax.alternative.codon_substitution_rate : 6,
							trait_relaxation_parameter: 7,
							trait_dependent_site_proportion: 8
						   };	
						   
	// start overall time measurement
	selection.io.startTimer (trait_relax.json [terms.json.timers], "Overall", 0);
}


/**
 * @name trait_relax.process_input
 * @usage requests input paths from the user and processes them into hyphy instances
 */
function trait_relax.process_input() 
{
	
	/* accept input paths from the user */
	trait_relax.trait_data_path = io.PromptUserForString("\n>Provide full path to the trait character alignment");
	trait_relax.tree_path =  io.PromptUserForString("\n>Provide full path to phylogenetic tree in newick format");
	trait_relax.seq_data_path = io.PromptUserForString("\n>Provide full path to codon sequence alignment");
	trait_relax.json_path = "debug";
	
	namespace trait_relax { 
		/* process the tree input files (sequence alignment, tree and trait data) */
		LoadFunctionLibrary ("modules/TraitRELAX_aux.bf"); // used for loading function the below function inside the scope of namespace "trait_relax" (-> all the instances generated within will be inside the namespace)
		load_input ({utility.getGlobalValue("terms.prefix"): "trait_relax", utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "trait_relax.set_mp_partition"}}, json_path, seq_data_path, tree_path, trait_data_path); // this function also validates that the species in the tree are exactly the same as the ones in the alignment
	}
	
	/* validate name compatibility in character data and codon data. The validation with the names in the tree is performed earlier (function: load_input) */
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


/**
 * @name trait_relax.declare_model_parameters
 * @usage sets the model parameters that extend the BS-REL model describing the test and reference instances: the selection intensity parameter k and the dependent site proportion p
 */
function trait_relax.declare_model_parameters() 
{
	
	// declare the selection intensity parameter k
	terms.trait_relax.k          = "relaxation or intensification parameter";
	terms.trait_relax.k_range    = {
			terms.lower_bound: "0",
			terms.upper_bound: "50"
		};
		
	// declare mixture parameter p for probability of each model
	terms.trait_relax.p = "trait dependent site proportion";
	parameters.DeclareGlobal(trait_dependent_site_proportion, None);
	trait_relax.p = 1;	// initial p is zero (same as null - all sites are independent of the character history)
	trait_relax.p :< 1; // hyphy bounds by default each parameter value to be in [0,1000], so p cannot receive a negative value 

	trait_relax.test_branches_name = "Test";
	trait_relax.reference_branches_name = "Reference";
}

/**
 * @name trait_relax.set_codon_model
 * @usage generates the codon model instance as a set of two sub-instances: test and reference, which are extensions of the BS-REL model
 */
function trait_relax.set_codon_model() 
{
	/* define the codon model as an extension of the BS-REL model, where each branch can have one of several (3) w values */

	// set the test model as an extension of the BS-REL model
	trait_relax.test =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
			"trait_relax.test", {
				"0": parameters.Quote(terms.global),
				"1": trait_relax.codon_data_info[terms.code],
				"2": parameters.Quote (trait_relax.rate_classes) // the number of rate classes
			},
			trait_relax.filter_names,
			None);
		
	// set the reference model as an extension of the BS-REL model
	trait_relax.reference =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
			"trait_relax.reference", {
				"0": parameters.Quote(terms.global),
				"1": trait_relax.codon_data_info[terms.code],
				"2": parameters.Quote (trait_relax.rate_classes) // the number of rate classes
			},
			trait_relax.filter_names,
			None);
			
	// restrict the model parameters to be the same for the test and reference
	trait_relax.bound_weights = models.BindGlobalParameters ({"0" : trait_relax.reference, "1" : trait_relax.test}, terms.mixture.mixture_aux_weight + ".+");
	models.BindGlobalParameters ({"0" : trait_relax.test, "1" : trait_relax.reference}, terms.nucleotideRate("[ACGT]","[ACGT]"));

	// add the relaxation parameter to the test model, and constrain the Q matrix such that the Q matrix of the test model is a function of the Q matrix of the reference model
	parameters.DeclareGlobalWithRanges (trait_relaxation_parameter, 1, terms.trait_relax.k_range[terms.lower_bound], terms.trait_relax.k_range[terms.upper_bound]); 
	model.generic.AddGlobal (trait_relax.test, trait_relaxation_parameter, terms.trait_relax.k); 

	// restrict the first two omegas (expected value is limited to be at most 1 --> use terms.range01)
	for (trait_relax.i = 1; trait_relax.i < trait_relax.rate_classes; trait_relax.i += 1) {
		parameters.SetRange (model.generic.GetGlobalParameter (trait_relax.reference , terms.AddCategory (terms.parameters.omega_ratio,trait_relax.i)), terms.range01);
		parameters.SetRange (model.generic.GetGlobalParameter (trait_relax.test , terms.AddCategory (terms.parameters.omega_ratio,trait_relax.i)), terms.range01);
		// restrict each omega of the test model to be equal to the corresponding omega in the reference model, raised to the power of the relaxation parameter
		parameters.SetConstraint (model.generic.GetGlobalParameter (trait_relax.test , terms.AddCategory (terms.parameters.omega_ratio,trait_relax.i)),
								  model.generic.GetGlobalParameter (trait_relax.reference , terms.AddCategory (terms.parameters.omega_ratio,trait_relax.i)) + "^" + trait_relaxation_parameter,
								  terms.global);
	}

	// restrict the last omega (expected value can be larger than 1 --> use terms.range_gte1) 
	parameters.SetRange (model.generic.GetGlobalParameter (trait_relax.reference , terms.AddCategory (terms.parameters.omega_ratio,trait_relax.rate_classes)), terms.range_gte1);
	parameters.SetRange (model.generic.GetGlobalParameter (trait_relax.test , terms.AddCategory (terms.parameters.omega_ratio,trait_relax.rate_classes)), terms.range_gte1);
	// restrict the last omega of the test model to be equal to the corresponding omega in the reference model, raised to the power of the relaxation parameter
	parameters.SetConstraint (model.generic.GetGlobalParameter (trait_relax.test , terms.AddCategory (terms.parameters.omega_ratio,trait_relax.i)),
							  model.generic.GetGlobalParameter (trait_relax.reference , terms.AddCategory (terms.parameters.omega_ratio,trait_relax.i)) + "^" + trait_relaxation_parameter,
							  terms.global);

							  
	// set a map of the two sub-models of which the RELAX model consists of
	trait_relax.model_object_map = { "trait_relax.reference" : trait_relax.reference,
									 "trait_relax.test" :       trait_relax.test };

}


/**
 * @name trait_relax.set_model_components
 * @usage creates the model instance and initialize its parameters to the values of some preliminary fitting in order to reduce optimization running time in practice
 * @The traitRalex model includes the following components: 
 * 1. parameters of the character model: mu and pie_0
 * 2. parameters of the sequence model. The sequence model consists of two components (test and reference). Each of these has the following parameters, whose values are identical in the two components: 
 * 2a. GTR parameters (5 transition parameters and 3 frequency parameters)
 * 2b. the ka/ks ratio parameters: 5 parameters in total. There are 3 omega's and 2 proportions parameters
 * 2c. The selection intensity parameter: k - this parameter is free to vary only in the test component. In the reference component it is fixed to one. 
 */
function trait_relax.set_model_components() 
{
	
	/* declare the model parameters that are additional to the standard BS-REL model */
	trait_relax.declare_model_parameters();
	
	/* set the codon model as two instances that are extensions of the BS-REL model: test and reference  */
	trait_relax.set_codon_model();
	
	// do preliminary fitting on the null model parameters
	trait_relax.final_partitioned_mg_results = None;
}


/** trait_relax.create_compute_delete - for debugging **/
function trait_relax.create_compute_delete() 
{
	
	iter = 0;
	mu = 1;
	pie_0 = 0.5;
	
	while (iter < 1000) { // fails in the second iteration - regardless if it is local or global function (verified)
	
		fprintf(stdout, "iteration: ", iter, "\n");
		iter += 1;
		
		// need to check what happens in the internal iteration of the grid (iterative process works) that doesn't happen in the external iteration od the grid (fails in the second iteration)
		
		results = trait_relax.compute_joint_likelihood(mu, pie_0, trait_relax.trait_filter_names, trait_relax.trait_trees, &trait_relax.codon_data, trait_relax.trees, 100, "Expected History",
			trait_relax.model_object_map, None, FALSE, maxSimulationsNum, {"return_extended_form": TRUE, utility.getGlobalValue("terms.run_options.retain_lf_object"): FALSE}, &base_model_map);
	}
	
}


/** trait_relax.create_optimize_delete - for debugging **/
function trait_relax.create_optimize_delete() 
{
	
	iter = 0;
	
	while (TRUE) {
		
		fprintf(stdout, "iteration: ", iter, "\n");
		iter += 1;
		
		mu = 1;
		pie_0 = 0.5;
	
		// compute the log likelihood value of the joint likelihood function given the corresponding character grid parameters 
		results = trait_relax.compute_joint_likelihood(mu, pie_0, trait_relax.trait_filter_names, trait_relax.trait_trees, &trait_relax.codon_data, trait_relax.trees, 100, "Exepcted history",
		trait_relax.model_object_map, None, TRUE, maxSimulationsNum, {"return_extended_form": TRUE, utility.getGlobalValue("terms.run_options.retain_lf_object"): FALSE}, &base_model_map);
	}
	
}



/* ____________________ SM AUXILIARY FUNCTIONS ____________________ */


/**
 * @name trait_relax.sample_mutations_given_ancestrals_per_branch
 * @param {Dict} branches_division 		- dictionary to fill in the time duration under 0 and 1 of the branch
 * @param {qDict} branches_history  	- dictionary to fill in the history along the branch
 * @param {String} sonName		   		- name of the node at the bottom of the branch
 * @param {String} fatherName	   		- name of the node at the top of the branch 
 * @param {String} sonState	       		- character state of the son node 
 * @param {String} fatherState	   		- character state of the father node 
 * @param {Float} branch_length	   		- the length of the branch
 * @param {Float} pie_0	  		   		- frequency of state 0
 * @param {Float} mu	  		   		- character substitution rate
 * @param {Integer} maxSimulationsNum	- maximal number of allowed attempt to simulate a branch history 
 * @usage generates a history of a given branch and the characters at its edges
 *	      produces a dictionary that maps to each branch the durations under state 0 and 1	   
 */
lfunction trait_relax.sample_mutations_given_ancestrals_per_branch (branches_division, branches_history, sonName, fatherName, sonState, fatherState, branch_length, pie_0, mu, maxSimulationsNum)
{	
	// reset a vector of transitions for the branch
	branch_history = {};
	branch_history["parent"] = fatherName;
	branch_history["sonState"] = sonState;
	branch_history["parentState"] = fatherState;
	branch_history["history"] = {};
	
	// reset a map to states 0 and 1 the durations under them along the branch
	branch_division = {}; 
	branch_division["parent"] = fatherName;
	branch_division["sonState"] = sonState;
	branch_division["parentState"] = fatherState;
	branch_division["durations"] = {};
	
	for (i=0; i<maxSimulationsNum;i=i+1) {                               // do not exceed the acceptable number of transitions along a branch
	
		branch_transitions_recorder = {};
		tansitionsCounter = 0;
		zero_duration = 0;
		one_duration = 0;
	
		disFromFather = 0;
		curState = fatherState;
		timeTillChange = 0;

		// if the states of the father and son are different, there has so be at least one transition
		// get the time in which it occurs, timeTillChange, according to Nielsen 2001, equations A1,A2
		// we sample timeTillChange conditional on it being smaller than branch_length
		
		if (fatherState != sonState) { 
			freq = pie_0;
			if (curState == 0) {
				freq = 1 - pie_0;
			}
			u = Random(0,1);
			tmp = u * (1 - Exp(-1*mu*freq*branch_length)); 		// no need to multiply rate_matrix[curState][curState] by -1 as its is already negative
			timeTillChange =  Log(1 - tmp) / (-1*mu*freq);      // unlike the demonstration in the article
		}
        
		// as long as the last jump didn't exceed the branch length -> add the current state and time to branch history and draw next state
		while (disFromFather + timeTillChange < branch_length) {
			
			// we now sample a new character for curState:
			// no need to sample the destination state over 2 states. simply choose the complement state to the current one
			if (timeTillChange > 0) // it will be 0 in the first iteration if father==son. Otherwise, not 0.
			{                                                  
				branch_transitions_recorder[tansitionsCounter] = timeTillChange; // add the transitions to the recorder
				tansitionsCounter = tansitionsCounter + 1;
				disFromFather = disFromFather + timeTillChange;
				// set the destination state
				if (curState == 0) {
					zero_duration += timeTillChange;
					curState = 1;
				} else {
					one_duration += timeTillChange;
					curState = 0;
				}
			}
			
			// Now sample the new timeTillChange:
			// the time to change is exponentially distributed with parameter lambda = sum of rates out 
			// of current state. The mean of this exponential distribution is 1/lambda.
			// ^RateMatrix[curState][curState] is negative, therefore, we multiply by -1.
			// sample timeTillChange from a exp(lambdaExpParam) distribution
			// with a Uniform helper. See here: https://en.wikipedia.org/wiki/Inverse_transform_sampling
			freq = pie_0;
			if (curState == 0) {
				freq = 1 - pie_0;
			}
			lambdaExpParam = -1.0 * (-1*mu*freq);         
			uniform_helper = Random(0,1);             				  		   // sample timeTillChange from a exp(lambdaExpParam) distribution
			timeTillChange = - (1.0 / lambdaExpParam) * Log(uniform_helper);   // with a Uniform helper. See here: https://en.wikipedia.org/wiki/Inverse_transform_sampling
		}
		
		// if the destination state of the last transition (which ends in the son) is the same state as the son's -> accept the history
		if (curState == sonState) {                                           // if the current state is the son's state -> the simulation succeeded -> record the last jump
			// record all branch history
			branch_transitions_recorder[tansitionsCounter] = branch_length - disFromFather;
			(branch_history["history"]) = branch_transitions_recorder;
			if (curState == 0) {
				zero_duration = zero_duration + (branch_length-disFromFather);
			} else {
				one_duration = one_duration + (branch_length-disFromFather);
			}
			tansitionsCounter = tansitionsCounter + 1;
			branch_history["transitionsNum"] = tansitionsCounter-1;
			(^branches_history)[sonName] = branch_history;
			(branch_division["durations"])[0] = zero_duration;
			(branch_division["durations"])[1] = one_duration;
			(^branches_division)[sonName] = branch_division;
			return 0;
			
		}
	}

	// if all simulations failed -> exit
	fprintf(stdout, "could not produce simulations with father = ", fatherName, " son = ", sonName, " branch length = ", branch_length, "\n\n");
	branch_history["transitionsNum"] = 0;
	(^branches_history)[sonName] = {};
	exit();
}



/**
 * @name trait_relax.generate_history
 * @param {Integer} history_num 		- number of the history that should be created (goes from 0 to num_of_histories)
 * @param {Dict} trees_duration  		- dictionary to fill in the durations under state 0 and state 1 for each branch in the tree
 * @param {Integer} lfId		   		- pointer of the character likelihood function
 * @param {Dict} tree		   			- dictionary holding the original user given tree
 * @param {Float} pie_0	  		   		- frequency of state 0
 * @param {Float} mu	  		   		- character substitution rate
 * @param {Integer} maxSimulationsNum	- maximal number of allowed attempt to simulate a branch history 
 * @return {Dict} historyTree			- a dictionary holding the sampled history in a tree format
 * @usage generates a history of a given branch and the characters at its edges
 *	      produces a dictionary that maps to each branch the durations under state 0 and 1	   
 */
lfunction trait_relax.generate_history (history_num, trees_duration, lfId, tree, pie_0, mu, maxSimulationsNum, run_options)
{
	// get the indicts of the tree branches from the tree instance
	branchNames = trees.BranchNames(tree); 		 		// branch names are ordered in post-traversal order
	ancestors = ancestral.build (lfId, 0, None); 		// last argument was initially None
	nodesNum = (ancestors["DIMENSIONS"])["BRANCHES"];			
	branchesHistory = {};                          		// histories along branches
	branches_division = {};						   		// overall 0 and 1 durations along branches
		
	for (sonID=1; sonID<nodesNum; sonID=sonID+1) { // disregard the last node, which is the root, as it is not a fatherless node
		
		// get the Id and state of the son node and its father, and the length of the branch connecting them
		sonName = ((ancestors["TREE_AVL"])[sonID])["Name"];
		fatherID = ((ancestors["TREE_AVL"])[sonID])["Parent"];
		fatherName = ((ancestors["TREE_AVL"])[fatherID])["Name"];
		sonState = (ancestors["MATRIX"])[sonID-1];
		fatherState = (ancestors["MATRIX"])[fatherID-1];
		branchLength = (tree[utility.getGlobalValue("terms.branch_length")])[branchNames[sonID-1]];
		
		// generate the history along the branch connecting the node to its father
		if (branchLength > 0) {
			// note: in the dictionary of history per branch, the last transition represent the length of the sub-branch connected to the son
			trait_relax.sample_mutations_given_ancestrals_per_branch(&branches_division, &branchesHistory, sonName, fatherName, sonState, fatherState, branchLength, pie_0, mu, maxSimulationsNum);
		}
	}
	
	// convert the SM info to tree with branch classification
	// represent the history as a tree with internal nodes, that have a single child
	// the trait states per branch in the new tree will dictate the label of the branch
	historyTree = trait_relax.generate_tree(tree, &branchesHistory, run_options, history_num);
	(^ trees_duration)[history_num] = branches_division; 
	
	// delete the ancestral character data
	DeleteObject(ancestors);
	
	return historyTree;
}


/**
 * @name trait_relax.average_histories
 * @param {Dict} tree		   					- dictionary holding the original user given tree
 * @param {Dict} histories_duration_analysis  	- dictionary that holds the durations under state 0 and state 1 for each branch in the tree
 * @return {Dict} expected_history			    - a dictionary holding the expected (averaged) history in a tree format
 * @usage generates a history of a given branch and the characters at its edges
 *	      produces a dictionary that maps to each branch the durations under state 0 and 1	   
 */
lfunction trait_relax.average_histories(tree, histories_duration_analysis, run_options) { 

	branch_names = trees.BranchNames(tree); // get the core branch names from the original tree
	branch_lengths = tree[utility.getGlobalValue("terms.branch_length")];
	histories_num = Columns(utility.Keys(^ histories_duration_analysis));
	
	// if documentation of history is required -> add "average_" prefix to the file
	if (None != run_options) {
		if (run_options["write_histories"] == TRUE) {
			run_options["average"] = TRUE;
		}
	}
	
	// initialize a dictionary that will hold the expected history
	expected_history_info = {};								  // the expected history of the tree branches
	expected_durations = {};							  	  // the expected durations of the tree branches
	for (b=0; b<Columns(branch_names)-1; b += 1) {		  	  // exclude the root
		sonName = branch_names[b];
		// fill in the expected history dictionary with the branch properties
		expected_branch_history = {};
		expected_branch_history["parent"] = (((^ histories_duration_analysis)[0])[sonName])["parent"];
		expected_branch_history["sonState"] = (((^ histories_duration_analysis)[0])[sonName])["sonState"];
		expected_branch_history["history"] ={}; // need to adjust according to son and parent state? no! need to set parent state to complement the son state?
		expected_history_info[sonName] = expected_branch_history;
		// set the initial durations to 0
		expected_durations[sonName] = {0:0, 1:0};
	}
	
	// set the ancestral states as the common states across all histories
	internal_to_statesDistribution = {};
	node_to_type = tree[utility.getGlobalValue("terms.trees.partitioned")];
	for (b=0; b<Columns(branch_names)-1; b += 1) {		
		if (node_to_type[branch_names[b]] % "internal") {
			internal_to_statesDistribution[branch_names[b]] = {0:0, 1:0};
		}
	}
	internal_to_statesDistribution["Node0"] = {0:0, 1:0};	// set an entry for the root as well

	// for each history, document the state of each internal node
	for (h=0; h<histories_num; h+=1) {
		root_set = 0;
		for (b=0; b<Columns(branch_names)-1; b += 1) {
			if (node_to_type[branch_names[b]] % "internal") {
				if ((((^ histories_duration_analysis)[h])[branch_names[b]])["sonState"] == 0) {
					(internal_to_statesDistribution[branch_names[b]])[0] = (internal_to_statesDistribution[branch_names[b]])[0] + 1;
				} else {
					(internal_to_statesDistribution[branch_names[b]])[1] = (internal_to_statesDistribution[branch_names[b]])[1] + 1;
				}
			}
			if ((((^ histories_duration_analysis)[0])[sonName])["parent"] % "Node0" && root_set == 0) {
				root_state = (((^ histories_duration_analysis)[0])[sonName])["parentState"];
				(internal_to_statesDistribution["Node0"])[root_state] = (internal_to_statesDistribution["Node0"])[root_state] + 1;
				root_set = 1;
			}
		}
	}

	// set the most common state in each internal node to be the its state in the expected history
	internal_to_state = {};
	internals = utility.Keys(internal_to_statesDistribution);
	for (k=0; k<Columns(internals); k += 1) {
		zero_count = (internal_to_statesDistribution[internals[k]])[0];
		one_count = (internal_to_statesDistribution[internals[k]])[1];
		if (zero_count > one_count) {
			internal_to_state[internals[k]] = 0;
		} else {
			internal_to_state[internals[k]] = 1;
		}
	}
	
	for (b=0; b<Columns(branch_names)-1; b += 1) {		  // exclude the root
		sonName = branch_names[b];
		(expected_history_info[sonName])["parentState"] = internal_to_state[(expected_history_info[sonName])["parent"]];
	}
	
	// set the histories along branches (up to 2 transitions, where one transition is of length 0)
	for (h=0; h<histories_num; h += 1) {
		durations_info = (^ histories_duration_analysis)[h];
		for (b=0; b<Columns(branch_names)-1; b += 1) {  // exclude the root
			(expected_durations[branch_names[b]])[0] = (expected_durations[branch_names[b]])[0] + ((durations_info[branch_names[b]])["durations"])[0];
			(expected_durations[branch_names[b]])[1] = (expected_durations[branch_names[b]])[1] + ((durations_info[branch_names[b]])["durations"])[1];	
		}
	}
	
	// exchange the sum with average
	for (b=0; b<Columns(branch_names)-1; b += 1) {		// exclude the root
		(expected_durations[branch_names[b]])[0] = (expected_durations[branch_names[b]])[0] / histories_num;
		(expected_durations[branch_names[b]])[1] = (expected_durations[branch_names[b]])[1] / histories_num;
	}
	
	// convert the expected durations dictionary to a history
	// no need in additional transition in case of parent_state==son_state because the expected history can have two consecutive branches with the same label
	for (b=0; b<Columns(branch_names)-1; b += 1) {		// exclude the root
		son_state = (expected_history_info[branch_names[b]])["sonState"];
		if ((expected_durations[branch_names[b]])[0] == 0 || (expected_durations[branch_names[b]])[1] == 0) {
			(expected_history_info[branch_names[b]])["transitionsNum"] = 0;
			((expected_history_info[branch_names[b]])["history"])[0] = (expected_durations[branch_names[b]])[0] + (expected_durations[branch_names[b]])[1]; // one term is 0 and the other is the branch length. Thus, it suffices to add them together
		} else { // history is read by trait_relax.generate_tree top to bottom, so the last transition denotes the duration between the son and the first internal node
			(expected_history_info[branch_names[b]])["transitionsNum"] = 1;
			((expected_history_info[branch_names[b]])["history"])[0] = (expected_durations[branch_names[b]])[1-son_state];
			((expected_history_info[branch_names[b]])["history"])[1] = (expected_durations[branch_names[b]])[son_state];
		}
	}
	expected_history = trait_relax.generate_tree(tree, &expected_history_info, run_options, 1);
	return expected_history;
}


/**
 * @name trait_relax.generate_tree
 * @param {Dict} tree		   	   - dictionary holding the original user given tree
 * @param {Dict} branches_history  - dictionary that holds the sampled character history
 * @return {Dict} history	       - a dictionary holding the history in a tree format
 * @usage convert the history in its dictionary form to a tree based on the original tree   
 */
lfunction trait_relax.generate_tree (tree, branches_history, run_options, history_num)
{
	branch_names = trees.BranchNames(tree); // get the core branch names from the original tree
	history_tree = {};
	internalsCounter = 1;
	tree_str = tree[utility.getGlobalValue("terms.trees.newick_with_lengths")]; 
	Topology historyTree = tree_str; // use utility.getGlobalValue("terms.trees.newick")] if you don't need branch lengths
	label_to_branch = {};
	
	// when not specifying a name, an internal node with a single child is expected to be generated
	for (i=0; i<Columns(branch_names)-1; i=i+1) {
		node = branch_names[i];
		branch_history = (^branches_history)[node];
		formerLabel = 1;
		// set the label of the node
		if (branch_history["sonState"] == 1) {
			label_to_branch[node] = 1;
		} else {
			label_to_branch[node] = 0;
			formerLabel = 0;
		}
		if (branch_history["transitionsNum"] > 0) {
			for (tc=branch_history["transitionsNum"]-1; tc>=0; tc=tc-1) { // in the branch history, the transitions are ordered from the parent node to the child node
																	      // therefore, they should be added to the tree in descending order
				parent = "I" + internalsCounter;
				if (formerLabel == 0) {
					label_to_branch[parent] = 1;
					formerLabel = 1;
				} else {
					label_to_branch[parent] = 0;
					formerLabel = 0;
				}
				historyTree + {"WHERE" : node , "PARENT" : parent, "PARENT_LENGTH" : (branch_history["history"])[tc]};
				node = parent;
				internalsCounter = internalsCounter + 1;
			}
		// in case of 0 transitions -> we don't add a new node. we only change the length of the branch connecting the last node with the core parent
		} else {  
			if (branch_history["parentState"] == 1) {
				label_to_branch[branch_history["parent"]] = 1;
			} else {
				label_to_branch[branch_history["parent"]] = 0;
			}
		}
	}
	// convert hostoryTree back to string historyTreeStr
	charHistoryTreeStr = Format(historyTree,1,1);
	Tree blankHistoryTree = charHistoryTreeStr;
	historyTreeStr = Format(blankHistoryTree,1,1);
	
	// edit the branch lengths of the code nodes
	branch_lengths = tree[utility.getGlobalValue("terms.branch_length")];
	for (i=0; i<Columns(branch_names); i=i+1) {
		node = branch_names[i];
		src_bl = branch_lengths[node];
		node_history = (^branches_history)[node];
		transitions_num = node_history["transitionsNum"];
		dst_bl = (node_history["history"])[transitions_num];
		src_expr = node + ":" + src_bl;
		dst_expr = node + ":" + dst_bl;
		historyTreeStr = custom_functions.strReplace(historyTreeStr, src_expr, dst_expr);		
	}

	// add labels to each node using string conversion and set the model map
	model_map = {};
	new_branches_names = utility.Keys(label_to_branch); // bug here due to the order of labeling
	for (i=0; i<Columns(new_branches_names); i=i+1) {
		node = new_branches_names[i];
		label = label_to_branch[node];
		model_map[node] = "";
		src_expr = node + ":";
		if (label) {
			dst_expr = node + "{Test}:";
			historyTreeStr = custom_functions.strReplace(historyTreeStr, src_expr, dst_expr);
			model_map[node] = "Test";
		}
	}
	
	// debug - write the last character histories to files
	if (None != run_options) {
		if (run_options["write_histories"] == TRUE) {
			history_output_path = run_options["histories_dir"] + "history_" + history_num;
			if (run_options["average"] == TRUE) {
				history_output_path =  run_options["histories_dir"] + "average_history";
			}
			historyTreeStr = custom_functions.strReplace(historyTreeStr, "{Test}", "{FG}");
			fprintf(history_output_path, historyTreeStr, ";");
		}
	}

	history = trees.ExtractTreeInfo(historyTreeStr);
	return history;
}


/* ________________ TRAIT RELAX AUXILIARY FUNCTIONS __________________ */


/**
 * @name trait_relax.generate_character_likelihood_function
 * @param {Float} mu    						- character substitution rate
 * @param {Float} pie_0    						- frequency of state 0 in the character model
 * @param {String} lf_prefix					- string that holds the prefix of the likelihood function in the global namespace
 * @param {Dict} filter_name  					- dictionary holding a pointer to the character data
 * @param {Dict} tree							- dictionary holding the tree provided by the user
 * @usage generates a likelihood function of the character model given initial parameters
 */
lfunction trait_relax.generate_character_likelihood_function(mu, pie_0, lf_prefix, data_filter, tree) {
	
	character_initial_parameters = {
	 "global":{
	   "Substitution rate from character 0 to character 1":
	   {utility.getGlobalValue("terms.fit.MLE"): mu}
	  },
	 utility.getGlobalValue("terms.efv_estimate"):{
	   utility.getGlobalValue("terms.model"):{
		{pie_0}
		{1-pie_0}
		}
	  }
	};
	
	trait_likelihood_function = custom_functions.CreateLFObject_FixedBLs (lf_prefix, data_filter, tree, "models.binaryML.ModelDescription", character_initial_parameters, {utility.getGlobalValue("terms.run_options.retain_lf_object"): TRUE, "set_freq": TRUE, "return_lf": TRUE});
	
	return trait_likelihood_function
}


/* _____________________________________________________________________*/
/* Generates a rule map that maps each branch in he history to is		*/ 
/* corresponding label 													*/
lfunction trait_relax.generate_history_rules_map(history) {
	test_dict = {};
	reference_dict = {};
	model_map_guide = history[utility.getGlobalValue("terms.trees.model_map")];
	branch_names = utility.Keys(model_map_guide);
	for (i=0; i<Columns(branch_names); i=i+1) {
		branch = branch_names[i];
		if (model_map_guide[branch] % "Test") {
			test_dict[branch] = "Test";
		} else {
			reference_dict[branch] = "Reference";
		}
	}
	model_map = { "trait_relax.test" : test_dict,
				  "trait_relax.reference" : reference_dict };	
				  
	return model_map
}


/**
 * @name trait_relax.wrapper_compute_joint_likelihood
 * @param {Dict} args_dict - dictonary holding the arguments to call trait_relax.compute_joint_likelihood with
 * @return {Float} log likelihood - the output from the call to trait_relax.compute_joint_likelihood with the parsed arguments
 */
lfunction trait_relax.wrapper_compute_joint_likelihood(args_dict) {
	
	mu = (^args_dict)[0];
	pie_0 = (^args_dict)[1];
	trait_filter_name = (^args_dict)[2];
	trait_tree = (^args_dict)[3];
	sequence_dataset = (^args_dict)[4];
	sequence_trees = (^args_dict)[5];
	histories_num = (^args_dict)[6];
	approximation_method = (^args_dict)[7];
	model_object_map = (^args_dict)[8];
	relax_initial_parameters = (^args_dict)[9];
	optimize_sequence_parameters = (^args_dict)[10];
	max_simulations_num = (^args_dict)[11];
	run_options = (^args_dict)[12];
	base_model_map = (^args_dict)[13];
	
	log_likelihood = trait_relax.compute_joint_likelihood(mu, pie_0, trait_filter_name, trait_tree, sequence_dataset, sequence_trees, histories_num, approximation_method, model_object_map, relax_initial_parameters, optimize_sequence_parameters, max_simulations_num, run_options, base_model_map);
	
	return log_likelihood;
} 



/**
 * @name trait_relax.compute_joint_likelihood
 * @param {Float} mu    						- character substitution rate
 * @param {Float} pie_0    						- frequency of state 0 in the character model
 * @param {Dict} filter_name  					- dictionary holding a pointer to the character data
 * @param {Dict} tree							- dictionary holding the tree provided by the user
 * @param {Integer} histories_num				- number of histories to consider in the character history approximation
 * @param {String} approximation_method 		- the approximation method that would be applied in the likelihood computation (integrating over all histories or averaging all to a single expected history)
 * @param {Dict} model_object_map  				- dictionary that holds the test and reference models that should be applied on the history's branches
 * @param {Dict} relax_initial_parameters		- dictionary the holds the initial sequence model parameters
 * @param {Bool} optimize_sequence_parameters	- TRUE if you wish to optimize the likelihood function over the sequence model parameters or not
 * @param {Integer} max_simulations_num			- the maximal number of simulations per branch before simulation failure is declared
 * @return {Float} log_likelihood	       		- the computed log likelihood
 * @usage simulates character histories, generates the joint likelihood function over these histories / the average history and computed / optimizes the likelihood function. must be a global function in order to have access to the joint likelihood function, which is generated in global namespace (solution based on estimators.bf with 'lfid' resulted in log likelihood -inf
 */
lfunction trait_relax.compute_joint_likelihood(mu, pie_0, trait_filter_name, trait_tree, sequence_dataset, sequence_trees, histories_num, approximation_method, model_object_map, relax_initial_parameters, optimize_sequence_parameters, max_simulations_num, run_options, base_model_map) {
	
	// create the codono model again, since it is deleted every time a likelihood function is deleted - it really slows down my code:( There must be another way... Must be able to work with the same likelihood function... 
	utility.ExecuteInGlobalNamespace("trait_relax.set_codon_model();");
	
	joint_sequence_filter_names = {}; 
	joint_sequence_trees = {}; 

	// initialize the input for the construction of the joint likelihood function
	lf_formula = "'Log((1-trait_relax.p)*SITE_LIKELIHOOD[0]+trait_relax.p*((";
	custom_functions.dup_data_filter("data_filter_0", sequence_dataset); // segmentation falut from here... did I delete sequence_dataset? possibly, when deleting the likelihood function:(
	joint_sequence_filter_names[0] = "data_filter_0";
	joint_sequence_trees[0] = sequence_trees[0];
	model_maps = {};
	model_maps[0] = ^ base_model_map;
	
	// create the likelihood function instance for the SM procedure
	trait_likelihood_function = trait_relax.generate_character_likelihood_function(mu, pie_0, "trait_relax.trait", trait_filter_name, trait_tree);
	trait_relax.char_log_likelihood = estimators.ComputeLF(trait_likelihood_function);

	// initialize 0 and 1 durations along the tree per generated history, and the histories themselves
	histories_durations_dict = {};
	histories_dict = {};
		
	// run the TraitRELAX pipeline on histories_num character (trait) histories 
	for (h = 0; h < histories_num; h = h+1) {  
		 
		// generate a stochastic mapping and put it in the path labelledtree_path	
		history = trait_relax.generate_history(h, &histories_durations_dict, trait_likelihood_function, trait_tree["0"], pie_0, mu, max_simulations_num, run_options);
		
		// if the user chose to fit the alternative model to multiple histories -> add the history to the likelihood function
		if (approximation_method == "Multiple Histories") {
			
			// apply the alternative model on the history
			histories_dict[h] = {"id": history_id, "tree": history};

			// add the history to the likelihood function formula
			if (h == 0) {
				lf_formula = lf_formula + "SITE_LIKELIHOOD[" + (h+1) + "]";
			} else {
				lf_formula = lf_formula + "+SITE_LIKELIHOOD[" + (h+1) + "]";
			}
			
			// add he history to the sequence dataset 
			joint_sequence_filter_names [h+1] = "data_filter_" + h+1;
			custom_functions.dup_data_filter(joint_sequence_filter_names [h+1], sequence_dataset);
			joint_sequence_trees[h+1] = history;
			model_maps[h+1] = trait_relax.generate_history_rules_map(history);
		}
	}
	
	// delete the trait likelihood function
	custom_functions.DeleteLikelihoodFunction_FixedBLs(trait_likelihood_function);
	
	// user chose to compute the likelihood over all the character histories --> include all of them in the likelihood function formula
	if (approximation_method == "Multiple Histories") {
		lf_formula = lf_formula + ")/" + histories_num + "))'";
	
	// user chose to average the histories into a single expected history    --> average the histories and use the expected history in the joint likelihood function formula
	} else { 
		expected_history_info = trait_relax.average_histories(sequence_trees[0], &histories_durations_dict, run_options); // compute the expected history
		custom_functions.dup_data_filter("data_filter_1", sequence_dataset) ;
		joint_sequence_filter_names[1] = "data_filter_1";
		joint_sequence_trees[1] = expected_history_info;
		model_maps[1] = trait_relax.generate_history_rules_map(expected_history_info);
		histories_dict[1] = {"id": expected_history_info_id, "tree": expected_history_info};
		lf_formula = lf_formula + "SITE_LIKELIHOOD[1])))'";
	}

	retain_lf = FALSE;
	if (run_options[utility.getGlobalValue("terms.run_options.retain_lf_object")]) {
		retain_lf = TRUE;
	}
	
	// generate and optimize the joint likelihood function (i.e, while optimizing the sequence model parameters)
	if (optimize_sequence_parameters) {
		trait_relax.alternative_model.initial_parameters = custom_functions.FitLF_FixedBLs(joint_sequence_filter_names, joint_sequence_trees, lf_formula, model_maps, relax_initial_parameters, model_object_map, {utility.getGlobalValue("terms.run_options.retain_lf_object"): retain_lf, terms.run_options.proportional_branch_length_scaler: trait_relax.scaler, "custom_lf_formula": TRUE});
	}
	
	// generate and compute the log likelihood of the joint likelihood function (i.e, while fixing the sequence model parameters)
	else {	
		trait_relax.alternative_model.initial_parameters = custom_functions.ComputeLF_FixedBLs(joint_sequence_filter_names, joint_sequence_trees, lf_formula, model_maps, relax_initial_parameters, model_object_map, {utility.getGlobalValue("terms.run_options.retain_lf_object"): retain_lf, terms.run_options.proportional_branch_length_scaler: trait_relax.scaler, "custom_lf_formula": TRUE});
	}
	
	trait_relax.seq_log_likelihood = trait_relax.alternative_model.initial_parameters[utility.getGlobalValue("terms.fit.log_likelihood")];

	trait_relax.joint_log_likelihood = trait_relax.char_log_likelihood + trait_relax.seq_log_likelihood;
	
	if (run_options["return_extended_form"]) { 
		trait_relax.alternative_model.initial_parameters [utility.getGlobalValue("terms.fit.log_likelihood")] = trait_relax.joint_log_likelihood;
		return &trait_relax.alternative_model.initial_parameters;
	} else {
		return trait_relax.joint_log_likelihood;
	}
}

		
/**
 * @name trait_relax.set_mp_partition
 * @param {Dict} partition_info  - dictionary that holds the user sequence input (codon alignment and tree)
 * @return {Dict} return_set	 - mapping of the tree branches to the test and reference models
 * @usage classifies all the branches except for the first one to reference (used as some meaningless partition for the null model optimization) 
 */		
lfunction trait_relax.set_mp_partition(partition_info, trait_data) {
	
	state_to_label = {"0": "Reference" , "1": "Test"};
	return_set = custom_functions.getMPTreePartition((partition_info["0"])[utility.getGlobalValue("terms.data.tree")], state_to_label, trait_data); // move to to he null section and change the branches selector accordingly
	return {"0": return_set }; 
	
}


/**
 * @name trait_relax.get_base_model_map
 * @param {Dict} tree  		   - dictionary that holds the user provided tree 
 * @usage creates a base model map that maps all branches into the reference category
 */
lfunction trait_relax.get_base_model_map(tree) {
	branch_names = trees.BranchNames(tree);
	model_map = {};
	reference_dict = {};
	for (i=0; i<Columns(branch_names); i=i+1) {
		branch = branch_names[i];
		reference_dict[branch] = "Reference";
	}
	model_map = { "trait_relax.test" : {},
				  "trait_relax.reference" : reference_dict };
	return model_map;
}	


/**
 * @name trait_relax.compute_bayes_measurements
 * @param {Matrix} sitewise_blockwise_likelihoods - matrix holding the site likelihood per codon site per history
 * @param {Integer} histories_num  		  		  - number of histories
 * @param {Float} dependent_propertion  		  - estimated proportion of sites whose evolution depend on the trait evolution
 * @return {Dict} site_bayes_measurements         - a dictionary that holds the bayes factor per site and and empirical bayes posterior per site vectors
 * @usage computed the byes factor and empirical bayes posterior measurement for each codon site
 */	
lfunction trait_relax.compute_bayes_measurements(sitewise_blockwise_likelihoods, histories_num, dependent_proportion) {
	sites_number = Columns(sitewise_blockwise_likelihoods) / (histories_num+1);   // the additional block is the independent case block (i.e, tree)
	sites_bayes_factors = {};
	sites_bayes_potesriors = {};
	for (i=0; i<sites_number; i+=1) {
		independent_case_likelihood = sitewise_blockwise_likelihoods[i*(histories_num+1)];
		dependent_case_likelihood = 0;
		for (j=1; j<=histories_num; j+=1) {
			dependent_case_likelihood = dependent_case_likelihood + sitewise_blockwise_likelihoods[i*histories_num+j]; // sitewise_blockwise_likelihoods is ordered by sites and then by blocks (=trees)
																													   // for example: say I have two sites and 3 trees: position 0 corresponds to (site 1, tree 1) 
																													   //												position 1 corresponds to (site 1, tree 2) 
																													   //												position 3 corresponds to (site 2, tree 1) 
		}
		dependent_case_likelihood = dependent_case_likelihood / histories_num;
		sites_bayes_factors[i] = dependent_case_likelihood / independent_case_likelihood; // compute the bayes factor of site i
		sites_bayes_potesriors[i] = (dependent_case_likelihood * dependent_proportion) / (dependent_case_likelihood * dependent_proportion + independent_case_likelihood * (1-dependent_proportion)); // bug here
	}
	site_bayes_measurements = {}; // first line is bayes factors, second line is empirical bayes posteriors
	site_bayes_measurements["bayes factor per site"] = sites_bayes_factors;
	site_bayes_measurements["empirical bayes posterior per site"] = sites_bayes_potesriors;
	return site_bayes_measurements;	
}
		
		
		
/* ______________ RELAX AUXILIARY FUNCTIONS (DEEP COPY) _______________ */

lfunction relax.extract.k(branch_info) {
    return (branch_info[utility.getGlobalValue("terms.trait_relax.k")])[utility.getGlobalValue("terms.fit.MLE")];
}

//------------------------------------------------------------------------------

lfunction relax.set.k (tree_name, node_name, model_description) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.trait_relax.k"), "String")) {
        k = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.trait_relax.k")];
        parameters.SetValue (tree_name + "." + node_name + "." + k, 1);
        parameters.SetRange (tree_name + "." + node_name + "." + k, utility.getGlobalValue ("terms.trait_relax.k_range"));
    }
    return tree_name + "." + node_name + "." + k;
}

//------------------------------------------------------------------------------

lfunction relax.init.k (lfId, components, data_filter, tree, model_map, initial_values, trait_relax.model_object_map) {
    parameter_set = estimators.TraverseLocalParameters (lfId, trait_relax.model_object_map, "relax.set.k");
    parameters.SetConstraint (model.generic.GetGlobalParameter (utility.getGlobalValue("trait_relax.ge") , terms.AddCategory (utility.getGlobalValue("terms.parameters.omega_ratio"),2)), utility.getGlobalValue("terms.parameters.one"), utility.getGlobalValue("terms.global"));
    /*parameters.SetConstraint (model.generic.GetGlobalParameter (utility.getGlobalValue("trait_relax.ge") , terms.AddCategory (utility.getGlobalValue("terms.parameters.omega_ratio"),utility.getGlobalValue ("trait_relax.rate_classes"))),
                             "1/(" +
                                Join ("*", utility.Map (
                                    utility.Range (utility.getGlobalValue ("trait_relax.rate_classes") - 1, 1, 1),
                                    "_value_",
                                    'model.generic.GetGlobalParameter (trait_relax.ge , terms.AddCategory (terms.parameters.omega_ratio,_value_))'
                                    ))
                             + ")",
                            "global");*/

    return 0;
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL.ModelDescription (type, code, components) {
    model = models.codon.BS_REL.ModelDescription(utility.getGlobalValue ('terms.global'), code, components);
    model [utility.getGlobalValue("terms.model.defineQ")] = "relax.BS_REL._DefineQ ";
    return model;
}

//------------------------------------------------------------------------------

lfunction relax.DistributionGuess (mean) {
    guess = {{0.05,0.7}{0.25,0.2}{10,0.1}};

    norm = + guess[-1][1];
    guess_mean = 1/(+(guess [-1][0] $ guess [-1][1]))/norm;
    return guess["_MATRIX_ELEMENT_VALUE_*(guess_mean*(_MATRIX_ELEMENT_COLUMN_==0)+(_MATRIX_ELEMENT_COLUMN_==1)*(1/norm))"];
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL._GenerateRate (fromChar, toChar, namespace, model_type, _tt, alpha, alpha_term, beta, beta_term, omega, omega_term) {

    p = {};
    diff = models.codon.diff(fromChar, toChar);

    if (None != diff) {
        p[model_type] = {};
        p[utility.getGlobalValue("terms.global")] = {};

        if (diff[utility.getGlobalValue("terms.diff.from")] > diff[utility.getGlobalValue("terms.diff.to")]) {
            nuc_rate = "theta_" + diff[utility.getGlobalValue("terms.diff.to")] + diff[utility.getGlobalValue("terms.diff.from")];
        } else {
            nuc_rate = "theta_" + diff[utility.getGlobalValue("terms.diff.from")] + diff[utility.getGlobalValue("terms.diff.to")];
        }
        nuc_rate = parameters.ApplyNameSpace(nuc_rate, namespace);
        (p[utility.getGlobalValue("terms.global")])[terms.nucleotideRate(diff[utility.getGlobalValue("terms.diff.from")], diff[utility.getGlobalValue("terms.diff.to")])] = nuc_rate;

        if (_tt[fromChar] != _tt[toChar]) {
            if (model_type == utility.getGlobalValue("terms.global")) {
                aa_rate = parameters.ApplyNameSpace(omega, namespace);
                (p[model_type])[omega_term] = aa_rate;
                utility.EnsureKey (p, utility.getGlobalValue("terms.local"));
                 (p[utility.getGlobalValue("terms.local")])[utility.getGlobalValue ("terms.trait_relax.k")] = "k";
                 aa_rate += "^k";
            } else {
                aa_rate = beta;
                (p[model_type])[beta_term] = aa_rate;
            }
            p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + aa_rate;
        } else {
            if (model_type == utility.getGlobalValue("terms.local")) {
                (p[model_type])[alpha_term] = alpha;
                p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate + "*" + alpha;
            } else {
                p[utility.getGlobalValue("terms.model.rate_entry")] = nuc_rate;
            }
        }
    }

    return p;
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL._DefineQ (bs_rel, namespace) {
    rate_matrices = {};

    bs_rel [utility.getGlobalValue("terms.model.q_ij")] = &rate_generator;
    bs_rel [utility.getGlobalValue("terms.mixture.mixture_components")] = {};

    _aux = parameters.GenerateSequentialNames (namespace + ".bsrel_mixture_aux", bs_rel[utility.getGlobalValue("terms.model.components")] - 1, "_");
    _wts = parameters.helper.stick_breaking (_aux, None);
    mixture = {};

    for (component = 1; component <= bs_rel[utility.getGlobalValue("terms.model.components")]; component += 1) {
       key = "component_" + component;
       ExecuteCommands ("
        function rate_generator (fromChar, toChar, namespace, model_type, _tt) {
           return relax.BS_REL._GenerateRate (fromChar, toChar, namespace, model_type, _tt,
                'alpha', utility.getGlobalValue('terms.parameters.synonymous_rate'),
                'beta_`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.nonsynonymous_rate'), component),
                'omega`component`', terms.AddCategory (utility.getGlobalValue('terms.parameters.omega_ratio'), component));
            }"
       );

       if ( component < bs_rel[utility.getGlobalValue("terms.model.components")]) {
            model.generic.AddGlobal ( bs_rel, _aux[component-1], terms.AddCategory (utility.getGlobalValue("terms.mixture.mixture_aux_weight"), component ));
            parameters.DeclareGlobalWithRanges (_aux[component-1], 0.5, 0, 1);
       }
       models.codon.generic.DefineQMatrix(bs_rel, namespace);
       rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rateMatrix")];
       (bs_rel [^'terms.mixture.mixture_components'])[key] = _wts [component-1];
    }


    bs_rel[utility.getGlobalValue("terms.model.rateMatrix")] = rate_matrices;
    parameters.SetConstraint(((bs_rel[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.nucleotideRate("A", "G")], "1", "");
    return bs_rel;
}





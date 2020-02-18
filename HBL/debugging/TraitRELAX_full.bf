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
LoadFunctionLibrary("libv3/models/binary.bf");
LoadFunctionLibrary("modules/custom_functions.bf");


/* ________________________ ENVIRONMENT VARIABLES ________________________ */

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);	// defined based on RELAX.bf
utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);	// defined based on RELAX.bf
utility.SetEnvVariable ("REPLACE_TREE_STRUCTURE", TRUE); 	// assures initialization of the tree every time it is re-loaded - required for loading histories repeatedly, with the same name
utility.SetEnvVariable("ACCEPT_ROOTED_TREES", TRUE);		// allows to preserve rooted trees in the generation of a HyPhy tree instance

/* ________________________ HARDCODED PARAMETERS ________________________ */

// character grid settings parameters
// KEREN 7.12.17: increased both from the original values 0.1 to reduce running time when debugging
mu_interval = 0.5;
pie_0_interval = 0.2;

// character histories sampling properties
maxSimulationsNum = 100000;

/* ________________________________ MAIN ________________________________ */

/* display model to the user */
trait_relax.display_model();


/* initiate script parameters */
trait_relax.init();

// memory capacity on dummy data is 96260

/* process user input */
trait_relax.process_input();

// memory capacity on dummy data is 96260

/* create the sequence model components (test and reference) and perform preliminary fitting for initialization of the model parameters */
trait_relax.set_model_components();

// memory capacity on dummy data is 173m

/* optimize the null model */
trait_relax.optimize_null_model();

// memory capacity on dummy data is 1932m - unexpected spike in trait_relax.optimize_null_model()

/* optimize the alternative model */
trait_relax.optimize_alternative_model();


/* ___________________ MAIN AUXILIARY FUNCTIONS ___________________ */


/**
 * @name trait_relax.display_model
 * @usage presents model into into the terminal
 */
function trait_relax.display_model() {
	
	trait_relax.analysis_description = {
								   terms.io.info : "TraitRELAX (A test of association between a binary phenotypic trait and the selection pattern (intensification or relaxation) that operates on the sequence data",
								   terms.io.version : "1.0",
								   terms.io.reference : "TraitRELAX: Detecting association between phenotypic traits and changes in selective pressure (2018)... hopefully?",
								   terms.io.authors : "Keren Halabi, Eli Levi Karin, Sergei L Kosakovsky Pond, Itay Mayrose, Tel Aviv University and Temple iGEM",
								   terms.io.contact : "halabikeren@mail.tau.ac.il",
								   terms.io.requirements : "trait character alignment, codon alignment and a phylogenetic tree"
								  };

	io.DisplayAnalysisBanner ( trait_relax.analysis_description );
	
}

		
/**
 * @name trait_relax.init
 * @usage sets the initial parameters for the pipeline
 */
function trait_relax.init() {
		
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
function trait_relax.process_input() {
	
	/* accept input paths from the user */
	trait_relax.trait_data_path = io.PromptUserForString("\n>Provide full path to the trait character alignment");
	trait_relax.tree_path =  io.PromptUserForString("\n>Provide full path to phylogenetic tree in newick format");
	trait_relax.seq_data_path = io.PromptUserForString("\n>Provide full path to codon sequence alignment");
	trait_relax.json_path = io.PromptUserForString("\n>Provide required output file suffix (will be created in the codon sequence alignment directory)");
	
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
function trait_relax.declare_model_parameters() {
	
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
function trait_relax.set_codon_model() {
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
 * @name trait_relax.do_perliminary_fitting
 * @usage perform preliminary fitting of MG-REV to the codon alignment to receive better initial values for the codon model parameters (should reduce running time)
 * The initial values for the sequence model parameters of TraitRELAX are obtained through optimization of simpler models to the data. 
 * This procedure contains two optimization steps. 
 * First, we fit a simple GTR model to nucleotide data. The obtained parameters will be used as initial parameters to MG-REV model. 
 * Second, we fit the MG-REV model. This fitting is divided into two steps. 
 * Third, 
 */
function trait_relax.do_perliminary_fitting() {
	
	/* duration measurement checkpoint */
	selection.io.startTimer (trait_relax.json [terms.json.timers], "Preliminary sequence model fitting", 1);
	
	/* get the initial GTR model parameters by optimizing based on the original tree */
	namespace trait_relax {
		doGTR_FixedBLs ("trait_relax");					   
	}
	// set the global MLEs of the GTR model according to the initial GTR parameters of the RELAX model (should reduce optimization duration)
	estimators.fixSubsetOfEstimates(trait_relax.gtr_results, trait_relax.gtr_results[terms.global]); 

	/* get the initial MG-REV model parameters by optimizing based on the original tree */
	namespace trait_relax {
		scaler_prefix = "trait_relax.scaler";
		doPartitionedMG_FixedBLs ("trait_relax", FALSE);
	}
	io.ReportProgressMessageMD ("TraitRELAX", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

	trait_relax.final_partitioned_mg_results = custom_functions.FitMGREV_FixedBLs (trait_relax.filter_names, trait_relax.trees, trait_relax.codon_data_info [terms.code], {
		terms.run_options.model_type: terms.local,
		terms.run_options.partitioned_omega: trait_relax.selected_branches,
	}, trait_relax.partitioned_mg_results);

	io.ReportProgressMessageMD("TraitRELAX", "codon-refit", "* " + selection.io.report_fit (trait_relax.final_partitioned_mg_results, 0, trait_relax.codon_data_info[terms.data.sample_size]));

	// extract a global omega value from the local syn-rate and nonsyn-rate estimated for each branch
	trait_relax.global_dnds = selection.io.extract_global_MLE_re (trait_relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);
	trait_relax.report_dnds = {};

	// document the initial values in the json
	utility.ForEach (trait_relax.global_dnds, "_value_", '
		io.ReportProgressMessageMD ("TraitRELAX", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));
		trait_relax.report_dnds [(regexp.FindSubexpressions (_value_[terms.description], "^" + terms.parameters.omega_ratio + ".+\\*(.+)\\*$"))[1]] = {"0" : {terms.json.omega_ratio : _value_[terms.fit.MLE], terms.json.proportion : 1}};
	');


	// store initial MG model parameters in the json file
	selection.io.json_store_lf_GTR_MG94 (trait_relax.json,
								trait_relax.MG94_name,
								trait_relax.final_partitioned_mg_results[terms.fit.log_likelihood],
								trait_relax.final_partitioned_mg_results[terms.parameters],
								trait_relax.sample_size,
								utility.ArrayToDict (utility.Map (trait_relax.global_dnds, "_value_", "{'key': _value_[terms.description], 'value' : Eval({{_value_ [terms.fit.MLE],1}})}")),
								(trait_relax.final_partitioned_mg_results[terms.efv_estimate])["VALUEINDEXORDER"][0],
								trait_relax.display_orders[trait_relax.MG94_name]);

	// store synonymous rate and non-synonymous rate for each branch (as these parameters are unique per branch)
	utility.ForEachPair (trait_relax.filter_specification, "_key_", "_value_",
		'selection.io.json_store_branch_attribute(trait_relax.json,trait_relax.MG94_name, terms.branch_length, trait_relax.display_orders[trait_relax.MG94_name],
												 _key_,
												 selection.io.extract_branch_info((trait_relax.final_partitioned_mg_results[terms.branch_length])[_key_], "selection.io.branch.length"));');
												 
	// set the omega parameters from the MG-REV model optimization to be the ones of the RELAX reference and test models
	trait_relax.ge_guess = relax.DistributionGuess(utility.Map (selection.io.extract_global_MLE_re (trait_relax.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio + ".+test.+"), "_value_", "_value_[terms.fit.MLE]"));
	trait_relax.distribution = models.codon.BS_REL.ExtractMixtureDistribution(trait_relax.reference);
	// assumption: this extracts 3 initial omegas from the single global dN/dS ratio evaluated when optimizing the MG-REV model
	parameters.SetStickBreakingDistribution (trait_relax.distribution, trait_relax.ge_guess);
	

	selection.io.stopTimer (trait_relax.json [terms.json.timers], "Preliminary sequence model fitting");
}


/**
 * @name trait_relax.set_model_components
 * @usage creates the model instance and initialize its parameters to the values of some preliminary fitting in order to reduce optimization running time in practice
 * @The traitRalex model includes the following components: 
 * 1. parameters of the character model: mu and pie_0
 * 2. parameters of the sequence model. The sequence model consists of two components (test and reference). Each of these has the following parameters, whose values are identical in the two components: 
 *2a. GTR parameters (5 transition parameters and 3 frequency parameters)
 *2b. the ka/ks ratio parameters: 5 parameters in total. There are 3 omega's and 2 proportions parameters
 *2c. The selection intensity parameter: k - this parameter is free to vary only in the test component. In the reference component it is fixed to one. 
 */
function trait_relax.set_model_components() {
	
	/* declare the model parameters that are additional to the standard BS-REL model */
	trait_relax.declare_model_parameters();
	
	/* set the codon model as two instances that are extensions of the BS-REL model: test and reference  */
	trait_relax.set_codon_model();
	
	// do preliminary fitting on the null model parameters
	trait_relax.do_perliminary_fitting();
}


/**
 * @name trait_relax.document_null_fit
 * @usage documents the null fitting result in the json output file
 */
function trait_relax.document_null_fit() {
	trait_relax.distribution_for_json = {trait_relax.test_branches_name : utility.Map (utility.Range (trait_relax.rate_classes, 0, 1),
														 "_index_",
														 "{terms.json.omega_ratio : trait_relax.inferred_distribution [_index_][0],
														   terms.json.proportion  : trait_relax.inferred_distribution [_index_][1]}")};

	trait_relax.distribution_for_json   [trait_relax.reference_branches_name] =   trait_relax.distribution_for_json   [trait_relax.test_branches_name];

	selection.io.json_store_lf (trait_relax.json,
								trait_relax.null_name ,
								trait_relax.null_model.fit[terms.fit.log_likelihood],
								trait_relax.null_model.fit[terms.parameters] + 9 , // +9 comes from CF3x4
								trait_relax.codon_data_info[terms.data.sample_size],
								trait_relax.distribution_for_json,
								trait_relax.display_orders[trait_relax.null_name]
							);

	selection.io.json_store_branch_attribute(trait_relax.json, trait_relax.null_name, terms.branch_length, trait_relax.display_orders[trait_relax.null_name],
												 0,
												 selection.io.extract_branch_info((trait_relax.null_model.fit[terms.branch_length])[0], "selection.io.branch.length"));	
}


/**
 * @name trait_relax.optimize_null_model
 * @usage requests input paths from the user and processes them into hyphy instances
 */
function trait_relax.optimize_null_model() {
	
	selection.io.startTimer (trait_relax.json [terms.json.timers], "TraitRELAX null model fitting", 2);
	
	/* optimize the character model parameters (independent of the sequence model in the null case) */
	// set the character grid
	character_grid = trait_relax.set_character_grid(&trait_relax.trait_dataset, trait_relax.trees["0"], mu_interval, pie_0_interval);
	
	// memory capacity on dummy data is 173m

	// optimize the joint model using a two step procedure until convergence
	best_char_fit = {};
	best_char_fit["mu"] = -1;
	best_char_fit["pie_0"] = -1;
	best_log_likelihood = -100000000000000000000000000000000000; //-inf didn't work
	
	for (m=0; m<Abs(character_grid["mu"]); m+=1) {
		for (p=0; p<Abs(character_grid["pie_0"]); p+=1) {
			
			// set initial character model parameters
			mu = (character_grid["mu"])[m];
			pie_0 = (character_grid["pie_0"])[p];
			
			trait_likelihood_function = trait_relax.generate_character_likelihood_function(mu, pie_0, "trait_relax.trait", trait_relax.trait_filter_names, trait_relax.trait_trees);
			log_likelihood = estimators.ComputeLF(&trait_likelihood_function);
			
			if (log_likelihood > best_log_likelihood) {
				best_log_likelihood = log_likelihood; // update the best achieved log likelihood
				best_char_fit["mu"] = mu;
				best_char_fit["pie_0"] = pie_0;
			}
			
			// delete trait_likelihood_function
			// need to clear parameters here as well!
			custom_functions.DeleteLikelihoodFunction_FixedBLs(trait_likelihood_function);
			
		}
	}
	
	// compute the log likelihood of the character model given the best parameters
	trait_likelihood_function = trait_relax.generate_character_likelihood_function(best_char_fit["mu"], best_char_fit["pie_0"], "trait_relax.trait", trait_relax.trait_filter_names, trait_relax.trait_trees);
	trait_relax.null_model.fit.character_log_likelihood = estimators.ComputeLF(&trait_likelihood_function);
	
	// delete trait_likelihood_function
	custom_functions.DeleteLikelihoodFunction_FixedBLs(trait_likelihood_function);

	/* optimize the sequence model (RELAX null), using a maximum parsimony based partition
	/* note the difference from RELAX.bf: since in RELAX.bf, the alternative model is optimized first, the null model is bing optimized by restricting k=1 after the optimization of the alternative model. in TraitRELAX.bf, we optimize the null model first, and since we don't have any prior division of the tree branches into test and reference */
				  
	// generate a maximum parsimony based partition of the tree branches into test and reference, using the user given tree and the character states of the species
	parsimonious_model_map = {"trait_relax.test" : utility.Filter (trait_relax.selected_branches, '_value_', '_value_ == trait_relax.test_branches_name'),
                    "trait_relax.reference" : utility.Filter (trait_relax.selected_branches, '_value_', '_value_ == trait_relax.reference_branches_name')};

	// optimize the null model while preserving the relations between branch lengths
	io.ReportProgressMessageMD ("TraitRELAX", "null", "Fitting the null (K := 1) model");
	

	// constraint k to be 1
	parameters.SetConstraint (model.generic.GetGlobalParameter (relax.test , terms.relax.k), terms.parameters.one, terms.global);
	
	// fit the RELAX null model to the data - here on dummy data memory capacity is 173m (same as RELAX)
	trait_relax.null_model.fit = custom_functions.FitLF_FixedBLs (trait_relax.filter_names, trait_relax.trees, None, { "0" : parsimonious_model_map}, 
																		trait_relax.final_partitioned_mg_results, trait_relax.model_object_map, 
																		{terms.run_options.retain_lf_object: TRUE, terms.run_options.proportional_branch_length_scaler: trait_relax.scaler, "custom_lf_formula": FALSE});
																		
	// extract the log likelihood of the sequence model (RELAX) given the best parameters
	trait_relax.null_model.fit.sequence_log_likelihood = trait_relax.null_model.fit[utility.getGlobalValue("terms.fit.log_likelihood")];
	
	// sum the log likelihood of the character model and log likelihood of the sequence model to achieve the log likelihood of the TraitRELAX null model
	trait_relax.null_model.fit.log_likelihood = trait_relax.null_model.fit.character_log_likelihood + trait_relax.null_model.fit.sequence_log_likelihood;
	
	// clear the constraint on k (for the purpose of the alternative model fitting)
	parameters.RemoveConstraint(model.generic.GetGlobalParameter (relax.test , terms.relax.k));
	
	// report null fitting results
	io.ReportProgressMessageMD ("TraitRELAX", "null", "* " + selection.io.report_fit (trait_relax.null_model.fit, 9, trait_relax.codon_data_info[terms.data.sample_size]));

	io.ReportProgressMessageMD("TraitRELAX", "null", "* The following rate distribution for test/reference branches was inferred");
	trait_relax.inferred_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistributionFromFit (trait_relax.test, trait_relax.null_model.fit)) % 0;
	selection.io.report_dnds (trait_relax.inferred_distribution);
	
	// document fitting results in the json file
	trait_relax.document_null_fit();
	
	//  here on dummy data memory capacity is 1357m (way bigger than RELAX), which is 247m)
	selection.io.stopTimer (trait_relax.json [terms.json.timers], "TraitRELAX null model fitting"); 
	
}


/**
* @name trait_relax.set_starting_point
* @usage sets the initial values of the sequence model parameters in the alternative model optimization
*/
function trait_relax.set_starting_point() {
	
	// if user chose to set the starting point according to the RELAX alternative MLEs -> perform RELAX alternative fitting to the data based on the maximum parsimony partition of the tree branches
	if (trait_relax.alternative_starting_point == "RELAX alternative") {
		
		// fit the RELAX alternative model to the data, based on the maximum parsimonious partition of the tree branches
		trait_relax.alternative_model.initial_parameters = custom_functions.FitLF_FixedBLs (trait_relax.filter_names, trait_relax.trees, None, { "0" : parsimonious_model_map}, 
																		trait_relax.final_partitioned_mg_results, trait_relax.model_object_map, 
																		{terms.run_options.retain_lf_object: TRUE, terms.run_options.proportional_branch_length_scaler: trait_relax.scaler, "custom_lf_formula": FALSE});
	
	// if the user chose to set the starting point according to the RELAX null model -> set the starting point to hold the MLEs achieved for the sequence model in the trait relax null model fitting 
	} else {
		trait_relax.alternative_model.initial_parameters = trait_relax.null_model.fit; 
	}
	
}


/**
 * @name trait_relax.optimize_by_grid
 * @usage performs a grid optimization of the alternative model with two layers: the first one fixes the sequence model and optimizes the character model, and the second fixes the character model and optimizes the sequence model
 */
function trait_relax.optimize_by_grid() {

	// set the character gird
	character_grid = trait_relax.set_character_grid(&trait_relax.trait_dataset, trait_relax.trees["0"], mu_interval, pie_0_interval);

	// set the initial values of the alternative model to be the same as the optimized values of the null model
	// trait_relax.alternative_model.initial_parameters will hold the initial parameters of the alternative sequence model, until convergence
	prev_log_likelihood = -100000000000000000000000000;
	curr_log_likelihood = -10000000000000000000000000;
	is_first = 1;

	// optimize the joint model using a two step procedure until convergence
	while  (curr_log_likelihood > prev_log_likelihood + min_log_likelihood_diff) {
		
		// update the log likelihood to be the best achieved to far
		prev_log_likelihood = curr_log_likelihood;
		
		/* step 1: optimize the character model parameters while fixing the sequence model parameters using grid optimization */
		best_char_fit = {};
		best_char_fit["mu"] = -1;
		best_char_fit["pie_0"] = -1;
		best_log_likelihood = -100000000000000000000000000;

		io.ReportProgressMessageMD ("TraitRELAX", "alt", "### Searching the grid for the optimized character model parameters");
		
		// iterate over all the possible character grid parameters
		for (m=0; m<Abs(character_grid["mu"]); m+=1) {
			for (p=0; p<Abs(character_grid["pie_0"]); p+=1) {
			
				// compute the log likelihood value of the joint likelihood function given the corresponding character grid parameters 
				results = trait_relax.compute_joint_likelihood((character_grid["mu"])[m], (character_grid["pie_0"])[p], trait_relax.trait_filter_names, trait_relax.trait_trees, &trait_relax.codon_data, trait_relax.trees, num_of_histories, trait_relax.approximation_method,
				trait_relax.model_object_map, &trait_relax.alternative_model.initial_parameters, trait_relax.model_object_map, FALSE, maxSimulationsNum, None, &parsimonious_model_map);
				log_likelihood = (^ results) [utility.getGlobalValue("terms.fit.log_likelihood")];
				if (log_likelihood > best_log_likelihood ) {
					best_log_likelihood = log_likelihood; // update the best achieved log likelihood
					best_char_fit["mu"] = (character_grid["mu"])[m];
					best_char_fit["pie_0"] = (character_grid["pie_0"])[p];
				}
				
			}
		}

		io.ReportProgressMessageMD ("TraitRELAX", "alt", "optimal character model parameters fixed");
		fprintf(stdout, "* mu = ", best_char_fit["mu"], "\n"); 
		fprintf(stdout, "* pie_0 = ", best_char_fit["pie_0"], "\n\n"); 
		
		
		/* step 2: optimize the sequence model (RELAX) parameters while fixing the character model parameters using conjugate gradient method */
		io.ReportProgressMessageMD ("TraitRELAX", "alt", "### optimizing RELAX parameters");

		// optimize the sequence model given the best character parameters
		results = trait_relax.compute_joint_likelihood(best_char_fit["mu"], best_char_fit["pie_0"], trait_relax.trait_filter_names, trait_relax.trait_trees, &trait_relax.codon_data, trait_relax.trees, num_of_histories, trait_relax.approximation_method, trait_relax.model_object_map, &trait_relax.alternative_model.initial_parameters, trait_relax.model_object_map, TRUE, maxSimulationsNum, None, &parsimonious_model_map);
		log_likelihood = (^ results) [utility.getGlobalValue("terms.fit.log_likelihood")];
		sitewise_blockwise_likelihoods = (^ results) ["sitewise_likelihoods"];

		io.ReportProgressMessageMD ("TraitRELAX", "alt", "RELAX parameters optimized");
		fprintf(stdout, "* Log likelihood = ", curr_log_likelihood, "\n\n"); // debug
		fprintf(stdout, "* Delta = ", curr_log_likelihood - prev_log_likelihood, "\n\n"); // debug
	
	}	
	// set the optimized alternative model parameters to the ones that triggered the convergence
	trait_relax.alternative_model.fit = trait_relax.alternative_model.initial_parameters;
	
	// compute the log likelihood of the character model given the best parameters
	trait_relax.alternative_model.fit.character_log_likelihood = best_log_likelihood;
	
	// extract the log likelihood of the sequence model given the best parameters
	trait_relax.alternative_model.fit.sequence_log_likelihood = trait_relax.alternative_model.fit[terms.fit.log_likelihood];
	
	// sum the log likelihood of the character model and the log likelihood of the sequence model to achieve the log likelihood of the TraitRELAX alternative model
	trait_relax.alternative_model.fit.log_likelihood = trait_relax.alternative_model.fit.character_log_likelihood + trait_relax.alternative_model.fit.sequence_log_likelihood;
	
}


/**
 * @name trait_relax.document_alternative_fit 
 * @usage documents the null fitting result in the json output file
 */
// WIERED INDENTATION - MAYBE A NOTEPAD++ ISSUE
function trait_relax.document_alternative_fit() {
	
	// TO DO: add here documentation of the trait_relax.p - need to add this parameter to the dictionary trait_relax.alternative_model.fit[terms.parameters] (should be done automatically
	selection.io.json_store_lf (trait_relax.json,
								trait_relax.alternative_name,
								trait_relax.alternative_model.fit[terms.fit.log_likelihood],
								trait_relax.alternative_model.fit[terms.parameters] + 9 , // +9 comes from CF3x4
								trait_relax.codon_data_info[terms.data.sample_size],
								trait_relax.distribution_for_json,
								trait_relax.display_orders[trait_relax.alternative_name]
							);

	selection.io.json_store_branch_attribute(trait_relax.json, trait_relax.alternative_name, terms.branch_length, trait_relax.display_orders[trait_relax.alternative_name],
												 0,
												 selection.io.extract_branch_info((trait_relax.alternative_model.fit[terms.branch_length])[0], "selection.io.branch.length"));


	selection.io.stopTimer (trait_relax.json [terms.json.timers], "TraitRELAX alternative model fitting");	
	
	trait_relax.distribution_for_json = {trait_relax.test_branches_name : utility.Map (utility.Range (trait_relax.rate_classes, 0, 1),
												 "_index_",
												 "{terms.json.omega_ratio : trait_relax.inferred_distribution [_index_][0],
												   terms.json.proportion  : trait_relax.inferred_distribution [_index_][1]}"),

							trait_relax.reference_branches_name : utility.Map (utility.Range (trait_relax.rate_classes, 0, 1),
												 "_index_",
												 "{terms.json.omega_ratio : trait_relax.inferred_distribution_ref [_index_][0],
												   terms.json.proportion  : trait_relax.inferred_distribution_ref [_index_][1]}")
						   };	
}


/**
 * @name trait_relax.optimize_alternative_model
 * @usage requests input paths from the user and processes them into hyphy instances
 */
function trait_relax.optimize_alternative_model() {
	
	// report to the user on alternative model fitting start time
	io.ReportProgressMessageMD ("TraitRELAX", "alt", "Fitting the alternative model to test K != 1");
	selection.io.startTimer (trait_relax.json [terms.json.timers], "TraitRELAX alternative model fitting", 3);	
	
	// select the number of character histories to generate per grid iteration
	num_of_histories = io.PromptUser("\n>Select the number of character histories to sample in the alternative model optimization", 3, 1, 10000, FALSE);	// get the required number of sampled character histories

	// request user to choose a starting point for the trait_relax alternative model optimization (either the results of relax null achieved earlier or the results of of relax alternative fitting to a partition based on the character maximum parsimony
	trait_relax.alternative_starting_point = io.SelectAnOption ({
											{"RELAX null", "[Default] initialize the sequence model parameters of the TraitRELAX alternative model to be the MLEs of the RELAX null model"}
											{"RELAX alternative", "initialize the sequence model parameters of the TraitRELAX alternative model to be the MLEs of the RELAX alternative model"}
										}, " Alternative model optimization starting point");	
	
	//request the user to provide chosen approximation method
	trait_relax.approximation_method = io.SelectAnOption ({
											{"Multiple Histories", "[Default] Fit the alternative model to multiple trees representing multiple character histories"}
											{"Expected History", "Fit the alternative model to a single expected history by averaging over multiple histories"}
										}, " TraitRELAX approximation method");
	
	// select convergence cutoff
	min_log_likelihood_diff = io.PromptUser("\n>Select the alternative model optimization convergence cutoff", 0.01, 0.000001, 0.1, FALSE);		 		// optimization minimum difference between previous LL value and current one (condition for convergence)
										
	// set trait_relax.p to be constrained to 1 or not according to the user's choice (first option is TraitPropRELAX - p<=1 and the second option in TraitRELAX - p=1)
	trait_relax.site_propertions_freedom = io.SelectAnOption ({
											{"Yes", "[Default] Allow trait dependent site proportion to be lower than 1"}
											{"No", "Fix trait dependent sites proportion to 1"}
										}, " Allow trait dependent site proportion to be lower than 1 or not");

	// apply a restriction according to the user's request
	if (trait_relax.site_propertions_freedom == "No") {
		trait_relax.p := 1;
	}
	
	// set the starting point of the sequence model in the alternative model
	trait_relax.set_starting_point();
																	
	// optimize the alternative model using the grid approach
	trait_relax.optimize_by_grid();
	
	// report the fitting results
	io.ReportProgressMessageMD("TraitRELAX", "alt", "* " + selection.io.report_fit (trait_relax.alternative_model.fit, 9, trait_relax.codon_data_info[terms.data.sample_size]));


	trait_relax.fitted.K = estimators.GetGlobalMLE (trait_relax.alternative_model.fit,terms.trait_relax.k);
	io.ReportProgressMessageMD("TraitRELAX", "alt", "* Relaxation/intensification parameter (K) = " + Format(trait_relax.fitted.K,8,2));
	io.ReportProgressMessageMD("TraitRELAX", "alt", "* The following rate distribution was inferred for branches under character 1");
	trait_relax.inferred_distribution = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (trait_relax.test)) % 0;
	selection.io.report_dnds (trait_relax.inferred_distribution);

	io.ReportProgressMessageMD("TraitRELAX", "alt", "* The following rate distribution was inferred for branches under character 0");
	trait_relax.inferred_distribution_ref = parameters.GetStickBreakingDistribution (models.codon.BS_REL.ExtractMixtureDistribution (trait_relax.reference)) % 0;
	selection.io.report_dnds (trait_relax.inferred_distribution_ref);

	// document the results in the json file
	trait_relax.document_alternative_fit();	
}


/**
 * @name trait_relax.optimize_alternative_model
 * @usage requests input paths from the user and processes them into hyphy instances
 */
function trait_relax.perform_test() {
	
	trait_relax.alpha = io.PromptUser("\n>Select the p-value threshold to use when testing for selection", 0.05, 0, 1, FALSE);

	// ask Itay: if not using prop - number of degrees of freedom should be 1 (no parameter p), right?
	parameters_num_diff = 1;
	if (trait_relax.site_propertions_freedom == "Yes") {
		parameters_num_diff = 2;
	}
	trait_relax.LRT = math.DoLRT (trait_relax.null_model.fit.log_likelihood, trait_relax.alternative_model.fit.log_likelihood, parameters_num_diff);

	// compute the bayes factors and the empirical bayes posteriors per site
	sites_bayes_measurements = trait_relax.compute_bayes_measurements(sitewise_blockwise_likelihoods, num_of_histories, trait_relax.p);

	console.log ("----\n## Test for relaxation (or intensification) of selection in association with a phenotypic trait [TraitRELAX]");
	console.log ( "Likelihood ratio test **p = " + Format (trait_relax.LRT[terms.p_value], 8, 4) + "**.");

	if (trait_relax.LRT[terms.p_value] <= trait_relax.alpha) {
		if (trait_relax.fitted.K > 1) {
			console.log (">Evidence for *intensification of selection* among phylogeny under character state 1 _relative_ to the phylogeny under character state 0 at P<="+ trait_relax.alpha);
		} else {
			console.log (">Evidence for *relaxation of selection* among among phylogeny under character state 1 _relative_ to the phylogeny under character state 0 at P<="+ trait_relax.alpha);
		}
	} else {
		console.log (">No significant evidence for relaxation (or intensification) of selection associated with character evolution at P<="+ trait_relax.alpha);
	}

	trait_relax.test_results = {};
	trait_relax.test_results["LRT"] = trait_relax.LRT;
	trait_relax.test_results[terms.trait_relax.k] = trait_relax.fitted.K;
	trait_relax.test_results[terms.trait_relax.p] = trait_relax.p;
	trait_relax.test_results["codon sites to trait association"] = sites_bayes_measurements;

	trait_relax.json [terms.json.test_results] = trait_relax.test_results;

	console.log ("----\n");
	
	selection.io.stopTimer (trait_relax.json [terms.json.timers], "Overall");

	io.SpoolJSON (trait_relax.json, trait_relax.codon_data_info [terms.json.json]);

	return trait_relax.json;
	
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
	
	trait_likelihood_function = custom_functions.CreateLFObject_FixedBLs (lf_prefix, data_filter, tree, "models.binary.ModelDescription", character_initial_parameters, {terms.run_options.retain_lf_object: TRUE, "set_freq": TRUE, "return_lf": TRUE});
	return trait_likelihood_function
}


/**
 * @name trait_relax.set_character_grid
 * @param {Dict} character_data_filter  		- dictionary holding a pointer to the character data
 * @param {Dict} tree							- dictionary holding the tree provided by the user
 * @return {Dict} character_grid	       		- dictionary holding the generated character grid in two dimensions: "mu" and "pie_0"
 * @usage sets the character grid as a matrix of mu and pie0 combinations for the character model
 */
lfunction trait_relax.set_character_grid(character_data_filter, tree, mu_interval, pie_0_interval) {
	
	character_grid = {};
	character_grid["mu"] = {};
	character_grid["pie_0"] = {};
	
	// calculate the maximum parsimony score of the character data along the tree to construct the character grid
	char_mp_score = custom_functions.getMPScore(tree, character_data_filter, FALSE); // FALSE -> return MP score and not events map
	// set the upper and lower limit of mu according to the mp score
	mu_lower_limit = char_mp_score;
	mu_upper_limit = 2*char_mp_score;
	mu_counter = 0;
	for (mu=mu_lower_limit; mu <=mu_upper_limit; mu += mu_interval) {
		(character_grid["mu"])[mu_counter] = mu;
		mu_counter += 1;
	}
	pie_0_lower_limit = 0.05; // setting above 0 to avoid having 1 being a absorbing state in the character transition matrix
	pie_0_upper_limit = 0.95; // setting below 1 to avoid having 0 being a absorbing state in the character transition matrix
	pie_0_counter = 0;
	for (pie_0=pie_0_lower_limit; pie_0 <=pie_0_upper_limit; pie_0 += pie_0_interval) {
		(character_grid["pie_0"])[pie_0_counter] = pie_0;
		pie_0_counter += 1;
	}
	
	return character_grid

	
	
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
lfunction trait_relax.compute_joint_likelihood(mu, pie_0, trait_filter_name, trait_tree, sequence_dataset, sequence_trees, histories_num, approximation_method, model_object_map, relax_initial_parameters, sequence_model_objects, optimize_sequence_parameters, max_simulations_num, run_options, base_model_map) {
	
	joint_sequence_filter_names = {}; 
	joint_sequence_trees = {}; 

	// initialize the input for the construction of the joint likelihood function
	lf_formula = "'Log((1-trait_relax.p)*SITE_LIKELIHOOD[0]+trait_relax.p*((";
	custom_functions.dup_data_filter("data_filter_0", sequence_dataset);
	joint_sequence_filter_names[0] = "data_filter_0";
	joint_sequence_trees[0] = sequence_trees[0];
	model_maps = {};
	model_maps[0] = ^ base_model_map;
	
	// create the likelihood function instance for the SM procedure
	trait_likelihood_function = trait_relax.generate_character_likelihood_function(mu, pie_0, "trait_relax.trait", trait_filter_name, trait_tree);

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
		}
		
		// add he history to the sequence dataset 
		joint_sequence_filter_names [h+1] = "data_filter_" + h+1;
		custom_functions.dup_data_filter(joint_sequence_filter_names [h+1], sequence_dataset);
		joint_sequence_trees[h+1] = history;
		model_maps[h+1] = trait_relax.generate_history_rules_map(history);
	}
	
	// delete the trait likelihood function
	custom_functions.DeleteLikelihoodFunction_FixedBLs(trait_likelihood_function);
	
	// user chose to compute the likelihood over all the character histories --> include all of them in the likelihood function formula
	if (approximation_method == "Multiple Histories") {
		lf_formula = lf_formula + ")/" + histories_num + "))'";
	
	// user chose to average the histories into a single expected history    --> average the histories and use the expected history in the joint likelihood function formula
	} else { 
		expected_history_info = trait_relax.average_histories(sequence_trees[0], &histories_durations_dict, run_options); // compute the expected history
		joint_sequence_filter_names[2] = "trait_relax.filter.default";
		joint_sequence_trees[2] = expected_history_info;
		model_maps[2] = trait_relax.generate_history_rules_map(expected_history_info);
		histories_dict[2] = {"id": expected_history_info_id, "tree": expected_history_info};
		lf_formula = lf_formula + "SITE_LIKELIHOOD[1])))'";
	}

	// generate and optimize the joint likelihood function (i.e, while optimizing the sequence model parameters)
	if (optimize_sequence_parameters == TRUE) {
		trait_relax.alternative_model.initial_parameters = custom_functions.FitLF_FixedBLs(joint_sequence_filter_names, joint_sequence_trees, lf_formula, model_maps, relax_initial_parameters, sequence_model_objects, {terms.run_options.retain_lf_object: FALSE, terms.run_options.proportional_branch_length_scaler: trait_relax.scaler, "custom_lf_formula": TRUE});
	}
	
	// generate and compute the log likelihood of the joint likelihood function (i.e, while fixing the sequence model parameters)
	else {
		trait_relax.alternative_model.initial_parameters = custom_functions.ComputeLF_FixedBLs(joint_sequence_filter_names, joint_sequence_trees, lf_formula, model_maps, relax_initial_parameters, sequence_model_objects, {terms.run_options.retain_lf_object: TRUE, terms.run_options.proportional_branch_length_scaler: trait_relax.scaler, "custom_lf_formula": TRUE});
	}
	
	return &trait_relax.alternative_model.initial_parameters;
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
	return return_set
	
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





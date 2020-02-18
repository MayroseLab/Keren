LoadFunctionLibrary("libv3/models/model_functions.bf");
LoadFunctionLibrary("libv3/convenience/regexp.bf");
LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/models/parameters.bf");
LoadFunctionLibrary ("libv3/models/frequencies.bf");
LoadFunctionLibrary ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");



/* _____________ MY OWN AUXILIARY FUNCTIONS ____________ */

/**
 * @name custom_functions_bl_opt.strReplace
 * @param {String} data_pointer	  - pointer to the data instance for which a data filter needs to be created
 * @param {String} filter_pointer - pointer that will hold the created data filter
 */
lfunction custom_functions_bl_opt.dup_data_filter(datafilter_name, data_pointer) {
	utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null"); 
	ExecuteCommands('ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "chooseGeneticCode.def");', {"0" : "Universal"});
	utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);
	DataSetFilter ^ datafilter_name = CreateFilter (^ data_pointer, 3, , , GeneticCodeExclusions);
	return	filter_name; 
}

/**
 * @name custom_functions_bl_opt.strReplace
 * @param {String} string - the string in which the regex substitution is required to occur
 * @param {String} src	  - the sub-string that needs to be replaced
 * @param {String} dst	  - the sub-string to replace src with 
 * @returns converted string
 */
lfunction custom_functions_bl_opt.strReplace (string, src, dst) {
	replacor = {2,1};
	replacor[0] = src;
	replacor[1] = dst;
	string = string^replacor;
	return string;
}


/**
 * @name getMPScore
 * @param {Dict} tree 				- dictionary representing the tree
 * @param {String} data				- name of the corresponding dataset
 * Calculated the maximum parsimony score of the tree given the data
 */
lfunction custom_functions_bl_opt.getMPScore (tree, data) {
	treeStr = tree[utility.getGlobalValue("terms.trees.newick_with_lengths")];
	Tree tree_instance = treeStr;
	nodes = trees.BranchNames(tree);
	tree_avl = tree_instance ^ 0; // get a list of the tree nodes in post-order
	
	// generate a map of node to its parent
	node_to_parent = {};
	for (i=0; i<Columns(nodes)-1; i+=1) {
		node_to_parent[nodes[i]] = (tree_avl[(tree_avl[i+1])["Parent"]])["Name"];
	}
	
	// generate a reverse map of node to its children (hould constain only internal nodes)
	node_to_children = {};
	for (i=0; i<Columns(nodes); i+=1) {
		
		if (nodes[i] % "Node0") {
			is_leaf = 0;
		} else {
			is_leaf = (tree[utility.getGlobalValue("terms.trees.partitioned")])[nodes[i]] % "leaf"
		}
		
		if (!is_leaf) {
			children = {};
			children_counter = 0;
			for (j=0; j<Columns(nodes); j+=1) {
				if (node_to_parent[nodes[j]] % nodes[i]) {
					children[children_counter] = nodes[j];
					children_counter += 1;
				}
			}
			node_to_children[nodes[i]] = {"childrenNum": children_counter, "childrenNames": children};
		}
	}
	
	// set the state of each node (could either be a single state or a union of multiple states), while counting the number of inersection events
	// currently, since the character sequence is consisted of a single site, we don't iterate over sites. in order to generalize this function, we would have to do so
	node_to_state = {};
	intersection_events_num = 0;
	for (i=0; i<Columns(nodes); i+=1) {
		if ((tree[utility.getGlobalValue("terms.trees.partitioned")])[nodes[i]] % "leaf") {
			node_to_state[nodes[i]] = (alignments.GetIthSequence(data, i))[utility.getGlobalValue("terms.data.sequence")]; 
		} else { // the node is internal and therefore has children
			children_num = (node_to_children[nodes[i]])["childrenNum"];
			children_names = (node_to_children[nodes[i]])["childrenNames"];
			children_states = {};
			states_counter = 0; // in our case, can only go up to 2
			for (j=0; j<children_num; j+=1) {
				child_name = children_names[j];
				child_state = node_to_state[child_name]; // since we traverse the tree in post order traversal, we are assured to visit the children before the parent
				is_new = 1;
				for (k=0; k<=states_counter; k+=1) {
					if (children_states[k] % child_state) {
						is_new = 0;
					}
				}
				if (is_new == 1) {
					is_both = child_state % "both";
					if (!is_both) {
						children_states[states_counter] = child_state;
						states_counter += 1;
					}
				}
			}
			if (states_counter > 1) {
				intersection_events_num += 1;
				node_to_state[nodes[i]] = "both";
			} else {
				node_to_state[nodes[i]] = child_state;
			}
		}
	}
	
	// return the number of intersection events, which is equal to the MP score
	return intersection_events_num;  
}


/**
 * @name custom_functions_bl_opt.fix_branch_lengths
 * @param {String} id
 * @param {Dictionary} model_info
 * @param {List} branch_names
 * @param {Dictionary} branch_lengths
 * Must be a function (and not lfunction) in order to apply a constraint on a model external to the function (doesn't apply otherwise)
 */
 // tried adjusting the ID to be non global: adjusted_id = custom_functions_bl_opt.strReplace(id, "([a-zA-Z][^.]+)\.", ""); - didn't work:(
function custom_functions_bl_opt.fix_branch_lengths (id, model_info, branch_names, branch_lengths) {
	// doing nothing
}


/**
 * @name custom_functions_bl_opt.clear_local_constraints
 * @param {Long} lf_id - pointer to a likelihood function
 * Must be a function (and not lfunction) in order to clear a constraint on a model external to the function (doesn't apply otherwise)
 */
function custom_functions_bl_opt.clear_local_constraints(lf_id) {
	// doing nothing
}


/**
 * @name custom_functions_bl_opt.DeleteLikelihoodFunction_FixedBLs
 * @param {Long} lfID - pointer to the a likelihood function with constrained branches that need to be cleared upon deletion
 */
function custom_functions_bl_opt.DeleteLikelihoodFunction_FixedBLs(lfId) {
	
	// get the information of the likelihood function
	GetString(lf_info, ^ lfId, -1);
	
	local_constraints = lf_info[utility.GetGlobalValue("terms.parameters.local_constrained")];
	constraints_num = Columns(local_constraints);
	for (c=0; c< constraints_num; c+=1) {
		ClearConstraints(local_constraints[C]);
	}
	
	// delete the likelihood function
	DeleteObject(^ lfId); 
}





/* ______________________ FUNCTIONS BASED ON FUNCTIONS FROM THE LIBV3 OF HYPHY ______________ */

/**
 * @name custom_functions_bl_opt.FitLF_FixedBLs
 * @param {Dict} data_filter  		- a vector of DataFilters
 * @param {Dict} tree  				- a vector of Trees
 * @param {String} lf_formula		- custom formula for the likelihood function computation
 * @param {Dict} model_map 			- map that holds the labeling of the trees branches into models 
 * @param {Dict} initial_values 	- dictionary that holds the initial values of the models parameters 
 * @param {Dict} model_objects 		- dictionary that holds the models instances
 * @param {Dict} run_options 		- dictionary that holds the function run settings 
 * @returns LF results after optimization
 */
lfunction custom_functions_bl_opt.FitLF_FixedBLs (data_filter, tree, lf_formula, model_map, initial_values, model_objects, run_options) {
	
	// memory capacity  on dummy data is 2858m
	/* console.log("breakpoint at line 223");
	while (TRUE) {
		// DO NOTHING
	} */
	
    if (Type(data_filter) == "String") {
        return custom_functions_bl_opt.FitLF_FixedBLs ({
            {
                data_filter__
            }
        }, {
            "0": tree
        },
        {
            "0" : model_map
        },
        initial_values, model_objects, run_options);
    }

    components = utility.Array1D(data_filter);
	
	if (run_options["custom_lf_formula"] == TRUE) {
		lf_components = {2 * components + 1, 1};
    } else {
		lf_components = {2 * components, 1};
	}
	
	// memory capacity  on dummy data is 2858m
		
    for (i = 0; i < components; i += 1) {
        lf_components[2 * i] = data_filter[i];
        lf_components[2 * i + 1] = &tree_id + "_" + i;
        custom_functions_bl_opt.ApplyModelToTree_FixedBLs (lf_components[2*i + 1], tree[i], model_objects, model_map[i]);
    }
	
	// memory capacity  on dummy data is 2858m
	
	
	if (run_options["custom_lf_formula"] == TRUE) {
		lf_components[2 * components] = lf_formula;
	}
		
    lf_id = &likelihoodFunction;
	num = io.PromptUser("\n>creating lf from FitLF. Press some number", 0, 1, 10, FALSE);
    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `lf_id` = (`&lf_components`)");
	num = io.PromptUser("\n>created lf from FitLF. Press some number", 0, 1, 10, FALSE);
	// memory capacity  on dummy data is 2858m
	
    df = 0;

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
            df = estimators.ApplyExistingEstimates("`&likelihoodFunction`", model_objects, initial_values, run_options[utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler")]);
    }
	
	

    if (utility.Has (run_options,utility.getGlobalValue("terms.run_options.apply_user_constraints"),"String")) {
        df += Call (run_options[utility.getGlobalValue("terms.run_options.apply_user_constraints")], lf_id, lf_components, data_filter, tree, model_map, initial_values, model_objects);
    }
	
	fprintf(stdout, "Optimize: start\n"); // debug
   	Optimize (mles, likelihoodFunction);
	fprintf(stdout, "Optimize: stop\n"); // debug
	
	// clear the constraints on the branch lengths
	custom_functions_bl_opt.clear_local_constraints(&likelihoodFunction);

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }

    results = estimators.ExtractMLEs( & likelihoodFunction, model_objects);

    results[utility.getGlobalValue ("terms.fit.log_likelihood")] = mles[1][0];
    results[utility.getGlobalValue ("terms.parameters")] = mles[1][1] + df;
	ConstructCategoryMatrix(sitewise_likelihoods ,likelihoodFunction ,SITE_LOG_LIKELIHOODS);
	results["sitewise_likelihoods"] = sitewise_likelihoods;

    results[utility.getGlobalValue ("terms.fit.filters")] = {
        1,
        components
    };

    for (i = 0; i < components; i += 1) {
        (results[utility.getGlobalValue ("terms.fit.filters")])[i] = lf_components[2 * i];

    }

    if (run_options[utility.getGlobalValue("terms.run_options.retain_lf_object")]) {
        results[utility.getGlobalValue("terms.likelihood_function")] = & likelihoodFunction;
    } else {
        DeleteObject(likelihoodFunction);
    }
    return results;
}


/**
 * @name custom_functions_bl_opt.FitLF_FixedBLs
 * @param {Dict} data_filter  		- a vector of DataFilters
 * @param {Dict} tree  				- a vector of Trees
 * @param {String} lf_formula		- custom formula for the likelihood function computation
 * @param {Dict} model_map 			- map that holds the labeling of the trees branches into models 
 * @param {Dict} initial_values 	- dictionary that holds the initial values of the models parameters 
 * @param {Dict} model_objects 		- dictionary that holds the models instances
 * @param {Dict} run_options 		- dictionary that holds the function run settings 
 * @returns LF results after computation
 */
lfunction custom_functions_bl_opt.ComputeLF_FixedBLs (data_filter, tree, lf_formula, model_map, initial_values, model_objects, run_options) {
    
	if (Type(data_filter) == "String") {
        return custom_functions_bl_opt.FitLF_FixedBLs ({
            {
                data_filter__
            }
        }, {
            "0": tree
        },
        {
            "0" : model_map
        },
        initial_values, model_objects, run_options);
    }

    components = utility.Array1D(data_filter);
	
	if (run_options["custom_lf_formula"] == TRUE) {
		lf_components = {2 * components + 1, 1};
    } else {
		lf_components = {2 * components, 1};
	}
		
    for (i = 0; i < components; i += 1) {
        lf_components[2 * i] = data_filter[i];
        lf_components[2 * i + 1] = &tree_id + "_" + i;
        custom_functions_bl_opt.ApplyModelToTree_FixedBLs (lf_components[2*i + 1], tree[i], model_objects, model_map[i]);
    }
	
	if (run_options["custom_lf_formula"] == TRUE) {
		lf_components[2 * components] = lf_formula;
	}
		
    lf_id = &likelihoodFunction;
	num = io.PromptUser("\n>creating lf from ComputeLF. Press some number to continue", 0, 1, 10, FALSE);
    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `lf_id` = (`&lf_components`)");
	num = io.PromptUser("\n>created lf from ComputeLF. Press some number to continue", 0, 1, 10, FALSE);
	
    df = 0;

	num = io.PromptUser("\n>applying existing estimates. Press some number to continue", 0, 1, 10, FALSE);
    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
            df = estimators.ApplyExistingEstimates("`&likelihoodFunction`", model_objects, initial_values, run_options[utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler")]);
    }
	num = io.PromptUser("\n>applied existing estimates. Press some number to continue", 0, 1, 10, FALSE);

	num = io.PromptUser("\n>applying user constraints. Press some number to continue", 0, 1, 10, FALSE);
    if (utility.Has (run_options,utility.getGlobalValue("terms.run_options.apply_user_constraints"),"String")) {
        df += Call (run_options[utility.getGlobalValue("terms.run_options.apply_user_constraints")], lf_id, lf_components, data_filter, tree, model_map, initial_values, model_objects);
    }
	num = io.PromptUser("\n>applied user constraints. Press some number to continue", 0, 1, 10, FALSE);
	
	num = io.PromptUser("\n>computing log likelihood. Press some number to continue", 0, 1, 10, FALSE);
   	log_likelihood = estimators.ComputeLF (& likelihoodFunction);
	num = io.PromptUser("\n>computed log likelihood. Press some number to continue", 0, 1, 10, FALSE);

	// clear the constraints on the branch lengths
	num = io.PromptUser("\n>clearing constraints. Press some number to continue", 0, 1, 10, FALSE);
	custom_functions_bl_opt.clear_local_constraints(& likelihoodFunction);
	num = io.PromptUser("\n>cleared constraints. Press some number to continue", 0, 1, 10, FALSE);

	
    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }
	
	num = io.PromptUser("\n>constructing results. Press some number to continue", 0, 1, 10, FALSE);
	results = {};

    results[utility.getGlobalValue ("terms.fit.log_likelihood")] = log_likelihood;
    results[utility.getGlobalValue ("terms.parameters")] = initial_values;
	ConstructCategoryMatrix(sitewise_likelihoods ,likelihoodFunction ,SITE_LOG_LIKELIHOODS);
	results["sitewise_likelihoods"] = sitewise_likelihoods;
	
    results[utility.getGlobalValue ("terms.fit.filters")] = {
        1,
        components
    };

    for (i = 0; i < components; i += 1) {
        (results[utility.getGlobalValue ("terms.fit.filters")])[i] = lf_components[2 * i];

    }
	num = io.PromptUser("\n>constructed results. Press some number to continue", 0, 1, 10, FALSE);

	num = io.PromptUser("\n>deleting lf. Press some number to continue", 0, 1, 10, FALSE);
    if (run_options[utility.getGlobalValue("terms.run_options.retain_lf_object")]) {
        results[utility.getGlobalValue("terms.likelihood_function")] = & likelihoodFunction;
    } else {
        DeleteObject(likelihoodFunction); // can only be deleted if the lf_compomnents contain the same instance just once
    }
	num = io.PromptUser("\n>deleted lf. Press some number to continue", 0, 1, 10, FALSE);
    return results;
}


/**
 * @name custom_functions_bl_opt.ApplyModelToTree_FixedBLs
 * @param id
 * @param tree
 * @param model_list
 * @param rules
 * @usage applies models on branches of the tree by the name 'id' (parameter id) according to a given mapping (parameter rules) and a list of models (parameter model_list). restricts branch lengths to be fixed if and only if the constraint is a linear combination with regards to the constrained parameter
 */
function custom_functions_bl_opt.ApplyModelToTree_FixedBLs (id, tree, model_list, rules) {
	
	// extract the tree branch names and lengths for the procedure of constraining branch lengths 
	branch_names = trees.BranchNames(tree);
	branch_lengths = tree[utility.getGlobalValue("terms.branch_length")];

	if (Type (rules) == "AssociativeList") {
	    // this has the form
	    // model id : list of branches to apply the model (as a string COLUMN matrix with branch names,
	    // or a dictionary where keys are the branch names)
	    // OR
	    // DEFAULT : model id
	    if (Abs (rules["DEFAULT"])) {
            ExecuteCommands ("UseModel (" + rules["DEFAULT"] + ");
                              Tree `id` = " + tree["string_with_lengths"] + ";
                              ");
	    } else {
            ExecuteCommands ("UseModel (USE_NO_MODEL);
                              Tree `id` = " + tree["string_with_lengths"] + ";
                              ");
	    }
		
		// in the handled data filter is not the default one, get the list of rules 
	    custom_functions_bl_opt.ApplyModelToTree_FixedBLs.ids = Rows (rules); // = {{"trait_relax.test", "trait_relax.reference"}}
	    for (custom_functions_bl_opt.ApplyModelToTree_FixedBLs.k = 0; custom_functions_bl_opt.ApplyModelToTree_FixedBLs.k < Abs (rules); custom_functions_bl_opt.ApplyModelToTree_FixedBLs.k += 1) {
	        custom_functions_bl_opt.ApplyModelToTree_FixedBLs.name = custom_functions_bl_opt.ApplyModelToTree_FixedBLs.ids[custom_functions_bl_opt.ApplyModelToTree_FixedBLs.k];
	        if ( custom_functions_bl_opt.ApplyModelToTree_FixedBLs.name != "DEFAULT") {
                custom_functions_bl_opt.ApplyModelToTree_FixedBLs.list = rules[custom_functions_bl_opt.ApplyModelToTree_FixedBLs.name];
                if (Type (custom_functions_bl_opt.ApplyModelToTree_FixedBLs.list) == "AssociativeList") {
                    custom_functions_bl_opt.ApplyModelToTree_FixedBLs.list = Rows (custom_functions_bl_opt.ApplyModelToTree_FixedBLs.list);
                }

			   // for each model in the list, apply it on the relevant branches and fix the branch length under the corresponding model (by constraining the local parameter of the branch)
                for (custom_functions_bl_opt.ApplyModelToTree_FixedBLs.b = 0; custom_functions_bl_opt.ApplyModelToTree_FixedBLs.b < Columns (custom_functions_bl_opt.ApplyModelToTree_FixedBLs.list); custom_functions_bl_opt.ApplyModelToTree_FixedBLs.b += 1) {
                    // assign to the branch in place b in the list of branches under model custom_functions_bl_opt.ApplyModelToTree_FixedBLs.apply_model the model custom_functions_bl_opt.ApplyModelToTree_FixedBLs.apply_model
					ExecuteCommands ("SetParameter (`id`." + custom_functions_bl_opt.ApplyModelToTree_FixedBLs.list[custom_functions_bl_opt.ApplyModelToTree_FixedBLs.b] + ",MODEL," + custom_functions_bl_opt.ApplyModelToTree_FixedBLs.name + ")");
				}
				// need to somehow get the info (not names!) of the model that goes by the name apply_model
				custom_functions_bl_opt.fix_branch_lengths (id, model_list[custom_functions_bl_opt.ApplyModelToTree_FixedBLs.name], branch_names, branch_lengths);
            }
	    }
	// if there is a single model that applies on all the trees, simply apply it on the tree and constrain all the local parameters of the branches under the same model
	} else {
	    // TO DO: REMOVE HARDCODING
		custom_functions_bl_opt.ApplyModelToTree_FixedBLs.modelID = model_list[model_list ["INDEXORDER"][0]]; // ASK STEPHANIE: IS custom_functions_bl_opt.ApplyModelToTree_FixedBLs.modelID[terms.id] THE MODEL INSTANCE? BECAUSE ABOVE, IN SetParameter(), custom_functions_bl_opt.ApplyModelToTree_FixedBLs.apply_model WASN'T A PARAMETER
		ExecuteCommands ("UseModel (" + custom_functions_bl_opt.ApplyModelToTree_FixedBLs.modelID[terms.id] + ");
						  Tree `id` = " + tree["string"] + ";
						  ");

		// fix the branches lengths
		custom_functions_bl_opt.fix_branch_lengths (id, custom_functions_bl_opt.ApplyModelToTree_FixedBLs.modelID, branch_names, branch_lengths);
	}
}


/**
 * @name custom_functions_bl_opt.FitSingleModel_Ext_FixedBLs
 * @param {DataFilter} data_filter
 * @param {Tree} tree
 * @param {Dict} model
 * @param {Matrix} initial_values
 * @param {Dict} run_options
 * @returns results (unlike estimators.FitSingleModel_Ext, this function also fixes the branch lengths before optimization)
 */
lfunction custom_functions_bl_opt.FitSingleModel_Ext_FixedBLs (data_filter, tree, model_template, initial_values, run_options) {
	
    this_namespace = (&_);
    this_namespace = this_namespace[0][Abs (this_namespace)-3];

    df = custom_functions_bl_opt.CreateLFObject_FixedBLs (this_namespace, data_filter, tree, model_template, initial_values, run_options);

	fprintf(stdout, "Optimize: start\n"); // debug
	Optimize(mles, likelihoodFunction);
	fprintf(stdout, "Optimize: stop\n"); // debug
	
	if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }

    model_id_to_object = {
        (this_namespace + ".model"): user_model
    };

	results = estimators.ExtractMLEs( & likelihoodFunction, model_id_to_object);


    results[utility.getGlobalValue("terms.fit.log_likelihood")] = mles[1][0];
    results[utility.getGlobalValue("terms.parameters")] = mles[1][1] + (user_model [utility.getGlobalValue("terms.parameters")]) [utility.getGlobalValue("terms.model.empirical")] + df;
	
    if (option[utility.getGlobalValue("terms.run_options.retain_model_object")]) {
		results[utility.getGlobalValue("terms.model")] = model_id_to_object;
    }

   if (run_options[utility.getGlobalValue("terms.run_options.retain_lf_object")]) {
        results[utility.getGlobalValue("terms.likelihood_function")] = & likelihoodFunction;
    } else {
        DeleteObject(likelihoodFunction);
    }

    return results;
}


/**
 * @name custom_functions_bl_opt.FitMGREV_FixedBLs
 * @param {DataFilter} codon_data
 * @param {Tree} tree
 * @param {String} genetic_code
 * @param {Dictionary} option
 * @param {Dictionary} initial_values
 * @returns MGREV results
 */
lfunction custom_functions_bl_opt.FitMGREV_FixedBLs (codon_data, tree, genetic_code, option, initial_values) {

    //TODO (TEAM): Where is data_filter being set?
    if (Type(data_filter) == "String") {
        return custom_functions_bl_opt.FitMGREV_FixedBLs({
                {
                    codon_data__
                }
            }, {
                "0": tree
            },
            genetic_code,
            option,
            initial_values)
    }

	// extract the number of likelihood function components ( component = data_filter + tree)
    components = utility.Array1D(codon_data);

    lf_components = {
        2 * components,
        1
    };

	// extract each component from the provided dictionary
    for (i = 0; i < components; i += 1) {
        GetDataInfo(fi, ^ (codon_data[i]), "PARAMETERS");
        DataSetFilter * ("filter_" + i) = CreateFilter( ^ (codon_data[i]), 3, '', '', fi["EXCLUSIONS"]);
        // need to do this for global references
        lf_components[2 * i] = "filter_" + i;
    }

    name_space = & model_MGREV;

	// define the MG-REV model
	// KEREN: note that the option of propertional branch scaler is not set to the defnie model function - maybe it's not a conventional parameter
    mg_rev = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        name_space, {
            "0": parameters.Quote(option[utility.getGlobalValue("terms.run_options.model_type")]),
            "1": genetic_code
        },
        codon_data,
        None);
	
    df = 0;
    model_assignment = {
        "default": mg_rev
    };
    rules = None;
    model_id_to_object = {
        name_space: mg_rev
    };

    // load the MG-REV model on the branches of each tree in the partition
	for (i = 0; i < components; i += 1) {
        lf_components[2 * i + 1] = "tree_" + i;
        custom_functions_bl_opt.ApplyModelToTree_FixedBLs(Eval("&`lf_components[2*i + 1]`"), tree[i], model_assignment, None); 
    }

	// set the components of the model according to the number of omega partitions (3 in case of MG-REV model)
    partition_omega = {};
    if (option[utility.getGlobalValue("terms.run_options.model_type")] == utility.getGlobalValue("terms.local") && Type(option[utility.getGlobalValue("terms.run_options.partitioned_omega")]) == "AssociativeList") {
        /**
            Assumes that option["partitioned-omega"] is a dictionary where each partition has
            an entry (0-index based), which itself is a dictionary of the form: "branch-name" : "branch-set"
        */
        utility.ForEach(option[utility.getGlobalValue("terms.run_options.partitioned_omega")], "_value_", "utility.AddToSet(`&partition_omega`,utility.Values(_value_))");
    }


    if (Abs(partition_omega)) {

        /**
            declare the global ratios for each branch set
            and add them to the model parameter set
        */

        new_globals = {};
        utility.ForEachPair(partition_omega, "_key_", "_value_",
            '`&new_globals` [_key_] = (`&name_space` + ".omega_" + Abs (`&new_globals`)); model.generic.AddGlobal (`&mg_rev`, `&new_globals` [_key_] , (utility.getGlobalValue("terms.parameters.omega_ratio")) + " for *" + _key_ + "*")');
        parameters.DeclareGlobal(new_globals, None);

        /**
            now replicate the local constraint for individual branches
        */

        alpha = model.generic.GetLocalParameter(mg_rev, utility.getGlobalValue("terms.parameters.synonymous_rate"));
        beta = model.generic.GetLocalParameter(mg_rev, utility.getGlobalValue("terms.parameters.nonsynonymous_rate"));
        io.CheckAssertion("None!=`&alpha` && None!=`&beta`", "Could not find expected local synonymous and non-synonymous rate parameters in \`custom_functions_bl_opt.FitMGREV_FixedBLs\`");

        apply_constraint: = component_tree + "." + node_name + "." + beta + ":=" + component_tree + "." + node_name + "." + alpha + "*" + new_globals[branch_map[node_name]];

        for (i = 0; i < components; i += 1) {
            component_tree = lf_components[2 * i + 1];
            ClearConstraints( * component_tree); // as a result of this global clearing, the constraints applied to fix the branch lengths in ApplyModelToTree_FixedBLs (line 312) are cleared as well
												 // therefore, I must re-apply them in lines 367-373
												 // when I re-apply, the model parameter omega_0 disappears
            branch_map = (option[utility.getGlobalValue("terms.run_options.partitioned_omega")])[i];
            component_branches = BranchName( * component_tree, -1);
            for (j = 0; j < Columns(component_branches) - 1; j += 1) {
                /**
                    -1 in the upper bound because we don't want to count the root node
                */

                node_name = (component_branches[j]);
                ExecuteCommands(apply_constraint);
            }
        }
		
    } else {}

    LikelihoodFunction likelihoodFunction = (lf_components);

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
        df += estimators.ApplyExistingEstimates("`&likelihoodFunction`", model_id_to_object, initial_values, option[utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler")]);
    }

	fprintf(stdout, "Optimize: start\n"); // debug
    Optimize(mles, likelihoodFunction);
	fprintf(stdout, "Optimize: stop\n"); // debug
	
	// clear the constraints on the branch lengths
	custom_functions_bl_opt.clear_local_constraints(& likelihoodFunction); 

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }

    results = estimators.ExtractMLEs( & likelihoodFunction, model_id_to_object);


    results[utility.getGlobalValue("terms.fit.log_likelihood")] = mles[1][0];
    results[utility.getGlobalValue("terms.parameters")] = mles[1][1] + (mg_rev [utility.getGlobalValue("terms.parameters")]) [utility.getGlobalValue("terms.model.empirical")] + df;


    if (option[utility.getGlobalValue("terms.run_options.retain_lf_object")]) {
        results[utility.getGlobalValue("terms.likelihood_function")] = & likelihoodFunction;
    } else {
        DeleteObject(likelihoodFunction);
    }

    if (option[utility.getGlobalValue("terms.run_options.retain_model_object")]) {
        results[utility.getGlobalValue("terms.model")] = model_id_to_object;
    }

    return results;
}


/**
 * @name custom_functions_bl_opt.CreateLFObject_FixedBLs
 * @param {String} context
 * @param {DataFilter} data_filter
 * @param {Tree} tree
 * @param {Dictionary} model_template
 * @param {Dictionary} initial_values
 * @param {Dictionary} run_options - additional option to fix the frequency to the one given in the initial values
 */
lfunction custom_functions_bl_opt.CreateLFObject_FixedBLs (context, data_filter, tree, model_template, initial_values, run_options) {
    if (Type(data_filter) == "String") {
        return estimators.FitSingleModel_Ext ({
            {
                data_filter__
            }
        }, {
            "0": tree
        }, model_template, initial_values, run_options)
    }

    components = utility.Array1D(data_filter);
	
    filters = utility.Map({
        components,
        1
    }["_MATRIX_ELEMENT_ROW_"], "_value_", "''+ '`context`.nuc_data_' + _value_");

    lf_components = {
        2 * components,
        1
    };

    for (i = 0; i < components; i += 1) {
        lf_components[2 * i] = filters[i];
        DataSetFilter ^ (filters[i]) = CreateFilter( ^ (data_filter[i]), 1);
    }

    user_model_id = context + ".user_model";
    utility.ExecuteInGlobalNamespace ("`user_model_id` = 0");
	^(user_model_id) = model.generic.DefineModel(model_template, context + ".model", {
            "0": "terms.global"
        }, filters, None);

	// run over the empirically estimated frequencies with the initial values
	// TO DO: keep track on hyphy features to see if I can set custom frequencies without running over empirical ones
	if (run_options["set_freq"] == TRUE) {
		intial_frequencies = (initial_values[utility.getGlobalValue("terms.efv_estimate")])[utility.getGlobalValue("terms.model")];
		states_num = Columns((^(user_model_id))[utility.getGlobalValue("terms.alphabet")]);
		for (i=0; i<states_num; i+=1) {
			((^(user_model_id))[utility.getGlobalValue("terms.efv_estimate")])[i][0] = intial_frequencies[i][0];
		}
	}
	
	for (i = 0; i < components; i += 1) {
        lf_components[2 * i + 1] = "`context`.tree_" + i;
        custom_functions_bl_opt.ApplyModelToTree_FixedBLs(lf_components[2 * i + 1], tree[i], {
            "default": ^(user_model_id)
        }, None);
    }
	
    lfid = context + ".likelihoodFunction";
    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `lfid` = (`&lf_components`)");
	
	df = 0;
    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
            df = estimators.ApplyExistingEstimates(lfid, {
                (^user_model_id)[utility.getGlobalValue ("terms.id")]: ^(user_model_id)
            }, initial_values, run_options[utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler")]);
    }
	
	if (run_options["return_lf"] == TRUE) {
		return lfid;
	} else {
		return df;
	}
}




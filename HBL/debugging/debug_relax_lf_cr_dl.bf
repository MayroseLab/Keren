RequireVersion("2.3.3");

LoadFunctionLibrary("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/libv3/all-terms.bf"); // must be loaded before CF3x4

// namespace 'utility' for convenience functions
LoadFunctionLibrary("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/libv3/UtilityFunctions.bf");

// namespace 'io' for interactive/datamonkey i/o functions
LoadFunctionLibrary("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/libv3/IOFunctions.bf");

LoadFunctionLibrary ("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/libv3/models/codon/MG_REV.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/libv3/tasks/estimators.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/libv3/tasks/alignments.bf");

// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/libv3/tasks/trees.bf");

LoadFunctionLibrary("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/modules/selection_lib.ibf");
LoadFunctionLibrary("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/libv3/models/codon/BS_REL.bf");
LoadFunctionLibrary("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/libv3/convenience/math.bf");


utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable ("ASSUME_REVERSIBLE_MODELS", TRUE);

/*------------------------------------------------------------------------------*/

fprintf(stdout, "checking initial memory status\n");
num = io.PromptUser(">Type a number to continue", 0, 1, 10000, FALSE);

relax.json    = { terms.json.input: {},
                  terms.json.fits : {},
                  terms.json.timers : {},
                  terms.json.test_results : {}
                  };

relax.relaxation_parameter        = "relax.K";
relax.rate_classes     = 3;

terms.relax.k          = "relaxation or intensification parameter";
terms.relax.k_range    = {
        terms.lower_bound: "0",
        terms.upper_bound: "50"
    };

relax.test_branches_name = "Test";
relax.reference_branches_name = "Reference";

/*------------------------------------------------------------------------------*/

namespace relax {
    LoadFunctionLibrary ("/groups/itay_mayrose/halabikeren/hyphy_orig/2.3.7/hyphy/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "relax", utility.getGlobalValue("terms.settings") : {utility.getGlobalValue("terms.settings.branch_selector") : "relax.select_branches"}});
}

/* now fit the two main models for RELAX */

relax.test.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "relax.test", {
            "0": parameters.Quote(terms.global),
            "1": relax.codon_data_info[terms.code],
            "2": parameters.Quote (relax.rate_classes) // the number of rate classes
        },
        relax.filter_names,
        None);



relax.reference.bsrel_model =  model.generic.DefineMixtureModel("models.codon.BS_REL.ModelDescription",
        "relax.reference", {
            "0": parameters.Quote(terms.global),
            "1": relax.codon_data_info[terms.code],
            "2": parameters.Quote (relax.rate_classes) // the number of rate classes
        },
        relax.filter_names,
        None);

relax.bound_weights = models.BindGlobalParameters ({"0" : relax.reference.bsrel_model, "1" : relax.test.bsrel_model}, terms.mixture.mixture_aux_weight + ".+");
models.BindGlobalParameters ({"0" : relax.test.bsrel_model, "1" : relax.reference.bsrel_model}, terms.nucleotideRate("[ACGT]","[ACGT]"));

parameters.DeclareGlobalWithRanges (relax.relaxation_parameter, 1, 0, 50);
model.generic.AddGlobal (relax.test.bsrel_model, relax.relaxation_parameter, terms.relax.k);

for (relax.i = 1; relax.i < relax.rate_classes; relax.i += 1) {
    parameters.SetRange (model.generic.GetGlobalParameter (relax.reference.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)), terms.range01);
    parameters.SetRange (model.generic.GetGlobalParameter (relax.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)), terms.range01);
    parameters.SetConstraint (model.generic.GetGlobalParameter (relax.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)),
                              model.generic.GetGlobalParameter (relax.reference.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)) + "^" + relax.relaxation_parameter,
                              terms.global);
}
parameters.SetRange (model.generic.GetGlobalParameter (relax.reference.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.rate_classes)), terms.range_gte1);
parameters.SetRange (model.generic.GetGlobalParameter (relax.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.rate_classes)), terms.range_gte1);
parameters.SetConstraint (model.generic.GetGlobalParameter (relax.test.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)),
                          model.generic.GetGlobalParameter (relax.reference.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,relax.i)) + "^" + relax.relaxation_parameter,
                          terms.global);

relax.model_map = {
                    "relax.test" : utility.Filter (relax.selected_branches[0], '_value_', '_value_ == relax.test_branches_name'),
                    "relax.reference" : utility.Filter (relax.selected_branches[0], '_value_', '_value_ == relax.reference_branches_name')
                  };


// constrain the proportions to be the same

relax.model_object_map = { "relax.reference" : relax.reference.bsrel_model,
                            "relax.test" :       relax.test.bsrel_model };

relax.alternative_model.fit = CreateComputeLF (relax.filter_names, relax.trees, { "0" : relax.model_map}, None, relax.model_object_map, {terms.run_options.retain_lf_object: FALSE});



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

/**
 * based on the function estimators.FitLF
 * modifications can be seen in lines 207-258
 * exaplanation: switched the call to Optimize() to 3 calls equivalent to the content of estimators.ComputeLF() 
 * the purpose was to monitor the memory uptake of each step in the computation of the likelihood function
 * the memory update was tracked using the command: free -m | grep Mem | awk '{print $3}' (retuned the used memory in mega units in a node) via the same node in which the HBL script is run
 */
lfunction CreateComputeLF(data_filter, tree, model_map, initial_values, model_objects, run_options) {

    if (Type(data_filter) == "String") {
        return estimators.FitLF ({
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


    lf_components = {
        2 * components,
        1
    };


    for (i = 0; i < components; i += 1) {
        lf_components[2 * i] = data_filter[i];
        lf_components[2 * i + 1] = &tree_id + "_" + i;
        model.ApplyModelToTree(lf_components[2*i + 1], tree[i], None, model_map[i]);
    }

	fprintf(stdout, "checking memory status before creation of a likelihood function\n");
	num = io.PromptUser(">Type 1 to continue", 0, 1, 10000, FALSE);

    lf_id = &likelihoodFunction;
    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `lf_id` = (`&lf_components`)");
	
	fprintf(stdout, "checking memory status after creation of the likelihood function\n");
	num = io.PromptUser(">Type 1 to continue", 0, 1, 10000, FALSE);

    df = 0;

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
            df = estimators.ApplyExistingEstimates("`&likelihoodFunction`", model_objects, initial_values, run_options[utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler")]);
    }

    if (utility.Has (run_options,utility.getGlobalValue("terms.run_options.apply_user_constraints"),"String")) {
        df += Call (run_options[utility.getGlobalValue("terms.run_options.apply_user_constraints")], lf_id, lf_components, data_filter, tree, model_map, initial_values, model_objects);
    }
	
	fprintf(stdout, "checking memory status before computation of the likelihood function\n");
	num = io.PromptUser(">Type 1 to continue", 0, 1, 10000, FALSE);

   	// copy of content of estimators.ComputeLF() - used in order to track the memory consumption step by step
	/* according to HBL documentation from SKP's HTML:
	1. `receptacle`: 
	-- LF_START_COMPUTE - must be called before any computations are done to set up the internals, 
	-- is the identifier of the variable which receives the result, 
	-- LF_DONE_COMPUTE - must be called after all the computations are done to clean up the internals. 2. `likelihood_funtion_id` -- the likelihood function to be evaluated.
	*/
	LFCompute (likelihoodFunction,LF_START_COMPUTE);
	
	fprintf(stdout, "checking memory status after estimators.ComputeLF() - step 1 - set up the internals \n");
	num = io.PromptUser(">Type 1 to continue", 0, 1, 10000, FALSE);
	
	LFCompute (likelihoodFunction,logl);
	
	fprintf(stdout, "checking memory status after estimators.ComputeLF() - step 2 - compute logl \n");
	num = io.PromptUser(">Type 1 to continue", 0, 1, 10000, FALSE);
	
	LFCompute (likelihoodFunction,LF_DONE_COMPUTE);
	
	fprintf(stdout, "checking memory status after estimators.ComputeLF() - step 3 - clean up the internals \n");
	num = io.PromptUser(">Type 1 to continue", 0, 1, 10000, FALSE);

    results = {};

    results[utility.getGlobalValue ("terms.fit.log_likelihood")] = logl;

    results[utility.getGlobalValue ("terms.fit.filters")] = {
        1,
        components
    };

    for (i = 0; i < components; i += 1) {
        (results[utility.getGlobalValue ("terms.fit.filters")])[i] = lf_components[2 * i];

    }
	
	fprintf(stdout, "checking memory status before deletion of the likelihood function\n");
	num = io.PromptUser(">Type 1 to continue", 0, 1, 10000, FALSE);

    if (run_options[utility.getGlobalValue("terms.run_options.retain_lf_object")]) {
        results[utility.getGlobalValue("terms.likelihood_function")] = & likelihoodFunction;
    } else {
        DeleteObject(likelihoodFunction);
    }
	
	fprintf(stdout, "checking memory status after deletion of the likelihood function\n");
	num = io.PromptUser(">Type 1 to continue", 0, 1, 10000, FALSE);

    return results;
}

//------------------------------------------------------------------------------

lfunction relax.extract.k(branch_info) {
    return (branch_info[utility.getGlobalValue("terms.relax.k")])[utility.getGlobalValue("terms.fit.MLE")];
}

//------------------------------------------------------------------------------

lfunction relax.set.k (tree_name, node_name, model_description) {
    if (utility.Has (model_description [utility.getGlobalValue ("terms.local")], utility.getGlobalValue ("terms.relax.k"), "String")) {
        k = (model_description [utility.getGlobalValue ("terms.local")])[utility.getGlobalValue ("terms.relax.k")];
        parameters.SetValue (tree_name + "." + node_name + "." + k, 1);
        parameters.SetRange (tree_name + "." + node_name + "." + k, utility.getGlobalValue ("terms.relax.k_range"));
    }
    return tree_name + "." + node_name + "." + k;
}

//------------------------------------------------------------------------------

lfunction relax.init.k (lf_id, components, data_filter, tree, model_map, initial_values, model_objects) {
    parameter_set = estimators.TraverseLocalParameters (lf_id, model_objects, "relax.set.k");
    parameters.SetConstraint (model.generic.GetGlobalParameter (utility.getGlobalValue("relax.ge.bsrel_model") , terms.AddCategory (utility.getGlobalValue("terms.parameters.omega_ratio"),2)), utility.getGlobalValue("terms.parameters.one"), utility.getGlobalValue("terms.global"));
    /*parameters.SetConstraint (model.generic.GetGlobalParameter (utility.getGlobalValue("relax.ge.bsrel_model") , terms.AddCategory (utility.getGlobalValue("terms.parameters.omega_ratio"),utility.getGlobalValue ("relax.rate_classes"))),
                             "1/(" +
                                Join ("*", utility.Map (
                                    utility.Range (utility.getGlobalValue ("relax.rate_classes") - 1, 1, 1),
                                    "_value_",
                                    'model.generic.GetGlobalParameter (relax.ge.bsrel_model , terms.AddCategory (terms.parameters.omega_ratio,_value_))'
                                    ))
                             + ")",
                            "global");*/

    return 0;
}

//------------------------------------------------------------------------------

lfunction relax.BS_REL.ModelDescription (type, code, components) {
    model = models.codon.BS_REL.ModelDescription(utility.getGlobalValue ('terms.global'), code, components);
    model [utility.getGlobalValue("terms.model.defineQ")] = "relax.BS_REL._DefineQ";
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
                 (p[utility.getGlobalValue("terms.local")])[utility.getGlobalValue ("terms.relax.k")] = "k";
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

lfunction relax.BS_REL._DefineQ(bs_rel, namespace) {
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
       rate_matrices [key] = bs_rel[utility.getGlobalValue("terms.model.rate_matrix")];
       (bs_rel [^'terms.mixture.mixture_components'])[key] = _wts [component-1];
    }


    bs_rel[utility.getGlobalValue("terms.model.rate_matrix")] = rate_matrices;
    parameters.SetConstraint(((bs_rel[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.global")])[terms.nucleotideRate("A", "G")], "1", "");
    return bs_rel;
}

//------------------------------------------------------------------------------

lfunction relax.select_branches(partition_info) {

    io.CheckAssertion("utility.Array1D (`&partition_info`) == 1", "RELAX only works on a single partition dataset");
    available_models = {};
    branch_set = {};


    tree_for_analysis = (partition_info[0])[utility.getGlobalValue("terms.data.tree")];
    utility.ForEach (tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")], "_value_", "`&available_models`[_value_] += 1");
    list_models   = utility.Keys   (available_models); // get keys
    branch_counts = utility.Values (available_models);
    option_count  = Abs (available_models);

    io.CheckAssertion("`&option_count` >= 2", "RELAX requires at least one designated set of branches in the tree.");

    selectTheseForTesting = {
        option_count, 2
    };

    for (k = 0; k < option_count; k += 1) {
        if (list_models[k] != "") {
            selectTheseForTesting[k][0] = list_models[k];
            selectTheseForTesting[k][1] = "Set " + list_models[k] + " with " + available_models[list_models[k]] + " branches";
        } else {
            selectTheseForTesting[k][0] = "Unlabeled branches";
            selectTheseForTesting[k][1] = "Set of " + available_models[list_models[k]] + " unlabeled branches";
        }
    }

    ChoiceList(testSet, "Choose the set of branches to use as the _test_ set", 1, NO_SKIP, selectTheseForTesting);
    io.CheckAssertion ("`&testSet` >= 0", "User cancelled branch selection; analysis terminating");
    if (option_count > 2) {
        ChoiceList(referenceSet, "Choose the set of branches to use as the _reference_ set", 1, testSet, selectTheseForTesting);
        io.CheckAssertion ("`&referenceSet` >= 0", "User cancelled branch selection; analysis terminating");
    } else {
        referenceSet = 1-testSet;
    }

    return_set = {};

    tree_configuration = {};
    tree_for_analysis = (partition_info[0])[utility.getGlobalValue("terms.data.tree")];

    tag_test = selectTheseForTesting [testSet][0];
    if (tag_test == "Unlabeled branches") {
        tag_test = "";
    }
    tag_reference = selectTheseForTesting [referenceSet][0];
    if (tag_reference == "Unlabeled branches") {
        tag_reference = "";
    }

    utility.ForEachPair (tree_for_analysis[utility.getGlobalValue("terms.trees.model_map")], "_key_", "_value_", "
        if (`&tag_test` == _value_ ) {
            `&tree_configuration`[_key_] = utility.getGlobalValue('relax.test_branches_name');
        } else {
            if (`&tag_reference` == _value_ ) {
                `&tree_configuration`[_key_] = utility.getGlobalValue('relax.reference_branches_name');
            } else {
                `&tree_configuration`[_key_] = utility.getGlobalValue('relax.unclassified_branches_name');
            }
        }
    ");

    return_set + tree_configuration;
    return return_set;
}

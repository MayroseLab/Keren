LoadFunctionLibrary("libv3/models/model_functions.bf");
LoadFunctionLibrary("libv3/convenience/regexp.bf");
LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("custom_functions.bf"); // custom wrapper functions for TraitRELAX implementation

/**
 * Loads file information into a namespace specified by prefix
 * Sets the following variables:
 * * <prefix>.codon_data.mapping
 * * <prefix>.codon_data.sites
 * * <prefix>.codon_data.species
 * * <prefix>.codon_data.unique_sites
 * * <prefix>.codon_data_info
 * * <prefix>.codon_filter.sequence_map
 * * <prefix>.codon_filter.site_freqs
 * * <prefix>.codon_filter.site_map
 * * <prefix>.codon_filter.sites
 * * <prefix>.codon_filter.species
 * * <prefix>.filter_names
 * * <prefix>.filter_specification
 * * <prefix>.name_mapping
 * * <prefix>.partition_count
 * * <prefix>.partitions_and_trees
 * * <prefix>.prefix
 * * <prefix>.sample_size
 * * <prefix>.selected_branches
 * * <prefix>.trees
 * @param prefix {String} : The namespace to prefix all file information variables with
 * @return nothing, the function sets variables within a namespace
 */

utility.SetEnvVariable ("MARKDOWN_OUTPUT", TRUE);

/**
 * @name process_char_data
 * @param {String} data_path - path to the trait data path 
 * @usage loads the character alignment
 */
function process_char_data(prefix, alignment_path, partitions_and_trees) {

	// read the character alignment
	trait_info = alignments.ReadNucleotideDataSet (prefix+".trait_dataset", alignment_path);					// ASK STEPHANIE: why read as nuc data?
	trait_name_mapping = relative_prot_rates.alignment_info[utility.getGlobalValue("terms.data.name_mapping")];	// ASK STEPHANIE: why call from relative_prot_rates
	if (None == trait_name_mapping) {  
		trait_name_mapping = {};
		utility.ForEach (alignments.GetSequenceNames (prefix+".trait_dataset"), "_value_", "`&trait_name_mapping`[_value_ && 1] = _value_ && 1"); // 20.12.17: bypass of HyPhy issue 725 - converting species names to upper case letters regardless of the format of the input using '&& 1' on the dataset species names
	}
	// apply name mapping on the character alignment
	trait_filter_specification = alignments.DefineFiltersForPartitions (partitions_and_trees, "`prefix`.trait_dataset" , "`prefix`.trait_filter.", trait_info);
	trait_filter_names = utility.Map (trait_filter_specification, "_value_", "_value_[terms.data.name]");
	
	// set the tree to be used by the character model to be the user input tree
	trait_trees = utility.Map (partitions_and_trees, "_value_", "_value_[terms.data.tree]");
	
}


/**
 * @name load_input
 * @param {String} prefix 				- prefix to the instances generated within the function
 * @param {String} codon_alignment_path	- full path to the codon sequences alignment file (should be in nexus format)
 * @param {String} tree_path			- full path to the codon sequences alignment file (should be in newick format)
 * @usage loads the alignment and tree into an instance by the name codon_data_info with the given prefix 
 * Can only be used after including TraitRELAX_aux
 */
function load_input (run_info, json_path, codon_alignment_path, tree_path, trait_data_path) {
	
    settings = run_info[utility.getGlobalValue("terms.settings")];
    prefix = run_info [utility.getGlobalValue("terms.prefix")];
    
	// process codon alignment
	utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
	ExecuteCommands ('codon_data_info = alignments.PromptForGeneticCodeAndAlignment(prefix+".codon_data", prefix+".codon_filter");',
						 {"0" : "Universal", "1": codon_alignment_path});				
	utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);

    codon_data_info[utility.getGlobalValue("terms.data.sample_size")] = codon_data_info[utility.getGlobalValue("terms.data.sites")] * codon_data_info[utility.getGlobalValue("terms.data.sequences")];
    upper_prefix = prefix && 1; // uppercase the prefix for json name
    codon_data_info[utility.getGlobalValue("terms.json.json")] = json_path + "."+upper_prefix+".json";

    name_mapping = codon_data_info[utility.getGlobalValue("terms.data.name_mapping")];
    if (None == name_mapping) { /** create a 1-1 mapping if nothing was done */
        name_mapping = {};
        utility.ForEach (alignments.GetSequenceNames (prefix+".codon_data"), "_value_", "`&name_mapping`[_value_ && 1] = _value_ && 1"); // 20.12.17: bypass of HyPhy issue 725 - converting species names to upper case letters regarless of the format of the input using '&& 1' on the dataset species names
    }
	
	// process tree
	utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
    ExecuteCommands('partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions (codon_data_info[utility.getGlobalValue("terms.data.partitions")], name_mapping);', {"0": tree_path}); // this triggers the request for a tree
    utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);
	

	partition_count = Abs (partitions_and_trees);

	// TODO: DE-HARDCODE "filter-string"
    utility.ForEachPair (partitions_and_trees,
                            "_key_",
                            "_value_",
                            '(`&partitions_and_trees`[_key_])[utility.getGlobalValue("terms.data.filter_string")] = selection.io.adjust_partition_string (_value_[utility.getGlobalValue("terms.data.filter_string")], 3*`&codon_data_info`[utility.getGlobalValue("terms.data.sites")])');

    io.ReportProgressMessage ("", ">Loaded a multiple sequence alignment with **" + codon_data_info[utility.getGlobalValue("terms.data.sequences")] + "** sequences, **" + codon_data_info[utility.getGlobalValue("terms.data.sites")] + "** codons **\n");

	// process character data
	process_char_data(prefix, trait_data_path, partitions_and_trees);
	
	// set the branches partition
	if (utility.Has (settings, utility.getGlobalValue("terms.settings.branch_selector"), "String")) {
        selected_branches =  Call (settings[utility.getGlobalValue("terms.settings.branch_selector")], partitions_and_trees, "`prefix`.trait_dataset");
    } else {
        selected_branches = selection.io.defineBranchSets(partitions_and_trees);
    }
	
	// Place in own attribute called `tested`
    selection.io.json_store_key_value_pair (json, None, utility.getGlobalValue("terms.json.tested"), selected_branches);

    /** Input attribute to JSON **/
    json[utility.getGlobalValue("terms.json.input")] = {};
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.file")] =  codon_data_info[utility.getGlobalValue("terms.data.file")];
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.sequences")] = codon_data_info[utility.getGlobalValue("terms.data.sequences")];
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.sites")] = codon_data_info[utility.getGlobalValue("terms.data.sites")];
    (json[utility.getGlobalValue("terms.json.input")])[utility.getGlobalValue("terms.json.partition_count")] = partition_count;

    // The trees should go into input as well and they should be w/ their branch lengths but ONLY if they have any.
    t = (partitions_and_trees["0"])[utility.getGlobalValue("terms.data.tree")];
    abs_branch_lengths = Abs(t[utility.getGlobalValue("terms.branch_length")]);  
	
    if (abs_branch_lengths == 0){
        selection.io.json_store_key_value_pair (json,
                                             utility.getGlobalValue("terms.json.input"), utility.getGlobalValue("terms.json.trees"),
                                             utility.Map (partitions_and_trees, "_pt_", '(_pt_[terms.data.tree])[terms.trees.newick]')
                                             );    
    } else {
        selection.io.json_store_key_value_pair (json,
                                             utility.getGlobalValue("terms.json.input"), utility.getGlobalValue("terms.json.trees"),
                                             utility.Map (partitions_and_trees, "_pt_", '(_pt_[terms.data.tree])[terms.trees.newick_with_lengths]')
                                             ); 
    }

     filter_specification = alignments.DefineFiltersForPartitions (partitions_and_trees, "`prefix`.codon_data" , "`prefix`.filter.", codon_data_info);


     selection.io.json_store_key_value_pair (json, None, utility.getGlobalValue("terms.json.partitions"),
                                                         filter_specification);
     trees = utility.Map (partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');
	 // adjust the tree to fit the selected branches
	 //adjust_trees_to_partition(trees, selected_branches);
	 
     filter_names = utility.Map (filter_specification, "_partition_", '_partition_[terms.data.name]');

     /* Store original name mapping */
    for (partition_index = 0; partition_index < partition_count; partition_index += 1) {
        selection.io.json_store_branch_attribute(json, utility.getGlobalValue ("terms.original_name"), utility.getGlobalValue ("terms.json.node_label"), display_orders[utility.getGlobalValue ("terms.original_name")],
                                         partition_index,
                                         name_mapping);
    }

}


/**
 * @name adjust_trees_to_partition
 * @param {Dict} trees 			   - the dictionary holding the trees that need to be adjusted
 * @param {Dict} selected_branches - the dictionary holding the partition of the branches by which the tree needs to be adjusted
*/
lfunction adjust_trees_to_partition(trees, selected_branches) {
	
	for (i=0; i<Abs(trees); i+=1) {
		partition = selected_branches[i];
		(trees[i])[utility.getGlobalValue("terms.trees.model_list")] = {"", "Test"};
		annotated_tree = (trees[i])[utility.getGlobalValue("terms.trees.newick_annotated")];
		fprintf(stdout, "annotated_tree before edit: ", annotated_tree, "\n");
		branches = utility.Keys(partition);
		for (b=0; b<Columns(branches); b+=1) {
			branch = branches[b];
			branch_label = partition[branch];
			if (branch_label % "Test") {
				src = branch + ":";
				dest = branch + "{" + branch_label + "}:";
				annotated_tree = custom_functions.strReplace (annotated_tree, src, dest);
			}
		}
		fprintf(stdout, "annotated_tree after edit : ", annotated_tree, "\n"); 
		(trees[i])[utility.getGlobalValue("terms.trees.newick_annotated")] = annotated_tree;
	}
	fprintf(stdout, "adjusted trees: ", trees, "\n"); // debug
	exit();
}


/**
 * @name doGTR_FixedBLs
 * @param {String} prefix 				- prefix to the intances generated witihn the function
 * @usage fits the GTR model to a given alignment and tree instances, defined under the given prefix, while fixing the branch lengths
 * Can only be used after including TraitRELAX_aux
 */
function doGTR_FixedBLs (prefix) {

    io.ReportProgressMessageMD (prefix, "nuc-fit", "Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model");

    gtr_results = utility.Extend (parameters.helper.tree_lengths_to_initial_values (trees, None),
                                  {
                                    utility.getGlobalValue ("terms.global") : {
                                        terms.nucleotideRate ("A","C") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("A","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("C","G") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("G","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25}
                                    }
                                 });

	gtr_results = custom_functions.FitSingleModel_Ext_FixedBLs (filter_names, trees, "models.DNA.GTR.ModelDescription", gtr_results, {});

    io.ReportProgressMessageMD (prefix, "nuc-fit", "* " + selection.io.report_fit (gtr_results, 0, 3*(^"`prefix`.codon_data_info")[utility.getGlobalValue ("terms.data.sample_size")]));



    /* Store nucleotide fit */
    gtr_rates = utility.Map( 
                    utility.Map (gtr_results[utility.getGlobalValue("terms.global")], "_value_", '   {terms.fit.MLE : _value_[terms.fit.MLE]}'),
                    "_value_",
                    "_value_[terms.fit.MLE]");
    efv = (gtr_results[utility.getGlobalValue("terms.efv_estimate")])["VALUEINDEXORDER"][0];
    selection.io.json_store_lf_GTR_MG94 (json,
                                utility.getGlobalValue ("terms.json.nucleotide_gtr"),
                                gtr_results[utility.getGlobalValue ("terms.fit.log_likelihood")],
                                gtr_results[utility.getGlobalValue ("terms.parameters")] ,
                                codon_data_info[utility.getGlobalValue ("terms.data.sample_size")],
                                gtr_rates, 
                                efv,
                                display_orders[utility.getGlobalValue ("terms.json.nucleotide_gtr")]);

    
    /* TODO: Why does this not work here? */
    /*
    utility.ForEachPair (filter_specification, "_key_", "_value_",
        'selection.io.json_store_branch_attribute(json, utility.getGlobalValue ("terms.json.nucleotide_gtr"), utility.getGlobalValue ("terms.branch_length"), display_orders[terms.json.nucleotide_gtr],
                                         _key_,
                                         selection.io.extract_branch_info((gtr_results[utility.getGlobalValue ("terms.branch_length")])[_key_], "selection.io.branch.length"));');
    */
    
    /* Store branch lengths */
    for (partition_index = 0; partition_index < Abs(filter_specification); partition_index += 1) {
        selection.io.json_store_branch_attribute(json, utility.getGlobalValue ("terms.json.nucleotide_gtr"), utility.getGlobalValue ("terms.branch_length"), display_orders[terms.json.nucleotide_gtr],
                                         partition_index,
                                         selection.io.extract_branch_info((gtr_results[utility.getGlobalValue ("terms.branch_length")])[partition_index], "selection.io.branch.length"));
        }

}


/**
 * @name doPartitionedMG_FixedBLs
 * @param {String} prefix 				- prefix to the intances generated witihn the function
 * @param {Bool}   keep_lf				- boolean indicating in the generated likelihood function instance should be retained or deleted
 * Can only be used after including TraitRELAX_aux
 */
function doPartitionedMG_FixedBLs (prefix, keep_lf) {
    io.ReportProgressMessageMD ("`prefix`", "codon-fit", "Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases");
    /**
        Declare per-partition branch length scalers
        slac.scaler_prefix_0
        slac.scaler_prefix_1
        etc
    */
    scaler_variables = utility.PopulateDict (0, partition_count, "`prefix`.scaler_prefix + '_' + _k_", "_k_");
    utility.ForEach (scaler_variables, "_value_", "parameters.DeclareGlobal(_value_, None);parameters.SetValue(_value_, 3);");

    partitioned_mg_results = custom_functions.FitMGREV_FixedBLs (filter_names, trees, codon_data_info [utility.getGlobalValue("terms.code")], {
        utility.getGlobalValue("terms.run_options.model_type"): utility.getGlobalValue("terms.local"), 
        utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler"): scaler_variables,
        utility.getGlobalValue("terms.run_options.partitioned_omega"): selected_branches,
        utility.getGlobalValue("terms.run_options.retain_lf_object"): keep_lf
    }, gtr_results);

    io.ReportProgressMessageMD("`prefix`", "codon-fit", "* " + selection.io.report_fit (partitioned_mg_results, 0, (^"`prefix`.codon_data_info")[utility.getGlobalValue ("terms.data.sample_size")]));
    global_dnds = selection.io.extract_global_MLE_re (partitioned_mg_results, "^" + utility.getGlobalValue("terms.parameters.omega_ratio"));
    utility.ForEach (global_dnds, "_value_", 'io.ReportProgressMessageMD ("`prefix`", "codon-fit", "* " + _value_[utility.getGlobalValue("terms.description")] + " = " + Format (_value_[utility.getGlobalValue("terms.fit.MLE")],8,4));');

    /** extract and report dN/dS estimates */
}


/* ______________________________________ TEST TO CHECK IF ANY CHANGE OF THE ORGINAL OPTIMIZATION FUNCTION DRASTICALLY BLUW UP MEMORY USAGE _______________________________________*/



function doGTR (prefix) {

    io.ReportProgressMessageMD (prefix, "nuc-fit", "Obtaining branch lengths and nucleotide substitution biases under the nucleotide GTR model");

    gtr_results = utility.Extend (parameters.helper.tree_lengths_to_initial_values (trees, None),
                                  {
                                    utility.getGlobalValue ("terms.global") : {
                                        terms.nucleotideRate ("A","C") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("A","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("C","G") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25} ,
                                        terms.nucleotideRate ("G","T") : { utility.getGlobalValue ("terms.fit.MLE") : 0.25}
                                    }
                                 });

	gtr_results = estimators.FitGTR(filter_names,
                                         trees,
                                         gtr_results);

    io.ReportProgressMessageMD (prefix, "nuc-fit", "* " + selection.io.report_fit (gtr_results, 0, 3*(^"`prefix`.codon_data_info")[utility.getGlobalValue ("terms.data.sample_size")]));



    /* Store nucleotide fit */
    gtr_rates = utility.Map( 
                    utility.Map (gtr_results[utility.getGlobalValue("terms.global")], "_value_", '   {terms.fit.MLE : _value_[terms.fit.MLE]}'),
                    "_value_",
                    "_value_[terms.fit.MLE]");
    efv = (gtr_results[utility.getGlobalValue("terms.efv_estimate")])["VALUEINDEXORDER"][0];
    selection.io.json_store_lf_GTR_MG94 (json,
                                utility.getGlobalValue ("terms.json.nucleotide_gtr"),
                                gtr_results[utility.getGlobalValue ("terms.fit.log_likelihood")],
                                gtr_results[utility.getGlobalValue ("terms.parameters")] ,
                                codon_data_info[utility.getGlobalValue ("terms.data.sample_size")],
                                gtr_rates, 
                                efv,
                                display_orders[utility.getGlobalValue ("terms.json.nucleotide_gtr")]);

    
    /* TODO: Why does this not work here? */
    /*
    utility.ForEachPair (filter_specification, "_key_", "_value_",
        'selection.io.json_store_branch_attribute(json, utility.getGlobalValue ("terms.json.nucleotide_gtr"), utility.getGlobalValue ("terms.branch_length"), display_orders[terms.json.nucleotide_gtr],
                                         _key_,
                                         selection.io.extract_branch_info((gtr_results[utility.getGlobalValue ("terms.branch_length")])[_key_], "selection.io.branch.length"));');
    */
    
    /* Store branch lengths */
    for (partition_index = 0; partition_index < Abs(filter_specification); partition_index += 1) {
        selection.io.json_store_branch_attribute(json, utility.getGlobalValue ("terms.json.nucleotide_gtr"), utility.getGlobalValue ("terms.branch_length"), display_orders[terms.json.nucleotide_gtr],
                                         partition_index,
                                         selection.io.extract_branch_info((gtr_results[utility.getGlobalValue ("terms.branch_length")])[partition_index], "selection.io.branch.length"));
        }

}


/**
 * @name doPartitionMG
 * Can only be used after including shared-load-file
 * @return
 */
function doPartitionedMG (prefix, keep_lf) {
    io.ReportProgressMessageMD ("`prefix`", "codon-fit", "Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases");


    /**
        Declare per-partition branch length scalers
        slac.scaler_prefix_0
        slac.scaler_prefix_1
        etc
    */
    scaler_variables = utility.PopulateDict (0, partition_count, "`prefix`.scaler_prefix + '_' + _k_", "_k_");

    utility.ForEach (scaler_variables, "_value_", "parameters.DeclareGlobal(_value_, None);parameters.SetValue(_value_, 3);");



    partitioned_mg_results = estimators.FitMGREV(filter_names, trees, codon_data_info [utility.getGlobalValue("terms.code")], {
        utility.getGlobalValue("terms.run_options.model_type"): utility.getGlobalValue("terms.local"), 
        utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler"): scaler_variables,
        utility.getGlobalValue("terms.run_options.partitioned_omega"): selected_branches,
        utility.getGlobalValue("terms.run_options.retain_lf_object"): keep_lf
    }, gtr_results);


    io.ReportProgressMessageMD("`prefix`", "codon-fit", "* " + selection.io.report_fit (partitioned_mg_results, 0, (^"`prefix`.codon_data_info")[utility.getGlobalValue ("terms.data.sample_size")]));
    global_dnds = selection.io.extract_global_MLE_re (partitioned_mg_results, "^" + utility.getGlobalValue("terms.parameters.omega_ratio"));
    utility.ForEach (global_dnds, "_value_", 'io.ReportProgressMessageMD ("`prefix`", "codon-fit", "* " + _value_[utility.getGlobalValue("terms.description")] + " = " + Format (_value_[utility.getGlobalValue("terms.fit.MLE")],8,4));');

    /** extract and report dN/dS estimates */
}



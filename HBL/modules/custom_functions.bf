LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/models/model_functions.bf");
LoadFunctionLibrary("libv3/convenience/regexp.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/models/parameters.bf");
LoadFunctionLibrary ("libv3/models/frequencies.bf");
LoadFunctionLibrary ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");

/* _____________ MY OWN AUXILIARY FUNCTIONS ____________ */


/* ____________ GENERAL AUXILIARY FUNCTIONS ____________ */

/**
 * @name custom_functions.strReplace
 * @param {String} string - the string in which the regex substitution is required to occur
 * @param {String} src	  - the sub-string that needs to be replaced
 * @param {String} dst	  - the sub-string to replace src with 
 * @returns converted string
 */
lfunction custom_functions.strReplace (string, src, dst) 
{
	replacor = {2,1};
	replacor[0] = src;
	replacor[1] = dst;
	string = string^replacor;
	return string;
}
 
 
/** 
 * @name custom_functions.compute_tree_size
 * @param {Dict} tree
 * @return the total branches lengths of the provided tree
 */
lfunction custom_functions.compute_tree_size(tree) 
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

 
 /**
 * @name custom_functions.dup_data_filter
 * @param {String} data_pointer	  - pointer to the data instance for which a data filter needs to be created
 * @param {String} filter_pointer - pointer that will hold the created data filter
 */
lfunction custom_functions.dup_data_filter(datafilter_name, data_pointer) 
{
	utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null"); 
	ExecuteCommands('ExecuteAFile(HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "chooseGeneticCode.def");', {"0" : "Universal"});
	utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);
	DataSetFilter ^ datafilter_name = CreateFilter (^ data_pointer, 3, , , GeneticCodeExclusions);
	return	filter_name; 
}



/* ______________ BRENT AUXILIARY FUNCTIONS ____________ */

/** 
 * @name custom_functions.sign
 * @param {Float} a
 * @param {Float} b
 * @return the magnitude of a times the sign of b 
 */
lfunction custom_functions.sign(a,b)
 {		
	if (b < 0) {
		return -a;
	}
	return a;
}


/** 
 * @name custom_functions.shift
 * @param {Integer} var1 - pointer to variable holding a float
 * @param {Integer} var2 - pointer to variable holding a float
 * @param {Integer} var3 - pointer to variable holding a float
 * @param {Integer} var4 - pointer to variable holding a float
 * @action converts: var1 -> var2, var2 -> var3, var3 -> var4
 */
lfunction custom_functions.shift(var1, var2, var3, var4) 
{
	^ var1 = ^ var2;
	^ var2 = ^ var3;
	^ var3 = ^ var4;	
}


/** 
 * @name custom_functions.get_args_str
 * @param {Float} varying_arg 		  - pointer to variable holding a float
 * @param {Integer} varying_arg_index - pointer to variable holding a float
 * @param {Dict} const_args_dict	  - pointer to variable holding a float
 * @return {String} args_str 		  - an arguments string representing the given dictionary + the varying argument
 */
lfunction custom_functions.get_args_str(varying_arg, varying_arg_index, const_args_dict) 
{
	args_str = "";
	args_num = Abs(const_args_dict) + 1;
	for (i=0; i<args_num; i+=1) {
		if (i > 0) {
			args_str = args_str + ", " ;
		}
		if (i+1 == varying_arg_index) {
			args_str = args_str + varying_arg;
		} else {
			args_str = args_str + const_args_dict[i+1];
		}

	}
	return args_str
}


/** 
 * @name custom_functions.prepare_args
 * @param {Float} varying_arg 		  - pointer to variable holding a float
 * @param {Integer} varying_arg_index - pointer to variable holding a float
 * @param {Dict} args_dict			  - pointer to variable holding a float
 */
lfunction custom_functions.prepare_args(varying_arg, varying_arg_index, args_dict) 
{
	(^args_dict)[varying_arg_index] = varying_arg;
}


/** 
 * @name custom_functions.brent_optimize
 * @param {Float} x_lower 		 - the left (lower) bound of values x can accept
 * @param {Float} x_start 	     - the initial value of x (must be between x_lower and x_upper)
 * @param {Float} x_upper 		 - the right (upper) bound of values x can accept
 * @param {String} func_id 		 - name of the function that needs to be optimized
 * @param {Integer} x_arg_index  - the number of argument in which x should be placed when calling func_id
 * @param {Dict} func_args		 - dictionary that holds the rest of arguments for func_id, which remain constant (key is the argument index and value is the argument value).
 * @param {Float} min_dist 		 - the minimal allowed distance between two assignments of x on which the func_id is evaluated
 * @param {Integer} x_max 		 - a pointer holding the value of x that maximizes func_id 
 * @action performs brent optimization on function ^func_id with regards to a given parameter x
 */
lfunction custom_functions.brent_optimize(x_lower, x_start, x_upper, func_id, x_arg_index, func_args,  min_dist, x_max) 
{	
	IterationsLimit = 100;
	almostZero = 0.1;
	goldenRatio = 0.3819660;
	
	final_dist=0.0;																	  // e - the distance moved on the step before last
	a = Min(x_lower, x_upper);														  // reset a and b according to the input
	b = Max(x_lower, x_upper);								
																					  // x_start is still between a and b because x_start is between x_lower and x_upper.
	x=x_start;																	 	  // reset x, w, and v (1st minimal = 2nd minimal = 3rd minimal = bx)
	w = x;
	v = x;
	custom_functions.prepare_args(x, x_arg_index, func_args); // convert the arguments dictionary to a string of arguments to run the function on
	fx_id = &fx;
	utility.ExecuteInGlobalNamespace ("`fx_id` = " + func_id + "(`&func_args`);");   // compute the likelihood function for the reset point
	fw = fx;
	fv = fx;
	
	for (iter=1;iter<=IterationsLimit;iter+=1) { 									  // as long as the number of iteration has not exceeded limit IterationsLimit -> continue
		
		xm=0.5*(a+b);																  // set xm to be the middle between a and b
		min_dist_1=min_dist*Abs(x)+almostZero;
		min_dist_2=2.0*min_dist_1;
		if (Abs(x-xm) <= (min_dist_2-0.5*(b-a))) {									  // check if the local maximum is reached, and if so, return it
			^ x_max = x;
			return fx;
		}
		if (Abs(e) > min_dist_1) { 												      // this part need to be modified for getting the max - construct a trial parabolic fit ?? see mnbrak() in page 400 (424 in the pdf)
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=Abs(q);
			temp_dist=final_dist;
			final_dist=d;
			if (Abs(p) >= Abs(0.5*q*temp_dist) || p <= q*(a-x) || p >= q*(b-x))	{ 	  // the above conditions determine the acceptability of the parabolic fit. Here we take the golden section step into the larger of the two segments.
				if (x >= xm) {
					final_dist = a-x;
				} else {
					final_dist = b-x;
				}
				d=goldenRatio*final_dist;
			}
			else {
				d=p/q;																 // take the parabolic step ?? see mnbrak() in page 400 (424 in the pdf)
				u=x+d;
				if (u-a < min_dist_2 || b-u < min_dist_2) {
					d=custom_functions.sign(min_dist_1,xm-x);
				}
			}
		} else {
			if (x >= xm) {
				final_dist = a - x;
			} else {
				final_dist = b - x;
			}
			d=goldenRatio*final_dist;
		}
		if (Abs(d) >= min_dist_1) { 
			u = x+d;
		} else {
			u = x+custom_functions.sign(min_dist_1,d);
		}
		custom_functions.prepare_args(u, x_arg_index, func_args); // evaluate the function with regards to u
		fu_id = &fu;
		utility.ExecuteInGlobalNamespace ("`fu_id` = " + func_id + "(`&func_args`);"); // 6.2.18: segmentation falut in the second exernal iteration of TraitRELAX
		if (fu >= fx) { 																   // if you improved in this step (i.e, the new computed value is lower than the last computed value) ->
			if (u >= x) {
				a=x;
			} else {
				b=x;																	   // reduce the search range to be (x,b) if u > x else (a,x)
			}
			custom_functions.shift(&v,&w,&x,&u);										   // update v(=w),w(=x) and x(=u) and their likelihood														
			custom_functions.shift(&fv,&fw,&fx,&fu);													
		} else {																		   // if you didn't improve (reduce function) in this step -> 
			if (u < x) {
				a=u;
			} else {
				b=u;																	  // reduce the search range to be (u,b) if x > u else (a,u)
			}
			if (fu >= fw || w == x) {													  // if f(3rd minimal) <= f(2nd minimal) or 2nd minimal == 1st minimal ->
				v=w;																	  // switch the 3rd minimal with the 2nd minimal
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu >= fv || v == x || v == w) {									  // if f(last computed) <= f(3rd minimal) or 3rd minimal == 1st minimal or 3rd minimal == 2nd minimal ->
				v=u;																	  // switch the 3rd minimal with the last computed
				fv=fu;
			}
		} 
	}
	fprintf(stdout, "Too many iterations in brent\n");
	^ x_max = x;
	return fx;
}


/* __________ STOCHASTIC MAPPING AUXILIARY FUNCTIONS __________ */

/**
 * @name custom_functions.sample_mutations_given_ancestrals_per_branch
 * @param {Dict} branches_division 		- dictionary to fill in the time duration under 0 and 1 of the branch
 * @param {qDict} branches_history  	- dictionary to fill in the history along the branch
 * @param {String} sonName		   		- name of the node at the bottom of the branch
 * @param {String} fatherName	   		- name of the node at the top of the branch 
 * @param {String} sonState	       		- character state of the son node 
 * @param {String} fatherState	   		- character state of the father node 
 * @param {Float} branch_length	   		- the length of the branch
 * @param {Float} pi_0	  		   		- frequency of state 0
 * @param {Float} mu	  		   		- character substitution rate
 * @param {Integer} maxSimulationsNum	- maximal number of allowed attempt to simulate a branch history 
 * @usage generates a history of a given branch and the characters at its edges
 *	      produces a dictionary that maps to each branch the durations under state 0 and 1	   
 */
lfunction custom_functions.sample_mutations_given_ancestrals_per_branch (branches_division, branches_history, sonName, fatherName, sonState, fatherState, branch_length, pi_0, mu, maxSimulationsNum)
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
			freq = pi_0;
			if (curState == 0) {
				freq = 1 - pi_0;
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
			freq = pi_0;
			if (curState == 0) {
				freq = 1 - pi_0;
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
 * @name custom_functions.generate_history
 * @param {Integer} history_num 		- number of the history that should be created (goes from 0 to num_of_histories)
 * @param {Dict} trees_duration  		- dictionary to fill in the durations under state 0 and state 1 for each branch in the tree
 * @param {Integer} lfId		   		- pointer of the character likelihood function
 * @param {Dict} tree		   			- dictionary holding the original user given tree
 * @param {Float} pi_0	  		   		- frequency of state 0
 * @param {Float} mu	  		   		- character substitution rate
 * @param {Integer} maxSimulationsNum	- maximal number of allowed attempt to simulate a branch history 
 * @return {Dict} historyTree			- a dictionary holding the sampled history in a tree format
 * @usage generates a history of a given branch and the characters at its edges
 *	      produces a dictionary that maps to each branch the durations under state 0 and 1	   
 */
lfunction custom_functions.generate_history (history_num, trees_duration, lfId, tree, pi_0, mu, maxSimulationsNum, run_options)
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
			custom_functions.sample_mutations_given_ancestrals_per_branch(&branches_division, &branchesHistory, sonName, fatherName, sonState, fatherState, branchLength, pi_0, mu, maxSimulationsNum);
		}
	}
	
	// convert the SM info to tree with branch classification
	// represent the history as a tree with internal nodes, that have a single child
	// the trait states per branch in the new tree will dictate the label of the branch
	historyTree = custom_functions.generate_tree(tree, &branchesHistory, run_options, history_num);
	(^ trees_duration)[history_num] = branches_division; 
	
	// delete the ancestral character data
	DeleteObject(ancestors);
	
	return historyTree;
}


/**
 * @name custom_functions.average_histories
 * @param {Dict} tree		   					- dictionary holding the original user given tree
 * @param {Dict} histories_duration_analysis  	- dictionary that holds the durations under state 0 and state 1 for each branch in the tree
 * @return {Dict} expected_history			    - a dictionary holding the expected (averaged) history in a tree format
 * @usage generates a history of a given branch and the characters at its edges
 *	      produces a dictionary that maps to each branch the durations under state 0 and 1	   
 */
lfunction custom_functions.average_histories(tree, histories_duration_analysis, run_options) 
{ 
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
		} else { // history is read by custom_functions.generate_tree top to bottom, so the last transition denotes the duration between the son and the first internal node
			(expected_history_info[branch_names[b]])["transitionsNum"] = 1;
			((expected_history_info[branch_names[b]])["history"])[0] = (expected_durations[branch_names[b]])[1-son_state];
			((expected_history_info[branch_names[b]])["history"])[1] = (expected_durations[branch_names[b]])[son_state];
		}
	}
	expected_history = custom_functions.generate_tree(tree, &expected_history_info, run_options, 1);
	return expected_history;
}


/**
 * @name custom_functions.generate_tree
 * @param {Dict} tree		   	   - dictionary holding the original user given tree
 * @param {Dict} branches_history  - dictionary that holds the sampled character history
 * @return {Dict} history	       - a dictionary holding the history in a tree format
 * @usage convert the history in its dictionary form to a tree based on the original tree   
 */
lfunction custom_functions.generate_tree (tree, branches_history, run_options, history_num)
{
	branch_names = trees.BranchNames(tree); // get the core branch names from the original tree
	history_tree = {};
	internalsCounter = 1;
	tree_str = tree[utility.getGlobalValue("terms.trees.newick_with_lengths")]; 
	Topology historyTree = tree_str; // use utility.getGlobalValue("terms.trees.newick")] if you don't need branch lengths
	label_to_branch = {};
	length_to_branch = {};
	
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
				length_to_branch[node] = (branch_history["history"])[tc]; // not the same as distances from the tree - need to figure out with ete3
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
			// change the length of the branch
			length_to_branch[branch_history["parent"]] = (branch_history["history"])[0];
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
		dst_bl = (node_history["history"])[transitions_num]; // here might need to use length_to_branch instead
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
		} else {
			dst_expr = node + "{Reference}:";
			historyTreeStr = custom_functions.strReplace(historyTreeStr, src_expr, dst_expr);
			model_map[node] = "Reference";	
		}
	}
	
	// debug - write the last character histories to files
	if (None != run_options) {
		if (run_options["write_histories"] == TRUE) {
			history_output_path = run_options["histories_dir"] + "history_" + history_num;
			if (run_options["average"] == TRUE) {
				history_output_path =  run_options["histories_dir"] + "average_history";
			}
			fprintf(history_output_path, historyTreeStr, ";");
		}
	}

	history = trees.ExtractTreeInfo(historyTreeStr);
	
	return history;
}
 
 
/* _________ MAXIMUM PARSIMONY AUXILIARY FUNCTIONS _________ */


/**
 * @name mapNodeToParent
 * @param {Dict} tree 				- dictionary representing the tree
 * @return {Dict} node_to_parent 	- a map of node to its parent along the tree
 */
 lfunction custom_functions.mapNodeToParent(tree) 
 {
	 // get a list of the tree nodes in post-order
	treeStr = tree[utility.getGlobalValue("terms.trees.newick_with_lengths")];
	Tree tree_instance = treeStr;
	nodes = trees.BranchNames(tree);
	tree_avl = tree_instance ^ 0;
	
	// generate a map of node to its parent
	node_to_parent = {};
	for (i=0; i<Columns(nodes)-1; i+=1) {
		node_to_parent[nodes[i]] = (tree_avl[(tree_avl[i+1])["Parent"]])["Name"];
	}
	
	return node_to_parent;	 
 }
 
 
 /**
 * @name mapNodeToChildren
 * @param {Dict} tree 				- dictionary representing the tree
 * @param {Dict} nodeto_parent		- map of node to its parent
 * @return {Dict} node_to_children 	- a map of node to its children along the tree
 */
 lfunction custom_functions.mapNodeToChildren(tree, node_to_parent) 
 {	 
	 // get a list of nodes
	 nodes = trees.BranchNames(tree);
	 
	 // generate a reverse map of node to its children (should contain only internal nodes)
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
	
	return node_to_children;

}

 
 /**
 * @name getMPScore
 * @param {Dict} tree 				- dictionary representing the tree
 * @param {String} data				- name of the corresponding dataset. Must consist of binary states
 * Calculates the maximum parsimony score of the tree given the data
 */
lfunction custom_functions.getMPScore (tree, data) 
{
	// get a list of the tree nodes
	nodes = trees.BranchNames(tree);
	
	// generate a map of node to its parent
	node_to_parent = custom_functions.mapNodeToParent(tree);
	
	// get a map of node to its children
	node_to_children = custom_functions.mapNodeToChildren(tree, node_to_parent);
	
	// set the state of each node (could either be a single state or a union of multiple states), while counting the number of intersection events
	// currently, since the character sequence is consisted of a single site, we don't iterate over sites. in order to generalize this function, we would have to do so
	node_to_state = {};
	intersection_events_num = 0;
	leafs = alignments.GetSequenceNames(data);
	for (i=0; i<Columns(leafs); i+=1) { // first get the statuses of the leafs
		node_to_state[leafs[i]] = (alignments.GetIthSequence(data, i))[utility.getGlobalValue("terms.data.sequence")]; 
		leaf_num += 1;
	}

	for (i=0; i<Columns(nodes); i+=1) { // how hanlde internal nodes
		is_leaf = (tree[utility.getGlobalValue("terms.trees.partitioned")])[nodes[i]] % "leaf";
		if (!is_leaf) {
			children_num = (node_to_children[nodes[i]])["childrenNum"];
			children_names = (node_to_children[nodes[i]])["childrenNames"];
			children_states = {};
			states_counter = 0; // in our case, can only go up to 2
			handled = False;
			for (j=0; j<children_num; j+=1) {
				child_name = children_names[j];
				child_state = node_to_state[child_name]; // since we traverse the tree in post order traversal, we are assured to visit the children before the parent
				is_both = child_state % "both";			 // check if the child already holds an intersection between the two states 
				if (is_both) {
					node_to_state[nodes[i]] = "both";
					handled = TRUE;
				} else {
					if (j == 0 && !is_both) {
						children_states[j] = child_state;
						states_counter += 1;
					} else {
						is_new = 1;
						for (k=0; k<states_counter; k+=1) {
							if (children_states[k] % child_state) {
								is_new = 0;
							}
						}
						if (is_new == 1 && !is_both) {
								children_states[states_counter] = child_state;
								states_counter += 1;
						}
					}
				}
			}
			if (!is_leaf && !handled) {
				if (states_counter > 1) {
					intersection_events_num += 1;
					node_to_state[nodes[i]] = "both";
				} else {
					node_to_state[nodes[i]] = child_state;
				}
			}
		}
	}
	
	// return the number of intersection events, which is equal to the MP score
	return {"0": intersection_events_num, "1": node_to_state};
}


/** 
 * @name getMPTreePartition
 * @param {Dict} tree 				- dictionary representing the tree
 * @param {Dict} state_to_label		- map of data state to label name
 * @param {String} data				- name of the corresponding dataset. Must consist of binary states
 * Creates a maximum parsimony based partition of the tree, based on the provided binary data
 */
 lfunction custom_functions.getMPTreePartition(tree, state_to_label, data) 
 { 
	 // get a list of possible states
	 states = utility.Keys(state_to_label);
	 
	 // get the state intersection and union events map
	 mp_data = custom_functions.getMPScore (tree, data);
	 mp_score = mp_data["0"];
	 node_to_state = mp_data["1"];
	 
     // get a list of the tree nodes in pre-order
	 treeStr = tree[utility.getGlobalValue("terms.trees.newick_with_lengths")];
	 Tree tree_instance = treeStr;
	 nodes = trees.BranchNames(tree);
	 tree_avl = tree_instance ^ 1; // root is names "Node0" and has no "Parent" field
	 for (i=0; i< Columns(nodes); i+=1) {
	 	hanlded_node = (tree_avl[i+1])["Name"];
	 	node_state = node_to_state[hanlded_node];
	 	if (i > 0) { // if node is not the root
	 		parent = (tree_avl[((tree_avl[i+1])["Parent"])])["Name"];
	 		parent_state = node_to_state[parent];
	 		if (node_state % "both") {
	 			node_to_state[hanlded_node] = parent_state;
	 		}
	 	} else { // the node is the root -> set the state arbitrarily if it can hold both
	 		if (node_state % "both") {
				if ((node % "Node0") == 0) { // if the node is not the root (the root shoudn't recieve any state in the partition since it doesn't represnt any branch)
					node_to_state[hanlded_node] = states[0];
				}
	 		}
	 	}
	 }
	
	 // convent the states mapping to a partition based on the state->label map
	 return_set = {};
	 for(i=0; i<Columns(states); i+=1) {
	 	state = states[i];
	 	for (j=0; j<Columns(nodes); j+=1) {
	 		node = nodes[j];
	 		node_state = node_to_state[node];
	 		if (node_state == state) {
	 			return_set[node] = state_to_label[state];
	 		}
	 	}
	 } 
	
	return {"0": mp_score, "1": return_set};
 } 


/* _________ CONSTRAIN BRANCH LENGTHS AUXILIARY FUNCTIONS _________ */


/**
 * @name custom_functions.fix_branch_lengths
 * @param {String} id
 * @param {Dictionary} model_info
 * @param {List} branch_names
 * @param {Dictionary} branch_lengths
 */
 // for improvement: is this implementation similar to the future HBL feature "constrain_branch_length"
lfunction custom_functions.fix_branch_lengths (id, model_info, branch_names, branch_lengths) 
{	
	branch_length_expr = model_info[utility.getGlobalValue("terms.model.branch_length_string")];

	for (bi=0; bi<Columns(branch_names); bi+=1) {   // ignore the last branch name that represents the root
		bl = branch_lengths[branch_names[bi]];
		local_params_dict = (model_info[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.local")];
		local_params_names = utility.Keys(local_params_dict);
		if (Columns(local_params_names) > 1) {
			num_to_constrain = Columns(local_params_names) - 1;
		} else {
			num_to_constrain = Columns(local_params_names);
		}
		for (i=0; i<num_to_constrain; i+=1) { // do not constrain the last parameter in the list as it would create a circular constraints issue
			branch_parameter_name = local_params_names[i];
			branch_parameter = ((model_info[utility.getGlobalValue("terms.parameters")])[utility.getGlobalValue("terms.local")])[branch_parameter_name];
			bl_expr_leftover = Simplify (branch_length_expr, {branch_parameter : 1});
			if (Columns(local_params_names) > 1) {
				for (j=0; j<Columns(local_params_names); j+=1) {
					parameter = local_params_dict[local_params_names[j]];
					bl_specific_parameter = "`id`." + branch_names[bi] + "." + parameter;
					bl_expr_leftover = custom_functions.strReplace(bl_expr_leftover, parameter, bl_specific_parameter); // the issue is probably here... look for slip of regex of parameter t that affects parameter k
				}
			}
			constraint := "`id`." + branch_names[bi] + "." + branch_parameter + ":=" + bl + "/(" + bl_expr_leftover + ");";
			utility.ExecuteInGlobalNamespace (constraint);
		}
	}
}


/**
 * @name custom_functions.clear_model_constraints
 * @param {Long} lf_id - pointer to a likelihood function
 */
lfunction custom_functions.clear_model_constraints(lf_id) 
{
	GetString(lf_info, ^ lf_id, -1);
	
	// clear local constraints
	local_constrained_parameters = lf_info[utility.getGlobalValue("terms.parameters.local_constrained")]; 
	for (i=0; i<Columns(local_constrained_parameters); i+=1) {
		utility.ExecuteInGlobalNamespace ("ClearConstraints(" + local_constrained_parameters[i] + ");");
	}

	// clear global constraints
	global_constrained_parameters = lf_info[utility.getGlobalValue("terms.parameters.global_constrained")]; 
	for (i=0; i<Columns(global_constrained_parameters); i+=1) {
		utility.ExecuteInGlobalNamespace ("ClearConstraints(" + global_constrained_parameters[i] + ");");
	}	
	
}


/**
 * @name custom_functions.DeleteLikelihoodFunction_FixedBLs
 * @param {Long} lfID - pointer to the a likelihood function with constrained branches that need to be cleared upon deletion
 */
 // for improvement: if the model and data filter wouldn't have been deleted each time a likelihood function is deleted, I wouldn't have to re-create them iteratively
lfunction custom_functions.DeleteLikelihoodFunction_FixedBLs(lfId) 
{
	// clear the constraints on the model parameters
	custom_functions.clear_model_constraints(lfId);
	
	// delete the likelihood function
	DeleteObject(^ lfId); 
}


/* ______________________ FUNCTIONS BASED ON FUNCTIONS FROM THE LIBV3 OF HYPHY ______________ */


/**
 * @name custom_functions.FitLF_FixedBLs
 * @param {Dict} data_filter  		- a vector of DataFilters
 * @param {Dict} tree  				- a vector of Trees
 * @param {String} lf_formula		- custom formula for the likelihood function computation
 * @param {Dict} model_map 			- map that holds the labeling of the trees branches into models 
 * @param {Dict} initial_values 	- dictionary that holds the initial values of the models parameters 
 * @param {Dict} model_objects 		- dictionary that holds the models instances
 * @param {Dict} run_options 		- dictionary that holds the function run settings 
 * @returns LF results after optimization
 */
lfunction custom_functions.FitLF_FixedBLs(data_filter, tree, lf_formula, model_map, initial_values, model_objects, run_options)
{
    if (Type(data_filter) == "String") {
        return custom_functions.FitLF_FixedBLs ({
            {
                data_filter__
            }
        }, {
            "0": tree
        }, lf_formula, 
        {
            "0" : model_map
        },
        initial_values, model_objects, run_options);
    }

    components = utility.Array1D(data_filter); 
	components_num = 2 * components;
	if (run_options["custom_lf_formula"]) { // add an additional slot from the likelihood function formula
		components_num += 1; 
	}
	lf_components = {
        components_num,
        1
    };
	
	/* for improvement: here, I need to chain for each tree the same data filter - but can't use a reference to the same instance of data filter 
	(else I get an error in the deletion of the likelihood function */
    for (i = 0; i < components; i += 1) {
        lf_components[2 * i] = data_filter[i];
        lf_components[2 * i + 1] = &tree_id + "_" + i;
		if (run_options["optimize_branch_lengths"]) {
			model.ApplyModelToTree (lf_components[2*i + 1], tree[i], None, model_map[i]);
		} else {
			custom_functions.ApplyModelToTree_FixedBLs (lf_components[2*i + 1], tree[i], model_objects, model_map[i]);
		}
    }
	
	if (run_options["custom_lf_formula"]) { // add the custom likelihood function formula as the last component of the input for the likelihood function constructor
		lf_components[components_num-1] = lf_formula;
	}

    lf_id = &likelihoodFunction;
    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `lf_id` = (`&lf_components`)");

    df = 0;

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
            df = estimators.ApplyExistingEstimates("`&likelihoodFunction`", model_objects, initial_values, run_options[utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler")]);
    }

    if (utility.Has (run_options,utility.getGlobalValue("terms.run_options.apply_user_constraints"),"String")) {
        df += Call (run_options[utility.getGlobalValue("terms.run_options.apply_user_constraints")], lf_id, lf_components, data_filter, tree, model_map, initial_values, model_objects);
    }
	
    //assert (0);

    //utility.SetEnvVariable ("VERBOSITY_LEVEL", 10);
   	Optimize (mles, likelihoodFunction);
	
	/* fprintf(stdout, "likelihoodFunction: ", likelihoodFunction, "\n"); // debug
	GetString(lf_info, likelihoodFunction, -1); // debug
	fprintf(stdout, "lf_info: ", lf_info, "\n\n"); // debug */
	
    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }

    results = estimators.ExtractMLEs( & likelihoodFunction, model_objects);

    results[utility.getGlobalValue ("terms.fit.log_likelihood")] = mles[1][0];
    results[utility.getGlobalValue ("terms.parameters")] = mles[1][1] + df;

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
        custom_functions.DeleteLikelihoodFunction_FixedBLs(lf_id);
    }

    return results;
}	


/**
 * @name custom_functions.FitLF_FixedBLs
 * @param {Dict} data_filter  		- a vector of DataFilters
 * @param {Dict} tree  				- a vector of Trees
 * @param {String} lf_formula		- custom formula for the likelihood function computation
 * @param {Dict} model_map 			- map that holds the labeling of the trees branches into models 
 * @param {Dict} initial_values 	- dictionary that holds the initial values of the models parameters 
 * @param {Dict} model_objects 		- dictionary that holds the models instances
 * @param {Dict} run_options 		- dictionary that holds the function run settings 
 * @returns LF results after computation
 */
lfunction custom_functions.ComputeLF_FixedBLs (data_filter, tree, lf_formula, model_map, initial_values, model_objects, run_options) 
{    
	if (Type(data_filter) == "String") {
        return custom_functions.FitLF_FixedBLs ({
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
        if (run_options["optimize_branch_lengths"]) {
			model.ApplyModelToTree (lf_components[2*i + 1], tree[i], None, model_map[i]);
		} else {
			custom_functions.ApplyModelToTree_FixedBLs (lf_components[2*i + 1], tree[i], model_objects, model_map[i]);
		}
    }
	
	if (run_options["custom_lf_formula"] == TRUE) {
		lf_components[2 * components] = lf_formula;
	}
		
    lf_id = &likelihoodFunction;
    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `lf_id` = (`&lf_components`)");
	
    df = 0;

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", 1);
        df = estimators.ApplyExistingEstimates("`&likelihoodFunction`", model_objects, initial_values, run_options[utility.getGlobalValue("terms.run_options.proportional_branch_length_scaler")]); 
    }

    if (utility.Has (run_options,utility.getGlobalValue("terms.run_options.apply_user_constraints"),"String")) {
        df += Call (run_options[utility.getGlobalValue("terms.run_options.apply_user_constraints")], lf_id, lf_components, data_filter, tree, model_map, initial_values, model_objects);
    }
	
	// for improvement: this is the only difference from FitLF: instead of optimizing, I am calling estimators.ComputeLF (which is very expensive in time)
   	log_likelihood = estimators.ComputeLF (& likelihoodFunction);
	
	//fprintf(stdout, "likelihood function after restricting parameters:\n", likelihoodFunction, "\n\n"); // debug
	
    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }
	
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

    if (run_options[utility.getGlobalValue("terms.run_options.retain_lf_object")]) {
        results[utility.getGlobalValue("terms.likelihood_function")] = & likelihoodFunction;
    } else {
        custom_functions.DeleteLikelihoodFunction_FixedBLs(lf_id);
    }
    return results;
}


/**
 * @name custom_functions.ApplyModelToTree_FixedBLs
 * @param id
 * @param tree
 * @param model_list
 * @param rules
 * @usage applies models on branches of the tree by the name 'id' (parameter id) according to a given mapping (parameter rules) and a list of models (parameter model_list). restricts branch lengths to be fixed if and only if the constraint is a linear combination with regards to the constrained parameter
 */
function custom_functions.ApplyModelToTree_FixedBLs (id, tree, model_list, rules) 
{	
	// extract the tree branch names and lengths for the procedure of constraining branch lengths 
	branch_names = trees.BranchNames(tree);
	branch_lengths = tree[utility.getGlobalValue("terms.branch_length")];

	if (Type (rules) == "AssociativeList") { // keren 31.1.18: always want to enter here when using trait relax
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
	    custom_functions.ApplyModelToTree_FixedBLs.ids = Rows (rules); // = {{"trait_relax.test", "trait_relax.reference"}}
	    for (custom_functions.ApplyModelToTree_FixedBLs.k = 0; custom_functions.ApplyModelToTree_FixedBLs.k < Abs (rules); custom_functions.ApplyModelToTree_FixedBLs.k += 1) {
	        custom_functions.ApplyModelToTree_FixedBLs.name = custom_functions.ApplyModelToTree_FixedBLs.ids[custom_functions.ApplyModelToTree_FixedBLs.k];
	        if ( custom_functions.ApplyModelToTree_FixedBLs.name != "DEFAULT") { 
                custom_functions.ApplyModelToTree_FixedBLs.list = rules[custom_functions.ApplyModelToTree_FixedBLs.name];
                if (Type (custom_functions.ApplyModelToTree_FixedBLs.list) == "AssociativeList") {
                    custom_functions.ApplyModelToTree_FixedBLs.list = Rows (custom_functions.ApplyModelToTree_FixedBLs.list);
                }

			   // for each model in the list, apply it on the relevant branches and fix the branch length under the corresponding model (by constraining the local parameter of the branch)
                for (custom_functions.ApplyModelToTree_FixedBLs.b = 0; custom_functions.ApplyModelToTree_FixedBLs.b < Columns (custom_functions.ApplyModelToTree_FixedBLs.list); custom_functions.ApplyModelToTree_FixedBLs.b += 1) {
                    // assign to the branch in place b in the list of branches under model custom_functions.ApplyModelToTree_FixedBLs.apply_model the model custom_functions.ApplyModelToTree_FixedBLs.apply_model
					ExecuteCommands ("SetParameter (`id`." + custom_functions.ApplyModelToTree_FixedBLs.list[custom_functions.ApplyModelToTree_FixedBLs.b] + ",MODEL," + custom_functions.ApplyModelToTree_FixedBLs.name + ")");
				}
				// for improvement: maybe my implementation of the fixation of branch lengths is not optimal
				custom_functions.fix_branch_lengths (id, model_list[custom_functions.ApplyModelToTree_FixedBLs.name], branch_names, branch_lengths);
            }
	    }
	// if there is a single model that applies on all the trees, simply apply it on the tree and constrain all the local parameters of the branches under the same model
	} else { // make sure not to reach here, as relax doesn't reach here either, and it's not working...
	    // TO DO: REMOVE HARDCODING
		custom_functions.ApplyModelToTree_FixedBLs.modelID = model_list[model_list ["INDEXORDER"][0]]; // ASK STEPHANIE: IS custom_functions.ApplyModelToTree_FixedBLs.modelID[terms.id] THE MODEL INSTANCE? BECAUSE ABOVE, IN SetParameter(), custom_functions.ApplyModelToTree_FixedBLs.apply_model WASN'T A PARAMETER
		ExecuteCommands ("UseModel (" + custom_functions.ApplyModelToTree_FixedBLs.modelID[terms.id] + ");
						  Tree `id` = " + tree["string"] + ";
						  ");
		// fix the branches lengths
		custom_functions.fix_branch_lengths (id, custom_functions.ApplyModelToTree_FixedBLs.modelID, branch_names, branch_lengths);
	}
}


/**
 * @name custom_functions.FitSingleModel_Ext_FixedBLs
 * @param {DataFilter} data_filter
 * @param {Tree} tree
 * @param {Dict} model
 * @param {Matrix} initial_values
 * @param {Dict} run_options
 * @returns results (unlike estimators.FitSingleModel_Ext, this function also fixes the branch lengths before optimization)
 */
lfunction custom_functions.FitSingleModel_Ext_FixedBLs (data_filter, tree, model_template, initial_values, run_options) 
{	
    this_namespace = (&_);
    this_namespace = this_namespace[0][Abs (this_namespace)-3];

    df = custom_functions.CreateLFObject_FixedBLs (this_namespace, data_filter, tree, model_template, initial_values, run_options);

	Optimize(mles, likelihoodFunction);
	
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
        custom_functions.DeleteLikelihoodFunction_FixedBLs(&likelihoodFunction);
    }

    return results;
}


/**
 * @name custom_functions.FitMGREV_FixedBLs
 * @param {DataFilter} codon_data
 * @param {Tree} tree
 * @param {String} genetic_code
 * @param {Dictionary} option
 * @param {Dictionary} initial_values
 * @returns MGREV results
 */
lfunction custom_functions.FitMGREV_FixedBLs (codon_data, tree, genetic_code, option, initial_values) 
{
    //TODO (TEAM): Where is data_filter being set?
    if (Type(data_filter) == "String") {
        return custom_functions.FitMGREV_FixedBLs({
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
        custom_functions.ApplyModelToTree_FixedBLs(Eval("&`lf_components[2*i + 1]`"), tree[i], model_assignment, None); 
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
        io.CheckAssertion("None!=`&alpha` && None!=`&beta`", "Could not find expected local synonymous and non-synonymous rate parameters in \`custom_functions.FitMGREV_FixedBLs\`");

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

    Optimize(mles, likelihoodFunction);
	
	// clear the constraints on the branch lengths
	custom_functions.clear_model_constraints(& likelihoodFunction); 

    if (Type(initial_values) == "AssociativeList") {
        utility.ToggleEnvVariable("USE_LAST_RESULTS", None);
    }

    results = estimators.ExtractMLEs( & likelihoodFunction, model_id_to_object);


    results[utility.getGlobalValue("terms.fit.log_likelihood")] = mles[1][0];
    results[utility.getGlobalValue("terms.parameters")] = mles[1][1] + (mg_rev [utility.getGlobalValue("terms.parameters")]) [utility.getGlobalValue("terms.model.empirical")] + df;


    if (option[utility.getGlobalValue("terms.run_options.retain_lf_object")]) {
        results[utility.getGlobalValue("terms.likelihood_function")] = & likelihoodFunction;
    } else {
        custom_functions.DeleteLikelihoodFunction_FixedBLs(&likelihoodFunction);
    }

    if (option[utility.getGlobalValue("terms.run_options.retain_model_object")]) {
        results[utility.getGlobalValue("terms.model")] = model_id_to_object;
    }

    return results;
}


/**
 * @name custom_functions.CreateLFObject_FixedBLs
 * @param {String} context
 * @param {DataFilter} data_filter
 * @param {Tree} tree
 * @param {Dictionary} model_template
 * @param {Dictionary} initial_values
 * @param {Dictionary} run_options - additional option to fix the frequency to the one given in the initial values
 */
lfunction custom_functions.CreateLFObject_FixedBLs (context, data_filter, tree, model_template, initial_values, run_options) 
{
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
	if (run_options["set_freq"] == TRUE) {
		intial_frequencies = (initial_values[utility.getGlobalValue("terms.efv_estimate")])[utility.getGlobalValue("terms.model")];
		states_num = Columns((^(user_model_id))[utility.getGlobalValue("terms.alphabet")]);
		for (i=0; i<states_num; i+=1) {
			((^(user_model_id))[utility.getGlobalValue("terms.efv_estimate")])[i][0] = intial_frequencies[i][0];
		}
	}
	
	for (i = 0; i < components; i += 1) {
        lf_components[2 * i + 1] = "`context`.tree_" + i;
        custom_functions.ApplyModelToTree_FixedBLs(lf_components[2 * i + 1], tree[i], {
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
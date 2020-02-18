/* SLKP */	
RequireVersion ("2.220141023");

/* imports */
LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("TreeTools");
// namespace 'estimators' for various estimator related functions
LoadFunctionLibrary("BranchSiteTemplate");

/* hard-coded parameters */
REPLACE_TREE_STRUCTURE = 1; // assures initialization of the tree every time it is re-loaded
AUTOMATICALLY_CONVERT_BRANCH_LENGTHS = 1;   
ACCEPT_ROOTED_TREES   = TRUE; 
scriptsDir = "/groups/itay_mayrose/halabikeren/HyPhy/";
chars = {{"0","1"}};
stateCharCount = 2;
numOfHistories = 1;
maxSimulationsNum = 1000;
OPTIMIZATION_METHOD     = 0;
USE_LAST_RESULTS	= 1;
VERBOSITY_LEVEL         = 1;



/* ====================================================================================	*/
/*		MAIN LOOP																		*/
/* ==================================================================================== */

/* process input from command line */

// ask for a tree file
/* fprintf(MESSAGE_LOG, "Please provide a path to tree file in newick format: ");
fscanf(stdin, "String", treePath);

// ask for characterData (i.e, traitData) file 
fprintf(MESSAGE_LOG, "Please provide a path to trait data file in fasta format: ");
fscanf(stdin, "String", traitDataPath);

// ask for sequenceData file
fprintf(MESSAGE_LOG, "Please provide a path to sequence data file in nexus format: ");
fscanf(stdin, "String", sequenceDatapath);
fprintf(MESSAGE_LOG, "\n"); */

// for debugging
treePath = "./data/tree.nwk";
traitDataPath = "./data/charData.fas";
sequenceDatapath = "./data/seqData.nex";

// process the character (trait) data
DataSet traitData  = ReadDataFile(traitDataPath);                 
DataSetFilter leavesCharData = CreateFilter(traitData,1);

// process the sequence data
ExecuteAFile (scriptsDir + "chooseGeneticCode.def");              // NEED TO REDUCE
DataSet codon_data = ReadDataFile (sequenceDatapath);
codon_filter = processCodonData (codon_data);                    // generate codon data and filter instances and assert their content

// set the character model
HarvestFrequencies(CharsFreq,leavesCharData,1,1,1); 
global mu = 1.0;                          // set the initial value of mu
ACCEPT_ROOTED_TREES   = TRUE; 
RateMatrix = {{*,mu*t}                                  
			  {mu*t,*}};
Model charModel = (RateMatrix, CharsFreq);
fscanf(treePath, "Tree", inTree);                               // extract the tree from the tree file. This must happen after defining the model 
LikelihoodFunction charModel_LF = (leavesCharData, inTree);     // set the likelihood function	
	
// set the sequence model - seqModel, which would apply to all the histories
HarvestFrequencies(nucFreqByPos, codon_filter, 3, 1, 1);                              
codon_frequencies = getCodonFrequencies(nucFreqByPos);          // this field  sets the codon frequencies according to the local CF3X4 model																						
                                                                // optimization is done here, in function CF3x4() (line 871)

// NEED TO COMPLETE - BE RELAX!

// initialize the input for the joint likelihood function
// jointLF = "" 

// run the traitRELAX pipeline on numOfHistories character (trait) histories                                              
for (this_rep = 0; this_rep < numOfHistories; this_rep = this_rep+1) {    
	
	// generate a stochastic mapping and put it in the path labelledTreePath
	labelledTreeStr = generateHistory(charModel_LF, inTree);
	fprintf(stdout, "labelledTreeStr: ", labelledTreeStr, "\n"); // debug
}
	
	// load the tree string into a string object
	/*ExecuteCommands ("Tree T_" + this_rep + " = labelledTreeStr;\n"); // I am using the wrapper ExecuteCommands in order to allow the trees to have names 
	                                                                  // corresponding to the history's number
	
	// add the tree loaded with the model to the input string of the joint likelihood function 
	if (this_rep) {
		jointLF = jointLF + ","
	}
	jointLF = jointLF + "codon_filter,T_" + this_rep;
	fprintf(MESSAGE_LOG, "jointLF: ", jointLF, "\n"); // debug!
}

// optimize the joint likelihood function - the character model expression is still missing!
ExecuteCommands ("LikelihoodFunction lf = (" + jointLF + ");\n"); // jointLF is a string that holds 
Optimize (res,lf); // optimizes the joint LF!!
USE_LAST_RESULTS	= 0;
OPTIMIZATION_METHOD     = 4; */

 
// Auxiliary functions
/*______________General Auxiliary Functions_____________________________*/
/* ____________________________________________________________________	*/

/* TRY RE-IMPLEMENTING TREE GENERATION AND ONLY IF FAILED ADD THIS AUXILIARY FUNCTION: */
/* auxiliary function to convert nodes names in the data */

/* auxiliary function to replace src sub-string with dest sub-string    */
/* in string                                                            */
function strReplace(string, src, dest)
{
	replacor = {2,1};
	replacor[0] = src;
	replacor[1] = dest;
	string = string^replacor;
	return string;

}

/* process the codon sequence data                                     */
function processCodonData (codon_data) 
{ 
    DataSetFilter codon_filter = CreateFilter (codon_data, 3, , , GeneticCodeExclusions);
	assert (codon_filter.sites*3==codon_data.sites, "The input alignment must not contain stop codons")
	// to do: assert no stop codon in the middle
    return codon_filter;
}

/* get the codon frequencies by CF3X4 model                            */
function getCodonFrequencies (nucFreq) 
{
    nucCF = CF3x4	(nucFreq, GeneticCodeExclusions); // CHECK: optimization is done here! what kind of optimization is it?
    codon3x4 = BuildCodonFrequencies (nucCF);
    return {"raw": nucFreq, 
            "nucleotide": nucCF,
            "codon": codon3x4};
}


/*______________SM Auxiliary Functions__________________________________*/
/* ____________________________________________________________________	*/

/* ____________________________________________________________________	*/
/*  map node indices to their father node indices                       */
function mapNodeIDtoFatherID(nodesNum, nodeIdToBranchIdMap, nodePOTIdToFatherPOTId)
{
	nodeIDToFatherID = {nodesNum, 1};
	for (sonID=0; sonID<nodesNum; sonID=sonID+1) {
		sonBranchID = nodeIdToBranchIdMap[sonID][0];
		sonPOTID = sonBranchID;
		fatherPOTID = nodePOTIdToFatherPOTId[sonPOTID];
		fatherID = -1;
		// if fatherPOTID==1 it means that the son is the root. We are going to keep -1 as the father index
		// otherwise we update:
		if (fatherPOTID != -1)
		{
			fatherBranchID = fatherPOTID;
			fatherID = nodeIdToBranchIdMap[fatherBranchID][1];
		}
		nodeIDToFatherID[sonID] = fatherID;
		
	}
	return nodeIDToFatherID;
}


/* ____________________________________________________________________	*/
/*  map node indices to the indices of their lest most leaf indices     */
function mapNodeIDtoLeftMostLeafID(nodesNum, nodeIDToFatherID, nodeIdToBranchId)
{
	// initialize the map
	NodeIDtoLeftMostLeafID = {nodesNum, 1};
	for (nodeID=0; nodeID<nodesNum; nodeID=nodeID+1) {
		NodeIDtoLeftMostLeafID[nodeID] = -1;
	}
	
	// fill the map by post order traversal update
	for (branchID=0; branchID<nodesNum; branchID=branchID+1) {
		nodeID = nodeIdToBranchId[branchID][1]; 
		if (nodeID<leavesNum) {                             // if the node is a leaf -> put itself in the map
			NodeIDtoLeftMostLeafID[nodeID] = nodeID;
		}
		fatherID = nodeIDToFatherID[nodeID];               // go to its father and update if required
		if (fatherID != -1) {
			if (NodeIDtoLeftMostLeafID[fatherID] == -1) { 
				NodeIDtoLeftMostLeafID[fatherID] = NodeIDtoLeftMostLeafID[nodeID];
			}
		}
	}
	
	return NodeIDtoLeftMostLeafID;
}     


/* ____________________________________________________________________	*/
/*  map node indices to their corresponding branch indices and vice  	*/
/*	versa.   															*/
function setupMapToTree (branchNames, leavesCharData, internalsCharData)
{	    
	nodesNum = Columns (branchNames);             // get total number of branches in tree 
	nodeIdToBranchIdMap = {nodesNum, 2};	      // initialize a 2xh matrix: the first column maps node indices to branch indices
	                                              //                          the second column maps branch indices to node indices
	
	// start with the mapping the leaves
	leavesNum = leavesCharData.species;
	for (k=0; k<leavesNum; k=k+1)	
	{
		GetString (seqName, leavesCharData, k);
		for (v=0; v<nodesNum; v=v+1)
		{
			if (branchNames[v] % seqName) {	 // find a branch with a matching string label to the one of the leaf node
				nodeIdToBranchIdMap[k][0] = v;	 // map leaf index to branch index
				nodeIdToBranchIdMap[v][1] = k;	 // map branch index to leaf index
				break;
			}
		}
	}

	// map the species number to the branches number and vice versa
	// map the root node (with index = leavesNum) to the last (demme) branch (with index = h-1)
	nodeIdToBranchIdMap[leavesNum][0] = nodesNum-1; 
	nodeIdToBranchIdMap[nodesNum-1][1] = leavesNum;
	
	// continue with internal nodes
	internalsNum = internalsCharData.species;
	for (k=1; k<internalsNum; k=k+1) {                // start k from 1 because 0+leavesNum is the index of the root
		GetString (seqName, internalsCharData, k); 
		nodeIdToBranchIdMap[leavesNum+k][0] = -1;         // initialize the mapped branch index
		for (v=0; v<nodesNum; v=v+1) {
			if (branchNames[v] % seqName) {           // find branch with a matching string label to the one of the node
				nodeIdToBranchIdMap[k+leavesNum][0] = v;  // map node index (given in increasing order from the last leaf index) to branch index 
				nodeIdToBranchIdMap[v][1] = k+leavesNum;  // map branch index to node index
				break;
			}
		}
	}

	GetDataInfo    (leavesDupInfo, leavesCharData);             // spill the leaves character data into a matrix - leavesDupInfo. 
	GetDataInfo	   (internalsDupInfo, internalsCharData);       // spill the internal nodes character data into a matrix - internalsDupInfo. 
	                                                                // the function GetDataInfo is given 3 arguments. the third argument has a default value - Blank
																	// the Blank parameter dictates that the resulting matrix will contain a site duplicate map
																	// meaning, a set of unique patterns and a map that matches each site to its pattern.
	traitIndices  = {1,stateCharCount};                             // initialize a vectors with a cell for each unique character index, which will map characters to indices

	for (h=Columns(traitIndices)-1; h>=0; h=h-1) {
		traitIndices  [h] = h;                                      // traitIndices holds a unique index for each unique character
	}
 
	charInfoLeafs  = {leavesNum, 1};	   // initialize a matrix with a row per leaf and a column per unique character pattern of the leaves
                                           // character data along the tree
	charInfoInternals = {internalsNum, 1}; // initialize a matrix with a row per internal node and a column per unique character pattern of 
	                                       // the internal nodes character data along the tree

	// translate the character data of the leaves into indices
	// for each combination of species, charInfoLeafs will hold the index corresponding to the character assigned to the given species
	for (h=0; h<leavesNum;h=h+1) {
		GetDataInfo (traitInfo, leavesCharData, h, 0); // this command sets a vector traitInfo of length corresponding to the number of unique characters, 
														   // matching the character data of the species in index h, in position v.
														   // this vector holds the value 1 in the index that matches to the character which is assigned to the 
														   // species in index h in position v, and 0 in all the other entries.
		trait = traitIndices * traitInfo;                  // this returns the index of the character that is assigned to species h in position v by a simple 
														   // vectors multiplication.
		charInfoLeafs[h][0] = trait[0];
	}

	// translate the character data of the internal nodes into indices
	// for each combination of node, charInfoLeafs will hold the index corresponding to the character assigned to the given node
	for (h=0; h<internalsNum;h=h+1) {
		GetDataInfo (traitInfo, internalsCharData, h, 0);
		traitInfo = traitIndices * traitInfo;
		charInfoInternals[h][v] = traitInfo[0];
	}
													  // of nodes. Each node has its traversal index (root is last, leftmost leaf is first. 0-based).
													  // The print is such that the traversal index of the father of the node is printed. For the root it
													  // is -1 (no father)
	GetInformation (seqStrings, leavesCharData);  // seqStrings is a vector of strings matching to the character data of the leaves
	                                                  // the order of the strings matches the order of the leaves in the dataSet
	return nodeIdToBranchIdMap;
}


/* ____________________________________________________________________	*/
/* generate a history of a given branch and the characters at its edges */
function sampleMutationsGivenAncestralsPerBranch(branchesHistory, branchesTransitionsCounter, sonName, fatherState, sonState, branch_length)
{
	// reset a vector of transitions for the branch
	transitionsRecorder = {};
	tansitionsCounter = 0;
	
	for (i=0; i<maxSimulationsNum;i=i+1) {                               // do not exceed the acceptable number of transitions along a branch
		disFromNode = 0;
		curState = fatherState;
		lambdaExpParam = -1.0 * RateMatrix[curState][curState]; 		 // the time to change is exponentially distributed with parameter lambda = sum of rates out 
																		 // of current state. The mean of this exponential distribution is 1/lambda.
																		 // RateMatrix[curState][curState] is negative, therefore, we multiply by -1.
		uniform_helper = Random(0,1);             				  		 // sample timeTillChange from a exp(lambdaExpParam) distribution
		timeTillChange = - (1.0 / lambdaExpParam) * Log(uniform_helper); // with a Uniform helper. See here: https://en.wikipedia.org/wiki/Inverse_transform_sampling

		// STOPPING POINT 17.7.17
		
		// if the states of the father and son are different, there has so be at least one transition
		// get the time in which it occurs, timeTillChange, according to Nielsen 2001, equations A1,A2
		if (fatherState != sonState) { 
			u = Random(0,1);
			tmp = u * (1 - Exp(RateMatrix[curState][curState] * branch_length)); // no need to multiply RateMatrix[curState][curState] by -1 as its is already negative
			timeTillChange =  Log(1 - tmp) / RateMatrix[curState][curState];     // unlike the demonstration in the article
		}
		
		// STOPPING POINT 7.8.17
        
		// as long as the last jump didn't exceed the branch length -> add the current state and time to branch history and draw next state
		while (disFromNode + timeTillChange < branch_length) {
		
			// NOTE TO ELI: I THIBK THAT SINCE THE SWITCH HERE IS BINARY, WE DON'T REALLY NEED THE JUMP PROBABILITIES (WE ALWAYS JUMP TO THE COMPLEMENT OF CurState)
			
			JumpProbabilities = {stateCharCount, 2};       // generate a vector of transition probabilities
													       // for each destination state i, the jump probability is rate(curState->i)/rate(leaving curState)
			for (r=0; r<stateCharCount; r=r+1) {
				JumpProbabilities[r][0] = r;
				JumpProbabilities[r][1] = 0;
				if (r != curState) {
					JumpProbabilities[r][1] = RateMatrix[curState][r] / (-1 * RateMatrix[curState][curState]);
				}
			}
			sample = Random(JumpProbabilities, {"PDF":"Multinomial","ARG0":1}); // sample a single destination state (ARG0=1) from a multinomial distribution with the jump probabilities
			for (m=0; m<stateCharCount; m=m+1) {
				if (sample[m][1] == 1) {
					former = curState;
					curState = sample[m][0];
				}
			}
			data = {2,1};                                                       // add the transitions to the recorder
			data[0] = timeTillChange;
			data[1] = curState;
			transitionsRecorder[tansitionsCounter] = data;
			tansitionsCounter = tansitionsCounter + 1;
			disFromNode = disFromNode + timeTillChange;
			lambdaExpParam = -1.0 * RateMatrix[curState][curState];         
			uniform_helper = Random(0,1);             				  		   // sample timeTillChange from a exp(lambdaExpParam) distribution
			timeTillChange = - (1.0 / lambdaExpParam) * Log(uniform_helper);   // with a Uniform helper. See here: https://en.wikipedia.org/wiki/Inverse_transform_sampling
       	}
		// the last jump exceeded the branch length
		if (curState == sonState) {                                           // if the current state is the son's state -> the simulation succeeded -> record the last jump
			data = {2,1};
			data[0] = branch_length - disFromNode;
			data[1] = curState;
			transitionsRecorder[tansitionsCounter] = data;
			tansitionsCounter = tansitionsCounter + 1;
			branchesTransitionsCounter[sonName] = tansitionsCounter;
			branchesHistory[sonName] = transitionsRecorder;
			return 0;
		}
	}
	// if all simulations failed -> exit
	fprintf(MESSAGE_LOG, "could not produce simulations with father = ", fatherState, " son = ", sonState, " branch length = ", bl + "\n\n");
	return 1;

}


/* ____________________________________________________________________	*/
/*	generate simulated history in internal nodes and long branches   	*/
function generateHistory (charModel_LF, origTree)
{
	// get the indices of the tree branches from the tree instance
	branchNames = BranchName (origTree, -1); // -1 means "retrieve for all branches in tree". the branch names are ordered by the post-order traversal.
	nodesNum = Columns (branchNames);
	
	// generate character data in the internal nodes, according to the likelihood function of the character model
	DataSet ancestralSeqs = SampleAncestors (charModel_LF);                                // reconstruct ancestors using the given likelihood function
	DataSetFilter	internalsCharData  = CreateFilter(ancestralSeqs, 1);                   // generate data set filter of character data from reconstructed ancestors 
	
	// get auxiliary maps for the history generation process
	nodeIdToBranchIdMap = setupMapToTree (branchNames, leavesCharData, internalsCharData); // generate a map that fits node indices to branch indices and vice versa
	                                                                                       // this function also grants internal nodes indices
	nodeIdToFatherTraversalId = Abs(origTree);                                             // Abs returns a post-order traversal of the tree. number of entries is the same as number 
	nodeIDToFatherID = mapNodeIDtoFatherID(nodesNum, nodeIdToBranchIdMap, nodeIdToFatherTraversalId);
	NodeIDtoLeftMostLeafID = mapNodeIDtoLeftMostLeafID(nodesNum, nodeIDToFatherID, nodeIdToBranchIdMap);
	
	// initialize the arrays that will hold the mapping history
	branchesHistory = {};                                       // histories along branches
	branchesTransitionsCounter = {};                            // counter of transitions per branch 
	
	// generate mapping history for the internal branches
	sonId = leavesNum + 1;		                              // the index of the first internal node is leavesNum + 1 (the right-most child of the root, whose index is leavesNum)
	for (h = 1; h < internalsNum; h = h + 1) {                // iterate over internal nodes - up to bottom
		
		// get the Id and the name of the node and its father
		branchId = nodeIdToBranchIdMap[sonId][0];	            // get the index of the branch corresponding to the node
		fatherId = nodeIDToFatherID[sonId];                     // get the index of the of the father of the node
		sonName = branchNames [branchId];
		fatherName = branchNames [fatherPOTId];              
		branchLength = BranchLength (origTree, sonName);          // get the length of the branch connecting the son to its father
		
		// get the character states of the father and the son from the internals data filter - charInfoInternals
		sonByInternalsId = h;                                    // h is the node's index according to the internal nodes data index
		fatherByInternalsId = fatherId - leavesNum;              // get the index of the father node according to charInfoInternals
		sonState = charInfoInternals[sonByInternalsId] [0];	     // get the character index of the state of the node
		fatherState = charInfoInternals[fatherByInternalsId][0]; // get the character index of the state of the father	
		
		// generate the history along the branch connecting the node to its father
		if (branchLength > 0) {
			sampleMutationsGivenAncestralsPerBranch(branchesHistory, branchesTransitionsCounter, sonName, fatherState, sonState, branchLength);
		}
		sonId = sonId + 1;		                                 // move on to the next internal node (up to bottom)
	}
	
	// handle the leaves and their corresponding branches
	for (h = 0; h < leavesNum; h = h + 1) {
		
		// get the Id and the name of the node and its father
		sonId = h;                                                // the leavesIds in the tree match their Ids in the data -> no need for reduction
		branchId = nodeIdToBranchIdMap[sonId][0];                
		fatherId = nodeIDToFatherID[sonId][1];                    // get the index of the node mapped to the branch that starts at the father
		sonName = branchNames [branchId];
		fatherName = branchNames [fatherPOTId];   
		branchLength = BranchLength (origTree, sonName);            // get the length of the branch connecting the son to its father
		
		// get the character states of the father and the son from the internals data filter - charInfoInternals
		fatherByInternalsId = fatherId - leavesNum;
		sonState = charInfoLeafs [sonId] [0];                     // get the character index of the state of the node
		fatherState = charInfoInternals[fatherByInternalsId][0];  // get the character index of the state of the father	
		
		// generate the history along the branch connecting the node to its father
		if (branchLength > 0) {
			sampleMutationsGivenAncestralsPerBranch(branchesHistory, branchesTransitionsCounter, sonName, fatherState, sonState, branchLength);
		}
	}
	
	// convert the SM info to tree with branch classification
	// represent the history as a tree with internal nodes, that have a single child
	// the trait states per branch in the new tree will dictate the label of the branch	
	labelledTreeString = generateTree(inTree, nodesNum, branchesHistory, branchesTransitionsCounter, branchNames, nodeIDToFatherID);
	return labelledTreeString;
}


/* ____________________________________________________________________	*/
/*	write the history as a tree and return a tree object                */
function generateTree(tree, nodesNum, branchToHistory, branchToTransitionsNum, branchToName, nodeToFather)
{
	GetString(TreeStr, inTree, 0);
	UseModel (USE_NO_MODEL);
	Tree labelled_tree = TreeStr;
	fprintf(stdout, "TreeStr: ", TreeStr, "\n"); // debug
	fprintf(stdout, "labelled_tree: ", labelled_tree, "\n"); // debug
	internalsCounter = 0;             // counter of added internal nodes (required for unique naming)
	for (nodeId=0; nodeId<nodesNum-1; nodeId=nodeId+1) {  // replace each node in the tree string by the string representing the path along the branch corresponding to it
		nodeName = branchToName[nodeId];
		branchTransitionsNumber = branchToTransitionsNum[nodeName];
		branchHistory = branchToHistory[nodeName];
		formerNode = branchToName[nodeToFather[nodeId]];
		for (j=0; j<branchTransitionsNumber; j=j+1) {
			
			// extract the transition information
			data = branchHistory[j];
			branchLength = data[0];
			classification = data[1]; // the classification is opposite to the assignment, as it is complement to the state in the previous internal node
			
			// convert the data into a corresponding sub-string
			if (j == 0) {                           // handle the special case in which the node is the first one (at the lower edge of the branch) 
				prefix = nodeName;
			} else {
				prefix = "internal" + internalsCounter;
				internalsCounter = internalsCounter + 1;
			}
			if (classification == 1) {              // if the branch state is 0 -> label as 'Donor' ("R set")
				lastNode = prefix + "_Donor_";
			} else {                                // if the branch state is 1 -> label as 'Recipient' ("T set")
				lastNode = prefix + "_Recipient_";
			}
			// add the new node to the tree
			fprintf(stdout, "\nformerNode: ", formerNode, "\n"); // debug
			fprintf(stdout, "lastNode: ", lastNode, "\n");     // debug
			labelled_tree + {"NAME": lastNode, "WHERE" : formerNode, "LENGTH" : branchLength};
			// T + {"NAME": "internal2", "WHERE" : "internal1", "LENGTH" : 0.472138904374846};
			fprintf(stdout, "labelled_tree: ", labelled_tree, "\n\n"); // debug
			formerNode = lastNode;
		}
	}
	GetString(labelledTreeStr, labelled_tree, 0);
	labelledTreeStr = strReplace(labelledTreeStr, "_:", "}:");
	labelledTreeStr = strReplace(labelledTreeStr, "_", "{");
	return labelledTreeStr;
}


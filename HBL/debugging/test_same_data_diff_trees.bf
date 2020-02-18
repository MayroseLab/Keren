dataPath = "./binData.2";  // 1 substitution
// REPLACE_TREE_STRUCTURE = 1; // assures initialization of the tree every time it is re-loaded

TreeStr1 = "((S1:0.1,S2:0.1),(S3:0.1,S4:0.1));";         // T1: all species are in the same distance from the others
TreeStr2 = "((S1:0.1,S2:0.1),(S3:0.1,S4:0.001));";       // T2: the different species (S4) is closer to all
TreeStr3 = "((S1:0.1,S2:0.1),(S3:0.1,S4:1000));";        // T3: the different species (S4) is further from all


// load the data
DataSet data = ReadDataFile (dataPath);

AUTOMATICALLY_CONVERT_BRANCH_LENGTHS = 1; 

DataSetFilter dataFilter = CreateFilter (data,1);



fprintf(stdout, "\n#####################################\n");

// set the model with low probability to leave the state you're in
fprintf(stdout, "\n\nattempting model with preserving branch lengths\n");


// // // // mu = 0.001:
fprintf(stdout, "low leaving rate: \n");
global mu = 0.001;
RateMatrix = {{*,mu*t}
			 {mu*t,*}};			 
Freqs = {{0.5},{0.5}};
Model M = (RateMatrix, Freqs);	

// load the model onto tree
Tree T1 = TreeStr1;
Tree T2 = TreeStr2;
Tree T3 = TreeStr3;

LikelihoodFunction L1 = (dataFilter, T1);
LikelihoodFunction L2 = (dataFilter, T2);
LikelihoodFunction L3 = (dataFilter, T3);
fprintf(stdout, "same dist: ", L1, "\n");
fprintf(stdout, "lower dist: ", L2, "\n");
fprintf(stdout, "higher dist: ", L3, "\n");

// // // // mu = 1:
fprintf(stdout, "medium leaving rate: \n");
global mu = 1;
RateMatrix = {{*,mu*t}
			 {mu*t,*}};			 
Freqs = {{0.5},{0.5}};
Model M = (RateMatrix, Freqs);	

// load the model onto tree
Tree T1 = TreeStr1;
Tree T2 = TreeStr2;
Tree T3 = TreeStr3;

LikelihoodFunction L1 = (dataFilter, T1);
LikelihoodFunction L2 = (dataFilter, T2);
LikelihoodFunction L3 = (dataFilter, T3);
fprintf(stdout, "same dist: ", L1, "\n");
fprintf(stdout, "lower dist: ", L2, "\n");
fprintf(stdout, "higher dist: ", L3, "\n");


// // // // mu = 1000:
fprintf(stdout, "high leaving rate: \n");
global mu = 1000;
RateMatrix = {{*,mu*t}
			 {mu*t,*}};			 
Freqs = {{0.5},{0.5}};
Model M = (RateMatrix, Freqs);	

// load the model onto tree
Tree T1 = TreeStr1;
Tree T2 = TreeStr2;
Tree T3 = TreeStr3;

LikelihoodFunction L1 = (dataFilter, T1);
LikelihoodFunction L2 = (dataFilter, T2);
LikelihoodFunction L3 = (dataFilter, T3);
fprintf(stdout, "same dist: ", L1, "\n");
fprintf(stdout, "lower dist: ", L2, "\n");
fprintf(stdout, "higher dist: ", L3, "\n");

dataPath1 = "./binData.1";  // no substitutions
dataPath2 = "./binData.2";  // 1 substitution
dataPath3 = "./binData.3";  // 2 substitutions
REPLACE_TREE_STRUCTURE = 1; // assures initialization of the tree every time it is re-loaded
AUTOMATICALLY_CONVERT_BRANCH_LENGTHS = 1; 
LIKELIHOOD_FUNCTION_OUTPUT = 0;

TreeStr = "((S1:0.1,S2:0.1),(S3:0.1,S4:0.1));";         // T1: all species are in the same distance from the others

// load the data
DataSet data1 = ReadDataFile (dataPath1);
DataSetFilter dataFilter1 = CreateFilter (data1,1);
DataSet data2 = ReadDataFile (dataPath2);
DataSetFilter dataFilter2 = CreateFilter (data2,1);
DataSet data3 = ReadDataFile (dataPath3);
DataSetFilter dataFilter3 = CreateFilter (data3,1);

// set the model with low probability to leave the state you're in
fprintf(stdout, "\n\nattempting model with no preserving branch lengths\n");

fprintf(stdout, "low leaving rate: \n");
global mu = 0.001;
RateMatrix = {{*,mu}
			 {mu,*}};			 
Freqs = {{0.5},{0.5}};
Model M = (RateMatrix, Freqs);	

// load the model onto tree
Tree T = TreeStr;

LikelihoodFunction L1 = (dataFilter1, T);
LikelihoodFunction L2 = (dataFilter2, T);
LikelihoodFunction L3 = (dataFilter3, T);
fprintf(stdout, "no substitutions: ", L1, "\n");
fprintf(stdout, "1 substitution: ", L2, "\n");
fprintf(stdout, "2 substitutions: ", L3, "\n");

fprintf(stdout, "medium leaving rate: \n");
global mu = 1;
RateMatrix = {{*,mu}
			 {mu,*}};			 
Freqs = {{0.5},{0.5}};
Model M = (RateMatrix, Freqs);	

// load the model onto tree
Tree T = TreeStr;

LikelihoodFunction L1 = (dataFilter1, T);
LikelihoodFunction L2 = (dataFilter2, T);
LikelihoodFunction L3 = (dataFilter3, T);
fprintf(stdout, "no substitutions: ", L1, "\n");
fprintf(stdout, "1 substitution: ", L2, "\n");
fprintf(stdout, "2 substitutions: ", L3, "\n");

fprintf(stdout, "high leaving rate: \n");
global mu = 1000;
RateMatrix = {{*,mu}
			 {mu,*}};			 
Freqs = {{0.5},{0.5}};
Model M = (RateMatrix, Freqs);	

// load the model onto tree
Tree T = TreeStr;

LikelihoodFunction L1 = (dataFilter1, T);
LikelihoodFunction L2 = (dataFilter2, T);
LikelihoodFunction L3 = (dataFilter3, T);
fprintf(stdout, "no substitutions: ", L1, "\n");
fprintf(stdout, "1 substitution: ", L2, "\n");
fprintf(stdout, "2 substitutions: ", L3, "\n");

fprintf(stdout, "\n#####################################\n");

// set the model with low probability to leave the state you're in
fprintf(stdout, "\n\nattempting model with preserving branch lengths\n");

fprintf(stdout, "low leaving rate: \n");
global mu = 0.001;
RateMatrix = {{*,mu*resFactor}
			 {mu*resFactor,*}};			 
Freqs = {{0.5},{0.5}};
Model M = (RateMatrix, Freqs);	

// load the model onto tree
Tree T = TreeStr;

LikelihoodFunction L1 = (dataFilter1, T);
LikelihoodFunction L2 = (dataFilter2, T);
LikelihoodFunction L3 = (dataFilter3, T);
fprintf(stdout, "no substitutions: ", L1, "\n");
fprintf(stdout, "1 substitution: ", L2, "\n");
fprintf(stdout, "2 substitutions: ", L3, "\n");

fprintf(stdout, "medium leaving rate: \n");
global mu = 1;
RateMatrix = {{*,mu*resFactor}
			 {mu*resFactor,*}};			 
Freqs = {{0.5},{0.5}};
Model M = (RateMatrix, Freqs);	

// load the model onto tree
Tree T = TreeStr;

LikelihoodFunction L1 = (dataFilter1, T);
LikelihoodFunction L2 = (dataFilter2, T);
LikelihoodFunction L3 = (dataFilter3, T);
fprintf(stdout, "no substitutions: ", L1, "\n");
fprintf(stdout, "1 substitution: ", L2, "\n");
fprintf(stdout, "2 substitutions: ", L3, "\n");

fprintf(stdout, "high leaving rate: \n");
global mu = 1000;
RateMatrix = {{*,mu*resFactor}
			 {mu*resFactor,*}};			 
Freqs = {{0.5},{0.5}};
Model M = (RateMatrix, Freqs);	

// load the model onto tree
Tree T = TreeStr;

LikelihoodFunction L1 = (dataFilter1, T);
LikelihoodFunction L2 = (dataFilter2, T);
LikelihoodFunction L3 = (dataFilter3, T);
fprintf(stdout, "no substitutions: ", L1, "\n");
fprintf(stdout, "1 substitution: ", L2, "\n");
fprintf(stdout, "2 substitutions: ", L3, "\n");

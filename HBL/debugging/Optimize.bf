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

// conserve branches lengths in thee the cases
AUTOMATICALLY_CONVERT_BRANCH_LENGTHS = 1;

// load the tree and the sequence data
simpleTree = "(A:1,B:1,(C:1,D:1):2);";
fprintf(stdout, "Original tree: ", simpleTree, "\n"); // debug
DataSet data = ReadFromString (
"
>A
A
>B
C
>C
A
>D
A
"
);
DataSetFilter dataFilter = CreateFilter (data,1);


// define JC69 model
fprintf(stdout, "\nJC69:\n");
Freqs = {{0.25},{0.25},{0.25},{0.25}};
global mu = 100;
RateMatrix = {{*,mu*t,mu*t,mu*t}
			  {mu*t,*,mu*t,mu*t}
			  {mu*t,mu*t,*,mu*t}
			  {mu*t,mu*t,mu*t,*}};	  
// load the model on the tree and 
Model M = (RateMatrix, Freqs);
GetString (branch_length, M, -1);
fprintf (stdout, "\nBranch length expression = ", branch_length, "\n");
Tree T = simpleTree;
LikelihoodFunction fit = (dataFilter, T);
fprintf(stdout, "tree after LFs definition:\n");
LIKELIHOOD_FUNCTION_OUTPUT = 2;
fprintf(stdout, "fit: ", fit, "\n");
fprintf(stdout, "Parameters after LFs setting:\n");
LIKELIHOOD_FUNCTION_OUTPUT = 4;
fprintf(stdout, "fit and tree properties: ", fit, "\n");
// constrain preservation of branches lengths
branchNames = BranchName(T, -1);
branchesNum = Columns(branchNames);
for (i=0; i<branchesNum-1; i=i+1) { // exclude the root
	BL = BranchLength(T, branchNames[i]);
	BL_expression_leftover = strReplace(branch_length, "\*t", "");
	constraint = "T." + branchNames[i] + ".t := " + BL + "/(" + BL_expression_leftover + ");";
	ExecuteCommands(constraint);
}
// optimize
Optimize(MLEs, fit); 
LIKELIHOOD_FUNCTION_OUTPUT = 2;
fprintf(stdout, "\nParameters after LFs optimization: ", fit, "\n");
LIKELIHOOD_FUNCTION_OUTPUT = 4;
fprintf(stdout, "fit and tree properties: ", fit, "\n");
fprintf(stdout, "T.A.t: ", T.A.t, "\n"); // debug





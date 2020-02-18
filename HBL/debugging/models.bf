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

//////////////////////////////////////////////////////////////

// JC69
// define model
fprintf(stdout, "\nJC69:\n");
Freqs = {{0.25},{0.25},{0.25},{0.25}};
global mu = 1;
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

//////////////////////////////////////////////////////////////

// K80
// define model
fprintf(stdout, "\nK80:\n");
Freqs = {{0.25},{0.25},{0.25},{0.25}};
global a = 1;
global b = 1;
RateMatrix = {{*,b*t,a*t,b*t}
			  {b*t,*,b*t,a*t}
			  {a*t,b*t,*,b*t}
			  {b*t,a*t,b*t,*}};	  
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

// F81
// define model
fprintf(stdout, "\nF81:\n");
global      eqFreqA = 0.25;
global      eqFreqC = 0.25;
global      eqFreqG = 0.25;
global      eqFreqT := 1.0 - eqFreqA - eqFreqC - eqFreqG;      eqFreqT :> 0;
Freqs = {{eqFreqA,eqFreqC,eqFreqG,eqFreqT}};
global mu = 1;
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

// HKY85
// define model
fprintf(stdout, "\nHKY85:\n");
global      eqFreqA = 0.25;
global      eqFreqC = 0.25;
global      eqFreqG = 0.25;
global      eqFreqT := 1.0 - eqFreqA - eqFreqC - eqFreqG;      eqFreqT :> 0;
Freqs = {{eqFreqA,eqFreqC,eqFreqG,eqFreqT}};
global a = 1;
global b = 1;
RateMatrix = {{*,b*t,a*t,b*t}
			  {b*t,*,b*t,a*t}
			  {a*t,b*t,*,b*t}
			  {b*t,a*t,b*t,*}};	    
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

// GTR
// define model
fprintf(stdout, "\nGTR:\n");
global      eqFreqA = 0.25;
global      eqFreqC = 0.25;
global      eqFreqG = 0.25;
global      eqFreqT := 1.0 - eqFreqA - eqFreqC - eqFreqG;      eqFreqT :> 0;
Freqs = {{eqFreqA,eqFreqC,eqFreqG,eqFreqT}};
global a = 1;
global b = 1;
global c = 1;
global d = 1;
global e = 1;
global f = 1;
RateMatrix = {{*,a*t,b*t,c*t}
              {a*t,*,d*t,e*t}
			  {b*t,d*t,*,f*t}
			  {c*t,e*t,f*t,*}};    
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

/* branch length expression in GTR explanation:
GTR:
Branch length expression = 
2*t*eqFreqG*eqFreqT*f+ (G->T AND T->G)
2*t*eqFreqC*eqFreqT*e+ (C->T AND T->C)
2*t*eqFreqC*eqFreqG*d+ (C->G AND G->C)
2*t*eqFreqA*eqFreqT*c+ (A->T AND T->A)
2*t*b*eqFreqA*eqFreqG+ (A->G AND G->A)
2*t*a*eqFreqA*eqFreqC  (A->C AND C->A) */



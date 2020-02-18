LoadFunctionLibrary("libv3/models/model_functions.bf");
LoadFunctionLibrary("libv3/models/DNA/GTR.bf");
LoadFunctionLibrary("libv3/convenience/regexp.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

fprintf(stdout, "checking initial memory status\n");

num = io.PromptUser("\n>Type a number to continue", 0, 1, 10000, FALSE);

/* tree_str = "((a1:1,a2:1):1,a3:1,a4:2);";
alignment = ">a1
A
>a2
C
>a3
A
>a4
A";
DataSet dataset = ReadFromString(alignment);
DataSetFilter df = CreateFilter(dataset, 1); */

// try with bigger data - need to annoy the virtual memory
fscanf("/groups/itay_mayrose/halabikeren/myScripts/HBL/my_data/real_data/A/unlabeled_tree.nwk", "String", tree_str);
DataSet data = ReadDataFile ("/groups/itay_mayrose/halabikeren/myScripts/HBL/my_data/real_data/A/seqData.nex");
DataSetFilter df = CreateFilter (data,1);
	
Freqs = {{0.25},{0.25},{0.25},{0.25}};
RateMatrix = {{*,a*t,b*t,c*t}
			  {a*t,*,d*t,e*t}
			  {b*t,d*t,*,f*t}
			  {c*t,e*t,f*t,*}};
Model model = (RateMatrix, Freqs); 
	
Tree tree = tree_str;

fprintf(stdout, "checking memory status before creation of a likelihood function\n");

num = io.PromptUser("\n>Type a number to continue", 0, 1, 10000, FALSE);

LikelihoodFunction lf1 = (df, tree); 

fprintf(stdout, "checking memory status after creation of the likelihood function\n");

num = io.PromptUser("\n>Type a number to continue", 0, 1, 10000, FALSE);

ComputeLF(lf1);

fprintf(stdout, "checking memory status after computation of the likelihood function\n");

num = io.PromptUser("\n>Type a number to continue", 0, 1, 10000, FALSE);

DeleteObject(lf1);

fprintf(stdout, "checking memory status after deletion of the likelihood function\n");

num = io.PromptUser("\n>Type a number to continue", 0, 1, 10000, FALSE);

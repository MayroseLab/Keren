LoadFunctionLibrary("libv3/models/model_functions.bf");
LoadFunctionLibrary("libv3/models/DNA/GTR.bf");
LoadFunctionLibrary("libv3/convenience/regexp.bf");
LoadFunctionLibrary("libv3/all-terms.bf");

tree_str = "((a1:1,a2:1):1,a3:1,a4:2);";
datafilter_name = "df";
tree_1_name = "T1";
tree_2_name = "T2";
alignment = ">a1
A
>a2
C
>a3
A
>a4
A";
DataSet dataset = ReadFromString(alignment);
DataSetFilter ^ datafilter_name = CreateFilter(dataset, 1);

lf_test(datafilter_name, tree_1_name, tree_2_name, tree_str);

lfunction lf_test(datafilter_name, tree_1_name, tree_2_name, tree_str) {
	
	Freqs = {{0.25},{0.25},{0.25},{0.25}};
	RateMatrix = {{*,a*t,b*t,c*t}
				  {a*t,*,d*t,e*t}
				  {b*t,d*t,*,f*t}
				  {c*t,e*t,f*t,*}};
	Model model = (RateMatrix, Freqs); 
		
	Tree ^ tree_1_name = tree_str;
	Tree ^ tree_2_name = tree_str;

	lf_components = {4,1};
	lf_components[0] = datafilter_name;
	lf_components[1] = tree_1_name;
	lf_components[2] = datafilter_name;
	lf_components[3] = tree_2_name;
	
	
	lf_id = &likelihoodFunction;
    utility.ExecuteInGlobalNamespace ("LikelihoodFunction `lf_id` = (`&lf_components`)"););
	DeleteObject(likelihoodFunction);
}

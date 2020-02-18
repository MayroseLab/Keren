seq_data_path = "/groups/itay_mayrose/halabikeren/myScripts/HBL/my_data/real_data/H/seqData.nex";
tree_path = "/groups/itay_mayrose/halabikeren/myScripts/HBL/my_data/real_data/H/labeled_tree.nwk";
replicates = 30;
num_codons = 200;
_outPrefix = "/groups/itay_mayrose/halabikeren/myScripts/HBL/my_data/simulations/H_aln_based/"; 

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

skipCodeSelectionStep 		= 0;
LoadFunctionLibrary("libv3/models/model_functions.bf");
LoadFunctionLibrary("libv3/convenience/regexp.bf");
LoadFunctionLibrary("libv3/all-terms.bf");
utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
ExecuteCommands ('LoadFunctionLibrary("chooseGeneticCode");',
					 {"0" : "Universal"});	
utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);	
LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");

/*---------------------------------------------------------------------------------------------------------------------------------------------*/

function BuildCodonFrequencies (obsF)
{
	PIStop = 1.0;
	result = {ModelMatrixDimension,1};
	hshift = 0;

	for (h=0; h<64; h=h+1)
	{
		first = h$16;
		second = h%16$4;
		third = h%4;
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			PIStop = PIStop-obsF[first][0]*obsF[second][1]*obsF[third][2];
			continue; 
		}
		result[h-hshift][0]=obsF[first][0]*obsF[second][1]*obsF[third][2];
	}
	return result*(1.0/PIStop);
}

/*---------------------------------------------------------------------------------------------------------------------------------------------*/


CODON_NUMBER = 500;

//define nucleotide transition rates and initialise
global AC = 1;
global AT = 1;
global CG = 1;
global CT = 1;
global GT = 1;

MGCustomRateBiasTerms = {{"AC*","","AT*","CG*","CT*","GT*"}};	

_nucBiasTerms = {4,4};
_nucBiasTerms[0][0] = "";

hv = 0;

for (h=0; h<4; h=h+1)
{
	for (v=h+1; v<4; v=v+1)
	{
		_nucBiasTerms[h][v] = MGCustomRateBiasTerms[hv];
		_nucBiasTerms[v][h] = MGCustomRateBiasTerms[hv];
		hv = hv + 1;	
	}
}

h=0;
v=0;

//populate a model matrix
function PopulateModelMatrix (ModelMatrixName&, EFV, synrateP, globalP)
{
	if (!ModelMatrixDimension)
	{
		ModelMatrixDimension = 64;
		for (h = 0 ;h<64; h=h+1)
		{
			if (_Genetic_Code[h]==10)
			{
				ModelMatrixDimension = ModelMatrixDimension-1;
			}
		}
	}
	
	ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension}; 

	hshift = 0;
	
	modelDefString = "";
	modelDefString*16384;
	
	catCounterAL = {};
	
	for (h=0; h<64; h=h+1)
	{
		if (_Genetic_Code[h]==10) 
		{
			hshift = hshift+1;
			continue; 
		}
		vshift = hshift;
		for (v = h+1; v<64; v=v+1)
		{
			diff = v-h;
			if (_Genetic_Code[v]==10) 
			{
				vshift = vshift+1;
				continue; 
			}
			nucPosInCodon = 2;
			if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0))
			{
				if (h$4==v$4)
				{
					transition = v%4;
					transition2= h%4;
				}
				else
				{
					if(diff%16==0)
					{
						transition = v$16;
						transition2= h$16;
						nucPosInCodon = 0;
					}
					else
					{
						transition = v%16$4;
						transition2= h%16$4;
						nucPosInCodon = 1;
					}
				}
				hs  = Format(h-hshift,0,0);
				vs  = Format(v-vshift,0,0);
				ts  = Format(transition,0,0);
				ts2 = Format(transition2,0,0);
				ps  = Format(nucPosInCodon,0,0);
				aa1 = _Genetic_Code[0][h];
				aa2 = _Genetic_Code[0][v];
				if (aa1==aa2) 
				{
					modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+_nucBiasTerms[transition][transition2]+""+synrateP+"*t*EFV__["+ts+"]["+ps+"];\n"+
									"ModelMatrixName["+vs+"]["+hs+"] := "+_nucBiasTerms[transition][transition2]+""+synrateP+"*t*EFV__["+ts2+"]["+ps+"];\n");
				}
				else
				{
					modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := "+_nucBiasTerms[transition][transition2]+""+synrateP+"*"+globalP+"*t*EFV__["+ts+"]["+ps+"];\n"+
									"ModelMatrixName["+vs+"]["+hs+"] := "+_nucBiasTerms[transition][transition2]+""+synrateP+"*"+globalP+"*t*EFV__["+ts2+"]["+ps+"];\n");	
				}
			}
	    }
    }		
	modelDefString*0;
	ExecuteCommands (modelDefString);
	return 0;
}

//prompt for file
DataSet 			ds 				= ReadDataFile(seq_data_path);
DataSetFilter 		dsf 			= CreateFilter(ds,3,"","",GeneticCodeExclusions);

//harvest frequencies from data
HarvestFrequencies	(nuc3, dsf, 3, 1, 1);
//fprintf(stdout, "nuc3:\n", nuc3, "\n"); // debug
//exit();

//define synonymous rate parameter (the "low" is leftovers)
global synlow := 1;

nucCF						= CF3x4	(nuc3, GeneticCodeExclusions);
//fprintf(stdout, "nucCF:\n", nucCF, "\n"); // debug
//exit();

//For getting initial branch lengths
global omegaInit = 0.4;
PopulateModelMatrix			  ("MGMatrixInit",  nucCF, "synlow", "omegaInit");

//constrain synonymous rate to be 1. Current version fixes tree in advance
//for speed, so this need to be unconstrained later
synlow := 1;

codon3x4					= BuildCodonFrequencies (nucCF);
//MUST CALL THE INITIALIZATION MODEL "FG" TO ENSURE THAT THE ALL BRANCHES GET A MODEL
Model		FG			= (MGMatrixInit, codon3x4, 0);

//prompt for tree - foreground branches must be labeled "{FG}"
utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
ExecuteCommands ('LoadFunctionLibrary			  ("queryTree");',
					 {"0" : tree_path});	
utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);					 

fprintf 					  (stdout, "[PHASE 0] Fitting the local MG94 (no site-to-site variation) to obtain initial parameter estimates\n");
VERBOSITY_LEVEL				 = 0;

//Optimize initial model with global omega
LikelihoodFunction	base_LF	 = (dsf, givenTree);
Optimize					  (res_base,base_LF);

localLL						 = res_base[1][0];
localParams					 = res_base[1][1] + 9;

fprintf(stdout, "\nLog L = ", localLL, " with ", localParams, " degrees of freedom\n");

//BACKGROUND
//define omega values for background
global omega1 = 0.2;
global omega2 = 1.0;
global omega3 = 2.0;
omega1 :< 1;
omega2 :< 1; //SHOULD THIS BE =1? DUNNO!
omega3 :> 1;

//set up model matrices for all permutations of synonymous rates and omega values
PopulateModelMatrix			  ("MGMatrixBG1",  nucCF, "synlow", "omega1");
PopulateModelMatrix			  ("MGMatrixBG2",  nucCF, "synlow", "omega2");	
PopulateModelMatrix			  ("MGMatrixBG3",  nucCF, "synlow", "omega3");

//set up proportion paramenters for omega values
global p1 = 0.55;
p1	:< 1;
global p2 = 0.89;
p2	:< 1;

//FOREGROUND
//define RELAXATION PARAMETER "relax"
global relax = 1;
//NULL CONSTRAINT

//populate model matrices (syn rate is called "synlow" because there used to be low and high)
PopulateModelMatrix			  ("MGMatrix1",  nucCF, "synlow", "omega1^relax");
PopulateModelMatrix			  ("MGMatrix2",  nucCF, "synlow", "omega2^relax");
PopulateModelMatrix			  ("MGMatrix3",  nucCF, "synlow", "omega3^relax");

//NOTE THE SHARED PROPORTION PARAMETERS
Model 		FG		=		  ("Exp(MGMatrix1)*p1+Exp(MGMatrix2)*(1-p1)*p2+Exp(MGMatrix3)*(1-p1)*(1-p2)",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL);
Model 		BG 		=		  ("Exp(MGMatrixBG1)*p1+Exp(MGMatrixBG2)*(1-p1)*p2+Exp(MGMatrixBG3)*(1-p1)*(1-p2)",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL);
//apply appropriate background model to both foreground models' trees
Tree						   fullTree = treeString;

//constrain all 4 trees to be the same as the input tree
ReplicateConstraint 		  ("this1.?.t:=this2.?.t__",fullTree,givenTree);
//ClearConstraints			  (fullTree);
ClearConstraints			  (synlow);

ASSUME_REVERSIBLE_MODELS	  = 1;

//define the likelihood function
LikelihoodFunction	lf	 = (dsf, fullTree);

//NULL MODEL
fprintf 					  (stdout, "[PHASE 1] Fitting the null\n");

omega_sets = {3,1};
omega_sets[0] = 0.001;
omega_sets[1] = 1.0;
omega_sets[2] = 8.0;


k_values = {11,1};
for (k = 0; k <11; k = k + 1) {
    k_values[k][0] = 1 - (k/10)^1.5;
}

for (k = 0; k <= 10; k = k + 1) {
	
	omega1 = omega_sets[0];
	omega2 = omega_sets[1];
	omega3 = omega_sets[2];
	
	relax = k_values[k][0];

	// Print LF
	LIKELIHOOD_FUNCTION_OUTPUT = 7;
	fprintf(_outPrefix + "." + k + ".null.fit", CLEAR_FILE, lf);
	LIKELIHOOD_FUNCTION_OUTPUT = 2;

	for (_i = 0; _i < replicates; _i += 1) {
		DataSet       sim           = SimulateDataSet (lf);
		DataSetFilter simFilter     = CreateFilter (sim, 1);
		_saveTo = _outPrefix + ".pre_k_" + k + "._replicate_" + _i;
		fprintf (_saveTo, simFilter);
		SetParameter (STATUS_BAR_STATUS_STRING, "Replicate " + (_i+1) + " of " + replicates + " generated", 0);
	}
}



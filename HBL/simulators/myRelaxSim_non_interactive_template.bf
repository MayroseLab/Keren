tree_path = "/scratch300/halabikeren/TraitRELAX/real_data/A/labeled_tree.nwk";
replicates = 10;
num_codons = 1000;
_outPrefix = "/scratch300/halabikeren/TraitRELAX/simulations/omegas_0.5_1_2/50_taxa_1000_codons_A_labels/simulated_dataset";

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


// codon universal properties - taken from codonSimulator.bf
charactersUniversalCode = {{"A","C","G","T"}
			  			   {"3",GeneticCodeExclusions,"",""}};

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


CODON_NUMBER = 300;

//define nucleotide transition rates and initialize - according to the fitting results from the SimulateDataSet() based simulation
global AT=0.5133264174461531;
global AC=0.60084650761655;
global CG=0.6357304411729526;
global GT=0.3454754695092268;
global CT=1.053827894165044;

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



nuc3 = {	 /* 1st codon pos*/ /* 2nd codon pos*/ /* 3rd codon pos*/
	/*A*/		{    0.355885780886,    0.377738927739,     0.47027972028} 
	/*C*/		{    0.188053613054,    0.197494172494,    0.115617715618}
	/*G*/		{    0.292657342657,    0.152097902098,    0.194230769231}
	/*T*/		{    0.163403263403,    0.272668997669,    0.219871794872}
				};

//define synonymous rate parameter (the "low" is leftovers)
global synlow := 1;

nucCF						= CF3x4	(nuc3, GeneticCodeExclusions);

//For getting initial branch lengths
global omegaInit = 0.4;
PopulateModelMatrix			  ("MGMatrixInit",  nucCF, "synlow", "omegaInit");

//constrain synonymous rate to be 1. Current version fixes tree in advance
//for speed, so this need to be unconstrained later
synlow := 1;

//codon3x4					= BuildCodonFrequencies (nucCF);
// based on fit results from the SimulateDataSet() based simulation
codon3x4 = {
    {0.02208405111263963}
    {0.02470694553694215}
    {0.02135305117967318}
    {0.02501285536067612}
    {0.0126488958838864}
    {0.01415118902377324}
    {0.01223020721141143}
    {0.01432640241603234}
    {0.01031757586431883}
    {0.01154298111577433}
    {0.009976055767905154}
    {0.01168590089973594}
    {0.01824926779965479}
    {0.02041670992860017}
    {0.01764520229237115}
    {0.02066950006512866}
    {0.01734654335365297}
    {0.01940677005803593}
    {0.01677235875481991}
    {0.01964705559227256}
    {0.009935433481228353}
    {0.01111545217200323}
    {0.009606562606417995}
    {0.0112530784928925}
    {0.008104232150258685}
    {0.009066761407764379}
    {0.007835975518883974}
    {0.009179021799480545}
    {0.01433440420167989}
    {0.0160368830025357}
    {0.01385992384220705}
    {0.01623544417413874}
    {0.02908326520637171}
    {0.03253744731096957}
    {0.02812058563241731}
    {0.03294031073885745}
    {0.01665777676761052}
    {0.01863619954802309}
    {0.01610639089921231}
    {0.01886694423923887}
    {0.01358757927240544}
    {0.01520135863433434}
    {0.01313781942143227}
    {0.01538957473467717}
    {0.02403310391433075}
    {0.0268874848399166}
    {0.0232375887590386}
    {0.02722039307965551}
    {0.01988068996685593}
    {0.02012684335535864}
    {0.01017806013756683}
    {0.01138689528510733}
    {0.009841158124419814}
    {0.01152788248744426}
    {0.009288174797360396}
    {0.008027332699513185}
    {0.009403176625928092}
    {0.01468445521542064}
    {0.01642850912618173}
    {0.01419838788459574}
    {0.01663191922895982}
};
//MUST CALL THE INITIALIZATION MODEL "FG" TO ENSURE THAT THE ALL BRANCHES GET A MODEL
Model		FG			= (MGMatrixInit, codon3x4, 0);

//prompt for tree - foreground branches must be labeled "{FG}"
utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", "/dev/null");
ExecuteCommands ('LoadFunctionLibrary			  ("queryTree");',
					 {"0" : tree_path});	
utility.ToggleEnvVariable ("GLOBAL_FPRINTF_REDIRECT", None);		

fprintf 					  (stdout, "[PHASE 0] Fitting the local MG94 (no site-to-site variation) to obtain initial parameter estimates\n");
VERBOSITY_LEVEL				 = 0;

fprintf (_outPrefix, CLOSE_FILE);

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

//NULL MODEL
fprintf 					  (stdout, "[PHASE 1] Fitting the null\n");
//LFCompute(lf,LF_START_COMPUTE);
//LFCompute(lf,logL);
//LFCompute(lf,LF_DONE_COMPUTE);
//fprintf(stdout, logL);

omega_sets = {3,1};
omega_sets[0] = 0.5;
omega_sets[1] = 1.0;
omega_sets[2] = 2.0;


k_values = {5,1};
k_values[0][0] = 0.2;
k_values[1][0] = 0.5;
k_values[2][0] = 1;
k_values[3][0] = 2;
k_values[4][0] = 5;


for (k = 0; k < 5; k = k + 1) {
	
	omega1 = omega_sets[0];
	omega2 = omega_sets[1];
	omega3 = omega_sets[2];
	
	relax = k_values[k][0];

	// Print LF
	LIKELIHOOD_FUNCTION_OUTPUT = 7;
	//fprintf(_outPrefix + "." + k + ".null.fit", CLEAR_FILE, lf);
	LIKELIHOOD_FUNCTION_OUTPUT = 2;

	for (_i = 0; _i < replicates; _i += 1) {
		DataSet       sim           = Simulate (fullTree,codon3x4,charactersUniversalCode,num_codons,0);
		DataSetFilter simFilter     = CreateFilter (sim, 1);
		_saveTo = _outPrefix + ".k_" + relax + "._replicate_" + _i;
		fprintf (_saveTo, simFilter);
		SetParameter (STATUS_BAR_STATUS_STRING, "Replicate " + (_i+1) + " of " + replicates + " generated", 0);
	}
}


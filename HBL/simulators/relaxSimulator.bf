skipCodeSelectionStep = 1; // skipping the general questions about the input alignment when launching chooseGeneticCode.def - as in this case, there's no such input

LoadFunctionLibrary("chooseGeneticCode");
LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");


/*---------------------------------------------------------------------------------------------------------------------------------------------*/

// simulation parameters


// set omega parameters for the simulations
omega_sets = {4,3};

omega_sets[0][0] = 0.1;
omega_sets[0][1] = 1.0;
omega_sets[0][2] = 2.0;

omega_sets[1][0] = 0.01;
omega_sets[1][1] = 1.0;
omega_sets[1][2] = 4.0;

omega_sets[2][0] = 0.001;
omega_sets[2][1] = 1.0;
omega_sets[2][2] = 8.0;

omega_sets[3][0] = 0.0001;
omega_sets[3][1] = 0.001;
omega_sets[3][2] = 2.0; // set with low proportion here

omega_frequencies = {4,2};

omega_frequencies[0][0] = 0.55;
omega_frequencies[0][1] = 0.89;

omega_frequencies[1][0] = 0.55;
omega_frequencies[1][1] = 0.89;

omega_frequencies[2][0] = 0.55;
omega_frequencies[2][1] = 0.89;

omega_frequencies[3][0] = 0.99;
omega_frequencies[3][1] = 0.01;

// set relaxation parameters for the simulations
k_values = {3,1};
k_values[0][0] = 0.5;
k_values[1][0] = 0.7;
k_values[2][0] = 1;


// path to the input tree to simulate data by
treeString = "(((((((((((((Buchnera_aphidicola_str_Sg{FG}:2.0,Buchnera_aphidicola_str_APS{FG}:1.0){FG}:1.0,(Buchnera_aphidicola_str_Bp{FG}:1.0,Buchnera_aphidicola_str_Cc{FG}:1.0){FG}:1.0){FG}:1.0,Ishikawaella_capsulata{FG}:1.0):1.0,Pantoea_ananatis:1.0):1.0,((((Erwinia_amylovora_CFBP1430:1.0,Erwinia_amylovora_ATCC_49946:1.0):1.0,Erwinia_pyrifoliae:1.0):1.0,Erwinia_tasmaniensis:1.0):1.0,Erwinia_billingiae:1.0):1.0):1.0,((((Klebsiella_variicola:1.0,Klebsiella_pneumoniae:1.0):1.0,(Enterobacter_sp_638:1.0,Enterobacter_cloacae:1.0):1.0):1.0,(((Salmonella_enterica:1.0,Citrobacter_koseri:1.0):1.0,Escherichia_coli:1.0):1.0,Citrobacter_rodentium:1.0):1.0):1.0,(Cronobacter_turicensis:1.0,Cronobacter_sakazakii:1.0):1.0):1.0):1.0,(Edwardsiella_tarda:1.0,Edwardsiella_ictaluri:1.0):1.0):1.0,((((Pectobacterium_wasabiae:1.0,Pectobacterium_atrosepticum:1.0):1.0,Pectobacterium_carotovorum:1.0):1.0,((Dickeya_zeae:1.0,Dickeya_dadantii_Ech586:1.0):1.0,Dickeya_dadantii_Ech703:1.0):1.0):1.0,((((Blochmannia_pennsylvanicus{FG}:1.0,Blochmannia_floridanus{FG}:1.0){FG}:1.0,Wigglesworthia_glossinidia{FG}:1.0){FG}:1.0,Baumannia_cicadellinicola{FG}:1.0){FG}:1.0,Sodalis_glossinidius:1.0):1.0):1.0):1.0,(((Regiella_insecticola:1.0,Hamiltonella_defensa:1.0):1.0,Yersinia_pestis:1.0):1.0,Serratia_proteamaculans:1.0):1.0):1.0,(((Xenorhabdus_nematophila:1.0,Xenorhabdus_bovienii:1.0):1.0,(Photorhabdus_luminescens:1.0,Photorhabdus_asymbiotica:1.0):1.0):1.0,((Riesia_pediculicola:1.0,Arsenophonus_nasoniae:1.0):1.0,Proteus_mirabilis:1.0):1.0):1.0):1.0,(Pasteurella_multocida:1.0,Haemophilus_influenzae:1.0):1.0):1.0,Vibrio_cholerae:1.0):1.0,Pseudomonas_aeruginosa:1.0,Xanthomonas_axonopodis:1.0);";
outputPath = "/groups/itay_mayrose/halabikeren/HyPhy/simulations/simulated_dataset";

// use universal code
geneticCodeID = 0; 
modelType = 0;
ApplyGeneticCodeTable (0); //use universal code

// number of simulated datasets to generate
replicates = 1;
num_codons = 50;

DATA_FILE_PRINT_FORMAT = 4;
/* what format the output should be saved as:
1. DATA_FILE_PRINT_FORMAT = 0:
   -- FASTA sequential;
2. DATA_FILE_PRINT_FORMAT = 1:
   -- FASTA interleaved;
3. DATA_FILE_PRINT_FORMAT = 2: 
   -- PHYLIP sequential;
4. DATA_FILE_PRINT_FORMAT = 3: 
   -- PHYLIP interleaved.
5. DATA_FILE_PRINT_FORMAT = 4:
  -- NEXUS sequential with sequence labels in the matrix;
6. DATA_FILE_PRINT_FORMAT = 5: 
  -- NEXUS interleaved with sequence labels in the matrix;
7. DATA_FILE_PRINT_FORMAT = 6:
  -- NEXUS sequential without sequence labels in the matrix;
8. DATA_FILE_PRINT_FORMAT = 7: 
  -- NEXUS interleaved without sequence labels in the matrix;
9. DATA_FILE_PRINT_FORMAT = 8;
  -- comma separated character data
  */ 

// nucleotides frequencies per position to calculate the frequencies of codons from
observedFreq = {	 /* 1st codon pos*/ /* 2nd codon pos*/ /* 3rd codon pos*/
	/*A*/		{    0.355885780886,    0.377738927739,     0.47027972028} 
	/*C*/		{    0.188053613054,    0.197494172494,    0.115617715618}
	/*G*/		{    0.292657342657,    0.152097902098,    0.194230769231}
	/*T*/		{    0.163403263403,    0.272668997669,    0.219871794872}
				};

// codon universal properties
charactersUniversalCode = {{"A","C","G","T"}
			  			   {"3",GeneticCodeExclusions,"",""}};

// nucleotide initial transition rates
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

//define synonymous rate parameter (the "low" is leftovers)
global synlow := 1;
ModelMatrixDimension = 61;
synlow := 1; //constrain synonymous rate to be 1. Current version fixes tree in advance - for speed, so this need to be unconstrained later

// define omega values and their constraints
global omega1 = 0.2;
global omega2 = 1.0;
global omega3 = 2.0;
omega1 :< 1;
omega2 :< 1; 
omega3 :> 1;

//set up proportion parameters for omega values
global p1 = 0.55;
p1	:< 1;
global p2 = 0.89;
p2	:< 1;
// set the relaxation parameter
global relax = 1; // null constraint: relax=1

ASSUME_REVERSIBLE_MODELS = 1;
		
/*---------------------------------------------------------------------------------------------------------------------------------------------*/

// auxiliary functions

// build from C3X4 nucleotides frequencies per codon position the codon frequencies vector
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

// populate a model matrix according to a given set of parameters
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
			
/*---------------------------------------------------------------------------------------------------------------------------------------------*/

// simulator
			
// extract codon frequencies from nucleotide frequencies per position			
nucCF = CF3x4	(observedFreq, GeneticCodeExclusions);
codon3x4 = BuildCodonFrequencies (nucCF);

//set up model matrices for the background (reference) branches
PopulateModelMatrix("MGMatrixBG1",  nucCF, "synlow", "omega1");
PopulateModelMatrix("MGMatrixBG2",  nucCF, "synlow", "omega2");	
PopulateModelMatrix("MGMatrixBG3",  nucCF, "synlow", "omega3");

//set up model matrices for the foreground (test) branches
PopulateModelMatrix("MGMatrix1",  nucCF, "synlow", "omega1^relax");
PopulateModelMatrix("MGMatrix2",  nucCF, "synlow", "omega2^relax");
PopulateModelMatrix("MGMatrix3",  nucCF, "synlow", "omega3^relax");

// set the FG and BG models accordingly
Model FG = ("Exp(MGMatrix1)*p1+Exp(MGMatrix2)*(1-p1)*p2+Exp(MGMatrix3)*(1-p1)*(1-p2)",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL);       // EXPLICIT_FORM_MATRIX_EXPONENTIAL = 0
Model BG = ("Exp(MGMatrixBG1)*p1+Exp(MGMatrixBG2)*(1-p1)*p2+Exp(MGMatrixBG3)*(1-p1)*(1-p2)",codon3x4,EXPLICIT_FORM_MATRIX_EXPONENTIAL); // the equilibrium frequencies should be integrated into the rate matrix explicitly.

//apply the BG model on the tree (same as null)
Tree fullTree = treeString;

// set constrains on the model before loading the likelihood function
ClearConstraints(synlow);

for (k = 0; k < 3; k = k + 1) {
	for (m = 0; m < 4; m = m + 1) {
		omega1 = omega_sets[m][0];
		omega2 = omega_sets[m][1];
		omega3 = omega_sets[m][2];
		p1 = omega_frequencies[m][0];
		p2 = omega_frequencies[m][1];
		relax = k_values[k][0]; // changing the value of k causes simulation under the alternative model

		for (_i = 0; _i < replicates; _i += 1) {
			DataSet       sim           = Simulate (fullTree,codon3x4,charactersUniversalCode,num_codons,0);
			DataSetFilter simFilter     = CreateFilter (sim, 3, "", "");
			_saveTo = outputPath + ".relax_" + relax + ".omega1_" + omega1 + ".omega2_" + omega2 + ".omega3_" + omega3 + ".p1_" + p1 + ".p2_" + p2 + ".AC_" + AC + ".AT_" + AT + ".CG_" + CG + ".CT_" + CT + ".GT_" + GT + ".replicate_" + _i;
			fprintf (outputPath, "\n", _i, ",", _saveTo);
			
			fprintf (_saveTo, simFilter);
			SetParameter (STATUS_BAR_STATUS_STRING, "Replicate " + (_i+1) + " of " + replicates + " generated", 0);
		}
	}
}

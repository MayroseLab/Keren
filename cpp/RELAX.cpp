// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>

using namespace std;

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/MixedSubstitutionModel.h>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.h>
#include <Bpp/Phyl/Model/FrequenciesSet/FrequenciesSet.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/BppOFrequenciesSetFormat.h>


using namespace bpp;
using namespace std;

#include <map>
#include <vector>
//#include <boost/math/distributions/chi_squared.hpp>

/******************************************************************************/
/**************************** Auxiliary functions *****************************/
/******************************************************************************/

double reportScalingFactor(TreeLikelihood* tl, double origTreeLength)
{
    const Tree& tree = tl->getTree();
    vector <const Node*> nodes = (dynamic_cast<const TreeTemplate<Node>&>(tree)).getNodes();
    double treeSize = 0;
    for (size_t b=0; b<nodes.size(); ++b)
    {
        if (nodes[b]->getId() != tree.getRootId())

        {
            treeSize = treeSize + nodes[b]->getDistanceToFather();
        }
    }
	cout << "\n\norigTreeLength: " << origTreeLength << endl; // debug
	cout << "treeSize: " << treeSize << endl; // debug
    double scalingFactor = treeSize / origTreeLength;
    cout << "The tree has been scaled by a sequence scaling factor of: " << scalingFactor << endl;
	return scalingFactor;
}

RNonHomogeneousMixedTreeLikelihood* getModelParameters(DiscreteRatesAcrossSitesTreeLikelihood* tl, MixedSubstitutionModelSet* modelSet, double origTreeLength)
{
  ParameterList parameters;
  for (size_t m = 0; m < modelSet->getNumberOfModels(); ++m) {
    ApplicationTools::displayMessage("\nmodel " + TextTools::toString(m+1) + "\n");
    TransitionModel* model = modelSet->getModel(m);
    parameters = model->getParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
  }
  parameters = tl->getRateDistributionParameters();
  for (size_t i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
  
  double scalingFactor = reportScalingFactor(tl, origTreeLength);
  return tl;
}

VectorSiteContainer* process_alignment(Alphabet* alphabet, BppApplication bppml)
{
    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppml.getParams()); // here, gaps will be converted to unknown characters
    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppml.getParams(), "", true, false); // convert the alignemnt to condensed format of unique sites
    delete allSites; // delete the non-condenced intance of the sequence data
    SiteContainerTools::changeGapsToUnknownCharacters(*sites); // convert gaps to unknown characters (as done in bppML.cpp)
    return(sites);

}

Tree* process_tree(BppApplication bppml)
{
    Tree* tree = PhylogeneticsApplicationTools::getTree(bppml.getParams());
    TreeTemplate<Node> ttree(*tree);
    vector<Node *> nodes = ttree.getNodes();
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
      if(nodes[i]->isLeaf())
        nodes[i]->setName("Leaf" + TextTools::toString(nodes[i]->getId()));
      else
        nodes[i]->setName("Internal" + TextTools::toString(nodes[i]->getId()));
      ApplicationTools::displayResult(nodes[i]->getName(), TextTools::toString(nodes[i]->getId())); 
    }
    // try to add father to Internal4 (parent of Leaf3) - no such function - need to use detailsSimulation instance

    string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", bppml.getParams(), "Input", "", true, 1); // process tree branches lengths
    string cmdName;
    map<string, string> cmdArgs;
    KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs); // this line process cmdName - which dictates weather the tree should be processed as is, or ultrameterized
    if (cmdName == "Clock") // if needed, turn the tree into an ultrametric tree
    {
      TreeTools::convertToClockTree(*tree, tree->getRootId(), true);
    }
    return(tree);
}

MixedSubstitutionModelSet* setRELAXModel(BppApplication* bppml, const VectorSiteContainer* codon_data, const CodonAlphabet* codonAlphabet, Tree* tree)
{
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml->getParams(), "Standard", "", true, true);
    GeneticCode* gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);
    
    // create the set of models
    MixedSubstitutionModelSet* modelSet = dynamic_cast<MixedSubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(codonAlphabet, gCode, codon_data, bppml->getParams()));
    return modelSet;
}

/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Test non-homogenous models    " << BPP_VERSION << "      *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  try
  { 
  
    // process input from params file
    BppApplication bppml(args, argv, "bppML");
	
	// set random seed
	// process seed from parameter file, if exists
	double seed;
	seed = ApplicationTools::getDoubleParameter("seed", bppml.getParams(), 1);
	if (seed == 1)
	{
		// else, choose a ransom seed
		seed = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
	}
	cout << "seed=" << seed << endl;
	RandomTools::setSeed(static_cast<long>(seed));

    map<string,string> parans = bppml.getParams(); // debug

    /* process alphabet type */
    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppml.getParams(), "", false);
    unique_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml.getParams(), "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    }

    /* process alignment */
    VectorSiteContainer* sites = process_alignment(alphabet, bppml);

    /* process tree */
    Tree* tree = process_tree(bppml);
	double origTreeLength = tree->getTotalLength();
  
    /* set a branch site model */
    MixedSubstitutionModelSet* model = setRELAXModel(&bppml, sites, codonAlphabet, tree);

    /* generate a likelihood function */
    DiscreteDistribution* rDist = new ConstantRateDistribution();
    bppml.startTimer();
    RNonHomogeneousMixedTreeLikelihood* tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true, true);
    tl->initialize();

    /*  save the inital value of k and then set RELAX.k_2 to 1 */
    // extract the user initial value of k for potential later use
    string FGModelInitialValues = ApplicationTools::getStringParameter("model2", bppml.getParams(), "RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1)", "", true, true);
    string modelName = "RELAX";
    map<string, string> cmdArgs2;
    KeyvalTools::parseProcedure(FGModelInitialValues, modelName, cmdArgs2);
    double userInitialValue = TextTools::toDouble(cmdArgs2["k"]);
    tl->setParameterValue("RELAX.k_2", 1);

    /* compute likelihood */
    cout << "\nComputing intial log likelihood" << endl;
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl->getValue(), 15));
    tl = getModelParameters(tl, model, origTreeLength); // debug - print model parameters
    bppml.done();

    /* fit the null */
    cout << "\nFitting the null model" << endl;
    bppml.startTimer();
	bppml.getParam("optimization.ignore_parameters") = "BrLen,RELAX.k_1,RELAX.k_2";
	int scaleTree = ApplicationTools::getIntParameter("optimization.scale.tree", bppml.getParams(), 0);
	double prevLogLikelihood = -tl->getValue();
	double currLogLikelihood = -tl->getValue();
	size_t index = 1;
	do
	{
		cout << "Optimization cycle: " << TextTools::toString(index) << endl;
		index = index + 1;
		if (scaleTree >= 1)
		{
			OptimizationTools::optimizeTreeScale(tl, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
		}
		reportScalingFactor(tl, origTreeLength);
		PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams());
		tl = getModelParameters(tl, model, origTreeLength); // debug - print model parameters
		ApplicationTools::displayResult("Current log likelihood", TextTools::toString(-tl->getValue(), 15));
		prevLogLikelihood = currLogLikelihood;
		currLogLikelihood = -tl->getValue();
		ApplicationTools::displayResult("Current diff", TextTools::toString((currLogLikelihood-prevLogLikelihood), 15));
	} while (currLogLikelihood - prevLogLikelihood > 0.01);
	double nullLogl = currLogLikelihood;
	cout << "iteraive optimzation complete" << endl;
	ApplicationTools::displayResult("Log likelihood", TextTools::toString(nullLogl, 15));
    bppml.done();

    /* fit the alternative */
    // set the initial value of k
    model->getModel(1)->setParameterValue("k", userInitialValue);
    cout << "\nFitting the alternative model" << endl;
	bppml.startTimer();
	bppml.getParam("optimization.ignore_parameters") = "BrLen,RELAX.k_1,RELAX.1_Full.theta_1,RELAX.1_Full.theta1_1,RELAX.1_Full.theta2_1,RELAX.2_Full.theta_1,RELAX.2_Full.theta1_1,RELAX.2_Full.theta2_1,RELAX.3_Full.theta_1,RELAX.3_Full.theta1_1,RELAX.3_Full.theta2_1"; // ignore frequency parameters to reduce optimization duration - results in one unit of ll reduction in optimality and 1 minutre reduction in duration
	prevLogLikelihood = -tl->getValue();
	currLogLikelihood = -tl->getValue();
	index = 1;
	do
	{
		cout << "Optimization cycle: " << TextTools::toString(index) << endl;
		index = index + 1;
		if (scaleTree >= 1)
		{
			OptimizationTools::optimizeTreeScale(tl, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
		}
		reportScalingFactor(tl, origTreeLength);
		PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams());
		tl = getModelParameters(tl, model, origTreeLength); // debug - print model parameters
		ApplicationTools::displayResult("Current log likelihood", TextTools::toString(-tl->getValue(), 15));
		prevLogLikelihood = currLogLikelihood;
		currLogLikelihood = -tl->getValue();
		ApplicationTools::displayResult("Current diff", TextTools::toString((currLogLikelihood-prevLogLikelihood), 15));
	} while (currLogLikelihood - prevLogLikelihood > 0.01);
	cout << "iteraive optimzation complete" << endl;
	double alternativeLogl = currLogLikelihood;
	ApplicationTools::displayResult("Log likelihood", TextTools::toString(alternativeLogl, 15));
    bppml.done();


    /* free resources */
    delete alphabet;
    delete sites;
    delete tree;
    if (model) delete model;
    delete rDist;
    delete tl;
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}
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


void printModelParameters(DiscreteRatesAcrossSitesTreeLikelihood* tl)
{
  ParameterList parameters = tl->getParameters();
  for (size_t p = 0; p < parameters.size(); ++p) 
  {
	ApplicationTools::displayResult(parameters[p].getName(), TextTools::toString(parameters[p].getValue()));
  }
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

MixedSubstitutionModelSet* setBSModel(BppApplication* bppml, const VectorSiteContainer* codon_data, const CodonAlphabet* codonAlphabet, Tree* tree)
{
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml->getParams(), "Standard", "", true, true);
    GeneticCode* gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);
    
    // create the set of models
    MixedSubstitutionModelSet* modelSet = dynamic_cast<MixedSubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(codonAlphabet, gCode, codon_data, bppml->getParams()));
    return modelSet;
}

MixedSubstitutionModel* setSModel(BppApplication* bppml, const VectorSiteContainer* codon_data, const CodonAlphabet* codonAlphabet, Tree* tree)
{
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml->getParams(), "Standard", "", true, true);
    GeneticCode* gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);
    MixedSubstitutionModel* model = dynamic_cast<MixedSubstitutionModel*>(PhylogeneticsApplicationTools::getTransitionModel(codonAlphabet, gCode, codon_data, bppml->getParams()));
	return model;
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

    /* set the model and likelihood function*/
	DiscreteDistribution* rDist = new ConstantRateDistribution();
    MixedSubstitutionModel* model    = 0;
    MixedSubstitutionModelSet* modelSet = 0;
	DiscreteRatesAcrossSitesTreeLikelihood* tl = 0;
	string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppml.getParams(), "no", "", true, 1);
	cout << "nhOpt: " << nhOpt << endl;
	if (nhOpt.compare("no") != 0)
	{
		modelSet = setBSModel(&bppml, sites, codonAlphabet, tree);
		tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, true, true);
	}
	else
	{
		model = setSModel(&bppml, sites, codonAlphabet, tree);
		tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true, true); 
	}
	tl->initialize();

    /* compute likelihood */
    cout << "\nComputing intial log likelihood" << endl;
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl->getValue(), 15));
    printModelParameters(tl); // debug - print model parameters
	
    /* first optimization call */
	OptimizationTools::optimizeTreeScale(tl, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
	tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams()));
	printModelParameters(tl); // debug - print model parameters
	ApplicationTools::displayResult("first log likelihood", TextTools::toString(-tl->getValue(), 15));

    /* second optimization call */
	OptimizationTools::optimizeTreeScale(tl, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
	tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams()));
	printModelParameters(tl); // debug - print model parameters
	ApplicationTools::displayResult("second log likelihood", TextTools::toString(-tl->getValue(), 15));

    /* free resources */
    delete alphabet;
    delete sites;
    delete tree;
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
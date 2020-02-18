// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Io/BppOSequenceReaderFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentReaderFormat.h>
#include <Bpp/Seq/Io/BppOSequenceWriterFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.h>
#include <Bpp/Phyl/Model/FrequenciesSet/FrequenciesSet.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/BppOFrequenciesSetFormat.h>

// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <vector>
//#include <boost/math/distributions/chi_squared.hpp>

using namespace bpp;
using namespace std;

typedef vector<vector<double>> VVDouble;
typedef vector<double> VDouble;
typedef unsigned int uint;

/******************************************************************************/
/**************************** Auxiliary functions *****************************/
/******************************************************************************/

const CodonAlphabet *getCodonAlphabet()
{
  map<string, string> alphabetSettings;
  alphabetSettings["alphabet"] = "Codon(letter=DNA)";
  const Alphabet *alphabet = SequenceApplicationTools::getAlphabet(alphabetSettings, "", false);
  const CodonAlphabet *codonAlphabet = dynamic_cast<const CodonAlphabet *>(alphabet);
  return codonAlphabet;
}

/******************************************************************************/

VectorSiteContainer *processCharacterData(BppApplication *bppml, const BinaryAlphabet *alphabet)
{
  string charDataFilePath = ApplicationTools::getAFilePath("input.character.file", bppml->getParams(), true, true, "", true, "none", 1);
  string sequenceFormat = ApplicationTools::getStringParameter("input.character.format", bppml->getParams(), "Fasta()", "", true, 1);
  BppOAlignmentReaderFormat bppoReader(1);
  unique_ptr<IAlignment> iAln(bppoReader.read(sequenceFormat));
  map<string, string> args(bppoReader.getUnparsedArguments());
  ApplicationTools::displayResult("character data file ", charDataFilePath);
  ApplicationTools::displayResult("chatacter data format ", iAln->getFormatName());
  const SequenceContainer *charCont = iAln->readAlignment(charDataFilePath, alphabet);
  VectorSiteContainer *sites = new VectorSiteContainer(*dynamic_cast<const OrderedSequenceContainer *>(charCont));
  delete charCont;
  return sites;
}

/******************************************************************************/

void giveNamesToInternalNodes(Tree *tree)
{
  TreeTemplate<Node> *ttree = dynamic_cast<TreeTemplate<Node> *>(tree);
  vector<Node *> nodes = ttree->getNodes();
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    if (!nodes[i]->hasName())
      nodes[i]->setName("_baseInternal_" + TextTools::toString(nodes[i]->getId()));
  }
}

/******************************************************************************/

Tree *processTree(BppApplication *bppml)
{
  Tree *tree = PhylogeneticsApplicationTools::getTree(bppml->getParams());
  giveNamesToInternalNodes(tree);

  string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", bppml->getParams(), "Input", "", true, 1); // process tree branches lengths
  string cmdName;
  map<string, string> cmdArgs;
  KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs); // this line process cmdName - which dictates weather the tree should be processed as is, or ultrameterized
  if (cmdName == "Clock")                                         // if needed, turn the tree into an ultrametric tree
  {
    TreeTools::convertToClockTree(*tree, tree->getRootId(), true);
  }
  return tree;
}

/******************************************************************************/

TransitionModel *setCharacterModel(BppApplication *bppml, VectorSiteContainer *charData, const BinaryAlphabet *alphabet, Tree *tree)
{
  // create the model
  // TO DO: ADD HERE PROCESSING OF INITIAL CHARACTER MODEL PARAMETERS AND PADD TO CONSTRUCTOR
  // extract the user initial value of k for potential later use
  double init_mu = ApplicationTools::getDoubleParameter("character_model.mu", bppml->getParams(), 1);
  double init_pi0 = ApplicationTools::getDoubleParameter("character_model.pi0", bppml->getParams(), 0.5);
  SubstitutionModel *model = new TwoParameterBinarySubstitutionModel(alphabet, init_mu, init_pi0);

  // compute the maximum parsimony score and set the lower and upper bounds on mu (=rate) as mp/tree_size, 2*mp/tree_size
  VDouble treeBranches = dynamic_cast<TreeTemplate<Node> *>(tree)->getBranchLengths();
  double treeSize = 0;
  for (size_t i = 0; i < treeBranches.size(); ++i)
  {
    treeSize += treeBranches[i];
  }

  model->setParameterValue(string("mu"), init_mu);
  map<int, double> frequencies;
  frequencies[0] = init_pi0;
  frequencies[1] = 1 - init_pi0;
  model->setFreq(frequencies);
 
  return dynamic_cast<TransitionModel *>(model);
}

/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/

int main(int args, char **argv)
{
  cout << "************************************************************************************************" << endl;
  cout << "*       Running character model optmization twice to test convergence in BrentOneDimension       *" << endl;
  cout << "************************************************************************************************" << endl;
  cout << endl;

  try
  {
	  
	// set random seed
	double seed = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
	cout << "seed=" << seed << endl;
	RandomTools::setSeed(static_cast<long>(seed));

    /* process input from params file */
    BppApplication bppml(args, argv, "bppml");

    /* process character data */
    const BinaryAlphabet *balpha = new BinaryAlphabet();
    VectorSiteContainer *charData = processCharacterData(&bppml, balpha);

    /* process tree */
    Tree *tree = processTree(&bppml);
    vector<Node *> nodes = (dynamic_cast<TreeTemplate<Node> *>(tree))->getNodes();

    /* set the character model */
    TransitionModel *charModel = setCharacterModel(&bppml, charData, balpha, tree);
	
	/* set the tree likelihood */
	DiscreteDistribution *rDist = new ConstantRateDistribution();
	RHomogeneousTreeLikelihood* tl = new RHomogeneousTreeLikelihood(*tree, *charData, charModel, rDist, false);
    tl->initialize();

    /* fit the null model: separate optimization of the character model and the sequence model, where the selection intensity parameter k is 1 */
    cout << "\n**** model fitting 1 ****" << endl;
	/*BrentOneDimension* characterParametersOptimizer = new BrentOneDimension(tl);
	characterParametersOptimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);
	characterParametersOptimizer->getStopCondition()->setTolerance(0.000001); // set the tolerance to be slighly less strict to account for the instability of the joint likelihood function
	characterParametersOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
	characterParametersOptimizer->setProfiler(0);
	characterParametersOptimizer->setMessageHandler(0);
	characterParametersOptimizer->setVerbose(1);
	ParameterList pi0, mu;
	mu.addParameter(tl->getParameter("TwoParameterBinary.mu"));
	const IntervalConstraint* muBounds = dynamic_cast<const IntervalConstraint*>(tl->getParameter("TwoParameterBinary.mu").getConstraint());
	pi0.addParameter(tl->getParameter("TwoParameterBinary.pi0"));
	const IntervalConstraint* pi0Bounds = dynamic_cast<const IntervalConstraint*>(tl->getParameter("TwoParameterBinary.pi0").getConstraint()); 
	
	// optimize the joint model with respect to pi0
	characterParametersOptimizer->setInitialInterval(pi0Bounds->getLowerBound(), pi0Bounds->getUpperBound()); // search within stricter bounds that the actual ones of pi0 to avoid failute of stochasitc mapping
	characterParametersOptimizer->init(pi0);
	characterParametersOptimizer->optimize();

	// optimize the joint model with respect to mu
	characterParametersOptimizer->setInitialInterval(muBounds->getLowerBound(), muBounds->getUpperBound()); // search within stricter bounds that the actual ones of pi0 to avoid failute of stochasitc mapping
	characterParametersOptimizer->init(mu);
	characterParametersOptimizer->optimize(); */

	PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams());
	
    ParameterList parameters = tl->getParameters();
	for (size_t p=0; p<parameters.size(); ++p)
	{
		ApplicationTools::displayResult(parameters[p].getName(), TextTools::toString(parameters[p].getValue()));
	}
	ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl->getValue()));

    /* fit the alternative model: sequencial optimization of the character model and then sequence model, given an expected history based on the character model */
    cout << "\n**** model fitting 2****" << endl;
	/*// optimize the joint model with respect to pi0
	characterParametersOptimizer->setInitialInterval(pi0Bounds->getLowerBound(), pi0Bounds->getUpperBound()); // search within stricter bounds that the actual ones of pi0 to avoid failute of stochasitc mapping
	characterParametersOptimizer->init(pi0);
	characterParametersOptimizer->optimize();

	// optimize the joint model with respect to mu
	characterParametersOptimizer->setInitialInterval(muBounds->getLowerBound(), muBounds->getUpperBound()); // search within stricter bounds that the actual ones of pi0 to avoid failute of stochasitc mapping
	characterParametersOptimizer->init(mu);
	characterParametersOptimizer->optimize(); */
	
	PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams());
	
	parameters = tl->getParameters();
	for (size_t p=0; p<parameters.size(); ++p)
	{
		ApplicationTools::displayResult(parameters[p].getName(), TextTools::toString(parameters[p].getValue()));
	}
	ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl->getValue()));

	
    // free parameters
    delete balpha;
    delete rDist;
    delete charData;
    delete charModel;
    delete tree;

    bppml.done();
  }
  catch (exception &e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}
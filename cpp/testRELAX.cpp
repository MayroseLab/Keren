// unit tests for RELAX

// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>

using namespace std;

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

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

/******************************************************************************/
/**************************** Auxiliary functions *****************************/
/******************************************************************************/


// if the initial likelihood is 0, either a parameter or some of the branch lengths are smaller than 0.000001
// to resolve this, change their values to 0.000001
DiscreteRatesAcrossSitesTreeLikelihood* correct_for_underflow(DiscreteRatesAcrossSitesTreeLikelihood* tl)
{
    double logL = tl->getValue();
    if (std::isinf(logL))
    {
      ParameterList pl = tl->getBranchLengthsParameters();
      for (unsigned int i = 0; i < pl.size(); i++)
      {
        if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
      }
      tl->matchParametersValues(pl); // load the new parameters values to the likelihood instance
    }
    return(tl);
}


void printModelParameters(DiscreteRatesAcrossSitesTreeLikelihood* tl)
  {
    double logL = tl->getValue();
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-logL, 15));
    ParameterList parameters = tl->getSubstitutionModelParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
    parameters = tl->getRateDistributionParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
    ParameterList pl = tl->getBranchLengthsParameters();
    for(unsigned int i = 0; i < pl.size(); i++)
    {
       ApplicationTools::displayResult(pl[i].getName(), TextTools::toString(pl[i].getValue())); 
    }
}

// function that optimizes the model and then prints to the screen the likelihood and the MLEs
DiscreteRatesAcrossSitesTreeLikelihood* optimize(DiscreteRatesAcrossSitesTreeLikelihood* tl, BppApplication bppml) 
{
    tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(
    PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams()));
    
    // write the optimized parameters to screen
    printModelParameters(tl);
    
    // check if the optimization procedure acheived valid parameters
    PhylogeneticsApplicationTools::checkEstimatedParameters(tl->getParameters());
    return(tl);

}

void doLRT(DiscreteRatesAcrossSitesTreeLikelihood* null_tl, DiscreteRatesAcrossSitesTreeLikelihood* alternative_tl)
{
  // get the logl and parameters number of the null model
  double null_logL = null_tl->getValue();
  unsigned int null_params_num = static_cast<unsigned int>(null_tl->getSubstitutionModelParameters().size() + null_tl->getRateDistributionParameters().size()); // the number of parameters does not include the sclaing parameter
  // get the logl and parameters number of the alternative model
  double alternative_logL = alternative_tl->getValue();
  unsigned int alternative_params_num = static_cast<unsigned int>(alternative_tl->getSubstitutionModelParameters().size() + alternative_tl->getRateDistributionParameters().size()); // the number of parameters does not include the sclaing parameter
  // compute the LRT and its p-value, given the chi-2 distribution with the appropriate degrees of freedom
  double stat = 2*alternative_logL - 2*null_logL;
  double pvalue = 1. - RandomTools::pChisq(stat, alternative_params_num - null_params_num);
  ApplicationTools::displayResult("LRT statistic", TextTools::toString(stat, 15));
  ApplicationTools::displayResult("LRT p-value", TextTools::toString(pvalue, 15));
}

VectorSiteContainer* process_alignment(Alphabet* alphabet, BppApplication bppml)
{
    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppml.getParams()); // here, gaps will be converted to unknown characters
    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppml.getParams(), "", true, false); // convert the alignemnt to condensed format of unique sites
    delete allSites; // delete the non-condenced intance of the sequence data
    return(sites);

}

Tree* process_tree(BppApplication bppml)
{
    Tree* tree = PhylogeneticsApplicationTools::getTree(bppml.getParams());
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



/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/

int main(int args, char** argv)
{
  cout << "**************************************************" << endl;
  cout << "*       Unit tests for RELAX implementation      *" << endl;
  cout << "**************************************************" << endl;
  cout << endl;

  try
  { 
    /* process input */
    
    // process input from params file
    BppApplication relaxTest(args, argv, "RELAXTest");
    relaxTest.startTimer();

    // process alphabet type
    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(relaxTest.getParams(), "", false);
    unique_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", relaxTest.getParams(), "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    }

    // process alignment
    VectorSiteContainer* sites = process_alignment(alphabet, relaxTest);

    // process tree
    Tree* tree = process_tree(relaxTest);

    // set the branch-site model
    SubstitutionModelSet* modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), sites, relaxTest.getParams());
    DiscreteDistribution* bs_rDist = PhylogeneticsApplicationTools::getRateDistribution(relaxTest.getParams());
    DiscreteRatesAcrossSitesTreeLikelihood* bs_tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, dynamic_cast<MixedSubstitutionModelSet*>(modelSet), bs_rDist, true, true);
    bs_tl->initialize();

    // compute the likelihood
    cout << "\nCompute log likelihood" << endl;
    printModelParameters(bs_tl); // debug - print model parameters

    /* free resources */
    delete alphabet;
    delete sites;
    if (modelSet) delete modelSet;
    delete bs_rDist;
    delete bs_tl;
    delete tree;
    relaxTest.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  return 0;
}